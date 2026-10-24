/* =======================================================================
Copyright 2026 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Threading;

namespace ilPSP.Utils {

    /// <summary>
    /// Implementations which can serve the ordinary BLAS calls, see <see cref="BLAS.Provider"/>.
    /// LAPACK uses the same provider and threading mode; PARDISO remains independent.
    /// </summary>
    public enum BLASProvider {
        /// <summary>Intel MKL (`libBoSSSnative_seq.so` resp. `libBoSSSnative_omp.so`), the baseline</summary>
        MKL = 0,

        /// <summary>OpenBLAS, bundled with the PETSc native wrapper, optional</summary>
        OpenBLAS = 1
    }

    static public partial class BLAS {

        // ==========================================================================================
        // provider state
        // ==========================================================================================

        static readonly object s_ProviderLock = new object();

        static BLASProvider s_Provider = BLASProvider.MKL;

        /// <summary>
        /// true, if the last call to <see cref="ActivateOMP"/>/<see cref="ActivateSEQ"/> selected threaded BLAS.
        /// </summary>
        static bool s_ThreadedMode = true; // consistent with the static ctor (m_BLAS = m_omp_BLAS)

        /// <summary>
        /// the implementation which serves the ordinary BLAS calls (<see cref="F77_BLAS"/>, <see cref="dgemm"/>, ...);
        /// can be changed by <see cref="TryActivateProvider"/>, independent of the sequential/threaded mode.
        /// </summary>
        public static BLASProvider Provider {
            get {
                return s_Provider;
            }
        }

        /// <summary>
        /// human-readable explanation why <see cref="Provider"/> was selected
        /// (explicit request, startup calibration result, fallback reason);
        /// </summary>
        public static string ProviderSelectionReason {
            get;
            private set;
        } = "default (startup selection not executed)";

        /// <summary>
        /// The provider requested via environment variable `BOSSS_BLAS_PROVIDER` (`mkl` (default), `auto` or `openblas`);
        /// null if the startup selection was not executed yet.
        /// </summary>
        public static string RequestedProvider {
            get;
            private set;
        }

        /// <summary>
        /// Selects the implementation for the ordinary BLAS calls;
        /// only at a quiescent boundary, i.e., not within a parallel section (<see cref="Environment.InParallelSection"/>).
        /// The sequential/threaded mode is preserved.
        /// </summary>
        /// <returns>
        /// false, if the provider is unavailable or cannot be switched now; then, the previous provider remains active
        /// and <paramref name="reason"/> explains why.
        /// </returns>
        public static bool TryActivateProvider(BLASProvider p, out string reason) {
            if(p != BLASProvider.MKL && p != BLASProvider.OpenBLAS)
                throw new ArgumentOutOfRangeException(nameof(p));
            if(ilPSP.Environment.InParallelSection) {
                reason = "provider cannot be changed within a parallel section";
                return false;
            }
            if(p == BLASProvider.OpenBLAS && !OpenBLASAvailable) {
                reason = "OpenBLAS unavailable: " + OpenBLASUnavailableReason;
                return false;
            }
            lock(s_ProviderLock) {
                s_Provider = p;
                ApplyProviderAndMode();
            }
            reason = null;
            return true;
        }

        /// <summary>
        /// Use Intel MKL for the ordinary BLAS calls (always available)
        /// </summary>
        public static void ActivateMKL() {
            if(!TryActivateProvider(BLASProvider.MKL, out string reason))
                throw new NotSupportedException(reason);
            ProviderSelectionReason = "explicit API call";
        }

        /// <summary>
        /// Use OpenBLAS for the ordinary BLAS calls, if available.
        /// </summary>
        /// <returns>false, if OpenBLAS could not be activated (<paramref name="reason"/>); then, the previous provider remains active</returns>
        public static bool ActivateOpenBLAS(out string reason) {
            bool r = TryActivateProvider(BLASProvider.OpenBLAS, out reason);
            if(r)
                ProviderSelectionReason = "explicit API call";
            return r;
        }

        static void SetThreadingMode(bool threaded) {
            lock(s_ProviderLock) {
                s_ThreadedMode = threaded;
                ApplyProviderAndMode();
            }
        }

        /// <summary>
        /// Number of threads that threaded OpenBLAS may use; the same budget as MKL, see <see cref="Environment.PinOMPthreads"/>;
        /// note that <see cref="Environment.NumThreads"/> is already clamped to the CPUs which are dedicated to this rank.
        /// </summary>
        static int OpenBLASThreadBudget {
            get {
                return Math.Max(1, Math.Min(ilPSP.Environment.MaxNumOpenMPthreads, ilPSP.Environment.NumThreads));
            }
        }

        /// <summary>
        /// must be called within <see cref="s_ProviderLock"/>
        /// </summary>
        internal static UnsafeLAPACK ActiveLAPACK => s_Provider == BLASProvider.OpenBLAS
            ? s_OpenBLAS.Value.Lapack
            : (s_ThreadedMode ? LAPACK.omp_F77_LAPACK : LAPACK.seq_F77_LAPACK);

        static void ApplyProviderAndMode() {
            if(s_Provider == BLASProvider.OpenBLAS) {
                var ob = s_OpenBLAS.Value;
                SetOpenBLASThreads(s_ThreadedMode ? OpenBLASThreadBudget : 1);
                m_BLAS = ob.Blas;
                s_ActiveOpenBLAS = ob.Blas;
            } else {
                s_ActiveOpenBLAS = null;
                m_BLAS = s_ThreadedMode ? m_omp_BLAS : m_seq_BLAS;
            }
        }

        // ==========================================================================================
        // OpenBLAS process-global thread count
        // ==========================================================================================
        //
        // OpenBLAS has only one, process-global thread count, i.e., there is no separate sequential library like for MKL.
        // * The count is only changed in `ApplyProviderAndMode`, i.e., at the coarse boundaries of parallel sections
        //   (`Environment.DisableOpenMP`/`EnableOpenMP`), resp. by provider switches, never per call.
        // * Within a parallel section the count is 1; it is raised only after the workers have joined.
        // * A tensor kernel (`SeqDGEMM`/`SeqDGEMV`) outside of a parallel section uses OpenBLAS only if the count is 1;
        //   it holds a lease (`s_SeqLeases`) for the duration of the call, and the count is never raised while a lease is held.
        //   Otherwise, sequential MKL is used: a tensor kernel never runs threaded BLAS.

        /// <summary>non-null iff OpenBLAS is the active provider</summary>
        static volatile UnsafeDBLAS s_ActiveOpenBLAS;

        /// <summary>current OpenBLAS thread count as set by us; 0 if never set</summary>
        static volatile int s_OpenBLASThreads = 0;

        /// <summary>set while the thread count is about to be raised; then, new tensor calls use MKL</summary>
        static volatile bool s_OpenBLASRaising = false;

        static int s_SeqLeases = 0;

        static long s_OpenBLASThreadChanges = 0;

        /// <summary>
        /// number of changes of the OpenBLAS thread count, for diagnostics
        /// </summary>
        public static long OpenBLASThreadCountChanges => Interlocked.Read(ref s_OpenBLASThreadChanges);

        /// <summary>
        /// current OpenBLAS thread count, as reported by OpenBLAS, or -1 if OpenBLAS is not loaded
        /// </summary>
        public static int OpenBLASNumThreads {
            get {
                if(!s_OpenBLAS.IsValueCreated || s_OpenBLAS.Value == null)
                    return -1;
                return s_OpenBLAS.Value.Control.GetNumThreads();
            }
        }

        static void SetOpenBLASThreads(int n) {
            if(n == s_OpenBLASThreads)
                return; // also, a redundant `set` would re-create a shut-down thread pool
            var ctrl = s_OpenBLAS.Value.Control;
            if(n > s_OpenBLASThreads && s_OpenBLASThreads != 0) {
                s_OpenBLASRaising = true;
                Interlocked.MemoryBarrier();
                var sw = new SpinWait();
                while(Volatile.Read(ref s_SeqLeases) > 0)
                    sw.SpinOnce();
                ctrl.SetNumThreads(n);
                s_OpenBLASThreads = n;
                s_OpenBLASRaising = false;
            } else {
                ctrl.SetNumThreads(n);
                s_OpenBLASThreads = n;
            }
            Interlocked.Increment(ref s_OpenBLASThreadChanges);
        }

        /// <summary>
        /// Returns the OpenBLAS instance if it may be used single-threaded now; then, <paramref name="leased"/> must be released by <see cref="ReleaseSeqLease"/>.
        /// </summary>
        static UnsafeDBLAS AcquireSeqOpenBLAS(out bool leased) {
            var ob = s_ActiveOpenBLAS;
            leased = false;
            if(ob == null)
                return null;
            if(ilPSP.Environment.InParallelSection) {
                // count is 1 for the entire section; it is raised only after the workers have joined
                return s_OpenBLASThreads == 1 ? ob : null;
            }
            Interlocked.Increment(ref s_SeqLeases);
            if(s_OpenBLASThreads == 1 && !s_OpenBLASRaising && s_ActiveOpenBLAS != null) {
                leased = true;
                return ob;
            }
            Interlocked.Decrement(ref s_SeqLeases);
            return null;
        }

        static void ReleaseSeqLease(bool leased) {
            if(leased)
                Interlocked.Decrement(ref s_SeqLeases);
        }

        /// <summary>
        /// Single-threaded DGEMM for tensor kernels (e.g. <see cref="MultidimensionalArray.Multiply(double, MultidimensionalArray, MultidimensionalArray, double, string)"/>);
        /// uses the active <see cref="Provider"/> if it can be guaranteed to run sequentially, otherwise sequential MKL.
        /// </summary>
        unsafe static public void SeqDGEMM(ref int TRANSA, ref int TRANSB,
                                           ref int M, ref int N, ref int K,
                                           ref double ALPHA,
                                           double* A, ref int LDA,
                                           double* B, ref int LDB,
                                           ref double BETA,
                                           double* C, ref int LDC) {
            var ob = AcquireSeqOpenBLAS(out bool leased);
            if(ob != null) {
                try {
                    ob.DGEMM(ref TRANSA, ref TRANSB, ref M, ref N, ref K, ref ALPHA, A, ref LDA, B, ref LDB, ref BETA, C, ref LDC);
                } finally {
                    ReleaseSeqLease(leased);
                }
            } else {
                m_seq_BLAS.DGEMM(ref TRANSA, ref TRANSB, ref M, ref N, ref K, ref ALPHA, A, ref LDA, B, ref LDB, ref BETA, C, ref LDC);
            }
        }

        /// <summary>
        /// Single-threaded DGEMV for tensor kernels, see <see cref="SeqDGEMM"/>.
        /// </summary>
        unsafe static public void SeqDGEMV(ref int TRANSA,
                                           ref int M, ref int N,
                                           ref double ALPHA,
                                           double* A, ref int LDA,
                                           double* X, ref int INCX,
                                           ref double BETA,
                                           double* Y, ref int INCY) {
            var ob = AcquireSeqOpenBLAS(out bool leased);
            if(ob != null) {
                try {
                    ob.DGEMV(ref TRANSA, ref M, ref N, ref ALPHA, A, ref LDA, X, ref INCX, ref BETA, Y, ref INCY);
                } finally {
                    ReleaseSeqLease(leased);
                }
            } else {
                m_seq_BLAS.DGEMV(ref TRANSA, ref M, ref N, ref ALPHA, A, ref LDA, X, ref INCX, ref BETA, Y, ref INCY);
            }
        }

    }
}

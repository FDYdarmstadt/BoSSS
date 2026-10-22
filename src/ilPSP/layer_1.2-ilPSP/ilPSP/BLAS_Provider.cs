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
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using ilPSP.Tracing;
using MPI.Wrappers.Utils;

namespace ilPSP.Utils {

    /// <summary>
    /// Implementations which can serve the ordinary BLAS calls, see <see cref="BLAS.Provider"/>.
    /// LAPACK, PARDISO etc. always use Intel MKL, independent of this choice.
    /// </summary>
    public enum BLASProvider {
        /// <summary>Intel MKL (`libBoSSSnative_seq.so` resp. `libBoSSSnative_omp.so`), the baseline</summary>
        MKL = 0,

        /// <summary>OpenBLAS, bundled with the PETSc native wrapper, optional</summary>
        OpenBLAS = 1
    }

    /// <summary>
    /// Thread controls of the OpenBLAS module;
    /// the prototypes were verified by disassembly of `libBoSSSnative_petsc_mpi_omp.so` (2026-10-22):
    /// the wrappers tail-call the OpenBLAS C-API `openblas_set_num_threads(int)`, `openblas_get_num_threads()`,
    /// `openblas_get_parallel()`, `blas_thread_shutdown_()`; `BoSSS_openblas_get_config(char*, int)` is a `snprintf` of corename and config.
    /// </summary>
    unsafe sealed class OpenBLASControl {

        readonly delegate* unmanaged<int, void> m_set_num_threads;
        readonly delegate* unmanaged<int> m_get_num_threads;
        readonly delegate* unmanaged<int> m_get_parallel;
        readonly delegate* unmanaged<byte*, int, void> m_get_config;

        OpenBLASControl(IntPtr set, IntPtr getN, IntPtr getPar, IntPtr getCfg) {
            m_set_num_threads = (delegate* unmanaged<int, void>)set;
            m_get_num_threads = (delegate* unmanaged<int>)getN;
            m_get_parallel = (delegate* unmanaged<int>)getPar;
            m_get_config = (delegate* unmanaged<byte*, int, void>)getCfg;
        }

        /// <summary>
        /// resolves the controls on the module handle of <paramref name="lib"/>; null if any is missing.
        /// </summary>
        internal static OpenBLASControl TryCreate(DynLibLoader lib, out string reason) {
            IntPtr Get(string name, ref string missing) {
                IntPtr p = DynamicLibraries.LoadSymbol(lib.LibHandle, name, out string err);
                if(p == IntPtr.Zero)
                    missing = (missing == null ? "" : missing + ", ") + name;
                return p;
            }
            string missing = null;
            var set = Get("BoSSS_openblas_set_num_threads", ref missing);
            var getN = Get("BoSSS_openblas_get_num_threads", ref missing);
            var getPar = Get("BoSSS_openblas_get_parallel", ref missing);
            var getCfg = Get("BoSSS_openblas_get_config", ref missing);
            if(missing != null) {
                reason = $"missing OpenBLAS thread-control export(s) {missing} in {lib.CurrentLibraryName}";
                return null;
            }
            reason = null;
            return new OpenBLASControl(set, getN, getPar, getCfg);
        }

        internal void SetNumThreads(int n) => m_set_num_threads(n);
        internal int GetNumThreads() => m_get_num_threads();

        /// <summary> 0: sequential, 1: pthreads, 2: OpenMP </summary>
        internal int GetParallel() => m_get_parallel();

        internal string GetConfig() {
            byte* buf = stackalloc byte[512];
            buf[0] = 0;
            m_get_config(buf, 512);
            return Marshal.PtrToStringAnsi((IntPtr)buf);
        }
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

        // ==========================================================================================
        // OpenBLAS discovery
        // ==========================================================================================

        sealed class OpenBLASModule {
            internal UnsafeDBLAS Blas;
            internal OpenBLASControl Control;
            internal string Config;
            internal string LoadInfo;
        }

        static string s_OpenBLASUnavailableReason = null;

        /// <summary>
        /// Lazy discovery; the module is never unloaded, so delegates remain valid for the lifetime of the process.
        /// </summary>
        static readonly Lazy<OpenBLASModule> s_OpenBLAS = new Lazy<OpenBLASModule>(DiscoverOpenBLAS, LazyThreadSafetyMode.ExecutionAndPublication);

        /// <summary>
        /// true, if OpenBLAS can be loaded (triggers the discovery on first call)
        /// </summary>
        public static bool OpenBLASAvailable {
            get {
                return s_OpenBLAS.Value != null;
            }
        }

        /// <summary>
        /// null if <see cref="OpenBLASAvailable"/>, otherwise the reason (missing library, dependency, symbol or platform)
        /// </summary>
        public static string OpenBLASUnavailableReason {
            get {
                _ = s_OpenBLAS.Value;
                return s_OpenBLASUnavailableReason;
            }
        }

        /// <summary>
        /// OpenBLAS core name and build configuration, null if unavailable
        /// </summary>
        public static string OpenBLASConfig {
            get {
                return s_OpenBLAS.Value?.Config;
            }
        }

        /// <summary>
        /// library path and environment settings applied when loading OpenBLAS, null if unavailable
        /// </summary>
        public static string OpenBLASLoadInfo {
            get {
                return s_OpenBLAS.Value?.LoadInfo;
            }
        }

        /// <summary>
        /// Full path of the loaded library for an MKL wrapper
        /// </summary>
        public static string MKLLibraryPath => m_seq_BLAS.CurrentLibraryName;

        [DllImport("libc.so.6", EntryPoint = "getenv")]
        static extern IntPtr libc_getenv(string name);

        [DllImport("libc.so.6", EntryPoint = "setenv")]
        static extern int libc_setenv(string name, string value, int overwrite);

        /// <summary>
        /// Sets a variable in the **native** environment (on Linux, <see cref="System.Environment.SetEnvironmentVariable(string, string)"/> does not reach `getenv` of native libraries),
        /// unless it is already defined there.
        /// </summary>
        /// <returns>the effective value</returns>
        static string SetNativeEnvDefault(string name, string value) {
            IntPtr cur = libc_getenv(name);
            if(cur != IntPtr.Zero)
                return Marshal.PtrToStringAnsi(cur);
            if(libc_setenv(name, value, 0) != 0)
                return null;
            return Marshal.PtrToStringAnsi(libc_getenv(name));
        }

        [DllImport("libc.so.6", EntryPoint = "dlopen")]
        static extern IntPtr libc_dlopen(string filename, int flags);

        [DllImport("libc.so.6", EntryPoint = "dlerror")]
        static extern IntPtr libc_dlerror();

        /// <summary>
        /// Human-readable reason why <paramref name="libName"/> could not be loaded:
        /// not found in the search path, missing dependency, or missing symbol.
        /// (The error string recorded by <see cref="DynLibLoader"/> is often empty, since the runtime may call `dlerror` in between.)
        /// </summary>
        static string DiagnoseLoadFailure(string libName) {
            var dirs = new List<string>();
            if(!string.IsNullOrWhiteSpace(DynLibLoader.PrimaryLibrarySearchPath))
                dirs.Add(DynLibLoader.PrimaryLibrarySearchPath);
            foreach(var d in (System.Environment.GetEnvironmentVariable("LD_LIBRARY_PATH") ?? "").Split(':', StringSplitOptions.RemoveEmptyEntries))
                if(!dirs.Contains(d))
                    dirs.Add(d);
            string path = dirs.Select(d => System.IO.Path.Combine(d, libName)).FirstOrDefault(System.IO.File.Exists);
            if(path == null)
                return $"{libName} not found in {dirs.ToConcatString("", ":", "")}";

            string err = null;
            try {
                libc_dlerror(); // resolve the P/Invoke and clear a pending error, before the dlopen
                IntPtr h = libc_dlopen(path, 2); // RTLD_NOW
                if(h != IntPtr.Zero) {
                    DynamicLibraries.dlclose(new DynamicLibraries.DynLibHandle() { val = h });
                    return $"{path} can be loaded, but some BLAS symbol is missing (incompatible or outdated native library)";
                }
                IntPtr e = libc_dlerror();
                err = e != IntPtr.Zero ? Marshal.PtrToStringAnsi(e) : null;
            } catch(Exception ex) when(ex is DllNotFoundException || ex is EntryPointNotFoundException) {
                // glibc < 2.34: dlopen is in libdl
            }
            return $"{path} cannot be loaded (missing dependency, e.g. MPI?): " + (err ?? "no error information");
        }

        /// <summary>
        /// Default for `OPENBLAS_THREAD_TIMEOUT` (OpenBLAS default is 28, i.e. idle workers spin ~2^28 cycles):
        /// lower values avoid that spinning workers slow down subsequent parallel sections (native handoff measurement).
        /// </summary>
        const string OpenBLASThreadTimeoutDefault = "14";

        static OpenBLASModule DiscoverOpenBLAS() {
            using(var tr = new FuncTrace()) {
                try {
                    if(System.Environment.OSVersion.Platform != PlatformID.Unix) {
                        s_OpenBLASUnavailableReason = "OpenBLAS is only supported on Linux";
                        return null;
                    }

                    // OpenBLAS reads its environment when it is loaded:
                    var info = new StringBuilder();
                    IntPtr already = DynamicLibraries.dlopen("libopenblas.so.0", 2 | 4); // RTLD_NOW | RTLD_NOLOAD
                    if(already != IntPtr.Zero) {
                        DynamicLibraries.dlclose(new DynamicLibraries.DynLibHandle() { val = already });
                        info.Append("OpenBLAS was already loaded (e.g. by PETSc), environment settings could not be applied; ");
                    } else {
                        string timeout = SetNativeEnvDefault("OPENBLAS_THREAD_TIMEOUT", OpenBLASThreadTimeoutDefault);
                        // avoids that OpenBLAS creates a thread for every CPU in the mask at load time; the count can be raised later
                        string numThreads = SetNativeEnvDefault("OPENBLAS_NUM_THREADS", OpenBLASThreadBudget.ToString(CultureInfo.InvariantCulture));
                        info.Append($"OPENBLAS_THREAD_TIMEOUT={timeout}, OPENBLAS_NUM_THREADS={numThreads}; ");
                    }

                    UnsafeDBLAS blas;
                    try {
                        blas = new UnsafeDBLAS(BLASProvider.OpenBLAS);
                    } catch(ApplicationException) {
                        // missing library, missing dependency (e.g. libmpi.so.40) or missing BLAS symbol
                        s_OpenBLASUnavailableReason = DiagnoseLoadFailure(BLAS_OpenBLAS_Libstuff.LinuxLibName);
                        tr.Info("OpenBLAS unavailable: " + s_OpenBLASUnavailableReason);
                        return null;
                    }

                    var ctrl = OpenBLASControl.TryCreate(blas, out string reason);
                    if(ctrl == null) {
                        // without thread control, single-thread execution in parallel sections cannot be guaranteed;
                        // the module stays loaded (it might be used by PETSc).
                        s_OpenBLASUnavailableReason = reason;
                        tr.Info("OpenBLAS unavailable: " + reason);
                        return null;
                    }

                    var r = new OpenBLASModule() {
                        Blas = blas,
                        Control = ctrl,
                        Config = ctrl.GetConfig() + $", threading model {ctrl.GetParallel()} (0=seq, 1=pthreads, 2=OpenMP)",
                    };
                    info.Insert(0, blas.CurrentLibraryName + "; ");
                    info.Append($"initial thread count {ctrl.GetNumThreads()}");
                    r.LoadInfo = info.ToString();
                    tr.Info("OpenBLAS loaded: " + r.LoadInfo + "; " + r.Config);
                    return r;
                } catch(DllNotFoundException e) {
                    s_OpenBLASUnavailableReason = "libc/libdl not found: " + e.Message;
                    return null;
                } catch(EntryPointNotFoundException e) {
                    s_OpenBLASUnavailableReason = "libc/libdl entry missing: " + e.Message;
                    return null;
                }
            }
        }

        // ==========================================================================================
        // startup selection
        // ==========================================================================================

        /// <summary>
        /// Environment variable for the choice of the BLAS provider: `mkl` (default), `auto` (calibration at startup) or `openblas`
        /// </summary>
        public const string ProviderEnvVar = "BOSSS_BLAS_PROVIDER";

        static bool s_StartupSelectionDone = false;

        /// <summary>
        /// timings of the last startup calibration; null if not executed
        /// </summary>
        public static ProviderCalibrationResult LastCalibration {
            get;
            private set;
        }

        /// <summary>
        /// Selects the BLAS provider, once per process:
        /// must be called after native library search paths and the CPU allocation (<see cref="Environment.InitThreading"/>) are set up,
        /// and before parallel numerical work.
        /// According to <see cref="ProviderEnvVar"/>: MKL by default (on AMD EPYC Zen 4, MKL with the fake-Intel probe outperformed OpenBLAS 0.3.34);
        /// for `auto`, a short calibration on representative shapes, which keeps MKL for ties, incomplete or inconclusive measurements and if OpenBLAS is unavailable.
        /// In the default case, the OpenBLAS library is not even loaded.
        /// </summary>
        public static void SelectProviderAtStartup() {
            if(s_StartupSelectionDone)
                return;
            using(var tr = new FuncTrace()) {
                s_StartupSelectionDone = true;
                string envVal = System.Environment.GetEnvironmentVariable(ProviderEnvVar);
                string req = (envVal ?? "mkl").Trim().ToLowerInvariant();
                if(req.Length == 0) {
                    req = "mkl";
                    envVal = null;
                }
                RequestedProvider = req;

                switch(req) {
                    case "mkl":
                    ProviderSelectionReason = envVal == null ? $"default ({ProviderEnvVar} not set)" : $"explicit request ({ProviderEnvVar}=mkl)";
                    break;

                    case "openblas":
                    if(TryActivateProvider(BLASProvider.OpenBLAS, out string reason)) {
                        ProviderSelectionReason = $"explicit request ({ProviderEnvVar}=openblas)";
                    } else {
                        ProviderSelectionReason = $"explicit request ({ProviderEnvVar}=openblas) failed, using MKL: {reason}";
                        tr.Warning(ProviderSelectionReason);
                    }
                    break;

                    case "auto":
                    if(!OpenBLASAvailable) {
                        ProviderSelectionReason = "auto: OpenBLAS unavailable, using MKL (" + OpenBLASUnavailableReason + ")";
                        break;
                    }
                    try {
                        var cal = CalibrateProviders();
                        LastCalibration = cal;
                        if(cal.OpenBLASWins) {
                            if(TryActivateProvider(BLASProvider.OpenBLAS, out string r2))
                                ProviderSelectionReason = "auto: " + cal.Summary;
                            else
                                ProviderSelectionReason = "auto: activation failed, using MKL: " + r2;
                        } else {
                            ProviderSelectionReason = "auto: " + cal.Summary;
                        }
                    } catch(Exception e) when(e is not ArithmeticException) {
                        tr.Warning("BLAS provider calibration failed: " + e);
                        ProviderSelectionReason = "auto: calibration failed, using MKL: " + e.Message;
                        lock(s_ProviderLock) {
                            s_Provider = BLASProvider.MKL;
                            ApplyProviderAndMode();
                        }
                    }
                    break;

                    default:
                    ProviderSelectionReason = $"unknown value {ProviderEnvVar}={req} (expected mkl, auto or openblas), using MKL";
                    tr.Warning(ProviderSelectionReason);
                    Console.Error.WriteLine("Warning: " + ProviderSelectionReason);
                    break;
                }
                tr.Info($"BLAS provider: {Provider}; {ProviderSelectionReason}");
            }
        }

        /// <summary>
        /// Result of <see cref="CalibrateProviders"/>
        /// </summary>
        public sealed class ProviderCalibrationResult {

            /// <summary>
            /// one measured case
            /// </summary>
            public sealed class Case {
                /// <summary>e.g. `dgemm nn 35x35x49`</summary>
                public string Name;
                /// <summary>seconds per call, best of several batches</summary>
                public double MKL, OpenBLAS;
                /// <summary>OpenBLAS/MKL; less than 1 means OpenBLAS is faster</summary>
                public double Ratio => OpenBLAS / MKL;
                /// <summary>maximum deviation, relative to the expected round-off</summary>
                public double ErrorRatio;
                /// <summary>threaded (true) or sequential (false) measurement</summary>
                public bool Threaded;
            }

            /// <summary> measured cases </summary>
            public List<Case> Cases = new List<Case>();

            /// <summary> geometric mean of the sequential <see cref="Case.Ratio"/> </summary>
            public double SeqGeoMeanRatio;

            /// <summary> worst sequential <see cref="Case.Ratio"/> </summary>
            public double SeqWorstRatio;

            /// <summary> threaded ratio, NaN if not measured </summary>
            public double ThreadedRatio = double.NaN;

            /// <summary> wall time of the calibration in seconds, including the first calls </summary>
            public double WallTime;

            /// <summary> decision </summary>
            public bool OpenBLASWins;

            /// <summary> decision and its reason </summary>
            public string Summary;

            /// <summary>
            /// table of all cases
            /// </summary>
            public override string ToString() {
                var stw = new StringBuilder();
                stw.AppendLine(Summary);
                foreach(var c in Cases)
                    stw.AppendLine(string.Format(CultureInfo.InvariantCulture, "  {0,-28} {1}  MKL {2,10:0.000} us  OpenBLAS {3,10:0.000} us  ratio {4:0.000}  err {5:0.00}",
                        c.Name, c.Threaded ? "thr" : "seq", c.MKL * 1e6, c.OpenBLAS * 1e6, c.Ratio, c.ErrorRatio));
                return stw.ToString();
            }
        }

        /// <summary>
        /// OpenBLAS is selected only if it is faster by at least this factor in the geometric mean over the sequential cases
        /// </summary>
        public const double CalibrationWinMargin = 1.10;

        /// <summary>
        /// OpenBLAS is rejected if any sequential case is slower than MKL by this factor (protects small shapes)
        /// </summary>
        public const double CalibrationWorstCaseLimit = 1.25;

        /// <summary>
        /// OpenBLAS is rejected if threaded OpenBLAS is slower than threaded MKL by this factor
        /// </summary>
        public const double CalibrationThreadedLimit = 1.10;

        /// <summary>
        /// Calibration is aborted (MKL is kept) if it exceeds this wall time (seconds)
        /// </summary>
        public const double CalibrationBudget = 0.5;

        /// <summary>
        /// Correctness-checks and times OpenBLAS against MKL on representative DGEMM/DGEMV shapes (from the Multiply census, `CENSUS-4.md`,
        /// plus some medium square cases), sequential and, if the thread budget is larger than one, threaded.
        /// Uses private buffers; the provider and threading mode are restored on exit.
        /// Must not be called within a parallel section.
        /// </summary>
        public static ProviderCalibrationResult CalibrateProviders() {
            if(ilPSP.Environment.InParallelSection)
                throw new NotSupportedException("calibration must not be called within a parallel section");
            if(!OpenBLASAvailable)
                throw new NotSupportedException("OpenBLAS unavailable: " + OpenBLASUnavailableReason);

            var ret = new ProviderCalibrationResult();
            var total = Stopwatch.StartNew();
            var ob = s_OpenBLAS.Value.Blas;

            // (M, N, K, transA, transB), Fortran order
            var gemmCases = new (int M, int N, int K, char ta, char tb)[] {
                (35, 35, 49, 'n', 'n'),    // LECEdgeQuadrature2, batch of cells
                (35, 35, 49, 't', 'n'),
                (35, 700, 1029, 'n', 'n'), // LECVolumeQuadrature2, 3D p4
                (56, 20, 1536, 'n', 't'),  // LECVolumeQuadrature2/NECQuadratureVolume
                (60, 35, 3, 'n', 'n'),     // short reduction
                (10, 10, 10, 'n', 'n'),    // tiny
                (20, 20, 20, 'n', 't'),
                (128, 128, 128, 'n', 'n'), // medium
                (300, 300, 300, 't', 'n'),
            };
            var gemvCases = new (int M, int N, char t)[] {
                (1029, 35, 't'),
                (35, 1029, 'n'),
                (200, 200, 't'),
                (16, 16, 't'),
            };

            int maxA = 0, maxB = 0, maxC = 0;
            foreach(var c in gemmCases) {
                maxA = Math.Max(maxA, c.M * c.K);
                maxB = Math.Max(maxB, c.K * c.N);
                maxC = Math.Max(maxC, c.M * c.N);
            }
            const int Nthr = 512;
            maxA = Math.Max(maxA, Nthr * Nthr); maxB = Math.Max(maxB, Nthr * Nthr); maxC = Math.Max(maxC, Nthr * Nthr);

            var rnd = new Random(4711);
            double[] A = new double[maxA], B = new double[maxB], C0 = new double[maxC], C1 = new double[maxC], C2 = new double[maxC];
            for(int i = 0; i < A.Length; i++) A[i] = rnd.NextDouble() - 0.5;
            for(int i = 0; i < B.Length; i++) B[i] = rnd.NextDouble() - 0.5;
            for(int i = 0; i < C0.Length; i++) C0[i] = rnd.NextDouble() - 0.5;

            bool budgetExceeded = false;
            string correctnessFailure = null;
            try {
                // sequential: OpenBLAS count 1, sequential MKL; 
                // the active provider (m_BLAS) is not changed, the OpenBLAS instance is called directly.
                lock(s_ProviderLock) {
                    SetOpenBLASThreads(1);
                }

                unsafe {
                    fixed(double* pA = A, pB = B, pC0 = C0, pC1 = C1, pC2 = C2) {
                        double* _pA = pA, _pB = pB, _pC0 = pC0, _pC1 = pC1, _pC2 = pC2;

                        foreach(var c in gemmCases) {
                            if(total.Elapsed.TotalSeconds > CalibrationBudget) {
                                budgetExceeded = true;
                                break;
                            }
                            int M = c.M, N = c.N, K = c.K, ta = c.ta, tb = c.tb;
                            int lda = c.ta == 'n' ? M : K, ldb = c.tb == 'n' ? K : N, ldc = M;
                            double alpha = 1.25, beta = 0.5;
                            int len = M * N;

                            void Run(UnsafeDBLAS lib, double* pC) {
                                lib.DGEMM(ref ta, ref tb, ref M, ref N, ref K, ref alpha, _pA, ref lda, _pB, ref ldb, ref beta, pC, ref ldc);
                            }

                            // correctness: both from the same initial C
                            Buffer.MemoryCopy(_pC0, _pC1, len * sizeof(double), len * sizeof(double));
                            Buffer.MemoryCopy(_pC0, _pC2, len * sizeof(double), len * sizeof(double));
                            Run(m_seq_BLAS, _pC1);
                            Run(ob, _pC2);
                            double err = ErrorRatio(_pC1, _pC2, len, K);

                            var cs = new ProviderCalibrationResult.Case() {
                                Name = $"dgemm {c.ta}{c.tb} {M}x{N}x{K}",
                                ErrorRatio = err,
                            };
                            ret.Cases.Add(cs);
                            if(!(err <= 1.0)) {
                                correctnessFailure = $"{cs.Name}: deviation {err:0.0e0} x round-off bound";
                                break;
                            }
                            // beta = 0.5 repeatedly applied keeps the values bounded
                            (cs.MKL, cs.OpenBLAS) = TimePair(() => Run(m_seq_BLAS, _pC1), () => Run(ob, _pC2), 2.0 * M * N * K);
                        }

                        if(!budgetExceeded && correctnessFailure == null) {
                            foreach(var c in gemvCases) {
                                if(total.Elapsed.TotalSeconds > CalibrationBudget) {
                                    budgetExceeded = true;
                                    break;
                                }
                                int M = c.M, N = c.N, t = c.t, lda = M, inc = 1;
                                double alpha = 1.25, beta = 0.5;
                                int lenY = c.t == 'n' ? M : N, K = c.t == 'n' ? N : M;
                                void Run(UnsafeDBLAS lib, double* pY) {
                                    lib.DGEMV(ref t, ref M, ref N, ref alpha, _pA, ref lda, _pB, ref inc, ref beta, pY, ref inc);
                                }
                                Buffer.MemoryCopy(_pC0, _pC1, lenY * sizeof(double), lenY * sizeof(double));
                                Buffer.MemoryCopy(_pC0, _pC2, lenY * sizeof(double), lenY * sizeof(double));
                                Run(m_seq_BLAS, _pC1);
                                Run(ob, _pC2);
                                double err = ErrorRatio(_pC1, _pC2, lenY, K);
                                var cs = new ProviderCalibrationResult.Case() {
                                    Name = $"dgemv {c.t} {M}x{N}",
                                    ErrorRatio = err,
                                };
                                ret.Cases.Add(cs);
                                if(!(err <= 1.0)) {
                                    correctnessFailure = $"{cs.Name}: deviation {err:0.0e0} x round-off bound";
                                    break;
                                }
                                (cs.MKL, cs.OpenBLAS) = TimePair(() => Run(m_seq_BLAS, _pC1), () => Run(ob, _pC2), 2.0 * M * N);
                            }
                        }

                        // threaded: only relevant if ordinary BLAS runs threaded
                        int budget = OpenBLASThreadBudget;
                        if(!budgetExceeded && correctnessFailure == null && budget > 1 && ilPSP.Environment.OpenMPenabled) {
                            lock(s_ProviderLock) {
                                SetOpenBLASThreads(budget);
                            }
                            int M = Nthr, N = Nthr, K = Nthr, ta = 'n', tb = 'n';
                            double alpha = 1.0, beta = 0.5;
                            void Run(UnsafeDBLAS lib, double* pC) {
                                lib.DGEMM(ref ta, ref tb, ref M, ref N, ref K, ref alpha, _pA, ref M, _pB, ref K, ref beta, pC, ref M);
                            }
                            var cs = new ProviderCalibrationResult.Case() {
                                Name = $"dgemm nn {M}x{N}x{K} ({budget} thr)",
                                Threaded = true,
                                ErrorRatio = double.NaN
                            };
                            (cs.MKL, cs.OpenBLAS) = TimePair(() => Run(m_omp_BLAS, _pC1), () => Run(ob, _pC2), 2.0 * M * N * K);
                            ret.Cases.Add(cs);
                            ret.ThreadedRatio = cs.Ratio;
                        }
                    }
                }
            } finally {
                lock(s_ProviderLock) {
                    if(s_Provider == BLASProvider.OpenBLAS)
                        ApplyProviderAndMode();
                    else
                        SetOpenBLASThreads(OpenBLASThreadBudget); // as if loaded with OPENBLAS_NUM_THREADS, e.g. for PETSc
                }
            }

            var seq = ret.Cases.Where(c => !c.Threaded && c.MKL > 0 && c.OpenBLAS > 0).ToArray();
            ret.SeqGeoMeanRatio = seq.Length > 0 ? Math.Exp(seq.Average(c => Math.Log(c.Ratio))) : double.NaN;
            ret.SeqWorstRatio = seq.Length > 0 ? seq.Max(c => c.Ratio) : double.NaN;
            ret.WallTime = total.Elapsed.TotalSeconds;

            string stats = string.Format(CultureInfo.InvariantCulture,
                "OpenBLAS/MKL time ratio: sequential geo. mean {0:0.000}, worst {1:0.000}, threaded {2:0.000}; {3} cases in {4:0.000} s",
                ret.SeqGeoMeanRatio, ret.SeqWorstRatio, ret.ThreadedRatio, ret.Cases.Count, ret.WallTime);
            if(correctnessFailure != null) {
                ret.OpenBLASWins = false;
                ret.Summary = "MKL kept, OpenBLAS failed correctness check (" + correctnessFailure + ")";
            } else if(budgetExceeded) {
                ret.OpenBLASWins = false;
                ret.Summary = $"MKL kept, calibration exceeded budget of {CalibrationBudget} s (incomplete); " + stats;
            } else if(!(ret.SeqGeoMeanRatio <= 1.0 / CalibrationWinMargin)) {
                ret.OpenBLASWins = false;
                ret.Summary = $"MKL kept, OpenBLAS not faster by margin {CalibrationWinMargin}; " + stats;
            } else if(!(ret.SeqWorstRatio <= CalibrationWorstCaseLimit)) {
                ret.OpenBLASWins = false;
                ret.Summary = $"MKL kept, OpenBLAS slower than limit {CalibrationWorstCaseLimit} on some case; " + stats;
            } else if(!double.IsNaN(ret.ThreadedRatio) && !(ret.ThreadedRatio <= CalibrationThreadedLimit)) {
                ret.OpenBLASWins = false;
                ret.Summary = $"MKL kept, threaded OpenBLAS slower than limit {CalibrationThreadedLimit}; " + stats;
            } else {
                ret.OpenBLASWins = true;
                ret.Summary = "OpenBLAS selected; " + stats;
            }
            return ret;
        }

        /// <summary>
        /// max. deviation between <paramref name="X"/> and <paramref name="Y"/>, relative to a round-off bound for a reduction of length <paramref name="K"/>:
        /// inputs are in [-0.5, 0.5] and alpha = 1.25, so the sum of absolute products is bounded by 0.3125*K, |beta*C| by 0.5;
        /// both results may deviate from the exact one by (K + 2)*eps times that.
        /// </summary>
        unsafe static double ErrorRatio(double* X, double* Y, int L, int K) {
            double maxDiff = 0;
            for(int i = 0; i < L; i++) {
                double d = Math.Abs(X[i] - Y[i]);
                if(double.IsNaN(d))
                    return double.PositiveInfinity;
                maxDiff = Math.Max(maxDiff, d);
            }
            double bound = 2.0 * (K + 2) * MachineEps * (0.3125 * K + 0.5);
            return maxDiff / bound;
        }

        /// <summary>
        /// best time per call (seconds) of each action; batches of both are alternated to cancel drift (e.g. CPU frequency)
        /// </summary>
        static (double tA, double tB) TimePair(Action a, Action b, double flops) {
            // warm-up
            a(); b(); a(); b();

            // batch size: about 0.2 ms of work, assuming 5 GFlop/s
            int reps = (int)Math.Max(1, Math.Min(1000, 0.2e-3 * 5e9 / flops));
            double bestA = double.MaxValue, bestB = double.MaxValue;
            var sw = new Stopwatch();
            for(int r = 0; r < 5; r++) {
                sw.Restart();
                for(int i = 0; i < reps; i++) a();
                bestA = Math.Min(bestA, sw.Elapsed.TotalSeconds / reps);

                sw.Restart();
                for(int i = 0; i < reps; i++) b();
                bestB = Math.Min(bestB, sw.Elapsed.TotalSeconds / reps);
            }
            return (bestA, bestB);
        }
    }
}

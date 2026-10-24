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
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading;
using MPI.Wrappers.Utils;

namespace ilPSP.Utils {

    // Optional native BLAS extensions of the BoSSS native libraries (PLAN-5 H4).
    // All of them are single-threaded; each capability is probed independently on the module handle which owns it,
    // and its absence never affects the mandatory BLAS binding (BLAS.m_seq_BLAS).
    // Modules are never unloaded, i.e. function pointers, JIT handles and kernels remain valid for the process lifetime.
    //
    // Prototypes, verified by disassembly of the installed libraries (2026-10-22; see FINDINGS.md), all C calling convention:
    //  - MKL (libBoSSSnative_seq.so), argument-forwarding wrappers with the CBLAS signatures, MKL_INT = int32:
    //    BoSSS_dgemm_batch_strided(layout, transa, transb, m, n, k, double alpha, A, lda, stridea, B, ldb, strideb, double beta, C, ldc, stridec, batch)
    //    BoSSS_dgemm_batch(layout, transa[], transb[], m[], n[], k[], alpha[], A[], lda[], B[], ldb[], beta[], C[], ldc[], group_count, group_size[])
    //    BoSSS_dgemv_batch_strided(layout, trans, m, n, double alpha, A, lda, stridea, x, incx, stridex, double beta, y, incy, stridey, batch)
    //    int BoSSS_mkl_jit_create_dgemm(void** jitter, layout, transa, transb, m, n, k, double alpha, lda, ldb, double beta, ldc)
    //    void* BoSSS_mkl_jit_get_dgemm_ptr(void* jitter) -> kernel void(void* jitter, double* a, double* b, double* c)
    //    int BoSSS_mkl_jit_destroy(void* jitter)
    //  - BLIS (embedded in libBoSSSnative_seq.so), typed API with int64 dimensions/strides, transpose flag 0/nonzero:
    //    BoSSS_bli_dgemm(int ta, int tb, int64 m, n, k, double alpha, a, int64 rsa, csa, b, int64 rsb, csb, double beta, c, int64 rsc, csc)
    //    BoSSS_bli_dgemv(int ta, int64 m, n, double alpha, a, int64 rsa, csa, x, int64 incx, double beta, y, int64 incy)
    //    (BLIS uses pthread mutexes for its memory pools, i.e. concurrent calls are safe.)
    //  - LIBXSMM 1.17 (libBoSSSnative_libxsmm.so), forwarding wrapper:
    //    kernel* BoSSS_libxsmm_dmmdispatch(int m, n, k, const int* lda, ldb, ldc, const double* alpha, beta, const int* flags, const int* prefetch)
    //    kernel: void(const double* a, const double* b, double* c), column-major, for prefetch = LIBXSMM_PREFETCH_NONE

    /// <summary>
    /// Intel MKL batched BLAS (sequential library); CBLAS conventions
    /// </summary>
    public static unsafe class MKLBatch {

        /// <summary> CBLAS layout </summary>
        public const int RowMajor = 101, ColMajor = 102;

        /// <summary> CBLAS transpose </summary>
        public const int NoTrans = 111, Trans = 112;

        static delegate* unmanaged<int, int, int, int, int, int, double, double*, int, int, double*, int, int, double, double*, int, int, int, void> s_dgemm_batch_strided;
        static delegate* unmanaged<int, int*, int*, int*, int*, int*, double*, double**, int*, double**, int*, double*, double**, int*, int, int*, void> s_dgemm_batch;
        static delegate* unmanaged<int, int, int, int, double, double*, int, int, double*, int, int, double, double*, int, int, int, void> s_dgemv_batch_strided;

        static readonly Lazy<string> s_Discovery = new Lazy<string>(delegate () {
            var h = BLAS.m_seq_BLAS.LibHandle;
            var p1 = NativeExtensionUtil.Sym(h, "BoSSS_dgemm_batch_strided");
            var p2 = NativeExtensionUtil.Sym(h, "BoSSS_dgemm_batch");
            var p3 = NativeExtensionUtil.Sym(h, "BoSSS_dgemv_batch_strided");
            string missing = NativeExtensionUtil.Missing(("BoSSS_dgemm_batch_strided", p1), ("BoSSS_dgemm_batch", p2), ("BoSSS_dgemv_batch_strided", p3));
            if(missing != null)
                return $"{missing} not exported by {BLAS.MKLLibraryPath}";
            s_dgemm_batch_strided = (delegate* unmanaged<int, int, int, int, int, int, double, double*, int, int, double*, int, int, double, double*, int, int, int, void>)p1;
            s_dgemm_batch = (delegate* unmanaged<int, int*, int*, int*, int*, int*, double*, double**, int*, double**, int*, double*, double**, int*, int, int*, void>)p2;
            s_dgemv_batch_strided = (delegate* unmanaged<int, int, int, int, double, double*, int, int, double*, int, int, double, double*, int, int, int, void>)p3;
            return null;
        }, LazyThreadSafetyMode.ExecutionAndPublication);

        /// <summary> true, if all batch routines are exported </summary>
        public static bool Available => s_Discovery.Value == null;

        /// <summary> null if <see cref="Available"/> </summary>
        public static string UnavailableReason => s_Discovery.Value;

        static void CheckAvail() {
            if(!Available)
                throw new NotSupportedException("MKL batch BLAS unavailable: " + UnavailableReason);
        }

        /// <summary> `cblas_dgemm_batch_strided` </summary>
        public static void DgemmBatchStrided(int layout, int transa, int transb, int m, int n, int k,
            double alpha, double* A, int lda, int strideA, double* B, int ldb, int strideB,
            double beta, double* C, int ldc, int strideC, int batch) {
            CheckAvail();
            s_dgemm_batch_strided(layout, transa, transb, m, n, k, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC, batch);
        }

        /// <summary> `cblas_dgemm_batch` (group API with pointer arrays) </summary>
        public static void DgemmBatch(int layout, int* transa, int* transb, int* m, int* n, int* k,
            double* alpha, double** A, int* lda, double** B, int* ldb, double* beta, double** C, int* ldc, int groupCount, int* groupSize) {
            CheckAvail();
            s_dgemm_batch(layout, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, groupCount, groupSize);
        }

        /// <summary> `cblas_dgemv_batch_strided` </summary>
        public static void DgemvBatchStrided(int layout, int trans, int m, int n, double alpha, double* A, int lda, int strideA,
            double* x, int incx, int strideX, double beta, double* y, int incy, int strideY, int batch) {
            CheckAvail();
            s_dgemv_batch_strided(layout, trans, m, n, alpha, A, lda, strideA, x, incx, strideX, beta, y, incy, strideY, batch);
        }
    }

    /// <summary>
    /// Intel MKL JIT DGEMM kernels (sequential library), cached per shape and scalars.
    /// </summary>
    public static unsafe class MKLJit {

        /// <summary> `MKL_JIT_SUCCESS`, `MKL_NO_JIT` (valid handle, but standard GEMM behind it), `MKL_JIT_ERROR` </summary>
        public const int JIT_SUCCESS = 0, NO_JIT = 1, JIT_ERROR = 2;

        /// <summary>
        /// maximum number of JIT kernels; beyond, no new kernels are created (existing ones remain valid)
        /// </summary>
        public const int MaxKernels = 4096;

        static delegate* unmanaged<IntPtr*, int, int, int, int, int, int, double, int, int, double, int, int> s_create;
        static delegate* unmanaged<IntPtr, IntPtr> s_get_ptr;

        static readonly Lazy<string> s_Discovery = new Lazy<string>(delegate () {
            var h = BLAS.m_seq_BLAS.LibHandle;
            var p1 = NativeExtensionUtil.Sym(h, "BoSSS_mkl_jit_create_dgemm");
            var p2 = NativeExtensionUtil.Sym(h, "BoSSS_mkl_jit_get_dgemm_ptr");
            var p3 = NativeExtensionUtil.Sym(h, "BoSSS_mkl_jit_destroy");
            string missing = NativeExtensionUtil.Missing(("BoSSS_mkl_jit_create_dgemm", p1), ("BoSSS_mkl_jit_get_dgemm_ptr", p2), ("BoSSS_mkl_jit_destroy", p3));
            if(missing != null)
                return $"{missing} not exported by {BLAS.MKLLibraryPath}";
            s_create = (delegate* unmanaged<IntPtr*, int, int, int, int, int, int, double, int, int, double, int, int>)p1;
            s_get_ptr = (delegate* unmanaged<IntPtr, IntPtr>)p2;
            return null;
        }, LazyThreadSafetyMode.ExecutionAndPublication);

        /// <summary> true, if the JIT routines are exported </summary>
        public static bool Available => s_Discovery.Value == null;

        /// <summary> null if <see cref="Available"/> </summary>
        public static string UnavailableReason => s_Discovery.Value;

        /// <summary>
        /// A JIT-compiled DGEMM; computes C = alpha*op(A)*op(B) + beta*C for the shape it was created for.
        /// </summary>
        public readonly struct Kernel {
            readonly IntPtr m_Jitter;
            readonly delegate* unmanaged<IntPtr, double*, double*, double*, void> m_Fn;

            internal Kernel(IntPtr jitter, IntPtr fn) {
                m_Jitter = jitter;
                m_Fn = (delegate* unmanaged<IntPtr, double*, double*, double*, void>)fn;
            }

            /// <summary> true for a usable kernel </summary>
            public bool IsValid => m_Fn != null;

            /// <summary> executes the kernel </summary>
            public void Invoke(double* A, double* B, double* C) {
                m_Fn(m_Jitter, A, B, C);
            }
        }

        readonly record struct Key(int Layout, int TransA, int TransB, int M, int N, int K, long AlphaBits, int Lda, int Ldb, long BetaBits, int Ldc);

        // value: kernel, or default (invalid) for shapes where MKL declined JIT (cached, so the decision is not repeated)
        static readonly ConcurrentDictionary<Key, Kernel> s_Cache = new ConcurrentDictionary<Key, Kernel>();
        static readonly object s_CreateLock = new object();
        static int s_NumKernels = 0;

        /// <summary> number of cached JIT kernels </summary>
        public static int NumKernels => s_NumKernels;

        /// <summary>
        /// Returns a cached JIT kernel or creates one; false if JIT is unavailable, declined by MKL (`MKL_NO_JIT`, error) or the cache is full.
        /// The exact values of <paramref name="alpha"/> and <paramref name="beta"/> are part of the key, since they are compiled into the kernel.
        /// </summary>
        public static bool TryGetKernel(int layout, int transa, int transb, int m, int n, int k, double alpha, int lda, int ldb, double beta, int ldc, out Kernel kernel) {
            var key = new Key(layout, transa, transb, m, n, k, BitConverter.DoubleToInt64Bits(alpha), lda, ldb, BitConverter.DoubleToInt64Bits(beta), ldc);
            if(s_Cache.TryGetValue(key, out kernel))
                return kernel.IsValid;
            if(!Available) {
                kernel = default;
                return false;
            }
            lock(s_CreateLock) {
                if(s_Cache.TryGetValue(key, out kernel))
                    return kernel.IsValid;
                if(s_Cache.Count >= MaxKernels) {
                    kernel = default; // full: skip new specializations, no eviction
                    return false;
                }
                IntPtr jitter = IntPtr.Zero;
                int status = s_create(&jitter, layout, transa, transb, m, n, k, alpha, lda, ldb, beta, ldc);
                if(status == JIT_SUCCESS && jitter != IntPtr.Zero) {
                    IntPtr fn = s_get_ptr(jitter);
                    kernel = fn != IntPtr.Zero ? new Kernel(jitter, fn) : default;
                    if(kernel.IsValid)
                        s_NumKernels++;
                } else {
                    // MKL_NO_JIT: a valid handle with the standard GEMM behind it, no benefit over the ordinary call;
                    // the handle is kept (never destroyed), but not used.
                    kernel = default;
                }
                s_Cache[key] = kernel;
                return kernel.IsValid;
            }
        }
    }

    /// <summary>
    /// BLIS (embedded in the sequential MKL library), single-threaded, arbitrary row and column strides
    /// </summary>
    public static unsafe class BLIS {

        static delegate* unmanaged<int, int, long, long, long, double, double*, long, long, double*, long, long, double, double*, long, long, void> s_dgemm;
        static delegate* unmanaged<int, long, long, double, double*, long, long, double*, long, double, double*, long, void> s_dgemv;

        static readonly Lazy<string> s_Discovery = new Lazy<string>(delegate () {
            var h = BLAS.m_seq_BLAS.LibHandle;
            var pi = NativeExtensionUtil.Sym(h, "BoSSS_bli_init");
            var p1 = NativeExtensionUtil.Sym(h, "BoSSS_bli_dgemm");
            var p2 = NativeExtensionUtil.Sym(h, "BoSSS_bli_dgemv");
            string missing = NativeExtensionUtil.Missing(("BoSSS_bli_init", pi), ("BoSSS_bli_dgemm", p1), ("BoSSS_bli_dgemv", p2));
            if(missing != null)
                return $"{missing} not exported by {BLAS.MKLLibraryPath}";
            ((delegate* unmanaged<void>)pi)(); // once; never finalized, since other users might hold BLIS objects
            s_dgemm = (delegate* unmanaged<int, int, long, long, long, double, double*, long, long, double*, long, long, double, double*, long, long, void>)p1;
            s_dgemv = (delegate* unmanaged<int, long, long, double, double*, long, long, double*, long, double, double*, long, void>)p2;
            return null;
        }, LazyThreadSafetyMode.ExecutionAndPublication);

        /// <summary> true, if BLIS is exported </summary>
        public static bool Available => s_Discovery.Value == null;

        /// <summary> null if <see cref="Available"/> </summary>
        public static string UnavailableReason => s_Discovery.Value;

        /// <summary>
        /// C = alpha*op(A)*op(B) + beta*C, where X(i,j) is located at X + i*rsX + j*csX;
        /// note: for beta = 0, BLIS overwrites C without reading it (NaN/Inf in C are not propagated).
        /// </summary>
        public static void Dgemm(bool transA, bool transB, long m, long n, long k, double alpha,
            double* A, long rsA, long csA, double* B, long rsB, long csB, double beta, double* C, long rsC, long csC) {
            if(!Available)
                throw new NotSupportedException("BLIS unavailable: " + UnavailableReason);
            s_dgemm(transA ? 1 : 0, transB ? 1 : 0, m, n, k, alpha, A, rsA, csA, B, rsB, csB, beta, C, rsC, csC);
        }

        /// <summary>
        /// y = alpha*op(A)*x + beta*y, A(i,j) at A + i*rsA + j*csA, A is m-by-n
        /// </summary>
        public static void Dgemv(bool transA, long m, long n, double alpha, double* A, long rsA, long csA,
            double* x, long incx, double beta, double* y, long incy) {
            if(!Available)
                throw new NotSupportedException("BLIS unavailable: " + UnavailableReason);
            s_dgemv(transA ? 1 : 0, m, n, alpha, A, rsA, csA, x, incx, beta, y, incy);
        }
    }

    /// <summary>
    /// LIBXSMM 1.17 small-matrix kernels (separate library `libBoSSSnative_libxsmm.so`), single-threaded, column-major,
    /// no transposition; only alpha = 1 and beta in {0, 1}.
    /// </summary>
    public static unsafe class LIBXSMM {

        /// <summary> library name </summary>
        public const string LinuxLibName = "libBoSSSnative_libxsmm.so";

        /// <summary> maximum number of cached dispatch results </summary>
        public const int MaxKernels = 4096;

        static delegate* unmanaged<int, int, int, int*, int*, int*, double*, double*, int*, int*, IntPtr> s_dispatch;

        static string s_LibPath;

        static readonly Lazy<string> s_Discovery = new Lazy<string>(delegate () {
            if(System.Environment.OSVersion.Platform != PlatformID.Unix)
                return "LIBXSMM is only supported on Linux";
            string path = NativeExtensionUtil.FindLibrary(LinuxLibName, out string notFound);
            if(path == null)
                return notFound;
            var h = DynamicLibraries.LoadDynLib(path, out string err, false); // never unloaded
            if(h.val == IntPtr.Zero)
                return $"{path} cannot be loaded: {err}";
            var pi = NativeExtensionUtil.Sym(h, "BoSSS_libxsmm_init");
            var pd = NativeExtensionUtil.Sym(h, "BoSSS_libxsmm_dmmdispatch");
            string missing = NativeExtensionUtil.Missing(("BoSSS_libxsmm_init", pi), ("BoSSS_libxsmm_dmmdispatch", pd));
            if(missing != null)
                return $"{missing} not exported by {path}";
            ((delegate* unmanaged<void>)pi)(); // once; never finalized while kernels may be alive
            s_dispatch = (delegate* unmanaged<int, int, int, int*, int*, int*, double*, double*, int*, int*, IntPtr>)pd;
            s_LibPath = path;
            return null;
        }, LazyThreadSafetyMode.ExecutionAndPublication);

        /// <summary> true, if LIBXSMM can be loaded </summary>
        public static bool Available => s_Discovery.Value == null;

        /// <summary> null if <see cref="Available"/> </summary>
        public static string UnavailableReason => s_Discovery.Value;

        /// <summary> full path of the loaded library </summary>
        public static string LibraryPath => Available ? s_LibPath : null;

        /// <summary>
        /// A LIBXSMM kernel; C = A*B + beta*C (column-major, leading dimensions as dispatched)
        /// </summary>
        public readonly struct Kernel {
            readonly delegate* unmanaged<double*, double*, double*, void> m_Fn;

            internal Kernel(IntPtr fn) {
                m_Fn = (delegate* unmanaged<double*, double*, double*, void>)fn;
            }

            /// <summary> true for a usable kernel </summary>
            public bool IsValid => m_Fn != null;

            /// <summary> executes the kernel </summary>
            public void Invoke(double* A, double* B, double* C) {
                m_Fn(A, B, C);
            }
        }

        readonly record struct Key(int M, int N, int K, int Lda, int Ldb, int Ldc, bool Beta0);

        static readonly ConcurrentDictionary<Key, Kernel> s_Cache = new ConcurrentDictionary<Key, Kernel>();

        /// <summary>
        /// Kernel for C = A*B + beta*C (column-major, no transposition, alpha = 1, beta = 0 or 1);
        /// false if LIBXSMM is unavailable, the dispatch declined the shape (null) or the cache is full.
        /// </summary>
        public static bool TryGetKernel(int m, int n, int k, int lda, int ldb, int ldc, bool beta0, out Kernel kernel) {
            var key = new Key(m, n, k, lda, ldb, ldc, beta0);
            if(s_Cache.TryGetValue(key, out kernel))
                return kernel.IsValid;
            if(!Available || s_Cache.Count >= MaxKernels) {
                kernel = default;
                return false;
            }
            // dispatch is thread-safe and returns the same kernel for the same descriptor
            double alpha = 1.0, beta = beta0 ? 0.0 : 1.0;
            int flags = 0;    // LIBXSMM_GEMM_FLAG_NONE: no transposition
            int prefetch = 0; // LIBXSMM_PREFETCH_NONE: kernel takes exactly (a, b, c)
            IntPtr fn = s_dispatch(m, n, k, &lda, &ldb, &ldc, &alpha, &beta, &flags, &prefetch);
            kernel = fn != IntPtr.Zero ? new Kernel(fn) : default;
            s_Cache.TryAdd(key, kernel);
            return kernel.IsValid;
        }
    }

    static class NativeExtensionUtil {

        internal static IntPtr Sym(DynamicLibraries.DynLibHandle h, string name) {
            if(h.val == IntPtr.Zero)
                return IntPtr.Zero;
            return DynamicLibraries.LoadSymbol(h, name, out _);
        }

        internal static string Missing(params (string name, IntPtr p)[] syms) {
            var m = syms.Where(s => s.p == IntPtr.Zero).Select(s => s.name).ToArray();
            return m.Length == 0 ? null : "missing export(s) " + string.Join(", ", m);
        }

        /// <summary>
        /// full path of <paramref name="libName"/> in the native library search path (<see cref="DynLibLoader.PrimaryLibrarySearchPath"/>, `LD_LIBRARY_PATH`)
        /// </summary>
        internal static string FindLibrary(string libName, out string reason) {
            var dirs = new List<string>();
            if(!string.IsNullOrWhiteSpace(DynLibLoader.PrimaryLibrarySearchPath))
                dirs.Add(DynLibLoader.PrimaryLibrarySearchPath);
            foreach(var d in (System.Environment.GetEnvironmentVariable("LD_LIBRARY_PATH") ?? "").Split(':', StringSplitOptions.RemoveEmptyEntries))
                if(!dirs.Contains(d))
                    dirs.Add(d);
            string path = dirs.Select(d => Path.Combine(d, libName)).FirstOrDefault(File.Exists);
            reason = path == null ? $"{libName} not found in {string.Join(":", dirs)}" : null;
            return path;
        }
    }
}

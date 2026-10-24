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
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using ilPSP.Tracing;
using MPI.Wrappers.Utils;

namespace ilPSP.Utils {
    /// <summary>Thread controls bound to the same module as the OpenBLAS imports.</summary>
    unsafe sealed class OpenBLASControl : DynLibLoader {
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        delegate void SetThreads(int count);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        delegate int GetInt();
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        delegate void GetConfiguration(byte* buffer, int capacity);
#pragma warning disable 649
        SetThreads openblas_set_num_threads;
        GetInt openblas_get_num_threads;
        GetInt openblas_get_parallel;
        GetConfiguration openblas_get_config;
#pragma warning restore 649

        internal OpenBLASControl(DynLibLoader owner) : base(owner, BoSSS_Prefix) { }
        internal void SetNumThreads(int n) => openblas_set_num_threads(n);
        internal int GetNumThreads() => openblas_get_num_threads();
        internal int GetParallel() => openblas_get_parallel();
        internal string GetConfig() {
            byte* buffer = stackalloc byte[512];
            buffer[0] = 0;
            openblas_get_config(buffer, 512);
            return Marshal.PtrToStringAnsi((IntPtr)buffer);
        }
    }

    public static partial class BLAS {
        // ==========================================================================================
        // OpenBLAS discovery
        // ==========================================================================================

        sealed class OpenBLASModule {
            internal UnsafeDBLAS Blas;
            internal UnsafeLAPACK Lapack;
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
                        string timeout = DynLibLoader.SetNativeEnvironmentDefault("OPENBLAS_THREAD_TIMEOUT", OpenBLASThreadTimeoutDefault);
                        // avoids that OpenBLAS creates a thread for every CPU in the mask at load time; the count can be raised later
                        string numThreads = DynLibLoader.SetNativeEnvironmentDefault("OPENBLAS_NUM_THREADS", OpenBLASThreadBudget.ToString(CultureInfo.InvariantCulture));
                        info.Append($"OPENBLAS_THREAD_TIMEOUT={timeout}, OPENBLAS_NUM_THREADS={numThreads}; ");
                    }

                    var blas = DynLibLoader.TryLoad(() => new UnsafeDBLAS(Parallelism.SEQ, BLASProvider.OpenBLAS), out string reason);
                    if(blas == null) {
                        s_OpenBLASUnavailableReason = reason;
                        return null;
                    }
                    var lapack = DynLibLoader.TryLoad(() => new UnsafeLAPACK(blas), out reason);
                    var ctrl = DynLibLoader.TryLoad(() => new OpenBLASControl(blas), out string controlReason);
                    if(lapack == null || ctrl == null) {
                        // BLAS and LAPACK must both be usable before publishing the provider.
                        // Keep the native module loaded: PETSc may already hold references to it.
                        s_OpenBLASUnavailableReason = reason ?? controlReason;
                        return null;
                    }

                    var r = new OpenBLASModule() {
                        Blas = blas,
                        Lapack = lapack,
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

    }
}

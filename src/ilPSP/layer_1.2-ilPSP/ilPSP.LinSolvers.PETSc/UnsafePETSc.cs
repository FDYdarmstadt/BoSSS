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
using System.Runtime.InteropServices;
using System.Security;
using System.Text;
using MPI.Wrappers.Utils;

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// Dynamic binding to the BoSSS-owned PETSc C ABI in a dedicated native library.
    /// </summary>
    [SuppressUnmanagedCodeSecurity]
    public sealed class UnsafePETSc : DynLibLoader {

        static (string[] _LibNames, string[][][] PrequesiteLibraries, GetNameMangling[] NameMangling, PlatformID[] OsFilter, int[] PointerSizeFilter) SelectLibrary(PETScNativeVariant variant) {
            string linuxName;
            string windowsName;
            bool includeLegacyLinuxName = false;
            switch(variant) {
                case PETScNativeVariant.Mpi:
                    linuxName = "libBoSSSnative_petsc_mpi.so";
                    windowsName = "BoSSSnative_petsc_mpi.dll";
                    includeLegacyLinuxName = true;
                    break;
                case PETScNativeVariant.OpenMp:
                    linuxName = "libBoSSSnative_petsc_omp.so";
                    windowsName = "BoSSSnative_petsc_omp.dll";
                    break;
                case PETScNativeVariant.MpiOpenMp:
                    linuxName = "libBoSSSnative_petsc_mpi_omp.so";
                    windowsName = "BoSSSnative_petsc_mpi_omp.dll";
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(variant), variant, "Unknown PETSc native variant.");
            }

            string[] libraryNames = includeLegacyLinuxName
                ? new[] { windowsName, linuxName, "libBoSSSnative_mpi.so" }
                : new[] { windowsName, linuxName };
            var nameMangling = new GetNameMangling[libraryNames.Length];
            var osFilter = new PlatformID[libraryNames.Length];
            var pointerSizeFilter = new int[libraryNames.Length];
            for(int i = 0; i < libraryNames.Length; i++) {
                nameMangling[i] = DynLibLoader.Identity;
                osFilter[i] = i == 0 ? PlatformID.Win32NT : PlatformID.Unix;
                pointerSizeFilter[i] = -1;
            }
            return (libraryNames, new string[libraryNames.Length][][], nameMangling, osFilter, pointerSizeFilter);
        }

        /// <summary>
        /// Loads the selected native library and resolves the PETSc entry points.
        /// </summary>
        public UnsafePETSc(PETScNativeVariant variant = PETScNativeVariant.Mpi) : base(SelectLibrary(variant)) {
            NativeVariant = variant;
        }

        /// <summary>
        /// Native runtime selected by this loader.
        /// </summary>
        public PETScNativeVariant NativeVariant { get; }

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_initialize();

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _BoSSS_PETSc_Initialize();

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _BoSSS_PETSc_GetVersion(byte[] buffer, int bufferLength);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_create_matrix_solver(
            int localRows,
            int globalRows,
            int globalRowStart,
            int[] rowPointers,
            int[] columnIndices,
            double[] values,
            out IntPtr solver);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_create_split_csr_matrix_solver(
            int localRows,
            int globalRows,
            int globalRowStart,
            int[] diagonalRowPointers,
            int[] diagonalColumnIndices,
            double[] diagonalValues,
            int[] offdiagonalRowPointers,
            int[] offdiagonalColumnIndices,
            double[] offdiagonalValues,
            out IntPtr solver);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_apply(
            IntPtr context,
            IntPtr ownedInput,
            int ownedLength,
            IntPtr haloInput,
            int haloLength,
            IntPtr output);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_create_shell_solver(
            int localRows,
            int globalRows,
            int globalRowStart,
            int ghostCount,
            int[] ghostGlobalIds,
            [MarshalAs(UnmanagedType.LPStr)] string petscOptions,
            _bosss_petsc_apply applyCallback,
            IntPtr matrixMultContext,
            out IntPtr solver);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_create_ksp_for_matrix(IntPtr solver);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_configure_ksp_solver(
            IntPtr solver,
            [MarshalAs(UnmanagedType.LPStr)] string petscOptions,
            int firstFieldRowsCount,
            int[] firstFieldRowIndices,
            int secondFieldRowsCount,
            int[] secondFieldRowIndices);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_setup_ksp_solver(IntPtr solver);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate void _bosss_petsc_ksp_monitor(IntPtr context, int iteration, double residualNorm);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_set_ksp_monitor(
            IntPtr solver,
            _bosss_petsc_ksp_monitor monitor,
            IntPtr context);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_solve(
            IntPtr solver,
            double[] solution,
            double[] rightHandSide,
            int useInitialGuess,
            out int iterations,
            out double residualNorm,
            out int convergenceReason);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _bosss_petsc_destroy_solver(IntPtr solver);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate IntPtr _bosss_petsc_last_error();

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _BoSSS_PETSc_Finalize();

#pragma warning disable 649
        _bosss_petsc_initialize bosss_petsc_initialize;
        _BoSSS_PETSc_Initialize BoSSS_PETSc_Initialize;
        _BoSSS_PETSc_GetVersion BoSSS_PETSc_GetVersion;
        _bosss_petsc_create_matrix_solver bosss_petsc_create_matrix_solver;
        _bosss_petsc_create_split_csr_matrix_solver bosss_petsc_create_split_csr_matrix_solver;
        _bosss_petsc_create_shell_solver bosss_petsc_create_shell_solver;
        _bosss_petsc_create_ksp_for_matrix bosss_petsc_create_ksp_for_matrix;
        _bosss_petsc_configure_ksp_solver bosss_petsc_configure_ksp_solver;
        _bosss_petsc_setup_ksp_solver bosss_petsc_setup_ksp_solver;
        _bosss_petsc_set_ksp_monitor bosss_petsc_set_ksp_monitor;
        _bosss_petsc_solve bosss_petsc_solve;
        _bosss_petsc_destroy_solver bosss_petsc_destroy_solver;
        _bosss_petsc_last_error bosss_petsc_last_error;
        _BoSSS_PETSc_Finalize BoSSS_PETSc_Finalize;
#pragma warning restore 649

        internal _bosss_petsc_initialize Initialize => bosss_petsc_initialize;

        internal _BoSSS_PETSc_Initialize InitializeRuntime => BoSSS_PETSc_Initialize;

        /// <summary>
        /// Gets the PETSc release reported by the loaded native PETSc wrapper.
        /// </summary>
        public string GetPETScVersion() {
            var buffer = new byte[256];
            int error = BoSSS_PETSc_GetVersion(buffer, buffer.Length);
            if(error != 0) {
                throw new InvalidOperationException($"PETSc version query failed with error {error}.");
            }
            int length = Array.IndexOf(buffer, (byte)0);
            if(length < 0) {
                length = buffer.Length;
            }
            return Encoding.UTF8.GetString(buffer, 0, length);
        }

        /// <summary>
        /// Returns a concise description of the loaded PETSc release and native capabilities.
        /// </summary>
        public string GetPETScRuntimeDescription() {
            string capabilities;
            string mpiRuntime = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                ? "Microsoft MPI 10"
                : "OpenMPI 4";
            switch(NativeVariant) {
                case PETScNativeVariant.Mpi:
                    capabilities = mpiRuntime;
                    break;
                case PETScNativeVariant.OpenMp:
                    capabilities = "MPIUNI, OpenMP";
                    break;
                case PETScNativeVariant.MpiOpenMp:
                    capabilities = mpiRuntime + ", OpenMP";
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(NativeVariant), NativeVariant, "Unknown PETSc native variant.");
            }
            return $"Loaded PETSc Version {GetPETScVersion()} with capabilities {capabilities}.";
        }

        internal _bosss_petsc_create_matrix_solver CreateMatrixSolver => bosss_petsc_create_matrix_solver;

        internal _bosss_petsc_create_split_csr_matrix_solver CreateSplitCsrMatrixSolver => bosss_petsc_create_split_csr_matrix_solver;

        internal _bosss_petsc_create_shell_solver CreateShellSolver => bosss_petsc_create_shell_solver;

        internal _bosss_petsc_create_ksp_for_matrix CreateKspForMatrix => bosss_petsc_create_ksp_for_matrix;

        internal _bosss_petsc_configure_ksp_solver ConfigureKspSolver => bosss_petsc_configure_ksp_solver;

        internal _bosss_petsc_setup_ksp_solver SetupKspSolver => bosss_petsc_setup_ksp_solver;

        internal _bosss_petsc_set_ksp_monitor SetKspMonitor => bosss_petsc_set_ksp_monitor;

        internal _bosss_petsc_solve Solve => bosss_petsc_solve;

        internal _bosss_petsc_destroy_solver DestroySolver => bosss_petsc_destroy_solver;

        internal _bosss_petsc_last_error LastError => bosss_petsc_last_error;

        /// <summary>
        /// Finalizes PETSc at application shutdown. All PETSc solvers must already be disposed.
        /// </summary>
        public int FinalizePETScRuntime() {
            return BoSSS_PETSc_Finalize();
        }
    }
}

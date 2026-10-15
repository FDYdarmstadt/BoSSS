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
using System.Threading;
using MPI.Wrappers.Utils;

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// Direct dynamic binding to the PETSc C ABI used by BoSSS PETSc solvers.
    /// </summary>
    /// <remarks>
    /// This class deliberately exposes only the PETSc functions needed by the
    /// existing BoSSS wrapper and direct-ABI solvers. It assumes the BoSSS PETSc ABI contract:
    /// real double scalars and 32-bit <c>PetscInt</c>. PETSc handles are opaque
    /// pointers; ownership follows the PETSc C API exactly.
    /// </remarks>
    [SuppressUnmanagedCodeSecurity]
    public sealed class UnsafePETSc : DynLibLoader {

        static (string[] _LibNames, string[][][] PrequesiteLibraries, GetNameMangling[] NameMangling, PlatformID[] OsFilter, int[] PointerSizeFilter) SelectLibrary(PETScNativeVariant variant) {
            string windowsName;
            string linuxName;
            switch(variant) {
                case PETScNativeVariant.Mpi:
                    windowsName = "libpetsc_mpi.dll";
                    linuxName = "libpetsc_mpi.so";
                    break;
                case PETScNativeVariant.OpenMp:
                    windowsName = "libpetsc_omp.dll";
                    linuxName = "libpetsc_omp.so";
                    break;
                case PETScNativeVariant.MpiOpenMp:
                    windowsName = "libpetsc_mpi_omp.dll";
                    linuxName = "libpetsc_mpi_omp.so";
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(variant), variant, "Unknown PETSc native variant.");
            }

            return (
                new[] { windowsName, linuxName },
                new[] { Array.Empty<string[]>(), Array.Empty<string[]>() },
                new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.Identity },
                new[] { PlatformID.Win32NT, PlatformID.Unix },
                new[] { 8, 8 });
        }

        /// <summary>
        /// Loads the shared PETSc ABI corresponding to <paramref name="variant"/>.
        /// </summary>
        /// <remarks>
        /// Windows and Linux shared PETSc artifacts are produced by the PETSc
        /// platform build scripts and staged beside BoSSS native libraries.
        /// </remarks>
        public UnsafePETSc(PETScNativeVariant variant = PETScNativeVariant.Mpi) : base(SelectLibrary(variant)) {
            NativeVariant = variant;
        }

        /// <summary>Native PETSc runtime selected by this binding.</summary>
        public PETScNativeVariant NativeVariant { get; }

        /// <summary>Current BoSSS PETSc ABI uses a 32-bit <c>PetscInt</c>.</summary>
        public const int PetscIntSize = 4;

        /// <summary>Current BoSSS PETSc ABI uses a real double <c>PetscScalar</c>.</summary>
        public const int PetscScalarSize = 8;

        /// <summary>Signature of a PETSc error handler.</summary>
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        public delegate int PetscErrorHandler(IntPtr communicator, int line, IntPtr function, IntPtr file, int errorCode, int errorType, IntPtr message, IntPtr context);

        /// <summary>Signature of a PETSc KSP convergence monitor.</summary>
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        public delegate int KspMonitor(IntPtr ksp, int iteration, double residualNorm, IntPtr context);

        /// <summary>Signature of a PETSc shell-preconditioner apply callback.</summary>
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        public delegate int PcShellApply(IntPtr pc, IntPtr input, IntPtr output);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscInitializeNoArguments();
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscInitialized(out int initialized);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscFinalize();
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscGetVersion([Out] byte[] version, UIntPtr length);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscErrorMessage(int errorCode, out IntPtr text, IntPtr specific);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscPushErrorHandler(PetscErrorHandler handler, IntPtr context);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscPopErrorHandler();
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscReturnErrorHandler(IntPtr communicator, int line, IntPtr function, IntPtr file, int errorCode, int errorType, IntPtr message, IntPtr context);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscOptionsCreate(out IntPtr options);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscOptionsDestroy(ref IntPtr options);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscOptionsInsertString(IntPtr options, [MarshalAs(UnmanagedType.LPStr)] string optionsString);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscOptionsClear(IntPtr options);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscObjectSetOptions(IntPtr petscObject, IntPtr options);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PetscObjectTypeCompare(IntPtr petscObject, [MarshalAs(UnmanagedType.LPStr)] string typeName, out int matches);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _ISCreateGeneral(IntPtr communicator, int length, IntPtr indices, int copyMode, out IntPtr indexSet);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _ISCreateStride(IntPtr communicator, int length, int first, int step, out IntPtr indexSet);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _ISDestroy(ref IntPtr indexSet);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecCreateMPI(IntPtr communicator, int localSize, int globalSize, out IntPtr vector);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecCreateSeq(IntPtr communicator, int size, out IntPtr vector);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecDestroy(ref IntPtr vector);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecSet(IntPtr vector, double value);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecGetArray(IntPtr vector, out IntPtr values);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecGetArrayRead(IntPtr vector, out IntPtr values);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecRestoreArray(IntPtr vector, ref IntPtr values);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecRestoreArrayRead(IntPtr vector, ref IntPtr values);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecScatterCreate(IntPtr source, IntPtr sourceIndices, IntPtr destination, IntPtr destinationIndices, out IntPtr scatter);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecScatterBegin(IntPtr scatter, IntPtr source, IntPtr destination, int insertMode, int scatterMode);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecScatterEnd(IntPtr scatter, IntPtr source, IntPtr destination, int insertMode, int scatterMode);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _VecScatterDestroy(ref IntPtr scatter);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatCreateAIJ(IntPtr communicator, int localRows, int localColumns, int globalRows, int globalColumns, int diagonalNonzeros, IntPtr diagonalNonzerosPerRow, int offDiagonalNonzeros, IntPtr offDiagonalNonzerosPerRow, out IntPtr matrix);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatSeqAIJSetPreallocationCSR(IntPtr matrix, IntPtr rowPointers, IntPtr columnIndices, IntPtr values);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatMPIAIJSetPreallocationCSR(IntPtr matrix, IntPtr rowPointers, IntPtr columnIndices, IntPtr values);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatCreateMPIAIJWithSplitArrays(IntPtr communicator, int localRows, int localColumns, int globalRows, int globalColumns, IntPtr diagonalRowPointers, IntPtr diagonalColumnIndices, IntPtr diagonalValues, IntPtr offDiagonalRowPointers, IntPtr offDiagonalColumnIndices, IntPtr offDiagonalValues, out IntPtr matrix);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatCreateShell(IntPtr communicator, int localRows, int localColumns, int globalRows, int globalColumns, IntPtr context, out IntPtr matrix);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatDestroy(ref IntPtr matrix);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatSetValues(IntPtr matrix, int rowCount, IntPtr rows, int columnCount, IntPtr columns, IntPtr values, int insertMode);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatAssemblyBegin(IntPtr matrix, int assemblyType);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatAssemblyEnd(IntPtr matrix, int assemblyType);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatGetOwnershipRange(IntPtr matrix, out int firstRow, out int lastRow);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatCreateVecs(IntPtr matrix, out IntPtr rightVector, out IntPtr leftVector);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatShellGetContext(IntPtr matrix, out IntPtr context);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _MatShellSetOperation(IntPtr matrix, int operation, IntPtr callback);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPCreate(IntPtr communicator, out IntPtr ksp);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPDestroy(ref IntPtr ksp);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPSetType(IntPtr ksp, [MarshalAs(UnmanagedType.LPStr)] string typeName);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPSetUp(IntPtr ksp);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPSolve(IntPtr ksp, IntPtr rightHandSide, IntPtr solution);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPSetTolerances(IntPtr ksp, double relativeTolerance, double absoluteTolerance, double divergenceTolerance, int maximumIterations);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPSetInitialGuessNonzero(IntPtr ksp, int initialGuessNonzero);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPGetResidualNorm(IntPtr ksp, out double residualNorm);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPGetIterationNumber(IntPtr ksp, out int iterations);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPGetConvergedReason(IntPtr ksp, out int reason);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPGetPC(IntPtr ksp, out IntPtr pc);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPMonitorSet(IntPtr ksp, KspMonitor monitor, IntPtr context, IntPtr destroy);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPGMRESSetRestart(IntPtr ksp, int restart);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPSetFromOptions(IntPtr ksp);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _KSPSetOperators(IntPtr ksp, IntPtr matrix, IntPtr preconditionerMatrix);

        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCSetType(IntPtr pc, [MarshalAs(UnmanagedType.LPStr)] string typeName);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCShellSetApply(IntPtr pc, PcShellApply apply);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCShellSetContext(IntPtr pc, IntPtr context);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCShellGetContext(IntPtr pc, out IntPtr context);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCShellSetName(IntPtr pc, [MarshalAs(UnmanagedType.LPStr)] string name);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCFactorSetMatSolverType(IntPtr pc, [MarshalAs(UnmanagedType.LPStr)] string solverType);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCFactorSetMatOrderingType(IntPtr pc, [MarshalAs(UnmanagedType.LPStr)] string orderingType);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCFactorSetLevels(IntPtr pc, int levels);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCASMSetOverlap(IntPtr pc, int overlap);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCASMSetType(IntPtr pc, int type);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCASMSetLocalType(IntPtr pc, int type);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCASMGetSubKSP(IntPtr pc, out int localSubdomains, out int firstSubdomain, out IntPtr subdomainKsps);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCFieldSplitSetIS(IntPtr pc, [MarshalAs(UnmanagedType.LPStr)] string splitName, IntPtr indexSet);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCMGSetLevels(IntPtr pc, int levels, IntPtr communicator);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCMGSetInterpolation(IntPtr pc, int level, IntPtr interpolation);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCMGSetRestriction(IntPtr pc, int level, IntPtr restriction);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCMGGetSmoother(IntPtr pc, int level, out IntPtr ksp);
        [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
        internal delegate int _PCMGSetGalerkin(IntPtr pc, int galerkinType);

#pragma warning disable 649
        _PetscInitializeNoArguments PetscInitializeNoArguments;
        _PetscInitialized PetscInitialized;
        _PetscFinalize PetscFinalize;
        _PetscGetVersion PetscGetVersion;
        _PetscErrorMessage PetscErrorMessage;
        _PetscPushErrorHandler PetscPushErrorHandler;
        _PetscPopErrorHandler PetscPopErrorHandler;
        _PetscReturnErrorHandler PetscReturnErrorHandler;
        _PetscOptionsCreate PetscOptionsCreate;
        _PetscOptionsDestroy PetscOptionsDestroy;
        _PetscOptionsInsertString PetscOptionsInsertString;
        _PetscOptionsClear PetscOptionsClear;
        _PetscObjectSetOptions PetscObjectSetOptions;
        _PetscObjectTypeCompare PetscObjectTypeCompare;
        _ISCreateGeneral ISCreateGeneral;
        _ISCreateStride ISCreateStride;
        _ISDestroy ISDestroy;
        _VecCreateMPI VecCreateMPI;
        _VecCreateSeq VecCreateSeq;
        _VecDestroy VecDestroy;
        _VecSet VecSet;
        _VecGetArray VecGetArray;
        _VecGetArrayRead VecGetArrayRead;
        _VecRestoreArray VecRestoreArray;
        _VecRestoreArrayRead VecRestoreArrayRead;
        _VecScatterCreate VecScatterCreate;
        _VecScatterBegin VecScatterBegin;
        _VecScatterEnd VecScatterEnd;
        _VecScatterDestroy VecScatterDestroy;
        _MatCreateAIJ MatCreateAIJ;
        _MatSeqAIJSetPreallocationCSR MatSeqAIJSetPreallocationCSR;
        _MatMPIAIJSetPreallocationCSR MatMPIAIJSetPreallocationCSR;
        _MatCreateMPIAIJWithSplitArrays MatCreateMPIAIJWithSplitArrays;
        _MatCreateShell MatCreateShell;
        _MatDestroy MatDestroy;
        _MatSetValues MatSetValues;
        _MatAssemblyBegin MatAssemblyBegin;
        _MatAssemblyEnd MatAssemblyEnd;
        _MatGetOwnershipRange MatGetOwnershipRange;
        _MatCreateVecs MatCreateVecs;
        _MatShellGetContext MatShellGetContext;
        _MatShellSetOperation MatShellSetOperation;
        _KSPCreate KSPCreate;
        _KSPDestroy KSPDestroy;
        _KSPSetType KSPSetType;
        _KSPSetUp KSPSetUp;
        _KSPSolve KSPSolve;
        _KSPSetTolerances KSPSetTolerances;
        _KSPSetInitialGuessNonzero KSPSetInitialGuessNonzero;
        _KSPGetResidualNorm KSPGetResidualNorm;
        _KSPGetIterationNumber KSPGetIterationNumber;
        _KSPGetConvergedReason KSPGetConvergedReason;
        _KSPGetPC KSPGetPC;
        _KSPMonitorSet KSPMonitorSet;
        _KSPGMRESSetRestart KSPGMRESSetRestart;
        _KSPSetFromOptions KSPSetFromOptions;
        _KSPSetOperators KSPSetOperators;
        _PCSetType PCSetType;
        _PCShellSetApply PCShellSetApply;
        _PCShellSetContext PCShellSetContext;
        _PCShellGetContext PCShellGetContext;
        _PCShellSetName PCShellSetName;
        _PCFactorSetMatSolverType PCFactorSetMatSolverType;
        _PCFactorSetMatOrderingType PCFactorSetMatOrderingType;
        _PCFactorSetLevels PCFactorSetLevels;
        _PCASMSetOverlap PCASMSetOverlap;
        _PCASMSetType PCASMSetType;
        _PCASMSetLocalType PCASMSetLocalType;
        _PCASMGetSubKSP PCASMGetSubKSP;
        _PCFieldSplitSetIS PCFieldSplitSetIS;
        _PCMGSetLevels PCMGSetLevels;
        _PCMGSetInterpolation PCMGSetInterpolation;
        _PCMGSetRestriction PCMGSetRestriction;
        _PCMGGetSmoother PCMGGetSmoother;
        _PCMGSetGalerkin PCMGSetGalerkin;
#pragma warning restore 649

        internal _PetscInitializeNoArguments InitializeNoArguments => PetscInitializeNoArguments;
        internal _PetscInitialized Initialized => PetscInitialized;
        internal _PetscFinalize Finalize => PetscFinalize;
        internal _PetscGetVersion GetVersion => PetscGetVersion;
        internal _PetscErrorMessage ErrorMessage => PetscErrorMessage;
        internal _PetscPushErrorHandler PushErrorHandler => PetscPushErrorHandler;
        internal _PetscPopErrorHandler PopErrorHandler => PetscPopErrorHandler;
        internal _PetscReturnErrorHandler ReturnErrorHandler => PetscReturnErrorHandler;
        internal _PetscOptionsCreate OptionsCreate => PetscOptionsCreate;
        internal _PetscOptionsDestroy OptionsDestroy => PetscOptionsDestroy;
        internal _PetscOptionsInsertString OptionsInsertString => PetscOptionsInsertString;
        internal _PetscOptionsClear OptionsClear => PetscOptionsClear;
        internal _PetscObjectSetOptions ObjectSetOptions => PetscObjectSetOptions;
        internal _PetscObjectTypeCompare ObjectTypeCompare => PetscObjectTypeCompare;
        internal _ISCreateGeneral CreateGeneralIndexSet => ISCreateGeneral;
        internal _ISCreateStride CreateStrideIndexSet => ISCreateStride;
        internal _ISDestroy DestroyIndexSet => ISDestroy;
        internal _VecCreateMPI CreateMpiVector => VecCreateMPI;
        internal _VecCreateSeq CreateSequentialVector => VecCreateSeq;
        internal _VecDestroy DestroyVector => VecDestroy;
        internal _VecSet SetVector => VecSet;
        internal _VecGetArray GetVectorArray => VecGetArray;
        internal _VecGetArrayRead GetVectorArrayRead => VecGetArrayRead;
        internal _VecRestoreArray RestoreVectorArray => VecRestoreArray;
        internal _VecRestoreArrayRead RestoreVectorArrayRead => VecRestoreArrayRead;
        internal _VecScatterCreate CreateVectorScatter => VecScatterCreate;
        internal _VecScatterBegin BeginVectorScatter => VecScatterBegin;
        internal _VecScatterEnd EndVectorScatter => VecScatterEnd;
        internal _VecScatterDestroy DestroyVectorScatter => VecScatterDestroy;
        internal _MatCreateAIJ CreateAijMatrix => MatCreateAIJ;
        internal _MatSeqAIJSetPreallocationCSR SetSequentialAijPreallocationCsr => MatSeqAIJSetPreallocationCSR;
        internal _MatMPIAIJSetPreallocationCSR SetMpiAijPreallocationCsr => MatMPIAIJSetPreallocationCSR;
        internal _MatCreateMPIAIJWithSplitArrays CreateMpiAijWithSplitArrays => MatCreateMPIAIJWithSplitArrays;
        internal _MatCreateShell CreateShellMatrix => MatCreateShell;
        internal _MatDestroy DestroyMatrix => MatDestroy;
        internal _MatSetValues SetMatrixValues => MatSetValues;
        internal _MatAssemblyBegin BeginMatrixAssembly => MatAssemblyBegin;
        internal _MatAssemblyEnd EndMatrixAssembly => MatAssemblyEnd;
        internal _MatGetOwnershipRange GetMatrixOwnershipRange => MatGetOwnershipRange;
        internal _MatCreateVecs CreateMatrixVectors => MatCreateVecs;
        internal _MatShellGetContext GetMatrixShellContext => MatShellGetContext;
        internal _MatShellSetOperation SetMatrixShellOperation => MatShellSetOperation;
        internal _KSPCreate CreateKsp => KSPCreate;
        internal _KSPDestroy DestroyKsp => KSPDestroy;
        internal _KSPSetType SetKspType => KSPSetType;
        internal _KSPSetUp SetUpKsp => KSPSetUp;
        internal _KSPSolve SolveKsp => KSPSolve;
        internal _KSPSetTolerances SetKspTolerances => KSPSetTolerances;
        internal _KSPSetInitialGuessNonzero SetKspInitialGuessNonzero => KSPSetInitialGuessNonzero;
        internal _KSPGetResidualNorm GetKspResidualNorm => KSPGetResidualNorm;
        internal _KSPGetIterationNumber GetKspIterationNumber => KSPGetIterationNumber;
        internal _KSPGetConvergedReason GetKspConvergedReason => KSPGetConvergedReason;
        internal _KSPGetPC GetKspPreconditioner => KSPGetPC;
        internal _KSPMonitorSet SetKspMonitor => KSPMonitorSet;
        internal _KSPGMRESSetRestart SetKspGmresRestart => KSPGMRESSetRestart;
        internal _KSPSetFromOptions SetKspFromOptions => KSPSetFromOptions;
        internal _KSPSetOperators SetKspOperators => KSPSetOperators;
        internal _PCSetType SetPcType => PCSetType;
        internal _PCShellSetApply SetPcShellApply => PCShellSetApply;
        internal _PCShellSetContext SetPcShellContext => PCShellSetContext;
        internal _PCShellGetContext GetPcShellContext => PCShellGetContext;
        internal _PCShellSetName SetPcShellName => PCShellSetName;
        internal _PCFactorSetMatSolverType SetPcFactorMatSolverType => PCFactorSetMatSolverType;
        internal _PCFactorSetMatOrderingType SetPcFactorMatOrderingType => PCFactorSetMatOrderingType;
        internal _PCFactorSetLevels SetPcFactorLevels => PCFactorSetLevels;
        internal _PCASMSetOverlap SetPcAsmOverlap => PCASMSetOverlap;
        internal _PCASMSetType SetPcAsmType => PCASMSetType;
        internal _PCASMSetLocalType SetPcAsmLocalType => PCASMSetLocalType;
        internal _PCASMGetSubKSP GetPcAsmSubKsps => PCASMGetSubKSP;
        internal _PCFieldSplitSetIS SetPcFieldSplitIndexSet => PCFieldSplitSetIS;
        internal _PCMGSetLevels SetPcMgLevels => PCMGSetLevels;
        internal _PCMGSetInterpolation SetPcMgInterpolation => PCMGSetInterpolation;
        internal _PCMGSetRestriction SetPcMgRestriction => PCMGSetRestriction;
        internal _PCMGGetSmoother GetPcMgSmoother => PCMGGetSmoother;
        internal _PCMGSetGalerkin SetPcMgGalerkin => PCMGSetGalerkin;

        /// <summary>
        /// Sets one matrix entry through <see cref="SetMatrixValues"/>.
        /// </summary>
        /// <remarks>
        /// PETSc implements <c>MatSetValue</c> as a C inline function rather
        /// than an exported shared-library symbol. The BoSSS C wrapper calls
        /// that inline helper; direct managed callers must use this equivalent
        /// one-entry <c>MatSetValues</c> call instead.
        /// </remarks>
        internal unsafe int SetMatrixValue(IntPtr matrix, int row, int column, double value, int insertMode) {
            return MatSetValues(matrix, 1, (IntPtr)(&row), 1, (IntPtr)(&column), (IntPtr)(&value), insertMode);
        }
    }

    /// <summary>
    /// Process-lifetime PETSc ABI bindings for every BoSSS PETSc runtime variant.
    /// </summary>
    /// <remarks>
    /// PETSc objects are still owned and destroyed by their individual solver.
    /// Only the dynamic-library handles and their function delegates are shared.
    /// Keeping these libraries loaded avoids unloading native code while PETSc or
    /// an OpenMP/MPI dependency may retain process-global state.
    /// </remarks>
    public static class PETScNativeLibraries {

        static readonly Lazy<UnsafePETSc> s_OpenMp = Create(PETScNativeVariant.OpenMp);
        static readonly Lazy<UnsafePETSc> s_Mpi = Create(PETScNativeVariant.Mpi);
        static readonly Lazy<UnsafePETSc> s_MpiOpenMp = Create(PETScNativeVariant.MpiOpenMp);

        static Lazy<UnsafePETSc> Create(PETScNativeVariant variant) {
            return new Lazy<UnsafePETSc>(() => new UnsafePETSc(variant), LazyThreadSafetyMode.ExecutionAndPublication);
        }

        /// <summary>
        /// Loads and binds all shipped PETSc ABI variants for the rest of this
        /// process. Call this after BoSSS native-library bootstrap when PETSc
        /// startup latency should not occur in the first solve.
        /// </summary>
        public static void PreloadAll() {
            _ = s_OpenMp.Value;
            _ = s_Mpi.Value;
            _ = s_MpiOpenMp.Value;
        }

        /// <summary>
        /// Attempts to load every PETSc variant without making PETSc a mandatory
        /// runtime dependency of the host application.
        /// </summary>
        /// <param name="failureMessage">
        /// A description of every unavailable variant, or <see langword="null"/>
        /// when all variants were loaded.
        /// </param>
        /// <returns><see langword="true"/> when every variant was loaded.</returns>
        public static bool TryPreloadAll(out string failureMessage) {
            string failures = null;
            TryPreload(s_OpenMp, PETScNativeVariant.OpenMp, ref failures);
            TryPreload(s_Mpi, PETScNativeVariant.Mpi, ref failures);
            TryPreload(s_MpiOpenMp, PETScNativeVariant.MpiOpenMp, ref failures);
            failureMessage = failures;
            return failures == null;
        }

        static void TryPreload(Lazy<UnsafePETSc> binding, PETScNativeVariant variant, ref string failures) {
            try {
                _ = binding.Value;
            } catch(Exception error) {
                string failure = $"{variant}: {error.Message}";
                failures = failures == null ? failure : failures + System.Environment.NewLine + failure;
            }
        }

        /// <summary>Returns the already process-owned ABI binding for <paramref name="variant"/>.</summary>
        internal static UnsafePETSc Get(PETScNativeVariant variant) {
            switch(variant) {
                case PETScNativeVariant.OpenMp: return s_OpenMp.Value;
                case PETScNativeVariant.Mpi: return s_Mpi.Value;
                case PETScNativeVariant.MpiOpenMp: return s_MpiOpenMp.Value;
                default: throw new ArgumentOutOfRangeException(nameof(variant), variant, "A concrete PETSc native variant is required.");
            }
        }
    }
}

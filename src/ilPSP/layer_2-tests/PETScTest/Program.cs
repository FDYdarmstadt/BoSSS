using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using BoSSS.Application.SipPoisson;
using BoSSS.Foundation.IO;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PETSc;
using ilPSP.Utils;
using MPI.Wrappers;

namespace PETScTest {

    static class Program {

        static int Main(string[] args) {
            if(Array.Exists(args, argument => string.Equals(argument, "--test-auto-selection", StringComparison.Ordinal)))
                return RunAutoSelectionTests();

            TestOptions options;
            try {
                options = ParseOptions(args);
            } catch(ArgumentException error) {
                Console.Error.WriteLine(error.Message);
                return 4;
            }
            if(options.ExpectOpenMpMpiRejection && options.NativeVariant != PETScNativeVariant.OpenMp) {
                Console.Error.WriteLine("--expect-openmp-mpi-rejection is valid only with --variant petsc_omp.");
                return 5;
            }
            if(options.ExpectUnsupportedDirectConfiguration && options.SolverMode == SolverMode.Krylov) {
                Console.Error.WriteLine("--expect-unsupported-direct-configuration requires a direct solver mode.");
                return 5;
            }

            bool mpiInitialized;
            ilPSP.Environment.Bootstrap(options.BootstrapArguments, Utils.GetNativeLibraryDir(), out mpiInitialized);
            ilPSP.Environment.InitThreading(true, null);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int mpiRank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int mpiSize);
            PETScNativeVariant effectiveNativeVariant = PETScSolverBase.ResolveNativeVariant(
                options.NativeVariant,
                mpiSize,
                ilPSP.Environment.NumThreads);
            int exitCode = 1;
            try {
                if(options.ProbeMumpsLifecycle || options.ProbeKrylovLifecycle)
                    exitCode = RunPetscLifecycleProbe(
                        options.NativeVariant,
                        options.ProbeKrylovLifecycle ? SolverMode.Krylov : SolverMode.DirectLu,
                        options.PinCsrBuffers,
                        mpiRank);
                else
                    exitCode = options.RunTutorial
                        ? RunSipPoissonTutorial(options.NativeVariant, options.SolverMode, options.DirectBackend)
                        : RunSolve(options.NativeVariant, options.SolverMode, options.DirectBackend, options.LocalRows);
                if(options.ExpectOpenMpMpiRejection) {
                    Console.Error.WriteLine("Expected the OpenMP-only variant to reject multiple MPI ranks, but the solve ran.");
                    exitCode = 6;
                }
                if(options.ExpectUnsupportedDirectConfiguration) {
                    Console.Error.WriteLine("Expected the unsupported PETSc direct configuration to be rejected, but the solve ran.");
                    exitCode = 6;
                }
            } catch(NotSupportedException error) when(
                options.ExpectOpenMpMpiRejection
                && error.Message.StartsWith("The PETSc OpenMP variant uses MPIUNI", StringComparison.Ordinal)) {
                Console.WriteLine("Expected OpenMP-only multi-rank rejection observed: " + error.Message);
                exitCode = 0;
            } catch(NotSupportedException error) when(options.ExpectUnsupportedDirectConfiguration
                && error.Message.StartsWith("PETSc direct backend", StringComparison.Ordinal)) {
                Console.WriteLine("Expected unsupported PETSc direct configuration observed: " + error.Message);
                exitCode = 0;
            } catch(InvalidOperationException error) when(options.AllowUnavailableDirectBackend
                && options.SolverMode != SolverMode.Krylov
                && error.Message.StartsWith("PETSc could not initialize the requested direct-solver backend", StringComparison.Ordinal)) {
                Console.WriteLine("Optional PETSc direct-solver backend unavailable; skipped: " + error.Message);
                exitCode = 0;
            } finally {
                if((options.ProbeMumpsLifecycle || options.ProbeKrylovLifecycle) && mpiRank == 0)
                    Console.WriteLine($"PETSc {LifecycleProbeName(options)} lifecycle probe: Finalize begin");
                using(var native = new UnsafePETSc(effectiveNativeVariant)) {
                    int petscError = native.FinalizePETScRuntime();
                    if(petscError != 0) {
                        Console.Error.WriteLine($"PETSc finalization failed with error {petscError}.");
                        exitCode = 3;
                    }
                }
                if((options.ProbeMumpsLifecycle || options.ProbeKrylovLifecycle) && mpiRank == 0)
                    Console.WriteLine($"PETSc {LifecycleProbeName(options)} lifecycle probe: Finalize complete");
                if(mpiInitialized)
                    csMPI.Raw.mpiFinalize();
            }
            return exitCode;
        }

        sealed class TestOptions {
            public PETScNativeVariant NativeVariant = PETScNativeVariant.Auto;
            public SolverMode SolverMode = SolverMode.Krylov;
            public PETScDirectBackend DirectBackend = PETScDirectBackend.Default;
            public int LocalRows = 4;
            public bool ExpectOpenMpMpiRejection;
            public bool ExpectUnsupportedDirectConfiguration;
            public bool AllowUnavailableDirectBackend;
            public bool RunTutorial;
            public bool ProbeMumpsLifecycle;
            public bool ProbeKrylovLifecycle;
            public bool PinCsrBuffers;
            public string[] BootstrapArguments;
        }

        static TestOptions ParseOptions(string[] args) {
            var options = new TestOptions();
            var bootstrapArguments = new List<string>();
            for(int i = 0; i < args.Length; i++) {
                switch(args[i]) {
                    case "--variant":
                        options.NativeVariant = ParseNativeVariant(NextOptionValue(args, ref i));
                        break;
                    case "--solver":
                        options.SolverMode = ParseSolverMode(NextOptionValue(args, ref i));
                        break;
                    case "--backend":
                        options.DirectBackend = ParseDirectBackend(NextOptionValue(args, ref i));
                        break;
                    case "--local-rows":
                        options.LocalRows = ParseLocalRows(NextOptionValue(args, ref i));
                        break;
                    case "--expect-openmp-mpi-rejection":
                        options.ExpectOpenMpMpiRejection = true;
                        break;
                    case "--expect-unsupported-direct-configuration":
                        options.ExpectUnsupportedDirectConfiguration = true;
                        break;
                    case "--allow-unavailable-direct-backend":
                        options.AllowUnavailableDirectBackend = true;
                        break;
                    case "--tutorial":
                        options.RunTutorial = true;
                        break;
                    case "--probe-mumps-lifecycle":
                        options.ProbeMumpsLifecycle = true;
                        options.SolverMode = SolverMode.DirectLu;
                        options.DirectBackend = PETScDirectBackend.Mumps;
                        break;
                    case "--probe-krylov-lifecycle":
                        options.ProbeKrylovLifecycle = true;
                        options.SolverMode = SolverMode.Krylov;
                        break;
                    case "--pin-csr":
                        options.PinCsrBuffers = true;
                        break;
                    default:
                        bootstrapArguments.Add(args[i]);
                        break;
                }
            }
            options.BootstrapArguments = bootstrapArguments.ToArray();
            return options;
        }

        static string LifecycleProbeName(TestOptions options) {
            return options.ProbeKrylovLifecycle ? "Krylov" : "MUMPS";
        }

        static string NextOptionValue(string[] args, ref int index) {
            if(index + 1 >= args.Length)
                throw new ArgumentException($"Missing value after {args[index]}.");
            return args[++index];
        }

        static PETScNativeVariant ParseNativeVariant(string value) {
            switch((value ?? "auto").Trim().ToLowerInvariant()) {
                case "auto":
                    return PETScNativeVariant.Auto;
                case "petsc_mpi":
                case "mpi":
                    return PETScNativeVariant.Mpi;
                case "petsc_omp":
                case "omp":
                    return PETScNativeVariant.OpenMp;
                case "petsc_mpi_omp":
                case "mpi_omp":
                    return PETScNativeVariant.MpiOpenMp;
                default:
                    throw new ArgumentException($"Unsupported PETSc variant '{value}'. Expected auto, petsc_mpi, petsc_omp, or petsc_mpi_omp.");
            }
        }

        static int RunAutoSelectionTests() {
            var cases = new[] {
                (MpiSize: 1, NumThreads: 1, Expected: PETScNativeVariant.OpenMp),
                (MpiSize: 1, NumThreads: 4, Expected: PETScNativeVariant.OpenMp),
                (MpiSize: 4, NumThreads: 1, Expected: PETScNativeVariant.Mpi),
                (MpiSize: 4, NumThreads: 4, Expected: PETScNativeVariant.MpiOpenMp)
            };
            foreach(var testCase in cases) {
                PETScNativeVariant actual = PETScSolverBase.ResolveNativeVariant(
                    PETScNativeVariant.Auto,
                    testCase.MpiSize,
                    testCase.NumThreads);
                Console.WriteLine($"PETSc auto selection: ranks={testCase.MpiSize}, threads={testCase.NumThreads}, variant={actual}");
                if(actual != testCase.Expected)
                    return 1;
            }
            return 0;
        }

        enum SolverMode {
            Krylov,
            DirectLu,
            DirectCholesky
        }

        static SolverMode ParseSolverMode(string value) {
            switch((value ?? "krylov").Trim().ToLowerInvariant()) {
                case "krylov":
                    return SolverMode.Krylov;
                case "direct_lu":
                    return SolverMode.DirectLu;
                case "direct_cholesky":
                    return SolverMode.DirectCholesky;
                default:
                    throw new ArgumentException($"Unsupported PETSc solver mode '{value}'. Expected krylov, direct_lu, or direct_cholesky.");
            }
        }

        static int ParseLocalRows(string value) {
            if(string.IsNullOrWhiteSpace(value))
                return 4;

            int localRows;
            if(!int.TryParse(value, out localRows) || localRows <= 0)
                throw new ArgumentException($"Unsupported local row count '{value}'. Expected a positive integer.");
            return localRows;
        }

        static PETScDirectBackend ParseDirectBackend(string value) {
            switch((value ?? "default").Trim().ToLowerInvariant()) {
                case "default":
                    return PETScDirectBackend.Default;
                case "mumps":
                    return PETScDirectBackend.Mumps;
                case "superlu_dist":
                case "superludist":
                    return PETScDirectBackend.SuperLuDist;
                case "umfpack":
                    return PETScDirectBackend.Umfpack;
                default:
                    throw new ArgumentException($"Unsupported PETSc direct backend '{value}'. Expected default, mumps, superlu_dist, or umfpack.");
            }
        }

        static int RunSipPoissonTutorial(PETScNativeVariant nativeVariant, SolverMode solverMode, PETScDirectBackend directBackend) {
            using(var solver = new SipPoissonMain()) {
                var control = SipHardcodedControl.TestCartesian1(xRes: 8, yRes: 4, pDG: 1);
                if(solverMode == SolverMode.Krylov) {
                    control.LinearSolver = new PETScKrylovSolverConfig {
                        NativeVariant = nativeVariant,
                        PETScOptions = "-ksp_type cg -pc_type jacobi",
                        ConvergenceCriterion = 1.0e-10,
                        MaxSolverIterations = 2000
                    };
                } else {
                    control.LinearSolver = new PETScDirectSolverConfig {
                        NativeVariant = nativeVariant,
                        Factorization = solverMode == SolverMode.DirectLu
                            ? PETScDirectFactorization.LU
                            : PETScDirectFactorization.Cholesky,
                        Backend = directBackend
                    };
                }

                solver.Init(control);
                solver.RunSolverMode();
                double error = (double)solver.QueryHandler.QueryResults["SolL2err"];
                Console.WriteLine($"PETSc SIP Poisson tutorial: solver={solverMode}, backend={directBackend}, l2_error={error:E3}");
                return error <= 5.0e-9 ? 0 : 1;
            }
        }

        static int RunPetscLifecycleProbe(PETScNativeVariant nativeVariant, SolverMode solverMode, bool pinCsrBuffers, int rank) {
            var partition = new Partitioning(1, csMPI.Raw._COMM.WORLD);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);
            var matrix = new MsrMatrix(partition);
            matrix[partition.i0, partition.i0] = 1.0;
            var solution = new double[1];
            var rightHandSide = new[] { 1.0 };
            string probeName = solverMode == SolverMode.Krylov ? "Krylov" : "MUMPS";
            PETScSolverBase solver = solverMode == SolverMode.Krylov
                ? (PETScSolverBase)new PETScKrylovSolver {
                    NativeVariant = nativeVariant,
                    Options = "-ksp_type cg -pc_type jacobi -ksp_rtol 1e-12",
                    PinManagedCsrBuffers = pinCsrBuffers
                }
                : new PETScDirectSolver {
                    NativeVariant = nativeVariant,
                    Factorization = PETScDirectFactorization.LU,
                    Backend = PETScDirectBackend.Mumps,
                    PinManagedCsrBuffers = pinCsrBuffers
                };

            try {
                Console.WriteLine(
                    $"PETSc {probeName} lifecycle probe input: rank={rank}, size={size}, localRows=1, "
                    + $"globalRows={partition.TotalLength}, globalRowStart={partition.i0}, "
                    + $"rowPointers=[0,1], columns=[{partition.i0}], values=[1], rhs=[1]");
                if(rank == 0)
                    Console.WriteLine($"PETSc {probeName} lifecycle probe: DefineMatrix begin (first PETSc solver in process; pinned_csr={pinCsrBuffers})");
                solver.DefineMatrix(matrix);
                if(rank == 0)
                    DumpLoadedNativeLibraries();
                if(rank == 0)
                    Console.WriteLine($"PETSc {probeName} lifecycle probe: DefineMatrix complete; KSPSolve begin");
                SolverResult result = solver.Solve(solution, rightHandSide);
                if(rank == 0)
                    Console.WriteLine($"PETSc {probeName} lifecycle probe: KSPSolve complete; converged={result.Converged}, solution={solution[0]:R}");
                return result.Converged && Math.Abs(solution[0] - 1.0) <= 1.0e-10 ? 0 : 1;
            } finally {
                if(rank == 0)
                    Console.WriteLine($"PETSc {probeName} lifecycle probe: Dispose begin");
                solver.Dispose();
                if(rank == 0)
                    Console.WriteLine($"PETSc {probeName} lifecycle probe: Dispose complete");
            }
        }

        static void DumpLoadedNativeLibraries() {
            const string mapsFile = "/proc/self/maps";
            if(!File.Exists(mapsFile)) {
                Console.WriteLine("PETSc MUMPS lifecycle probe: /proc/self/maps is unavailable.");
                return;
            }

            string[] markers = { "mpi", "petsc", "mumps", "scalapack", "blas", "lapack", "gfortran", "quadmath" };
            string[] libraries = File.ReadLines(mapsFile)
                .Select(line => line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries).LastOrDefault())
                .Where(path => path != null && path.StartsWith("/", StringComparison.Ordinal))
                .Where(path => markers.Any(marker => path.IndexOf(marker, StringComparison.OrdinalIgnoreCase) >= 0))
                .Distinct()
                .OrderBy(path => path, StringComparer.Ordinal)
                .ToArray();

            Console.WriteLine($"PETSc MUMPS lifecycle probe: loaded native libraries ({libraries.Length}):");
            foreach(string library in libraries)
                Console.WriteLine("  " + library);
        }

        static int RunSolve(PETScNativeVariant nativeVariant, SolverMode solverMode, PETScDirectBackend directBackend, int localRows) {
            var partition = new Partitioning(localRows, csMPI.Raw._COMM.WORLD);
            var matrix = new MsrMatrix(partition);
            long first = partition.i0;
            long globalRows = partition.TotalLength;

            for(int localRow = 0; localRow < localRows; localRow++) {
                long row = first + localRow;
                matrix[row, row] = 4.0;
                if(row > 0)
                    matrix[row, row - 1] = -1.0;
                if(row + 1 < globalRows)
                    matrix[row, row + 1] = -1.0;
            }

            var expected = new double[localRows];
            for(int i = 0; i < localRows; i++)
                expected[i] = 1.0;
            var rightHandSide = new double[localRows];
            for(int i = 0; i < localRows; i++) {
                long globalRow = first + i;
                rightHandSide[i] = globalRow == 0 || globalRow + 1 == globalRows ? 3.0 : 2.0;
            }
            var solution = new double[localRows];

            SolverResult result;
            double residual;
            int reason;
            if(solverMode == SolverMode.Krylov) {
                using(var solver = new PETScKrylovSolver {
                    NativeVariant = nativeVariant,
                    Options = $"-ksp_type cg -pc_type jacobi -ksp_rtol 1e-12 -ksp_atol 1e-50 -ksp_max_it {Math.Max(1000L, globalRows + 10)}"
                }) {
                    solver.DefineMatrix(matrix);
                    result = solver.Solve(solution, rightHandSide);
                    residual = solver.ResidualNorm;
                    reason = solver.ConvergenceReason;
                }
            } else {
                using(var solver = new PETScDirectSolver {
                    NativeVariant = nativeVariant,
                    Factorization = solverMode == SolverMode.DirectLu
                        ? PETScDirectFactorization.LU
                        : PETScDirectFactorization.Cholesky,
                    Backend = directBackend
                }) {
                    solver.DefineMatrix(matrix);
                    result = solver.Solve(solution, rightHandSide);
                    residual = solver.ResidualNorm;
                    reason = solver.ConvergenceReason;
                }
            }

            double localError = 0.0;
            for(int i = 0; i < localRows; i++)
                localError = Math.Max(localError, Math.Abs(solution[i] - expected[i]));
            double globalError = localError.MPIMax(csMPI.Raw._COMM.WORLD);

            int rank;
            int size;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
            PETScNativeVariant effectiveNativeVariant = PETScSolverBase.ResolveNativeVariant(
                nativeVariant,
                size,
                ilPSP.Environment.NumThreads);
            if(rank == 0)
                Console.WriteLine($"PETSc managed smoke test: requested_variant={nativeVariant}, variant={effectiveNativeVariant}, solver={solverMode}, backend={directBackend}, local_rows={localRows}, global_rows={globalRows}, ranks={size}, threads={ilPSP.Environment.NumThreads}, converged={result.Converged}, reason={reason}, iterations={result.NoOfIterations}, residual={residual:E3}, max_error={globalError:E3}");

            return result.Converged && globalError < 1.0e-10 ? 0 : 1;
        }
    }
}

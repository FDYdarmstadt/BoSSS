using System;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PETSc;
using ilPSP.Utils;
using MPI.Wrappers;

namespace PETScTest {

    static class Program {

        static int Main(string[] args) {
            PETScNativeVariant nativeVariant;
            SolverMode solverMode;
            try {
                nativeVariant = ParseNativeVariant(System.Environment.GetEnvironmentVariable("BOSSS_PETSC_VARIANT"));
                solverMode = ParseSolverMode(System.Environment.GetEnvironmentVariable("BOSSS_PETSC_SOLVER_MODE"));
            } catch(ArgumentException error) {
                Console.Error.WriteLine(error.Message);
                return 4;
            }
            bool expectOpenMpMpiRejection = string.Equals(
                System.Environment.GetEnvironmentVariable("BOSSS_PETSC_EXPECT_OPENMP_MPI_REJECTION"),
                "1",
                StringComparison.Ordinal);
            if(expectOpenMpMpiRejection && nativeVariant != PETScNativeVariant.OpenMp) {
                Console.Error.WriteLine("BOSSS_PETSC_EXPECT_OPENMP_MPI_REJECTION is valid only for petsc_omp.");
                return 5;
            }
            string nativeDirectory = System.Environment.GetEnvironmentVariable("BOSSS_NATIVE_DIR");
            if(string.IsNullOrWhiteSpace(nativeDirectory)) {
                Console.Error.WriteLine("BOSSS_NATIVE_DIR must point to the directory containing the selected libBoSSSnative_petsc_*.so.");
                return 2;
            }

            bool mpiInitialized;
            ilPSP.Environment.Bootstrap(args, nativeDirectory, out mpiInitialized);
            int exitCode = 1;
            try {
                exitCode = RunSolve(nativeVariant, solverMode);
                if(expectOpenMpMpiRejection) {
                    Console.Error.WriteLine("Expected the OpenMP-only variant to reject multiple MPI ranks, but the solve ran.");
                    exitCode = 6;
                }
            } catch(NotSupportedException error) when(
                expectOpenMpMpiRejection
                && error.Message.StartsWith("The PETSc OpenMP variant uses MPIUNI", StringComparison.Ordinal)) {
                Console.WriteLine("Expected OpenMP-only multi-rank rejection observed: " + error.Message);
                exitCode = 0;
            } finally {
                using(var native = new UnsafePETSc(nativeVariant)) {
                    int petscError = native.FinalizePETScRuntime();
                    if(petscError != 0) {
                        Console.Error.WriteLine($"PETSc finalization failed with error {petscError}.");
                        exitCode = 3;
                    }
                }
                if(mpiInitialized)
                    csMPI.Raw.mpiFinalize();
            }
            return exitCode;
        }

        static PETScNativeVariant ParseNativeVariant(string value) {
            switch((value ?? "petsc_mpi").Trim().ToLowerInvariant()) {
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
                    throw new ArgumentException($"Unsupported BOSSS_PETSC_VARIANT '{value}'. Expected petsc_mpi, petsc_omp, or petsc_mpi_omp.");
            }
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
                    throw new ArgumentException($"Unsupported BOSSS_PETSC_SOLVER_MODE '{value}'. Expected krylov, direct_lu, or direct_cholesky.");
            }
        }

        static int RunSolve(PETScNativeVariant nativeVariant, SolverMode solverMode) {
            const int localRows = 4;
            var partition = new Partitioning(localRows, csMPI.Raw._COMM.WORLD);
            var matrix = new MsrMatrix(partition);
            long first = partition.i0;
            long globalRows = partition.TotalLength;

            for(int localRow = 0; localRow < localRows; localRow++) {
                long row = first + localRow;
                matrix[row, row] = 2.0;
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
                rightHandSide[i] = globalRow == 0 || globalRow + 1 == globalRows ? 1.0 : 0.0;
            }
            var solution = new double[localRows];

            SolverResult result;
            double residual;
            int reason;
            if(solverMode == SolverMode.Krylov) {
                using(var solver = new PETScSolver {
                    NativeVariant = nativeVariant,
                    Options = "-ksp_type cg -pc_type jacobi -ksp_rtol 1e-12 -ksp_atol 1e-50 -ksp_max_it 1000"
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
                    Backend = PETScDirectBackend.Default
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
            if(rank == 0)
                Console.WriteLine($"PETSc managed smoke test: variant={nativeVariant}, solver={solverMode}, ranks={size}, OMP_NUM_THREADS={System.Environment.GetEnvironmentVariable("OMP_NUM_THREADS") ?? "unset"}, converged={result.Converged}, reason={reason}, iterations={result.NoOfIterations}, residual={residual:E3}, max_error={globalError:E3}");

            return result.Converged && globalError < 1.0e-10 ? 0 : 1;
        }
    }
}

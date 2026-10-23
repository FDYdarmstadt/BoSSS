/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.GridImport;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using ilPSP.Kraypis;
using BoSSS.Solution.Tecplot;
using ilPSP.LinSolvers.MUMPS;
using ilPSP.LinSolvers.PETSc;
using MPI.Wrappers;
using MPI.Wrappers.Utils;
using System;
using System.IO;
using System.Linq;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Checks the BoSSS installation, paricularly by checking of
    /// required/optional libraries are present
    /// </summary>
    public class InstallationChecker {

        // Only rank 0 writes user-facing status and the shared diagnostic log.
        private static bool isReportingRank;

        /// <summary>
        /// Performs a check for
        /// - MPI
        /// - Required libraries
        /// - Optional libraries
        /// </summary>
        public static void CheckSetup(bool runDetailedChecks = false) {
            Console.WriteLine(
                "Congratulations! Your basic BoSSS installation seems to be working");


            int numberOfBits = IntPtr.Size * 8;
            Console.WriteLine();
            Console.WriteLine(
                "You are running {0} bit BoSSS on a {1} operating system",
                numberOfBits,
                System.Environment.OSVersion.Platform);

            string installdir = BoSSS.Foundation.IO.Utils.GetBoSSSInstallDir();
            Console.WriteLine("Installation directory: " + installdir);



            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);
            isReportingRank = rank == 0;

            Console.Out.Flush();
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            Console.WriteLine();
            Console.WriteLine("BoSSS has been started with {0} MPI processes", size);

            ilPSP.Environment.StdoutOnlyOnRank0 = false;
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            Console.Out.Flush();
            Console.WriteLine(
                "Hello from rank {0}!", rank);
            ilPSP.Environment.StdoutOnlyOnRank0 = true;

            
            Console.WriteLine();
            Console.WriteLine("Status of REQUIRED external libraries:");
            {
                // MPI should already be loaded
                CheckDynamicLibrary("MPI", () => (DynLibLoader)csMPI.Raw);
                CheckDynamicLibrary("BLAS (sequential)", () => new UnsafeDBLAS(Parallelism.SEQ));
                CheckDynamicLibrary("LAPACK (sequential)", () => new LAPACK(Parallelism.SEQ));
                CheckDynamicLibrary("BLAS (OpenMP-parallel)", () => new UnsafeDBLAS(Parallelism.OMP));
                CheckDynamicLibrary("LAPACK (OpenMPparallel)", () => new LAPACK(Parallelism.OMP));
            }


            Console.WriteLine();
            Console.WriteLine("Status of OPTIONAL external libraries:");
            {
                string prepend = String.Format("- gnuplot:").PadRight(23);
                try {
                    string path = Gnuplot.GetGnuplotPath();
                    string file = Path.GetFileName(path);
                    path = Path.GetDirectoryName(path);

                    if(isReportingRank)
                        ReportSuccess("gnuplot", file, path);
                } catch (Exception e) {
                    if(isReportingRank)
                        ReportError("gnuplot", isOptional: true, detailedMessage: e.Message);
                }

                CheckDynamicLibrary("METIS", () => new UnsafeMETIS(), isOptional: true);
                CheckDynamicLibrary("ParMETIS", () => new ParMetis(), isOptional: true);

                CheckDynamicLibrary("CGNS", () => new CgnsDriver(false), isOptional: true);
                CheckDynamicLibrary("CGNS (HDF5)", () => new CgnsDriver(true), isOptional: true);

                CheckDynamicLibrary("Tecplot", () => new UnsafeTECIO(), isOptional: true);

                CheckDynamicLibrary("MUMPS (sequential)", () => new UnsafeMUMPS(ilPSP.Parallelism.SEQ), isOptional: true);
                CheckDynamicLibrary("MUMPS (MPI)", () => new UnsafeMUMPS(ilPSP.Parallelism.MPI), isOptional: true);

                bool petscLoaded = CheckDynamicLibrary("PETSc MPI+OMP",
                    () => new UnsafePETSc(), isOptional: true);
                if(isReportingRank && petscLoaded)
                    ReportPETScLoadCoverage();

                CheckDynamicLibrary("PARDISO (Intel MKL, sequential)", () => new ilPSP.LinSolvers.PARDISO.Wrapper_MKL(ilPSP.Parallelism.SEQ), isOptional: true);
                CheckDynamicLibrary("PARDISO (Intel MKL, OMP)", () => new ilPSP.LinSolvers.PARDISO.Wrapper_MKL(ilPSP.Parallelism.OMP), isOptional: true);
                
                CheckDynamicLibrary("UMFPACK (sequential)", () => new ilPSP.LinSolvers.UMFPACK.UnsafeUMFPACK(ilPSP.Parallelism.SEQ), isOptional: true);
                CheckDynamicLibrary("UMFPACK (OMP)", () => new ilPSP.LinSolvers.UMFPACK.UnsafeUMFPACK(ilPSP.Parallelism.OMP), isOptional: true);

                CheckDynamicLibrary("PARDISO (v5)", () => new ilPSP.LinSolvers.PARDISO.Wrapper_v5(), isOptional: true);

                CheckDynamicLibrary("OpenCL", () => new ilPSP.LinSolvers.monkey.CL.cl(), isOptional: true);
                CheckDynamicLibrary("CUDA", () => new ilPSP.LinSolvers.monkey.CUDA.cu(), isOptional: true);
            }

            if(isReportingRank)
                ReportBLASProvider();

            if (log != null) {
                log.Close();
            }

            // Detailed checks use an interactive prompt and include single-rank direct solvers.
            Console.WriteLine();
            Console.WriteLine("Basic setup check completed!");
            if(runDetailedChecks) {
                Console.WriteLine("Starting detailed solver checks.");
                DetailCheckSetup();
                return;
            }
            if(size > 1) {
                return;
            }
            if(Console.IsInputRedirected) {
                Console.WriteLine("Detailed checks require an interactive console.");
                return;
            }
            Console.WriteLine("Do you wish to perform a detailed check of functionality? y/n");
            char performCheck = Console.ReadKey().KeyChar;
            if (performCheck == 'y')
            {
                Console.Write("\n");
                DetailCheckSetup();
            }

        }

        /// <summary>
        /// Helper method to check if the library returned by
        /// <paramref name="libraryFactory"/> can be loaded without problems
        /// </summary>
        /// <param name="libraryName"></param>
        /// <param name="libraryFactory"></param>
        /// <param name="isOptional"></param>
        private static bool CheckDynamicLibrary(string libraryName, Func<DynLibLoader> libraryFactory, bool isOptional = false) {
            try {
                DynLibLoader library = libraryFactory();
                string file = Path.GetFileName(library.CurrentLibraryName);

                string path = Path.GetDirectoryName(library.CurrentLibraryName);
                if (path.IsNullOrEmpty()) {
                    path = "local search path";
                }

                if(isReportingRank)
                    ReportSuccess(libraryName, file, path);
                return true;
            } catch (Exception e) {
                if(isReportingRank)
                    ReportError(libraryName, isOptional, e.Message);
                return false;
            }
        }
        
        /// <summary>
        /// Writes a success message to the console
        /// </summary>
        /// <param name="libraryName"></param>
        /// <param name="fileName"></param>
        /// <param name="path"></param>
        private static void ReportSuccess(string libraryName, string fileName, string path) {
            string prepend = String.Format("- {0}:", libraryName).PadRight(38);
            Console.Write(prepend);

            //Console.ForegroundColor = ConsoleColor.Green;
            Console.Write("Good", prepend, fileName, path);
            //Console.ResetColor();

            Console.WriteLine(" (using {0} from {1})", fileName, path);
        }

        /// <summary>
        /// Reports an optional component which is deliberately not applicable to
        /// the current execution configuration.
        /// </summary>
        private static void ReportSkipped(string libraryName, string reason) {
            string prepend = String.Format("- {0}:", libraryName).PadRight(38);
            Console.WriteLine(prepend + "Skipped (" + reason + ")");
        }

        /// <summary>
        /// BLAS provider selection (MKL/OpenBLAS) and optional native BLAS extensions;
        /// reuses the runtime discovery of <see cref="BLAS"/>, i.e. no re-calibration and no PETSc initialization.
        /// </summary>
        private static void ReportBLASProvider() {
            Console.WriteLine();
            Console.WriteLine("BLAS provider:");
            string Row(string name) => String.Format("- {0}:", name).PadRight(38);

            Console.WriteLine(Row("requested") + (BLAS.RequestedProvider ?? "(startup selection not executed)") + $" (environment variable {BLAS.ProviderEnvVar}: mkl (default), auto or openblas)");
            Console.WriteLine(Row("effective") + BLAS.Provider + " (" + BLAS.ProviderSelectionReason + ")");
            if(BLAS.LastCalibration != null) {
                foreach(var c in BLAS.LastCalibration.Cases)
                    Console.WriteLine(Row("  " + c.Name) + string.Format(System.Globalization.CultureInfo.InvariantCulture,
                        "OpenBLAS/MKL time {0:0.000} ({1})", c.Ratio, c.Threaded ? "threaded" : "sequential"));
            }
            Console.WriteLine(Row("threads") + $"{ilPSP.Environment.NumThreads} TPL thread(s), max. {ilPSP.Environment.MaxNumOpenMPthreads} OpenMP thread(s), OpenMP {(ilPSP.Environment.OpenMPenabled ? "enabled" : "disabled")}; tensor kernels always single-threaded");
            Console.WriteLine(Row("CPU affinity (this thread)") + CPUAffinity.GetCurrentThreadAffinity().ToConcatString("[", ",", "]") + $" of {CPUAffinity.TotalNumberOfCPUs} CPU(s)");
            Console.WriteLine(Row("MKL") + "Good (using " + BLAS.MKLLibraryPath + "; also used for LAPACK and PARDISO)");
            if(BLAS.OpenBLASAvailable) {
                Console.WriteLine(Row("OpenBLAS") + "Good (" + BLAS.OpenBLASLoadInfo + "; " + BLAS.OpenBLASConfig + $"; current thread count {BLAS.OpenBLASNumThreads})");
            } else {
                Console.WriteLine(Row("OpenBLAS") + "Unavailable (" + BLAS.OpenBLASUnavailableReason + ")");
            }

            // optional extensions: presence of the exports only; not used by this BoSSS version yet
            string Exports(DynamicLibraries.DynLibHandle h, params string[] names) {
                var missing = names.Where(n => DynamicLibraries.LoadSymbol(h, n, out _) == IntPtr.Zero).ToArray();
                return missing.Length == 0 ? "exports present" : "missing export(s) " + missing.ToConcatString("", ", ", "");
            }
            var mkl = BLAS.m_seq_BLAS.LibHandle;
            string Usable(bool avail, string reason) => avail ? "Good" : "Unavailable (" + reason + ")";
            Console.WriteLine(Row("MKL batch GEMM/GEMV") + Usable(MKLBatch.Available, MKLBatch.UnavailableReason) + "; not used by Multiply (no measured gain over a loop of DGEMMs)");
            Console.WriteLine(Row("MKL JIT GEMM") + Usable(MKLJit.Available, MKLJit.UnavailableReason) + (MKLJit.Available ? $"; used by Multiply for small GEMMs ({MKLJit.NumKernels} kernel(s) cached in this process)" : "; Multiply uses the blocked SIMD path / DGEMM instead"));
            Console.WriteLine(Row("BLIS (in MKL library)") + Usable(BLIS.Available, BLIS.UnavailableReason) + "; general-stride Multiply path available, but off by default (no consistent gain measured)");
            Console.WriteLine(Row("MKL fake-Intel probe") + Exports(mkl, "BOSSS_use_mkl_fakes", "BOSSS_mkl_fakes_were_called") + "; " + FakeIntelIndicator(mkl, "embedded"));
            foreach(var pre in BLAS.m_seq_BLAS.PreReqLibHandles.Where(h => h.val != IntPtr.Zero))
                Console.WriteLine(Row("  libfakeintel (prerequisite)") + FakeIntelIndicator(pre, "separate library"));
            Console.WriteLine(Row("LIBXSMM") + Usable(LIBXSMM.Available, LIBXSMM.UnavailableReason) + (LIBXSMM.Available ? " (" + LIBXSMM.LibraryPath + ")" : "") + "; not used by Multiply (MKL JIT is faster)");
        }

        /// <summary>
        /// Reads the sticky indicator `int BOSSS_mkl_fakes_were_called(void)` (prototype verified by disassembly);
        /// a zero does not prove that the fake was not used, since MKL caches its CPU dispatch.
        /// </summary>
        [System.Runtime.InteropServices.UnmanagedFunctionPointer(System.Runtime.InteropServices.CallingConvention.Cdecl)]
        private delegate int IntFunc();

        private static string FakeIntelIndicator(DynamicLibraries.DynLibHandle h, string what) {
            IntPtr p = DynamicLibraries.LoadSymbol(h, "BOSSS_mkl_fakes_were_called", out _);
            if(p == IntPtr.Zero)
                return what + ": no indicator";
            int v = System.Runtime.InteropServices.Marshal.GetDelegateForFunctionPointer<IntFunc>(p)();
            return what + (v != 0 ? $": fake Intel CPU check was called ({v})" : ": indicator 0 (not called, or MKL dispatch cached)");
        }

        private static void ReportPETScLoadCoverage() {
            Console.WriteLine("  PETSc load check (MPI+OMP): native bridge, PETSc runtime, and its compiler/OpenMP runtime DLLs.");
        }

        /// <summary>
        /// Reports that <paramref name="libraryName"/> could not be loaded and
        /// writes <paramref name="detailedMessage"/> to the log file
        /// (see <see cref="log"/>). The option <paramref name="isOptional"/> is
        /// used to display whether this is a serious issue
        /// </summary>
        /// <param name="libraryName"></param>
        /// <param name="isOptional"></param>
        /// <param name="detailedMessage"></param>
        private static void ReportError(string libraryName, bool isOptional, string detailedMessage) {
            string prepend = String.Format("- {0}:", libraryName).PadRight(38);
            Console.Write(prepend);

            //Console.ForegroundColor = isOptional ? ConsoleColor.Yellow : ConsoleColor.Red;
            Console.Write("Failed");
            //Console.ResetColor();

            Console.WriteLine(" (detailed error messages have been written {0})", LOG_FILE_NAME);
            WriteToErrorLog(libraryName, detailedMessage);
        }

        /// <summary>
        /// Name of the log file <see cref="log"/>
        /// </summary>
        private const string LOG_FILE_NAME = "BoSSSPadCheck.log";

        /// <summary>
        /// A log file where error/warning messages are written
        /// </summary>
        private static TextWriter log;

        /// <summary>
        /// Writes to <see cref="log"/>
        /// </summary>
        /// <param name="libraryName"></param>
        /// <param name="message"></param>
        private static void WriteToErrorLog(string libraryName, string message) {
            if(!isReportingRank)
                return;
            if (log == null) {
                log = new StreamWriter(LOG_FILE_NAME);
            }

            log.WriteLine("Loading {0} with the following message:", libraryName);
            log.WriteLine(message);
            log.WriteLine();
            log.Flush();
        }

        private static void DetailCheckSetup() {
            int n = 50;
            int spacing = 1;
            double[] x = new double[n];
            double[] b = new double[n];
            MultidimensionalArray M = MultidimensionalArray.Create(n, n);

            // create some array
            for(int i = 0; i< n; i++)
            {
                b[i] = 1;
                for(int j = 0; j < n; j++){
                    M[i, j] = i+Math.Pow(j,i);
                }
            }

            RunDetailedCheck("LAPACK Function DGETRF", () => {
                M.Solve(x, b);
            });

            RunDetailedCheck("BLAS Function DDOT", () => {
                BLAS.DDOT(ref n, x, ref spacing, b,ref spacing);
            });

            bool singleMpiRank = IsSingleMpiRank();
            RunDetailedCheck(singleMpiRank ? "MUMPS sequential solver" : "MUMPS MPI solver", () => {
                using(var solver = new MUMPSSolver { Parallelism = singleMpiRank ? Parallelism.SEQ : Parallelism.MPI }) {
                    SolveAndVerify(solver);
                }
            });

            RunDetailedCheck(singleMpiRank ? "PARDISO OpenMP solver" : "PARDISO OpenMP solver (per rank)", () => {
                using(var solver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver { Parallelism = Parallelism.OMP }) {
                    SolveAndVerify(solver, csMPI.Raw._COMM.SELF);
                }
            });

            RunDetailedCheck(singleMpiRank ? "UMFPACK OpenMP solver" : "UMFPACK OpenMP solver (per rank)", () => {
                using(var solver = new ilPSP.LinSolvers.UMFPACK.UMFPACKSolver { Parallelism = Parallelism.OMP }) {
                    SolveAndVerify(solver, csMPI.Raw._COMM.SELF);
                }
            });

            RunDetailedCheck("PETSc MPI+OMP Krylov solver", () => {
                Console.WriteLine("  Components: PETSc Mat/Vec/KSP; KSP=CG; PC=Jacobi; 3-by-3 diagonal solve.");
                using(var solver = new PETScKrylovSolver {
                    Options = "-ksp_type cg -pc_type jacobi -ksp_rtol 1e-12 -ksp_atol 1e-50 -ksp_max_it 100"
                }) {
                    SolveAndVerify(solver);
                }
            });
        }

        private static bool IsSingleMpiRank() {
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);
            return size == 1;
        }

        private static void RunDetailedCheck(string name, Action test) {
            Console.Write("Starting detailed test of " + name + "\n");
            try {
                test();
                Console.Write("Successfully completed " + name + "\n\n");
            } catch(Exception exception) {
                Console.Write("FAILED " + name + ": " + exception.Message + "\n\n");
            }
        }

        private static void SolveAndVerify(ISparseSolver solver) {
            SolveAndVerify(solver, csMPI.Raw._COMM.WORLD);
        }

        private static void SolveAndVerify(ISparseSolver solver, MPI_Comm mpiComm) {
            var partitioning = new Partitioning(3, mpiComm);
            var matrix = new MsrMatrix(partitioning);
            for(int i = 0; i < partitioning.LocalLength; i++) {
                long row = partitioning.i0 + i;
                matrix[row, row] = 2.0;
            }

            var solution = new double[partitioning.LocalLength];
            var rightHandSide = new double[partitioning.LocalLength];
            for(int i = 0; i < rightHandSide.Length; i++)
                rightHandSide[i] = 2.0;

            solver.DefineMatrix(matrix);
            SolverResult result = solver.Solve(solution, rightHandSide);
            if(!result.Converged)
                throw new ArithmeticException("solver reported non-convergence.");
            for(int i = 0; i < solution.Length; i++) {
                if(Math.Abs(solution[i] - 1.0) > 1.0e-10)
                    throw new ArithmeticException($"incorrect solution at local index {i}: {solution[i]:R}.");
            }
        }


        // Deprecated
        class Metis : DynLibLoader {

            public Metis() :
                base(new string[] { "metis.dll", "libBoSSSnative_seq.so", "libmetis.so" },
                     new string[3][][],
                     new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.BoSSS_Prefix, StandardMangling },
                     new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix },
                     new int[] { -1, -1, -1 }) {
            }

            static string StandardMangling(string _name) {
                return "METIS_" + _name;
            }
        }


        class ParMetis : DynLibLoader {

            public ParMetis() :
                base(new string[] { "parmetis.dll", "libparmetis.so" },
                     new string[2][][],
                     new GetNameMangling[] { StandardMangling, StandardMangling },
                     new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
                     new int[] { -1, -1 }) {
            }

            static string StandardMangling(string _name) {
                return "ParMETIS_" + _name;
            }
        }

        // Deprecated
        class Tecplot : DynLibLoader {

            public Tecplot() :
                base(new string[] { "tecio.dll","libBoSSSnative_seq.so", "libtecio.so" },
                     new string[3][][],
                     new GetNameMangling[] { NoMangling, DynLibLoader.BoSSS_Prefix, NoMangling },
                     new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix },
                     new int[] { -1, -1, -1 }) {
            }

            static string NoMangling(string _name) {
                return _name;
            }
        }

        // Deprecated
        class MUMPS : DynLibLoader {

            public MUMPS() :
                base(new string[] { "dmumps-mpi.dll","libBoSSSnative_mpi.so", "libdmumps.so" },
                     new string[3][][],
                     new GetNameMangling[] { NoMangling, DynLibLoader.BoSSS_Prefix, NoMangling },
                     new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix },
                     new int[] { -1, -1, -1 }) {
            }

            static string NoMangling(string _name) {
                return _name;
            }
        }
    }
}

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
using ilPSP.Utils;
using MPI.Wrappers;
using MPI.Wrappers.Utils;
using System;
using System.IO;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Checks the BoSSS installation, paricularly by checking of
    /// required/optional libraries are present
    /// </summary>
    public class InstallationChecker {

        /// <summary>
        /// Performs a check for
        /// - MPI
        /// - Required libraries
        /// - Optional libraries
        /// </summary>
        public static void CheckSetup() {
            Console.WriteLine(
                "Congratulations! Your basic BoSSS installation seems to be working");

            int numberOfBits = IntPtr.Size * 8;
            Console.WriteLine();
            Console.WriteLine(
                "You are running {0} bit BoSSS on a {1} operating system",
                numberOfBits,
                System.Environment.OSVersion.Platform);

            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int size);

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
                CheckDynamicLibrary("BLAS", () => new UnsafeDBLAS());
                CheckDynamicLibrary("LAPACK", () => new LAPACK());
            }


            Console.WriteLine();
            Console.WriteLine("Status of OPTIONAL external libraries:");
            {
                string prepend = String.Format("- gnuplot:").PadRight(23);
                try {
                    string path = Gnuplot.GetGnuplotPath();
                    string file = Path.GetFileName(path);
                    path = Path.GetDirectoryName(path);

                    ReportSuccess("gnuplot", file, path);
                } catch (Exception e) {
                    ReportError("gnuplot", isOptional: true, detailedMessage: e.Message);
                }

                CheckDynamicLibrary("METIS", () => new Metis(), isOptional: true);
                CheckDynamicLibrary("ParMETIS", () => new ParMetis(), isOptional: true);

                CheckDynamicLibrary("CGNS", () => new CgnsDriver(false), isOptional: true);
                CheckDynamicLibrary("CGNS (HDF5)", () => new CgnsDriver(true), isOptional: true);

                CheckDynamicLibrary("Tecplot", () => new Tecplot(), isOptional: true);

                CheckDynamicLibrary("MUMPS", () => new MUMPS(), isOptional: true);

                CheckDynamicLibrary("PARDISO (Intek MKL)", () => new ilPSP.LinSolvers.PARDISO.Wrapper_MKL(), isOptional: true);
                CheckDynamicLibrary("PARDISO (v4)", () => new ilPSP.LinSolvers.PARDISO.Wrapper_v4(), isOptional: true);
                CheckDynamicLibrary("PARDISO (v5)", () => new ilPSP.LinSolvers.PARDISO.Wrapper_v5(), isOptional: true);

                CheckDynamicLibrary("OpenCL", () => new ilPSP.LinSolvers.monkey.CL.cl(), isOptional: true);
                CheckDynamicLibrary("CUDA", () => new ilPSP.LinSolvers.monkey.CUDA.cu(), isOptional: true);
            }

            if (log != null) {
                log.Close();
            }
        }

        /// <summary>
        /// Helper method to check if the library returned by
        /// <paramref name="libraryFactory"/> can be loaded without problems
        /// </summary>
        /// <param name="libraryName"></param>
        /// <param name="libraryFactory"></param>
        /// <param name="isOptional"></param>
        private static void CheckDynamicLibrary(string libraryName, Func<DynLibLoader> libraryFactory, bool isOptional = false) {
            try {
                DynLibLoader library = libraryFactory();
                string file = Path.GetFileName(library.CurrentLibraryName);

                string path = Path.GetDirectoryName(library.CurrentLibraryName);
                if (path.IsNullOrEmpty()) {
                    path = "local search path";
                }

                ReportSuccess(libraryName, file, path);
            } catch (Exception e) {
                ReportError(libraryName, isOptional, e.Message);
            }
        }
        
        /// <summary>
        /// Writes a success message to the console
        /// </summary>
        /// <param name="libraryName"></param>
        /// <param name="fileName"></param>
        /// <param name="path"></param>
        private static void ReportSuccess(string libraryName, string fileName, string path) {
            string prepend = String.Format("- {0}:", libraryName).PadRight(23);
            Console.Write(prepend);

            Console.ForegroundColor = ConsoleColor.Green;
            Console.Write("Good", prepend, fileName, path);
            Console.ResetColor();

            Console.WriteLine(" (using {0} from {1})", fileName, path);
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
            string prepend = String.Format("- {0}:", libraryName).PadRight(23);
            Console.Write(prepend);

            Console.ForegroundColor = isOptional ? ConsoleColor.Yellow : ConsoleColor.Red;
            Console.Write("Failed");
            Console.ResetColor();

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
            if (log == null) {
                log = new StreamWriter(LOG_FILE_NAME);
            }

            log.WriteLine("Loading {0} with the following message:", libraryName);
            log.WriteLine(message);
            log.WriteLine();
            log.Flush();
        }

        class Metis : DynLibLoader {

            public Metis() :
                base(new string[] { "metis.dll", "libmetis.so" },
                     new string[2][][],
                     new GetNameMangling[] { StandardMangling, StandardMangling },
                     new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
                     new int[] { -1, -1 }) {
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


        class Tecplot : DynLibLoader {

            public Tecplot() :
                base(new string[] { "tecio.dll", "libtecio.so" },
                     new string[2][][],
                     new GetNameMangling[] { NoMangling, NoMangling },
                     new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
                     new int[] { -1, -1 }) {
            }

            static string NoMangling(string _name) {
                return _name;
            }
        }


        class MUMPS : DynLibLoader {

            public MUMPS() :
                base(new string[] { "dmumps.dll", "libdmumps.so" },
                     new string[2][][],
                     new GetNameMangling[] { NoMangling, NoMangling },
                     new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
                     new int[] { -1, -1 }) {
            }

            static string NoMangling(string _name) {
                return _name;
            }
        }
    }
}

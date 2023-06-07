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

using System;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading;
using ilPSP.Utils;
using MPI.Wrappers;

namespace ilPSP {

    /// <summary>
    /// some basic entry points into the application
    /// </summary>
    public static class Environment {

        static private bool m_BootStrapDone = false;

        [DllImport("kernel32.dll", CharSet = CharSet.Auto)]
        private static extern void SetDllDirectory(string lpPathName);


        /// <summary>
        /// checks for the native library perquisites, inits MPI
        /// </summary>
        /// <remarks>
        /// this function performs the MPI-init
        /// </remarks>
        /// <param name="CommandLineArgs">
        /// startup arguments of the app, passes to
        /// <see cref="IMPIdriver.Init"/>
        /// </param>
        /// <param name="mpiInitialized">
        /// on exit, true, if <see cref="IMPIdriver.Init"/> was called within
        /// this function; (i.e. this is the first call to this method in the
        /// whole application).
        /// </param>
        /// <param name="__nativeDir">
        /// primary directory to search for native libraries. 
        /// </param>
        /// <returns>
        /// File directory for native files.
        /// </returns>
        public static string Bootstrap(string[] CommandLineArgs, string __nativeDir, out bool mpiInitialized) {
            //be forgiving on multiple calls
            mpiInitialized = false;
            string ret = "";
            if (m_BootStrapDone == true) {
                return ret;
            }


            var nativeDir = __nativeDir.IsEmptyOrWhite() ? default : new DirectoryInfo(__nativeDir);

            StdOut = new DuplicatingTextWriter(new StreamWriter(Console.OpenStandardOutput()), 25, false);
            Console.SetOut(StdOut);
            StdErr = new DuplicatingTextWriter(new StreamWriter(Console.OpenStandardError()), 1, false);
            Console.SetError(StdErr);
            

            //Console.WriteLine("bootstrapping necessary.");
            if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                // ++++++++++
                // MS windows 
                // ++++++++++

                
                // error
                // =====

                if (nativeDir == null || !nativeDir.Exists)
                    throw new ApplicationException("unable to do native library bootstrapping: missing directory 'x86' or 'amd64'");


                // set location for native ilPSP libraries
                // =======================================
                
                // search for the right Dll's (either 64 or 32 Bit)
                SetDllDirectory(nativeDir.FullName);
                MPI.Wrappers.Utils.DynLibLoader.PrimaryLibrarySearchPath = nativeDir.FullName;
                ret = nativeDir.FullName;

                m_BootStrapDone = true;

            } else if (System.Environment.OSVersion.Platform == PlatformID.Unix || System.Environment.OSVersion.Platform == PlatformID.MacOSX) {
               

                if (nativeDir == null || !nativeDir.Exists)
                    throw new ApplicationException ("unable to do native library bootstrapping: missing directory 'x86' or 'amd64-openmpi'");


                // set location for native ilPSP libraries
                // =======================================
                // append the directory to LD_LIBRARY_PATH
                var ld_library_path = System.Environment.GetEnvironmentVariable ("LD_LIBRARY_PATH");
                if (ld_library_path.IsEmptyOrWhite ())
                    ld_library_path = "";
                ld_library_path = nativeDir.FullName + ":" + ld_library_path;
                System.Environment.SetEnvironmentVariable ("LD_LIBRARY_PATH", ld_library_path);
                MPI.Wrappers.Utils.DynLibLoader.PrimaryLibrarySearchPath = nativeDir.FullName;


                ret = nativeDir.FullName;

                m_BootStrapDone = true;
            } else {
                Console.WriteLine("WARNING: Unable to determine os type (MS Windows or Unix?).");
                Console.WriteLine("WARNING: No bootstrapping performed");
            }

            // MPI init
            // ========
            if (!csMPI.Raw.Initialized()) {
                csMPI.Raw.Init(CommandLineArgs);
                mpiInitialized = true;
            }
     

            // init MPI enviroment
            // ===================
            m_MpiEnv = new MPIEnviroment();
            //System.Threading.Thread.Sleep(10000);
            //Console.WriteLine("StdoutOnlyOnRank0 set to false");
            StdoutOnlyOnRank0 = true;
            NativeLibraryDir = ret;
            return ret;
        }

        /// <summary>
        /// Directory where native libraries are located
        /// </summary>
        public static string NativeLibraryDir {
            get;
            private set;
        }

        /// <summary>
        /// This text writer id hooked into the standard output stream (<see cref="Console.Out"/>),
        /// in order to provide e.g. capturing to text files.
        /// </summary>
        public static DuplicatingTextWriter StdOut {
            get;
            private set;
        }

        /// <summary>
        /// This text writer id hooked into the standard error stream (<see cref="Console.Error"/>),
        /// in order to provide e.g. capturing to text files.
        /// </summary>
        public static DuplicatingTextWriter StdErr {
            get;
            private set;
        }

        static bool m_StdoutOnlyOnRank0 = false;

        /// <summary>
        /// if true, the standard - output stream will not be visible on screen on processer with MPI rank 
        /// unequal to 0.
        /// </summary>
        public static bool StdoutOnlyOnRank0 {
            get {
                return m_StdoutOnlyOnRank0;
            }
            set {
                if(StdOut != null) {
                    m_StdoutOnlyOnRank0 = value;
                    if(m_StdoutOnlyOnRank0) {
                        StdOut.surpressStream0 = (MPIEnv.MPI_Rank != 0);
                    } else {
                        StdOut.surpressStream0 = false;
                    }
                }
            }
        }


        static bool FileExistsSafe(FileInfo fi) {
            bool exists;
            try {
                exists = File.Exists(fi.FullName);
            } catch (IOException) {
                exists = true;
            }
            return exists;
        }

        static MPIEnviroment m_MpiEnv;

        /// <summary>
        /// environment of the world communicator
        /// </summary>
        public static MPIEnviroment MPIEnv {
            get {
                return m_MpiEnv;
            }
        }


        /*
        /// <summary>
        /// (tries to) do a recursive copy of a directory
        /// </summary>
        /// <param name="srcDir"></param>
        /// <param name="dstDir"></param>
        static void CopyDirectoryRec(DirectoryInfo srcDir, DirectoryInfo dstDir) {
            FileInfo[] srcFiles = srcDir.GetFiles();


            foreach (FileInfo srcFile in srcFiles) {
                TryCopy(srcFile.FullName, Path.Combine(dstDir.FullName, srcFile.Name));
            }

            foreach (DirectoryInfo srcSubDir in srcDir.GetDirectories()) {
                DirectoryInfo dstSubDir = new DirectoryInfo(Path.Combine(dstDir.FullName, srcSubDir.Name));
                if (!dstSubDir.Exists)
                    dstSubDir.Create();
                CopyDirectoryRec(srcSubDir, dstSubDir);
            }
        }
        */
        /// <summary>
        /// Utility function which tries to copy a file from
        /// <paramref name="sourceFileName"/> to
        /// <paramref name="destFileName"/> overwriting existing files if
        /// required. Issues a warning (but proceeds as normal) if the copy
        /// process fails.
        /// </summary>
        /// <param name="sourceFileName">
        /// The path to the file to be copied
        /// </param>
        /// <param name="destFileName">The path to the destination</param>
        private static void TryCopy(string sourceFileName, string destFileName) {
            try {
                File.Copy(sourceFileName, destFileName, true);
                //Console.WriteLine("Copy: " + sourceFileName + " -> " + destFileName);
            } catch (Exception e) {
                Console.WriteLine("WARNING: Unable to copy to: '"
                    + destFileName + "': " + e.GetType().Name + " says:'" + e.Message + "'");
            }
        }

    }
}

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
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;

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
            m_MpiEnv = new MPIEnvironment();
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

        /// <summary>
        /// Number of threads used in multi-thread-parallelization;
        /// - This variable refers to the number of threads used for the C#-parts of BoSSS, e.g. the quadrature kernel.
        /// - OpenMP threading is controlled by 
        /// </summary>
        public static int NumThreads {
            get;
            set;
        } = 4;

        public static int MaxNumOpenMPthreads { 
            get; 
            private set; 
        }

        /// <summary>
        /// Can be turned on and off via <see cref="EnableOpenMP"/> and <see cref="DisableOpenMP"/>, respectively.
        /// </summary>
        public static bool OpenMPenabled {
            get {
                return !OpenMPdisabled;
            }
        }


        static bool OpenMPdisabled = false;
        static int backup_MaxNumOpenMPthreads = -1;

        /// <summary>
        /// Disable the use of OpenMP in external libraries
        /// </summary>
        public static void DisableOpenMP() {
            if(OpenMPdisabled == false) {
                OpenMPdisabled = true;
                backup_MaxNumOpenMPthreads = MaxNumOpenMPthreads;
                MaxNumOpenMPthreads = 1;
                BLAS.ActivateSEQ();
                LAPACK.ActivateSEQ();
            }
        }

        /// <summary>
        /// Enable/Re-enable the use of OpenMP in external libraries (mostly Intel MKL, which provides BLAS, LAPACK and PARDISO)
        /// </summary>
        public static void EnableOpenMP() {
            if(OpenMPdisabled) {
                MaxNumOpenMPthreads = backup_MaxNumOpenMPthreads;
                OpenMPdisabled = false;
            }

            BLAS.ActivateOMP();
            LAPACK.ActivateOMP();
        }


        public static ParallelLoopResult ParallelFor(int fromInclusive, int toExclusive, Action<int, ParallelLoopState> body, bool enablePar = false) {
            if (InParallelSection) {
                throw new ApplicationException("trying to call a ParallelFor inside of a ParallelFor");
            }

            int __Numthreads = enablePar ? NumThreads : 1;


            var options = new ParallelOptions {
                MaxDegreeOfParallelism = __Numthreads,
            };
            //ThreadPool.SetMinThreads(__Numthreads, 1);
            //ThreadPool.SetMaxThreads(__Numthreads, 2);

            try {
                InParallelSection = true;
                BLAS.ActivateSEQ();
                LAPACK.ActivateSEQ();

                return Parallel.For(fromInclusive, toExclusive, options, body);
            } finally { 
                InParallelSection = false;
                BLAS.ActivateOMP();
                LAPACK.ActivateOMP();
                //SetOMPbinding();
            }
        }

        public static void ParallelFor(int fromInclusive, int toExclusive, Action<int> body, bool enablePar = true) {
            if (InParallelSection == true) {
                for (int i = 0; i < toExclusive; i++) {
                    body(i);
                }
            } else {

                int __Numthreads = enablePar ? NumThreads : 1;


                var options = new ParallelOptions {
                    MaxDegreeOfParallelism = __Numthreads,
                };
                //ThreadPool.SetMinThreads(__Numthreads, 1);
                //ThreadPool.SetMaxThreads(__Numthreads, 2);

                try {
                    InParallelSection = true;
                    BLAS.ActivateSEQ(); // within a parallel section, we don't want BLAS/LAPACK to spawn into further threads
                    LAPACK.ActivateSEQ();

                    Parallel.For(fromInclusive, toExclusive, options, body);
                } finally {
                    InParallelSection = false;
                    BLAS.ActivateOMP(); // restore parallel 
                    LAPACK.ActivateOMP();
                    //SetOMPbinding();
                }
            }
        }

        public static ParallelLoopResult ParallelFor<TLocal>(int fromInclusive, int toExclusive, Func<TLocal> localInit, Func<int, ParallelLoopState, TLocal, TLocal> body, Action<TLocal> localFinally, bool enablePar = true) {
            if (InParallelSection) {
                throw new ApplicationException("trying to call a ParallelFor inside of a ParallelFor");
            }

            int __Numthreads = enablePar ? NumThreads : 1;

            var options = new ParallelOptions {
                MaxDegreeOfParallelism = __Numthreads,
            };
            //ThreadPool.SetMinThreads(__Numthreads, 1);
            //ThreadPool.SetMaxThreads(__Numthreads, 2);

            try {
                InParallelSection = true;
                BLAS.ActivateSEQ();
                LAPACK.ActivateSEQ();

                return Parallel.For(fromInclusive, toExclusive, options, localInit, body, localFinally);
            } finally {
                InParallelSection = false;
                BLAS.ActivateOMP();
                LAPACK.ActivateOMP();
                //SetOMPbinding();
            }
        }

        /// <summary>
        /// before we start messing with OpenMP affinity
        /// </summary>
        static System.Collections.Generic.IReadOnlyList<int> ReservedCPUsInitially = null;

        static IEnumerable<int> ReservedCPUsForThisRank = null;

        static IEnumerable<int> ReservedCPUsOnSMP = null;

        public static bool MpiJobOwnsEntireComputer => ReservedCPUsOnSMP.Count() == CPUAffinity.TotalNumberOfCPUs;
        
        public static bool MpiRnkOwnsEntireComputer => ReservedCPUsInitially.Count() == CPUAffinity.TotalNumberOfCPUs;


        public static void InitThreading(bool LookAtEnvVar, int? NumThreadsOverride) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                bool bkup = StdoutOnlyOnRank0;
                StdoutOnlyOnRank0 = false;
                tr.Info($"MPI Rank {MPIEnv.MPI_Rank}: Value for OMP_PLACES: {System.Environment.GetEnvironmentVariable("OMP_PLACES")}");
                tr.Info($"MPI Rank {MPIEnv.MPI_Rank}: Value for OMP_PROC_BIND: {System.Environment.GetEnvironmentVariable("OMP_PROC_BIND")}");


                // ===========================
                // Determine Number of Threads
                // ===========================


                if (NumThreadsOverride != null) {
                    tr.Info("API override of number of threads: " + NumThreadsOverride.Value + " (ignoring OMP_NUM_THREADS, etc.)");
                    NumThreads = NumThreadsOverride.Value;
                } else {
                    if (LookAtEnvVar) {
                        int? omp_num_threads = null;
                        try {
                            string _omp_num_treads = System.Environment.GetEnvironmentVariable("OMP_NUM_THREADS");
                            if (_omp_num_treads != null) {
                                omp_num_threads = Int32.Parse(_omp_num_treads);
                                tr.Info("OMP_NUM_THREADS = " + omp_num_threads);
                            } else {
                                tr.Info("not defined: OMP_NUM_THREADS");
                            }
                        } catch (Exception e) {
                            tr.Error("Exception parsing OMP_NUM_THREADS: " + e);
                        }

                        if (omp_num_threads != null) {
                            NumThreads = omp_num_threads.Value;
                        } else {

                            int MPIranksOnNode = int.MaxValue;
                            for (int iSMP = 0; iSMP < MPIEnv.NoOfSMPs; iSMP++) {
                                MPIranksOnNode = Math.Min(MPIranksOnNode, MPIEnv.MPIProcessesPerSMP(iSMP));
                            }

                            int num_procs_tot = System.Environment.ProcessorCount;
                            int num_procs = Math.Max(1, num_procs_tot - 2); // leave some cores for the system.
                            int num_procs_per_smp = Math.Max(1, num_procs / MPIranksOnNode);
                            tr.Info($"Failed to determine user wish for number of threads; trying to use all! System reports {num_procs_tot} CPUs, will use all but 2 for BoSSS ({num_procs} total, {num_procs_per_smp} per MPI rank, MPI ranks on current node is {MPIranksOnNode}).");

                            NumThreads = num_procs_per_smp;

                        }
                    } else {
                        tr.Info($"Using default value for number of threads ({NumThreads})");
                    }
                }

                tr.Info("Finally, setting number of OpenMP and Parallel Task Library threads to " + NumThreads);

                if (NumThreads <= 0)
                    throw new NotSupportedException($"Number of threads must be at least 1; set to {NumThreads}");


                // ===========================
                // OpenMP configuration
                // ===========================
                if(ReservedCPUsInitially == null)
                    ReservedCPUsInitially = CPUAffinity.GetAffinity().ToList().AsReadOnly();
                IEnumerable<int> ReservedCPUs = ReservedCPUsInitially.ToArray();

                //if(ReservedCPUs.Count() == 1) {
                //Debugger.Launch();
                //ReservedCPUs = CPUAffinity.GetAffinity();
                //ReservedCPUs = CPUAffinity.GetAffinity().ToConcatString("[", ", ", "]")
                //}
                if (System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                    if (System.Environment.GetEnvironmentVariable("CCP_AFFINITY").IsNonEmpty()) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // Running on MS HPC Cluster, which defines the `CCP_AFFINITY` variable
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                        var _ReservedCPUs = CPUAffinityWindows.Decode_CCP_AFFINITY();
                        bool eqalAff = _ReservedCPUs.SetEquals(ReservedCPUs);
                        string listdiffs;
                        if (!eqalAff)
                            listdiffs = " (From Win32: " + ReservedCPUs.ToConcatString("[", ", ", "]") + " from CCP_AFFINITY: " + _ReservedCPUs.ToConcatString("[", ", ", "]") + ")";
                        else
                            listdiffs = "";
                        if (eqalAff == false) {
                            tr.Error("Mismatch in CPU affinity! " + listdiffs);
                        }
                        tr.Info("Win32 reports same affinity as CPUs from CCP_AFFINITY? " + eqalAff);
                        ReservedCPUs = _ReservedCPUs;
                    }
                }
                tr.Info($"R{MPIEnv.MPI_Rank}: reserved CPUs: {ReservedCPUs.ToConcatString("[", ", ", "]")}, C# reports mask {Process.GetCurrentProcess().ProcessorAffinity:X}");

                ReservedCPUsOnSMP = CPUAffinity.CpuListOnSMP(ReservedCPUs, out bool disjoint, out bool allequal);
                if (disjoint == true && allequal == true) {
                    throw new ApplicationException("Error in algorithm.");
                }

                

                tr.Info($"MpiJobOwnsEntireComputer = {MpiJobOwnsEntireComputer}, RnkJobOwnsEntireComputer = {MpiRnkOwnsEntireComputer}");

                if (allequal) {
                    if (ReservedCPUsOnSMP.Count() >= NumThreads * MPIEnv.ProcessesOnMySMP) {
                        //
                        // Sufficient CPUs to give each MPI rank `NumThreads` CPUs
                        //


                        ReservedCPUsForThisRank = ReservedCPUsOnSMP.ToArray().GetSubVector(MPIEnv.ProcessRankOnSMP * NumThreads, NumThreads);
                        MaxNumOpenMPthreads = ReservedCPUsOnSMP.Count();
                    }

                } else if (disjoint) {

                    ReservedCPUsForThisRank = ReservedCPUs.ToArray();
                    if (ReservedCPUsForThisRank.Count() < NumThreads) {
                        tr.Error($"R{MPIEnv.MPI_Rank}: Insufficient number of CPUs: NumThreads = {NumThreads}, but got affinity to {ReservedCPUsForThisRank.ToConcatString("[", ", ", "]")}");
                    }

                    MaxNumOpenMPthreads = Math.Min(ReservedCPUsForThisRank.Count(), NumThreads);

                } else {
                    // just hope for the best
                    MKLservice.Dynamic = true;
                }

                if(ReservedCPUsForThisRank != null) {
                    tr.Info($"R{MPIEnv.MPI_Rank}: using CPUs {ReservedCPUsForThisRank.ToConcatString("[", ",", "]")} for OpenMP.");
                } else {
                    tr.Info($"R{MPIEnv.MPI_Rank}: using dynamic OpenMP tread placement.");
                }
                
               
                BLAS.ActivateOMP();
                LAPACK.ActivateOMP();
                SetOMPbinding();
                tr.Info($"R{MPIEnv.MPI_Rank}: CPU affinity after OpenMP binding: " + CPUAffinity.GetAffinity().ToConcatString("[", ", ", "]"));
                CheckOMPThreading();
                StdoutOnlyOnRank0 = bkup;
            }
        }

        private static void SetOMPbinding() {
            using (var tr = new FuncTrace("SetOMPbinding")) {
                tr.InfoToConsole = true;
                if (ReservedCPUsForThisRank == null
                    || MpiRnkOwnsEntireComputer) {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // In these cases, we might just let the OpenMP threads float
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    tr.Info($"Floating OpenMP configuration ({ReservedCPUsForThisRank?.ToConcatString("[", ", ", "]") ?? "NULL"}, MpiRnkOwnsEntireComputer = {MpiRnkOwnsEntireComputer})");

                    // just hope that dynamic thread will avoid the deadlocks.
                    MKLservice.Dynamic = true;
                    MKLservice.SetNumThreads(Math.Min(MaxNumOpenMPthreads, NumThreads));
                } else {
                    var OpenMPcpuIdx = CPUAffinity.ToOpenMpCPUindices(ReservedCPUsForThisRank.Take(Math.Min(NumThreads, MaxNumOpenMPthreads))).ToArray();
                    tr.Info($"Binding to CPUs {OpenMPcpuIdx.ToConcatString("[", ", ", "]")} configuration ({ReservedCPUsForThisRank?.ToConcatString("[", ", ", "]") ?? "NULL"}, MpiRnkOwnsEntireComputer = {MpiRnkOwnsEntireComputer})");

                    MKLservice.BindOMPthreads(OpenMPcpuIdx);
                }
            }
        }

        /// <summary>
        /// We are trying to identify if OpenMP-thread from different MPI ranks dead-lock each other;
        /// Therefore:
        /// 1. A reference measurement of a GEMM operation is performed on rank 0
        /// 2. The same operation is then performed on rank [0], [0,1], [0,1,2], ...
        /// 3. The runtime measurements are then compared to the reference measurement
        /// 4. an exception is thrown if the parallel runs take much longer than the reference run
        /// </summary>
        /// <exception cref="ApplicationException"></exception>
        static void CheckOMPThreading() {
            using (var tr = new FuncTrace()) {
                /* SOME TESTS ON LICHTENBERG: 
                 * 
                 * there seems to be no clear advantage in setting OMP_PLACES;
                 * 
                 * without any OMP_PLACES: -------------------------------------------

                01/25/2024 13:57:37  Running with 4 MPI processes 
                Working path: /work/home/fk69umer/XNSE-25jan24
                Binary path: /work/home/fk69umer/XNSE-25jan24
                Current commit hash: 1ffd4df84aec90b393ee393933b8853b88c85b7a
                User: fk69umer
                Node: mpsc0536 (ranks 0, 1, 2, 3)

                MPI Rank 2: Value for OMP_PLACES: 
                MPI Rank 3: Value for OMP_PLACES: 
                MPI Rank 3: Value for OMP_PROC_BIND: 
                MPI Rank 1: Value for OMP_PLACES: 
                MPI Rank 2: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 0: Value for OMP_PLACES: 
                MPI Rank 1: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 0: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 0: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 1: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 2: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 3: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                Ref run: (min|avg|max) : (	7.879E-02 |	8.43E-02 |	1.328E-01)
                Now, doing work on 1 ranks ...
                R0: 1 workers: (min|avg|max) : (	7.878E-02 |	7.9E-02   |	7.924E-02)  --- 		( 1E00     |	9.371E-01 |	5.97E-01)
                Now, doing work on 2 ranks ...
                R0: 2 workers: (min|avg|max) : (	7.579E-02 |	7.853E-02 |	8.01E-02)  --- 		    ( 9.62E-01 |	9.316E-01 |	6.03E-01)
                R1: 2 workers: (min|avg|max) : (	3.419E-02 |	8.897E-02 |	4.743E-01)  --- 		( 4.34E-01 |	1.055E00  |	3.57E00)
                Now, doing work on 3 ranks ...
                R1: 3 workers: (min|avg|max) : (	5.337E-02 |	7.226E-02 |	8.145E-02)  --- 		( 6.77E-01 |	8.572E-01 |	6.14E-01)
                R0: 3 workers: (min|avg|max) : (	7.728E-02 |	8.364E-02 |	9.387E-02)  --- 		( 9.81E-01 |	9.922E-01 |	7.07E-01)
                R2: 3 workers: (min|avg|max) : (	7.862E-02 |	8.826E-02 |	1.322E-01)  --- 		( 9.98E-01 |	1.047E00  |	9.96E-01)
                Now, doing work on 4 ranks ...
                R2: 4 workers: (min|avg|max) : (	7.767E-02 |	8.378E-02 |	9.692E-02)  --- 		( 9.86E-01 |	9.939E-01 |	7.3E-01)
                R0: 4 workers: (min|avg|max) : (	8.063E-02 |	1.092E-01 |	1.83E-01)  --- 		    ( 1.02E00  |	1.296E00  |	1.38E00)
                R1: 4 workers: (min|avg|max) : (	6.777E-02 |	1.142E-01 |	3.435E-01)  --- 		( 8.6E-01  |	1.354E00  |	2.59E00)
                R3: 4 workers: (min|avg|max) : (	3.407E-02 |	1.279E-01 |	6.426E-01)  --- 		( 4.32E-01 |	1.517E00  |	4.84E00)
                            -------------------------------------

                and now, with OMP_PLACES set: -----------------------------------
                01/25/2024 14:05:01  Running with 4 MPI processes 
                Working path: /work/home/fk69umer/XNSE-25jan24
                Binary path: /work/home/fk69umer/XNSE-25jan24
                Current commit hash: 1ffd4df84aec90b393ee393933b8853b88c85b7a
                User: fk69umer
                Node: mpsd0114 (ranks 0, 1, 2, 3)


                MPI Rank 3: Value for OMP_PLACES: 
                MPI Rank 1: Value for OMP_PLACES: 
                MPI Rank 2: Value for OMP_PLACES: 
                MPI Rank 3: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 1: Value for OMP_PROC_BIND: 
                MPI Rank 2: Value for OMP_PROC_BIND: 
                MPI Rank 0: Value for OMP_PLACES: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 0: Value for OMP_PROC_BIND: 
                OMP_NUM_THREADS = 4
                Finally, setting number of OpenMP and Parallel Task Library threads to 4
                MPI Rank 2: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 0: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 3: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                MPI Rank 1: assigned to CPUs: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
                R0, SMP rank 0: setting OMP_PLACES = {0,1,2,3}
                R1, SMP rank 1: setting OMP_PLACES = {4,5,6,7}
                R3, SMP rank 3: setting OMP_PLACES = {12,13,14,15}
                R2, SMP rank 2: setting OMP_PLACES = {8,9,10,11}
                Ref run: (min|avg|max) : (	3.31E-02 |	5.071E-02 |	1.917E-01)
                Now, doing work on 1 ranks ...
                R0: 1 workers: (min|avg|max) : (	3.297E-02 |	3.804E-02 |	5.443E-02)  --- 		( 9.96E-01 |	7.502E-01 |	2.84E-01)
                Now, doing work on 2 ranks ...
                R0: 2 workers: (min|avg|max) : (	3.623E-02 |	4.84E-02 |	6.778E-02)  --- 		( 1.09E00 |	9.545E-01 |	3.54E-01)
                Now, doing work on 3 ranks ...
                R1: 2 workers: (min|avg|max) : (	4.857E-02 |	7.535E-02 |	2.967E-01)  --- 		( 1.47E00 |	1.486E00 |	1.55E00)
                R1: 3 workers: (min|avg|max) : (	5.121E-02 |	5.61E-02 |	7.05E-02)  --- 		( 1.55E00 |	1.106E00 |	3.68E-01)
                R0: 3 workers: (min|avg|max) : (	4.167E-02 |	6.363E-02 |	7.679E-02)  --- 		( 1.26E00 |	1.255E00 |	4E-01)
                R2: 3 workers: (min|avg|max) : (	4.825E-02 |	7.676E-02 |	2.685E-01)  --- 		( 1.46E00 |	1.514E00 |	1.4E00)
                Now, doing work on 4 ranks ...
                R1: 4 workers: (min|avg|max) : (	5.134E-02 |	5.854E-02 |	9.231E-02)  --- 		( 1.55E00 |	1.154E00 |	4.81E-01)
                R2: 4 workers: (min|avg|max) : (	5.128E-02 |	6.41E-02 |	8.576E-02)  --- 		( 1.55E00 |	1.264E00 |	4.47E-01)
                R3: 4 workers: (min|avg|max) : (	5.178E-02 |	8.39E-02 |	3.045E-01)  --- 		( 1.56E00 |	1.655E00 |	1.59E00)
                R0: 4 workers: (min|avg|max) : (	3.749E-02 |	8.64E-02 |	1.526E-01)  --- 		( 1.13E00 |	1.704E00 |	7.96E-01)
                -------------------------------------

                */


                const int N = 2048;
                const int Runs = 5;


                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MpiSz);
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int Rank);


                var A = MultidimensionalArray.Create(N, N);
                var B = MultidimensionalArray.Create(N, N);
                var C = MultidimensionalArray.Create(N, N);

                A.Storage.FillRandom();
                B.Storage.FillRandom();
                C.Storage.FillRandom();



                (double minTime, double avgTime, double maxTime) GEMMbench() {
                    double mintime = double.MaxValue;
                    double maxtime = 0.0;

                    var RunTimes = new System.Collections.Generic.List<double>();

                    for (int i = 0; i < Runs; i++) {
                        var start = DateTime.Now;
                        A.GEMM(1.0, B, C, 0.1);
                        var end = DateTime.Now;

                        double secs = (end - start).TotalSeconds;

                        RunTimes.Add(secs);
                        mintime = Math.Min(mintime, secs);
                        maxtime = Math.Max(maxtime, secs);
                    }

                    // remove the outliers before computing the average:
                    RunTimes.Sort();
                    RunTimes.RemoveAt(0);
                    RunTimes.RemoveAt(RunTimes.Count - 1);
                    double avgTime = RunTimes.Sum() / RunTimes.Count;

                    avgTime /= Runs;

                    return (mintime, avgTime, maxtime);
                }

                (double minTime, double avgTime, double maxTime) TimeRef0 = (0, 0, 0);
                if (Rank == 0) {
                    tr.Info("Now, doing reference run on rank 0...");
                    TimeRef0 = GEMMbench();
                    tr.Info($"Ref run: (min|avg|max) : (\t{TimeRef0.minTime:0.###E-00} |\t{TimeRef0.avgTime:0.###E-00} |\t{TimeRef0.maxTime:0.###E-00})");
                }

                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                var TimeRef = TimeRef0.MPIBroadcast(0);

                for (int ranksToBench = 0; ranksToBench < MpiSz; ranksToBench++) {
                    tr.Info("Now, doing work on " + (ranksToBench + 1) + " ranks ...");

                    (double minTime, double avgTime, double maxTime) TimeX = (BLAS.MachineEps, BLAS.MachineEps, BLAS.MachineEps);
                    if (Rank <= ranksToBench) {
                        TimeX = GEMMbench();
                    }


                    double minFactor = TimeX.minTime / TimeRef.minTime;
                    double avgFactor = TimeX.avgTime / TimeRef.avgTime;
                    double maxFactor = TimeX.maxTime / TimeRef.maxTime;


                    if (minFactor.MPIMax() > 10 || maxFactor.MPIMax() > 5 || avgFactor.MPIMax() > 10) {

                        string scaling = $"R{Rank}: {ranksToBench + 1} workers: (min|avg|max) : (\t{TimeX.minTime:0.###E-00} |\t{TimeX.avgTime:0.###E-00} |\t{TimeX.maxTime:0.###E-00})  --- \t\t( {minFactor:0.##E-00} |\t{avgFactor:0.###E-00} |\t{maxFactor:0.##E-00})";
                        tr.Info("Scaling involving " + (ranksToBench + 1) + " ranks: " + scaling);
                        Console.WriteLine("Suspicious OpenMP runtime behavior: " + scaling);

                        //if(avgFactor > 7)
                        //    throw new ApplicationException("Some very slow processor detected -- maybe some OpenMP locking: " + scaling);
                    }



                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                }


            }
        }

        /// <summary>
        /// true, if the code currently runs in multi-threaded; then, further spawning into sub-threads should not occur.
        /// </summary>
        public static bool InParallelSection {
            get;
            private set;
        } = false;





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

        static MPIEnvironment m_MpiEnv;

        /// <summary>
        /// environment of the world communicator
        /// </summary>
        public static MPIEnvironment MPIEnv {
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

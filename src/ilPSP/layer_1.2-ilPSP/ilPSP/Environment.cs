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
            if(m_BootStrapDone == true) {
                return ret;
            }


            var nativeDir = __nativeDir.IsEmptyOrWhite() ? default : new DirectoryInfo(__nativeDir);

            StdOut = new DuplicatingTextWriter(new StreamWriter(Console.OpenStandardOutput()), 25, false);
            Console.SetOut(StdOut);
            StdErr = new DuplicatingTextWriter(new StreamWriter(Console.OpenStandardError()), 1, false);
            Console.SetError(StdErr);


            //Console.WriteLine("bootstrapping necessary.");
            if(System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                // ++++++++++
                // MS windows 
                // ++++++++++


                // error
                // =====

                if(nativeDir == null || !nativeDir.Exists)
                    throw new ApplicationException("unable to do native library bootstrapping: missing directory 'x86' or 'amd64'");


                // set location for native ilPSP libraries
                // =======================================

                // search for the right Dll's (either 64 or 32 Bit)
                SetDllDirectory(nativeDir.FullName);
                MPI.Wrappers.Utils.DynLibLoader.PrimaryLibrarySearchPath = nativeDir.FullName;
                ret = nativeDir.FullName;

                m_BootStrapDone = true;

            } else if(System.Environment.OSVersion.Platform == PlatformID.Unix || System.Environment.OSVersion.Platform == PlatformID.MacOSX) {


                if(nativeDir == null || !nativeDir.Exists)
                    throw new ApplicationException("unable to do native library bootstrapping: missing directory 'x86' or 'amd64-openmpi'");


                // set location for native ilPSP libraries
                // =======================================
                // append the directory to LD_LIBRARY_PATH
                var ld_library_path = System.Environment.GetEnvironmentVariable("LD_LIBRARY_PATH");
                if(ld_library_path.IsEmptyOrWhite())
                    ld_library_path = "";
                ld_library_path = nativeDir.FullName + ":" + ld_library_path;
                System.Environment.SetEnvironmentVariable("LD_LIBRARY_PATH", ld_library_path);
                MPI.Wrappers.Utils.DynLibLoader.PrimaryLibrarySearchPath = nativeDir.FullName;


                ret = nativeDir.FullName;

                m_BootStrapDone = true;
            } else {
                Console.WriteLine("WARNING: Unable to determine os type (MS Windows or Unix?).");
                Console.WriteLine("WARNING: No bootstrapping performed");
            }

            // MPI init
            // ========
            if(!csMPI.Raw.Initialized()) {
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

        /// <summary>
        /// The maximum number of OpenMP threads on the entire computer (aka. Symmetric Multi Processing Node, SMP Node), among all MPI Ranks;
        /// This is important if, e.g. in parallel PARDISO solves, the matrix is gathered on one MPI rank and all other Ranks pause;
        /// </summary>
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
            if(DisableOpenMP_becauseIsSlow)
                return;

            if(OpenMPdisabled) {
                MaxNumOpenMPthreads = backup_MaxNumOpenMPthreads;
                OpenMPdisabled = false;
            }

            BLAS.ActivateOMP();
            LAPACK.ActivateOMP();
            PinOMPthreads();
        }


        public static ParallelLoopResult ParallelFor(int fromInclusive, int toExclusive, Action<int, ParallelLoopState> body, bool enablePar = false) {
            if(InParallelSection) {
                throw new ApplicationException("trying to call a ParallelFor inside of a ParallelFor");
            }

            int __Numthreads = enablePar ? NumThreads : 1;
            PinTPLThreads();

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
                PinOMPthreads();
            }
        }

        public static void ParallelFor(int fromInclusive, int toExclusive, Action<int> body, bool enablePar = true) {
            if(InParallelSection == true || !enablePar || NumThreads <= 1) {
                for(int i = fromInclusive; i < toExclusive; i++) {
                    body(i);
                }
            } else {

                int __Numthreads = enablePar ? NumThreads : 1;
                PinTPLThreads();

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
                    PinOMPthreads();
                }
            }
        }

        public static void ParallelFor(int fromInclusive, int toExclusive, Action<int, int> body, bool enablePar = true) {
            if(InParallelSection == true || !enablePar || NumThreads <= 1) {
                body(fromInclusive, toExclusive);
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

                    void _body(int ithread) {
                        PinTPLThread(ithread);
                        int L = toExclusive - fromInclusive;
                        int i0 = (L * ithread) / __Numthreads;
                        int iE = (L * (ithread + 1)) / __Numthreads;

                        body(i0, iE);
                    }


                    Parallel.For(0, __Numthreads, options, _body);
                } finally {
                    InParallelSection = false;
                    BLAS.ActivateOMP(); // restore parallel 
                    LAPACK.ActivateOMP();
                    PinOMPthreads();
                }
            }
        }

        static bool PerformTPLthreadPinning = false;

        static void PinTPLThreads() {
            if(PerformTPLthreadPinning) {
                Parallel.For(0, ilPSP.Environment.NumThreads,
                    new ParallelOptions { MaxDegreeOfParallelism = ilPSP.Environment.NumThreads },
                    PinTPLThread);
            }
        }

        static void PinTPLThread(int ithread) {
            if(PerformTPLthreadPinning) {
                int L = DedicatedCPUsForThisRank.Length;
                if(ilPSP.Environment.NumThreads > L)
                    throw new ApplicationException("Configuration error: more threads than CPUs available");
                int skip = L - NumThreads;
                int iCpu = DedicatedCPUsForThisRank[skip + ithread]; // use the left-over CPUs **at the beginning** for spare; I assume that background threads rather grab those, resp. mpiexec is forcing them to do so.
                CPUAffinity.SetCurrentThreadAffinity(iCpu);
                //CPUAffinity.SetCurrentThreadAffinity(DedicatedCPUsForThisRank);
            }
        }

        public static void ParallelFor(int fromInclusive, int toExclusive, Action<int, int, int> body, bool enablePar = true) {
            if(InParallelSection == true || !enablePar || NumThreads <= 1) {
                body(0, fromInclusive, toExclusive);
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

                    void _body(int ithread) {
                        PinTPLThread(ithread);

                        int L = toExclusive - fromInclusive;
                        int i0 = (L * ithread) / __Numthreads;
                        int iE = (L * (ithread + 1)) / __Numthreads;

                        body(ithread, i0, iE);
                    }


                    Parallel.For(0, __Numthreads, options, _body);
                } finally {
                    InParallelSection = false;
                    BLAS.ActivateOMP(); // restore parallel 
                    LAPACK.ActivateOMP();
                    PinOMPthreads();
                }
            }
        }


        public static ParallelLoopResult ParallelFor<TLocal>(int fromInclusive, int toExclusive, Func<TLocal> localInit, Func<int, ParallelLoopState, TLocal, TLocal> body, Action<TLocal> localFinally, bool enablePar = true) {
            if(InParallelSection) {
                throw new ApplicationException("trying to call a ParallelFor inside of a ParallelFor");
            }

            int __Numthreads = enablePar ? NumThreads : 1;
            PinTPLThreads();

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
                PinOMPthreads();
            }
        }


        /// <summary>
        /// Turns of Multithread-Parallelizatin in a certain subsection
        /// </summary>
        public class SerialSection : IDisposable {
            
            int numThreadsBackup;

            
            public SerialSection() {
                if(InParallelSection) {
                    throw new ApplicationException("already in parallel section");
                }
                numThreadsBackup = NumThreads;
                NumThreads = 1;
            }

            public void Dispose() {
                NumThreads = numThreadsBackup;
            }
        }



        /// <summary>
        /// before we start messing with OpenMP affinity
        /// </summary>
        static System.Collections.Generic.IReadOnlyList<int> ReservedCPUsInitially = null;

        static int[] DedicatedCPUsForThisRank = null;

        static IEnumerable<int> ReservedCPUsOnSMP = null;

        public static bool MpiJobOwnsEntireComputer => ReservedCPUsOnSMP.Count() == CPUAffinity.TotalNumberOfCPUs;

        public static bool MpiRnkOwnsEntireComputer => MpiJobOwnsEntireComputer && MPIEnv.ProcessesOnMySMP == 1;


        public static void InitThreading(bool LookAtEnvVar, int? NumThreadsOverride) {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                StdoutOnlyOnRank0 = false;
                //tr.StdoutOnAllRanks();


                tr.Info($"MPI Rank {MPIEnv.MPI_Rank}: Value for OMP_PLACES: {System.Environment.GetEnvironmentVariable("OMP_PLACES")}");
                tr.Info($"MPI Rank {MPIEnv.MPI_Rank}: Value for OMP_PROC_BIND: {System.Environment.GetEnvironmentVariable("OMP_PROC_BIND")}");
                tr.Info($"Number of CPUs in system: {CPUAffinity.TotalNumberOfCPUs}");

                // ===========================
                // Determine Number of Threads
                // ===========================


                if(NumThreadsOverride != null) {
                    tr.Info("API override of number of threads: " + NumThreadsOverride.Value + " (ignoring OMP_NUM_THREADS, etc.)");
                    NumThreads = NumThreadsOverride.Value;
                } else {
                    int? omp_num_threads = null;
                    if(LookAtEnvVar) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // // try to infer number of threads from envvar OMP_NUM_THREADS
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        try {
                            string _omp_num_treads = System.Environment.GetEnvironmentVariable("OMP_NUM_THREADS");
                            if(_omp_num_treads != null) {
                                omp_num_threads = Int32.Parse(_omp_num_treads);
                                tr.Info("OMP_NUM_THREADS = " + omp_num_threads);
                            } else {
                                tr.Info("not defined: OMP_NUM_THREADS");
                            }
                        } catch(Exception e) {
                            tr.Error("Exception parsing OMP_NUM_THREADS: " + e);
                        }
                    }

                    if(omp_num_threads != null) {
                        NumThreads = omp_num_threads.Value;
                    } else {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++
                        // nothing defined, so use everything on the computer
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++


                        int MPIranksOnNode = int.MaxValue;
                        for(int iSMP = 0; iSMP < MPIEnv.NoOfSMPs; iSMP++) {
                            MPIranksOnNode = Math.Min(MPIranksOnNode, MPIEnv.MPIProcessesPerSMP(iSMP));
                        }

                        int num_procs_tot = System.Environment.ProcessorCount;
                        tr.Info($"System reports {num_procs_tot} CPUs, {MPIranksOnNode} MPI ranks on current node.");
                        int num_procs = Math.Max(1, num_procs_tot - 2); // leave some cores for the system.
                        int num_procs_per_process = Math.Max(1, num_procs / MPIranksOnNode);
                        if(num_procs_per_process > 1 && num_procs_per_process % 2 != 0)
                            num_procs_per_process--; // choose an even number
                        tr.Info($"Failed to determine user wish for number of threads; trying to use all! System reports {num_procs_tot} CPUs, will use all but 2 for BoSSS ({num_procs} total, {num_procs_per_process} per MPI rank, MPI ranks on current node is {MPIranksOnNode}).");

                        NumThreads = num_procs_per_process;

                    }
                }
                MaxNumOpenMPthreads = NumThreads;

                
                //else {
                //    tr.Info($"Using default value for number of threads ({NumThreads})");
                //}


                tr.Info("Finally, setting number of OpenMP and Parallel Task Library threads to " + NumThreads);

                if(NumThreads <= 0)
                    throw new NotSupportedException($"Number of threads must be at least 1; set to {NumThreads}");
                // ===========================
                // OpenMP configuration
                // ===========================
                if(ReservedCPUsInitially == null)
                    ReservedCPUsInitially = CPUAffinity.GetCurrentThreadAffinity().ToList().AsReadOnly();
                IEnumerable<int> ReservedCPUs = ReservedCPUsInitially.ToArray();

                //if(ReservedCPUs.Count() == 1) {
                //Debugger.Launch();
                //ReservedCPUs = CPUAffinity.GetAffinity();
                //ReservedCPUs = CPUAffinity.GetAffinity().ToConcatString("[", ", ", "]")
                //}
                if(System.Environment.OSVersion.Platform == PlatformID.Win32NT) {
                    if(System.Environment.GetEnvironmentVariable("CCP_AFFINITY").IsNonEmpty()) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // Running on MS HPC Cluster, which defines the `CCP_AFFINITY` variable
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        tr.Info($"CCP_AFFINITY is set as '{System.Environment.GetEnvironmentVariable("CCP_AFFINITY")}'");

                        var _ReservedCPUs = CPUAffinityWindows.Decode_CCP_AFFINITY();
                        bool eqalAff = _ReservedCPUs.SetEquals(ReservedCPUs);
                        string listdiffs;
                        if(!eqalAff)
                            listdiffs = " (From Win32: " + ReservedCPUs.ToConcatString("[", ",", "]") + " from CCP_AFFINITY: " + _ReservedCPUs.ToConcatString("[", ",", "]") + ")";
                        else
                            listdiffs = "";
                        if(eqalAff == false) {
                            tr.Error("Mismatch in CPU affinity (" + MPIEnv.MPI_Rank + "of" + MPIEnv.MPI_Size + ")! " + listdiffs);
                        }
                        tr.Info("Win32 reports same affinity as CPUs from CCP_AFFINITY? " + eqalAff);
                        //ReservedCPUs = _ReservedCPUs;
                    } else {
                        tr.Info($"CCP_AFFINITY not set");
                    }
                }
                tr.Info($"R{MPIEnv.MPI_Rank}: reserved CPUs: {ReservedCPUs.ToConcatString("[", ",", "]")}, C# reports mask {Process.GetCurrentProcess().ProcessorAffinity:X}");

                ReservedCPUsOnSMP = CPUAffinity.CpuListOnSMP(ReservedCPUs, out bool disjoint, out bool allequal);
                if(disjoint == true && allequal == true) {
                    throw new ApplicationException("Error in algorithm.");
                }
                tr.Info($"R{MPIEnv.MPI_Rank}: reserved CPUs for this MPI rank: {ReservedCPUs.ToConcatString("[", ",", "]")}, on entire SMP: {ReservedCPUsOnSMP.ToConcatString("[", ",", "]")}, disjount CPU affinities? {disjoint}");


                tr.Info($"MpiJobOwnsEntireComputer = {MpiJobOwnsEntireComputer}, RnkJobOwnsEntireComputer = {MpiRnkOwnsEntireComputer}");

                if(allequal) {
                    int l = ReservedCPUsOnSMP.Count();
                    int r = MPIEnv.ProcessRankOnSMP;
                    int s = MPIEnv.ProcessesOnMySMP;

                    int i0 = (l * r) / s, iE = (l*(r+1)) / s;

                    DedicatedCPUsForThisRank = ReservedCPUsOnSMP.ToArray().GetSubVector(i0, iE - i0);

                    /*if(ReservedCPUsOnSMP.Count() >= NumThreads * MPIEnv.ProcessesOnMySMP) {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // Sufficient CPUs to give each MPI rank `NumThreads` CPUs
                        //
                        // We might have more CPUs at hand than `NumThreads`;
                        // But, despite havening more, we only want to use `NumThreads` CPUs, 
                        // since the user only requested `NumThreads` CPUs
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        int MaxSkip = ReservedCPUsOnSMP.Count() - NumThreads * MPIEnv.ProcessesOnMySMP;
                        int skip;
                        if(MaxSkip.MPIMax() > 0) {
                            //
                            // If we have more CPUs than requested, still do some shifting, to prevent using the same CPU over and over
                            //

                            if(MPIEnv.ProcessRankOnSMP == 0) {
                                MaxSkip = rnd.Next(0, MaxSkip + 1);
                            } else {
                                MaxSkip = -1;
                            }

                            int[] _MaxSkip = new int[MPIEnv.NoOfSMPs];
                            _MaxSkip[MPIEnv.SMPrank] = MaxSkip;
                            _MaxSkip = _MaxSkip.MPIMax();
                            skip = _MaxSkip[MPIEnv.SMPrank];
                        } else {
                            skip = 0;
                        }

                        DedicatedCPUsForThisRank = ReservedCPUsOnSMP.ToArray().GetSubVector(skip + MPIEnv.ProcessRankOnSMP * NumThreads, NumThreads);
                        MaxNumOpenMPthreads = Math.Min(ReservedCPUsOnSMP.Count(), MPIEnv.ProcessesOnMySMP * NumThreads);
                    }*/

                } else if(disjoint) {

                    DedicatedCPUsForThisRank = ReservedCPUs.ToArray();
                    
                } else {
                    DedicatedCPUsForThisRank = ReservedCPUs.ToArray();
                    // just hope for the best
                    MKLservice.Dynamic = true;
                    MKLservice.SetNumThreads(NumThreads);
                    tr.Warning("Reserved CPUs for all ranks are neither equal nor disjoint; some CPUs are owned by multiple ranks, some are exclusive; Pinning will be disabled.");
                }
                if(NumThreads > DedicatedCPUsForThisRank.Count()) {
                    NumThreads = Math.Max(1, DedicatedCPUsForThisRank.Count());
                    MaxNumOpenMPthreads = NumThreads; 
                    tr.Error($"R{MPIEnv.MPI_Rank}: Insufficient number of CPUs: NumThreads = {NumThreads}, but got affinity to {DedicatedCPUsForThisRank.ToConcatString("[", ",", "]")}; reducing to {MaxNumOpenMPthreads} threads");
                }


                if(DedicatedCPUsForThisRank != null) {
                    tr.Info($"R{MPIEnv.MPI_Rank}: using CPUs {DedicatedCPUsForThisRank.ToConcatString("[", ",", "]")} for OpenMP.");
                } else {
                    tr.Info($"R{MPIEnv.MPI_Rank}: using dynamic OpenMP tread placement.");
                }

                PerformOMPthreadPinning = MPIEnv.MPI_Size > 1 && (allequal != disjoint);
                PerformTPLthreadPinning = MPIEnv.MPI_Size > 1 && (allequal != disjoint);

/*
                if(NumThreads*2 < DedicatedCPUsForThisRank.Length) {
                    PerformOMPthreadPinning = false;
                    PerformTPLthreadPinning = false;
                }
*/

                tr.Info($"R{MPIEnv.MPI_Rank}: TPL thread pinning: {PerformTPLthreadPinning}, OMP thread pinning: {PerformOMPthreadPinning}");

                BLAS.ActivateOMP();
                LAPACK.ActivateOMP();
                PinTPLThreads();
                PinOMPthreads();
                StdoutOnlyOnRank0 = true;
                tr.Info($"R{MPIEnv.MPI_Rank}: CPU affinity after OpenMP binding: " + CPUAffinity.GetCurrentThreadAffinity().ToConcatString("[", ",", "]"));
            }
        }

        static bool PerformOMPthreadPinning = false;


        private static void PinOMPthreads() {



            if(PerformOMPthreadPinning) {
                //var cpus = DedicatedCPUsForThisRank.GetSubVector(DedicatedCPUsForThisRank.Length - NumThreads, NumThreads); // use the left-over CPUs **at the beginning** for spare; I assume that background threads rather grab those.
                //MKLservice.BindOMPthreads_1To1(cpus);
            }
            
            /*{
                //tr.InfoToConsole = true;
                if (DedicatedCPUsForThisRank == null || MpiRnkOwnsEntireComputer) {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // In these cases, we might just let the OpenMP threads float
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    
                    // just hope that dynamic thread will avoid the deadlocks.
                    MKLservice.SetNumThreads(Math.Min(MaxNumOpenMPthreads, NumThreads));
                    MKLservice.Dynamic = true;
                } else {
                    int[] OpenMPcpuIdx;

                    int L = DedicatedCPUsForThisRank.Count();
                    int Nt = Math.Min(NumThreads, MaxNumOpenMPthreads);

                    if (L > Nt) {

                        int skip = 0; // rnd.Next(0, L - Nt + 1);
                        if (skip + Nt > L) {
                            throw new ApplicationException("skipping done wrong");
                        }

                        OpenMPcpuIdx = CPUAffinity.ToOpenMpCPUindices(DedicatedCPUsForThisRank.Skip(skip).Take(Math.Min(NumThreads, MaxNumOpenMPthreads))).ToArray();
                    } else {
                        OpenMPcpuIdx = CPUAffinity.ToOpenMpCPUindices(DedicatedCPUsForThisRank).ToArray();
                    }

                    if(OMPbindingStrategy == null) {
                        //OMPbindingStrategy = OnlinePerformanceMeasurement.FindBestOMPstrategy(OpenMPcpuIdx, out DisableOpenMP_becauseIsSlow);
                        OMPbindingStrategy = ilPSP.OMPbindingStrategy.CloneFromMain;
                        OnlinePerformanceMeasurement.Log.OMPbindingStrategy = OMPbindingStrategy.Value;
                    }

                    if (DisableOpenMP_becauseIsSlow) {
                        DisableOpenMP();
                        CPUAffinity.SetCurrentThreadAffinity(DedicatedCPUsForThisRank);
                    } else {
                        if(last_OpenMPcpuIdx.SetEquals(OpenMPcpuIdx) == false) {
                            last_OpenMPcpuIdx = OpenMPcpuIdx;
                            MKLservice.BindOMPthreads(OpenMPcpuIdx, OMPbindingStrategy.Value);
                        }
                    }
                }
            }*/

        }

        //static int[] last_OpenMPcpuIdx; 
        //static OMPbindingStrategy? OMPbindingStrategy;

        static bool DisableOpenMP_becauseIsSlow = false;

        /// <summary>
        /// true, if the code currently runs in multi-threaded; then, further spawning into sub-threads should not occur.
        /// </summary>
        public static bool InParallelSection {
            get;
            private set;
        } = false;

        static bool m_StdoutOnlyOnRank0 = false;

        /// <summary>
        /// if true, the standard - output stream will not be visible on screen on processor with MPI rank 
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

        static MPIEnvironment m_MpiEnv;

        /// <summary>
        /// environment of the world communicator
        /// </summary>
        public static MPIEnvironment MPIEnv {
            get {
                return m_MpiEnv;
            }
        }
    }
}

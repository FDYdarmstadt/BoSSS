using BoSSS.Application.BoSSSpad;
using BoSSS.Solution;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using log4net.Appender;
using log4net.Config;
using log4net.Layout;
using MPI.Wrappers;
using NUnit.Framework;
using NUnitLite;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading;

namespace PublicTestRunner {

    /// <summary>
    /// Redirection of fails in <see cref="System.Diagnostics.Debug"/> to NUnit assertions 
    /// </summary>
    class MyListener : System.Diagnostics.TraceListener {
        public override void Write(string message) => Console.Error.Write(message);

        public override void WriteLine(string message) => Console.Error.WriteLine(message);

        public override void Fail(string message) {
            Console.Error.WriteLine("Assertion Fail: >>>" + message + "<<<");
            Assert.Fail(message);
        }

        public override void Fail(string message, string detailMessage) {
            Console.Error.WriteLine("Assertion Fail: >>>" + message + ": " + ", detailed message is: " + detailMessage );
            Assert.Fail(message + ", details: " + detailMessage);
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public interface ITestTypeProvider {

        /// <summary>
        /// List of serial tests (one MPI core) that should be executed in DEBUG and RELEASE mode.
        /// Referencing any type of the assembly will do.
        /// </summary>
        Type[] FullTest { get; }

        /// <summary>
        /// List of serial tests (one MPI core) that should be executed only in RELEASE mode.
        /// Referencing any type of the assembly will do.
        /// </summary>
        Type[] ReleaseOnlyTests { get; }

        /// <summary>
        /// List of parallel tests that should be executed in DEBUG and RELEASE mode.
        /// Referencing any type of the assembly will do.
        /// </summary>
        (Type type, int NoOfProcs)[] MpiFullTests { get; }

        /// <summary>
        /// List of parallel tests that should be executed only in RELEASE mode.
        /// Referencing any type of the assembly will do.
        /// </summary>
        (Type type, int NoOfProcs)[] MpiReleaseOnlyTests { get; }

        /// <summary>
        /// Root for data searches
        /// </summary>
        DirectoryInfo GetRepositoryBaseDir();

        /// <summary>
        /// If true, the managed assemblies are not copied for every job; there is only a single copy to reduce IO load.
        /// Uses the <see cref="Job.EntryAssemblyRedirection"/> - hack;
        /// </summary>
        bool CopyManagedAssembliesCentrally { get; }

        /// <summary>
        /// If true, the deployment directories are deleted for jobs which finished successfully
        /// </summary>
        bool DeleteSuccessfulTestFiles { get; }

        /// <summary>
        /// Number of tries when a job fails
        /// </summary>
        int RetryCount { get; }
    }

    /// <summary>
    /// 
    /// </summary>
    public class PublicTests : ITestTypeProvider {

        /// <summary>
        /// Note: for better use of parallel resources, try to put the expensive tests in front
        /// </summary>
        public virtual Type[] FullTest => new Type[] {
                        typeof(ilPSP.MultidimensionalArray_Tests),
                        typeof(BoSSS.Application.SipPoisson.SipPoissonMain),
                        typeof(AdvancedSolverTests.TestsMain),
                        typeof(BoSSS.Application.CDG_ProjectionTest.AllUpTest),
                        typeof(BoSSS.Application.Matrix_MPItest.AllUpTest),
                        typeof(BoSSS.Application.XdgPoisson3.XdgPoisson3Main),
                        typeof(CDG_Projection_MPI.ConstrainedDGField_Tests),
                        typeof(BoSSS.Application.ElementTests.ElementTests),
                        typeof(BoSSS.Application.MultigridTest.MultigridMain),
                        typeof(BoSSS.Application.DatabaseTests.DatabaseTestsProgram),
                        typeof(BoSSS.Application.XDGTest.UnitTest),
                        typeof(BoSSS.Application.DerivativeTest.DerivativeTestMain),
                        typeof(CutCellQuadrature.CutCellQuadratureMain),
                        typeof(BoSSS.Application.SpecFEM.AllUpTest),
                        typeof(BoSSS.Application.ipViscosity.TestSolution),
                        //typeof(BoSSS.Application.LevelSetTestBench.LevelSetTestBenchMain),
                        //typeof(BoSSS.Application.AdaptiveMeshRefinementTest.AllUpTest),
                        typeof(BoSSS.Application.ExternalBinding.CodeGen.Test),
                        typeof(BoSSS.Application.ExternalBinding.Initializer),
                        typeof(BoSSS.Application.XNSE_Solver.XNSE),
                        typeof(MPITest.Program),
                        typeof(BoSSS.Application.CutCellQuadratureScaling.AllTests),
                        typeof(BoSSS.Application.TraceDGtest.UnitTests)
                    };

        /// <summary>
        /// Note: for better use of parallel resources, try to put the expensive tests in front
        /// </summary>
        virtual public Type[] ReleaseOnlyTests => new Type[] {
                    typeof(BoSSS.Application.LsTest.SolverWithLevelSetUpdaterTestCenter),
                    typeof(BoSSS.Application.XNSERO_Solver.XNSERO),
                    //typeof(BoSSS.Application.XNSE_Solver.XNSE),
                    typeof(BoSSS.Application.XNSFE_Solver.XNSFE),
                    typeof(BoSSS.Application.XdgTimesteppingTest.XdgTimesteppingMain),
                    typeof(CNS.CNSProgram),
                    typeof(NSE_SIMPLE.SIMPLESolver),
                    typeof(BoSSS.Application.ZwoLsTest.AllUpTest),
                    typeof(QuadratureAndProjectionTest.QuadratueAndProjectionTest),
                    typeof(BoSSS.Application.XdgNastyLevsetLocationTest.AllUpTest),
                    typeof(LTSTests.Program),
                    typeof(BoSSS.Application.TutorialTests.AllUpTest),
                    typeof(BoSSS.Application.XNSEC.XNSEC),
                    typeof(BoSSS.Application.ExternalBinding.CahnHilliardTest),
                    //typeof(BoSSS.Application.XNSE_ViscosityAgglomerationTest.XNSE_ViscosityAgglomerationTestMain),
                    typeof(ALTSTests.Program),
                    typeof(ZwoLevelSetSolver.ZLS),
                    typeof(HFSISolver.HFSI),
                    typeof(HangingNodesTests.HangingNodesTestMain),
                    typeof(BoSSS.Application.CahnHilliard.CahnHilliardMain),
                    typeof(IntersectingLevelSetTest.AllUpTest),
                    typeof(BUIDT.Tests.BUIDTTestProgram),
                    typeof(SAIDT.Tests.SAIDTTestProgram),
                    typeof(XESF.Tests.XESFTestProgram),
                    //typeof(XNSE_ParallelTests.XNSE_ParallelTests),
                    typeof(StokesHelical_Ak.HelicalMain)
                };
            
        

        public virtual (Type type, int NoOfProcs)[] MpiFullTests => new (Type type, int NoOfProcs)[] {
                        (typeof(CDG_Projection_MPI.ConstrainedDGField_Tests), 4),
                        (typeof(CDG_Projection_MPI.ConstrainedDGField_Tests), 2),
                        (typeof(AdvancedSolverTests.TestsMain),4),
                        (typeof(AdvancedSolverTests.MPITests),4),
                        (typeof(MPITest.Program), 4),
                        (typeof(MPITest.Program), 3),
                        (typeof(MPITest.Program), 2),
                        (typeof(BoSSS.Application.SpecFEM.AllUpTest), 4),
                    };

        public virtual (Type type, int NoOfProcs)[] MpiReleaseOnlyTests => new (Type type, int NoOfProcs)[] {
                        (typeof(BoSSS.Application.XNSE_Solver.XNSE_Solver_LargeMPItest), 8),
                        (typeof(BoSSS.Application.XNSE_Solver.XNSE_Solver_MPItest), 3),
                        (typeof(BoSSS.Application.XNSE_Solver.XNSE_Solver_MPItest), 4),
                        (typeof(BoSSS.Application.XdgPoisson3.XdgPoisson3Main), 4),
                        (typeof(MPITest.Program), 4),
                        //(typeof(HangingNodesTests.HangingNodesTestMain), 2), // fk, 29mar22: parallel runs executed directly in `release.yml` to allow serial-parallel comparison
                        //(typeof(HangingNodesTests.HangingNodesTestMain), 4), // fk, 29mar22: parallel runs executed directly in `release.yml` to allow serial-parallel comparison
                        (typeof(BoSSS.Application.SpecFEM.AllUpTest), 4),
                        (typeof(BoSSS.Application.Matrix_MPItest.AllUpTest), 4),
                        (typeof(BoSSS.Application.LoadBalancingTest.LoadBalancingTestMain), 4),
                        (typeof(ALTSTests.Program), 4),
                        (typeof(CNS_MPITests.Tests.LoadBalancing.ShockTubeLoadBalancingTests), 4),
                        (typeof(HilbertTest.HilbertTest), 4),
                    };

        public virtual DirectoryInfo GetRepositoryBaseDir() {
            DirectoryInfo repoRoot;
            try {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                repoRoot = dir?.Parent?.Parent?.Parent?.Parent?.Parent?.Parent;
                if(repoRoot == null)
                    return null;

                DirectoryInfo src = repoRoot.GetDirectories("src").SingleOrDefault();
                DirectoryInfo libs = repoRoot.GetDirectories("libs").SingleOrDefault();
                DirectoryInfo doc = repoRoot.GetDirectories("doc").SingleOrDefault();
                if ( src == null || !src.Exists ) {
                    return null;
                }
                //throw new Exception();
                if ( libs == null || !libs.Exists ) {
                    return null;
                }
                //throw new Exception();
                if ( doc == null || !doc.Exists ) {
                    return null;
                }
                //throw new Exception();

            } catch ( Exception ) {
                // must be null, for the nunit3 sub-program when running from somewhere.
                return null;
                //throw new IOException("Unable to find repository root. 'runjobmanger' must be invoked from its default location within the BoSSS git repository. You might want to use 'runjobmanger-ignore_tests_w_deps'");
            }

            return repoRoot;
        }

        public virtual bool CopyManagedAssembliesCentrally => true;

        public virtual int RetryCount => 3;

        public virtual bool DeleteSuccessfulTestFiles => true;
    }

    /// <summary>
    /// Subroutines of the test runner
    /// </summary>
    public static class PublicTestRunnerMain {

        static ITestTypeProvider TestTypeProvider;

        public static bool discoverRelease = true;

        /// <summary>
        /// Timout for job-manager run
        /// </summary>
        public static double TimeOutSec = 6*60*60; // 6 hours

        /// <summary>
        /// supposed to ignore tests depending on files in the source code repo;
        /// thereby, we can run the test runner from outside the source repository.
        /// </summary>
        static bool ignore_tests_w_deps = false;

        private static List<Job> AllJobs = new();

        /// <summary>
        /// finds all assemblies which potentially contain tests.
        /// </summary>
        static Assembly[] GetAllAssembliesForTests() {
            using ( new FuncTrace() ) {
                var R = new HashSet<Assembly>();

                if ( TestTypeProvider.FullTest != null ) {
                    foreach ( Type t in TestTypeProvider.FullTest ) {
                        //Console.WriteLine("test type: " + t.FullName);
                        Assembly a = t.Assembly;

                        //Console.WriteLine("  assembly: " + a.FullName + " @ " + a.Location);
                        _ = R.Add(a);
                        //Console.WriteLine("  added? " + added);
                    }
                }

                if ( discoverRelease ) {
                    if ( TestTypeProvider.ReleaseOnlyTests != null ) {
                        foreach ( Type t in TestTypeProvider.ReleaseOnlyTests ) {
                            _ = R.Add(t.Assembly);
                        }
                    }
                }

                return R.ToArray();
            }
        }

        //static bool IsReleaseOnlyAssembly

        static (Assembly Asbly, int NoOfProcs)[] GetAllMpiAssemblies() {
            var R = new List<(Assembly Asbly, int NoOfProcs)>();

            if ( TestTypeProvider.MpiFullTests != null ) {
                foreach ( (Type type, int NoOfProcs) t in TestTypeProvider.MpiFullTests ) {
                    //Console.WriteLine("test type: " + t.type.FullName + " (" + t.NoOfProcs + " procs).");
                    //Console.WriteLine("  assembly: " + t.type.Assembly.FullName + " @ " + t.type.Assembly.Location);
                    bool contains = R.Contains(t, (itm1, itm2) => (itm1.NoOfProcs == itm2.NoOfProcs) && itm1.Asbly.Equals(itm2.type.Assembly));
                    if ( !contains ) {
                        R.Add((t.type.Assembly, t.NoOfProcs));
                    }
                    //Console.WriteLine("  added? " + (!contains));
                }
            }

            if ( discoverRelease ) {
                if ( TestTypeProvider.MpiReleaseOnlyTests != null ) {
                    foreach ( (Type type, int NoOfProcs) t in TestTypeProvider.MpiReleaseOnlyTests ) {
                        bool contains = R.Contains(t, (itm1, itm2) => (itm1.NoOfProcs == itm2.NoOfProcs) && itm1.Asbly.Equals(itm2.type.Assembly));
                        if ( !contains ) {
                            R.Add((t.type.Assembly, t.NoOfProcs));
                        }
                    }
                }
            }

            return R.ToArray();
        }

        static string[] LocateFile(string PartialPath) {
            DirectoryInfo repoRoot = TestTypeProvider.GetRepositoryBaseDir();
            if ( repoRoot == null ) {
                return null;
            }

            // if we get here, we probably have access to the repository root directory.
            string[] r = LocateFileRecursive("", repoRoot, PartialPath);
            if ( r == null || r.Length <= 0 ) {
                throw new IOException("unable to find file '" + PartialPath + "'");
            }

            //if (r.Length > 1) {
            //    throw new IOException("The path '" + PartialPath + "' has been found several times: " + r.ToConcatString("", ", ", ";") + " Did you specify the path correctly?");
            //}

            return r;
        }

        static string[] LocateFileRecursive(string RelPath, DirectoryInfo absPath, string SomeFileName) {
            var ret = new List<string>();

            string _SomeFileName = "*" + SomeFileName;

            foreach ( FileInfo f in absPath.GetFiles() ) {
                string RelName = RelPath + f.Name;

                if ( RelName.EndsWith(SomeFileName) ) {
                    ret.Add(f.FullName);
                } else if ( SomeFileName.WildcardMatch(RelName) ) {
                    ret.Add(f.FullName);
                } else if ( _SomeFileName.WildcardMatch(RelName) ) {
                    ret.Add(f.FullName);
                }
            }

            foreach ( DirectoryInfo d in absPath.GetDirectories() ) {
                ret.AddRange(LocateFileRecursive(RelPath + d.Name + "/", d, SomeFileName));
            }

            return ret.ToArray();
        }

        static (int NoOfTests, string[] tests, string[] shortnames, IDictionary<string, string[]> RequiredFiles4Test, int?[] NumThreads) GetTestsInAssembly(Assembly a, string AssemblyFilter) {
            var r = new List<string>(); // full test names
            var n = new List<int?>(); // number of threads for test
            var l = new List<string>(); // short names 
            var d = new Dictionary<string, string[]>(); // pairs of (full test name | files for the test)

            var g = new HashSet<string>(); // global files for all tests

            Type[] ttt = a.GetTypes();
            foreach ( Type t in ttt ) { // loop over types in assembly...
                if ( t.IsClass ) {
                    MethodInfo[] mmm = t.GetMethods();

                    // first, seach for num_threads defined for the entire class...
                    int? num_threads_class = null;
                    if ( t.GetCustomAttribute(typeof(NUnitNumThreads)) != null ) {
                        var nta = t.GetCustomAttribute(typeof(NUnitNumThreads)) as NUnitNumThreads;
                        //if (dc != null)
                        //    throw new NotSupportedException("Palacing a `NUnitFileToCopyHackAttribute` together with `SetUpAttribute` or `OneTimeSetUpAttribute` is not supported (anymore).");

                        num_threads_class = nta.NumThreads;
                    }

                    foreach ( MethodInfo m in mmm ) { // loop over methods in type...
                        if ( t.IsAbstract && !m.IsStatic ) {
                            continue;
                        }

                        if ( !FilterTestMethod(m, AssemblyFilter) ) {
                            continue;
                        }

                        var s = new HashSet<string>();

                        if ( m.GetCustomAttribute(typeof(SetUpAttribute)) != null
                         || m.GetCustomAttribute(typeof(OneTimeSetUpAttribute)) != null ) {
                            //if (dc != null)
                            //    throw new NotSupportedException("Palacing a `NUnitFileToCopyHackAttribute` together with `SetUpAttribute` or `OneTimeSetUpAttribute` is not supported (anymore).");

                            if ( m.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) is NUnitFileToCopyHackAttribute dc ) {
                                foreach ( string someFile in dc.SomeFileNames ) {
                                    g.AddRange(LocateFile(someFile));
                                }
                            }
                        }

                        // search int num_threads is "overridden" for the respective test
                        int? num_threads = num_threads_class;
                        if ( m.GetCustomAttribute(typeof(NUnitNumThreads)) != null ) {
                            var nta = m.GetCustomAttribute(typeof(NUnitNumThreads)) as NUnitNumThreads;
                            //if (dc != null)
                            //    throw new NotSupportedException("Palacing a `NUnitFileToCopyHackAttribute` together with `SetUpAttribute` or `OneTimeSetUpAttribute` is not supported (anymore).");

                            num_threads = nta.NumThreads;
                        }

                        //bool testAdded = false;
                        if ( m.GetCustomAttribute(typeof(TestAttribute)) != null ) {
                            r.Add(t.FullName + "." + m.Name);
                            n.Add(num_threads);
                            l.Add(Path.GetFileNameWithoutExtension(a.ManifestModule.Name) + "#" + m.Name);
                            //Console.WriteLine("Added: " + r.Last());
                            //testAdded = true;
                            //}

                            //if (m.GetCustomAttribute(typeof(TestAttribute)) != null
                            //   || m.GetCustomAttribute(typeof(SetUpAttribute)) != null
                            //   || m.GetCustomAttribute(typeof(OneTimeSetUpAttribute)) != null) {

                            if ( m.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) is NUnitFileToCopyHackAttribute dc ) {
                                //Console.WriteLine("Added: " + r.Last() + " depends on " + dc.SomeFileNames[0]);
                                //if(ignore_tests_w_deps && testAdded) {
                                //    // supposed to ignore tests depending on files in the source code repo
                                //    r.RemoveAt(r.Count - 1);
                                //    l.RemoveAt(l.Count - 1);
                                //    continue; // skip this test
                                //}

                                foreach ( string someFile in dc.SomeFileNames ) {
                                    s.AddRange(LocateFile(someFile));
                                }
                            }

                            //if (System.Environment.GetEnvironmentVariable("BOSSS_RUNTESTFROMBACKUP").IsEmptyOrWhite() == false) {
                            //    s.Add("BOSSS_RUNTESTFROMBACKUP.txt");
                            //}

                            d.Add(r.Last(), s.ToArray());
                        }
                    }
                }
            }

            // add assembly-global files for test, check uniqueness of filenames
            foreach ( string testname in d.Keys ) {
                var s = new List<string>();
                s.AddRange(d[testname]);
                s.AddRange(g);

                List<string> FileNamesOnly = new();
                foreach ( string filePath in s ) {
                    string fileName = Path.GetFileName(filePath);
                    if ( FileNamesOnly.Contains(fileName, (string a, string b) => a.Equals(b, StringComparison.InvariantCultureIgnoreCase)) ) {
                        throw new IOException($"Dependent Filename {fileName} is not unique for test assembly {a}. (full Path {filePath}).");
                    }

                    FileNamesOnly.Add(fileName);
                }

                d[testname] = s.ToArray();
            }

            return n.Count != r.Count ? throw new ApplicationException("internal error") : ((int NoOfTests, string[] tests, string[] shortnames, IDictionary<string, string[]> RequiredFiles4Test, int?[] NumThreads)) (r.Count, r.ToArray(), l.ToArray(), d, n.ToArray());
        }

        /*
        class MyTestFilter : NUnit.Framework.Internal.TestFilter {
            public override TNode AddToXml(TNode parentNode, bool recursive) {
                throw new NotImplementedException();
            }

            public override bool Match(ITest test) {
                throw new NotImplementedException();
            }
        }
        */

        static TextWriterAppender logger_output = null;
        static Stream tracerfile;
        static TextWriter tracertxt;

        static void InitTraceFile(string basename) {

            if ( logger_output != null ) {
                throw new ApplicationException("Already called."); // is seems this object is designed so that it stores at max one session per lifetime
            }

            tracerfile = new FileStream($"trace_{basename}.txt", FileMode.Create, FileAccess.Write, FileShare.Read);
            tracertxt = new StreamWriter(tracerfile);

            TextWriterAppender fa = new TextWriterAppender();
            fa.ImmediateFlush = true;
            //fa.Writer = Console.Out;
            fa.Writer = tracertxt;
            fa.Layout = new PatternLayout("%date %-5level From__%logger: %message%newline");
            fa.ActivateOptions();
            _ = BasicConfigurator.Configure(fa);
            logger_output = fa;

        }

        static void CloseTracing() {
            try {
                logger_output.Close();
                logger_output = null;

                tracertxt.Flush();
                tracertxt.Close();
                tracertxt.Dispose();
                tracertxt = null;

                tracerfile.Flush();
                tracerfile.Close();
                tracerfile.Dispose();
            } catch ( Exception ) {
                //Console.Error.WriteLine(e.GetType().Name + " during closing of tracing: " + e.Message);
            }
        }

        /// <summary>
        /// Decides whether an assembly <paramref name="a"/> matches the filter (wildcard <paramref name="AssemblyFilter"/>) specified by the user
        /// </summary>
        static bool FilterTestAssembly(Assembly a, string AssemblyFilter) {
            if ( AssemblyFilter.IsEmptyOrWhite() ) {
                return true;
            }

            string[] sFilters = AssemblyFilter.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
            foreach ( string filter in sFilters ) {
                string modFilter;
                bool expect = false;
                if ( filter.StartsWith("!") ) {
                    modFilter = filter[1..];
                    expect = false;
                } else {
                    modFilter = filter;
                    expect = true;
                }

                modFilter = modFilter.Split("\\", StringSplitOptions.RemoveEmptyEntries)[0];

                if ( modFilter.WildcardMatch(Path.GetFileNameWithoutExtension(a.Location)) == expect ) {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Decides whether an assembly <paramref name="m"/> matches the filter (wildcard <paramref name="AssemblyFilter"/>) specified by the user
        /// </summary>
        static bool FilterTestMethod(MethodInfo m, string AssemblyFilter) {

            if ( AssemblyFilter.IsEmptyOrWhite() ) {
                return true;
            }

            string FullName = m.DeclaringType.Name + "." + m.Name;
            Assembly a = m.DeclaringType.Assembly;

            string[] sFilters = AssemblyFilter.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
            foreach ( string filter in sFilters ) {
                string modFilter;
                bool expect = false;
                if ( filter.StartsWith("!") ) {
                    modFilter = filter[1..];
                    expect = false;
                } else {
                    modFilter = filter;
                    expect = true;
                }

                string[] AssemblyNMethodFilter = modFilter.Split("\\", StringSplitOptions.RemoveEmptyEntries);
                string aFilter = AssemblyNMethodFilter[0];

                if ( aFilter.WildcardMatch(Path.GetFileNameWithoutExtension(a.Location)) == expect ) {

                    if ( AssemblyNMethodFilter.Length < 2 ) {
                        return true;
                    }

                    string mfilter = AssemblyNMethodFilter[1];

                    if ( mfilter.WildcardMatch(FullName) == expect ) {
                        return true;
                    }
                }
            }

            return false;

        }


        /// <summary>
        /// to avoid IO collisions for concurrent runs of the job manager on the same machine (e.g. DEBUG and RELEASE);
        /// appending of the user name avoids "unauthorized access"-exceptions
        /// </summary>
        static readonly Mutex IOsyncMutex;

        static PublicTestRunnerMain() {
            var rnd = new Random(Guid.NewGuid().GetHashCode());
            string mutex_name = "_BoSSS_test_runner_IOmutex_" + System.Environment.UserName;
            for ( int iRetry = 0; iRetry < 10; iRetry++ ) {
                try {
                    IOsyncMutex = new Mutex(false, mutex_name);
                    break;
                } catch ( Exception ) {
                    Thread.Sleep(rnd.Next(1000));
                    IOsyncMutex = null;
                }
            }

            if ( IOsyncMutex == null ) {
                Console.WriteLine("Unable to create mutex: " + mutex_name);
            }
        }

        /// <summary>
        /// to distinct the internalTestRunner
        /// </summary>
        public static string RunnerPrefix = "Pub";

        public static void CancelAllJobsOnExit(object sender, EventArgs args) {
            var running = AllJobs.Where(j => !(j.Status == JobStatus.FinishedSuccessful || j.Status == JobStatus.FailedOrCanceled)).ToArray();
            if(running.Length < 0) {
                Console.WriteLine("No running jobs.");
                return;
            }

            Console.WriteLine($"Forcefully Cancelling {running.Length} Jobs");
            foreach( Job job in running) {
                Console.WriteLine($"Canceling {job.Name}...");
                job.LatestDeployment.Cancel("Forcefully Cancelled");
            }
            
        }

        public static int JobManagerRun(string AssemblyFilter, int ExecutionQueueNo) {
            AppDomain.CurrentDomain.ProcessExit += CancelAllJobsOnExit;
            Console.CancelKeyPress += CancelAllJobsOnExit;

            // ===================================
            // phase 0: setup
            // ===================================

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MpiSize);
            if ( MpiSize != 1 ) {
                throw new NotSupportedException("runjobmanager subprogram must be executed serially");
            }

            BoSSSshell.ReloadExecutionQueues();

            if ( ExecutionQueueNo >= BoSSSshell.ExecutionQueues.Count ) {
                throw new ApplicationException($"Execution queue #{ExecutionQueueNo} does not exist on this machine/account (see configuration file ~/.BoSSS/etc/BatchProcessorConfig.json).");
            }

            BatchProcessorClient bpc = BoSSSshell.ExecutionQueues[ExecutionQueueNo];
            Console.WriteLine($"Using batch queue {ExecutionQueueNo}: {bpc}");

            FileStream ServerMutex;
            string DateNtime = null;
            try {
                _ = (IOsyncMutex?.WaitOne());

                var rnd = new Random(DateTime.Now.Millisecond + typeof(PublicTestRunnerMain).Assembly.Location.GetHashCode() + Directory.GetCurrentDirectory().GetHashCode());
                Thread.Sleep(rnd.Next(10000)); // sleep for a random amount of time to avoid 
                do {
                    DateNtime = DateTime.Now.ToString("MMMdd_HHmmss");
                    string MutexFileName = Path.Combine(bpc.DeploymentBaseDirectory, RunnerPrefix + DebugOrReleaseSuffix + "_" + DateNtime + ".lock");
                    try {
                        ServerMutex = File.Open(MutexFileName, FileMode.Create, FileAccess.Write, FileShare.None);
                        using var wrt = new StreamWriter(ServerMutex);
                        wrt.WriteLine("Locked by BoSSS test runner at " + DateNtime);
                    } catch ( Exception ee ) {
                        Console.Error.WriteLine($"Unable to get lock on {MutexFileName}: {ee}");
                        Console.Error.WriteLine("wait and try again...");
                        ServerMutex = null;
                        Thread.Sleep(rnd.Next(10000));
                    }
                } while ( ServerMutex == null );
                Console.WriteLine($"Using prefix'{DateNtime}' for all jobs.");
                if ( DateNtime == null ) {
                    throw new ApplicationException("internal error");
                }
            } catch ( Exception e ) {
                Console.Error.WriteLine("UNRECOVERABLE EXCEPTION DURING CREATION OF OUTPUT DIRECTORY");
                Console.Error.WriteLine(e.GetType() + ":  " + e.Message);
                Console.Error.WriteLine(e.StackTrace);
                Console.Error.WriteLine("TERMINATING APPLICATION");
                System.Environment.Exit(-666);
            } finally {
                IOsyncMutex?.ReleaseMutex();
            }

            Tracer.NamespacesToLog = new string[] { "" };
            InitTraceFile("JobManagerRun-" + DateNtime);

            bool I_StartedMinibatch = false;
            if ( bpc is MiniBatchProcessorClient ) {
                I_StartedMinibatch = MiniBatchProcessor.Server.StartIfNotRunning(RunExternal: true);
            }

            int returnCode = 0;
            using ( var tr = new FuncTrace() ) {

                // ===================================
                // phase 1: discover tests
                // ===================================
                //Debugger.Launch();
                BoSSSshell.WorkflowMgm.Init("BoSSStst" + DateNtime, bpc);

                // deployment of assemblies
                string NativeOverride;
                string RelManagedPath;
                if ( TestTypeProvider.CopyManagedAssembliesCentrally ) {
                    string mngdir = RunnerPrefix + DebugOrReleaseSuffix + "_" + DateNtime + "_managed";
                    var ManagedOverride = new DirectoryInfo(Path.Combine(bpc.DeploymentBaseDirectory, mngdir));
                    ManagedOverride.Create();
                    TestTypeProvider.GetType().Assembly.DeployAt(ManagedOverride);

                    RelManagedPath = "../" + mngdir + "/" + Path.GetFileName(TestTypeProvider.GetType().Assembly.Location);

                    if ( !bpc.RuntimeLocation.IsEmptyOrWhite() ) {

                        // 
                        // since we deploy the **managed** assemblies to an alternative location,
                        // we should also deploy the **native binaries**.
                        //

                        string suffix = bpc.RuntimeLocation.IsEmptyOrWhite() ? "amd64" : bpc.RuntimeLocation?.Split(new char[] { '/', '\\' }, StringSplitOptions.RemoveEmptyEntries)?.Last();
                        var _NativeOverride = new DirectoryInfo(Path.Combine(bpc.DeploymentBaseDirectory, RunnerPrefix + DebugOrReleaseSuffix + "_" + DateNtime + "_" + suffix));
                        _NativeOverride.Create();

                        string BosssInstall = BoSSS.Foundation.IO.Utils.GetBoSSSInstallDir();
                        string BosssBinNative = Path.Combine(BosssInstall, "bin", "native", bpc.RuntimeLocation);
                        MetaJobMgrIO.CopyDirectoryRec(BosssBinNative, _NativeOverride.FullName, null);

                        NativeOverride = bpc is SlurmClient slurm ? slurm.DeploymentDirectoryAtRemote(_NativeOverride.FullName) : _NativeOverride.FullName;
                    } else {

                        // not enough info to deploy the binaries;
                        // we just can hope that on the cluster the BOSSS_INSTALL var is correctly set up.

                        NativeOverride = null;
                    }
                } else {
                    NativeOverride = null;
                    RelManagedPath = null; //The job manager will deploy the assemblies
                }

                // collection for all tests:
                var allTests = new List<(Assembly ass, string testname, string shortname, string[] depfiles, int NoOfProcs, int? NumThreads)>();

                // Find all serial tests:
                {
                    Assembly[] assln = GetAllAssembliesForTests();
                    if ( assln != null ) {
                        foreach ( Assembly a in assln ) {
                            if ( FilterTestAssembly(a, AssemblyFilter) ) {
                                (int NoOfTests, string[] tests, string[] shortnames, IDictionary<string, string[]> RequiredFiles4Test, int?[] NumThreads) = GetTestsInAssembly(a, AssemblyFilter);
                                for ( int iTest = 0; iTest < NoOfTests; iTest++ ) {
                                    if ( ignore_tests_w_deps && RequiredFiles4Test[tests[iTest]].Length > 0 ) {
                                        Console.WriteLine($"Skipping all in {a} due to external test dependencies.");
                                        break;
                                    }

                                    allTests.Add((a, tests[iTest], shortnames[iTest], RequiredFiles4Test[tests[iTest]], 1, NumThreads[iTest]));
                                }
                            }
                        }
                    }
                }

                // Find all MPI-parallel tests:
                {
                    (Assembly Asbly, int NoOfProcs)[] ParAssln = GetAllMpiAssemblies();
                    if ( ParAssln != null ) {
                        foreach ( (Assembly Asbly, int NoOfProcs) in ParAssln ) {
                            if ( FilterTestAssembly(Asbly, AssemblyFilter) ) {

                                Assembly a = Asbly;
                                (int NoOfTests, string[] tests, string[] shortnames, IDictionary<string, string[]> RequiredFiles4Test, int?[] NumThreads) = GetTestsInAssembly(a, AssemblyFilter);

                                for ( int iTest = 0; iTest < NoOfTests; iTest++ ) {
                                    if ( ignore_tests_w_deps && RequiredFiles4Test[tests[iTest]].Length > 0 ) {
                                        Console.WriteLine($"Skipping all in {a} due to external test dependencies.");
                                        break;
                                    }

                                    allTests.Add((a, tests[iTest], shortnames[iTest], RequiredFiles4Test[tests[iTest]], NoOfProcs, NumThreads[iTest]));
                                }
                            }
                        }
                    }
                }

                Console.WriteLine($"Found {allTests.Count} individual tests ({DebugOrReleaseSuffix}):");
                int cnt = 0;
                foreach ( (Assembly ass, string testname, string shortname, string[] depfiles, int NoOfProcs, int? NumThreads) in allTests ) {
                    cnt++;
                    Console.WriteLine($"  #{cnt}: {testname}");
                    Console.WriteLine($"     {shortname}");
                    Console.WriteLine($"     {NoOfProcs} MPI processors.");
                }

                {
                    string BOSSS_TEST_RUNNER_GODMODE = Path.Combine(BoSSS.Foundation.IO.Utils.GetBoSSSUserSettingsPath(), "godmode.txt");
                    try {
                        string s = File.ReadAllText(BOSSS_TEST_RUNNER_GODMODE);
                        int godval = int.Parse(s);
                        if ( godval != 0 ) {
                            Console.WriteLine("Detected Godmode-Cheatfile. Setting all tests to success.");
                            return 0;
                        }
                    } catch ( Exception ) { }
                }

                // ===================================
                // phase 2: submit jobs
                // ===================================

                var monitor = new JobDeadlineMonitor(Path.Combine("..", "..", ".."));

                Console.WriteLine($"******* Starting job/test deployment/submission ({DateTime.Now}) *******");

                DateTime start = DateTime.Now;

                cnt = 0;
                var AllOpenJobs = new List<(Job job, string ResFile, string testname)>();
                using ( new BlockTrace("DEPLOYMENT", tr) ) {
                    var checkResFileName = new HashSet<string>();

                    foreach ( (Assembly ass, string testname, string shortname, string[] depfiles, int NoOfProcs, int? NumThreads) in allTests ) {
                        try {
                            cnt++;
                            Console.WriteLine($"Submitting {cnt} of {allTests.Count} ({shortname})...");
                            (Job job, string resultFile, string name) j = SubmitJob(ass, testname, shortname, TestTypeProvider.RetryCount, bpc, depfiles, DateNtime, NoOfProcs, NumThreads, NativeOverride, RelManagedPath, cnt);

                            AllJobs.Add(j.job);

                            // Check if there exists timing information for this job
                            if ( !monitor.JobExists(j.job, DateNtime) ) {
                                Console.ForegroundColor = ConsoleColor.Yellow;
                                Console.WriteLine($"Warning: There is no timing information for {j.job.Name}, please add the information manually to TimeRecords.json");
                                Console.ResetColor();
                            }

                            if ( checkResFileName.Add(j.resultFile) == false ) {
                                throw new IOException($"Result file name {j.resultFile} is used multiple times.");
                            }

                            AllOpenJobs.Add(j);
                            Console.WriteLine($"Successfully submitted {j.job.Name}. \n");
                        } catch ( Exception e ) {
                            Console.Error.WriteLine($"{e.GetType().Name} during job submission: {e.Message}.");
                            returnCode--;
                        }
                    }
                }

                if ( returnCode == 0 ) {
                    Console.WriteLine($"******* All jobs/tests deployed ({DateTime.Now}) *******");
                } else {
                    Console.WriteLine($"******* Deployed finished ({DateTime.Now}) -- SOME DEPLOYMENT(S) FAILED *******");
                }
                // ===================================
                // phase 3: wait until complete...
                // ===================================

                var AllFinishedJobs = new List<(Job job, string ResFile, string testname, JobStatus LastStatus, OnlineProfiling[] profilings)>();

                using ( var ot = new StreamWriter("allout-" + DateNtime + "-" + DebugOrReleaseSuffix + ".txt") ) {

                    (Job job, string ResFile, string testname, JobStatus LastStatus, OnlineProfiling[] profilings)[] UpdateFinishedJobs() {
                        using var trr = new FuncTrace("UpdateFinishedJobs");
                        string CurrentDir = Path.GetDirectoryName(typeof(PublicTestRunnerMain).Assembly.Location);

                        var RecentlyFinished = new List<(Job job, string ResFile, string testname, JobStatus LastStatus, OnlineProfiling[] profilings)>();


                        for ( int iJob = 0; iJob < AllOpenJobs.Count; iJob++ ) {
                            (Job job, string ResFile, string testname) = AllOpenJobs[iJob];
                            JobStatus s = job.Status;


                            if ( monitor.Overdue(job, DateNtime) ) {
                                Console.Error.WriteLine($" ------------------- Overdue: {job.Name}");
                                job.LatestDeployment.Cancel("Job is running unusually long!");
                                continue;
                            }

                            {
                                string resultArg = "--result=";
                                string resArg = job.EnvironmentVars.Values.Single(arg => arg.StartsWith(resultArg));
                                string _resFile = resArg.Replace(resultArg, "");
                                if ( _resFile != ResFile ) {
                                    throw new ApplicationException("internal mismatch in result file name");

                                }
                            }

                            if ( s is JobStatus.FinishedSuccessful ) {
                                monitor.UpdateEntry(job, DateNtime);
                            }

                            if ( s is JobStatus.FailedOrCanceled or JobStatus.Unknown ) {
                                Console.WriteLine($" ------------------- Job Failed ({job.Name}) reason: {s}, Exit Code {job?.LatestDeployment.ExitCode}, {job?.LatestDeployment.DeploymentDirectory}");
                                JobStatus s1 = job.GetStatus(true);
                                if ( s1 != s ) {
                                    Console.WriteLine("changed its mind to: " + s1);
                                    s = s1;
                                }

                                if ( job.SubmitCount < job.RetryCount ) {
                                    Console.WriteLine("Trying once again with failed job...");
                                    job.Reactivate();
                                    continue;
                                }
                            }

                            if ( s is JobStatus.FailedOrCanceled or JobStatus.FinishedSuccessful ) {

                                bool reallyDelete = true;

                                // message:
                                if ( s == JobStatus.FinishedSuccessful ) {
                                    Console.WriteLine(s + ": " + job.Name + " (" + job.LatestDeployment.RunTime + ") // " + testname + " (" + DateTime.Now + ")");
                                } else {
                                    Console.WriteLine(s + ": " + job.Name + " (" + job.LatestDeployment.RunTime + ") // " + testname + " at " + job.LatestDeployment.DeploymentDirectory.FullName + " (" + DateTime.Now + ")");
                                }

                                // copy stdout and stderr to logfile
                                LogResultFile(ot, job, testname, ResFile);

                                // copy xml result file
                                using ( new BlockTrace("copy_nunit_xml_result", trr) ) {
                                    try {
                                        string[] sourceFiles = Directory.GetFiles(job.LatestDeployment.DeploymentDirectory.FullName, "result-*.xml");
                                        sourceFiles = sourceFiles.Cat(Directory.GetFiles(job.LatestDeployment.DeploymentDirectory.FullName, "*.html")); // output of jupyter notebooks

                                        foreach ( string orig in sourceFiles ) {
                                            string n = Path.GetFileName(orig);
                                            string dest = Path.Combine(CurrentDir, n);
                                            File.Copy(orig, dest, true);
                                        }
                                    } catch ( IOException ioe ) {
                                        Console.Error.WriteLine(ioe.GetType().Name + ": " + ioe.Message);
                                        returnCode--;
                                    }
                                }

                                OnlineProfiling[] profilings = null;
                                if ( s == JobStatus.FinishedSuccessful ) {
                                    using ( new BlockTrace("load_profilings", trr) ) {
                                        try {
                                            string[] sourceFiles = Directory.GetFiles(job.LatestDeployment.DeploymentDirectory.FullName, "profiling_bin.*.txt");
                                            //if(sourceFiles.Length <= 0) {
                                            //    Console.WriteLine(" ------------------ No benchmark file");
                                            //}

                                            profilings = new OnlineProfiling[sourceFiles.Length];
                                            int rnk = 0;
                                            foreach ( string orig in sourceFiles ) {
                                                profilings[rnk] = OnlineProfiling.Deserialize(File.ReadAllText(orig));
                                                rnk++;
                                            }

                                            if ( DetectSlowBenchmark(profilings, out _) ) {
                                                reallyDelete = false;
                                            }
                                        } catch ( Exception ioe ) {
                                            Console.Error.WriteLine(ioe.GetType().Name + ": " + ioe.Message);
                                        }
                                    }
                                }

                                // delete deploy directory
                                if ( TestTypeProvider.DeleteSuccessfulTestFiles ) {
                                    using ( new BlockTrace("delete_deploy_dir", trr) ) {
                                        if ( s == JobStatus.FinishedSuccessful && reallyDelete ) {
                                            try {
                                                Directory.Delete(job.LatestDeployment.DeploymentDirectory.FullName, true);
                                            } catch ( Exception e ) {
                                                Console.Error.WriteLine($"{e.GetType().Name}: {e.Message}");
                                            }
                                        }
                                    }
                                }

                                // move job to 'finished' list
                                (Job job, string ResFile, string testname, JobStatus s, OnlineProfiling[] profilings) X = (job, ResFile, testname, s, profilings);
                                AllFinishedJobs.Add(X);
                                RecentlyFinished.Add(X);

                                AllOpenJobs.RemoveAt(iJob);
                                iJob--;
                            }
                        }

                        return RecentlyFinished.ToArray();
                    }

                    _ = UpdateFinishedJobs();
                    double RestTime = Math.Max(1, TimeOutSec - (DateTime.Now - start).TotalSeconds);
                    while ( RestTime > 1.0 && AllOpenJobs.Count > 0 ) {
                        using ( new BlockTrace("Sleeping", tr) ) {
                            Thread.Sleep(2 * 60 * 1000); // sleep for 2 minutes
                        }

                        (Job job, string ResFile, string testname, JobStatus LastStatus, OnlineProfiling[] profilings)[] ll = UpdateFinishedJobs();
                        RestTime = Math.Max(1, TimeOutSec - (DateTime.Now - start).TotalSeconds);
                        Console.WriteLine("Remaining minutes until timeout: " + Math.Round(RestTime / 60.0));
                        Console.Write("     Waiting for: ");
                        int i = 0;
                        foreach ( (Job job, string ResFile, string testname) in AllOpenJobs ) {
                            Console.Write(job.Name);
                            Console.Write(" ");
                            i++;
                            if ( i >= 8 ) {
                                break;
                            }
                        }

                        if ( i < AllOpenJobs.Count ) {
                            Console.WriteLine($"... (and {AllOpenJobs.Count - i} more.)");
                        } else {
                            Console.WriteLine();
                        }
                    }

                    // Stop all remaining Jobs
                    AllOpenJobs.ForEach(a => a.job.LatestDeployment.Cancel("Timeout"));
                    _ = UpdateFinishedJobs();
                }

                // ===================================
                // phase 4: summary
                // ===================================
                monitor.Save();
                GitlabResultsTable.WriteResultsTable(AllFinishedJobs.Select(t => t.job).ToList(), "TestRunner");
                if ( AllOpenJobs.Count == 0 ) {
                    Console.WriteLine("----------------------------------------------------------------------------------------------------");
                    Console.WriteLine($"All jobs/tests finished ({DateTime.Now}) - Summary:");
                    Console.WriteLine("----------------------------------------------------------------------------------------------------");
                } else {
                    Console.WriteLine("----------------------------------------------------------------------------------------------------");
                    Console.WriteLine($"Reached Time-out of {TimeOutSec / 60} Minutes; {AllOpenJobs.Count} jobs/tests still running; ({DateTime.Now}) - Summary:");
                    Console.WriteLine("----------------------------------------------------------------------------------------------------");

                    returnCode -= 1;
                }



                Console.WriteLine("--- Test Results -----------------------------------------------------------------------------------");

                int SuccessfulFinishedCount = 0;
                foreach ( (Job job, string ResFile, string testname, JobStatus LastStatus, OnlineProfiling[] profilings) in AllFinishedJobs.Where(ttt => ttt.LastStatus == JobStatus.FinishedSuccessful) ) { // all success 
                    Console.WriteLine($"{job.Status}: {job.Name} ({job.LatestDeployment.RunTime}) // {testname}");
                    SuccessfulFinishedCount++;
                }

                int OtherStatCount = 0;
                foreach ( (Job job, string ResFile, string testname, JobStatus LastStatus, OnlineProfiling[] profilings) in AllFinishedJobs.Where(ttt => ttt.LastStatus != JobStatus.FinishedSuccessful) ) { // all failed
                    Console.WriteLine($"{job.Status}: {job.Name} ({job.LatestDeployment.RunTime}) // {testname} at {job.LatestDeployment.DeploymentDirectory.FullName}");
                    returnCode -= 1;
                    OtherStatCount++;
                }

                foreach ( (Job job, string ResFile, string testname) in AllOpenJobs ) { // all still running
                    Console.WriteLine($"{job.Status}: {job.Name} ({job.LatestDeployment.RunTime}) // {testname} at {job.Status}");
                    returnCode -= 1;
                    OtherStatCount++;
                }

                // summary table
                var totalTime = new TimeSpan();
                try {
                    var resTable = new Dictionary<string, IList<object>>();
                    const string jn = "JobName";
                    const string tn = "TestName";
                    const string st = "Status";
                    const string rt = "RunTime";
                    const string nd = "Ndeploy";
                    const string os = "OMPstrat";

                    resTable.Add(jn, new List<object>());
                    resTable.Add(tn, new List<object>());
                    resTable.Add(st, new List<object>());
                    resTable.Add(rt, new List<object>());
                    resTable.Add(nd, new List<object>());
                    resTable.Add(os, new List<object>());

                    foreach ( (Job job, string ResFile, string testname, JobStatus LastStatus, OnlineProfiling[] profilings) in AllFinishedJobs ) {
                        resTable[jn].Add(job.Name);
                        resTable[tn].Add(testname);
                        resTable[st].Add(job.Status);
                        resTable[rt].Add(job.LatestDeployment.RunTime.TotalSeconds);
                        resTable[nd].Add(job.AllDeployments.Count);
                        totalTime += job.LatestDeployment.RunTime;

                        string osString;
                        try {
                            osString = profilings?.Select(prf => prf?.OnlinePerformanceLog?.OMPbindingStrategy).ToConcatString("", "-", "");
                            osString ??= "NIX";
                        } catch ( Exception e ) {
                            Console.WriteLine(e);
                            osString = "NIX";
                        }

                        resTable[os].Add(osString);
                    }

                    Console.WriteLine($"Saving Results table to: {Path.GetFullPath("ResultTable.csv")}");
                    resTable.SaveToCSVFile("ResultTable.csv");
                    //File.Copy("ResultTable.csv", Path.Combine("C:\\tmp", "ResultTable-" + DateTime.Now.ToString("MMMdd_HHmmss") + ".csv"));
                } catch ( Exception e ) {
                    Console.WriteLine($"{e.Message}, {e.StackTrace}");
                }

                // very final message:
                if ( SuccessfulFinishedCount == (AllOpenJobs.Count + AllFinishedJobs.Count) ) {
                    Console.WriteLine($"All tests/jobs finished successfully. Sum of runtime of all jobs: {totalTime}");
                    Console.WriteLine("SUCCESS.");

                    if ( returnCode != 0 ) {
                        Console.Error.WriteLine("Ignoring some other error occurred (maybe IO error after successful test run) -- check output log;");
                        //Console.Error.WriteLine("FAILURE.");
                    } else {
                    }

                    returnCode = 0;

                } else {
                    Console.Error.WriteLine($"Only {SuccessfulFinishedCount} tests/jobs finished successfully -- {OtherStatCount} have other states. Sum of runtime of all jobs: {totalTime}");
                    Console.Error.WriteLine("FAILURE.");
                }

                Console.WriteLine($"{DateTime.Now}");
            }

            if ( I_StartedMinibatch ) {
                MiniBatchProcessor.Server.SendTerminationSignal(WaitForOtherJobstoFinish: false);
            }

            CloseTracing();

            return returnCode;
        }

        private static void CurrentDomain_ProcessExit(object sender, EventArgs e) => throw new NotImplementedException();

        /// <summary>
        /// Checks all profiling's in all MPI ranks for something fishy.
        /// </summary>
        /// <param name="profilings">
        /// profiling for each MPI rank
        /// - index: enumeration of MPI ranks
        /// </param>
        /// <param name="worstBenchmarks">
        /// On exit, the worst (i.e. slowest) result of the respective benchmark, over all taken measurements in all MPI processes
        /// - keys: benchmark name in <see cref="ilPSP.OnlinePerformanceLog.BenchResults"/>
        /// - values: minimum 
        /// </param>
        /// <returns>
        /// - true: if the relative speed of any benchmark is below 0.1
        /// - false: everything is fine!
        /// </returns>
        private static bool DetectSlowBenchmark(OnlineProfiling[] profilings, out Dictionary<string, double> worstBenchmarks) {
            bool veryBadResultDetected = false;

            worstBenchmarks = new Dictionary<string, double>();
            foreach ( OnlineProfiling profiling in profilings ) {
                Dictionary<string, List<double>> BenchResults = profiling?.OnlinePerformanceLog?.BenchResults;

                if ( BenchResults != null ) {
                    foreach ( KeyValuePair<string, List<double>> kv in BenchResults ) {
                        IEnumerable<double> validMeasurements = kv.Value.Where(val => val > 0.0); // non-positive values indicate skipped benchmarks
                        if ( validMeasurements.Count() > 0 ) {
                            double worstRes = validMeasurements.Min();
                            if ( worstBenchmarks.ContainsKey(kv.Key) ) {
                                worstRes = Math.Min(worstRes, worstBenchmarks[kv.Key]);
                            }

                            if ( worstRes is > 0 and < 0.1 ) {
                                veryBadResultDetected = true;
                            }

                            worstBenchmarks[kv.Key] = worstRes;
                        }
                    }
                }
            }

            return veryBadResultDetected;
        }

        /// <summary>
        /// Writes a one-line summary of all profiling's in all MPI ranks for something fishy.
        /// </summary>
        /// <param name="profilings">
        /// profiling for each MPI rank
        /// - index: enumeration of MPI ranks
        /// </param>        /// <returns></returns>
        private static string SummarizeProfilings(OnlineProfiling[] profilings) {
            string BenchmarkSummary;
            using ( var stw = new StringWriter() ) {

                bool veryBadResultDetected = DetectSlowBenchmark(profilings, out Dictionary<string, double> worstBenchmarks);

                stw.Write(profilings?.Select(prof => prof?.OnlinePerformanceLog?.OMPbindingStrategy).ToConcatString("[", " | ", " ]"));

                if ( veryBadResultDetected ) {
                    stw.Write("!!! SLOW BENCHMARK RESULT !!! ");
                }

                int L = worstBenchmarks.Count;
                int l = 0;
                foreach ( KeyValuePair<string, double> worstRes in worstBenchmarks ) {
                    l++;
                    stw.Write(worstRes.Key + ": " + worstRes.Value.ToString("g4"));
                    if ( worstRes.Value < 0.1 ) {
                        stw.Write(" !!!");
                    }

                    if ( l < L ) {
                        stw.Write(", ");
                    } else {
                        stw.Write(";");
                    }
                }

                if ( veryBadResultDetected ) {
                    stw.Write(" ");
                    stw.Write(profilings?.Select(p => p?.OnlinePerformanceLog).ToConcatString("[", " | ", " ]"));
                }

                BenchmarkSummary = stw.ToString();
            }

            return BenchmarkSummary;
        }

        private static void LogResultFile(StreamWriter ot, Job j, string testname, string ResFile) {
            using ( new FuncTrace("LogResultFile") ) {
                int sz = j.NumberOfMPIProcs;
                ot.WriteLine("##fdhgjegf763748trfhe8hurdsinf598ugf498jvhsn*hbbvc#####!################");
                ot.WriteLine("########################################################################");
                ot.WriteLine("########################################################################");
                ot.WriteLine("####  " + testname);
                ot.WriteLine("########################################################################");
                ot.WriteLine("########################################################################");
                ot.WriteLine("########################################################################");
                ot.WriteLine("#### Deploy directory: " + j.LatestDeployment.DeploymentDirectory.FullName);
                ot.WriteLine("#### Full test name:   " + j.Name);
                ot.WriteLine("#### Number of procs:  " + j.NumberOfMPIProcs);
                ot.WriteLine("#### Status:           " + j.Status);
                ot.WriteLine("#### Job ID:           " + j.LatestDeployment.BatchProcessorIdentifierToken);
                //                                   +   +   +   +
                if ( j.NumberOfMPIProcs <= 1 ) {
                    ot.WriteLine("#### Result File:      " + ResFile);
                    ot.WriteLine("####    exists?        " + File.Exists(Path.Combine(j.LatestDeployment.DeploymentDirectory.FullName, ResFile)));
                } else {
                    for ( int i = 0; i < sz; i++ ) {
                        if ( i == 0 ) {
                            ot.WriteLine("#### Result Files:     " + MpiResFileNameMod(i, sz, ResFile));
                        } else {
                            ot.WriteLine("####                   " + MpiResFileNameMod(i, sz, ResFile));
                        }
                    }

                    for ( int i = 0; i < sz; i++ ) {
                        if ( i == 0 ) {
                            ot.WriteLine("####    exists?        " + File.Exists(Path.Combine(j.LatestDeployment.DeploymentDirectory.FullName, MpiResFileNameMod(i, sz, ResFile))));
                        } else {
                            ot.WriteLine("####                   " + File.Exists(Path.Combine(j.LatestDeployment.DeploymentDirectory.FullName, MpiResFileNameMod(i, sz, ResFile))));
                        }
                    }
                }

                ot.WriteLine("########################################################################");
                ot.WriteLine("########################################################################");
                ot.WriteLine("########################################################################");

                ot.WriteLine("[[[Stdout: ");
                string stdout = null;
                try {
                    stdout = j.Stdout;
                    ot.WriteLine(j.Stdout);
                } catch ( Exception e ) {
                    ot.WriteLine($"{e.GetType().Name} during reading of stdout stream: {e.Message}");
                    ot.WriteLine(e.StackTrace);
                    stdout = "";
                }

                ot.WriteLine("]]]");

                stdout ??= "";

                using ( var str = new StringReader(stdout) ) {
                    string magic = "arg #3 override from environment variable 'BOSSS_ARG_3': --result=";
                    for ( string line = str.ReadLine(); line != null; line = str.ReadLine() ) {
                        if ( line.StartsWith(magic) ) {
                            string file = line.Replace(magic, "");
                            if ( file != ResFile ) {
                                throw new ArgumentException("Internal result file mismatch: " + file + " vs. " + ResFile + " on job " + j.LatestDeployment.BatchProcessorIdentifierToken);
                            }

                            break;
                        }
                    }
                }

                try {
                    string stderr = j.Stderr;

                    if ( stderr.IsEmptyOrWhite() ) {
                        ot.WriteLine("[[[Empty Error Stream: this is good!]]]");
                    } else {
                        ot.WriteLine("[[[Stderr:");
                        ot.WriteLine(stderr);
                        ot.WriteLine("]]]");
                    }
                } catch ( Exception e ) {

                    ot.WriteLine($"{e.GetType().Name} during reading of stderr stream: {e.Message}");
                    ot.WriteLine(e.StackTrace);
                }

                ot.WriteLine();
                ot.WriteLine();
                ot.WriteLine();
                ot.Flush();
            }
        }

        static string DebugOrReleaseSuffix {
            get {
                string dor;
#if DEBUG
                dor = "DBG";
#else
                dor = "REL";
#endif
                return dor;
            }
        }

        public static (Job j, string resultFile, string name) SubmitJob(
            Assembly a,
            string TestName, string Shortname,
            int RetryCount,
            BatchProcessorClient bpc,
            string[] AdditionalFiles,
            string prefix,
            int NoOfMpiProcs,
            int? NumThreads,
            string nativeOverride,
            string TestRunnerRelPath,
            int cnt) {
            using ( new FuncTrace() ) {

                // define unique name (not to long) for the job
                string dor = DebugOrReleaseSuffix;
                string final_jName;
                {
                    string jName = NoOfMpiProcs <= 1 ? $"{prefix}-{dor}-{Shortname}" : $"{prefix}-{dor}-p{NoOfMpiProcs}-{Shortname}";

                    if ( jName.Length > 127 ) {
                        // Name length limit set by MS HPC cluster
                        jName = jName[..127];
                    }

                    int counter = 2;
                    final_jName = jName;
                    while ( BoSSSshell.WorkflowMgm.AllJobs.ContainsKey(final_jName) ) {
                        string suffix = "_" + counter;
                        counter++;
                        final_jName = jName.Length + suffix.Length > 127 ? jName[..(127 - suffix.Length)] : jName;
                        final_jName += suffix;
                    }
                }

                // create job
                var j = new Job(final_jName, TestTypeProvider.GetType()) {
                    SessionReqForSuccess = false,
                    RetryCount = RetryCount
                };
                string resultFile = $"result-{dor}-{cnt}.xml";
                j.MySetCommandLineArguments("nunit3", Path.GetFileNameWithoutExtension(a.Location), $"--test={TestName}", $"--result={resultFile}");
                foreach ( string f in AdditionalFiles ) {
                    j.AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(File.ReadAllBytes(f), Path.GetFileName(f)));
                }

                if ( BOSSS_RUNTESTFROMBACKUP_ENVVAR ) {
                    j.AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(File.ReadAllBytes("BOSSS_RUNTESTFROMBACKUP.txt"), "BOSSS_RUNTESTFROMBACKUP.txt"));
                }

                if ( BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER_ENVVAR ) {
                    j.AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(File.ReadAllBytes("BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt"), "BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt"));
                }

                if ( nativeOverride != null ) {
                    j.EnvironmentVars.Add(BoSSS.Foundation.IO.Utils.BOSSS_NATIVE_OVERRIDE, nativeOverride);
                }

                j.NumberOfMPIProcs = NoOfMpiProcs;
                if ( NumThreads != null ) {
                    j.NumberOfThreads = NumThreads.Value;
                } else {
                    j.NumberOfThreads = 1; // reduced no of threads to improve cluster throughput.
                }

                if ( TestRunnerRelPath != null ) {
                    j.EntryAssemblyRedirection = TestRunnerRelPath;
                }

                j.Activate(bpc);
                return (j, resultFile, TestName);
            }
        }

        public static void DeleteResultFiles() {
            using var tr = new FuncTrace("DeleteResultFiles");
            string CurrentDir = Path.GetDirectoryName(typeof(PublicTestRunnerMain).Assembly.Location);
            string[] FilesToDel = Directory.GetFiles(CurrentDir, "result-*.xml");

            foreach ( string f in FilesToDel ) {

                File.Delete(f);

            }
        }

        /// <summary>
        /// Copies additional files required for some test;
        /// these files are identified via the <see cref="NUnitFileToCopyHackAttribute"/>.
        /// </summary>
        /// <remarks>
        /// - This method only performs some file copy if the test runner is executed from within the source repository;
        /// - otherwise (e.g. if tests are deployed to a HPC cluster) it does not copy anything and we hope that all required files are in place.
        /// </remarks>
        static void MegaMurxPlusPlus(Assembly a, string filter) {
            using var tr = new FuncTrace();
            (int NoOfTests, string[] tests, string[] shortnames, IDictionary<string, string[]> RequiredFiles4Test, int?[] NumThreads) = GetTestsInAssembly(a, filter);

            string dir = Directory.GetCurrentDirectory();
            tr.Info("Current dir: " + dir);

            foreach ( string t in tests ) {
                foreach ( string fOrigin in RequiredFiles4Test[t] ) {
                    tr.Info("Origin file: " + fOrigin);
                    if ( File.Exists(fOrigin) ) {
                        string fDest = Path.Combine(dir, Path.GetFileName(fOrigin));
                        tr.Info("Destination file: " + fDest);
                        File.Copy(fOrigin, fDest, true);
                    } else {
                        tr.Info($"Origin file {fOrigin} NOT FOUND!");
                    }
                }
            }
        }


        const string TelemetryFolder = "c:\\telemetry";

        /// <summary>
        /// Runs all tests serially
        /// </summary>
        static int RunNunit3Tests(string AssemblyFilter, string[] args) {

            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MpiRank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MpiSize);
            ilPSP.Tracing.Tracer.NamespacesToLog = new string[] { "" };
            InitTraceFile($"Nunit3.{DateTime.Now:MMMdd_HHmmss}.{MpiRank}of{MpiSize}");

            string CCP_AFFINITY = System.Environment.GetEnvironmentVariable("CCP_AFFINITY");
            string TelemetryFile = Guid.NewGuid().ToString() + $".{MpiRank}of{MpiSize}";
            var start = DateTime.Now;
            if(CCP_AFFINITY.IsNonEmpty()) {
                if(!Directory.Exists(TelemetryFolder)) {
                    Thread.Sleep(Math.Abs(TelemetryFile.GetHashCode())%10000);
                    Directory.CreateDirectory(TelemetryFolder);
                }

                string TelemetryPath = Path.Combine(TelemetryFolder, TelemetryFile);
                using(var tele = new StreamWriter(TelemetryPath)) {

                    tele.WriteLine(start.ToFileTimeUtc());
                    tele.WriteLine(CPUAffinityWindows.TotalNumberOfCPUs + " " + CPUAffinityWindows.NumberOfCPUsPerGroup);
                    tele.WriteLine(CCP_AFFINITY);
                    tele.WriteLine(CPUAffinityWindows.GetCurrentThreadAffinity().ToConcatString("", ", ", ""));
                    tele.WriteLine(CPUAffinityWindows.Decode_CCP_AFFINITY().ToConcatString("", ", ", ""));
                    tele.WriteLine(start);

                }
            }


            Console.WriteLine($"Running an NUnit test on {MpiSize} MPI processes ...");
            Tracer.NamespacesToLog_EverythingOverrideTestRunner = true;

            {
                string BOSSS_TEST_RUNNER_GODMODE = Path.Combine(BoSSS.Foundation.IO.Utils.GetBoSSSUserSettingsPath(), "godmode.txt");
                try {
                    string s = File.ReadAllText(BOSSS_TEST_RUNNER_GODMODE);
                    int godval = int.Parse(s);
                    if ( godval != 0 ) {
                        Console.WriteLine("Detected Godmode-Cheatfile. Setting all tests to success.");
                        return 0;
                    }
                } catch ( Exception ) { }
            }

            using var ftr = new FuncTrace();
            Assembly[] assln = GetAllAssembliesForTests();

            if ( MpiSize != 1 ) {
                // this seems some parallel run
                // we have to fix the result argument

                if ( args.Where(a => a.StartsWith("--result=")).Count() != 1 ) {
                    throw new ArgumentException("MPI-parallel NUnit runs require the '--result' argument.");
                }

                int i = args.IndexWhere(a => a.StartsWith("--result="));
                string arg_i = args[i];
                string resFileName = arg_i.Replace("--result=", "");
                args[i] = "--result=" + MpiResFileNameMod(MpiRank, MpiSize, resFileName);

                (Assembly Asbly, int NoOfProcs)[] parAssis = GetAllMpiAssemblies();
                foreach ( (Assembly Asbly, int NoOfProcs) in parAssis ) {
                    Assembly a = Asbly;
                    if ( !assln.Contains(a) ) {
                        a.AddToArray(ref assln);
                    }
                }

                int ii = 0;
                foreach ( Assembly a in assln ) {
                    ftr.Info("Assembly #" + ii + ": " + a.ToString());
                    ii++;
                }
            }

            int count = 0;
            bool ret = false;
            foreach ( Assembly a in assln ) {
                if ( !FilterTestAssembly(a, AssemblyFilter) ) {
                    continue;
                }

                if ( GetTestsInAssembly(a, AssemblyFilter).NoOfTests <= 0 ) {
                    Console.WriteLine("Matching Assembly search string, but none of the methods match. (wrong wildcard?)");
                    continue;
                }

                Console.WriteLine("Matching assembly: " + a.Location);
                ftr.Info("found Assembly #" + count + ": " + a.Location);
                count++;

                if ( MpiRank == 0 ) {
                    MegaMurxPlusPlus(a, AssemblyFilter);
                }

                Console.WriteLine("Waiting for all processors to catch up BEFORE starting test(s)...");
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                Console.WriteLine("All Here.");

                int r;
                using ( var bt = new BlockTrace("RUNNING_TEST", ftr) ) {
                    var tr = new TextRunner(a);
                    r = tr.Execute(args);
                    //var summary = tr.Summary;
                    //File.WriteAllText("a.txt", summary.ToString());
                    //ftr.Info("Test summary: " + summary.ToString());
                }

                using ( var bt = new BlockTrace("StdOut/StdErr reset", ftr) ) {
                    Console.SetOut(new StreamWriter(Console.OpenStandardOutput()));
                    Console.SetError(new StreamWriter(Console.OpenStandardError()));

                    bt.Info("Waiting for all processors to catch up AFTER running test(s)...");
                    //csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                    MPICollectiveWatchDog.WatchAtRelease(csMPI.Raw._COMM.WORLD, 10 * 60, true); // 10 Minute timeout for all processors to arrive.
                    bt.Info("All Here.");

                    //var ar = new AutoRun(a);
                    //int r = ar.Execute(args);

                    int[] all_rS = r.MPIAllGather();
                    for ( int rnk = 0; rnk < all_rS.Length; rnk++ ) {
                        bt.Info($"Rank {rnk}: NUnit returned code " + r);
                    }

                    {
                        string currentDirectory = Directory.GetCurrentDirectory();
                        string[] pltFiles = Directory.GetFiles(currentDirectory, "*.plt");

                        long totalSizeBytes = 0;

                        foreach ( string pltFile in pltFiles ) {
                            var fileInfo = new FileInfo(pltFile);
                            totalSizeBytes += fileInfo.Length;
                        }

                        double totalSizeGigabytes = totalSizeBytes / (1024.0 * 1024 * 1024);

                        if ( totalSizeGigabytes > 2.0 ) {
                            bt.Error("Test produced more than 2 Gigabyte of plt-files -- please check!");
                            //throw new IOException("Test produced more than 2 Gigabyte of plt-files -- please check!");
                        }
                    }
                }

                ftr.Info($"failstate before most recent code {ret} (false means OK, r = {r})");
                ret |= (r != 0);
                ftr.Info($"failstate after most recent code {ret} (false means OK, r = {r})");
            }

            {
                ftr.Info("Found  " + count + " assemblies in total");

                if ( count <= 0 ) {
                    Console.WriteLine("Found no assembly matching: " + AssemblyFilter + " (hint: don't provide a filename extension, e.g. '.dll' or '.exe'; assembly names are compared without using an extension, e.g. 'XNSE_Solver', not 'XNSE_Solver.exe' or 'XNSE_Solver.dll'.)");
                    return -1;
                }
            }

            Console.WriteLine();
            ftr.Info($"failstate all tests: {ret} (false means OK)");
            return ret ? -1 : 0;
        }

        private static string MpiResFileNameMod(int MpiRank, int MpiSize, string resFileName) {
            if ( MpiRank > MpiSize ) {
                throw new ArgumentException();
            }

            if ( MpiSize == 1 ) {
                return resFileName;
            }

            resFileName = Path.GetFileNameWithoutExtension(resFileName);

            string newResFileName = resFileName + "." + MpiRank + "of" + MpiSize + ".xml";
            return newResFileName;
        }

        static void PrintMainUsage() {
            //                 0         1         2         3         4         5         6         7         8
            //                 01234567890123456789012345678901234567890123456789012345678901234567890123456789
            Console.WriteLine("");
            Console.WriteLine("Usage is: PublicTestRunner.exe SUBPROGRAM FILTER [additional arguments]");
            Console.WriteLine("  where SUBPROGRAM is one of:");
            Console.WriteLine("         nunit3        : execute all tests in FILTER within this process.");
            Console.WriteLine("                         additional arguments are passed to NUnit.");
            Console.WriteLine("         runjobmanager         : submit tests to the job manger.");
            Console.WriteLine("         runjobmanager-release : submit ALL tests to the job manger.");
            Console.WriteLine("         runjobmanager-debug   : submit only DEBUG tests to the job manger.");
            Console.WriteLine("         runjobmanager-ignore_tests_w_deps : ignore all tests which have the");
            Console.WriteLine("                                 'NUnitFileToCopyHack' attribute; these depend");
            Console.WriteLine("                                 on files within the source repo; skipping those");
            Console.WriteLine("                                 allows to run without source code in place.");
            Console.WriteLine("         list          : list discovered tests in FILTER (no execution).");
            Console.WriteLine("         help          : prints this message.");
            Console.WriteLine("         yaml          : write 'jobs.yml' file for Gitlab.");
            Console.WriteLine("  and FILTER selects some assembly, i.e. DerivativeTests.exe; it can be ");
            Console.WriteLine("  a wildcard, i.e. use * for submitting all tests.");

        }

        public static bool BOSSS_RUNTESTFROMBACKUP_ENVVAR = false;
        public static bool BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER_ENVVAR = false;

        /// <summary>
        /// the real main-function
        /// </summary>
        /// <param name="args"></param>
        /// <param name="ttp">
        /// A hook to find tests within the entire heap of assemblies.
        /// </param>
        /// <returns></returns>
        public static int _Main(string[] args, ITestTypeProvider ttp) {

            if ( args.Length > 5 ) {
                Console.WriteLine($"Warning: got {args.Length} arguments -- are you using this right?");

                Console.WriteLine("No of args: " + args.Length);
                //int i = 0;
                foreach ( string arg in args ) {
                    Console.WriteLine("    " + arg);
                }

                Console.WriteLine($"Remember: if you are using a bash shell, a wildcard * might get interpreted -- use quoted wildcard '*' !");
                return -args.Length;
            }

            TestTypeProvider = ttp;
            Console.WriteLine("BoSSS NUnit test runner.");

            args = BoSSS.Solution.Application.ArgsFromEnvironmentVars(args);

            System.Diagnostics.Trace.Listeners.Clear();
            System.Diagnostics.Trace.Listeners.Add(new MyListener());

            if ( args.Length < 1 ) {
                Console.WriteLine("Insufficient number of arguments.");
                PrintMainUsage();
                return -7777;
            }

            if ( System.Environment.GetEnvironmentVariable("BOSSS_RUNTESTFROMBACKUP").IsNonEmpty() ) {
                BOSSS_RUNTESTFROMBACKUP_ENVVAR = true;
                File.WriteAllText("BOSSS_RUNTESTFROMBACKUP.txt", "Hello, Suckers!");
                Console.WriteLine("trying to forward the BOSSS_RUNTESTFROMBACKUP hack via additional deployment files...");
            }

            if ( System.Environment.GetEnvironmentVariable("BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER").IsNonEmpty() ) {
                BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER_ENVVAR = true;
                File.WriteAllText("BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER.txt", "Hello, Suckers!");
                Console.WriteLine("trying to forward the BOSSS_DELTETE_OLD_DEPLOYMENTS_DATABASES_MASTER hack via additional deployment files...");
            }

            BoSSS.Solution.Application.InitMPI();

            //Console.WriteLine("De-activation of OpenMP parallelization in external libraries; (Multiple processes using OpenMP at the same tine in Windows seem to cause deadlocks with our current MKL version.) ");
            //ilPSP.Environment.DisableOpenMP();

            int ret = -1;
            switch ( args[0] ) {
                case "nunit3":

                    if ( args.Length < 2 ) {
                        Console.WriteLine("Insufficient number of arguments.");
                        PrintMainUsage();
                        return -7777;
                    }
                    ret = RunNunit3Tests(args[1], args.Skip(2).ToArray());
                    break;


                case "list":
                    if ( args.Length < 2 ) {
                        Console.WriteLine("Insufficient number of arguments.");
                        PrintMainUsage();
                        return -7777;
                    } {
                        string myFilter = args[1];
                        int total = 0;
                        Assembly[] assln = GetAllAssembliesForTests();
                        if ( assln == null || assln.Length == 0 ) {
                            Console.WriteLine("No assemblies available from the test provider.");
                            ret = 0;
                            break;
                        }

                        foreach ( var a in assln ) {
                            if ( !FilterTestAssembly(a, myFilter) )
                                continue;
                            var tinfo = GetTestsInAssembly(a, myFilter);
                            if ( tinfo.NoOfTests <= 0 )
                                continue;
                            Console.WriteLine($"Assembly: {Path.GetFileName(a.Location)} ({a.Location})");
                            for ( int i = 0; i < tinfo.NoOfTests; i++ ) {
                                Console.WriteLine($"  {i + 1}  : {tinfo.tests[i]}");
                            }
                            total = tinfo.NoOfTests;
                        }
                        Console.WriteLine($"Found {total} tests matching '{myFilter}'.");
                        ret = 0;
                    }
                    break;

                case "runjobmanager":
                case "runjobmanager-debug":
                case "runjobmanager-release":
                case "runjobmanager-ignore_tests_w_deps":
                    DeleteResultFiles();

                    bool runRelease;
#if DEBUG
                    runRelease = false;
#else
                runRelease = true;
#endif

                    if ( args[0].EndsWith("release") ) {
                        runRelease = true;
                    }

                    if ( args[0].EndsWith("debug") ) {
                        runRelease = false;
                    }

                    discoverRelease = runRelease;

                    if ( args[0].EndsWith("ignore_tests_w_deps") ) {
                        ignore_tests_w_deps = true;
                    }

                    int iQueue = 1;
                    string filter = args.Length > 1 ? args[1] : "*";
                    if ( args.Length == 3 ) {
                        Console.WriteLine("arg 2 is:" + args[2]);
                        if ( args[2].StartsWith("queue#") ) {
                            try {
                                iQueue = int.Parse(args[2].Split(new[] { '#' }, StringSplitOptions.RemoveEmptyEntries)[1]);
                            } catch ( Exception ) {
                                ret = -1;
                                Console.Error.WriteLine("Unable to parse queue number from " + args[1]);
                                PrintMainUsage();
                                break;
                            }
                        } else {
                            ret = -1;
                            PrintMainUsage();
                            break;
                        }
                    } else {

                    }

                    ret = JobManagerRun(filter, iQueue);
                    break;

                case "help":
                    PrintMainUsage();
                    ret = 0;
                    break;

                default:
                    PrintMainUsage();
                    ret = -1000;
                    break;
            }

            MPICollectiveWatchDog.WatchAtRelease(csMPI.Raw._COMM.WORLD);
            Console.WriteLine("ret b4 finalize = " + ret);
            csMPI.Raw.mpiFinalize();

            System.Environment.Exit(ret);

            return ret;
        }

        static int Main(string[] args) {
            Console.WriteLine($"received {args.Length} arguments.");
            //for(int i = 0; i < args.Length; i++) {
            //    Console.WriteLine($"  arg#{i}  >>>>>>{args[i]}<<<<<<");
            //}
            //Debugger.Launch();

            try {
                return _Main(args, new PublicTests());
            } catch ( Exception e ) {
                // note: this seemingly useless try-catch is here since our test runner server (FDYGITRUNNER)
                // seems to silently fail on all exceptions thrown after MPI init.

                int rank, size;
                if ( csMPI.Raw.Initialized() ) {
                    csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                    csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                } else {
                    rank = 0;
                    size = 0;
                }

                Console.WriteLine("Got some exception: " + e + "(" + e.StackTrace + ")");

                using var stw = new StreamWriter("Exception-" + DateTime.Now.ToString("MMMdd_HHmmss") + "." + rank + "of" + size + ".txt");
                stw.WriteLine("Got some exception: " + e);
                stw.WriteLine(e.StackTrace);
                stw.Flush();
                stw.Close();
                System.Environment.Exit(-667);

                return -667;
            }
        }
    }
}

using BoSSS.Application.BoSSSpad;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using log4net.Appender;
using log4net.Config;
using log4net.Layout;
using MPI.Wrappers;
using NUnit.Framework;
using NUnit.Framework.Interfaces;
using NUnitLite;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading;
using System.Threading.Tasks;


namespace PublicTestRunner {

    /// <summary>
    /// Redirection of fails in <see cref="System.Diagnostics.Debug"/> to NUnit assertions 
    /// </summary>
    class MyListener : System.Diagnostics.TraceListener {
        public override void Write(string message) {
            Console.Error.Write(message);
        }

        public override void WriteLine(string message) {
            Console.Error.WriteLine(message);
        }

        public override void Fail(string message) {
            Console.Error.WriteLine("Assertion Fail: >>>" + message + "<<<");
            Assert.Fail(message);
        }

        public override void Fail(string message, string detailMessage) {
            Console.Error.WriteLine("Assertion Fail: >>>" + message + ": " + ", detailed message is: " + detailMessage);
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
        bool CopyManagedAssembliesCentraly { get; } 
    }

    /// <summary>
    /// 
    /// </summary>
    public class PublicTests : ITestTypeProvider {


        /// <summary>
        /// Note: for better use of parallel resources, try to put the expensive tests in front
        /// </summary>
        virtual public Type[] FullTest {
            get {
                return new Type[] {
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
                        typeof(CutCellQuadrature.Program),
                        typeof(BoSSS.Application.SpecFEM.AllUpTest),
                        typeof(BoSSS.Application.ipViscosity.TestSolution),
                        //typeof(BoSSS.Application.LevelSetTestBench.LevelSetTestBenchMain),
                        //typeof(BoSSS.Application.AdaptiveMeshRefinementTest.AllUpTest),
                        typeof(BoSSS.Application.ExternalBinding.CodeGen.Test),
                        typeof(BoSSS.Application.ExternalBinding.Initializer),
                        //typeof(BoSSS.Application.XNSE_Solver.XNSE), // to expensive for debug
                        typeof(MPITest.Program)
                    };
            }
        }

        /// <summary>
        /// Note: for better use of parallel resources, try to put the expensive tests in front
        /// </summary>
        virtual public Type[] ReleaseOnlyTests {
            get {
                return new Type[] {
                        typeof(BoSSS.Application.XNSERO_Solver.XNSERO),
                        typeof(BoSSS.Application.XNSE_Solver.XNSE),
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
                        typeof(HangingNodesTests.HangingNodesTestMain),
                        typeof(BoSSS.Application.CahnHilliard.CahnHilliardMain),
                        typeof(IntersectingLevelSetTest.AllUpTest),
                    };
            }
        }

        virtual public (Type type, int NoOfProcs)[] MpiFullTests {
            get {
                return new (Type type, int NoOfProcs)[] {
                        (typeof(CDG_Projection_MPI.ConstrainedDGField_Tests), 4),
                        (typeof(CDG_Projection_MPI.ConstrainedDGField_Tests), 2),
                        (typeof(AdvancedSolverTests.TestsMain),4),
                        (typeof(AdvancedSolverTests.MPITests),4),
                        (typeof(MPITest.Program), 4),
                        (typeof(MPITest.Program), 3),
                        (typeof(MPITest.Program), 2),
                        (typeof(BoSSS.Application.SpecFEM.AllUpTest), 4),
                    };
            }
        }

        virtual public (Type type, int NoOfProcs)[] MpiReleaseOnlyTests {
            get {
                return new (Type type, int NoOfProcs)[] {
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
            }
        }

        virtual public DirectoryInfo GetRepositoryBaseDir() {
            DirectoryInfo repoRoot;
            try {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                repoRoot = dir.Parent.Parent.Parent.Parent.Parent.Parent;

                var src = repoRoot.GetDirectories("src").SingleOrDefault();
                var libs = repoRoot.GetDirectories("libs").SingleOrDefault();
                var doc = repoRoot.GetDirectories("doc").SingleOrDefault();
                if (src == null || !src.Exists)
                    return null;
                //throw new Exception();
                if (libs == null || !libs.Exists)
                    return null;
                //throw new Exception();
                if (doc == null || !doc.Exists)
                    return null;
                //throw new Exception();

            } catch (Exception) {
                // must be null, for the nunit3 sub-program when running from somewhere.
                return null;
                //throw new IOException("Unable to find repository root. 'runjobmanger' must be invoked from its default location within the BoSSS git repository. You might want to use 'runjobmanger-ignore_tests_w_deps'");
            }

            return repoRoot;
        }

        virtual public bool CopyManagedAssembliesCentraly => true;
    }

    /// <summary>
    /// Subroutines of the test runner
    /// </summary>
    public static class PublicTestRunnerMain {

        static ITestTypeProvider TestTypeProvider;

        static public bool discoverRelease = true;

        /// <summary>
        /// Timout for job-manager run
        /// </summary>
        public static double TimeOutSec = 24 * 60 * 60; // 24 hours;


        /// <summary>
        /// supposed to ignore tests depending on files in the source code repo;
        /// thereby, we can run the test runner from outside the source repositiory.
        /// </summary>
        static bool ignore_tests_w_deps = false;


        /// <summary>
        /// finds all assemblies which potentially contain tests.
        /// </summary>
        static Assembly[] GetAllAssembliesForTests() {
            using(new FuncTrace()) {
                var R = new HashSet<Assembly>();

                if(TestTypeProvider.FullTest != null) {
                    foreach(var t in TestTypeProvider.FullTest) {
                        //Console.WriteLine("test type: " + t.FullName);
                        var a = t.Assembly;
                        //Console.WriteLine("  assembly: " + a.FullName + " @ " + a.Location);
                        bool added = R.Add(a);
                        //Console.WriteLine("  added? " + added);
                    }
                }



                if(discoverRelease) {
                    if(TestTypeProvider.ReleaseOnlyTests != null) {
                        foreach(var t in TestTypeProvider.ReleaseOnlyTests) {
                            R.Add(t.Assembly);
                        }
                    }
                }

                return R.ToArray();
            }
        }


        

        

        

        //static bool IsReleaseOnlyAssembly


        static (Assembly Asbly, int NoOfProcs)[] GetAllMpiAssemblies() {
            var R = new List<(Assembly Asbly, int NoOfProcs)>();

            if (TestTypeProvider.MpiFullTests != null) {
                foreach (var t in TestTypeProvider.MpiFullTests) {
                    //Console.WriteLine("test type: " + t.type.FullName + " (" + t.NoOfProcs + " procs).");
                    //Console.WriteLine("  assembly: " + t.type.Assembly.FullName + " @ " + t.type.Assembly.Location);
                    bool contains = R.Contains(t, (itm1, itm2) => ((itm1.NoOfProcs == itm2.NoOfProcs) && itm1.Asbly.Equals(itm2.type.Assembly)));
                    if (!contains) {
                        R.Add((t.type.Assembly, t.NoOfProcs));
                    }
                    //Console.WriteLine("  added? " + (!contains));
                }
            }

            if(discoverRelease) {
                if(TestTypeProvider.MpiReleaseOnlyTests != null) {
                    foreach(var t in TestTypeProvider.MpiReleaseOnlyTests) {
                        bool contains = R.Contains(t, (itm1, itm2) => ((itm1.NoOfProcs == itm2.NoOfProcs) && itm1.Asbly.Equals(itm2.type.Assembly)));
                        if(!contains) {
                            R.Add((t.type.Assembly, t.NoOfProcs));
                        }
                    }
                }
            }

            return R.ToArray();
        }


        static string[] LocateFile(string PartialPath) {
            DirectoryInfo repoRoot = TestTypeProvider.GetRepositoryBaseDir();
            if (repoRoot == null)
                return null;

            // if we get here, we probably have access to the repository root directory.
            string[] r = LocateFileRecursive("", repoRoot, PartialPath);
            if (r == null || r.Length <= 0) {
                throw new IOException("unable to find file '" + PartialPath + "'");
            }

            //if (r.Length > 1) {
            //    throw new IOException("The path '" + PartialPath + "' has been found several times: " + r.ToConcatString("", ", ", ";") + " Did you specify the path correctly?");
            //}

            return r;
        }


        static string[] LocateFileRecursive(string RelPath, DirectoryInfo absPath, string SomeFileName) {
            List<string> ret = new List<string>();

            string _SomeFileName = "*" + SomeFileName;

            foreach (var f in absPath.GetFiles()) {
                string RelName = RelPath + f.Name;

                if (RelName.EndsWith(SomeFileName))
                    ret.Add(f.FullName);
                else if (SomeFileName.WildcardMatch(RelName))
                    ret.Add(f.FullName);
                else if (_SomeFileName.WildcardMatch(RelName))
                    ret.Add(f.FullName);

            }

            foreach (var d in absPath.GetDirectories()) {
                ret.AddRange(LocateFileRecursive(RelPath + d.Name + "/", d, SomeFileName));
            }


            return ret.ToArray();
        }

        static (int NoOfTests, string[] tests, string[] shortnames, IDictionary<string,string[]> RequiredFiles4Test, int?[] NumThreads) GetTestsInAssembly(Assembly a, string AssemblyFilter) {
            var r = new List<string>(); // full test names
            var n = new List<int?>(); // number of threads for test
            var l = new List<string>(); // short names 
            var d = new Dictionary<string, string[]>(); // pairs of (full test name | files for the test)

            var g = new HashSet<string>(); // global files for all tests


            var ttt = a.GetTypes();
            foreach (var t in ttt) { // loop over types in assembly...
                if (t.IsClass) {
                    var mmm = t.GetMethods();

                    // first, seach for num_threads defined for the entire class...
                    int? num_threads_class = null;
                    if (t.GetCustomAttribute(typeof(NUnitNumThreads)) != null) {
                        var nta = t.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) as NUnitNumThreads;
                        //if (dc != null)
                        //    throw new NotSupportedException("Palacing a `NUnitFileToCopyHackAttribute` together with `SetUpAttribute` or `OneTimeSetUpAttribute` is not supported (anymore).");

                        num_threads_class = nta.NumThreads;
                    }



                    foreach (var m in mmm) { // loop over methods in type...
                        if (t.IsAbstract && !m.IsStatic)
                            continue;

                        if (!FilterTestMethod(m, AssemblyFilter))
                            continue;

                        var s = new HashSet<string>();



                        if (m.GetCustomAttribute(typeof(SetUpAttribute)) != null
                         || m.GetCustomAttribute(typeof(OneTimeSetUpAttribute)) != null) {
                            var dc = m.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) as NUnitFileToCopyHackAttribute;
                            //if (dc != null)
                            //    throw new NotSupportedException("Palacing a `NUnitFileToCopyHackAttribute` together with `SetUpAttribute` or `OneTimeSetUpAttribute` is not supported (anymore).");

                            if (dc != null) {
                                foreach (string someFile in dc.SomeFileNames) {
                                    g.AddRange(LocateFile(someFile));
                                }
                            }
                        }

                        // search int num_threads is "overridden" for the respective test
                        int? num_threads = num_threads_class;
                        if (m.GetCustomAttribute(typeof(NUnitNumThreads)) != null) {
                            var nta = m.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) as NUnitNumThreads;
                            //if (dc != null)
                            //    throw new NotSupportedException("Palacing a `NUnitFileToCopyHackAttribute` together with `SetUpAttribute` or `OneTimeSetUpAttribute` is not supported (anymore).");

                            num_threads = nta.NumThreads;
                        }


                        //bool testAdded = false;
                        if (m.GetCustomAttribute(typeof(TestAttribute)) != null) {
                            r.Add(t.FullName + "." + m.Name);
                            n.Add(num_threads);
                            l.Add(Path.GetFileNameWithoutExtension(a.ManifestModule.Name) + "#" + m.Name);
                            //Console.WriteLine("Added: " + r.Last());
                            //testAdded = true;
                            //}

                            //if (m.GetCustomAttribute(typeof(TestAttribute)) != null
                            //   || m.GetCustomAttribute(typeof(SetUpAttribute)) != null
                            //   || m.GetCustomAttribute(typeof(OneTimeSetUpAttribute)) != null) {
                            var dc = m.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) as NUnitFileToCopyHackAttribute;


                            if (dc != null) {
                                //Console.WriteLine("Added: " + r.Last() + " depends on " + dc.SomeFileNames[0]);
                                //if(ignore_tests_w_deps && testAdded) {
                                //    // supposed to ignore tests depending on files in the source code repo
                                //    r.RemoveAt(r.Count - 1);
                                //    l.RemoveAt(l.Count - 1);
                                //    continue; // skip this test
                                //}

                                foreach (string someFile in dc.SomeFileNames) {
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
            foreach(var testname in d.Keys) {
                var s = new List<string>();
                s.AddRange(d[testname]);
                s.AddRange(g);

                List<string> FileNamesOnly = new();
                foreach (string filePath in s) {
                    string fileName = Path.GetFileName(filePath);
                    if (FileNamesOnly.Contains(fileName, (string a, string b) => a.Equals(b, StringComparison.InvariantCultureIgnoreCase)))
                        throw new IOException($"Dependent Filename {fileName} is not unique for test assembly {a}. (full Path {filePath}).");
                    FileNamesOnly.Add(fileName);
                }

                d[testname] = s.ToArray();
            }

            if (n.Count != r.Count)
                throw new ApplicationException("internal error");

            return (r.Count, r.ToArray(), l.ToArray(), d, n.ToArray());
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

            if (logger_output != null)
                throw new ApplicationException("Already called."); // is seems this object is designed so that it stores at max one session per lifetime


            tracerfile = new FileStream($"trace_{basename}.txt", FileMode.Create, FileAccess.Write, FileShare.Read);
            tracertxt = new StreamWriter(tracerfile);

            TextWriterAppender fa = new TextWriterAppender();
            fa.ImmediateFlush = true;
            //fa.Writer = Console.Out;
            fa.Writer = tracertxt;
            fa.Layout = new PatternLayout("%date %-5level %logger: %message%newline");
            fa.ActivateOptions();
            BasicConfigurator.Configure(fa);
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
            } catch (Exception) {
                //Console.Error.WriteLine(e.GetType().Name + " during closing of tracing: " + e.Message);
            }
        }

        
        /// <summary>
        /// Decides whether an assembly <paramref name="a"/> matches the filter (wildcard <paramref name="AssemblyFilter"/>) specified by the user
        /// </summary>
        static bool FilterTestAssembly(Assembly a, string AssemblyFilter) {
            if (AssemblyFilter.IsEmptyOrWhite())
                return true;
            string[] sFilters = AssemblyFilter.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var filter in sFilters) {
                string modFilter;
                bool expect = false;
                if(filter.StartsWith("!")) {
                    modFilter = filter.Substring(1);
                    expect = false;
                } else {
                    modFilter = filter;
                    expect = true;
                }

                modFilter = modFilter.Split("\\", StringSplitOptions.RemoveEmptyEntries)[0];

                if (modFilter.WildcardMatch(Path.GetFileNameWithoutExtension(a.Location)) == expect)
                    return true;
            }
            return false;
        }
        

         /// <summary>
        /// Decides whether an assembly <paramref name="m"/> matches the filter (wildcard <paramref name="AssemblyFilter"/>) specified by the user
        /// </summary>
        static bool FilterTestMethod(MethodInfo m, string AssemblyFilter) {

            if(AssemblyFilter.IsEmptyOrWhite())
                return true;
            
            string FullName = m.DeclaringType.Name + "." + m.Name;
            Assembly a = m.DeclaringType.Assembly;

            string[] sFilters = AssemblyFilter.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var filter in sFilters) {
                string modFilter;
                bool expect = false;
                if(filter.StartsWith("!")) {
                    modFilter = filter.Substring(1);
                    expect = false;
                } else {
                    modFilter = filter;
                    expect = true;
                }

                string[] AssemblyNMethodFilter = modFilter.Split("\\", StringSplitOptions.RemoveEmptyEntries);
                string aFilter = AssemblyNMethodFilter[0];

                if(aFilter.WildcardMatch(Path.GetFileNameWithoutExtension(a.Location)) == expect) {

                    if(AssemblyNMethodFilter.Length < 2)
                        return true;

                    string mfilter = AssemblyNMethodFilter[1];

                    
                    if (mfilter.WildcardMatch(FullName) == expect)
                        return true;

                }
            }
            return false;

        }

        static public int BuildYaml() {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out var MpiSize);
            if (MpiSize != 1) {
                throw new NotSupportedException("yaml subprogram must be executed serially");
            }

            using(var tr = new FuncTrace()) {

                // ===================================
                // phase 1: submit jobs
                // ===================================


           
                var allTests = new List<(Assembly ass, string testname, string shortname, string[] depfiles, int NoOfProcs)>();
                {
                    var assln = GetAllAssembliesForTests();
                    if(assln != null) {
                        foreach(var a in assln) {
                            {
                                var allTst4Assi = GetTestsInAssembly(a, null);
                                for(int iTest = 0; iTest < allTst4Assi.NoOfTests; iTest++) {
                                    allTests.Add((a, allTst4Assi.tests[iTest], allTst4Assi.shortnames[iTest], allTst4Assi.RequiredFiles4Test[allTst4Assi.tests[iTest]], 1));
                                }
                            }
                        }
                    }
                }
                {
                    var ParAssln = GetAllMpiAssemblies();
                    if(ParAssln != null) {
                        foreach(var TT in ParAssln) {
                            {

                                var a = TT.Asbly;
                                var allTst4Assi = GetTestsInAssembly(a, null);

                                for(int iTest = 0; iTest < allTst4Assi.NoOfTests; iTest++) {
                                    allTests.Add((a, allTst4Assi.tests[iTest], allTst4Assi.shortnames[iTest], allTst4Assi.RequiredFiles4Test[allTst4Assi.tests[iTest]], TT.NoOfProcs));
                                }
                            }
                        }
                    }
                }

                Console.WriteLine($"Found {allTests.Count} individual tests ({DebugOrReleaseSuffix}):");
                int cnt = 0;
                foreach(var t in allTests) {
                    cnt++;
                    Console.WriteLine($"  #{cnt}: {t.testname}");
                    Console.WriteLine($"     {t.shortname}");
                    Console.WriteLine($"     {t.NoOfProcs} MPI processors.");
                }

                Console.WriteLine($"******* Writing new yaml file ({DateTime.Now}) *******");

                string yamlName;
                 
                yamlName = "jobs.yml";

                using(var YAML = new StreamWriter(yamlName)) {
                    YAML.WriteLine("################################################################################");
                    YAML.WriteLine($"# this is an auto-generated file by {TestTypeProvider.GetType().Assembly.FullName}.");
                    YAML.WriteLine("# any modification might get over-written");
                    YAML.WriteLine($"# created: {DateTime.Now}");
                    YAML.WriteLine($"# user:    {System.Environment.UserName}");
                    YAML.WriteLine($"# system:  {System.Environment.MachineName}");
                    YAML.WriteLine("################################################################################");
                    YAML.WriteLine();

                    //Set Workflow
                    YAML.WriteLine("workflow:");
                    YAML.WriteLine("  rules:");
                    YAML.WriteLine("    - when: always");
                    YAML.WriteLine();

                    //Set Stages
                    YAML.WriteLine("stages:");
                    YAML.WriteLine("  - test");
                    YAML.WriteLine("  - test parallel");
                    YAML.WriteLine();

                    if (allTests.Count == 0)
                    {
                        YAML.WriteLine("EmptyTest:");
                        YAML.WriteLine("  stage: test");
                        YAML.WriteLine("  script:");
                        YAML.WriteLine("    - echo \"Empty\"");
                    }
                    else
                    {
                        //Set job class
                        // ======================================================================
                        //Gitlab yaml sets RUNNER_PATH, RUNNER_EXE, BUILD_DEPENDENCY, ARTIFACT_REF_PATH
                        //Gitlab automatically sets CI_PROJECT_PATH, CI_MERGE_REQUEST_REF_PATH

                        YAML.WriteLine(".Test:");
                        YAML.WriteLine("  before_script:");
                        YAML.WriteLine("    - cd $RUNNER_PATH");
                        YAML.WriteLine("    - bash -c \"chmod +x $RUNNER_EXE\"");
                        YAML.WriteLine("  artifacts:");
                        YAML.WriteLine("    reports:");
                        YAML.WriteLine("      junit: $RUNNER_PATH/TestResult.*");
                        YAML.WriteLine("    expire_in: 2 days");
                        YAML.WriteLine("  needs:");
                        YAML.WriteLine("    - project: $CI_PROJECT_PATH");
                        YAML.WriteLine("      job: $BUILD_DEPENDENCY");
                        YAML.WriteLine("      ref: $ARTIFACT_REF_PATH");
                        YAML.WriteLine("      artifacts: true");
                        YAML.WriteLine();

                        cnt = 0;
                        var checkResFileName = new HashSet<string>();
                        foreach (var t in allTests)
                        {

                            if (t.testname.Contains("TutorialTest"))
                            {
                                Console.WriteLine("skipping: " + t.testname);
                                continue;
                            }


                            if (t.NoOfProcs == 1)
                                YAML.WriteLine(DebugOrReleaseSuffix + "#" + t.shortname + ":" + t.testname + ":");
                            else 
                                YAML.WriteLine(DebugOrReleaseSuffix + "#p" + t.NoOfProcs + "#" + t.shortname + ":" + t.testname + ":");
                            YAML.WriteLine("   extends: .Test");

                            if (t.NoOfProcs == 1)
                                YAML.WriteLine("   stage: test");
                            else
                                YAML.WriteLine("   stage: test parallel");
                            YAML.WriteLine("   script:");
                            if (t.NoOfProcs == 1)
                                YAML.WriteLine($"     - '& ./$RUNNER_EXE nunit3 {Path.GetFileName(t.ass.Location)} --test={t.testname} --result=TestResult.xml'");
                            else
                                YAML.WriteLine($"     - mpiexec -n {t.NoOfProcs} ./$RUNNER_EXE nunit3 {Path.GetFileName(t.ass.Location)} --test={t.testname} --result=TestResult.xml");
                            if (t.NoOfProcs > 1)
                            {
                                YAML.WriteLine("   tags:");
                                YAML.WriteLine($"    - {t.NoOfProcs}cores");
                            }
                            YAML.WriteLine();
                        }
                    }
                }
            }
            return 0;
        }

        /// <summary>
        /// to avoid IO collisions for concurrent runs of the job manager on the same machine (e.g. DEBUG and RELEASE);
        /// appending of the user name avoids "unauthorized access"-exceptions
        /// </summary>
        static Mutex IOsyncMutex = new Mutex(false, "_BoSSS_test_runner_IOmutex_" + System.Environment.UserName);

        /// <summary>
        /// to distinct the internalTestRunner
        /// </summary>
        public static string RunnerPrefix = "Pub";

        

        static public int JobManagerRun(string AssemblyFilter, int ExecutionQueueNo) {

            // ===================================
            // phase 0: setup
            // ===================================


            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out var MpiSize);
            if(MpiSize != 1) {
                throw new NotSupportedException("runjobmanager subprogram must be executed serially");
            }

            BoSSSshell.ReloadExecutionQueues();

            if(ExecutionQueueNo >= BoSSSshell.ExecutionQueues.Count)
                throw new ApplicationException($"Execution queue #{ExecutionQueueNo} does not exist on this machine/account (see configuration file ~/.BoSSS/etc/BatchProcessorConfig.json).");
            BatchProcessorClient bpc = BoSSSshell.ExecutionQueues[ExecutionQueueNo];
            Console.WriteLine($"Using batch queue {ExecutionQueueNo}: {bpc.ToString()}");

            FileStream ServerMutex;
            string DateNtime = null;
            try {
                IOsyncMutex.WaitOne();

                var rnd = new Random(DateTime.Now.Millisecond + typeof(PublicTestRunnerMain).Assembly.Location.GetHashCode() + Directory.GetCurrentDirectory().GetHashCode());
                Thread.Sleep(rnd.Next(10000)); // sleep for a random amount of time to avoid 
                do {
                    DateNtime = DateTime.Now.ToString("MMMdd_HHmmss");
                    string MutexFileName = Path.Combine(bpc.DeploymentBaseDirectory, RunnerPrefix + DebugOrReleaseSuffix + "_" +  DateNtime + ".lock");
                    try {
                        ServerMutex = File.Open(MutexFileName, FileMode.Create, FileAccess.Write, FileShare.None);
                        using(var wrt = new StreamWriter(ServerMutex)) {
                            wrt.WriteLine("Locked by BoSSS test runner at " + DateNtime);
                        }
                    } catch(Exception ee) {
                        Console.Error.WriteLine($"Unable to get lock on {MutexFileName}: {ee}");
                        Console.Error.WriteLine("wait and try again...");
                        ServerMutex = null;
                        Thread.Sleep(rnd.Next(10000));
                    }
                } while(ServerMutex == null);
                Console.WriteLine($"Using prefix'{DateNtime}' for all jobs.");
                if(DateNtime == null)
                    throw new ApplicationException("internal error");
            } catch(Exception e) {
                Console.Error.WriteLine("UNRECOVERABLE EXCEPTION DURING CREATION OF OUTPUT DIRECTORY");
                Console.Error.WriteLine(e.GetType() + ":  " + e.Message);
                Console.Error.WriteLine(e.StackTrace);
                Console.Error.WriteLine("TERMINATING APPLICATION");
                System.Environment.Exit(-666);
            } finally {
                IOsyncMutex.ReleaseMutex();
            }
            Tracer.NamespacesToLog = new string[] { "" };
            InitTraceFile("JobManagerRun-" + DateNtime);

            bool I_StartedMinibatch = false;
            if(bpc is MiniBatchProcessorClient) {
                I_StartedMinibatch = MiniBatchProcessor.Server.StartIfNotRunning(RunExternal: true);
            }

            
            int returnCode = 0;
            using(var tr = new FuncTrace()) {

                // ===================================
                // phase 1: discover tests
                // ===================================

                BoSSSshell.WorkflowMgm.Init("BoSSStst" + DateNtime, bpc);

                // deployment of native libraries
                string NativeOverride;
                if(bpc.DeployRuntime == false) {
                    //
                    // DeployRuntime is false: 
                    // this means that no copy (of the native libraries) occurs for the **individual** jobs
                    // The TestRunner, however copies it centrally, at once, to ensure that it is running using the most recent binaries
                    //
                    var _NativeOverride = new DirectoryInfo(Path.Combine(bpc.DeploymentBaseDirectory, RunnerPrefix + DebugOrReleaseSuffix + "_" + DateNtime + "_amd64"));
                    _NativeOverride.Create();

                    if (bpc.RuntimeLocation != null) {
                        string BosssInstall = BoSSS.Foundation.IO.Utils.GetBoSSSInstallDir();
                        string BosssBinNative = Path.Combine(BosssInstall, "bin", "native", bpc.RuntimeLocation);
                        MetaJobMgrIO.CopyDirectoryRec(BosssBinNative, _NativeOverride.FullName, null);

                        if (bpc is SlurmClient slurm) {
                            NativeOverride = slurm.DeploymentDirectoryAtRemote(_NativeOverride.FullName);
                        } else {
                            NativeOverride = _NativeOverride.FullName;
                        }
                    } else {
                        NativeOverride = null;    
                    }
                } else {
                    NativeOverride = null;
                }
                
                // deployment of assemblies
                string RelManagedPath;
                if(TestTypeProvider.CopyManagedAssembliesCentraly) {
                    string mngdir = RunnerPrefix + DebugOrReleaseSuffix + "_" + DateNtime + "_managed";
                    DirectoryInfo ManagedOverride = new DirectoryInfo(Path.Combine(bpc.DeploymentBaseDirectory, mngdir));
                    ManagedOverride.Create();
                    TestTypeProvider.GetType().Assembly.DeployAt(ManagedOverride);

                    RelManagedPath = "../" + mngdir + "/" + Path.GetFileName(TestTypeProvider.GetType().Assembly.Location);
                } else {
                    RelManagedPath = null;
                }
                


                // collection for all tests:
                var allTests = new List<(Assembly ass, string testname, string shortname, string[] depfiles, int NoOfProcs, int? NumThreads)>();
                
                // Find all serial tests:
                {
                    var assln = GetAllAssembliesForTests();
                    if(assln != null) {
                        foreach(var a in assln) {
                            if(FilterTestAssembly(a, AssemblyFilter)) {
                                var allTst4Assi = GetTestsInAssembly(a, AssemblyFilter);
                                for(int iTest = 0; iTest < allTst4Assi.NoOfTests; iTest++) {
                                    if(ignore_tests_w_deps && allTst4Assi.RequiredFiles4Test[allTst4Assi.tests[iTest]].Length > 0) {
                                        Console.WriteLine($"Skipping all in {a} due to external test dependencies.");
                                        break;
                                    }
                                    allTests.Add((a, allTst4Assi.tests[iTest], allTst4Assi.shortnames[iTest], allTst4Assi.RequiredFiles4Test[allTst4Assi.tests[iTest]], 1, allTst4Assi.NumThreads[iTest]));
                                }
                            }
                        }
                    }
                }

               
                // Find all MPI-parallel tests:
                {
                    var ParAssln = GetAllMpiAssemblies();
                    if(ParAssln != null) {
                        foreach(var TT in ParAssln) {
                            if(FilterTestAssembly(TT.Asbly, AssemblyFilter)) {

                                var a = TT.Asbly;
                                var allTst4Assi = GetTestsInAssembly(a, AssemblyFilter);

                                for(int iTest = 0; iTest < allTst4Assi.NoOfTests; iTest++) {
                                    if(ignore_tests_w_deps && allTst4Assi.RequiredFiles4Test[allTst4Assi.tests[iTest]].Length > 0) {
                                        Console.WriteLine($"Skipping all in {a} due to external test dependencies.");
                                        break;
                                    }

                                    allTests.Add((a, allTst4Assi.tests[iTest], allTst4Assi.shortnames[iTest], allTst4Assi.RequiredFiles4Test[allTst4Assi.tests[iTest]], TT.NoOfProcs, allTst4Assi.NumThreads[iTest]));
                                }
                            }
                        }
                    }
                }

                Console.WriteLine($"Found {allTests.Count} individual tests ({DebugOrReleaseSuffix}):");
                int cnt = 0;
                foreach(var t in allTests) {
                    cnt++;
                    Console.WriteLine($"  #{cnt}: {t.testname}");
                    Console.WriteLine($"     {t.shortname}");
                    Console.WriteLine($"     {t.NoOfProcs} MPI processors.");
                }

                {
                    const string BOSSS_TEST_RUNNER_GODMODE = "c:\\tmp\\godmode.txt";
                    try {
                        var s = File.ReadAllText(BOSSS_TEST_RUNNER_GODMODE);
                        int godval = int.Parse(s);
                        if(godval != 0) {
                            Console.WriteLine("Detected Godmode-Cheatfile. Setting all tests to success.");
                            return 0;
                        }
                    } catch(Exception) { }
              
                }

                // ===================================
                // phase 2: submit jobs
                // ===================================

                Console.WriteLine($"******* Starting job/test deployment/submission ({DateTime.Now}) *******");

                DateTime start = DateTime.Now;

                cnt = 0;
                var AllOpenJobs = new List<(Job job, string ResFile, string testname)>();
                using(new BlockTrace("DEPLOYMENT", tr)) {
                    var checkResFileName = new HashSet<string>();

                    foreach(var t in allTests) {
                        try {
                            cnt++;
                            Console.WriteLine($"Submitting {cnt} of {allTests.Count} ({t.shortname})...");
                            var j = SubmitJob(t.ass, t.testname, t.shortname, bpc, t.depfiles, DateNtime, t.NoOfProcs, t.NumThreads, NativeOverride, RelManagedPath, cnt);
                            if(checkResFileName.Add(j.resultFile) == false) {
                                throw new IOException($"Result file name {j.resultFile} is used multiple times.");
                            }


                            Console.WriteLine($"Successfully submitted {j.j.Name}. \n");
                            AllOpenJobs.Add(j);
                        } catch(Exception e) {
                            Console.Error.WriteLine($"{e.GetType().Name} during job submission: {e.Message}.");
                            returnCode--;
                        }
                    }
                }

                if(returnCode == 0) {
                    Console.WriteLine($"******* All jobs/tests deployed ({DateTime.Now}) *******");
                } else {
                    Console.WriteLine($"******* Deployed finished ({DateTime.Now}) -- SOME DEPLOYMENT(S) FAILED *******");
                }
                // ===================================
                // phase 3: wait until complete...
                // ===================================


                var AllFinishedJobs = new List<(Job job, string ResFile, string testname, JobStatus LastStatus)>();


                
                using(var ot = new StreamWriter("allout-" + DateNtime + "-" + DebugOrReleaseSuffix + ".txt")) {

                    (Job job, string ResFile, string testname, JobStatus LastStatus)[] UpdateFinishedJobs() {
                        using(var trr = new FuncTrace("UpdateFinishedJobs")) {
                            string CurrentDir = Path.GetDirectoryName(typeof(PublicTestRunnerMain).Assembly.Location);

                            var RecentlyFinished = new List<(Job job, string ResFile, string testname, JobStatus LastStatus)>();

                            for(int iJob = 0; iJob < AllOpenJobs.Count; iJob++) {
                                var jj = AllOpenJobs[iJob];
                                var s = jj.job.Status;

                                {
                                    string resultArg = "--result=";
                                    string resArg = jj.job.EnvironmentVars.Values.Single(arg => arg.StartsWith(resultArg));
                                    string _resFile = resArg.Replace(resultArg, "");
                                    if(_resFile != jj.ResFile) {
                                        throw new ApplicationException("internal mismatch in result file name");
                                    }
                                }


                                if (s == JobStatus.FailedOrCanceled) {
                                    Console.WriteLine(" ------------------- Job Failed reason:");
                                    var s1 = jj.job.GetStatus(true);
                                    if (s1 != s) {
                                        Console.WriteLine("changed its mind to: " + s1);
                                        s = s1;
                                    }

                                    if (jj.job.SubmitCount < jj.job.RetryCount) {
                                        Console.WriteLine("Trying once again with failed job...");
                                        jj.job.Reactivate();
                                        continue;
                                    }
                                }


                                if(s == JobStatus.FailedOrCanceled || s == JobStatus.FinishedSuccessful) {
                                    


                                    // message:
                                    if(s == JobStatus.FinishedSuccessful)
                                        Console.WriteLine(s + ": " + jj.job.Name + " // " + jj.testname + " (" + DateTime.Now + ")");
                                    else
                                        Console.WriteLine(s + ": " + jj.job.Name + " // " + jj.testname + " at " + jj.job.LatestDeployment.DeploymentDirectory.FullName + " (" + DateTime.Now + ")");



                                    // copy stdout and stderr to logfile
                                    LogResultFile(ot, jj.job, jj.testname, jj.ResFile);

                                    // copy xml result file
                                    using(new BlockTrace("copy_nunit_xml_result", trr)) {
                                        try {
                                            string[] sourceFiles = Directory.GetFiles(jj.job.LatestDeployment.DeploymentDirectory.FullName, "result-*.xml");
                                            sourceFiles = sourceFiles.Cat(Directory.GetFiles(jj.job.LatestDeployment.DeploymentDirectory.FullName, "*.html")); // output of jupyter notebooks

                                            foreach(var orig in sourceFiles) {
                                                string n = Path.GetFileName(orig);
                                                string dest = Path.Combine(CurrentDir, n);
                                                File.Copy(orig, dest, true);
                                            }
                                        } catch(IOException ioe) {
                                            Console.Error.WriteLine(ioe.GetType().Name + ": " + ioe.Message);
                                            returnCode--;
                                        }
                                    }
                                    // delete deploy directory
                                    using(new BlockTrace("delete_deploy_dir", trr)) {
                                        if(s == JobStatus.FinishedSuccessful) {
                                            try {
                                                Directory.Delete(jj.job.LatestDeployment.DeploymentDirectory.FullName, true);
                                            } catch(Exception e) {
                                                Console.Error.WriteLine($"{e.GetType().Name}: {e.Message}");
                                            }
                                        }
                                    }


                                    // move job to 'finished' list
                                    var X = (jj.job, jj.ResFile, jj.testname, s);
                                    AllFinishedJobs.Add(X);
                                    RecentlyFinished.Add(X);

                                    AllOpenJobs.RemoveAt(iJob);
                                    iJob--;
                                }
                            }

                            return RecentlyFinished.ToArray();
                        }
                    }

                    UpdateFinishedJobs();
                    double RestTime = Math.Max(1, TimeOutSec - (DateTime.Now - start).TotalSeconds);
                    while(RestTime > 1.0 && AllOpenJobs.Count > 0) {
                        using(new BlockTrace("Sleeping", tr)) {
                            Thread.Sleep(2 * 60 * 1000); // sleep for 2 minutes
                        }
                        var ll = UpdateFinishedJobs();
                        RestTime = Math.Max(1, TimeOutSec - (DateTime.Now - start).TotalSeconds);
                        Console.WriteLine("Remaining minutes until timeout: " + Math.Round(RestTime / 60.0));
                        Console.Write("     Waiting for: ");
                        int i = 0;
                        foreach(var j in AllOpenJobs) {
                            Console.Write(j.job.Name);
                            Console.Write(" ");
                            i++;
                            if(i >= 8)
                                break;
                        }
                        if(i < AllOpenJobs.Count)
                            Console.WriteLine($"... (and {AllOpenJobs.Count - i} more.)");
                        else
                            Console.WriteLine();
                    }
                }

                // ===================================
                // phase 4: summary
                // ===================================


                if(AllOpenJobs.Count == 0) {
                    Console.WriteLine("----------------------------------------------------------------------------------------------------");
                    Console.WriteLine($"All jobs/tests finished ({DateTime.Now}) - Summary:");
                    Console.WriteLine("----------------------------------------------------------------------------------------------------");
                } else {
                    Console.WriteLine("----------------------------------------------------------------------------------------------------");
                    Console.WriteLine($"Reached Time-out of {TimeOutSec / 60} Minutes; {AllOpenJobs.Count} jobs/tests still running; ({DateTime.Now}) - Summary:");
                    Console.WriteLine("----------------------------------------------------------------------------------------------------");

                    returnCode -= 1;
                }

                int SuccessfulFinishedCount = 0;
                foreach(var jj in AllFinishedJobs.Where(ttt => ttt.LastStatus == JobStatus.FinishedSuccessful)) {
                    Console.WriteLine($"{jj.job.Status}: {jj.job.Name} // {jj.testname}");
                    SuccessfulFinishedCount++;
                }

                int OtherStatCount = 0;
                foreach(var jj in AllFinishedJobs.Where(ttt => ttt.LastStatus != JobStatus.FinishedSuccessful)) {
                    Console.WriteLine($"{jj.job.Status}: {jj.job.Name} // {jj.testname} at {jj.job.LatestDeployment.DeploymentDirectory.FullName}");
                    returnCode -= 1;
                    OtherStatCount++;
                }
                foreach(var jj in AllOpenJobs) {
                    Console.WriteLine($"{jj.job.Status}: {jj.job.Name} // {jj.testname} at {jj.job.Status}");
                    returnCode -= 1;
                    OtherStatCount++;
                }

                // very final message:
                if(SuccessfulFinishedCount == (AllOpenJobs.Count + AllFinishedJobs.Count)) {
                    Console.WriteLine("All tests/jobs finished successfully.");
                    Console.WriteLine("SUCCESS.");

                    if(returnCode != 0) {
                        Console.Error.WriteLine("Ignoring some other error occurred (maybe IO error after successful test run) -- check output log;");
                        //Console.Error.WriteLine("FAILURE.");
                    } else {
                    }

                    returnCode = 0;

                } else {
                    Console.Error.WriteLine($"Only {SuccessfulFinishedCount} tests/jobs finished successfully -- {OtherStatCount} have other states.");
                    Console.Error.WriteLine("FAILURE.");
                }
                Console.WriteLine($"{DateTime.Now}");
            }


            if(I_StartedMinibatch)
                MiniBatchProcessor.Server.SendTerminationSignal(WaitForOtherJobstoFinish: false);

            CloseTracing();

            return returnCode;
        }

        private static void LogResultFile(StreamWriter ot, Job j, string testname, string ResFile) {
            using (new FuncTrace("LogResultFile")) {
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
                if (j.NumberOfMPIProcs <= 1) {
                    ot.WriteLine("#### Result File:      " + ResFile);
                    ot.WriteLine("####    exists?        " + File.Exists(Path.Combine(j.LatestDeployment.DeploymentDirectory.FullName, ResFile)));
                } else {
                    for (int i = 0; i < sz; i++) {
                        if (i == 0)
                            ot.WriteLine("#### Result Files:     " + MpiResFileNameMod(i, sz, ResFile));
                        else
                            ot.WriteLine("####                   " + MpiResFileNameMod(i, sz, ResFile));
                    }
                    for (int i = 0; i < sz; i++) {
                        if (i == 0)
                            ot.WriteLine("####    exists?        " + File.Exists(Path.Combine(j.LatestDeployment.DeploymentDirectory.FullName, MpiResFileNameMod(i, sz, ResFile))));
                        else
                            ot.WriteLine("####                   " + File.Exists(Path.Combine(j.LatestDeployment.DeploymentDirectory.FullName, MpiResFileNameMod(i, sz, ResFile))));
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
                } catch (Exception e) {
                    ot.WriteLine($"{e.GetType().Name} during reading of stdout stream: {e.Message}");
                    ot.WriteLine(e.StackTrace);
                    stdout = "";
                }
                ot.WriteLine("]]]");

                if(stdout == null)
                    stdout = "";

                using (var str = new StringReader(stdout)) {
                    string magic = "arg #3 override from environment variable 'BOSSS_ARG_3': --result=";
                    for (string line = str.ReadLine(); line != null; line = str.ReadLine()) {
                        if (line.StartsWith(magic)) {
                            string file = line.Replace(magic, "");
                            if (file != ResFile) {
                                throw new ArgumentException("Internal result file mismatch: " + file + " vs. " + ResFile + " on job " + j.LatestDeployment.BatchProcessorIdentifierToken);
                            }

                            break;
                        }
                    }
                }

                try {
                    string stderr = j.Stderr;

                    if (stderr.IsEmptyOrWhite()) {
                        ot.WriteLine("[[[Empty Error Stream: this is good!]]]");
                    } else {
                        ot.WriteLine("[[[Stderr:");
                        ot.WriteLine(stderr);
                        ot.WriteLine("]]]");
                    }
                } catch (Exception e) {

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

        static public (Job j, string resultFile, string name) SubmitJob(
            Assembly a,
            string TestName, string Shortname,
            BatchProcessorClient bpc,
            string[] AdditionalFiles,
            string prefix,
            int NoOfMpiProcs,
            int? NumThreads,
            string nativeOverride,
            string TestRunnerRelPath,
            int cnt) {
            using (new FuncTrace()) {

                // define unique name (not to long) for the job
                string dor = DebugOrReleaseSuffix;
                string final_jName;
                {
                    string jName;
                    if (NoOfMpiProcs <= 1)
                        jName = $"{prefix}-{dor}-{Shortname}";
                    else
                        jName = $"{prefix}-{dor}-p{NoOfMpiProcs}-{Shortname}";

                    if (jName.Length > 127) {
                        // Name length limit set by MS HPC cluster
                        jName = jName.Substring(0, 127);
                    }
                    int counter = 2;
                    final_jName = jName;
                    while (BoSSSshell.WorkflowMgm.AllJobs.ContainsKey(final_jName)) {
                        string suffix = "_" + counter;
                        counter++;
                        if (jName.Length + suffix.Length > 127) {
                            final_jName = jName.Substring(0, 127 - suffix.Length);
                        } else {
                            final_jName = jName;
                        }
                        final_jName = final_jName + suffix;
                    }
                }

                // create job
                Job j = new Job(final_jName, TestTypeProvider.GetType());
                j.SessionReqForSuccess = false;
                j.RetryCount = 3;
                string resultFile = $"result-{dor}-{cnt}.xml";
                j.MySetCommandLineArguments("nunit3", Path.GetFileNameWithoutExtension(a.Location), $"--test={TestName}", $"--result={resultFile}");
                foreach (var f in AdditionalFiles) {
                    j.AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(File.ReadAllBytes(f), Path.GetFileName(f)));
                }
                if(BOSSS_RUNTESTFROMBACKUP_ENVVAR) {
                    j.AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(File.ReadAllBytes("BOSSS_RUNTESTFROMBACKUP.txt"), "BOSSS_RUNTESTFROMBACKUP.txt"));
                }
                if (nativeOverride != null) {
                    j.EnvironmentVars.Add(BoSSS.Foundation.IO.Utils.BOSSS_NATIVE_OVERRIDE, nativeOverride);
                }
                j.NumberOfMPIProcs = NoOfMpiProcs;
                if(NumThreads != null) {
                    j.NumberOfThreads = NumThreads.Value;
                }
                if(TestRunnerRelPath != null)
                    j.EntryAssemblyRedirection = TestRunnerRelPath;
                j.Activate(bpc);
                return (j, resultFile, TestName);
            }
        }

        public static void DeleteResultFiles() {
            using (var tr = new FuncTrace("DeleteResultFiles")) {
                string CurrentDir = Path.GetDirectoryName(typeof(PublicTestRunnerMain).Assembly.Location);
                string[] FilesToDel = Directory.GetFiles(CurrentDir, "result-*.xml");



                foreach (var f in FilesToDel) {

                    File.Delete(f);

                }
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
            using (var tr = new FuncTrace()) {
                var r = GetTestsInAssembly(a, filter);

                var dir = Directory.GetCurrentDirectory();
                tr.Info("Current dir: " + dir);

                foreach (var t in r.tests) { 
                    foreach (var fOrigin in r.RequiredFiles4Test[t]) {
                        tr.Info("Origin file: " + fOrigin);
                        if (File.Exists(fOrigin)) {
                            string fDest = Path.Combine(dir, Path.GetFileName(fOrigin));
                            tr.Info("Destination file: " + fDest);
                            File.Copy(fOrigin, fDest, true);
                        } else {
                            tr.Info($"Origin file {fOrigin} NOT FOUND!");
                        }
                    }
                }
            }
        }



        /// <summary>
        /// Runs all tests serially
        /// </summary>
        static int RunNunit3Tests(string AssemblyFilter, string[] args) {
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out var MpiRank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out var MpiSize);
            ilPSP.Tracing.Tracer.NamespacesToLog = new string[] { "" };
            InitTraceFile($"Nunit3.{DateTime.Now.ToString("MMMdd_HHmmss")}.{MpiRank}of{MpiSize}");

            Console.WriteLine($"Running an NUnit test on {MpiSize} MPI processes ...");
            Tracer.NamespacesToLog_EverythingOverrideTestRunner = true;

            using (var ftr = new FuncTrace()) {
                Assembly[] assln = GetAllAssembliesForTests();

                if(MpiSize != 1) {
                    // this seems some parallel run
                    // we have to fix the result argument

                    if(args.Where(a => a.StartsWith("--result=")).Count() != 1) {
                        throw new ArgumentException("MPI-parallel NUnit runs require the '--result' argument.");
                    }

                    int i = args.IndexWhere(a => a.StartsWith("--result="));
                    string arg_i = args[i];
                    string resFileName = arg_i.Replace("--result=", "");
                    args[i] = "--result=" + MpiResFileNameMod(MpiRank, MpiSize, resFileName);


                    var parAssis = GetAllMpiAssemblies();
                    foreach(var t in parAssis) {
                        Assembly a = t.Asbly;
                        if(!assln.Contains(a))
                            a.AddToArray(ref assln);
                    }


                    int ii = 0;
                    foreach(var a in assln) {
                        ftr.Info("Assembly #" + ii + ": " + a.ToString());
                        ii++;
                    }
                }

                int count = 0;
                bool ret = false;
                foreach(var a in assln) {
                    if(!FilterTestAssembly(a, AssemblyFilter)) {
                        continue;
                    }

                    if(GetTestsInAssembly(a, AssemblyFilter).NoOfTests <= 0) {
                        Console.WriteLine("Matching Assembly search string, but none of the methods match. (wrong wildcard?)");
                        continue;
                    }

                    Console.WriteLine("Matching assembly: " + a.Location);
                    ftr.Info("found Assembly #" + count + ": " + a.Location);
                    count++;

                    if(MpiRank == 0) {
                        MegaMurxPlusPlus(a, AssemblyFilter);
                    }




                    Console.WriteLine("Waiting for all processors to catch up BEFORE starting test(s)...");
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                    Console.WriteLine("All Here.");

                    int r;
                    using(var bt = new BlockTrace("RUNNING_TEST", ftr)) {
                        var tr = new TextRunner(a);
                        r = tr.Execute(args);
                        //var summary = tr.Summary;
                        //File.WriteAllText("a.txt", summary.ToString());
                        //ftr.Info("Test summary: " + summary.ToString());
                    }

                    using(var bt = new BlockTrace("StdOut/StdErr reset", ftr)) {
                        Console.SetOut(new StreamWriter(Console.OpenStandardOutput()));
                        Console.SetError(new StreamWriter(Console.OpenStandardError()));

                        

                        bt.Info("Waiting for all processors to catch up AFTER running test(s)...");
                        //csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                        MPICollectiveWatchDog.WatchAtRelease(csMPI.Raw._COMM.WORLD, 10*60, true); // 10 Minute timeout for all processors to arrive.
                        bt.Info("All Here.");

                        
                        //var ar = new AutoRun(a);
                        //int r = ar.Execute(args);

                        int[] all_rS = r.MPIAllGather();
                        for(int rnk = 0; rnk < all_rS.Length; rnk++) {
                            bt.Info($"Rank {rnk}: NUnit returned code " + r);
                        }


                        {
                            string currentDirectory = Directory.GetCurrentDirectory();
                            string[] pltFiles = Directory.GetFiles(currentDirectory, "*.plt");

                            long totalSizeBytes = 0;

                            foreach (string pltFile in pltFiles) {
                                FileInfo fileInfo = new FileInfo(pltFile);
                                totalSizeBytes += fileInfo.Length;
                            }

                            double totalSizeGigabytes = totalSizeBytes / (1024.0 * 1024 * 1024);

                            if (totalSizeGigabytes > 2.0) {
                                bt.Error("Test produced more than 2 Gigabyte of plt-files -- please check!");
                                //throw new IOException("Test produced more than 2 Gigabyte of plt-files -- please check!");
                            }
                        }

                    }

                    ftr.Info($"failstate before most recent code {ret} (false means OK, r = {r})");
                    ret = ret | (r != 0);
                    ftr.Info($"failstate after most recent code {ret} (false means OK, r = {r})");
                }

                {
                    ftr.Info("Found  " + count + " assemblies in total");
                    
                    if(count <= 0) {
                        Console.WriteLine("Found no assembly matching: " + AssemblyFilter + " (hint: don't provide a filename extension, e.g. '.dll' or '.exe'; assembly names are compared without using an extension, e.g. 'XNSE_Solver', not 'XNSE_Solver.exe' or 'XNSE_Solver.dll'.)");
                        return -1;
                    }
                }


                Console.WriteLine();
                ftr.Info($"failstate all tests: {ret} (false means OK)");
                return ret ? -1 : 0;
            }
        }

        private static string MpiResFileNameMod(int MpiRank, int MpiSize, string resFileName) {
            if (MpiRank > MpiSize)
                throw new ArgumentException();

            if (MpiSize == 1)
                return resFileName;

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
            Console.WriteLine("         help          : prints this message.");
            Console.WriteLine("         yaml          : write 'jobs.yml' file for Gitlab.");
            Console.WriteLine("  and FILTER selects some assembly, i.e. DerivativeTests.exe; it can be ");
            Console.WriteLine("  a wildcard, i.e. use * for submitting all tests.");

        }



        static public bool BOSSS_RUNTESTFROMBACKUP_ENVVAR = false;

        /// <summary>
        /// the real main-function
        /// </summary>
        /// <param name="args"></param>
        /// <param name="ttp">
        /// A hook to find tests within the entire heap of assemblies.
        /// </param>
        /// <returns></returns>
        public static int _Main(string[] args, ITestTypeProvider ttp) {
            
            if (args.Length > 5) {
                Console.WriteLine($"Warning: got {args.Length} arguments -- are you using this right?");

                Console.WriteLine("No of args: " + args.Length);
                //int i = 0;
                foreach (string arg in args) {
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
           
            if (args.Length < 1) {
                Console.WriteLine("Insufficient number of arguments.");
                PrintMainUsage();
                return -7777;
            }

            if (System.Environment.GetEnvironmentVariable("BOSSS_RUNTESTFROMBACKUP").IsEmptyOrWhite() == false) {
                BOSSS_RUNTESTFROMBACKUP_ENVVAR = true;
                File.WriteAllText("BOSSS_RUNTESTFROMBACKUP.txt", "Helo, Suckers!");
                Console.WriteLine("trying to forward the BOSSS_RUNTESTFROMBACKUP hack via additional deployment files...");
            }

            BoSSS.Solution.Application.InitMPI();
            ilPSP.Environment.InitThreading(true, null);
            
            int ret = -1;
            switch (args[0]) {
                case "nunit3":
                if (args.Length < 2) {
                    Console.WriteLine("Insufficient number of arguments.");
                    PrintMainUsage();
                    return -7777;
                }
                ret = RunNunit3Tests(args[1], args.Skip(2).ToArray());
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

                if (args[0].EndsWith("release"))
                    runRelease = true;
                if (args[0].EndsWith("debug"))
                    runRelease = false;
                discoverRelease = runRelease;

                if (args[0].EndsWith("ignore_tests_w_deps"))
                    ignore_tests_w_deps = true;

                int iQueue = 1;
                string filter = args.Length > 1 ? args[1] : "*";
                if (args.Length == 3) {
                    Console.WriteLine("arg 2 is:" + args[2]);
                    if (args[2].StartsWith("queue#")) {
                        try {
                            iQueue = int.Parse(args[2].Split(new[] { '#' }, StringSplitOptions.RemoveEmptyEntries)[1]);
                        } catch (Exception) {
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

                case "yaml":
#if DEBUG
                discoverRelease = false;
#else
                discoverRelease = true;
#endif
                ret = BuildYaml();
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


            Console.WriteLine("ret b4 finalize = " + ret);
            csMPI.Raw.mpiFinalize();


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
            } catch(Exception e) {
                // note: this seemingly useless try-catch is here since our test runner server (FDYGITRUNNER)
                // seems to silently fail on all exceptions thrown after MPI init.

                int rank, size;
                if (csMPI.Raw.Initialized()) {
                    csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                    csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                } else {
                    rank = 0;
                    size = 0;
                }

                Console.WriteLine("Got some exception: " + e);

                using (var stw = new StreamWriter("Exception-" + DateTime.Now.ToString("MMMdd_HHmmss") + "." + rank + "of" + size + ".txt")) {
                    stw.WriteLine("Got some exception: " + e);
                    stw.WriteLine(e.StackTrace);
                    stw.Flush();
                    stw.Close();
                }


                return -667;
            }
        }
    }
}

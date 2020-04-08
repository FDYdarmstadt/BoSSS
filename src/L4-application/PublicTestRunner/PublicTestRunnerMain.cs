using BoSSS.Application.BoSSSpad;
using ilPSP;
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
            //Assert.Fail(message);
        }

        public override void Fail(string message, string detailMessage) {
            //Assert.Fail(message + ", details: " + detailMessage);
        }
    }


    static class PublicTestRunnerMain {

        /// <summary>
        /// List of tests that should be executed in DEBUG and RELEASE; referencing any type of the assembly will do.
        /// </summary>
        static Type[] FullTestTypes = new Type[] {
            typeof(BoSSS.Application.DerivativeTest.DerivativeTestMain),
            typeof(BoSSS.Application.SipPoisson.SipPoissonMain),
            typeof(BoSSS.Application.Matrix_MPItest.AllUpTest),
            typeof(BoSSS.Application.ElementTests.ElementTests),
            typeof(BoSSS.Application.DatabaseTests.DatabaseTestsProgram),
            typeof(CutCellQuadrature.Program),
            typeof(BoSSS.Application.XDGTest.UnitTest),
            typeof(BoSSS.Application.SpecFEM.AllUpTest),
            typeof(BoSSS.Application.ipViscosity.TestSolution),
            typeof(BoSSS.Application.MultigridTest.MultigridMain),
            typeof(BoSSS.Application.ZwoLsTest.AllUpTest),
            typeof(BoSSS.Application.XdgTimesteppingTest.XdgTimesteppingMain),
            //typeof(BoSSS.Application.LevelSetTestBench.LevelSetTestBenchMain),
            typeof(BoSSS.Application.XdgPoisson3.XdgPoisson3Main),
            typeof(BoSSS.Application.AdaptiveMeshRefinementTest.AllUpTest),
            typeof(BoSSS.Application.ExternalBinding.CodeGen.Test),
            typeof(BoSSS.Application.ExternalBinding.Initializer),
            typeof(MPITest.Program)
        };

        static Type[] ReleaseOnlyTests = new Type[] {
            typeof(BoSSS.Application.TutorialTests.AllUpTest),
            typeof(CNS.Program),
            typeof(BoSSS.Application.TutorialTests.AllUpTest),
            typeof(QuadratureAndProjectionTest.QuadratueAndProjectionTest),
            typeof(BoSSS.Application.XdgNastyLevsetLocationTest.AllUpTest),
            typeof(LTSTests.Program),
            //typeof(BoSSS.Application.XNSE_ViscosityAgglomerationTest.XNSE_ViscosityAgglomerationTestMain),
            typeof(NSE_SIMPLE.SIMPLESolver),
            typeof(BoSSS.Application.IBM_Solver.IBM_SolverMain),
            typeof(ALTSTests.Program),
            typeof(BoSSS.Application.XNSE_Solver.XNSE_SolverMain)
        };


        static Assembly[] GetAllAssemblies() {
            var R = new HashSet<Assembly>();

            foreach(var t in FullTestTypes) {
                R.Add(t.Assembly);
            }
#if !DEBUG
            foreach (var t in ReleaseOnlyTests) {
                R.Add(t.Assembly);
            }
#endif
            return R.ToArray();
        }

        static string LocateFile(string SomeFileName) {
            DirectoryInfo repoRoot;
            try {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                repoRoot = dir.Parent.Parent.Parent.Parent.Parent;

                var src = repoRoot.GetDirectories("src").SingleOrDefault();
                var libs = repoRoot.GetDirectories("libs").SingleOrDefault();
                var doc = repoRoot.GetDirectories("doc").SingleOrDefault();

                if (src == null || !src.Exists)
                    return null;
                if (libs == null || !libs.Exists)
                    return null;
                if (doc == null || !doc.Exists)
                    return null;

            } catch (Exception) {
                return null; // unable to find file
            }

            // if we get here, we probably have access to the repository root directory.
            string[] r = LocateFileRecursive("", repoRoot, SomeFileName);
            if(r == null || r.Length <= 0) {
                throw new IOException("unable to find file '" + SomeFileName  + "'"); 
            }
            if(r.Length > 1) {
                throw new IOException("found multiple matches for '" + SomeFileName + "'");
            }

            return r[0];
        }


        static string[] LocateFileRecursive(string RelPath, DirectoryInfo absPath, string SomeFileName) {
            List<string> ret = new List<string>();


            foreach(var f in absPath.GetFiles()) {
                string RelName = RelPath + f.Name;

                if (RelName.EndsWith(SomeFileName))
                    ret.Add(f.FullName);
                else if (SomeFileName.WildcardMatch(RelName))
                    ret.Add(f.FullName);

            }

            foreach(var d in absPath.GetDirectories()) {
                ret.AddRange(LocateFileRecursive(RelPath + d.Name + "/", d, SomeFileName));
            }


            return ret.ToArray();
        }



        static (string[] tests, string[] RequiredFiles) GetTestsInAssembly(Assembly a) {
            var r = new List<string>();
            var s = new HashSet<string>();

            var ttt = a.GetTypes();
            foreach(var t in ttt) {
                if(t.IsClass) {
                    var mmm = t.GetMethods();

                    foreach(var m in mmm) {
                        if(m.GetCustomAttribute(typeof(TestAttribute)) != null) {
                            r.Add(t.FullName + "." +  m.Name);
                        }

                        if(m.GetCustomAttribute(typeof(TestAttribute)) != null
                           || m.GetCustomAttribute(typeof(SetUpAttribute)) != null
                           || m.GetCustomAttribute(typeof(OneTimeSetUpAttribute)) != null) {
                            var dc = m.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) as NUnitFileToCopyHackAttribute;

                            if(dc != null) {
                                foreach(string someFile in dc.SomeFileNames) {
                                    string filepath = LocateFile(someFile);
                                    s.Add(filepath);
                                }
                            }
                        }
                    }
                }
            }

            return (r.ToArray(), s.ToArray());
        }


        class MyTestFilter : NUnit.Framework.Internal.TestFilter {
            public override TNode AddToXml(TNode parentNode, bool recursive) {
                throw new NotImplementedException();
            }

            public override bool Match(ITest test) {
                throw new NotImplementedException();
            }
        }


        static public int JobManagerRun(string AssemblyFilter) {

            // ===================================
            // phase 1: submit jobs
            // ===================================

            string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss");

            InteractiveShell.ReloadExecutionQueues();
            InteractiveShell.WorkflowMgm.Init("BoSSStst" + DateNtime);

            BatchProcessorClient bpc = InteractiveShell.ExecutionQueues[1];

            var allTests = new List<(Assembly ass, string testname, string[] depfiles)>();
            {
                var assln = GetAllAssemblies();
                foreach (var a in assln) {
                    if (!AssemblyFilter.IsEmptyOrWhite()) {
                        if (!AssemblyFilter.WildcardMatch(Path.GetFileName(a.Location)))
                            continue;
                    }


                    var (allTest, depfiles) = GetTestsInAssembly(a);

                    foreach (var t in allTest) {
                        allTests.Add((a, t, depfiles));
                    }
                }
            }

            var allJobs = new List<(Job job, string ResFile)>();
            foreach(var t in allTests) {
                var j = JobManagerRun(t.ass, t.testname, bpc, t.depfiles, DateNtime);
                allJobs.Add(j);
            }

            // ===================================
            // phase 2: wait until complete...
            // ===================================

            while (InteractiveShell.WorkflowMgm.BlockUntilAnyJobTerminate(out var job, PollingIntervallSeconds: 120) > 0) {

                if(job != null) {
                    Console.WriteLine("just finished: " + job.Name + ": " + job.Status);
                }
            }
            Thread.Sleep(10000);
            Console.WriteLine("----------------------------------");
            Console.WriteLine("All jobs finished - Summary:");
            Console.WriteLine("----------------------------------");
            foreach (var j in allJobs) {
                Console.WriteLine(j.job.ToString());
            }

            // ===================================
            // phase 3: collect files
            // ===================================

            int returnCode = 0;

            string CurrentDir = Path.GetDirectoryName(typeof(PublicTestRunnerMain).Assembly.Location);

            foreach (var j in allJobs) {
                //Console.WriteLine(j.ToString());

                if (j.job.Status != JobStatus.FinishedSuccessful)
                    returnCode--;

                try {
                    string[] sourceFiles = Directory.GetFiles(j.job.DeploymentDirectory, "result-*.xml");

                    foreach (var orig in sourceFiles) {
                        string n = Path.GetFileName(orig);
                        string dest = Path.Combine(CurrentDir, n);
                        File.Copy(orig, dest);
                    }
                } catch(IOException ioe) {
                    Console.Error.WriteLine(ioe.GetType().Name + ": " + ioe.Message);
                    returnCode--;
                }
            }

            using (var ot = new StreamWriter("allout-" + DateNtime + ".txt")) {
                foreach (var jj in allJobs) {
                    var j = jj.job;
                    var FullResultFile = Path.Combine(j.DeploymentDirectory, jj.ResFile);
                    ot.WriteLine("##fdhgjegf763748trfhe8hurdsinf598ugf498jvhsn*hbbvc#####!################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("####  " + j.Name);
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("#### Deploy directory: " + j.DeploymentDirectory);
                    ot.WriteLine("#### Status:           " + j.Status);
                    ot.WriteLine("#### ID:               " + j.BatchProcessorIdentifierToken);
                    ot.WriteLine("#### Result File:      " + jj.ResFile);
                    ot.WriteLine("####    exists?        " + File.Exists(FullResultFile));
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("Stdout: ");
                    ot.WriteLine(j.Stdout);

                    string stderr = j.Stderr;
                    if (stderr.IsEmptyOrWhite()) {
                        ot.WriteLine("[[[Empty Error Stream: this is good!]]]");
                    } else {
                        ot.WriteLine("[[[Stderr:");
                        ot.WriteLine(stderr);
                        ot.WriteLine("]]]");
                    }

                    ot.WriteLine();
                    ot.WriteLine();
                    ot.WriteLine();
                }
            }

            return returnCode;
        }

        static string DebugOrReleaseSuffix {
            get {
                string dor;
#if DEBUG
                dor = "DEBUG";
#else
                dor = "RELEASE";
#endif
                return dor;
            }
        }

        static public (Job j, string resFileName) JobManagerRun(Assembly a, string TestName, BatchProcessorClient bpc, string[] AdditionalFiles, string prefix) {
            string dor = DebugOrReleaseSuffix;
            Job j = new Job($"{prefix}-{TestName}-{dor}", typeof(PublicTestRunnerMain));
            string resultFile = $"result-{TestName}-{dor}.xml";
            j.MySetCommandLineArguments("--nunit3", Path.GetFileName(a.Location), $"--test={TestName}", $"--result={resultFile}");

            foreach (var f in AdditionalFiles) {
                j.AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(File.ReadAllBytes(f), Path.GetFileName(f)));
            }
            j.Activate(bpc);
            return (j, resultFile);
        }


        public static void DeleteResultFiles() {
            string CurrentDir = Path.GetDirectoryName(typeof(PublicTestRunnerMain).Assembly.Location);
            string[] FilesToDel = Directory.GetFiles(CurrentDir, "result-*.xml");



            foreach (var f in FilesToDel) {

                File.Delete(f);

            }
        }


        static void MegaMurxPlusPlus(Assembly a) {
            var r = GetTestsInAssembly(a);

            var dir = Directory.GetCurrentDirectory();

            foreach (var fOrigin in r.RequiredFiles) {
                if(File.Exists(fOrigin)) {
                    string fDest = Path.Combine(dir, Path.GetFileName(fOrigin));

                    File.Copy(fOrigin, fDest, true);


                }

            }

        }



        /// <summary>
        /// Runs all tests serially
        /// </summary>
        static int RunSerial(string AssemblyFilter, string[] args) {
            var assln = GetAllAssemblies();


            int count = 0;
            bool ret = false;
            foreach(var a in assln) {
                if(!AssemblyFilter.IsEmptyOrWhite()) {
                    if(!AssemblyFilter.WildcardMatch(Path.GetFileName(a.Location)))
                        continue;

                    Console.WriteLine("Matching assembly: " + a.Location);
                }
                count++;

                MegaMurxPlusPlus(a);

                

                var tr = new TextRunner(a);
                int r = tr.Execute(args);

                Console.SetOut(new StreamWriter(Console.OpenStandardOutput()));
                Console.SetError(new StreamWriter(Console.OpenStandardError()));

                //var ar = new AutoRun(a);
                //int r = ar.Execute(args);

                Console.WriteLine("Nunit returend code " + r);

                ret = ret | (r != 0);
            }

            if (!AssemblyFilter.IsEmptyOrWhite()) {
                if(count <= 0) {
                    Console.WriteLine("Found no assembly matching: " + AssemblyFilter);
                }
            }


            return ret ? -1 : 0;
        }

        static int Main(string[] args) {
            args = BoSSS.Solution.Application.ArgsFromEnvironmentVars(args);

            var ll = System.Diagnostics.Debug.Listeners;
            ll.Clear();
            ll.Add(new MyListener());

            BoSSS.Solution.Application.InitMPI();


            int ret = -1; 
            switch (args[0]) {
                case "--nunit3":
                Console.WriteLine("Assembly filter: " + (args[1] != null ? args[1] : "NULL"));
                ret = RunSerial(args[1], args.Skip(2).ToArray());
                break;
                               

                case "--runjobmanager":
                DeleteResultFiles();
                ret = JobManagerRun(args.Length > 1 ? args[1] : null);
                break;

                default:
                throw new NotSupportedException("unknown subprogram.");
            }

            csMPI.Raw.mpiFinalize();
            csMPI.Raw.mpiFinalize();

            return ret;
            
        }
    }
}

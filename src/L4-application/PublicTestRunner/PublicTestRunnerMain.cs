using BoSSS.Application.BoSSSpad;
using ilPSP;
using NUnit.Framework;
using NUnit.Framework.Interfaces;
using NUnitLite;
using System;
using System.Collections.Generic;
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
            typeof(BoSSS.Application.EllipticReInitTest.EllipticReInitMain),
            typeof(BoSSS.Application.LevelSetTestBench.LevelSetTestBenchMain),
            typeof(BoSSS.Application.XdgPoisson3.XdgPoisson3Main),
            typeof(BoSSS.Application.AdaptiveMeshRefinementTest.AllUpTest),
            typeof(BoSSS.Application.ExternalBinding.CodeGen.Test),
            typeof(BoSSS.Application.ExternalBinding.Initializer)
        };

        static Type[] ReleaseOnlyTests = new Type[] {
            typeof(BoSSS.Application.TutorialTests.AllUpTest),
            typeof(CNS.Program)//,
            //typeof(LowMachCombustionNSE
            //TutorialTests.exe
            // QuadratureAndProjectionTest.exe XdgNastyLevsetLocationTest.exe LTSTests.exe XNSE_ViscosityAgglomerationTest.exe NSE_SIMPLE/bin/*/NSE_SIMPLE.exe EllipticReInit.exe IBM_Solver/bin/*/IBM_Solver.exe FSI_Solver/bin/*/FSI_Solver.exe ALTSTests.exe XNSE_Solver/bin/*/XNSE_Solver.exe XDGShock/bin/*/XDGShock.exe
        };


        static Assembly[] GetAllAssemblies() {
            var R = new HashSet<Assembly>();

            foreach(var t in FullTestTypes) {
                R.Add(t.Assembly);
            }
#if !DEBUG
            foreach (var t in FullTestTypes) {
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
                                string filepath = LocateFile(dc.SomeFileName);
                                s.Add(filepath);
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
            InteractiveShell.ReloadExecutionQueues();
            InteractiveShell.WorkflowMgm.Init("BoSSSTesting");

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

            List<Job> allJobs = new List<Job>();
            foreach(var t in allTests) {
                var j = JobManagerRun(t.ass, t.testname, bpc, t.depfiles);
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
                Console.WriteLine(j.ToString());
            }

            // ===================================
            // phase 3: collect files
            // ===================================

            int returnCode = 0;

            string CurrentDir = Path.GetDirectoryName(typeof(PublicTestRunnerMain).Assembly.Location);

            foreach (var j in allJobs) {
                //Console.WriteLine(j.ToString());

                if (j.Status != JobStatus.FinishedSuccessful)
                    returnCode--;

                try {
                    string[] sourceFiles = Directory.GetFiles(j.DeploymentDirectory, "result-*.xml");

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

            using (var ot = new StreamWriter("allout.txt")) {
                foreach (var j in allJobs) {
                    ot.WriteLine("##fdhgjegf763748trfhe8hurdsinf598ugf498jvhsn*hbbvc#####!################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("####  " + j.Name);
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

        static public Job JobManagerRun(Assembly a, string TestName, BatchProcessorClient bpc, string[] AdditionalFiles) {
            string dor = DebugOrReleaseSuffix;
            Job j = new Job($"test-{TestName}-{dor}", typeof(PublicTestRunnerMain));
            j.MySetCommandLineArguments("--nunit3", Path.GetFileName(a.Location), $"--test={TestName}", $"--result=result-{TestName}-{dor}.xml");

            foreach (var f in AdditionalFiles) {
                j.AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(File.ReadAllBytes(f), Path.GetFileName(f)));
            }
            j.Activate(bpc);
            return j;
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



            bool ret = false;
            foreach(var a in assln) {
                if(!AssemblyFilter.IsEmptyOrWhite()) {
                    if(!AssemblyFilter.WildcardMatch(Path.GetFileName(a.Location)))
                        continue;
                }

                MegaMurxPlusPlus(a);

                

                //var tr = new TextRunner(a);
                //int r = tr.Execute(args);

                var ar = new AutoRun(a);
                int r = ar.Execute(args);

                Console.WriteLine("Nunit returend code " + r);

                ret = ret | (r != 0);
            }

            
            return ret ? -1 : 0;
        }

        static int Main(string[] args) {
            args = BoSSS.Solution.Application.ArgsFromEnvironmentVars(args);

            var ll = System.Diagnostics.Debug.Listeners;
            ll.Clear();
            ll.Add(new MyListener());

            switch(args[0]) {
                case "--nunit3":
                Console.WriteLine("Assembly filter: " + (args[1] != null ? args[1] : "NULL"));
                return RunSerial(args[1], args.Skip(2).ToArray());
                               

                case "--runjobmanager":
                DeleteResultFiles();
                return JobManagerRun(args.Length > 1 ? args[1] : null);
            }

            throw new NotSupportedException("unknown subprogram.");
        }
    }
}

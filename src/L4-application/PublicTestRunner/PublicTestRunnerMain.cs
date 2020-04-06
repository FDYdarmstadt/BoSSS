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
            Assert.Fail(message);
        }

        public override void Fail(string message, string detailMessage) {
            Assert.Fail(message + ", details: " + detailMessage);
        }
    }


    static class PublicTestRunnerMain {

        /// <summary>
        /// List of tests that should be executed in DEBUG and RELEASE; referencing any type of the assembly will do.
        /// </summary>
        static Type[] FullTestTypes = new Type[] {
            //typeof(BoSSS.Application.DerivativeTest.DerivativeTestMain),
            //typeof(BoSSS.Application.SipPoisson.SipPoissonMain),
            //typeof(BoSSS.Application.Matrix_MPItest.AllUpTest),
            typeof(BoSSS.Application.ElementTests.ElementTests),
            //typeof(BoSSS.Application.DatabaseTests.DatabaseTestsProgram),
            typeof(CutCellQuadrature.Program),
            typeof(BoSSS.Application.XDGTest.UnitTest),
            //typeof(BoSSS.Application.SpecFEM.AllUpTest),
            //typeof(BoSSS.Application.ipViscosity.TestSolution),
            //typeof(BoSSS.Application.MultigridTest.MultigridMain),
            //typeof(BoSSS.Application.ZwoLsTest.AllUpTest),
            //typeof(BoSSS.Application.XdgTimesteppingTest.XdgTimesteppingMain),
            //typeof(BoSSS.Application.EllipticReInitTest.EllipticReInitMain),
            //typeof(BoSSS.Application.LevelSetTestBench.LevelSetTestBenchMain),
            //typeof(BoSSS.Application.XdgPoisson3.XdgPoisson3Main),
            //typeof(BoSSS.Application.AdaptiveMeshRefinementTest.AllUpTest),
            //typeof(BoSSS.Application.ExternalBinding.CodeGen.Test),
            typeof(BoSSS.Application.ExternalBinding.Initializer)
        };

        static Type[] ReleaseOnlyTests = new Type[] {
            typeof(BoSSS.Application.TutorialTests.AllUpTest)
        };


        static Assembly[] GetAllAssemblies() {
            var R = new HashSet<Assembly>();

            foreach(var t in FullTestTypes) {
                R.Add(t.Assembly);
            }

            return R.ToArray();
        }

        static string[] GetTestsInAssembly(Assembly a) {
            var r = new List<string>();

            var ttt = a.GetTypes();
            foreach(var t in ttt) {
                if(t.IsClass) {
                    var mmm = t.GetMethods();

                    foreach(var m in mmm) {
                        if(m.GetCustomAttribute(typeof(TestAttribute)) != null) {
                            r.Add(t.FullName + "." +  m.Name);
                        }
                    }
                }
            }

            return r.ToArray();
        }


        class MyTestFilter : NUnit.Framework.Internal.TestFilter {
            public override TNode AddToXml(TNode parentNode, bool recursive) {
                throw new NotImplementedException();
            }

            public override bool Match(ITest test) {
                throw new NotImplementedException();
            }
        }


        static public void JobManagerRun(string AssemblyFilter) {
            InteractiveShell.ReloadExecutionQueues();
            InteractiveShell.WorkflowMgm.Init("BoSSSTesting");

            BatchProcessorClient bpc = InteractiveShell.ExecutionQueues[1];

            var allTests = new List<ValueTuple<Assembly, string>>();
            var assln = GetAllAssemblies();
            foreach(var a in assln) {
                if(!AssemblyFilter.IsEmptyOrWhite()) {
                    if(!AssemblyFilter.WildcardMatch(Path.GetFileName(a.Location)))
                        continue;
                }


                string[] allTest = GetTestsInAssembly(a);

                foreach(var t in allTest) {
                    allTests.Add((a, t));
                }

            }

            List<Job> allJobs = new List<Job>();
            foreach(var t in allTests) {
                JobManagerRun(t.Item1, t.Item2, bpc);
            }

            while(InteractiveShell.WorkflowMgm.BlockUntilAnyJobTerminate(out var job, PollingIntervallSeconds: 120) > 0) {

                if(job != null) {
                    Console.WriteLine(job.Name + ": " + job.Status);
                }
            }
        }



        static public void JobManagerRun(Assembly a, string TestName, BatchProcessorClient bpc) {
            string dor;
#if DEBUG
            dor = "DEBUG";
#else
            dor = "RELEASE";
#endif
            Job j = new Job($"test-{TestName}-{dor}", typeof(PublicTestRunnerMain));


            j.MySetCommandLineArguments("--nunit3", Path.GetFileName(a.Location), $"--test={TestName}", $"--result=result-{TestName}-{dor}.xml");
            
            j.Activate(bpc);
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


     

                var tr = new TextRunner(a);
                int r = tr.Execute(args);
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
                JobManagerRun(args.Length > 1 ? args[1] : null);
                break;
            }



            return 0;
        }
    }
}

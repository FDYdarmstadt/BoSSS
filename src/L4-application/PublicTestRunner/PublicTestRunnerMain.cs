using BoSSS.Application.BoSSSpad;
using ilPSP;
using ilPSP.Utils;
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
        static Type[] FullTest = new Type[] {
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
            //typeof(BoSSS.Application.LevelSetTestBench.LevelSetTestBenchMain),
            typeof(BoSSS.Application.XdgPoisson3.XdgPoisson3Main),
            //typeof(BoSSS.Application.AdaptiveMeshRefinementTest.AllUpTest),
            typeof(BoSSS.Application.ExternalBinding.CodeGen.Test),
            typeof(BoSSS.Application.ExternalBinding.Initializer),
            //typeof(BoSSS.Application.TutorialTests.AllUpTest),
            typeof(MPITest.Program)
        };

        static Type[] ReleaseOnlyTests = new Type[] {
            //typeof(BoSSS.Application.TutorialTests.AllUpTest),
            typeof(BoSSS.Application.XdgTimesteppingTest.XdgTimesteppingMain),
            typeof(CNS.Program),
            typeof(QuadratureAndProjectionTest.QuadratueAndProjectionTest),
            typeof(BoSSS.Application.XdgNastyLevsetLocationTest.AllUpTest),
            typeof(LTSTests.Program),
            //typeof(BoSSS.Application.XNSE_ViscosityAgglomerationTest.XNSE_ViscosityAgglomerationTestMain),
            typeof(NSE_SIMPLE.SIMPLESolver),
            typeof(BoSSS.Application.IBM_Solver.IBM_SolverMain),
            typeof(ALTSTests.Program),
            typeof(BoSSS.Application.XNSE_Solver.XNSE_SolverMain)
        };

        static (Type type, int NoOfProcs)[] MpiFullTests = new (Type type, int NoOfProcs)[] {
            (typeof(MPITest.Program), 4),
            (typeof(MPITest.Program), 3),
            (typeof(MPITest.Program), 2),
            (typeof(BoSSS.Application.SpecFEM.AllUpTest), 4),
            (typeof(AdvancedSolverTests.SubBlocking.LocalTests),4),
            (typeof(AdvancedSolverTests.SubBlocking.ExternalTests),4)
        };

        static (Type type, int NoOfProcs)[] MpiReleaseOnlyTests = new (Type type, int NoOfProcs)[] {
            (typeof(MPITest.Program), 4),
            (typeof(BoSSS.Application.SpecFEM.AllUpTest), 4)
            //(typeof(BoSSS.Application.XNSE_Solver.XNSE_Solver_MPItest), 4),
            //(typeof(BoSSS.Application.Matrix_MPItest.AllUpTest), 4),
            //(typeof(BoSSS.Application.LoadBalancingTest.LoadBalancingTestMain), 4),
            //(typeof(ALTSTests.Program), 4),
            //(typeof(CNS_MPITests.Tests.LoadBalancing.ShockTubeLoadBalancingTests), 4),
            //(typeof(HilbertTest.HilbertTest), 4),
            //(typeof(BoSSS.Application.XdgPoisson3.XdgPoisson3Main), 4)
        };



        static Assembly[] GetAllAssemblies() {
            var R = new HashSet<Assembly>();

            foreach(var t in FullTest) {
                //Console.WriteLine("test type: " + t.FullName);
                var a = t.Assembly;
                //Console.WriteLine("  assembly: " + a.FullName + " @ " + a.Location);
                bool added = R.Add(a);
                //Console.WriteLine("  added? " + added);
            }
#if !DEBUG
            foreach (var t in ReleaseOnlyTests) {
                R.Add(t.Assembly);
            }
#endif

            

            return R.ToArray();
        }

        static (Assembly Asbly, int NoOfProcs)[] GetAllMpiAssemblies() {
            var R = new List<(Assembly Asbly, int NoOfProcs)>();

            foreach (var t in MpiFullTests) {
                //Console.WriteLine("test type: " + t.type.FullName + " (" + t.NoOfProcs + " procs).");
                //Console.WriteLine("  assembly: " + t.type.Assembly.FullName + " @ " + t.type.Assembly.Location);
                bool contains = R.Contains(t, (itm1, itm2) => ((itm1.NoOfProcs == itm2.NoOfProcs) && itm1.Asbly.Equals(itm2.type.Assembly)));
                if(!contains) {
                    R.Add((t.type.Assembly, t.NoOfProcs));
                }
                //Console.WriteLine("  added? " + (!contains));
            }
#if !DEBUG
            foreach (var t in MpiReleaseOnlyTests) {
                bool contains = R.Contains(t, (itm1, itm2) => ((itm1.NoOfProcs == itm2.NoOfProcs) && itm1.Asbly.Equals(itm2.type.Assembly)));
                if (!contains) {
                    R.Add((t.type.Assembly, t.NoOfProcs));
                }
            }
#endif
            return R.ToArray();
        }

        static string[] LocateFile(string PartialPath) {
            DirectoryInfo repoRoot;
            try {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                repoRoot = dir.Parent.Parent.Parent.Parent.Parent;

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
                return null;
                //throw new IOException("Unable to find repository root. 'runjobmanger' must be invoked from its default location within the BoSSS git repository.");
            }

            // if we get here, we probably have access to the repository root directory.
            string[] r = LocateFileRecursive("", repoRoot, PartialPath);
            if(r == null || r.Length <= 0) {
                throw new IOException("unable to find file '" + PartialPath  + "'"); 
            }
            
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

            foreach(var d in absPath.GetDirectories()) {
                ret.AddRange(LocateFileRecursive(RelPath + d.Name + "/", d, SomeFileName));
            }


            return ret.ToArray();
        }



        static (int NoOfTests, string[] tests, string[] shortnames, string[] RequiredFiles) GetTestsInAssembly(Assembly a) {
            var r = new List<string>(); // full test names
            var l = new List<string>(); // short names 
            var s = new HashSet<string>();

            var ttt = a.GetTypes();
            foreach(var t in ttt) {
                if(t.IsClass) {
                    var mmm = t.GetMethods();

                    foreach(var m in mmm) {
                        if (t.IsAbstract && !m.IsStatic)
                            continue;
                        
                        if(m.GetCustomAttribute(typeof(TestAttribute)) != null) {
                            r.Add(t.FullName + "." +  m.Name);
                            l.Add(Path.GetFileNameWithoutExtension(a.ManifestModule.Name) + "#" + m.Name);
                        }

                        if(m.GetCustomAttribute(typeof(TestAttribute)) != null
                           || m.GetCustomAttribute(typeof(SetUpAttribute)) != null
                           || m.GetCustomAttribute(typeof(OneTimeSetUpAttribute)) != null) {
                            var dc = m.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) as NUnitFileToCopyHackAttribute;

                            if(dc != null) {
                                foreach(string someFile in dc.SomeFileNames) {
                                    s.AddRange(LocateFile(someFile));
                                }
                            }
                        }
                    }
                }
            }

            return (r.Count, r.ToArray(), l.ToArray(), s.ToArray());
        }


        class MyTestFilter : NUnit.Framework.Internal.TestFilter {
            public override TNode AddToXml(TNode parentNode, bool recursive) {
                throw new NotImplementedException();
            }

            public override bool Match(ITest test) {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// (tries to) do a recursive copy of a directory
        /// </summary>
        static void CopyDirectoryRec(string srcDir, string dstDir) {

            void TryCopy(string sourceFileName, string destFileName) {
                try {
                    File.Copy(sourceFileName, destFileName, true);
                } catch (Exception e) {
                    Console.WriteLine("WARNING: Unable to copy to: '"
                        + destFileName + "': " + e.GetType().Name + " says:'" + e.Message + "'");
                }
            }

            string[] srcFiles = Directory.GetFiles(srcDir);

            foreach (string srcFile in srcFiles) {
                TryCopy(srcFile, Path.Combine(dstDir, Path.GetFileName(srcFile)));
            }

            string[] subDirs = Directory.GetDirectories(srcDir);
            foreach (string srcAbsDir in subDirs) {
                string srcRelDir = Path.GetFileName(srcAbsDir);
                string dstSubDir = Path.Combine(dstDir, srcRelDir);
                if (!Directory.Exists(dstSubDir))
                    Directory.CreateDirectory(dstSubDir);
                CopyDirectoryRec(srcAbsDir, dstSubDir);
            }
        }


        static public int JobManagerRun(string AssemblyFilter) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out var MpiSize);
            if(MpiSize != 1) {
                throw new NotSupportedException("runjobmanager subprogram must be executed serially"); 
            }

            // ===================================
            // phase 1: submit jobs
            // ===================================

            string DateNtime = DateTime.Now.ToString("MMMdd_HHmm");
            Console.WriteLine($"Using prefix'{DateNtime}' for all jobs.");

            InteractiveShell.ReloadExecutionQueues();
            InteractiveShell.WorkflowMgm.Init("BoSSStst" + DateNtime);

            BatchProcessorClient bpc = InteractiveShell.ExecutionQueues[1];

            DirectoryInfo NativeOverride;
            if (!bpc.DeployRuntime) {
                NativeOverride = new DirectoryInfo(Path.Combine(bpc.DeploymentBaseDirectory, DateNtime + "_amd64"));
                NativeOverride.Create();
                CopyDirectoryRec(ilPSP.Environment.NativeLibraryDir, NativeOverride.FullName);
            } else {
                NativeOverride = null;
            }

            var allTests = new List<(Assembly ass, string testname, string shortname, string[] depfiles, int NoOfProcs)>();
            {
                var assln = GetAllAssemblies();
                foreach (var a in assln) {
                    if (!AssemblyFilter.IsEmptyOrWhite()) {
                        if (!AssemblyFilter.WildcardMatch(Path.GetFileName(a.Location)))
                            continue;
                    }
                    
                    var allTst4Assi = GetTestsInAssembly(a);

                    for (int iTest = 0; iTest < allTst4Assi.NoOfTests; iTest++) { 
                        allTests.Add((a, allTst4Assi.tests[iTest], allTst4Assi.shortnames[iTest], allTst4Assi.RequiredFiles, 1));
                    }
                }
            }
            {
                var ParAssln = GetAllMpiAssemblies();
                foreach (var TT in ParAssln) {
                    if (!AssemblyFilter.IsEmptyOrWhite()) {
                        if (!AssemblyFilter.WildcardMatch(Path.GetFileName(TT.Asbly.Location)))
                            continue;
                    }

                    var a = TT.Asbly;
                    var allTst4Assi = GetTestsInAssembly(a);

                    for (int iTest = 0; iTest < allTst4Assi.NoOfTests; iTest++) {
                        allTests.Add((a, allTst4Assi.tests[iTest], allTst4Assi.shortnames[iTest], allTst4Assi.RequiredFiles, TT.NoOfProcs));
                    }
                }
            }

            Console.WriteLine($"Found {allTests.Count} individual tests ({DebugOrReleaseSuffix}):");
            int cnt = 0;
            foreach (var t in allTests) {
                cnt++;
                Console.WriteLine($"  #{cnt}: {t.testname}");
                Console.WriteLine($"     {t.shortname}");
                Console.WriteLine($"     {t.NoOfProcs} MPI processors.");
            }

            Console.WriteLine("******* Starting job deployment/submission *******");

            DateTime start = DateTime.Now;

            cnt = 0;
            var allJobs = new List<(Job job, string ResFile, string testname)>();
            foreach(var t in allTests) {
                cnt++;
                Console.WriteLine($"Submitting {cnt} of {allTests.Count} ({t.shortname})...");
                var j = JobManagerRun(t.ass, t.testname, t.shortname, bpc, t.depfiles, DateNtime, t.NoOfProcs, NativeOverride);
                Console.WriteLine($"Successfully submitted {j.j.Name}.");
                allJobs.Add(j);
            }

            Console.WriteLine("******* All jobs deployed *******");

            // ===================================
            // phase 2: wait until complete...
            // ===================================

            const double TimeOutSec = 200*60;

            var alreadyFinished = InteractiveShell.WorkflowMgm.AllJobs.Select(kv => kv.Value).Where(delegate (Job j) {
                var s = j.Status;
                if (s == JobStatus.Failed)
                    return true;
                if (s == JobStatus.FinishedSuccessful)
                    return true;
                return false;
            }).ToArray();

            foreach ( var t in alreadyFinished) {
                Console.WriteLine("already finished: " + t.Name + ": " + t.Status);
            }

            double RestTime = Math.Max(1, TimeOutSec - (DateTime.Now - start).TotalSeconds);

            while (InteractiveShell.WorkflowMgm.BlockUntilAnyJobTerminate(out var job, PollingIntervallSeconds: 120, TimeOutSeconds: RestTime) > 0) {

                if(job != null) {
                    Console.WriteLine("just finished: " + job.Name + ": " + job.Status);
                }

                RestTime = TimeOutSec - (DateTime.Now - start).TotalSeconds;
                if(RestTime <= 1.0) {
                    Console.Error.WriteLine("timeout.");
                    break;
                }
                RestTime = Math.Max(1, RestTime);
            }
            Thread.Sleep(10000);
            Console.WriteLine("----------------------------------");
            Console.WriteLine("All jobs finished - Summary:");
            Console.WriteLine("----------------------------------");
            int SuccessfulFinishedCount = 0;
            foreach (var j in allJobs) {
                if (j.job.Status == JobStatus.FinishedSuccessful)
                    SuccessfulFinishedCount++;
                Console.WriteLine(j.testname + ": " + j.job.ToString());
            }
            if(SuccessfulFinishedCount == allJobs.Count) {
                Console.WriteLine("All tests/jobs finished successfully.");
            } else {
                Console.WriteLine($"Only {SuccessfulFinishedCount} tests/jobs finished successfully -- {allJobs.Count - SuccessfulFinishedCount} have other states.");

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



            using (var ot = new StreamWriter("allout-" + DateNtime + "-" + DebugOrReleaseSuffix + ".txt")) {

                var allNotSuccessful = allJobs.Where(jj => jj.job.Status != JobStatus.FinishedSuccessful).ToArray();
                var allSuccessful = allJobs.Where(jj => jj.job.Status == JobStatus.FinishedSuccessful).ToArray();

                if (allNotSuccessful.Length > 0) {
                    ot.WriteLine("All jobs not finished successful:");
                    foreach (var jj in allJobs) {
                        Console.WriteLine($"{jj.job.Status}: {jj.job.Name} ({jj.testname}, at {jj.job.DeploymentDirectory})");
                    }
                }

                var _allJobs = ArrayTools.Cat(allNotSuccessful,allSuccessful);
                foreach (var jj in _allJobs) {
                    var j = jj.job;
                    int sz = j.NumberOfMPIProcs;
                    ot.WriteLine("##fdhgjegf763748trfhe8hurdsinf598ugf498jvhsn*hbbvc#####!################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("####  " + jj.testname);
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("#### Deploy directory: " + j.DeploymentDirectory);
                    ot.WriteLine("#### Full test name:   " + j.Name);
                    ot.WriteLine("#### Number of procs:  " + j.NumberOfMPIProcs);
                    ot.WriteLine("#### Status:           " + j.Status);
                    ot.WriteLine("#### Job ID:           " + j.BatchProcessorIdentifierToken);
                    //                                   +   +   +   +
                    if (j.NumberOfMPIProcs <= 1) {
                        ot.WriteLine("#### Result File:      " + jj.ResFile);
                        ot.WriteLine("####    exists?        " + File.Exists(Path.Combine(j.DeploymentDirectory, jj.ResFile)));
                    } else {
                        for (int i = 0; i < sz; i++) {
                            if (i == 0)
                                ot.WriteLine("#### Result Files:     " + MpiResFileNameMod(i, sz, jj.ResFile));
                            else
                                ot.WriteLine("####                   " + MpiResFileNameMod(i, sz, jj.ResFile));
                        }
                        for (int i = 0; i < sz; i++) {
                            if (i == 0)
                                ot.WriteLine("####    exists?        " + File.Exists(Path.Combine(j.DeploymentDirectory, MpiResFileNameMod(i, sz, jj.ResFile))));
                            else
                                ot.WriteLine("####                   " + File.Exists(Path.Combine(j.DeploymentDirectory, MpiResFileNameMod(i, sz, jj.ResFile))));
                        }
                    }
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");
                    ot.WriteLine("########################################################################");

                    ot.WriteLine("[[[Stdout: ");
                    try {
                        ot.WriteLine(j.Stdout);
                    } catch (Exception e) {
                        ot.WriteLine($"{e.GetType().Name} during reading of stdout stream: {e.Message}");
                        ot.WriteLine(e.StackTrace);
                    }
                    ot.WriteLine("]]]");

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

                    if(jj.job.Status == JobStatus.FinishedSuccessful) {
                        try {
                            Directory.Delete(jj.job.DeploymentDirectory, true);
                        } catch(Exception e) {
                            Console.Error.WriteLine($"{e.GetType().Name}: {e.Message}");
                        }
                    }
                }
            }

            return returnCode;
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

        static public (Job j, string resultFile, string name) JobManagerRun(
            Assembly a,
            string TestName, string Shortname, 
            BatchProcessorClient bpc, 
            string[] AdditionalFiles, 
            string prefix, 
            int NoOfMpiProcs,
            DirectoryInfo nativeOverride) {

            // define unique name (not to long) for the job
            string dor = DebugOrReleaseSuffix;
            string final_jName;
            {
                string jName;
                if (NoOfMpiProcs <= 1)
                    jName = $"{prefix}-{dor}-{Shortname}";
                else
                    jName = $"{prefix}p{NoOfMpiProcs}-{dor}-{Shortname}";

                if (jName.Length > 127) {
                    // Name length limit set by MS HPC cluster
                    jName = jName.Substring(0, 127);
                }
                int counter = 2;
                final_jName = jName;
                while (InteractiveShell.WorkflowMgm.AllJobs.ContainsKey(final_jName)) {
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
            Job j = new Job(final_jName, typeof(PublicTestRunnerMain));
            string resultFile = $"result-{dor}-{TestName}.xml";
            j.MySetCommandLineArguments("nunit3", Path.GetFileName(a.Location), $"--test={TestName}", $"--result={resultFile}");
            foreach (var f in AdditionalFiles) {
                j.AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(File.ReadAllBytes(f), Path.GetFileName(f)));
            }
            if(nativeOverride != null) {
                j.EnvironmentVars.Add(BoSSS.Foundation.IO.Utils.BOSSS_NATIVE_OVERRIDE, nativeOverride.FullName);
            }
            j.NumberOfMPIProcs = NoOfMpiProcs;
            j.Activate(bpc);
            return (j, resultFile, TestName);
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
        static int RunNunit3Tests(string AssemblyFilter, string[] args) {
            var assln = GetAllAssemblies();

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out var MpiSize);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out var MpiRank);

            if(MpiSize != 1) {
                // this seems some parallel run
                // we have to fix the result argument

                if (args.Where(a => a.StartsWith("--result=")).Count() != 1) {
                    throw new ArgumentException("MPI-parallel NUnit runs require the '--result' argument.");
                }

                int i = args.IndexWhere(a => a.StartsWith("--result="));
                string arg_i = args[i];
                string resFileName = Path.GetFileNameWithoutExtension(arg_i.Replace("--result=", ""));
                args[i] = "--result=" + MpiResFileNameMod( MpiRank, MpiSize, resFileName);

            }


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

                Console.WriteLine("NUnit returned code " + r);

                ret = ret | (r != 0);
            }

            if (!AssemblyFilter.IsEmptyOrWhite()) {
                if(count <= 0) {
                    Console.WriteLine("Found no assembly matching: " + AssemblyFilter);
                }
            }


            return ret ? -1 : 0;
        }

        private static string MpiResFileNameMod( int MpiRank, int MpiSize, string resFileName) {
            if (MpiRank > MpiSize)
                throw new ArgumentException(); 

            if (MpiSize == 1)
                return resFileName;

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
            Console.WriteLine("         runjobmanager : submit tests to the job manger.");
            Console.WriteLine("  and FILTER selects some assembly, i.e. DerivativeTests.exe; it can be ");
            Console.WriteLine("  a wildcard, i.e. use * for submitting all tests.");

        }


        static int Main(string[] args) {
            //Debugger.Launch();
            Console.WriteLine("BoSSS NUnit test runner.");

            args = BoSSS.Solution.Application.ArgsFromEnvironmentVars(args);

            var ll = System.Diagnostics.Debug.Listeners;
            ll.Clear();
            ll.Add(new MyListener());


            if(args.Length < 2) {
                Console.WriteLine("Insufficient number of arguments.");
                PrintMainUsage();
                return -7777;
            }

            BoSSS.Solution.Application.InitMPI();


            int ret = -1; 
            switch (args[0]) {
                case "nunit3":
                ret = RunNunit3Tests(args[1], args.Skip(2).ToArray());
                break;
                               

                case "runjobmanager":
                DeleteResultFiles();
                ret = JobManagerRun(args[1]);
                break;

                default:
                throw new NotSupportedException("unknown subprogram.");
            }

            csMPI.Raw.mpiFinalize();
            

            return ret;
            
        }
    }
}

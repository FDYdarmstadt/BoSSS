using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using BoSSS.Solution.Control;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using ilPSP;
using System.Diagnostics;
using ilPSP.Tracing;
using System.IO;
using System.Xml.Linq;
using System.Reflection;
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.BoSSSpad {
    /// <summary>
    /// Extensions for <see cref="BoSSS.Solution.Control.AppControl"/> objects.
    /// </summary>
    public static class AppControlExtensions {

        /// <summary>
        /// Runs the solver described by the control object <paramref name="ctrl"/> in the current process.
        /// The method blocks until the solver is finished.
        /// </summary>
        /// <param name="ctrl"></param>
        /// <returns>
        /// The Session information after the solver is finished.
        /// </returns>
        public static SessionInfo Run(this AppControl ctrl) {

            var solverClass = ctrl.GetSolverType();
            object solver = Activator.CreateInstance(solverClass);

            var app = (BoSSS.Solution.IApplication)solver;

            app.Init(ctrl);
            app.RunSolverMode();

            var S = app.CurrentSessionInfo;
            
            if(solver is IDisposable)
                ((IDisposable)solver).Dispose();

            return S;
        }

        /// <summary>
        /// Creates a job for the control object <paramref name="ctrl"/>.
        /// The method returns immediately.
        /// This job can still be configured (e.g. setting number of MPI processors) and must be activated (<see cref="Job.Activate(BatchProcessorClient)"/>)
        /// to run on a batch system.
        /// </summary>
        /// <param name="ctrl"></param>
        /// <returns></returns>
        public static Job CreateJob(this AppControl ctrl) {
            ctrl.ProjectName = BoSSSshell.WorkflowMgm.CurrentProject;

            string JobName = ctrl.SessionName;
            int ctrl_idx = BoSSSshell.WorkflowMgm.RegisterControl(ctrl);
            if (JobName.IsEmptyOrWhite()) {
                JobName = "UnnamedJob_" + ctrl_idx;
                ctrl.SessionName = JobName;
            }

            Type solverClass = ctrl.GetSolverType();
            Job job = new Job(ctrl.SessionName, solverClass);
            job.SetControlObject(ctrl);


            return job;
        }

        /// <summary>
        /// Runs the solver described by the control object <paramref name="ctrl"/> on a batch system.
        /// The method returns immediately.
        /// </summary>
        /// <param name="ctrl"></param>
        /// <param name="BatchSys"></param>
        /// <returns></returns>
        public static Job RunBatch(this AppControl ctrl, BatchProcessorClient BatchSys) {
            var job = ctrl.CreateJob();
            
            job.Activate(BatchSys);
            
            return job;
        }

        /// <summary>
        /// Runs the solver described by the control object <paramref name="ctrl"/> on a batch system from the currently defined queues (<see cref="BoSSSshell.ExecutionQueues"/>).
        /// The method returns immediately.
        /// </summary>
        /// <param name="ctrl"></param>
        /// <param name="queueIdx">
        /// Index int <see cref="BoSSSshell.ExecutionQueues"/>
        /// </param>
        public static Job RunBatch(this AppControl ctrl, int queueIdx) {
            var b = BoSSSshell.ExecutionQueues[queueIdx];
            return RunBatch(ctrl, b);
        }


        /// <summary>
        /// Runs the solver described by the control object <paramref name="ctrl"/> 
        /// on the default batch system 
        /// The method returns immediately after the job is deployed., i.e. it does not wait for the job to finish.
        /// </summary>
        public static Job RunBatch(this AppControl ctrl) {
            var b = BoSSSshell.GetDefaultQueue();
            return RunBatch(ctrl, b);
        }

        /// <summary>
        /// Runs the solver described by the control objects <paramref name="ctrls"/> on a batch system.
        /// The method returns immediately.
        /// </summary>
        /// <param name="ctrls"></param>
        /// <param name="BatchSys"></param>
        /// <param name="DeployAssembliesOnce">If true deploy all necessary native and managed libraries only once</param>
        /// <returns></returns>
        public static IEnumerable<Job> RunBatch(this IEnumerable<AppControl> ctrls, BatchProcessorClient BatchSys, bool DeployAssembliesOnce = true) {
            var jobs = ctrls.Select(ctrl => ctrl.CreateJob()).ToList();

            jobs.Activate(BatchSys, DeployAssembliesOnce);

            return jobs;
        }

        /// <summary>
        /// Runs the solver described by the control objects <paramref name="ctrls"/> on a batch system from the currently defined queues (<see cref="BoSSSshell.ExecutionQueues"/>).
        /// The method returns immediately.
        /// </summary>
        /// <param name="ctrls"></param>
        /// <param name="queueIdx">
        /// Index int <see cref="BoSSSshell.ExecutionQueues"/>
        /// </param>
        /// <param name="DeployAssembliesOnce">If true deploy all necessary native and managed libraries only once</param>
        public static IEnumerable<Job> RunBatch(this IEnumerable<AppControl> ctrls, int queueIdx, bool DeployAssembliesOnce = true) {
            var b = BoSSSshell.ExecutionQueues[queueIdx];
            return RunBatch(ctrls, b, DeployAssembliesOnce);
        }


        /// <summary>
        /// Runs the solver described by the control objects <paramref name="ctrls"/> 
        /// on the default batch system 
        /// The method returns immediately after the job is deployed., i.e. it does not wait for the job to finish.
        /// </summary>
        public static IEnumerable<Job> RunBatch(this IEnumerable<AppControl> ctrls, bool DeployAssembliesOnce = true) {
            var b = BoSSSshell.GetDefaultQueue();
            return RunBatch(ctrls, b, DeployAssembliesOnce);
        }

        /// <summary>
        /// Same functionality as <see cref="Job.Activate()"/>, but for a collection of jobs.
        /// </summary>
        /// <param name="jobs"></param>
        /// <param name="DeployAssembliesOnce">If true deploy all necessary native and managed libraries only once</param>
        public static void Activate(this IEnumerable<Job> jobs, bool DeployAssembliesOnce = true) {
            var b = BoSSSshell.GetDefaultQueue();
            Activate(jobs, b, DeployAssembliesOnce);
        }

        /// <summary>
        /// Same functionality as <see cref="Job.Activate()"/>, but for a collection of jobs.
        /// </summary>
        /// <param name="jobs"></param>
        /// <param name="queueIdx"></param>
        /// <param name="DeployAssembliesOnce">If true deploy all necessary native and managed libraries only once</param>
        public static void Activate(this IEnumerable<Job> jobs, int queueIdx, bool DeployAssembliesOnce = true) {
            var b = BoSSSshell.ExecutionQueues[queueIdx];
            Activate(jobs, b, DeployAssembliesOnce);
        }

        /// <summary>
        /// Same functionality as <see cref="Job.Activate()"/>, but for a collection of jobs.
        /// If activated, native and managed assemblies are copied only once and the directory overrides set to a relative location of the job deployment directory.
        /// </summary>
        /// <param name="jobs"></param>
        /// <param name="BatchSys"></param>
        /// <param name="DeployAssembliesOnce">If true deploy all necessary native and managed libraries only once</param>
        public static void Activate(this IEnumerable<Job> jobs, BatchProcessorClient BatchSys, bool DeployAssembliesOnce = true) {

            // business as usual
            if (!DeployAssembliesOnce) {
                jobs.ForEach(job => job.Activate(BatchSys));
                return;
            }

            Console.WriteLine(" Using the DeployAssembliesOnce option, this is experimental and untested if all necessary files are copied in all cases!");

            // some checks
            Assembly EntryAssembly = jobs.First().EntryAssembly;
            foreach (var job in jobs) {
                if (job.EntryAssembly.FullName != EntryAssembly.FullName) 
                    throw new ApplicationException("All jobs must use the same solver!");
            }


            string JobDirectoryBaseName() {
                string Exe = Path.GetFileNameWithoutExtension(EntryAssembly.Location);
                string Proj = BoSSSshell.WorkflowMgm.CurrentProject;

                return Proj
                    //+ "-" + Sess 
                    + "-" + Exe + "-binaries-";
            }

            string GetNewDeploymentDir() {
                if (BatchSys == null)
                    throw new NotSupportedException("Job is not activated yet.");

                if (!Path.IsPathRooted(BatchSys.DeploymentBaseDirectory))
                    throw new IOException($"Deployment base directory for {BatchSys.ToString()} must be rooted/absolute, but '{BatchSys.DeploymentBaseDirectory}' is not.");

                string ShortName = JobDirectoryBaseName();
                string DeployDir;
                int Counter = 0;
                do {
                    string Suffix = Counter > 0 ? "-" + Counter : "";
                    string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss");
                    DeployDir = Path.Combine(BatchSys.DeploymentBaseDirectory, ShortName + DateNtime + Suffix);
                    Counter++;
                } while (Directory.Exists(DeployDir) == true);

                return DeployDir;
            }

            // set directory for managed and native assemblies
            string DeployDirPath = GetNewDeploymentDir();
            DirectoryInfo  DeployDir = new DirectoryInfo(DeployDirPath);
            DeployDir.Create();

            // Deploy managed libs ...
            Console.WriteLine("Deploying executables and additional files ... once");
            EntryAssembly.DeployAt(DeployDir);


            // Deploy native libs ...
            if (BatchSys.DeployRuntime == true) {
                string BosssInstall = BoSSS.Foundation.IO.Utils.GetBoSSSInstallDir();
                var BosssBinNative = new DirectoryInfo(Path.Combine(BosssInstall, "bin", "native", BatchSys.RuntimeLocation));
                MetaJobMgrIO.CopyDirectoryRec(BosssBinNative.FullName, DeployDirPath, null);
                Console.WriteLine("   copied '" + BatchSys.RuntimeLocation + "' runtime.");
                // ... but only once
                BatchSys.DeployRuntime = false;
            }

            foreach(var job in jobs) {
                job.EnvironmentVars.Add(BoSSS.Foundation.IO.Utils.BOSSS_NATIVE_OVERRIDE, Path.Combine("..", DeployDir.Name));
                job.EntryAssemblyRedirection = Path.Combine("..", DeployDir.Name, Path.GetFileName(EntryAssembly.Location));
                job.Activate(BatchSys);
            }

        }

        /// <summary>
        /// Returns the job correlated to a control object
        /// </summary>
        public static Job GetJob(this AppControl ctrl) {
            var ret = new List<Job>();
            foreach (var j in BoSSSshell.WorkflowMgm.AllJobs.Values) {
                var cj = j.GetControl();
                if (cj == null)
                    continue;

                if (cj.Equals(ctrl))
                    ret.Add(j);
            }
            if (ret.Count <= 0) {
                Console.WriteLine("No Job assigned for given control object yet.");
                return null;
            } else {
                if(ret.Count > 1) {
                    string messeage = $"Unable to find a 1:1 correlation between control object and jobs: matching jobs {ret.ToConcatString("", ", ", "")};";
                    throw new ApplicationException(messeage);
                }

                return ret[0];
            }
        }
        
        /// <summary>
        /// Returns the job correlated to a session
        /// </summary>
        public static Job GetJob(this ISessionInfo s) {
            var ctrl = s.GetControl();
            return ctrl.GetJob();
        }

        /// <summary>
        /// Returns all sessions which can be correlated to a specific control object
        /// </summary>
        public static ISessionInfo[] GetAllSessions(this AppControl ctrl) {
            var AllCandidates = BoSSSshell.WorkflowMgm.Sessions.Where(
                    sinf => BoSSSshell.WorkflowMgm.SessionInfoAppControlCorrelation(sinf, ctrl));

            var cnt = AllCandidates.Count();

            if(cnt <= 0)
                return new ISessionInfo[0];

            return AllCandidates.ToArray();
        }


        /// <summary>
        /// Creates a restart-job from the latest session (<see cref="Job.LatestSession"/>) in a given job
        /// </summary>
        public static Job CreateRestartJob(this Job job2rest) {
            var stat = job2rest.Status;
            if (stat == JobStatus.InProgress)
                throw new NotSupportedException("Unable to restart: job is currently in progress");
            if (stat == JobStatus.PreActivation)
                throw new NotSupportedException("Unable to restart: job is not activated yet");
            if (stat == JobStatus.InProgress)
                throw new NotSupportedException("Unable to restart: job is currently in progress");

            var job = job2rest.LatestSession.CreateRestartJob();
            job.NumberOfMPIProcs = job2rest.NumberOfMPIProcs;
            job.RetryCount = job2rest.RetryCount;
            //job.MemPerCPU = job2rest.MemPerCPU;
            //job.ExecutionTime = job2rest.ExecutionTime;
            return job;
        }

        /// <summary>
        /// Creates a restart-job from the latest session (<see cref="Job.LatestSession"/>) in a given job
        /// </summary>
        public static Job Restart(this Job job2rest) {
            var newJob = job2rest.CreateRestartJob();
            newJob.Activate(job2rest.AssignedBatchProc);
            return newJob;
        }

        /// <summary>
        /// Creates a restart-job control object from a given session.
        /// </summary>
        /// <param name="sess"></param>
        /// <returns></returns>
        public static AppControl CreateRestartControl(this ISessionInfo sess) {
            var ctrl = sess.GetControl();

            Guid rst_ID = sess.ID;
            TimestepNumber rst_ts = sess.Timesteps.Last().TimeStepNumber;

            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();

            ctrl.RestartInfo = Tuple.Create(rst_ID, rst_ts);

            ctrl.SessionName = ctrl.SessionName + "_Restart";

            return ctrl;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sess"></param>
        /// <returns></returns>
        public static Job CreateRestartJob(this ISessionInfo sess) {
            var ctrl = sess.CreateRestartControl();
            var job = ctrl.CreateJob();
            return job;
        }


        
        /// <summary>
        /// Verifies a control object, especially if it is suitable for serialization.
        /// </summary>
        /// <param name="ctrl"></param>
        public static void VerifyEx(this AppControl ctrl) {
           
            // call basic verification
            ctrl.Verify();



            // see is legacy-features are used, which don't support serialization.
            if (ctrl.GridFunc != null)
                throw new ArgumentException("'GridFunc' is not supported - cannot be serialized.");
            //if(ctrl.DynamicLoadBalancing_CellCostEstimators.Count != 0)
            //    throw new ArgumentException("'DynamicLoadBalancing_CellCostEstimatorFactories' is not supported - cannot be serialized.");

            // try serialization/deserialization
            AppControl ctrlBack;
            if (ctrl.GeneratedFromCode) {
                string code = ctrl.ControlFileText;
                
                AppControl.FromCode(code, ctrl.GetType(), out AppControl c, out AppControl[] cS);
                if(cS != null) {
                    ctrlBack = cS[ctrl.ControlFileText_Index];
                } else {
                    ctrlBack = c;
                }
            } else {
                string JSON = ctrl.Serialize();
                ctrlBack = AppControl.Deserialize(JSON);//, ctrl.GetType());
            }
            ctrlBack.Verify();

            // compare original and de-serialized object
            if (!ctrl.Tags.SetEquals(ctrlBack.Tags))
                throw new ArgumentException("Unable to serialize/deserialize tags correctly.");

            if(!ctrl.InitialValues.Keys.SetEquals(ctrlBack.InitialValues.Keys))
                throw new ArgumentException("Unable to serialize/deserialize initial values correctly.");

            foreach(string ivk in ctrl.InitialValues.Keys) {
                var f1 = ctrl.InitialValues[ivk];
                var f2 = ctrlBack.InitialValues[ivk];

                if(!f1.Equals(f2))
                    throw new ArgumentException($"Unable to serialize/deserialize initial values for '{ivk}' correctly. Note that you cannot use delegates with the workflow management, use either formulas as text or 'GetFormulaObject(...)' in the BoSSSpad -- or use other implementations of IBoundaryAndInitialData.");
            }

            if (!ctrl.FieldOptions.Keys.SetEquals(ctrlBack.FieldOptions.Keys))
                throw new ArgumentException("Unable to serialize/deserialize field options correctly.");
            foreach(var fok in ctrl.FieldOptions.Keys) {
                var o1 = ctrl.FieldOptions[fok];
                var o2 = ctrlBack.FieldOptions[fok];

                if(!o1.Equals(o2))
                    throw new ArgumentException("Unable to serialize/deserialize field options correctly.");
            }

            if (!ctrl.InitialValues_Evaluators.Keys.SetEquals(ctrlBack.InitialValues_Evaluators.Keys))
                throw new ArgumentException("Unable to serialize/deserialize initial values correctly.");

            string errBnyd(string n) {
                return $"Unable to serialize/deserialize boundary values for '{n}' correctly. Note that you cannot use delegates with the workflow management, use either formulas as text or 'GetFormulaObject(...)' in the BoSSSpad.";
            }

            if(!ctrl.BoundaryValues.Keys.SetEquals(ctrlBack.BoundaryValues.Keys))
                throw new ArgumentException("Unable to serialize/deserialize boundary values correctly.");

            foreach(var bvk in ctrl.BoundaryValues.Keys) {
                var bvc = ctrl.BoundaryValues[bvk];
                var bvd = ctrlBack.BoundaryValues[bvk];

                if(!bvc.Evaluators.Keys.SetEquals(bvd.Evaluators.Keys))
                    throw new ArgumentException(errBnyd(bvk));

                if(!bvc.Values.Keys.SetEquals(bvd.Values.Keys))
                    throw new ArgumentException(errBnyd(bvk));

                foreach(string s in bvc.Values.Keys) {
                    var f1 = bvc.Values[s];
                    var f2 = bvd.Values[s];

                    if(!f1.Equals(f2))
                        throw new ArgumentException(errBnyd(bvk + "," + s));
                }
            }

            // reflection comparison
            var D1 = new Dictionary<string, object>();
            var D2 = new Dictionary<string, object>();
            BoSSS.Solution.Application.FindKeys(D1, ctrl);
            BoSSS.Solution.Application.FindKeys(D2, ctrlBack);
            foreach(var kv1 in D1) {
                string name = kv1.Key;
                var o1 = kv1.Value;

                if((!D2.ContainsKey(name)) || (!KeyComparison(o1, D2[name]))) {
                    throw new ArgumentException("Unable to serialize/deserialize '" + name + "' correctly. (Missing a DataMemberAtribute in control class?)");
                }

            }


            bool isInt(object a) {

                switch(Type.GetTypeCode(a.GetType())) {
                    case TypeCode.Byte:
                    case TypeCode.SByte:
                    case TypeCode.UInt16:
                    case TypeCode.UInt32:
                    case TypeCode.UInt64:
                    case TypeCode.Int16:
                    case TypeCode.Int32:
                    case TypeCode.Int64:
                    return true;

                    default:
                    return false;
                }
            }

            bool isFloat(object a) {

                switch(Type.GetTypeCode(a.GetType())) {
                    case TypeCode.Double:
                    case TypeCode.Single:
                    return true;

                    default:
                    return false;
                }
            }

            bool KeyComparison(object a, object b) {
                if(a == null && b == null)
                    return true;

                if(a == null || b == null)
                    return false;

                if(a.Equals(b))
                    return true;

                if(isInt(a) && isInt(b)) {
                    long _a = Convert.ToInt64(a);
                    long _b = Convert.ToInt64(b);

                    return _a == _b;
                }

                if(isFloat(a) && isFloat(b)) {
                    double _a = Convert.ToDouble(a);
                    double _b = Convert.ToDouble(b);

                    return _a == _b;
                }

                return false;
            }

        }

        /// <summary>
        /// Creates the <see cref="IDatabaseInfo"/> for the current settings of the control file #
        /// (<see cref="AppControl.DbPath"/>, <see cref="AppControl.AlternateDbPaths"/>),
        /// if accessible on the current computer.
        /// </summary>
        static public IDatabaseInfo GetDatabase(this AppControl Control) {
            List<ValueTuple<string, string>> allPaths = new List<(string, string)>();
            if(!Control.DbPath.IsNullOrEmpty())
                allPaths.Add((Control.DbPath, null));
            if(Control.AlternateDbPaths != null)
                allPaths.AddRange(Control.AlternateDbPaths);

            string mName = System.Environment.MachineName.ToLowerInvariant();

            string dbPath = null;
            foreach(var t in allPaths) {
                string path = t.Item1;
                string filter = t.Item2;

                if(path.IsNullOrEmpty() || path.IsEmptyOrWhite())
                    continue;

                if(!filter.IsNullOrEmpty() && !filter.IsEmptyOrWhite()) {
                    if(!mName.Contains(filter)) {
                        continue;
                    }
                }

                if(System.IO.Directory.Exists(path)) {
                    dbPath = path;
                    break;
                }
            }

            if(dbPath == null) {
                return null;
            }

            // try to match it with on of the already known databases
            foreach(var db in BoSSSshell.databases) {
                if(db.PathMatch(dbPath))
                    return db;
            }

            // otherwise, create new db
            return DatabaseInfo.Open(dbPath);
        }
    }
}

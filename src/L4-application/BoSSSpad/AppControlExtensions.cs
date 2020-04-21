using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.Control;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using ilPSP;

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
        /// Runs the solver described by the control object <paramref name="ctrl"/> on a batch system.
        /// The method returns immediately.
        /// </summary>
        /// <param name="ctrl"></param>
        /// <param name="BatchSys"></param>
        /// <returns></returns>
        public static Job RunBatch(this AppControl ctrl, BatchProcessorClient BatchSys) {
            ctrl.ProjectName = InteractiveShell.WorkflowMgm.CurrentProject;


            string JobName = ctrl.SessionName;
            int ctrl_idx = InteractiveShell.WorkflowMgm.RegisterControl(ctrl);
            if(JobName.IsEmptyOrWhite()) {
                JobName = "UnnamedJob_" + ctrl_idx;
                ctrl.SessionName = JobName;
            }

            Type solverClass = ctrl.GetSolverType();
            Job job = new Job(JobName, solverClass);

            //job.ExecutionTime = executionTime;
            //job.NumberOfMPIProcs = NumberOfMPIProcs;
            //job.UseComputeNodesExclusive = UseComputeNodesExclusive;
            job.SetControlObject(ctrl);
            job.Activate(BatchSys);
            
            return job;
        }

        /// <summary>
        /// Runs the solver described by the control object <paramref name="ctrl"/> on a batch system from the currently defined queues (<see cref="InteractiveShell.ExecutionQueues"/>).
        /// The method returns immediately.
        /// </summary>
        /// <param name="ctrl"></param>
        /// <param name="queueIdx">
        /// Index int <see cref="InteractiveShell.ExecutionQueues"/>
        /// </param>
        public static Job RunBatch(this AppControl ctrl, int queueIdx = 0) {
            var b = InteractiveShell.ExecutionQueues[queueIdx];
            return RunBatch(ctrl, queueIdx);
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
            ctrl.ProjectName = InteractiveShell.WorkflowMgm.CurrentProject;

            Type solverClass = ctrl.GetSolverType();
            Job job = new Job(ctrl.SessionName , solverClass);
            job.SetControlObject(ctrl);

            
            return job;
        }

        /// <summary>
        /// Returns the job correlated to a control object
        /// </summary>
        public static Job GetJob(this AppControl ctrl) {
            foreach(var j in InteractiveShell.WorkflowMgm.AllJobs.Values) {
                var cj = j.GetControl();
                if(cj == null)
                    continue;

                if(cj.Equals(ctrl))
                    return j;
            }
            Console.WriteLine("No Job assigned for given control object yet.");
            return null;
        }

        /// <summary>
        /// Returns all sessions which can be correlated to a specific control object
        /// </summary>
        public static ISessionInfo[] GetAllSessions(this AppControl ctrl) {
            var AllCandidates = InteractiveShell.WorkflowMgm.Sessions.Where(
                    sinf => InteractiveShell.WorkflowMgm.SessionInfoAppControlCorrelation(sinf, ctrl));

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
            job.MemPerCPU = job2rest.MemPerCPU;
            job.ExecutionTime = job2rest.ExecutionTime;
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
            if(ctrl.DynamicLoadBalancing_CellCostEstimatorFactories.Count != 0)
                throw new ArgumentException("'DynamicLoadBalancing_CellCostEstimatorFactories' is not supported - cannot be serialized.");

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
                    throw new ArgumentException($"Unable to serialize/deserialize initial values for '{ivk}' correctly. Note that you cannot use delegates with the workflow management, use either formulas as text or 'GetFormulaObject(...)' in the BoSSSpad.");
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

                if((!D2.ContainsKey(name)) || (!o1.Equals(D2[name]))) {
                    throw new ArgumentException("Unable to serialize/deserialize '" + name + "' correctly. (Missing a DataMemberAtribute in control class?)");
                }

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
            foreach(var db in InteractiveShell.databases) {
                if(db.PathMatch(dbPath))
                    return db;
            }

            // otherwise, create new db
            return new DatabaseInfo(dbPath);
        }
    }
}

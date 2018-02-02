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
            ((IDisposable)solver).Dispose();

            return S;
        }

        /// <summary>
        /// Runs the solver described by the control object <paramref name="ctrl"/> on a batch system.
        /// The method returns immediately.
        /// </summary>
        /// <param name="ctrl"></param>
        /// <param name="BatchSys"></param>
        /// <param name="NumberOfMPIProcs"></param>
        /// <param name="UseComputeNodesExclusive"></param>
        /// <returns></returns>
        public static Job RunBatch(this AppControl ctrl, BatchProcessorClient BatchSys, int NumberOfMPIProcs = 1, bool UseComputeNodesExclusive = false) {
            ctrl.ProjectName = InteractiveShell.WorkflowMgm.CurrentProject;

            Type solverClass = ctrl.GetSolverType();
            Job job = new Job(ctrl.SessionName , solverClass);
            job.NumberOfMPIProcs = NumberOfMPIProcs;
            job.UseComputeNodesExclusive = UseComputeNodesExclusive;
            job.SetControlObject(ctrl);
            job.Activate(BatchSys);
            
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
            string JSON = ctrl.Serialize();
            AppControl ctrlBack = AppControl.Deserialize(JSON);//, ctrl.GetType());
            ctrlBack.Verify();

            // compare original and de-serialized object
            if (!ctrl.Tags.SetEquals(ctrlBack.Tags))
                throw new ArgumentException("Unable to serialize/deserialize tags correctly.");

            if(!ctrl.InitialValues.Keys.SetEquals(ctrlBack.InitialValues.Keys))
                throw new ArgumentException("Unable to serialize/deserialize initial values correctly.");

            foreach(var ivk in ctrl.InitialValues.Keys) {
                var f1 = ctrl.InitialValues[ivk];
                var f2 = ctrlBack.InitialValues[ivk];

                if(!f1.Equals(f2))
                    throw new ArgumentException("Unable to serialize/deserialize initial values correctly.");
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

            if(!ctrl.BoundaryValues.Keys.SetEquals(ctrlBack.BoundaryValues.Keys))
                throw new ArgumentException("Unable to serialize/deserialize boundary values correctly.");

            foreach(var bvk in ctrl.BoundaryValues.Keys) {
                var bvc = ctrl.BoundaryValues[bvk];
                var bvd = ctrlBack.BoundaryValues[bvk];

                if(!bvc.Evaluators.Keys.SetEquals(bvd.Evaluators.Keys))
                    throw new ArgumentException("Unable to serialize/deserialize boundary values correctly.");

                if(!bvc.Values.Keys.SetEquals(bvd.Values.Keys))
                    throw new ArgumentException("Unable to serialize/deserialize boundary values correctly.");

                foreach(string s in bvc.Values.Keys) {
                    var f1 = bvc.Values[s];
                    var f2 = bvd.Values[s];

                    if(!f1.Equals(f2))
                        throw new ArgumentException("Unable to serialize/deserialize boundary values correctly.");
                }
            }

            // reflection comparison
            var D1 = new Dictionary<string, object>();
            var D2 = new Dictionary<string, object>();
            BoSSS.Solution.Application.FindKeys(D1, ctrl);
            BoSSS.Solution.Application.FindKeys(D2, ctrl);
            foreach(var kv1 in D1) {
                string name = kv1.Key;
                var o1 = kv1.Value;

                if((!D2.ContainsKey(name)) || (!o1.Equals(D2[name]))) {
                    throw new ArgumentException("Unable to serialize/deserialize '" + name + "' correctly. (Missing a DataMemberAtribute in control class?)");
                }

            }

        }
    }
}

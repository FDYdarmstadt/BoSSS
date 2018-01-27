using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.Control;
using BoSSS.Foundation.IO;
using BoSSS.Solution;

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
        /// <returns></returns>
        public static Job RunBatch(this AppControl ctrl, BatchProcessorClient BatchSys, int NumberOfMPIProcs = 1) {
            ctrl.ProjectName = InteractiveShell.WorkflowMgm.CurrentProject;

            Type solverClass = ctrl.GetSolverType();
            Job job = new Job(ctrl.SessionName , solverClass);
            job.NumberOfMPIProcs = NumberOfMPIProcs;
            job.SetControlObject(ctrl);
            job.Activate(BatchSys);
            
            return job;
        }

    }
}

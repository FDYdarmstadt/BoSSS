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

        //public static void Verify_(this AppControl ctrl, bool TestForSerilaization = true) {
        //    var AllProblems = new List<string>();

        //    if(ctrl.GridFunc != null) 

        //}


        public static SessionInfo Run(this AppControl ctrl) {

            var solverClass = ctrl.GetSolverType();
            object solver = Activator.CreateInstance(solverClass);

            //var app = (BoSSS.Solution.Application)solver;

            //IApplication<AppControl> app2 = (IApplication<AppControl>)app;
            //
            //app2.Init(ctrl, null);

            IApplication app = (IApplication)solver;

            var S = app.CurrentSessionInfo;
            ((IDisposable)solver).Dispose();

            return S;
        }


    }
}

using System;
using System.Runtime.Serialization;
using BoSSS.Solution;
using BoSSS.Solution.Control;


namespace BoSSS.Application.TraceDGtest {
    public class TraceDGtestControl : AppControlSolver {
        public TraceDGtestControl() {
            SetDGdegree(1);
        }

        public override void SetDGdegree(int p) {
            base.FieldOptions.Clear();
            AddFieldOption("Velocity*", p);
            AddFieldOption("SurfaceConcentration", p);
            AddFieldOption("Phi", Math.Min(2, p));
        }

        override public Type GetSolverType() {
            return typeof(TraceDGtestMain);
        }


       

    }

}

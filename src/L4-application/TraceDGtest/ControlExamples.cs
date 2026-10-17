using System;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;

namespace BoSSS.Application.TraceDGtest {
    public class ControlExamples {

        /// <summary>
        /// Example: 2D Cartesian grid with AMR at the Level Set.
        /// </summary>
        public static TraceDGtestControl SteadystateWithLevelSetAMR(int degree = 2, bool useAMR = false) {
            var c = new TraceDGtestControl();

            // Project name and basic field options
            c.ProjectName = "TraceDG_Cartesian_LevelSetAMR";
            c.SetDGdegree(degree);

            // Grid generation: 2D Cartesian grid
            c.GridFunc = delegate {
                var xNodes = GenericBlas.Linspace(-1, 1, 11);
                var yNodes = GenericBlas.Linspace(-1, 1, 11);
                return Grid2D.Cartesian2DGrid(xNodes, yNodes);
            };

            // Initial values
            c.InitialValues_Evaluators.Add("VelocityX", X => 1.0);
            c.InitialValues_Evaluators.Add("SurfaceConcentration", (double[] X) => X[0].Pow2() + X[1].Pow2());
            c.InitialValues_Evaluators.Add("RHS", (double[] X) => X[0].Pow2() + X[1].Pow2());

            // Example: Level Set function (circle)
            c.InitialValues_Evaluators.Add("Phi", X => 0.33 - (X[0].Pow2() + X[1].Pow2()));

           
            // turn on AMR
            c.AdaptiveMeshRefinement = useAMR;
            c.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 2 });
            c.AMR_startUpSweeps = 2;

            // steady state
            c.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            // return
            c.SkipSolveAndEvaluateResidual = true;
            return c;
        }

    }
}
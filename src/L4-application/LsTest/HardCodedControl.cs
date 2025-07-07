using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.LsTest {
    public class HardCodedControl {

        /// <summary>
        /// Initialize quadratic Levelset and test if an artificial interface is created at the inflow.
        /// If the algorithm works correct, the single planar initial interface
        /// is transported along the channel with a steady velocity.
        /// If it is not working correct, due to the non monotonic initial level set, 
        /// an artificial second interface is created at the inflow. 
        /// And/or the interface moves with a speed unequal to the advection velocity.
        /// </summary>
        /// <returns></returns>
        public static SolverWithLevelSetUpdaterTestControl BoundaryConditionTest() {
            var C = new SolverWithLevelSetUpdaterTestControl();

            C.DegreeOfLevelSets = 2;
            C.SetDGdegree(2);

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, 5, 5 + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 1, 2);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                
                grd.EdgeTagNames.Add(1, "wall");                

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 1;
                    return et;
                });

                return grd;
            };

            C.AddInitialValue("Phi", $"(X, t) => -X[0] * X[0] + 1.0", true);
            //C.AddInitialValue("Phi", $"(X, t) => -X[0] + 1.0", true);

            // initial values and exact solution
            // =================================
            C.SetAdvectionVelocity(0, new Func<double[], double, double>[] { (X, t) => 1.0, (X, t) => 1.0 });


            // timestepping and solver
            // =======================
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.LieSplitting;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler; // should not really matter, we are only projecting the exact underlying advection velocity in each timestep

            C.dtFixed = 0.01;
            C.NoOfTimesteps = 500;
            C.Endtime = 5.0;

            return C;
        }
    }
}

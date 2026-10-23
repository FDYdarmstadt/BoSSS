using System;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using ilPSP.Utils;

namespace BoSSS.Application.XNSE_Solver {

    /// <summary>
    /// Small control objects for the tests of the workflow management (BoSSSpad);
    /// they are located in the solver assembly, so that they can be used as control statements (`cs:...`)
    /// on a batch system: a deployed job only contains the assemblies which the solver depends on.
    /// </summary>
    public static class WorkflowTestControls {

        /// <summary>
        /// Session name set by <see cref="ChannelControl"/>; for name-based correlation of sessions and jobs, 
        /// the job must have the same name.
        /// </summary>
        public const string ChannelSessionName = "Channel_40x8k1";

        /// <summary>
        /// Grid for <see cref="ChannelControl"/>
        /// </summary>
        public static GridCommons ChannelGrid() {
            var xNodes = GenericBlas.Linspace(0, 10, 41);
            var yNodes = GenericBlas.Linspace(-1, 1, 9);
            GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes);
            IGrid_Extensions.DefineEdgeTags(grid, delegate (double[] X) {
                double x = X[0];
                double y = X[1];
                if(Math.Abs(y - (-1)) <= 1.0e-8)
                    return "wall"; // lower wall
                if(Math.Abs(y - (+1)) <= 1.0e-8)
                    return "wall"; // upper wall
                if(Math.Abs(x - (0.0)) <= 1.0e-8)
                    return "Velocity_Inlet"; // inlet
                if(Math.Abs(x - (+10.0)) <= 1.0e-8)
                    return "Pressure_Outlet"; // outlet
                throw new ArgumentOutOfRangeException("unknown domain");
            });
            return grid;
        }

        /// <summary>
        /// Steady-state channel flow, which runs within a few seconds.
        /// </summary>
        /// <param name="solverCreatesGrid">
        /// If true, the grid is created (and saved to the database) by the solver; this is required if the control is used as a control statement 
        /// (`cs:BoSSS.Application.XNSE_Solver.WorkflowTestControls.ChannelControl(true)`, see `Job.SetControlStatement`),
        /// since then, the job management cannot save a grid before the job starts.
        /// The database is not set here, but passed to the solver by the job (`--dbguid`).
        /// </param>
        public static XNSE_Control ChannelControl(bool solverCreatesGrid = false) {
            var c = new XNSE_Control();

            // general description:
            int k = 1;
            string desc = "Steady state, channel, k" + k;
            c.ProjectDescription = desc;
            c.savetodb = true;
            c.Tags.Add("k" + k);
            // DG polynomial degree
            c.SetDGdegree(k);
            // Physical parameters:
            double reynolds = 20;
            c.PhysicalParameters.rho_A = 1;
            c.PhysicalParameters.mu_A = 1.0 / reynolds;
            // Timestepping properties:
            c.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            c.SessionName = ChannelSessionName;

            c.AddBoundaryValue("Velocity_Inlet", "VelocityX", new Formula("X => 1 - X[0]*X[0]", false));

            // setting the grid:
            if(solverCreatesGrid) {
                c.GridFunc = ChannelGrid;
            } else {
                c.SetGrid(ChannelGrid());
            }

            return c;
        }
    }
}

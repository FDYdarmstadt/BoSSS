using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.LevelSetTools;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution;

namespace XNSE_ParallelTets {

    /// <summary>
    /// Tests 
    /// </summary>
    public class Controls {

        static public XNSE_Control Test_ChannelFlow2D(bool nonlinear = true, bool useAMR = true) {

            var C = new XNSE_Control();

            C.SetDGdegree(2);

            C.PhysicalParameters = new PhysicalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,

                mu_A = 1.0,
                mu_B = 1.0,

                Sigma = 0.1,

                IncludeConvection = nonlinear,
                Material = true
            };


            // grid generation
            // ===============
            #region grid

            double H = 2;
            int Lscale = 1;
            double L = Lscale * H;
            int kelem = 8;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, (Lscale * kelem) + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double U = 0.125;

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");

            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => ((-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0*Math.PI*(t/T)));
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", (X, t) => ((-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0 * Math.PI * (t / T)));

            C.AddBoundaryValue("pressure_outlet_right");

            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("VelocityX#A", X => (-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => (-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U);

            //C.InitialValues_Evaluators.Add("GravityX#A", X => 5.0);
            //C.InitialValues_Evaluators.Add("GravityX#B", X => 5.0);


            double[] center = new double[] { (H / 2.0) + 0.0, H / 2.0 };
            double radius = 0.4;

            C.InitialValues_Evaluators.Add("Phi",
                (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius)  // signed-distance form
                );

            #endregion


            // solver options
            // ==============
            #region solver

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            C.ReInitPeriod = 4;

            // C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            // C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.MaxSolverIterations = 50;

            // C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.5;
            C.NoOfTimesteps = 4;

            {
                C.AdaptiveMeshRefinement = useAMR;
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
                //C.activeAMRlevelIndicators.Add(new AMRonNarrowbandAtBoundary(new byte[] { 1 }) { maxRefinementLevel = AMRlevel_dropBL });
                C.activeAMRlevelIndicators.Add(new AMRLevelIndicatorLibrary.AMRonBoundary(new byte[] { 1 }) { maxRefinementLevel = 1 });
                C.AMR_startUpSweeps = 1;
            }

            #endregion


            return C;
        
        }

        static public XNSE_Control Test_ChannelFlow3D(bool nonlinear = true, bool useAMR = true) {

            var C = new XNSE_Control();

            C.SetDGdegree(2);

            C.PhysicalParameters = new PhysicalParameters() {
                rho_A = 1.0,
                rho_B = 1.0,

                mu_A = 1.0,
                mu_B = 1.0,

                Sigma = 0.1,

                IncludeConvection = nonlinear,
                Material = true
            };


            // grid generation
            // ===============
            #region grid

            double H = 2;
            double W = H;
            int Lscale = 1;
            double L = Lscale * H;
            int kelem = 4;

            C.GridFunc = delegate () {
                double[] Ynodes = GenericBlas.Linspace(0, L, (Lscale * kelem) + 1);
                double[] Xnodes = GenericBlas.Linspace(0, W, kelem + 1);
                double[] Znodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes, periodicY: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                grd.EdgeTagNames.Add(3, "velocity_inlet_front");
                grd.EdgeTagNames.Add(4, "velocity_inlet_back");

                grd.EdgeTagNames.Add(5, "velocity_inlet_left");
                grd.EdgeTagNames.Add(6, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[2]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[2] - H) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - W) <= 1.0e-8)
                        et = 4;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 5;
                    if (Math.Abs(X[1] - L) <= 1.0e-8)
                        et = 6;

                    return et;
                });

                return grd;
            };

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double U = 0.125;

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");

            C.AddBoundaryValue("velocity_inlet_front", "VelocityX#A", (X, t) => ((-4.0 * U / H.Pow2()) * (X[2] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0*Math.PI*(t/T)));
            C.AddBoundaryValue("velocity_inlet_front", "VelocityX#B", (X, t) => ((-4.0 * U / H.Pow2()) * (X[2] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0 * Math.PI * (t / T)));

            C.AddBoundaryValue("velocity_inlet_back", "VelocityX#A", (X, t) => ((-4.0 * U / H.Pow2()) * (X[2] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0*Math.PI*(t/T)));
            C.AddBoundaryValue("velocity_inlet_back", "VelocityX#B", (X, t) => ((-4.0 * U / H.Pow2()) * (X[2] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0 * Math.PI * (t / T)));

            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => ((-4.0 * U / H.Pow2()) * (X[2] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0*Math.PI*(t/T)));
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", (X, t) => ((-4.0 * U / H.Pow2()) * (X[2] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0 * Math.PI * (t / T)));

            C.AddBoundaryValue("pressure_outlet_right");

            #endregion


            // Initial Values
            // ==============
            #region init

            C.InitialValues_Evaluators.Add("VelocityX#A", X => (-4.0 * U / H.Pow2()) * (X[2] - H / 2.0).Pow2() + U);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => (-4.0 * U / H.Pow2()) * (X[2] - H / 2.0).Pow2() + U);

            //C.InitialValues_Evaluators.Add("GravityX#A", X => 5.0);
            //C.InitialValues_Evaluators.Add("GravityX#B", X => 5.0);


            double[] center = new double[] { W / 2.0, (H / 2.0) + 0.0, L / 2.0 };
            double radius = 0.4;

            C.InitialValues_Evaluators.Add("Phi",
                (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() + (X[2] - center[2]).Pow2()).Sqrt() - radius)  // signed-distance form
                );

            #endregion


            // solver options
            // ==============
            #region solver

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            C.ReInitPeriod = 4;

            // C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            // C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.MaxSolverIterations = 50;

            // C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 0.5;
            C.NoOfTimesteps = 4;

            {
                C.AdaptiveMeshRefinement = useAMR;
                C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
                //C.activeAMRlevelIndicators.Add(new AMRonNarrowbandAtBoundary(new byte[] { 1 }) { maxRefinementLevel = AMRlevel_dropBL });
                C.activeAMRlevelIndicators.Add(new AMRLevelIndicatorLibrary.AMRonBoundary(new byte[] { 1 }) { maxRefinementLevel = 1 });
                C.AMR_startUpSweeps = 1;
            }

            #endregion


            return C;

        }
    }
}

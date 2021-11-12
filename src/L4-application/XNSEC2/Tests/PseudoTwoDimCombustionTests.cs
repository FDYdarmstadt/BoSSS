using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using ilPSP;
using ilPSP.Utils;
using System;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    /// </summary>
    static public partial class FullNSEControlExamples {

        public static XNSEC_Control PseudoTwoDimensionalTwoPhaseFlow(int p = 1, int kelem = 5) {
            XNSEC_Control C = new XNSEC_Control();

            string _DbPath = null;

            int D = 2;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.ImmediatePlotPeriod = 1;
            // basic database options
            // ======================

            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSEC/Channel";
            C.ProjectDescription = "XNSEC testing";
            #endregion db

            // DG degrees
            // ==========
            C.SetDGdegree(p);

            // Physical Parameters
            // ===================

            #region physics

            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            double sigma = 0.0; // No surface tension
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;
            C.rhoOne = true;
            C.ChemicalReactionActive = false;
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };
            C.EnableTemperature = false;
            C.EnableMassFractions = false;
            C.UseSelfMadeTemporalOperator = false;
            #endregion physics

            // grid generation
            // ===============

            #region grid

            double L = 5;
            double H = 2;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 4);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicY: true);
                grd.EdgeTagNames.Add(1, "velocity_inlet_left");
                grd.EdgeTagNames.Add(2, "velocity_inlet_right");
                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;         
                    if(Math.Abs(X[0]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                return grd;
            };

            #endregion grid

            // Initial Values
            // ==============

            #region init
            Func<double[], double> PhiFunc = (X => X[0] - 2.1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);
            double U = 0.0;
            C.InitialValues_Evaluators.Add("VelocityX#A", X => 1.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 1.0);
            #endregion init

            // exact solution
            // ==============

            #region exact
            C.Phi = ((X, t) => PhiFunc(X));
            #endregion exact

            // boundary conditions
            // ===================

            #region BC
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => 1.0);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", (X, t) => 1.0);
            C.AddBoundaryValue("velocity_inlet_left", "Temperature#A", (X, t) => 1.0);
            C.AddBoundaryValue("velocity_inlet_left", "Temperature#B", (X, t) => 1.0);


            C.AddBoundaryValue("velocity_inlet_right", "VelocityX#A", (X, t) => 1.0);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityX#B", (X, t) => 1.0);
            C.AddBoundaryValue("velocity_inlet_right", "Temperature#A", (X, t) => 1.0);
            C.AddBoundaryValue("velocity_inlet_right", "Temperature#B", (X, t) => 1.0);

            #endregion BC


            #region solver

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;


            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.None;



   

            #endregion solver

            // Timestepping
            // ============

            #region time

            C.SkipSolveAndEvaluateResidual = false;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.MaxSolverIterations = 5;
            C.saveperiod = 10;

            #endregion time

            return C;
        }



        public static XNSEC_Control PseudoTwoDimensionalTwoPhaseFlow_DifferentDensities(int p = 1, int kelem = 15, int wallBC = 0) {
            XNSEC_Control C = new XNSEC_Control();

            string _DbPath = null;

            int D = 2;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.ImmediatePlotPeriod = 1;
            // basic database options
            // ======================

            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSEC/Channel";
            C.ProjectDescription = "Channel flow testing";

            //C.ContinueOnIoError = false;
            //C.LogValues = XNSE_Control.LoggingValues.ChannelFlow;

            #endregion db

            // DG degrees
            // ==========
            C.SetDGdegree(p);

            // Physical Parameters
            // ===================

            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 2;
            C.PhysicalParameters.mu_A = 1*0;
            C.PhysicalParameters.mu_B = 1*0;
            double sigma = 0.0; // No surface tension
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;
            C.rhoOne = true;
            C.ChemicalReactionActive = false;
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };
            C.EnableTemperature = false;
            C.EnableMassFractions = false;
            C.UseSelfMadeTemporalOperator = false;
            #endregion physics

            // grid generation
            // ===============

            #region grid

            double L = 5;
            double H = 2;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, 4);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicY: true);
                //var grd = Grid2D.UnstructuredTriangleGrid(Xnodes, Ynodes);


                grd.EdgeTagNames.Add(1, "velocity_inlet_left");
                grd.EdgeTagNames.Add(2, "velocity_inlet_right");


                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                return grd;
            };

            #endregion grid

            // Initial Values
            // ==============

            #region init

            Func<double[], double> PhiFunc = (X => X[0] - 2.0);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double U = +1.0;
            C.InitialValues_Evaluators.Add("VelocityX#A", X => U);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => U*2);

            C.InitialValues_Evaluators.Add("Temperature#A", X => 1.0);
            C.InitialValues_Evaluators.Add("Temperature#B", X => 1.0);

            C.InitialValues_Evaluators.Add("MassFraction0#A", X => 1.0);
            C.InitialValues_Evaluators.Add("MassFraction0#B", X => 1.0);

            #endregion init

            // exact solution
            // ==============

            #region exact

            C.Phi = ((X, t) => PhiFunc(X));



            #endregion exact

            // boundary conditions
            // ===================

            #region BC


            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => U);
            //C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", (X, t) => U);
            C.AddBoundaryValue("velocity_inlet_left", "Temperature#A", (X, t) => 1.0);
            C.AddBoundaryValue("velocity_inlet_left", "Temperature#B", (X, t) => 1.0);


            //C.AddBoundaryValue("velocity_inlet_right", "VelocityX#A", (X, t) => U);
            C.AddBoundaryValue("velocity_inlet_right", "VelocityX#B", (X, t) => U*2);
            C.AddBoundaryValue("velocity_inlet_right", "Temperature#A", (X, t) => 1.0);
            C.AddBoundaryValue("velocity_inlet_right", "Temperature#B", (X, t) => 1.0);

            #endregion BC

            // misc. solver options
            // ====================

            #region solver

            //C.ComputeEnergyProperties = true;
            //C.solveKineticEnergyEquation = true;
            ////C.CheckJumpConditions = true;
            //C.kinEViscousDiscretization = Solution.EnergyCommon.KineticEnergyViscousSourceTerms.laplaceKinE;
            //C.kinEPressureDiscretization = Solution.EnergyCommon.KineticEnergyPressureSourceTerms.divergence;
            //C.withDissipativePressure = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.0;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.Phi = (X,t) => ((X[0] - (center[0]+U*t)).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius;

            //C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            //C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;
            //C.useFiltLevSetGradientForEvolution = true;
            //C.ReInitPeriod = 1;
            //C.ReInitOnRestart = true;

            C.Option_LevelSetEvolution = LevelSetEvolution.None;
            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.None;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.SkipSolveAndEvaluateResidual = true;

            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
            C.AMR_startUpSweeps = 1;
            //C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            //C.BaseRefinementLevel = 2;
            //C.RefinementLevel = 2;

            C.InitSignedDistance = false;
            C.adaptiveReInit = false;

            #endregion solver

            // Timestepping
            // ============

            #region time

            C.SkipSolveAndEvaluateResidual = true;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            //C.dtFixed = 0.1;

            //C.NonLinearSolver.MinSolverIterations = 1;
            //C.NonLinearSolver.MaxSolverIterations = 1;
            //C.NoOfTimesteps = 100; // 500;
            C.saveperiod = 10;

            #endregion time

            return C;
        }

    }
}
using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using ZwoLevelSetSolver.SolidPhase;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Utils;

namespace ZwoLevelSetSolver.ControlFiles {

    public static class HardCodedControl {

        public static ZLS_Control QuasiStationaryDroplet(int p = 1, double agglomerationTreshold = 0.0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.AgglomerationThreshold = agglomerationTreshold;
            C.NoOfMultigridLevels = 1;
            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            // basic database options
            // ======================
            #region db
            C.savetodb = false;
            C.ProjectName = "XNSE/Droplet";

            C.ContinueOnIoError = false;
            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Reusken");
            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;
            double sigma = 0.5;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 3.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            int kelem = 15;
            double xSize = 0.25;
            double yTop = 0.25;
            int overlap = 1;
            double yBottom = yTop / (1 - (kelem / 2 )/ (overlap + 0.5));

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem / 2 + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "freeslip_lower");
                grd.EdgeTagNames.Add(2, "freeslip_upper");
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.1;
            double Theta_e = Math.PI / 2.0;

            double[] center = new double[] { 0, 0 };

            Func<double[], double> PhiFunc = (X => ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double pJump = sigma / R;
            C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //Func<double[], double> Phi1Func = (X => -(X[1] - 0.02 + 0.4 * X[0] * X[0]));
            Func<double[], double> Phi1Func = (X => -(X[1] -0.000));
            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            double xVelocity = 0.2;
            C.InitialValues_Evaluators.Add("VelocityX#A", X => xVelocity);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => xVelocity);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("freeslip_lower");
            C.AddBoundaryValue("freeslip_upper");
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => xVelocity);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => xVelocity);
            C.AddBoundaryValue("pressure_outlet_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            //C.PhysicalParameters.sliplength = 0.001;

            #endregion


            // exact solution
            // ==============

            // misc. solver options
            // ====================
            #region solver
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;



            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = compMode;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.NoOfTimesteps = 300;
            C.Endtime = 0.3;
            C.saveperiod = 1;

            #endregion


            return C;
        }

        public static ZLS_Control RecedingDroplet_OnPlate (int p = 3, int kelem = 24, int AMRlvl = 0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            
            C.SuperSampling = 3;

            C.AgglomerationThreshold = 0.1;
            C.NoOfMultigridLevels = 1;

            int D = 2;
            
            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb =false;
            C.ProjectName = "XNSE/Droplet";
            //C.ProjectDescription = "Static droplet on plate";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Reusken");
            C.PhysicalParameters.rho_A = 0.1;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.1;
            C.PhysicalParameters.mu_B = 0.1;
            double sigma = 0.01;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 0.25;
            double yTop = 0.2;
            int overlap = 5;
            double yBottom = yTop / (1 - (kelem/2 ) / (overlap + 0.5)) + 0.001;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem/2 + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] -yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.1;
            double Theta_e = Math.PI / 1.6;
            double s = 2 * R * Math.Sin(Theta_e);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0, -h};

            Func<double[], double> PhiFunc = Phi;

            double Phi(double[] X) {
                return ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2();
                //return ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R;
            }


            double pJump = sigma / R;

            C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);


            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            //Func<double[], double> Phi1Func = (X => -(X[1] -0.02 + 0.4 * X[0] * X[0]));
            Func<double[], double> Phi1Func = (X => -(X[1] -0.000));
            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("wall_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            //C.PhysicalParameters.sliplength = 0.001;

            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 80;
            C.NonLinearSolver.MinSolverIterations = 60;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-10;
            C.LinearSolver.ConvergenceCriterion = 1e-10;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;

            

            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            //C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = 1 });
            //C.AMR_startUpSweeps = 0;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Iterative;

            C.TimesteppingMode = compMode;
            double dt = 5e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static ZLS_Control Droplet_OnPlate(int p = 2, int kelem = 6, int AMRlvl = 0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;

            C.SuperSampling = 3;

            C.AgglomerationThreshold = 0.1;
            C.NoOfMultigridLevels = 1;

            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            C.ProjectName = "XNSE/Droplet";
            //C.ProjectDescription = "Static droplet on plate";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.Tags.Add("Reusken");
            C.PhysicalParameters.rho_A = 0.1;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.1;
            C.PhysicalParameters.mu_B = 0.1;
            double sigma = 0.01;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 0.5;
            double yTop = 0.5 - 0.001;
            double yBottom = -0.5 - 0.001;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 0.06;
            double Theta_e = Math.PI / 2;
            double s = 2 * R * Math.Sin(Theta_e);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0, -h };

            Func<double[], double> PhiFunc = Phi;

            double Phi(double[] X) {
                //return 3 *  (((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2());
                return ((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - R;
            }


            double pJump = sigma / R;

            C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);


            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            //Func<double[], double> Phi1Func = (X => -(X[1] -0.02 + 0.4 * X[0] * X[0]));
            Func<double[], double> Phi1Func = (X => -( X[1] - 0.000));
            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("wall_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            //C.PhysicalParameters.sliplength = 0.001;

            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 80;
            C.NonLinearSolver.MinSolverIterations = 2;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-10;
            C.LinearSolver.ConvergenceCriterion = 1e-10;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;



            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = BoSSS.Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = 3 });
            C.AMR_startUpSweeps = 2;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;

            C.TimesteppingMode = compMode;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static XNSE_Control StaticDroplet_Free(int p = 3, int kelem = 20, int AMRlvl = 0) {

            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.0;
            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;
            bool steadyInterface = true;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_test_db";
            string _DbPath = null; // @"\\HPCCLUSTER\hpccluster-scratch\smuda\XNSE_studyDB";
            //string _DbPath = @"\\terminal03\Users\smuda\local\terminal03_XNSE_studyDB";

            // basic database options
            // ======================
            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "StaticDroplet";
            //C.ProjectDescription = "Static droplet";
            //C.SessionName = "SD_meshStudy_Hysing_mesh" + kelem; // "_AMR"+AMRlvl;

            C.ContinueOnIoError = false;
            //C.LogValues = XNSE_Control.LoggingValues.Dropletlike;

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            //C.Tags.Add("Hysing");
            C.Tags.Add("La = 5000");
            C.PhysicalParameters.rho_A = 1e4;
            C.PhysicalParameters.rho_B = 1e4;
            C.PhysicalParameters.mu_A = 1.0;
            C.PhysicalParameters.mu_B = 1.0;
            double sigma = 1.0;
            C.PhysicalParameters.Sigma = sigma;

            //C.Tags.Add("La = 0.005");
            //C.PhysicalParameters.rho_A = 1;
            //C.PhysicalParameters.rho_B = 1;
            //C.PhysicalParameters.mu_A = 10;
            //C.PhysicalParameters.mu_B = 10;
            //double sigma = 1;
            //C.PhysicalParameters.Sigma = sigma;

            //Air - Water(lenght scale == centimeters, 3D space)
            //C.PhysicalParameters.rho_A = 1e3;      // kg / cm^3
            //C.PhysicalParameters.rho_B = 1.2;    // kg / cm^3
            //C.PhysicalParameters.mu_A = 1e-3;       // kg / cm * sec
            //C.PhysicalParameters.mu_B = 17.1e-6;    // kg / cm * sec
            //double sigma = 72.75e-3;                // kg / sec^2 
            //C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double Lscl = 1.0;
            double xSize = Lscl * 1.0;
            double ySize = Lscl * 1.0;
            double zSize = Lscl * 1.0;

            if(D == 2) {
                C.GridFunc = delegate () {
                    double[] Xnodes = GenericBlas.Linspace(-xSize / 2.0, xSize / 2.0, kelem + 0);
                    double[] Ynodes = GenericBlas.Linspace(-ySize / 2.0, ySize / 2.0, kelem + 0);
                    Ynodes[0] = -ySize / 2.0 - Lscl* 0.1;
                    var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
                    //var grd = Grid2D.UnstructuredTriangleGrid(Xnodes, Ynodes);

                    grd.EdgeTagNames.Add(1, "wall_lower");
                    grd.EdgeTagNames.Add(2, "wall_upper");
                    grd.EdgeTagNames.Add(3, "wall_left");
                    grd.EdgeTagNames.Add(4, "wall_right");

                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if(Math.Abs(X[1] + ySize / 2.0 + Lscl * 0.1) <= 1.0e-8)
                            et = 1;
                        if(Math.Abs(X[1] - ySize / 2.0) <= 1.0e-8)
                            et = 2;
                        if(Math.Abs(X[0] + xSize / 2.0) <= 1.0e-8)
                            et = 3;
                        if(Math.Abs(X[0] - xSize / 2.0) <= 1.0e-8)
                            et = 4;

                        return et;
                    });

                    return grd;
                };
            }

            #endregion


            // Initial Values
            // ==============
            #region init

            double r = Lscl * 0.25;

            Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()).Sqrt() - r);         // signed distance
            //Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()) - r.Pow2());         // quadratic
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = (X => -(X[1] + r + 0.1 ));
            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0.0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0.0);

            double Pjump = sigma / r;
            if(D == 3)
                Pjump *= 2.0;

            C.InitialValues_Evaluators.Add("Pressure#A", X => Pjump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            //C.InitialValues_Evaluators.Add("GravityY#A", X => 0.0);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => 0.0);


            // restart
            //var database = new DatabaseInfo(_DbPath);
            //Guid restartID = new Guid("c95b2833-288e-4014-ba08-affa65a2398e");
            //C.RestartInfo = new Tuple<Guid, Foundation.IO.TimestepNumber>(restartID, null);

            #endregion


            // exact solution
            // ==============
            #region exact

            C.Phi = ((X, t) => PhiFunc(X));

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            if(D == 2) {
                C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
                C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 });
            }

            if(D == 3) {
                C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0, (X, t) => 0.0 });
                C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0, (X, t) => 0.0 });
            }

            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionPressure.Add("A", (X, t) => Pjump);
            C.ExactSolutionPressure.Add("B", (X, t) => 0.0);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            C.AddBoundaryValue("wall_left");
            C.AddBoundaryValue("wall_right");
            if(D == 3) {
                C.AddBoundaryValue("wall_front");
                C.AddBoundaryValue("wall_back");
            }


            #endregion


            // misc. solver options
            // ====================
            #region solver
            C.AdaptiveMeshRefinement = false;
            C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            C.BaseRefinementLevel = 0;
            C.RefinementLevel = 0;
            //C.AMR_startUpSweeps = 2;

            //C.InitSignedDistance = false;
            C.adaptiveReInit = false;

            //C.LinearSolver.SolverCode = LinearSolverCode.exp_Kcycle_schwarz;
            C.LinearSolver.SolverCode = LinearSolverCode.automatic;
            //C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.LinearSolver.NoOfMultigridLevels = 1;
            C.NonLinearSolver.MaxSolverIterations = 50;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 80;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-9;
            C.LevelSet_ConvergenceCriterion = 1e-7;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            //C.LinearSolver = new DirectSolver() { WhichSolver = DirectSolver._whichSolver.PARDISO };


            C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.Option_LevelSetEvolution = (steadyInterface) ? LevelSetEvolution.None : LevelSetEvolution.Fourier;
            C.FastMarchingPenaltyTerms = BoSSS.Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.None;

            C.AdvancedDiscretizationOptions.SurfStressTensor = SurfaceSressTensor.Isotropic;
            //C.PhysicalParameters.mu_I = 1.0 * sigma;
            //C.PhysicalParameters.lambda_I = 2.0 * sigma;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.STFstabilization = DoNotTouchParameters.SurfaceTensionForceStabilization.None;



            #endregion

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;


            //C.Timestepper_LevelSetHandling = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            //C.LSunderrelax = 0.05;
            C.Timestepper_LevelSetHandling = (steadyInterface) ? LevelSetHandling.None : LevelSetHandling.Coupled_Once;

            C.TimesteppingMode = compMode;
            //C.CompMode = AppControl._CompMode.Transient; 

            double dt = 0.1; //0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1; // 12500; // (int)(125.0 / dt);
            C.saveperiod = 10;


            return C;

        }

        public static ZLS_Control Test_VerticalBeamInChannel(int p = 2, int kelem = 4, int AMRlvl = 0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            C.ProjectName = "ZLSTestVerticalBeam";
            //C.ProjectDescription = "Vertical Beam in channel";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 0.05;
            C.PhysicalParameters.mu_B = 0.05;
            
            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            C.Material = new Beam();

            #endregion


            // grid generation
            // ===============
            #region grid

            double xLeft = -3;
            double xRight = 8;
            double ySize = 2;

            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(xLeft, xRight, 11 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2*kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "freeslip_upper");
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;

                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] - xLeft) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xRight) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double power = 2;
            double a = 0.1;
            double b = 1.01;

            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                return - Math.Pow((Math.Pow(Math.Abs(X[0] / a), power) + Math.Pow(Math.Abs(X[1] / b), power)), 1.0 / 1.0) + 1;
                //return -((X[0]).Pow2() + (X[1] - 0.5).Pow2()).Sqrt() + 0.3;
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);


            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0);


            double R(double t) {
                if(t < 1) {
                    return (35 + (-84 + (70 - 20 * t) * t) * t) * t * t * t * t;
                } else {
                    return 1;
                }
            }

            double vmax = 0.01;
            double d = 0.2;
            double inflow(double[] x, double t) {
                if(x[1] >= d)
                    return vmax * R(t);
                else
                    return vmax * (x[1] / d) * (x[1] / d) * R(t);
            }

            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("freeslip_upper");
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", inflow);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", inflow);
            C.AddBoundaryValue("pressure_outlet_right");

            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 80;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-10;
            C.LinearSolver.ConvergenceCriterion = 1e-10;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = 2 });
            C.AMR_startUpSweeps = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.TimesteppingMode = compMode;
            double dt = 5e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static ZLS_Control BallInChannel(int p = 2, int kelem = 3, int AMRlvl = 0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            C.ProjectName = "ZLSTestVerticalBeam";
            //C.ProjectDescription = "Vertical Beam in channel";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 0.05;
            C.PhysicalParameters.mu_B = 0.05;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            C.Material = new HardSiliconeRubber();

            #endregion


            // grid generation
            // ===============
            #region grid

            double xLeft = -3;
            double xRight = 8;
            double ySize = 2;

            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(xLeft, xRight, 11 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2 * kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "freeslip_lower");
                grd.EdgeTagNames.Add(2, "freeslip_upper");
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;

                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] - xLeft) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xRight) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double power = 2;
            double a = 0.1;
            double b = 1.01;

            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                //return -Math.Pow((Math.Pow(Math.Abs(X[0] / a), power) + Math.Pow(Math.Abs(X[1] / b), power)), 1.0 / 1.0) + 1;
                return -((X[0]).Pow2() + (X[1] - 1).Pow2()).Sqrt() + 0.3;
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            double v0 = 0.1;
            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0);
            C.InitialValues_Evaluators.Add("VelocityX#C", X => 0);
            #endregion


            // boundary conditions
            // ===================
            #region BC
            double R(double t) {
                if (t < 1) {
                    return (35 + (-84 + (70 - 20 * t) * t) * t) * t * t * t * t;
                } else {
                    return 1;
                }
            }

            double vmax = 0.01;
            double inflow(double[] x, double t) {
                    return vmax * R(t);
            }


            C.AddBoundaryValue("freeslip_lower");
            C.AddBoundaryValue("freeslip_upper");
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", inflow);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", inflow);
            C.AddBoundaryValue("pressure_outlet_right");

            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 80;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.MaxSolverIterations = 50;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-10;
            C.LinearSolver.ConvergenceCriterion = 1e-10;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;



            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = 2 });
            C.AMR_startUpSweeps = 1;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;

            C.TimesteppingMode = compMode;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion

            return C;
        }
    }
}

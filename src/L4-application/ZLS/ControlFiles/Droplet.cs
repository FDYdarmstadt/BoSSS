using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
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

namespace ZwoLevelSetSolver.ControlFiles {

    public static class Droplet {

        public static ZLS_Control Aland( int p = 1, int AMRlvl = 0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;

            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 3;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            //C.DbPath = @"C:\Users\miao\Documents\Database\Droplet-EE";
            C.ProjectName = "Droplet";
            C.SessionName = "Droplet-check";
            C.ProjectDescription = "Droplet running on pc";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            //if(D == 3) {
            //    C.FieldOptions.Add("DisplacementZ", new FieldOpts() {
            //        Degree = p,
            //        SaveToDB = FieldOpts.SaveToDBOpt.FALSE
            //    });
            //}

            // Physical Parameters
            // ===================
            #region physics

            double scale = 1e-4;
            C.PhysicalParameters.rho_A = 1 * scale * scale * scale;
            C.PhysicalParameters.rho_B = 1260 * scale * scale * scale;
            C.PhysicalParameters.mu_A = 0.1 * scale ;
            C.PhysicalParameters.mu_B = 1.41 * scale;
            double sigma = 0.046;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid() {
                Density = 1000 * scale * scale * scale,
                Lame2 = 1000 * scale,
                Viscosity = 1 * scale
                //Viscosity = 0 0.05e-4 * scale,
            };

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 6;
            double yTop = (3.5 - 0.00625);
            double yBottom = (- 2.5 - 0.00625);
            int kelem = 12;
            //int kelem = 36;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem / 2 + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                //grd.EdgeTagNames.Add(3, "wall_left");
                //grd.EdgeTagNames.Add(4, "wall_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    //if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                    //    et = 3;
                    //if(Math.Abs(X[0] - xSize) <= 1.0e-8)
                    //   et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = (1.778);
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
            Func<double[], double> Phi1Func = (X => -( X[1] - 0.001));
            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("pressure_outlet_upper");
            //C.AddBoundaryValue("wall_left");
            //C.AddBoundaryValue("wall_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            C.PhysicalParameters.sliplength = 0.001;

            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 1;
            C.NonLinearSolver.MinSolverIterations = 0;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 0;// 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;
            C.NonLinearSolver.Globalization = BoSSS.Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;


            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;
            //C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = 6 });
            //C.AMR_startUpSweeps = 5;
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = 3 });
            //C.AMR_startUpSweeps = 2;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            

            C.TimesteppingMode = compMode;
            //double dt = 5e-7;
            double dt = 5e-4;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 10;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 1;

            #endregion

            return C;
        }


        public static ZLS_Control Wiki(int p = 4, int AMRlvl = 2, double SlipLength = 0, bool useGravity = true) {
            // --control "cs:ZwoLevelSetSolver.ControlFiles.Droplet.Wiki()"
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;

            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 3;

            //AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Steady;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            //C.DbPath = @"C:\Users\miao\Documents\Database\Droplet-EE";
            C.ProjectName = "Droplet";
            C.SessionName = "Droplet-check";
            C.ProjectDescription = "Droplet running on pc";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            //if(D == 3) {
            //    C.FieldOptions.Add("DisplacementZ", new FieldOpts() {
            //        Degree = p,
            //        SaveToDB = FieldOpts.SaveToDBOpt.FALSE
            //    });
            //}

            // Physical Parameters
            // ===================
            #region physics

            double scale = 1;// 4.4175e-4; // For a droplet with radius r = 176.7 micrometre
            //double SL = 0;// 5e-6 / scale;

            C.PhysicalParameters.rho_A = 1 * scale * scale * scale; //Liquid
            C.PhysicalParameters.rho_B = 1 * scale * scale * scale; //Gas
            //C.PhysicalParameters.rho_B = 10 * scale * scale * scale; //Gas
            C.PhysicalParameters.mu_A = 1.41 * scale;
            C.PhysicalParameters.mu_B = 1.41 * scale;
            //C.PhysicalParameters.mu_B = 0.1 * scale;
            double sigma = 1.0 * (0.2 * 0.2 + 0.5 * 0.1 * 0.1 * Math.PI) * 9.8;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = C.PhysicalParameters.mu_A / SlipLength;
            C.PhysicalParameters.betaS_B = C.PhysicalParameters.mu_B / SlipLength;

            //C.PhysicalParameters.betaL = (C.PhysicalParameters.mu_A + C.PhysicalParameters.mu_B) / SL / 2.0;
            C.PhysicalParameters.betaL = 0;
            C.PhysicalParameters.theta_e = 0.5 * Math.PI;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            C.Material = new Solid() {
                Density = 1 * scale * scale * scale,
                Lame2 = 1e7 * scale,
                Viscosity = 1.41 * scale
            };

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 1;
            double yTop = 1;
            double yBottom = -0.5;
            int kelem = 8;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0.0, xSize, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem / 2 * 3 + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                //grd.EdgeTagNames.Add(1, "NavierSlip_Linear_lower");
                grd.EdgeTagNames.Add(1, "wall_lower");
                //grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(2, "freeslip_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + 0.0) <= 1.0e-8)
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

            double R = 0.4;
            //double Theta_e = Math.PI / 2;
            //double s = 2 * R * Math.Sin(Theta_e);
            //double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            //double[] center = new double[] { 0.0, -1.5 };

            Func<double[], double> PhiFunc = Phi;

            double Phi(double[] X) {
                double x = X[0];
                double y = X[1];
                double val;
                double R = 0.05;
                if(y >= 0.2)
                    val = x * x - R * R;
                else
                    val = x * x + (y - 0.2) * (y - 0.2) - R * R;
                return val;
            }


            double pJump = sigma / R;

            //C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            //C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);
            //C.AddInitialValue("GravityY#A", new Formula($"X => -9.8", false));
            //C.AddInitialValue("GravityY#B", new Formula($"X => -9.8", false));
            if(useGravity)
                C.AddInitialValue("GravityY#C", new Formula($"X => -9.8", false));


            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            Func<double[], double> Phi1Func = Phi1;

            double Phi1(double[] X) {
                double x = X[0];
                double y = X[1];
                double val;
                double R = 0.1;
                if(x <= 2 * R)
                    val = (y - 0.2) * (y - 0.2) - R * R;
                else
                    val = (x - 2 * R) * (x - 2 * R) + (y - 0.2) * (y - 0.2) - R * R;
                return -val;
            }
            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("wall_lower");
            //C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("freeslip_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            if(SlipLength != 0.0) {
                C.AdvancedDiscretizationOptions.IBM_BoundaryType = IBM_BoundaryType.NavierSlip;
            } else {
                C.AdvancedDiscretizationOptions.IBM_BoundaryType = IBM_BoundaryType.NoSlip;
            }
            //C.AdvancedDiscretizationOptions.IBM_BoundaryType = IBM_BoundaryType.NavierSlip_Linear;
            //C.PhysicalParameters.sliplength = 0.001;

            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MinSolverIterations = 2;
            //C.Solver_MaxIterations = 50;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;
            C.NonLinearSolver.Globalization = BoSSS.Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;


            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 0, maxRefinementLevel = AMRlvl - 2});
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 1, maxRefinementLevel = AMRlvl });
            //C.AMR_startUpSweeps = AMRlvl;

            #endregion

            C.DynamicLoadBalancing_On = false;
            C.DynamicLoadBalancing_Period = 1000;
            C.DynamicLoadBalancing_RedistributeAtStartup = false;
            C.GridPartType = GridPartType.METIS;

            //C.ReInitPeriod = 10;
            C.FailOnSolverFail = false;

            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 1e-1;
            C.dtMax = dt;
            C.dtMin = dt;
            //C.Endtime = 100;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion

            return C;
        }



        public static ZLS_Control OnDB(int p = 2, int AMRlvl = 2) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;

            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 3;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            //C.DbPath = @"C:\Users\miao\Documents\Database\Droplet-EE";
            C.ProjectName = "Droplet";
            C.SessionName = "Droplet-check";
            C.ProjectDescription = "Droplet running on pc";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            //if(D == 3) {
            //    C.FieldOptions.Add("DisplacementZ", new FieldOpts() {
            //        Degree = p,
            //        SaveToDB = FieldOpts.SaveToDBOpt.FALSE
            //    });
            //}

            // Physical Parameters
            // ===================
            #region physics

            double scale = 1e-4;
            C.PhysicalParameters.rho_A = 1000.0 * scale * scale * scale;
            C.PhysicalParameters.rho_B = C.PhysicalParameters.rho_A;// 1260 * scale * scale * scale;
            C.PhysicalParameters.mu_A = 1000.0 * scale;
            C.PhysicalParameters.mu_B = C.PhysicalParameters.mu_A;// 1.41 * scale;
            double sigma = 0.1; //0.046;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 1e-12;
            C.PhysicalParameters.betaS_B = 1e-12;

            C.PhysicalParameters.betaL = 1e-12;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid() {
                Density = 1000 * scale * scale * scale,
                Lame2 = 1e5 * scale,
                Viscosity = 1 * scale
                //Viscosity = 0 0.05e-4 * scale,
            };

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 6;
            double yTop = 3;//(3.5 - 0.00625);
            double yBottom = -3;// (-2.5 - 0.00625);
            //int kelem = 12;
            int kelem = 5;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, 0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                //grd.EdgeTagNames.Add(1, "wall_lower");
                //grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(1, "NavierSlip_Linear_lower");
                grd.EdgeTagNames.Add(2, "freeslip_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 1.778;
            double Theta_e = Math.PI / 2;
            double s = 2 * R * Math.Sin(Theta_e);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0.0, -2.0 };

            Func<double[], double> PhiFunc = Phi;

            double Phi(double[] X) {
                //return 3 *  (((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2());
                return (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - 4.0 * R * R;
            }


            double pJump = sigma / R;

            C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);


            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            Func<double[], double> Phi1Func = (X => -(X[1] + 2.0));
            //C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            //C.AddBoundaryValue("wall_lower");
            //C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("NavierSlip_Linear_lower");
            C.AddBoundaryValue("freeslip_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

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

            C.NonLinearSolver.MaxSolverIterations = 160;
            C.NonLinearSolver.MinSolverIterations = 2;
            //C.Solver_MaxIterations = 50;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;
            C.NonLinearSolver.Globalization = BoSSS.Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;


            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = AMRlvl });
            //C.AMR_startUpSweeps = AMRlvl;
            //C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = 6 });
            //C.AMR_startUpSweeps = 5;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 0, maxRefinementLevel = AMRlvl });
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 1, maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


            C.TimesteppingMode = compMode;
            double dt = 1;
            //double dt = 5e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static ZLS_Control OnIBM(int p = 2, int AMRlvl = 2) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;

            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 3;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            //C.DbPath = @"C:\Users\miao\Documents\Database\Droplet-EE";
            C.ProjectName = "Droplet";
            C.SessionName = "Droplet-check";
            C.ProjectDescription = "Droplet running on pc";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            //if(D == 3) {
            //    C.FieldOptions.Add("DisplacementZ", new FieldOpts() {
            //        Degree = p,
            //        SaveToDB = FieldOpts.SaveToDBOpt.FALSE
            //    });
            //}

            // Physical Parameters
            // ===================
            #region physics

            double scale = 1e-4;
            C.PhysicalParameters.rho_A = 10.0 * scale * scale * scale;
            C.PhysicalParameters.rho_B = 1000.0 * scale * scale * scale;
            C.PhysicalParameters.mu_A = 100.0 * scale;
            C.PhysicalParameters.mu_B = 1000.0 * scale;
            double sigma = 1; //0.1; //0.046;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.00001;
            C.PhysicalParameters.betaS_B = 0.00001;

            C.PhysicalParameters.betaL = 0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid() {
                Density = 1000 * scale * scale * scale,
                Lame2 = 1e3 * scale,
                Viscosity = 1000 * scale
                //Viscosity = 0 0.05e-4 * scale,
            };

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 6 * 2;
            double yTop = 3 * 2;//(3.5 - 0.00625);
            double yBottom = -3 * 2;// (-2.5 - 0.00625);
            int kelem = 5;
            //int kelem = 5 * 8;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, 0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 1.778;
            double Theta_e = Math.PI / 2;
            double s = 2 * R * Math.Sin(Theta_e);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0.0, -1.5 };

            Func<double[], double> PhiFunc = Phi;

            //double Phi(double[] X) {
            //    //return 3 *  (((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2());
            //    return -((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - R * R);
            //}

            //double Phi(double[] X) {
            //    //return 3 *  (((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2());
            //    return -(X[0] + 3.0 - 1e-5);
            //}

            double Phi(double[] X) {
                double x = X[0];
                double y = X[1];
                double dx = x - center[0];
                double dy = y - center[1];

                if(y > -2) {
                    return -((dx * dx + dy * dy) - R * R);
                } else {
                    double dyEdge = -2 - center[1];
                    return -((dx * dx + dyEdge * dyEdge) - R * R);
                }
            }


            double pJump = sigma / R;

            C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);
            //C.AddInitialValue("GravityY#A", new Formula($"X => -9.8", false));
            //C.AddInitialValue("GravityY#B", new Formula($"X => -9.8", false));
            //C.AddInitialValue("GravityY#C", new Formula($"X => -9.8",false));


            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            Func<double[], double> Phi1Func = (X => -(X[1] + 2.0));
            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            C.AdvancedDiscretizationOptions.IBM_BoundaryType = IBM_BoundaryType.NavierSlip;
            C.PhysicalParameters.sliplength = 0.001;

            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 160;
            C.NonLinearSolver.MinSolverIterations = 2;
            //C.Solver_MaxIterations = 50;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;
            C.NonLinearSolver.Globalization = BoSSS.Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;


            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 0, maxRefinementLevel = AMRlvl });
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 1, maxRefinementLevel = AMRlvl });
            //C.AMR_startUpSweeps = AMRlvl;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


            C.TimesteppingMode = compMode;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static ZLS_Control OnIBMBug(int p = 2, int AMRlvl = 3) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;

            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 3;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"\\fdyprime\userspace\smuda\cluster\cluster_db";
            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            //C.DbPath = @"C:\Users\miao\Documents\Database\Droplet-EE";
            C.ProjectName = "Droplet";
            C.SessionName = "Droplet-check";
            C.ProjectDescription = "Droplet running on pc";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion

            //if(D == 3) {
            //    C.FieldOptions.Add("DisplacementZ", new FieldOpts() {
            //        Degree = p,
            //        SaveToDB = FieldOpts.SaveToDBOpt.FALSE
            //    });
            //}

            // Physical Parameters
            // ===================
            #region physics

            double scale = 1e-4;
            C.PhysicalParameters.rho_A = 10.0 * scale * scale * scale;
            C.PhysicalParameters.rho_B = 1000.0 * scale * scale * scale;
            C.PhysicalParameters.mu_A = 100.0 * scale;
            C.PhysicalParameters.mu_B = 1000.0 * scale;
            double sigma = 0.1; //0.046;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.00001;
            C.PhysicalParameters.betaS_B = 0.00001;

            C.PhysicalParameters.betaL = 0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid() {
                Density = 1000 * scale * scale * scale,
                Lame2 = 1e3 * scale,
                Viscosity = 1000 * scale
                //Viscosity = 0 0.05e-4 * scale,
            };

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 6 * 2;
            double yTop = 3 * 2;//(3.5 - 0.00625);
            double yBottom = -3 * 2;// (-2.5 - 0.00625);
            int kelem = 5;
            //int kelem = 5 * 8;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-xSize, 0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "pressure_outlet_upper");
                grd.EdgeTagNames.Add(3, "freeslip_left");
                grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] + xSize) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - 0.0) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // Initial Values
            // ==============
            #region init

            double R = 1.778;
            double Theta_e = Math.PI / 2;
            double s = 2 * R * Math.Sin(Theta_e);
            double h = Math.Sqrt(R.Pow2() - (0.25 * s.Pow2()));

            double[] center = new double[] { 0.0, -2.0 }; //bug of quad?
            //double[] center = new double[] { 0.0, -1.5 };

            Func<double[], double> PhiFunc = Phi;

            double Phi(double[] X) {
                //return 3 *  (((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2());
                return -((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - R * R);
            }

            //double Phi(double[] X) {
            //    //return 3 *  (((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2());
            //    return -(X[0] + 3.0 - 1e-5);
            //}


            double pJump = sigma / R;

            C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);
            //C.AddInitialValue("GravityY#A", new Formula($"X => -9.8", false));
            //C.AddInitialValue("GravityY#B", new Formula($"X => -9.8", false));
            //C.AddInitialValue("GravityY#C", new Formula($"X => -9.8",false));


            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            Func<double[], double> Phi1Func = (X => -(X[1] + 2.0));
            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("pressure_outlet_upper");
            C.AddBoundaryValue("freeslip_left");
            C.AddBoundaryValue("freeslip_right");

            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            C.AdvancedDiscretizationOptions.IBM_BoundaryType = IBM_BoundaryType.NavierSlip;
            C.PhysicalParameters.sliplength = 0.001;

            #endregion

            // misc. solver options
            // ====================
            #region solver


            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 160;
            C.NonLinearSolver.MinSolverIterations = 2;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;
            C.NonLinearSolver.Globalization = BoSSS.Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;


            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.FullySymmetric;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 0, maxRefinementLevel = AMRlvl });
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { levelSet = 1, maxRefinementLevel = AMRlvl });
            //C.AMR_startUpSweeps = AMRlvl;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


            C.TimesteppingMode = compMode;
            double dt = 1e-3;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 10000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static ZLS_Control DropletOnPlate(int p = 2, int kelem = 6, int AMRlvl = 0) {
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
            double sigma = 0.5;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new HardSiliconeRubber();

            #endregion


            // grid generation
            // ===============
            #region grid

            double xSize = 0.3;
            double yTop = 0.3 - 0.00625;
            double yBottom = -0.3 - 0.00625;

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

            double[] center = new double[] { 0, -0.28 };

            Func<double[], double> PhiFunc = Phi;

            double Phi(double[] X) {
                //return 3 *  (((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - R.Pow2());
                return (X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2() - R * R;
            }


            double pJump = sigma / R;

            C.InitialValues_Evaluators.Add("Pressure#A", X => pJump);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);


            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            //Func<double[], double> Phi1Func = (X => -(X[1] -0.02 + 0.4 * X[0] * X[0]));
            Func<double[], double> Phi1Func = (X => -(X[1] + 0.26));
            //C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

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
            //C.Solver_MaxIterations = 50;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-10;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;



            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new ContactPointRefiner { maxRefinementLevel = 5 });
            C.AMR_startUpSweeps = 4;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;


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
    }
}

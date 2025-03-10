using BoSSS.Application.XNSE_Solver;
using BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases.PrintingNip;
using BoSSS.Application.XNSE_Solver.Tests;
using HFSISolver.Tests;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using HFSISolver.SolidPhase;
using static BoSSS.Solution.Control.AppControl;
using MathNet.Numerics.Distributions;
using System.Drawing;

namespace HFSISolver.Tests {
    public static class FSTC {

        public static HFSI_Control Circle(int p = 2, int h = 3, int AMRlvl = 0, double BeamDensity = 0.1) {
            HFSI_Control C = new HFSI_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.2;
            C.NoOfMultigridLevels = 1;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            // basic database options
            // ======================
            #region db
            //string DbPath = @"C:\Databases\ZLS_GummiLippe";
            C.savetodb = false;
            //C.DbPath = DbPath;
            C.ProjectName = "Fluid-Solid-Taylor-Couette";
            //C.ProjectDescription = "Fluid-Solid-Taylor-Couette";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            // restart
            //Guid restartID = new Guid("1e7248ea-dd9e-49ee-b7c7-ee33c61c2dda");
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restartID, "83");

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1.3;
            //C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.2;
            //C.PhysicalParameters.mu_B = 0.05;
            C.PhysicalParameters.Sigma = 0.9;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_C = BeamDensity;

            C.ThermalParameters.c_A = 1;
            C.ThermalParameters.c_C = 1;

            C.ThermalParameters.k_A = 3;
            C.ThermalParameters.k_C = 40;

            C.ThermalParameters.hVap = 0;
            C.ThermalParameters.T_sat = 100;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid {
                Viscosity = 0.2,
                Lame2 = 10,
                Density = BeamDensity
            };
            #endregion

            // grid generation
            // ===============
            #region grid

            double xLeft = -2;
            double xRight = 2;
            double yTop = 2;
            double yBottom = -2;
            int kelem = Convert.ToInt32(Math.Pow(2, h));

            double[] CutOut1Point1 = new double[2] { 0.5, -0.5 };
            double[] CutOut1Point2 = new double[2] { -0.5, 0.5 };

            var CutOut1 = new BoSSS.Platform.Utils.Geom.BoundingBox(2);
            CutOut1.AddPoint(CutOut1Point1);
            CutOut1.AddPoint(CutOut1Point2);

            //*
            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(xLeft, xRight, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, CutOuts: CutOut1);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_lower");
                grd.EdgeTagNames.Add(2, "wall_ConstantTemperature_upper");
                grd.EdgeTagNames.Add(3, "wall_ConstantTemperature_left");
                grd.EdgeTagNames.Add(4, "wall_ConstantTemperature_right");
                grd.EdgeTagNames.Add(5, "wall_ConstantTemperature_inside");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 5;

                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] - xLeft) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xRight) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };
            //*/


            #endregion

            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                return -(X[1] * X[1] + X[0] * X[0] - (1 + 0.25 * Math.Sqrt(2)) * (1 + 0.25 * Math.Sqrt(2)));
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            double A1 = 0.922554663967561;
            double B1 = -1.690218655870246;
            double A2 = 0.067608746234810;
            double B2 = -0.033804373117405;

            double InletVelocityX(double[] X, double t, double A, double B) {
                return -(A * Math.Sqrt(X[0] * X[0] + X[1] * X[1]) + B / Math.Sqrt(X[0] * X[0] + X[1] * X[1])) * X[1] / Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
            }
            double InletVelocityY(double[] X, double t, double A, double B) {
                return (A * Math.Sqrt(X[0] * X[0] + X[1] * X[1]) + B / Math.Sqrt(X[0] * X[0] + X[1] * X[1])) * X[0] / Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
            }

            double PreDisplacementX(double[] X, double t, double A, double B) {
                if(t < 5) {
                    t = 0;
                } else if(t < 15) {
                    t = (t - 5) * 0.1;
                } else { 
                    t = 1; }
                return t * (-(A * Math.Sqrt(X[0] * X[0] + X[1] * X[1]) + B / Math.Sqrt(X[0] * X[0] + X[1] * X[1])) * X[1] / Math.Sqrt(X[0] * X[0] + X[1] * X[1]));
            }
            double PreDisplacementY(double[] X, double t, double A, double B) {
                if(t < 5) {
                    t = 0;
                } else if(t < 15) {
                    t = (t - 5) * 0.1;
                } else {
                    t = 1;
                }
                return t * ((A * Math.Sqrt(X[0] * X[0] + X[1] * X[1]) + B / Math.Sqrt(X[0] * X[0] + X[1] * X[1])) * X[0] / Math.Sqrt(X[0] * X[0] + X[1] * X[1]));
            }

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => InletVelocityX(X));
            //C.InitialValues_Evaluators.Add("VelocityY#A", X => InletVelocityY(X));
            //C.InitialValues_Evaluators.Add("DisplacementX#C", X => PreDisplacementX(X));
            //C.InitialValues_Evaluators.Add("DisplacementY#C", X => PreDisplacementY(X));
            //C.InitialValues_Evaluators.Add("Pressure#C", X => 0);
            //C.InitialValues_Evaluators.Add("Temperature#A", X => 5);
            //C.InitialValues_Evaluators.Add("Temperature#C", X => 1);

            #endregion

            // boundary conditions
            // ===================
            #region BC

            double T_out = 5;
            double T_in = 3;
            C.AddBoundaryValue("wall_ConstantTemperature_lower", "VelocityX", (X, t) => InletVelocityX(X, t, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_lower", "VelocityY", (X, t) => InletVelocityY(X, t, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#A", "X => " + T_out, false);

            C.AddBoundaryValue("wall_ConstantTemperature_upper", "VelocityX", (X, t) => InletVelocityX(X, t, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_upper", "VelocityY", (X, t) => InletVelocityY(X, t, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_upper", "Temperature#A", "X => " + T_out, false);

            C.AddBoundaryValue("wall_ConstantTemperature_left", "VelocityX", (X, t) => InletVelocityX(X, t, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_left", "VelocityY", (X, t) => InletVelocityY(X, t, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#A", "X => " + T_out, false);

            C.AddBoundaryValue("wall_ConstantTemperature_right", "VelocityX", (X, t) => InletVelocityX(X, t, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_right", "VelocityY", (X, t) => InletVelocityY(X, t, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#A", "X => " + T_out, false);

            C.AddBoundaryValue("wall_ConstantTemperature_inside", "DisplacementX", (X, t) => PreDisplacementX(X, t, A2, B2));
            C.AddBoundaryValue("wall_ConstantTemperature_inside", "DisplacementY", (X, t) => PreDisplacementY(X, t, A2, B2));
            C.AddBoundaryValue("wall_ConstantTemperature_inside", "Temperature#C", "X => " + T_in, false);

            #endregion

            // misc. solver options
            // ====================
            #region solver
            C.AdvancedDiscretizationOptions.GNBC_Localization = NavierSlip_Localization.Bulk;
            C.AdvancedDiscretizationOptions.GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;
            C.PhysicalParameters.sliplength = 0.001;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 80;
            C.NonLinearSolver.MinSolverIterations = 2;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;
            C.DynamicLoadBalancing_On = false;
            C.DynamicLoadBalancing_RedistributeAtStartup = false;
            C.GridPartType = GridPartType.METIS;
            #endregion

            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.TimesteppingMode = compMode;
            double dt = 1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.NoOfTimesteps = 20;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static HFSI_Control SmallCircle(IHFSITest tst, int p = 2, int h = 3, int AMRlvl = 0, double BeamDensity = 0.1) {
            HFSI_Control C = new HFSI_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.2;
            C.NoOfMultigridLevels = 1;

            // basic database options
            // ======================
            #region db
            //string DbPath = @"C:\Databases\ZLS_GummiLippe";
            C.savetodb = false;
            //C.DbPath = DbPath;
            C.ProjectName = "Fluid-Solid-Taylor-Couette";
            //C.ProjectDescription = "Fluid-Solid-Taylor-Couette";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            // restart
            //Guid restartID = new Guid("1e7248ea-dd9e-49ee-b7c7-ee33c61c2dda");
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restartID, "83");

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1.3;
            //C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.2;
            //C.PhysicalParameters.mu_B = 0.05;
            C.PhysicalParameters.Sigma = 0.9;

            C.ThermalParameters.rho_A = C.PhysicalParameters.rho_A;
            C.ThermalParameters.rho_C = BeamDensity;

            C.ThermalParameters.c_A = 1;
            C.ThermalParameters.c_C = 1;

            C.ThermalParameters.k_A = 3;
            C.ThermalParameters.k_C = 40;

            C.ThermalParameters.hVap = 0;
            C.ThermalParameters.T_sat = 100;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid {
                Viscosity = 0.2,
                Lame2 = 10,
                Density = BeamDensity
            };
            #endregion

            // grid generation
            // ===============
            #region grid

            double xLeft = -2;
            double xRight = 2;
            double yTop = 2;
            double yBottom = -2;
            int kelem = Convert.ToInt32(Math.Pow(2, h));

            double[] CutOut1Point1 = new double[2] { 0.5, -0.5 };
            double[] CutOut1Point2 = new double[2] { -0.5, 0.5 };

            var CutOut1 = new BoSSS.Platform.Utils.Geom.BoundingBox(2);
            CutOut1.AddPoint(CutOut1Point1);
            CutOut1.AddPoint(CutOut1Point2);

            //*
            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(xLeft, xRight, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, CutOuts: CutOut1);

                grd.EdgeTagNames.Add(1, "wall_ConstantTemperature_lower");
                grd.EdgeTagNames.Add(2, "wall_ConstantTemperature_upper");
                grd.EdgeTagNames.Add(3, "wall_ConstantTemperature_left");
                grd.EdgeTagNames.Add(4, "wall_ConstantTemperature_right");
                grd.EdgeTagNames.Add(5, "wall_ConstantTemperature_inside");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 5;

                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] - xLeft) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xRight) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };
            //*/


            #endregion

            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                return -(X[1] * X[1] + X[0] * X[0] - (1 + 0.25 * Math.Sqrt(2)) * (1 + 0.25 * Math.Sqrt(2)));
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            //Displacement equation for region 1(R1 <= r <= R2):
            //u_theta1(r) = 0.135217492469620 * r + -0.033804373117405 / r
            //Velocity equation for region 2(R2 <= r <= R3):
            //v_theta2(r) = 0.922554663967561 * r + -1.690218655870246 / r
            double A1 = 0.922554663967561;
            double B1 = -1.690218655870246;
            double A2 = 0.135217492469620;
            double B2 = -0.033804373117405;
            //The temperature distribution in region 1, T1(r), is:
            //T1(r) = 0.645015726004179 * log(r) + 1.447090831896623
            //The temperature distribution in region 2, T2(r), is:
            //T2(r) = 8.600209680055716 * log(r) + -0.961211091954969
            double A3 = 0.645015726004179;
            double B3 = 1.447090831896623;
            double A4 = 8.600209680055716;
            double B4 = -0.961211091954969;


            double InletVelocityX(double[] X, double A, double B) {
                return -(A * Math.Sqrt(X[0] * X[0] + X[1] * X[1]) + B / Math.Sqrt(X[0] * X[0] + X[1] * X[1])) * X[1] / Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
            }
            double InletVelocityY(double[] X, double A, double B) {
                return (A * Math.Sqrt(X[0] * X[0] + X[1] * X[1]) + B / Math.Sqrt(X[0] * X[0] + X[1] * X[1])) * X[0] / Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
            }
            double Temperature(double[] X, double A, double B) {
                return A * Math.Log(Math.Sqrt(X[0] * X[0] + X[1] * X[1])) + B;
            }

            //C.InitialValues_Evaluators.Add("VelocityX#A", X => InletVelocityX(X, A1, B1));
            //C.InitialValues_Evaluators.Add("VelocityY#A", X => InletVelocityY(X, A1, B1));
            //C.InitialValues_Evaluators.Add("DisplacementX#C", X => InletVelocityX(X, A2, B2));
            //C.InitialValues_Evaluators.Add("DisplacementY#C", X => InletVelocityY(X, A2, B2));
            //C.InitialValues_Evaluators.Add("Pressure#C", X => 0);


            #endregion

            // boundary conditions
            // ===================
            #region BC
            C.AddBoundaryValue("wall_ConstantTemperature_lower", "VelocityX", X => InletVelocityX(X, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_lower", "VelocityY", X => InletVelocityY(X, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_lower", "Temperature#A", X => Temperature(X, A4, B4));

            C.AddBoundaryValue("wall_ConstantTemperature_upper", "VelocityX", X => InletVelocityX(X, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_upper", "VelocityY", X => InletVelocityY(X, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_upper", "Temperature#A", X => Temperature(X, A4, B4));

            C.AddBoundaryValue("wall_ConstantTemperature_left", "VelocityX", X => InletVelocityX(X, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_left", "VelocityY", X => InletVelocityY(X, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_left", "Temperature#A", X => Temperature(X, A4, B4));

            C.AddBoundaryValue("wall_ConstantTemperature_right", "VelocityX", X => InletVelocityX(X, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_right", "VelocityY", X => InletVelocityY(X, A1, B1));
            C.AddBoundaryValue("wall_ConstantTemperature_right", "Temperature#A", X => Temperature(X, A4, B4));

            C.AddBoundaryValue("wall_ConstantTemperature_inside", "DisplacementX", X => InletVelocityX(X, A2, B2));
            C.AddBoundaryValue("wall_ConstantTemperature_inside", "DisplacementY", X => InletVelocityY(X, A2, B2));
            C.AddBoundaryValue("wall_ConstantTemperature_inside", "Temperature#C", X => Temperature(X, A3, B3));


            // initial values and exact solution
            // =================================

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionDisplacement = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionTemperature = new Dictionary<string, Func<double[], double, double>>();


            int D = tst.SpatialDimension;

            foreach(var spc in new[] { "A", "B", "C" }) {
                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));
                C.ExactSolutionDisplacement.Add(spc, D.ForLoop(d => tst.GetDis(spc, d)));
                C.ExactSolutionTemperature.Add(spc, tst.GetTemperature(spc));

                for(int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d) + "#" + spc, tst.GetU(spc, d).Convert_Xt2X(0.0));
                    C.InitialValues_Evaluators.Add(VariableNames.DisplacementComponent(d) + "#" + spc, tst.GetDis(spc, d).Convert_Xt2X(0.0));
                    //var Gravity_d = tst.GetF(spc, d).Convert_X2Xt();
                    //C.SetGravity(spc, d, Gravity_d);
                }

                C.InitialValues_Evaluators.Add(BoSSS.Solution.NSECommon.VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));
                C.InitialValues_Evaluators.Add(BoSSS.Solution.NSECommon.VariableNames.Temperature + "#" + spc, tst.GetTemperature(spc).Convert_Xt2X(0.0));
            }

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
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;
            C.DynamicLoadBalancing_On = false;
            C.DynamicLoadBalancing_RedistributeAtStartup = false;
            C.GridPartType = GridPartType.METIS;
            #endregion

            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.NoOfTimesteps = 20;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static HFSI_Control CircleCircle(int p = 2, int h = 3, int AMRlvl = 0, double BeamDensity = 0.1) {
            HFSI_Control C = new HFSI_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.2;
            C.NoOfMultigridLevels = 1;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            // basic database options
            // ======================
            #region db
            //string DbPath = @"C:\Databases\ZLS_GummiLippe";
            C.savetodb = false;
            //C.DbPath = DbPath;
            C.ProjectName = "Fluid-Solid-Taylor-Couette";
            //C.ProjectDescription = "Fluid-Solid-Taylor-Couette";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            // restart
            //Guid restartID = new Guid("1e7248ea-dd9e-49ee-b7c7-ee33c61c2dda");
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(restartID, "83");

            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1.3;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.2;
            C.PhysicalParameters.mu_B = 0.05;
            C.PhysicalParameters.Sigma = 0.9;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid {
                Viscosity = 0.2,
                Lame2 = 10,
                Density = BeamDensity
            };
            #endregion

            // grid generation
            // ===============
            #region grid

            double xLeft = -3;
            double xRight = 3;
            double yTop = 3;
            double yBottom = -3;
            int kelem = Convert.ToInt32(Math.Pow(2, h));

            double[] CutOut1Point1 = new double[2] { 0.5, -0.5 };
            double[] CutOut1Point2 = new double[2] { -0.5, 0.5 };

            var CutOut1 = new BoSSS.Platform.Utils.Geom.BoundingBox(2);
            CutOut1.AddPoint(CutOut1Point1);
            CutOut1.AddPoint(CutOut1Point2);

            //*
            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(xLeft, xRight, 3 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(yBottom, yTop, 3 * kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, CutOuts: CutOut1);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");
                grd.EdgeTagNames.Add(3, "wall_left");
                grd.EdgeTagNames.Add(4, "wall_right");
                grd.EdgeTagNames.Add(5, "wall_inside");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 5;

                    if(Math.Abs(X[1] - yBottom) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - yTop) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0] - xLeft) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - xRight) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };
            //*/


            #endregion

            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = delegate (double[] X) {
                return (X[1] * X[1] + X[0] * X[0] - 2 * 2);
            };

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                return -(X[1] * X[1] + X[0] * X[0] - (1 + 0.25 * Math.Sqrt(2)) * (1 + 0.25 * Math.Sqrt(2)));
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            double A1 = 0.005872850876260;
            double B1 = -0.002936425438130;
            double A2 = 0.320551778780949;
            double B2 = -0.587285087625986;
            double A3 = 0.761015594500438;
            double B3 = -2.349140350503944;

            double InletVelocityX(double[] X, double A, double B) {
                return -(A * Math.Sqrt(X[0] * X[0] + X[1] * X[1]) + B / Math.Sqrt(X[0] * X[0] + X[1] * X[1])) * X[1] / Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
            }
            double InletVelocityY(double[] X, double A, double B) {
                return (A * Math.Sqrt(X[0] * X[0] + X[1] * X[1]) + B / Math.Sqrt(X[0] * X[0] + X[1] * X[1])) * X[0] / Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
            }


            C.InitialValues_Evaluators.Add("VelocityX#A", X => InletVelocityX(X, A2, B2));
            C.InitialValues_Evaluators.Add("VelocityY#A", X => InletVelocityY(X, A2, B2));
            C.InitialValues_Evaluators.Add("VelocityX#B", X => InletVelocityX(X, A3, B3));
            C.InitialValues_Evaluators.Add("VelocityY#B", X => InletVelocityY(X, A3, B3));
            C.InitialValues_Evaluators.Add("DisplacementX#C", X => InletVelocityX(X, A1, B1));
            C.InitialValues_Evaluators.Add("DisplacementY#C", X => InletVelocityY(X, A1, B1));
            C.InitialValues_Evaluators.Add("Pressure#C", X => 0);

            #endregion

            // boundary conditions
            // ===================
            #region BC
            C.AddBoundaryValue("wall_lower", "VelocityX", X => InletVelocityX(X, A3, B3));
            C.AddBoundaryValue("wall_lower", "VelocityY", X => InletVelocityY(X, A3, B3));
            C.AddBoundaryValue("wall_upper", "VelocityX", X => InletVelocityX(X, A3, B3));
            C.AddBoundaryValue("wall_upper", "VelocityY", X => InletVelocityY(X, A3, B3));
            C.AddBoundaryValue("wall_left", "VelocityX", X => InletVelocityX(X, A3, B3));
            C.AddBoundaryValue("wall_left", "VelocityY", X => InletVelocityY(X, A3, B3));
            C.AddBoundaryValue("wall_right", "VelocityX", X => InletVelocityX(X, A3, B3));
            C.AddBoundaryValue("wall_right", "VelocityY", X => InletVelocityY(X, A3, B3));
            C.AddBoundaryValue("wall_inside", "DisplacementX", X => InletVelocityX(X, A1, B1));
            C.AddBoundaryValue("wall_inside", "DisplacementY", X => InletVelocityY(X, A1, B1));
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
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;
            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = AMRlvl });
            C.AMR_startUpSweeps = AMRlvl;
            C.DynamicLoadBalancing_On = false;
            C.DynamicLoadBalancing_RedistributeAtStartup = false;
            C.GridPartType = GridPartType.METIS;
            #endregion

            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.TimesteppingMode = compMode;
            double dt = 1;
            C.dtMax = dt;
            C.dtMin = dt;
            C.NoOfTimesteps = 20;
            C.saveperiod = 1;

            #endregion

            return C;
        }

    }
}

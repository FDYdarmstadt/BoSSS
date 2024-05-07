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
    public static class Beam {

        public static ZLS_Control Channel(int p = 2, int kelem = 5, int AMRlvl = 0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.2;
            C.NoOfMultigridLevels = 1;
            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;



            // basic database options
            // ======================
            #region db
            string DbPath = @"C:\Databases\ZLS_GummiLippe";
            C.savetodb = false;
            //C.DbPath = DbPath;
            C.ProjectName = "ZLSTestVerticalBeam";
            //C.ProjectDescription = "Vertical Beam in channel";

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

            C.PhysicalParameters.rho_A = 0.1;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 0.05;
            C.PhysicalParameters.mu_B = 0.05;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid {
                Viscosity = 0.1,
                Lame2 = 1000,
                Density = 0.1
            };
            #endregion

            // grid generation
            // ===============
            #region grid

            double xLeft = -3;
            double xRight = 8;
            double ySize = 2;

            //*
            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(xLeft, xRight, ((int)(xRight - xLeft) * kelem) + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2 * kelem + 1);

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
            //*/


            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                if(X[0] < 0 && X[1] < 0.9) {
                    //y = 0.11
                    return -(X[0]).Pow2() + 0.1 * 0.1;
                }
                if(X[0] > 0 && X[1] < 0.9) {
                    return -(X[0]).Pow2() + 0.1 * 0.1;
                }
                return -((X[0]).Pow2() + (X[1] - 0.9).Pow2()) + 0.1 * 0.1;
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);



            double R(double t) {
                if(t <= 1) {
                    return (35 + (-84 + (70 - 20 * t) * t) * t) * t * t * t * t;
                } else {
                    return 1;
                }
            }

            double d = 0.2;


            double G(double y, double t) {
                if(y >= d)
                    return 1;
                else
                    return (1 - ((y - d) / d) * ((y - d) / d));
            }


            double vmax = 1;
            double inflowX(double[] x, double t) {
                return vmax * R(t) * G(x[1], t);
            }


            #endregion

            C.InitialValues_Evaluators.Add("VelocityX#A", x => 0);
            C.InitialValues_Evaluators.Add("VelocityX#B", x => 0);

            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("freeslip_upper");
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", inflowX);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", inflowX);
            //C.AddBoundaryValue("velocity_inlet_left", "VelocityY#A", inflowY);
            //dtC.AddBoundaryValue("velocity_inlet_left", "VelocityY#B", inflowY);
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
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = 2 });
            C.AMR_startUpSweeps = 2;

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
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 20;
            C.NoOfTimesteps = 2000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

    }
}

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

    public static class HardCodedControl {

        public static ZLS_Control Channel(int p = 3, int kelem = 5, int AMRlvl = 0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            // basic database options
            // ======================
            #region db
            C.savetodb = false;
            C.ContinueOnIoError = false;
            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 0.1;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid() {
                Lame2 = 2,
                Viscosity = 1,
                Density = 1
            };
            #endregion


            // grid generation
            // ===============
            #region grid


            //*
            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(-1, 1, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-1, 1, 2 * kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;

                    if(Math.Abs(X[1] + 1) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - 1) <= 1.0e-8)
                        et = 2;
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
                return -(0.07-X[1]);
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion

            // boundary conditions
            // ===================
            #region BC
            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            double vel = 0.01;
            C.AddBoundaryValue("wall_lower", "VelocityX#A", X => vel);
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
            C.NonLinearSolver.ConvergenceCriterion = 1e-10;
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

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            #endregion

            return C;
        }
        
        public static ZLS_Control DynamicChannel(int p = 2, int kelem = 5, int AMRlvl = 0) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            // basic database options
            // ======================
            #region db
            C.savetodb = false;
            C.ContinueOnIoError = false;
            #endregion

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 0.1;
            C.PhysicalParameters.rho_B = 0.1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 1;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid() {
                Lame2 = 2,
                Viscosity = 1,
                Density = 1
            };

            #endregion


            // grid generation
            // ===============
            #region grid


            //*
            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(-1, 1, 2 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-1, 1, 2 * kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "wall_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;

                    if(Math.Abs(X[1] + 1) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - 1) <= 1.0e-8)
                        et = 2;
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
                return -(0.07-X[1]);
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            #endregion

            // boundary conditions
            // ===================
            #region BC
            C.AddBoundaryValue("wall_lower");
            C.AddBoundaryValue("wall_upper");
            double vel = 0.5;
            C.AddBoundaryValue("wall_lower", "VelocityX#A", X => vel * (X[0] + 1).Pow2() * (X[0] - 1).Pow2());
            C.AddBoundaryValue("wall_lower", "VelocityX#B", X => vel);
            C.AddBoundaryValue("wall_lower", "VelocityX#C", X => vel);


            C.AddBoundaryValue("wall_lower", "DisplacementX#A", (X,t) => t * vel * (X[0] + 1).Pow2() * (X[0] - 1).Pow2());

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
            C.NonLinearSolver.ConvergenceCriterion = 1e-10;
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
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100;
            C.NoOfTimesteps = 1000;


            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            #endregion

            return C;
        }

        public static ZLS_Control AcceleratedBallInChannel(int p = 2, int kelem = 4) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;
            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            C.ProjectName = "ZLSTestVerticalBeam";

            C.ContinueOnIoError = false;

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

            C.Material = new AllOne();
            C.Material.Viscosity = 0;
            C.Material.Density = 2;
            C.Material.Lame2 = 100;

            #endregion


            // grid generation
            // ===============
            #region grid

            double xLeft = -3;
            double xRight = 3;
            double ySize = 2.1;

            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(xLeft, xRight, 6 * kelem + 1);
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

            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                return -((X[0]).Pow2() + (X[1] - 1).Pow2()) + 0.3.Pow2();
                //return -((X[0]).Pow2() + (X[1] - 1).Pow2()).Sqrt() + 0.3;
            };

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            double v0 = 0.00;
            C.InitialValues_Evaluators.Add("VelocityX#A", X => v0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => v0);
            C.InitialValues_Evaluators.Add("VelocityX#C", X => v0);
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

            double vmax = 1;
            double inflow(double[] x, double t) {
                    return v0 + R(t) * vmax;
            }

            double RInt(double t) {
                double integralR(double t) {
                    return -2.5 * t.Pow(8) + 10 * t.Pow(7) - 14 * t.Pow(6) + 7 * t.Pow(5);
                }
                if(t < 1) {
                    return integralR(t);
                } else {
                    return 1 * (t-1) + integralR(1);
                }
            }

            double integratedInflow(double[] x, double t) {
                return v0 * t + RInt(t) * vmax;
            }

            C.AddBoundaryValue("freeslip_lower");
            C.AddBoundaryValue("freeslip_upper");
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", inflow);
            C.AddBoundaryValue("velocity_inlet_left", "VelocityY#A", X => 0);
            C.AddBoundaryValue("pressure_outlet_right");

            C.AddBoundaryValue("velocity_inlet_left", "DisplacementX#A", (X,t) => integratedInflow(X, t));
            C.AddBoundaryValue("velocity_inlet_left", "DisplacementY#A", X => 0.0);

            #endregion

            // misc. solver options
            // ====================
            #region solver

            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.MinSolverIterations = 1;
            //C.Solver_MaxIterations = 50;
            
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
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

            C.TimeSteppingScheme = TimeSteppingScheme.BDF2;
            C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            C.TimesteppingMode = compMode;
            double dt = 1e-2;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 100;
            C.NoOfTimesteps = 300;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static ZLS_Control BallInChannel(int p = 2, int kelem = 4)
        {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.AgglomerationThreshold = 0.1;
            C.NoOfMultigridLevels = 1;

            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;
            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            C.ProjectName = "ZLSTestVerticalBeam";

            C.ContinueOnIoError = false;

            #endregion

            // Physical Parameters
            // ===================
            #region physics


            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 0.1;
            C.PhysicalParameters.mu_B = 0.1;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new Solid();
            C.Material.Viscosity = 0;
            C.Material.Density = 2;
            C.Material.Lame2 = 1000;

            #endregion


            // grid generation
            // ===============
            #region grid

            double xLeft = -3;
            double xRight = 3;
            double ySize = 2;

            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(xLeft, xRight, 6 * kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, ySize, 2 * kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "freeslip_lower");
                grd.EdgeTagNames.Add(2, "freeslip_upper");
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;

                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - ySize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] - xLeft) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[0] - xRight) <= 1.0e-8)
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

            

            Func<double[], double> Phi1Func = delegate (double[] X) {
                return -((X[0]).Pow2() + (X[1] - 1).Pow2()) + 0.3.Pow2();
                //return -((X[0]).Pow2() + (X[1] - 1).Pow2()).Sqrt() + 0.3;
            };
            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);

            double v0 = 1;
            C.InitialValues_Evaluators.Add("VelocityX#A", X => v0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => v0);
            C.InitialValues_Evaluators.Add("VelocityX#C", X => v0);
            #endregion


            // boundary conditions
            // ===================
            #region BC


            double inflow(double[] x, double t)
            {
                return v0;
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

            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.MinSolverIterations = 2;
            //C.Solver_MaxIterations = 50;

            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
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

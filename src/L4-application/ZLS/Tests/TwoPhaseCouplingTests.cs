using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
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

namespace ZwoLevelSetSolver.Tests
{
    class TwoPhaseCouplingTests
    {
        public static ZLS_Control BallInChannel(int p = 2, int kelem = 4)
        {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.AgglomerationThreshold = 0.1;
            C.NoOfMultigridLevels = 1;

            
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
            double R(double t)
            {
                if (t < 1)
                {
                    return (35 + (-84 + (70 - 20 * t) * t) * t) * t * t * t * t;
                }
                else
                {
                    return 1;
                }
            }

            double vmax = 0.01;
            double inflow(double[] x, double t)
            {
                return v0 + R(t) * vmax;
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
            C.LinearSolver.MaxSolverIterations = 20;
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

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 5e-3;
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

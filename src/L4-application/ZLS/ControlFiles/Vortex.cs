using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.ControlFiles {
    public static class Vortex {
        public static ZLS_Control SteadyVortex(int p = 2, int kelem = 16) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;

            int D = 2;

            // basic database options
            // ======================
            #region db

            C.savetodb = false;
            C.ProjectName = "ConvergenceTests";
            //C.ProjectDescription = "Test for Convergence";

            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            #endregion


            // DG degrees
            // ==========
            #region degrees
            //C.SetDGdegree(p);

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            C.Material = new ConvergenceTest();

            #endregion


            // grid generation
            // ===============
            #region grid

            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(-1, 1, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-1, 1, kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;

                    if (Math.Abs(X[1] + 1) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - 1) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + 1) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - 1) <= 1.0e-8)
                        et = 1;

                    return et;
                });

                return grd;

            };
            C.AddBoundaryValue("wall");

            #endregion


            // Initial Values
            // ==============
            #region init

            int K = 0;
            ForceX Vx = new ForceX();
            ForceY Vy = new ForceY();

            C.AddInitialValue("Phi", new Formula("X => -1"));
            C.AddInitialValue(VariableNames.SolidLevelSetCG, new Formula("X => 1"));
            C.AddInitialValue("GravityX#C", Vx);
            C.AddInitialValue("GravityY#C", Vy);

            #endregion

            // misc. solver options
            // ====================
            #region solver

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.ConvergenceCriterion = 0.0; // solve as accurate as possible

            C.LevelSet_ConvergenceCriterion = 1e-12;



            //C.Option_LevelSetEvolution = (compMode == AppControl._TimesteppingMode.Steady) ? LevelSetEvolution.None : LevelSetEvolution.FastMarching;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit = false; 

            C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;
            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

            //C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband { maxRefinementLevel = 3 });
            //C.AMR_startUpSweeps = 2;

            #endregion


            // Timestepping
            // ============
            #region time

            //C.CheckJumpConditions = true;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            #endregion

            return C;
        }

        public class ForceX : IBoundaryAndInitialData {

            double amplitude;
            public ForceX() {
                this.amplitude = 0.1;
            }

            public double Evaluate(double[] X, double t) {
                double r = Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
                double sum = -X[1] / r;
                return amplitude * sum * A(r);
            }

            public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
                NonVectorizedScalarFunction.Vectorize(this.Evaluate, time)(input, output);
            }
        };

        public class ForceY : IBoundaryAndInitialData {
            double amplitude;

            public ForceY() {
                this.amplitude = 0.1;
            }

            public double Evaluate(double[] X, double t) {
                double r = Math.Sqrt(X[0] * X[0] + X[1] * X[1]);
                double sum = X[0] / r;
                return amplitude * sum * A(r);

            }

            public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
                NonVectorizedScalarFunction.Vectorize(this.Evaluate, time)(input, output);
            }
        }

        static double A(double r) {
            return Math.Sin(Math.PI / 2 * r) / Math.Pow(1 + r, 6);
        }
    }
}

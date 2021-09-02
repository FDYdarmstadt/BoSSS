using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.SolidPhase;

namespace ZwoLevelSetSolver.ControlFiles {
    public static class ConvergenceTests {
        public static ZLS_Control RycroftPaper(int p = 2, int kelem = 18) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;
            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";



            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 0.001;
            C.PhysicalParameters.mu_B = 0.001;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new ConvergenceTest();

            #endregion


            // grid generation
            // ===============
            #region grid
            double xSize = 2;
            double ySize = 2;

            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(-(xSize / 2), (xSize / 2), kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-(ySize / 2), (ySize / 2), kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true, periodicY: true);

                return grd;

            };


            #endregion


            // basic database options
            // ======================
            #region db
            C.ProjectName = "ConvergenceTests";
            //C.ProjectDescription = "Test for Convergence";
            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                // only solid
                return 1;
            };

            Func<Vector, double, double> Vvorx = delegate (Vector X, double lamda) {
                return -Math.Sin(Math.PI * X[1]) * Math.Exp(-lamda * (2 - Math.Cos(Math.PI * X[0]) - Math.Cos(Math.PI * X[1])));
            };

            Func<Vector, double, double> Vvory = delegate (Vector X, double lamda) {
                return Math.Cos(Math.PI * X[0]) * Math.Exp(-lamda * (2 - Math.Cos(Math.PI * X[0]) - Math.Cos(Math.PI * X[1])));
            };

            int K = 0;

            Func<double[], double> Vx = delegate (double[] X) {
                double sum = 0;
                for (int k = 0; k <= K; k++) {
                    Vector x = new Vector(X[0] - (5.0 + 2.0 * k) / 6.0, X[1] - (5.0 + 2.0 * k) / 6.0);
                    sum = sum + (Math.Pow(-1, k) * Vvorx(x, 2 * (k + 1)));
                }
                return sum;
            };

            Func<double[], double> Vy = delegate (double[] X) {
                double sum = 0;
                for (int k = 0; k <= K; k++) {
                    Vector x = new Vector(X[0] - (5.0 + 2.0 * k) / 6.0, X[1] - (5.0 + 2.0 * k) / 6.0);
                    sum = sum + (Math.Pow(-1, k) * Vvory(x, 2 * (k + 1)));
                }
                return sum;
            };



            //C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);
            C.AddInitialValue(VariableNames.SolidLevelSetCG, new Formula("X => 1"));


            C.AddInitialValue("VelocityX#A", new VelocityX(K));
            C.AddInitialValue("VelocityX#B", new VelocityX(K));
            C.AddInitialValue("VelocityX#C", new VelocityX(K));

            C.InitialValues_Evaluators.Add("VelocityY#A", Vy);
            C.InitialValues_Evaluators.Add("VelocityY#B", Vy);
            C.InitialValues_Evaluators.Add("VelocityY#C", Vy);
            #endregion

            // boundary conditions
            // ===================
            #region BC

            #endregion

            // misc. solver options
            // ====================
            #region solver

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.MaxSolverIterations = 10;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-10;
            C.LinearSolver.ConvergenceCriterion = 1e-10;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

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
            double dt = 0.001;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 0.5;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static ZLS_Control SimpleTransport(int p = 2, int kelem = 18) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;
            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";



            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 0.001;
            C.PhysicalParameters.mu_B = 0.001;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new ConvergenceTest();

            #endregion


            // grid generation
            // ===============
            #region grid
            double xSize = 2;
            double ySize = 2;

            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(-(xSize / 2), (xSize / 2), kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-(ySize / 2), (ySize / 2), kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true, periodicY: true);

                return grd;

            };


            #endregion


            // basic database options
            // ======================
            #region db
            C.ProjectName = "ConvergenceTests";
            //C.ProjectDescription = "Test for Convergence";
            #endregion


            // Initial Values
            // ==============
            #region init

            //C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);
            C.AddInitialValue(VariableNames.SolidLevelSetCG, new Formula("X => 1"));


            C.InitialValues_Evaluators.Add("VelocityX#A", X => 0);
            C.InitialValues_Evaluators.Add("VelocityX#B", X => 0);
            C.InitialValues_Evaluators.Add("VelocityX#C", X => 0);

            C.InitialValues_Evaluators.Add("VelocityY#A", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY#B", X => 0);
            C.InitialValues_Evaluators.Add("VelocityY#C", X => 0.1);
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

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

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
            double dt = 0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 0.5;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion

            return C;
        }

        public static ZLS_Control SimpleConvergence(int p = 2, int kelem = 18) {
            ZLS_Control C = new ZLS_Control(p);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.AgglomerationThreshold = 0.3;
            C.NoOfMultigridLevels = 1;
            int D = 2;

            AppControl._TimesteppingMode compMode = AppControl._TimesteppingMode.Transient;

            //_DbPath = @"D:\local\local_Testcase_databases\Testcase_ContactLine";
            //_DbPath = @"D:\local\local_spatialConvStudy\StaticDropletOnPlateConvergence\SDoPConvDB";



            C.ContinueOnIoError = false;

            //C.LogValues = XNSE_Control.LoggingValues.MovingContactLine;
            //C.PostprocessingModules.Add(new MovingContactLineLogging());

            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 0.001;
            C.PhysicalParameters.mu_B = 0.001;

            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = false;

            C.Material = new ConvergenceTest();

            #endregion


            // grid generation
            // ===============
            #region grid
            double xSize = 2;
            double ySize = 2;

            C.GridFunc = delegate () {

                double[] Xnodes = GenericBlas.Linspace(-(xSize / 2), (xSize / 2), kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(-(ySize / 2), (ySize / 2), kelem + 1);

                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true, periodicY: true);

                return grd;

            };


            #endregion


            // basic database options
            // ======================
            #region db
            C.ProjectName = "ConvergenceTests";
            //C.ProjectDescription = "Test for Convergence";
            #endregion


            // Initial Values
            // ==============
            #region init

            Func<double[], double> PhiFunc = (X => -1);
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            Func<double[], double> Phi1Func = delegate (double[] X) {
                // only solid
                return 1;
            };


            Func<double[], double> Vx = delegate (double[] X) {
                double sum = 0;
                //sum = Math.Sin(Math.PI * X[1]);
                sum = 0.1;
                return sum;
            };

            Func<double[], double> Vy = delegate (double[] X) {
                double sum = 0;
                return sum;
            };



            //C.InitialValues_Evaluators.Add(VariableNames.SolidLevelSetCG, Phi1Func);
            C.AddInitialValue(VariableNames.SolidLevelSetCG, new Formula("X => 1"));


            C.InitialValues_Evaluators.Add("VelocityX#A", Vx);
            C.InitialValues_Evaluators.Add("VelocityX#B", Vx);
            C.InitialValues_Evaluators.Add("VelocityX#C", Vx);

            C.InitialValues_Evaluators.Add("VelocityY#A", Vy);
            C.InitialValues_Evaluators.Add("VelocityY#B", Vy);
            C.InitialValues_Evaluators.Add("VelocityY#C", Vy);
            #endregion

            // boundary conditions
            // ===================
            #region BC

            #endregion

            // misc. solver options
            // ====================
            #region solver

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.2;
            //C.AdvancedDiscretizationOptions.PenaltySafety = 40;
            //C.AdvancedDiscretizationOptions.UseGhostPenalties = true;

            C.NonLinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MinSolverIterations = 1;
            C.LinearSolver.MaxSolverIterations = 10;
            //C.Solver_MaxIterations = 50;
            C.NonLinearSolver.ConvergenceCriterion = 1e-10;
            C.LinearSolver.ConvergenceCriterion = 1e-10;
            //C.Solver_ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-12;

            C.AdvancedDiscretizationOptions.ViscosityMode = ViscosityMode.Standard;

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
            double dt = 0.01;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 0.5;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion

            return C;
        }
    } 

    class VelocityX : IBoundaryAndInitialData {
        int K;

        public VelocityX(int K) {
            this.K = K;
        }

        public double Evaluate(double[] X, double t) {
            double sum = 0;
            for (int k = 0; k <= K; k++) {
                Vector x = new Vector(X[0] - (5.0 + 2.0 * k) / 6.0, X[1] - (5.0 + 2.0 * k) / 6.0);
                sum = sum + (Math.Pow(-1, k) * Vvorx(x, 2 * (k + 1)));
            }
            return sum;
        }

        double Vvorx(double[] X, double lamda) {
            return -Math.Sin(Math.PI * X[1]) * Math.Exp(-lamda * (2 - Math.Cos(Math.PI * X[0]) - Math.Cos(Math.PI * X[1])));
        }

        public void Evaluate(MultidimensionalArray input, double time, MultidimensionalArray output) {
            NonVectorizedScalarFunction.Vectorize(this.Evaluate, time)(input, output);
        }
    }
}

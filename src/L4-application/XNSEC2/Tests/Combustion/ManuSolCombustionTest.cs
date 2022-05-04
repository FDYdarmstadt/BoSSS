using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;

namespace BoSSS.Application.XNSEC {

    public static partial class FullNSEControlExamples {

        public static XNSEC_Control NUnitTestManuSol_3() {
            //var C = ControlManuSolLowMachCombustion(DGp: 2, SizeFactor: 3);
            //C.DbPath = "c:\\BoSSS_DB";
            var C = ControlManuSolLowMachCombustion(DGp: 2, SizeFactor: 3);

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.ReactionRateConstants = new double[] { 1e4, 3, 1, 1 };
            C.HeatRelease = 1;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.None;

            C.VariableOneStepParameters = false;
            C.NumberOfChemicalSpecies = 4; // number of chemical species, without inert

            //C.myThermalWallType = SIPDiffusionTemperature.ThermalWallType.fixedTemperature;

            return C;
        }

        /// <summary>
        /// Manufactued solution described in the fdy annual report 2015 by Martin Karcher.
        /// </summary>
        public static XNSEC_Control ControlManuSolLowMachCombustion(int DGp = 3, int SizeFactor = 5) {
            XNSEC_Control C = new XNSEC_Control();
            // Solver configuration
            // ==============

            C.ManufacturedSolutionSwitch = true;
            C.savetodb = false;
            C.ProjectName = "Manufactured Solution 2D tests or NUnitTest";

            // Parameters
            // ==============
            C.AnalyticsolutionSwitch = true;
            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;
            C.NumberOfChemicalSpecies = 4;
            C.SetDGdegree(DGp);
            C.physicsMode = PhysicsMode.Combustion;
            C.Reynolds = 1.0;
            C.Froude = 1.0;
            C.Schmidt = 1.0;
            C.MatParamsMode = MaterialParamsMode.Constant;
            C.Prandtl = 1.0;
            //C.ReactionRateConstants = new double[] { 1e1, 1.0, 1.0, 1.0 };
            C.ReactionRateConstants = new double[] { 1e2, 2, 1, 1 };
            C.MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.LinearSolver = new Solution.AdvancedSolvers.OrthoMGSchwarzConfig() {
                NoOfMultigridLevels = 4
            };
            C.HeatRelease = 5;

            C.PhysicalParameters.IncludeConvection = true;

            // Grid declaration
            // ===============
            double h = Math.Pow(2, -SizeFactor + 1); // cell length
            double cells = 1 / h;
            int cells2 = (int)cells;
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                if (Math.Abs(x - 0) < 1e-8)
                    return 1;
                if (Math.Abs(y - 0) < 1e-8)
                    return 2;
                if (Math.Abs(x - 1) < 1e-8)
                    return 3;
                if (Math.Abs(y - 1) < 1e-8)
                    return 4;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 1, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(0, 1, cells2 + 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                grd.EdgeTagNames.Add(1, "Pressure_Outlet_left");
                grd.EdgeTagNames.Add(2, "Pressure_Outlet_bottom");
                grd.EdgeTagNames.Add(3, "Velocity_Inlet_right");
                grd.EdgeTagNames.Add(4, "Velocity_Inlet_top");
                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            double[] alphas = new double[] { 0.3, 0.4, 0.1, 0.2 };

            C.AnalyticVelocityX = X => -Math.Cos(X[0]);
            C.AnalyticVelocityY = X => -Math.Cos(X[1]);
            C.AnalyticPressure = X => Math.Sin(X[0] * X[1]);
            C.AnalyticTemperature = X => Math.Cos(X[0] * X[1]);
            C.AnalyticY0 = X => (0.3 * Math.Cos(X[0] * X[1]));
            C.AnalyticY1 = X => 0.4 * Math.Cos(X[0] * X[1]);
            C.AnalyticY2 = X => 0.1 * Math.Cos(X[0] * X[1]);
            C.AnalyticY3 = X => 0.2 * Math.Cos(X[0] * X[1]);
            C.AnalyticY4 = X => (1.0 - Math.Cos(X[0] * X[1]));

            double mult1 = 0.8; // should be bigger than 0
            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => C.AnalyticPressure(X) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => C.AnalyticVelocityX(X) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => C.AnalyticVelocityY(X) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => C.AnalyticTemperature(X) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => C.AnalyticY0(X) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction1, X => C.AnalyticY1(X) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction2, X => C.AnalyticY2(X) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction3, X => C.AnalyticY3(X) * mult1);

            // Analytical / manufactured solutions
            C.InitialValues_Evaluators.Add(VariableNames.Pressure + "_an", X => C.AnalyticPressure(X));
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0) + "_an", X => C.AnalyticVelocityX(X));
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1) + "_an", X => C.AnalyticVelocityY(X));
            C.InitialValues_Evaluators.Add(VariableNames.Temperature + "_an", X => C.AnalyticTemperature(X));
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "_an", X => C.AnalyticY0(X));
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction1 + "_an", X => C.AnalyticY1(X));
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction2 + "_an", X => C.AnalyticY2(X));
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction3 + "_an", X => C.AnalyticY3(X));

            // boundary conditions
            // ===================
            C.AddBoundaryValue("Velocity_Inlet_right", VariableNames.Velocity_d(0), C.AnalyticVelocityX);
            C.AddBoundaryValue("Velocity_Inlet_right", VariableNames.Velocity_d(1), C.AnalyticVelocityY);
            C.AddBoundaryValue("Velocity_Inlet_right", VariableNames.Temperature, C.AnalyticTemperature);
            C.AddBoundaryValue("Velocity_Inlet_right", VariableNames.MassFraction0, C.AnalyticY0);
            C.AddBoundaryValue("Velocity_Inlet_right", VariableNames.MassFraction1, C.AnalyticY1);
            C.AddBoundaryValue("Velocity_Inlet_right", VariableNames.MassFraction2, C.AnalyticY2);
            C.AddBoundaryValue("Velocity_Inlet_right", VariableNames.MassFraction3, C.AnalyticY3);
            C.AddBoundaryValue("Velocity_Inlet_right", VariableNames.MassFraction4, C.AnalyticY4);

            C.AddBoundaryValue("Velocity_Inlet_top", VariableNames.Velocity_d(0), C.AnalyticVelocityX);
            C.AddBoundaryValue("Velocity_Inlet_top", VariableNames.Velocity_d(1), C.AnalyticVelocityY);
            C.AddBoundaryValue("Velocity_Inlet_top", VariableNames.Temperature, C.AnalyticTemperature);
            C.AddBoundaryValue("Velocity_Inlet_top", VariableNames.MassFraction0, C.AnalyticY0);
            C.AddBoundaryValue("Velocity_Inlet_top", VariableNames.MassFraction1, C.AnalyticY1);
            C.AddBoundaryValue("Velocity_Inlet_top", VariableNames.MassFraction2, C.AnalyticY2);
            C.AddBoundaryValue("Velocity_Inlet_top", VariableNames.MassFraction3, C.AnalyticY3);
            C.AddBoundaryValue("Velocity_Inlet_top", VariableNames.MassFraction4, C.AnalyticY4);

            C.AddBoundaryValue("Pressure_Outlet_left", VariableNames.Pressure, C.AnalyticPressure);
            C.AddBoundaryValue("Pressure_Outlet_bottom", VariableNames.Pressure, C.AnalyticPressure);

            //Manufactured Solutions
            C.ManufacturedSolution_MomentumX = delegate (double[] X, double t) {
                double[] x = X;
                double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

                double p0 = 1.0;
                double x_ = x[0];
                double y_ = x[1];

                double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
                double alpha1 = 0.3;
                double alpha2 = 0.4;
                double alpha3 = 0.1;
                double alpha4 = 0.2;

                double res = 0;

                double ConvectionTerm = p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Pow(Math.Cos(x_), 0.2e1) * y_ * Math.Sin(x_ * y_) - p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Pow(Math.Cos(x_), 0.2e1) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) - 0.2e1 * p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(x_) * Math.Sin(x_) + p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(x_) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) - p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) - p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(x_) * Math.Sin(y_);

                double ViscTerm = -0.4e1 / 0.3e1 * Math.Cos(x_) / C.Reynolds; // OK
                double PressureGradientTerm = y_ * Math.Cos(x_ * y_); // OK
                double BouyancyTerm = 0.0;
                res = ConvectionTerm + ViscTerm + PressureGradientTerm - BouyancyTerm;

                return -res;
            };

            C.ManufacturedSolution_MomentumY = delegate (double[] X, double t) {
                double[] x = X;
                double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

                double p0 = 1.0;
                double x_ = x[0];
                double y_ = x[1];

                double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
                double alpha1 = 0.3;
                double alpha2 = 0.4;
                double alpha3 = 0.1;
                double alpha4 = 0.2;

                double res = 0;

                double ConvectionTerm = p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(x_) * Math.Cos(y_) * y_ * Math.Sin(x_ * y_) - p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * Math.Cos(y_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) - p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_) * Math.Cos(y_) + p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Pow(Math.Cos(y_), 0.2e1) * x_ * Math.Sin(x_ * y_) - p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Pow(Math.Cos(y_), 0.2e1) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) - 0.2e1 * p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(y_) * Math.Sin(y_);
                double BouyancyTerm = -Math.Pow(C.Froude, -0.2e1) * p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5);

                double ViscTerm = -0.4e1 / 0.3e1 * Math.Cos(y_) / C.Reynolds; //
                double PressureGradientTerm = x_ * Math.Cos(x_ * y_);
                res = ConvectionTerm + ViscTerm + PressureGradientTerm - BouyancyTerm;

                return -res;
            };

            C.ManufacturedSolution_Continuity = delegate (double[] X, double t) {
                double x_ = X[0];
                double y_ = X[1];
                double p0 = 1.0;
                double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

                double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
                double alpha1 = 0.3;
                double alpha2 = 0.4;
                double alpha3 = 0.1;
                double alpha4 = 0.2;

                double dRhoUdx = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(x_) * y_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_);
                double dRhoVdy = -p0 * Math.Pow(Math.Cos(x_ * y_), -0.2e1) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Cos(y_) * x_ * Math.Sin(x_ * y_) + p0 / Math.Cos(x_ * y_) * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / Math.Cos(x_ * y_) / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                double man1 = -1 * (dRhoUdx + dRhoVdy);
                return man1;
            };

            C.ManufacturedSolution_Energy = delegate (double[] X, double t) {
                double[] x = X;
                double x_ = x[0];
                double y_ = x[1];
                double p0 = 1.0;
                double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };
                double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
                double alpha1 = 0.3;
                double alpha2 = 0.4;
                double alpha3 = 0.1;
                double alpha4 = 0.2;
                double[] Coefficients = new double[] { alpha1, alpha2, alpha3, alpha4 };

                double ReactionRate = 0;
                double ConvectionTermSwitch = 1.0;
                double DiffussionTermSwitch = 1.0;
                double Da = C.ReactionRateConstants[0];
                double Ta = C.ReactionRateConstants[1];
                double a = C.ReactionRateConstants[2];
                double b = C.ReactionRateConstants[3];
                double res = 0.0;
                double alpha;

                double ConvectionTerm = p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_) + p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha1 / M1, a) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha2 / M2, b);

                if (!C.ChemicalReactionActive)
                    ReactionRate = 0.0;
                double DiffussionTerm = 1.0 / (C.Reynolds * C.Prandtl) * Math.Cos(x_ * y_) * (Math.Pow(x_, 2) + Math.Pow(y_, 2)); // OK for lowmach AND combustion
                double SourceTerm = -C.HeatRelease * ReactionRate * M1;
                res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                return -res;
            };

            C.ManufacturedSolution_Species0 = delegate (double[] X, double t) {
                int SpeciesIndex = 0;
                double[] x = X;
                double x_ = x[0];
                double y_ = x[1];
                double p0 = 1.0;
                double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

                double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
                double alpha1 = 0.3;
                double alpha2 = 0.4;
                double alpha3 = 0.1;
                double alpha4 = 0.2;
                double[] Coefficients = new double[] { alpha1, alpha2, alpha3, alpha4 };
                double Da = C.ReactionRateConstants[0];
                double Ta = C.ReactionRateConstants[1];
                double a = C.ReactionRateConstants[2];
                double b = C.ReactionRateConstants[3];

                double ConvectionTerm = 0;
                double ReactionRate = 0;
                double SourceTerm = 0;
                double DiffussionTerm = 0;

                double ConvectionTermSwitch = 1.0;
                double DiffussionTermSwitch = 1.0;
                double res = 0.0;
                double[] alphas = new double[] { alpha1, alpha2, alpha3, alpha4 };
                double alpha = alphas[SpeciesIndex];
                double StoichiometricCoeff = C.StoichiometricCoefficients[SpeciesIndex];
                ConvectionTerm = p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_) + p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha1 / M1, a) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha2 / M2, b);

                if (!C.ChemicalReactionActive)
                    ReactionRate = 0.0;
                ConvectionTerm *= Coefficients[SpeciesIndex];

                DiffussionTerm = 1.0 / (C.Reynolds * C.Schmidt) * alpha * y_ * y_ * Math.Cos(x_ * y_) + alpha * x_ * x_ * Math.Cos(x_ * y_);

                SourceTerm = -MolarMasses[SpeciesIndex] * StoichiometricCoeff * ReactionRate;
                res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                return -(res);
            };

            C.ManufacturedSolution_Species1 = delegate (double[] X, double t) {
                int SpeciesIndex = 1;
                double[] x = X;
                double x_ = x[0];
                double y_ = x[1];
                double p0 = 1.0;
                double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

                double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
                double alpha1 = 0.3;
                double alpha2 = 0.4;
                double alpha3 = 0.1;
                double alpha4 = 0.2;
                double[] Coefficients = new double[] { alpha1, alpha2, alpha3, alpha4 };
                double Da = C.ReactionRateConstants[0];
                double Ta = C.ReactionRateConstants[1];
                double a = C.ReactionRateConstants[2];
                double b = C.ReactionRateConstants[3];

                double ConvectionTerm = 0;
                double ReactionRate = 0;
                double SourceTerm = 0;
                double DiffussionTerm = 0;

                double ConvectionTermSwitch = 1.0;
                double DiffussionTermSwitch = 1.0;
                double res = 0.0;
                double[] alphas = new double[] { alpha1, alpha2, alpha3, alpha4 };
                double alpha = alphas[SpeciesIndex];
                double StoichiometricCoeff = C.StoichiometricCoefficients[SpeciesIndex];
                ConvectionTerm = p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_) + p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha1 / M1, a) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha2 / M2, b);

                if (!C.ChemicalReactionActive)
                    ReactionRate = 0.0;

                ConvectionTerm *= Coefficients[SpeciesIndex];

                DiffussionTerm = 1.0 / (C.Reynolds * C.Schmidt) * alpha * y_ * y_ * Math.Cos(x_ * y_) + alpha * x_ * x_ * Math.Cos(x_ * y_);

                SourceTerm = -MolarMasses[SpeciesIndex] * StoichiometricCoeff * ReactionRate;
                res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                return -(res);
            };

            C.ManufacturedSolution_Species2 = delegate (double[] X, double t) {
                int SpeciesIndex = 2;
                double[] x = X;
                double x_ = x[0];
                double y_ = x[1];
                double p0 = 1.0;
                double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

                double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
                double alpha1 = 0.3;
                double alpha2 = 0.4;
                double alpha3 = 0.1;
                double alpha4 = 0.2;
                double[] Coefficients = new double[] { alpha1, alpha2, alpha3, alpha4 };
                double Da = C.ReactionRateConstants[0];
                double Ta = C.ReactionRateConstants[1];
                double a = C.ReactionRateConstants[2];
                double b = C.ReactionRateConstants[3];

                double ConvectionTerm = 0;
                double ReactionRate = 0;
                double SourceTerm = 0;
                double DiffussionTerm = 0;

                double ConvectionTermSwitch = 1.0;
                double DiffussionTermSwitch = 1.0;
                double res = 0.0;
                double[] alphas = new double[] { alpha1, alpha2, alpha3, alpha4 };
                double alpha = alphas[SpeciesIndex];
                double StoichiometricCoeff = C.StoichiometricCoefficients[SpeciesIndex];
                ConvectionTerm = p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_) + p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha1 / M1, a) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha2 / M2, b);

                if (!C.ChemicalReactionActive)
                    ReactionRate = 0.0;

                ConvectionTerm *= Coefficients[SpeciesIndex];

                DiffussionTerm = 1.0 / (C.Reynolds * C.Schmidt) * alpha * y_ * y_ * Math.Cos(x_ * y_) + alpha * x_ * x_ * Math.Cos(x_ * y_);

                SourceTerm = -MolarMasses[SpeciesIndex] * StoichiometricCoeff * ReactionRate;
                res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                return -(res);
            };

            C.ManufacturedSolution_Species3 = delegate (double[] X, double t) {
                int SpeciesIndex = 3;
                double[] x = X;
                double x_ = x[0];
                double y_ = x[1];
                double p0 = 1.0;
                double[] MolarMasses = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0 };

                double M1 = MolarMasses[0]; double M2 = MolarMasses[1]; double M3 = MolarMasses[2]; double M4 = MolarMasses[3]; double M5 = MolarMasses[4];
                double alpha1 = 0.3;
                double alpha2 = 0.4;
                double alpha3 = 0.1;
                double alpha4 = 0.2;
                double[] Coefficients = new double[] { alpha1, alpha2, alpha3, alpha4 };
                double Da = C.ReactionRateConstants[0];
                double Ta = C.ReactionRateConstants[1];
                double a = C.ReactionRateConstants[2];
                double b = C.ReactionRateConstants[3];

                double ConvectionTerm = 0;
                double ReactionRate = 0;
                double SourceTerm = 0;
                double DiffussionTerm = 0;

                double ConvectionTermSwitch = 1.0;
                double DiffussionTermSwitch = 1.0;
                double res = 0.0;
                double[] alphas = new double[] { alpha1, alpha2, alpha3, alpha4 };
                double alpha = alphas[SpeciesIndex];
                double StoichiometricCoeff = C.StoichiometricCoefficients[SpeciesIndex];
                ConvectionTerm = p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(x_) * (-alpha1 * y_ * Math.Sin(x_ * y_) / M1 - alpha2 * y_ * Math.Sin(x_ * y_) / M2 - alpha3 * y_ * Math.Sin(x_ * y_) / M3 - alpha4 * y_ * Math.Sin(x_ * y_) / M4 + (alpha1 * y_ * Math.Sin(x_ * y_) + alpha2 * y_ * Math.Sin(x_ * y_) + alpha3 * y_ * Math.Sin(x_ * y_) + alpha4 * y_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(x_) + p0 * Math.Pow(alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5, -0.2e1) * Math.Cos(y_) * (-alpha1 * x_ * Math.Sin(x_ * y_) / M1 - alpha2 * x_ * Math.Sin(x_ * y_) / M2 - alpha3 * x_ * Math.Sin(x_ * y_) / M3 - alpha4 * x_ * Math.Sin(x_ * y_) / M4 + (alpha1 * x_ * Math.Sin(x_ * y_) + alpha2 * x_ * Math.Sin(x_ * y_) + alpha3 * x_ * Math.Sin(x_ * y_) + alpha4 * x_ * Math.Sin(x_ * y_)) / M5) + p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * Math.Sin(y_);

                ReactionRate = Da * Math.Exp(-Ta / Math.Cos(x_ * y_)) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha1 / M1, a) * Math.Pow(p0 / (alpha1 * Math.Cos(x_ * y_) / M1 + alpha2 * Math.Cos(x_ * y_) / M2 + alpha3 * Math.Cos(x_ * y_) / M3 + alpha4 * Math.Cos(x_ * y_) / M4 + (0.10e1 - alpha1 * Math.Cos(x_ * y_) - alpha2 * Math.Cos(x_ * y_) - alpha3 * Math.Cos(x_ * y_) - alpha4 * Math.Cos(x_ * y_)) / M5) * alpha2 / M2, b);

                if (!C.ChemicalReactionActive)
                    ReactionRate = 0.0;

                ConvectionTerm *= Coefficients[SpeciesIndex];

                DiffussionTerm = 1.0 / (C.Reynolds * C.Schmidt) * alpha * y_ * y_ * Math.Cos(x_ * y_) + alpha * x_ * x_ * Math.Cos(x_ * y_);

                SourceTerm = -MolarMasses[SpeciesIndex] * StoichiometricCoeff * ReactionRate;
                res = ConvectionTerm * ConvectionTermSwitch + DiffussionTerm * DiffussionTermSwitch + SourceTerm;
                return -(res);
            };

            return C;
        }
    }
}
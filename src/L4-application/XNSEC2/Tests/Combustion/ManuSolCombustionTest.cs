using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;

namespace BoSSS.Application.XNSEC {

    static public partial class FullNSEControlExamples {

        static public XNSEC_Control NUnitTestManuSol_3() {
            //var C = ControlManuSolLowMachCombustion(DGp: 1, SizeFactor: 2);
            //C.DbPath = "c:\\BoSSS_DB";
            var C = ControlManuSolLowMachCombustion(DGp: 2, SizeFactor:4);
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.ReactionRateConstants = new double[] { 1e4, 3, 1, 1 };
            C.HeatRelease = 1;
            C.ImmediatePlotPeriod = 1;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;



            C.VariableOneStepParameters = false;
            C.NumberOfChemicalSpecies = 4; // number of chemical species, without inert
            //C.ImmediatePlotPeriod = 1;

            //C.myThermalWallType = SIPDiffusionTemperature.ThermalWallType.fixedTemperature;


            return C;
        }

        /// <summary>
        /// Manufactued solution described in the fdy annual report 2015 by Martin Karcher.
        /// </summary>
        static public XNSEC_Control ControlManuSolLowMachCombustion(int DGp = 3, int SizeFactor = 5) {
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
            C.LinearSolver.SolverCode = LinearSolverCode.exp_Kcycle_schwarz;
            C.LinearSolver.NoOfMultigridLevels = 4;
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

                if(Math.Abs(x - 0) < 1e-8)
                    return 1;
                if(Math.Abs(y - 0) < 1e-8)
                    return 2;
                if(Math.Abs(x - 1) < 1e-8)
                    return 3;
                if(Math.Abs(y - 1) < 1e-8)
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
            //C.InitialValues_Evaluators.Add("p_an", C.AnalyticPressure);
            //C.InitialValues_Evaluators.Add("u_an", C.AnalyticVelocityX);
            //C.InitialValues_Evaluators.Add("v_an", C.AnalyticVelocityY);
            //C.InitialValues_Evaluators.Add("T_an", C.AnalyticTemperature);
            //C.InitialValues_Evaluators.Add("Y0_an", C.AnalyticY0);
            //C.InitialValues_Evaluators.Add("Y1_an", C.AnalyticY1);
            //C.InitialValues_Evaluators.Add("Y2_an", C.AnalyticY2);
            //C.InitialValues_Evaluators.Add("Y3_an", C.AnalyticY3);

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
            return C;
        }
    }
}
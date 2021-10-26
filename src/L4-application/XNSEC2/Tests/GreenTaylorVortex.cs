using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;
using static System.Math;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    /// </summary>
    static public partial class FullNSEControlExamples {

        /// <summary>
        /// NUnit test
        /// </summary>
        static public XNSEC_Control NUnitUnsteadyTaylorVortex() {
            var C = UnsteadyTaylorVortex(3, 7);
            C.rhoOne = false;
            C.timeDerivativeConti_OK = false;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-7;
            //C.NonLinearSolver.MaxSolverIterations = 20;
            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.DbPath = null;// @"C:\Databases\BoSSS_DB";
            C.savetodb = false;
            C.ChemicalReactionActive = false;
            C.ImmediatePlotPeriod = 1;
            C.UseSelfMadeTemporalOperator = false;
            C.timeDerivativeConti_OK = false;
            C.timeDerivativeEnergyp0_OK = false;
            C.PhysicalParameters.IncludeConvection = true;
            //C.NonLinearSolver.Globalization = Solution.AdvancedSolvers.Newton.GlobalizationOption.LineSearch;
            return C;
        }

        // Valid for incompressible
        static public XNSEC_Control UnsteadyTaylorVortex(int DGp = 3, int ncells = 10) {
            XNSEC_Control C = new XNSEC_Control();

            C.savetodb = false;
            C.ProjectName = "Unsteady_LowMach Taylor vortex";
            C.ProjectDescription = "NUnitTest for UnsteadyLowMach";

            // Solver configuration
            // ==============

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.physicsMode = PhysicsMode.Combustion;
            C.Endtime = 0.5;
            C.dtFixed = 0.05;
            C.NumberOfChemicalSpecies = 1;
            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

            C.PhysicalParameters.IncludeConvection = false;

            C.EnableTemperature = false; // T=1
            C.EnableMassFractions = false; // Y_0 = 1
            //C.MomentumConvection_OK = false; // Only Stokes
            C.rhoOne = true;

            C.NoOfTimesteps = int.MaxValue;
            C.MatParamsMode = MaterialParamsMode.Constant;
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.Constant;
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient; // Unsteady simulation...

            // Parameters
            // ==============
            C.AnalyticsolutionSwitch = true;
            C.Reynolds = 100.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;

            //Grid declaration
            // ===============

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];
                //Left
                if(Math.Abs(x + 1) < 1e-8)
                    return 1;

                // bottom
                if(Math.Abs(y + 1) < 1e-8)
                    return 1;

                //right
                if(Math.Abs(x - 1) < 1e-8)
                    return 1;

                //top
                if(Math.Abs(y - 1) < 1e-8)
                    return 1;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1.0, 1.0, ncells + 1);
                var _yNodes = GenericBlas.Linspace(-1.0, 1.0, ncells + 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, periodicX: true, periodicY: true);
                //  grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // Analytic solutions
            // ===================
            double time = 0.5;
            double reynolds = C.Reynolds;
            Func<double[], double> _AnalyticVelocityX = X => -Cos(PI * X[0]) * Sin(PI * X[1]) * Exp(-2.0 * PI * PI * time / reynolds);
            Func<double[], double> _AnalyticVelocityY = X => Sin(PI * X[0]) * Cos(PI * X[1]) * Exp(-2.0 * PI * PI * time / reynolds);
            Func<double[], double> _AnalyticPressure = X => -0.25 * (Cos(2.0 * PI * X[0]) + Cos(2.0 * PI * X[1])) * Exp(-4.0 * PI * PI * time / reynolds);
            Func<double[], double> _AnalyticTemperature = X => 1.0;
            Func<double[], double> _AnalyticY0 = X => 1.0;

            C.AnalyticVelocityX = _AnalyticVelocityX;
            C.AnalyticVelocityY = _AnalyticVelocityY;
            C.AnalyticPressure = _AnalyticPressure;
            C.AnalyticTemperature = _AnalyticTemperature;
            C.AnalyticY0 = _AnalyticY0;

            // Analytical / manufactured solutions
            // ==============
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0) + "_an", _AnalyticVelocityX);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1) + "_an", _AnalyticVelocityY);
            C.InitialValues_Evaluators.Add(VariableNames.Pressure + "_an", _AnalyticPressure);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature + "_an", _AnalyticTemperature);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "_an", _AnalyticY0);

            // initial values
            // ==============
            C.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => -Cos(PI * X[0]) * Sin(PI * X[1]));
            C.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => Sin(PI * X[0]) * Cos(PI * X[1]));
            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => -0.25 * (Cos(2.0 * PI * X[0]) + Cos(2.0 * PI * X[1])));
            C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);
            return C;
        }
    }
}
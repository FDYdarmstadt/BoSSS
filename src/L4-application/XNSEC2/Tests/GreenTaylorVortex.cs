using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Queries;
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
            var C = UnsteadyTaylorVortex(3, 10,0.1,3,null,1.0);
            //var C = UnsteadyTaylorVortex(1, 5);

            C.rhoOne = false;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.DbPath = null;// @"C:\Databases\BoSSS_DB";
            C.savetodb = false;
         
            //C.NonLinearSolver.Globalization = Solution.AdvancedSolvers.Newton.GlobalizationOption.LineSearch;
            return C;
        }

        // Valid for incompressible
        static public XNSEC_Control UnsteadyTaylorVortex(int DGp = 3, int ncells = 10, double dt = 1e-1, int bdfOrder = 1, string dbPath = null, double finaltime = 1.0) {
            XNSEC_Control C = new XNSEC_Control();
            C.savetodb = false;
            
            if (dbPath != null) {
                C.DbPath = dbPath;
                C.savetodb = true;
            }
            
            C.ProjectName = "Unsteady_LowMach Taylor vortex";
            C.ProjectDescription = "NUnitTest for UnsteadyLowMach";

            // Solver configuration
            // ==============

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.physicsMode = PhysicsMode.Combustion;
            C.Endtime = finaltime;
            C.dtFixed = dt;
            C.NumberOfChemicalSpecies = 1;
            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };
            C.PhysicalParameters.IncludeConvection = true;

            switch (bdfOrder) {
                case 1:
                    C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
                    break;
                case 2:
                    C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.BDF2;
                    break;
                case 3:
                    C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.BDF3;
                    break;
                case 4:
                    C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.BDF4;
                    break;
                default:
                    C.TimeSteppingScheme = Solution.XdgTimestepping.TimeSteppingScheme.BDF3;
                    break;

                    
            }

            C.NoOfTimesteps = int.MaxValue;
            C.MatParamsMode = MaterialParamsMode.Constant;
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.Constant;
            C.TimesteppingMode = AppControl._TimesteppingMode.Transient; // Unsteady simulation...

            C.ChemicalReactionActive = false;
            C.UseSelfMadeTemporalOperator = false;
            C.timeDerivativeConti_OK = false;
            C.timeDerivativeEnergyp0_OK = false;
            C.PhysicalParameters.IncludeConvection = true;
            C.EnableTemperature = false;
            C.EnableMassFractions = false;

            // Parameters
            // ==============
            C.AnalyticsolutionSwitch = true;
            C.Reynolds = 10.0;
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
            double time = C.Endtime;
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

            int queryQuadOrder = 20;
            C.Queries.Add(
                "SolL2err_u",
                QueryLibrary.L2Error(VariableNames.VelocityX, C.AnalyticVelocityX, queryQuadOrder));
            C.Queries.Add(
                "SolL2err_v",
                QueryLibrary.L2Error(VariableNames.VelocityY, C.AnalyticVelocityY, queryQuadOrder));
            C.Queries.Add(
                "SolL2err_p",
                QueryLibrary.L2ErrorNoMean(VariableNames.Pressure, C.AnalyticPressure, queryQuadOrder));
            C.Queries.Add(
                "SolL2err_T",
                QueryLibrary.L2Error(VariableNames.Temperature, C.AnalyticTemperature, queryQuadOrder));


            return C;
        }
    }
}
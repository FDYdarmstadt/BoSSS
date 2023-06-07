using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;

namespace BoSSS.Application.XNSEC {

    static public partial class FullNSEControlExamples {

        static public XNSEC_Control NonReactingMixingLayerTests(int DGp = 1, int meshScaling = 8, double dt = 1.0, string dbPath = @"C:\Databases\BoSSS_DB") {
            XNSEC_Control C = new XNSEC_Control();

            // Solver configuration
            // ==============

            C.physicsMode = PhysicsMode.Combustion;
            C.MatParamsMode = MaterialParamsMode.Sutherland;
            C.rhoOne = false;
            C.ChemicalReactionActive = false;
            C.GravityDirection[1] = 0.0; // no gravity

            C.saveperiod = 1;
            C.SetSaveOptions(dbPath, 1);

            //if(dbPath != null) {
            //    C.DbPath = dbPath;
            //    C.savetodb = true;
            //} else {
            //    C.DbPath = null;
            //    C.savetodb = false;
            //}

            C.ProjectName = "Steady Mixing Layer Test";

            C.NumberOfChemicalSpecies = 2;
            C.SetDGdegree(DGp);

            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.verbose = true;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            if (dt < 0) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                C.dtFixed = dt;
                C.Endtime = 400.0;
                C.NoOfTimesteps = int.MaxValue;
                C.UseSelfMadeTemporalOperator = true;
            }

            // Grid declaration
            // ===============
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet O_2
                //2: Wall upper
                //3: Velocity inlet CH_4
                //4: Pressure outlet
                //5: Wall lower

                //Inlet oxidizer
                if (Math.Abs(x - 0) < 1e-8)
                    return 1;

                //upper Wall
                if (Math.Abs(y - 0.1) < 1e-8 && (Math.Abs(x - 0.1) < 0.1 + 1e-8 || Math.Abs(x - 0.65) < 0.35 + 1e-8))
                    return 2;

                //Inlet fuel
                if (Math.Abs(y - 0.1) < 1e-8 && Math.Abs(x - 0.25) < 0.05 + 1e-8)
                    return 3;

                //Outlet
                if (Math.Abs(x - 1) < 1e-8)
                    return 4;

                //lower Wall
                if (Math.Abs(y - (-0.1)) < 1e-8)
                    return 2;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 1, 2 * meshScaling + 1);
                var _yNodes = GenericBlas.Linspace(-0.1, 0.1, meshScaling + 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                grd.EdgeTagNames.Add(1, "Velocity_Inlet_O2");
                //grd.EdgeTagNames.Add(2, "Wall_upper");
                grd.EdgeTagNames.Add(2, "NoSlipNeumann");
                grd.EdgeTagNames.Add(3, "Velocity_Inlet_CH4");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");
                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // initial values
            // ==============

            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction1, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction2, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction3, X => 0.0);

            // boundary conditions
            // ===================
            // Species 0: Fuel
            // Species 1: Oxidizer
            // Species 2: Product1
            // Species 3: Product2

            double airProportion = 1.0;
            double velMult = 0.3;

            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Velocity_d(0), (X, t) => 1.0 - 100 * Math.Pow(X[1], 2));
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Velocity_d(1), (X, t) => 0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Temperature, (X, t) => 1.0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction0, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction1, (X, t) => 1.0 * airProportion);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction2, (X, t) => 1.0 * (1.0 - airProportion));
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction3, (X, t) => 0.0);

            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Velocity_d(0), (X, t) => 0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Velocity_d(1), (X, t) => -1.0 * velMult * (1 - 1.0 / Math.Pow(0.05, 2) * Math.Pow(X[0] - 0.25, 2)));
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Temperature, (X, t) => 1.5);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction0, (X, t) => 1.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction1, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction2, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction3, (X, t) => 0.0);

            //C.AddBoundaryValue("Wall_upper", VariableNames.Temperature, (X, t) => 1.0);
            //C.AddBoundaryValue("Wall_lower", VariableNames.Temperature, (X, t) => 1.0);
            C.AddBoundaryValue("NoSlipNeumann");

            C.AddBoundaryValue("Pressure_Outlet");

            // Parameters
            // ==============

            C.Reynolds = 1.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;

            C.T_ref_Sutherland = 600;

            C.ReactionRateConstants = new double[] { 1, 3, 1.0, 2.0 };
            C.StoichiometricCoefficients = new double[] { -1, -2, 1, 2 };
            C.M_ref = 16;
            C.PartialHeatCapacities = new double[] { 1.0, 1.0, 1.0, 1.0 };
            return C;
        }

        static public XNSEC_MF_Control MixingLayerTests_MixtureFraction(int DGp = 1, int meshScaling = 10, double dt = -1.0, string dbPath = @"C:\Databases\BoSSS_DB") {
            XNSEC_MF_Control C = new XNSEC_MF_Control();

            // Solver configuration
            // ==============

            C.physicsMode = PhysicsMode.MixtureFraction;
            C.MatParamsMode = MaterialParamsMode.Constant;
            C.rhoOne = true;
            C.ChemicalReactionActive = false;
            C.GravityDirection[1] = 0.0; // no gravity

            C.saveperiod = 1;
            C.SetSaveOptions(dbPath, 1);
            C.ImmediatePlotPeriod = 1;
            //if(dbPath != null) {
            //    C.DbPath = dbPath;
            //    C.savetodb = true;
            //} else {
            //    C.DbPath = null;
            //    C.savetodb = false;
            //}

            C.ProjectName = "Steady Mixing Layer Test";

            C.NumberOfChemicalSpecies = 2;
            C.SetDGdegree(DGp);
            C.HeatRelease = 10;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.verbose = true;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            //C.NonLinearSolver.MaxSolverIterations = 1; C.NonLinearSolver.MinSolverIterations= 1;
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            if (dt < 0) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                C.dtFixed = dt;
                C.Endtime = 400.0;
                C.NoOfTimesteps = int.MaxValue;
                C.UseSelfMadeTemporalOperator = true;
                C.timeDerivativeConti_OK = true;
                //C.PlotAdditionalParameters = true;
            }

            // Grid declaration
            // ===============
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet O_2
                //2: Wall upper
                //3: Velocity inlet CH_4
                //4: Pressure outlet
                //5: Wall lower

                //Inlet oxidizer
                if (Math.Abs(x - 0) <=1e-8 && (y - 0.0 >= 1e-8) )
                    return 1;

                //upper Wall
                if (Math.Abs(y - 0.1) <= 1e-8)
                    return 2;
                //lower Wall
                if (Math.Abs(y - (-0.1)) <= 1e-8)
                    return 2;

                //Inlet fuel
                if (Math.Abs(x - 0) <= 1e-8 && (y - 0.0 <= 1e-8))
                    return 3;

                //Outlet
                if (Math.Abs(x - 0.5) <= 1e-8)
                    return 4;

           
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 0.5, 2 * meshScaling + 1);
                var _yNodes = GenericBlas.Linspace(-0.1, 0.1, meshScaling + 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                grd.EdgeTagNames.Add(1, "Velocity_Inlet_O2");
                //grd.EdgeTagNames.Add(2, "Wall_upper");
                grd.EdgeTagNames.Add(2, "wall");
                grd.EdgeTagNames.Add(3, "Velocity_Inlet_CH4");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");
                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // initial values
            // ==============



            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.MixtureFraction, X => 1.0);

        



            // boundary conditions
            // ===================
            // Species 0: Fuel
            // Species 1: Oxidizer
            // Species 2: Product1
            // Species 3: Product2

            double airProportion = 1.0;
            double velMult = 0.3;

            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Velocity_d(0), (X, t) => 1.0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Velocity_d(1), (X, t) => 0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MixtureFraction, (X, t) => 0.0);


            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Velocity_d(0), (X, t) => 1.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Velocity_d(1), (X, t) => 0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MixtureFraction, (X, t) => 1.0);


            C.AddBoundaryValue("wall", VariableNames.Velocity_d(0), (X, t) => 1.0);
            C.AddBoundaryValue("wall", VariableNames.Velocity_d(1), (X, t) => 0);
            //C.AddBoundaryValue("wall", VariableNames.MixtureFraction, (X, t) => 1.0);


            //C.AddBoundaryValue("Wall_upper", VariableNames.Temperature, (X, t) => 1.0);
            //C.AddBoundaryValue("Wall_lower", VariableNames.Temperature, (X, t) => 1.0);

            C.AddBoundaryValue("Pressure_Outlet");

            // Parameters
            // ==============

            C.Reynolds = 1.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;

            C.T_ref_Sutherland = 600;

            C.ReactionRateConstants = new double[] { 1, 3, 1.0, 2.0 };
            C.StoichiometricCoefficients = new double[] { -1, -2, 1, 2 };
            C.M_ref = 16;
            C.PartialHeatCapacities = new double[] { 1.0, 1.0, 1.0, 1.0 };
            return C;
        }




    }
}
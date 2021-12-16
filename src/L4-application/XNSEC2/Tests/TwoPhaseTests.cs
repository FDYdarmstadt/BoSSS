using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    /// </summary>
    static public partial class FullNSEControlExamples {

        /// <summary>
        /// NUnit test
        /// </summary>
        static public XNSEC_Control StillBoxTest_NUnit() {
            var C = StillBoxTest(2, 3);
            C.NoOfMultigridLevels = 1;
            //C.savetodb = false;
            //C.DbPath = null;
            C.savetodb = true;
            C.DbPath = @"C:\Databases\BoSSS_DB";
            C.ChemicalReactionActive = false;
            C.ImmediatePlotPeriod = 1;
            return C;
        }

        /// <summary>
        /// A boring box where nothing is happening
        /// Just two inmiscible fluids without any movement...
        /// </summary>
        static public XNSEC_Control StillBoxTest(int DGp = 2, int nCellsMult = 3, string dbPath = @"C:\Databases\BoSSS_DB") {
            XNSEC_Control C = new XNSEC_Control();

            // Solver configuration
            // ==============

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

            C.DbPath = dbPath;
            C.ProjectName = "StillBox";

            C.physicsMode = PhysicsMode.Combustion;

            C.rhoOne = true;
            C.ChemicalReactionActive = false;

            C.NumberOfChemicalSpecies = 2;
            C.ChemicalReactionActive = false;
            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

            C.PhysicalParameters.IncludeConvection = true;
            C.HeatRelease = 0.0;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;

            // Parameters
            // ==============

            C.Reynolds = 1.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 4.0;
            C.PenaltyHeatConduction = 4.0;

            C.ImmediatePlotPeriod = 1;
            // Grid declaration
            // ===============
            double h = Math.Pow(2, -nCellsMult + 1); // cell length
            double cells = 1 / h;
            int cells2 = (int)cells;

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                return 1; // everywhere wall
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // initial values
            // ==============
            //C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
            //C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.0);
            //C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
            //C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
            //C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);

            // Phase A
            C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#A", X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0) + "#A", X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1) + "#A", X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#A", X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "#A", X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction1 + "#A", X => 0.0);
            //Phase B
            C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#B", X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0) + "#B", X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1) + "#B", X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#B", X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "#B", X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction1 + "#B", X => 1.0);

            C.InitialValues_Evaluators.Add("Phi", X => X[1]);

            // boundary conditions
            // ===================

            C.AddBoundaryValue("wall", VariableNames.Temperature + "#A", (X, t) => 1.0);
            C.AddBoundaryValue("wall", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
            C.AddBoundaryValue("wall", VariableNames.MassFraction1 + "#B", (X, t) => 1.0);

            return C;
        }
    }
}
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    /// </summary>
    public static partial class FullNSEControlExamples {

        ///// <summary>
        ///// NUnit test
        ///// </summary>
        //public static XNSEC_Control ImmersedBoundaryTest_(bool immersedBoundary) {
        //    var C = ImmersedBoundaryTest(DGp: 1, nCellsMult: 4, checkConsistency: false, immersedBoundary);
        //    C.NoOfMultigridLevels = 1;
        //    C.savetodb = false;
        //    C.DbPath = null;
        //    //C.savetodb = true;
        //    //C.DbPath = @"C:\Databases\BoSSS_DB";
        //    C.ChemicalReactionActive = false;
        //    return C;
        //}

        //public static XNSEC_Control ImmersedBoundaryTest(int DGp = 1, int nCellsMult = 3, bool checkConsistency = false, bool immersedBoundary = false) {
        //    XNSEC_Control C = new XNSEC_Control();
        //    Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");
        //    Console.WriteLine("ChannelFlowTest");
        //    Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");

        //    // Solver configuration
        //    // ==============

        //    C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
        //    C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
        //    C.NonLinearSolver.verbose = true;
        //    C.DbPath = null;
        //    C.ProjectName = "ChannelFlowTest";

        //    //C.physicsMode = PhysicsMode.LowMach;
        //    C.physicsMode = PhysicsMode.Combustion;

        //    C.rhoOne = true;
        //    C.AnalyticsolutionSwitch = true;
        //    C.PhysicalParameters.IncludeConvection = false;

        //    C.NumberOfChemicalSpecies = 2;
        //    C.ChemicalReactionActive = false;
        //    C.SetDGdegree(DGp);
        //    C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

        //    C.HeatRelease = 0.0;
        //    // Parameters
        //    // ==============

        //    C.Reynolds = 20.0;
        //    C.Prandtl = 1.0;
        //    C.Schmidt = 1.0;
        //    C.PenaltyViscMomentum = 1.0;
        //    C.PenaltyHeatConduction = 1.0;

        //    C.UseImmersedBoundary = immersedBoundary;
        //    C.ImmediatePlotPeriod = 1;
        //    double L = 1.0;

        //    if (immersedBoundary) {
        //        L = 0.75;
        //        Func<double[], double, double> PhiFunc2 = (X, t) => (X[1] - L);
        //        C.InitialValues_Evaluators_TimeDep.Add("Phi2", PhiFunc2);
        //        C.ThermalParameters.T_sat = 1.3; // boundary temperature
        //    }

        //    //C.ImmediatePlotPeriod = 1;
        //    // Grid declaration
        //    // ===============
        //    double h = Math.Pow(2, -nCellsMult + 1); // cell length
        //    double cells = 1 / h;
        //    int cells2 = (int)cells;

        //    Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
        //        double x = X[0];
        //        double y = X[1];

        //        //Edge tags
        //        //1: Velocity inlet
        //        //2: Wall upper
        //        //3: Pressure outlet
        //        //4: Wall lower

        //        //Inlet oxidizer
        //        if (Math.Abs(x - 0) < 1e-8)
        //            return 1;

        //        //upper Wall
        //        if (Math.Abs(y - 1) < 1e-8)
        //            return 2;

        //        //Outlet
        //        if (Math.Abs(x - 1) < 1e-8)
        //            return 3;

        //        //lower Wall
        //        if (Math.Abs(y + 0) < 1e-8)
        //            return 2;
        //        else throw new ArgumentOutOfRangeException();
        //    };

        //    bool periodic = true;

        //    C.GridFunc = delegate {
        //        var _xNodes = GenericBlas.Linspace(0, 1, cells2 + 1);
        //        var _yNodes = GenericBlas.Linspace(0, 1, (cells2) + 1);
        //        //var _xNodes = GenericBlas.Linspace(0, 10, 2+ 1);
        //        //var _yNodes = GenericBlas.Linspace(-1, 1, 2+ 1);
        //        var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, periodicX: periodic);

        //        if (periodic) {
        //            grd.EdgeTagNames.Add(2, "Velocity_Inlet");
        //        } else {
        //            grd.EdgeTagNames.Add(1, "Velocity_Inlet");
        //            grd.EdgeTagNames.Add(2, "wall");
        //            grd.EdgeTagNames.Add(3, "Pressure_Outlet");
        //        }
        //        grd.DefineEdgeTags(GridEdgeTagFunc);
        //        return grd;
        //    };

        //    C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
        //    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.5);
        //    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
        //    C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 0.5);
        //    C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 0.4);
        //    C.InitialValues_Evaluators.Add(VariableNames.MassFraction1, X => 0.6);

        //    C.InitialValues_Evaluators.Add("Phi", X => -1);
        //    // boundary conditions
        //    // ===================

        //    //C.AddBoundaryValue("wall", VariableNames.Temperature + "#A", (X, t) => 1.0);
        //    //C.AddBoundaryValue("wall", VariableNames.MassFraction0 + "#A", (X, t) => 0.0);

        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(0) + "#A", (X, t) => 0.0);
        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);
        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.Temperature + "#A", (X, t) => 1.0);
        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.Temperature + "#B", (X, t) => 1.0);
        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction0 + "#A", (X, t) => 0.0);
        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction1 + "#A", (X, t) => 1.0);

        //    return C;
        //}

        //public static XNSEC_Control ImmersedBoundaryTest_MixtureFraction(int DGp = 3, int nCellsMult = 3, bool MF = true) {
        //    XNSEC_Control C = new XNSEC_MF_Control();
        //    C.ImmediatePlotPeriod = 1;

        //    //C.physicsMode = PhysicsMode.LowMach;
        //    C.physicsMode = MF ? PhysicsMode.MixtureFraction : PhysicsMode.Combustion;

        //    // Solver configuration
        //    // ==============

        //    C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
        //    C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
        //    C.NonLinearSolver.verbose = true;
        //    C.DbPath = null;

        //    C.savetodb = false;
        //    C.DbPath = null;
        //    //C.savetodb = true;
        //    //C.DbPath = @"C:\Databases\BoSSS_DB";
        //    C.ChemicalReactionActive = false;
        //    //C.physicsMode = PhysicsMode.LowMach;
        //    C.physicsMode = PhysicsMode.MixtureFraction;

        //    C.rhoOne = true;
        //    C.AnalyticsolutionSwitch = true;
        //    C.PhysicalParameters.IncludeConvection = false;

        //    C.NumberOfChemicalSpecies = 2;
        //    C.ChemicalReactionActive = false;
        //    C.SetDGdegree(DGp);
        //    C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };
        //    bool reallyPseudo2D = true;
        //    C.HeatRelease = 0.0;
        //    // Parameters
        //    // ==============

        //    C.Reynolds = 20.0;
        //    C.Prandtl = 1.0;
        //    C.Schmidt = 1.0;
        //    C.PenaltyViscMomentum = 1.0;
        //    C.PenaltyHeatConduction = 1.0;

        //    C.UseImmersedBoundary = true;
        //    C.ImmediatePlotPeriod = 1;
        //    C.HeatRelease = 1;
        //    C.smoothingFactor = 10 * 0;
        //    C.zSt = 0.71;
        //    double L = 3.0;
        //    C.AdaptiveMeshRefinement = true;

        //    C.activeAMRlevelIndicators.Add(new BoSSS.Application.XNSEC.AMR_onFlameSheet(C.zSt, 5));
        //    L = 2.65;
        //    Func<double[], double, double> PhiFunc2 = (X, t) => (X[1] - L);
        //    C.InitialValues_Evaluators_TimeDep.Add("Phi2", PhiFunc2);
        //    C.ThermalParameters.T_sat = 1.3; // boundary temperature

        //    C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
        //    C.dtFixed = double.MaxValue / 1e4;
        //    C.NoOfTimesteps = 5;

        //    //C.ImmediatePlotPeriod = 1;
        //    // Grid declaration
        //    // ===============
        //    double h = Math.Pow(2, -nCellsMult + 1); // cell length
        //    double cells = 1 / h;
        //    int cells2 = (int)cells;

        //    Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
        //        double x = X[0];
        //        double y = X[1];

        //        //Edge tags
        //        //1: Velocity inlet
        //        //2: Wall upper
        //        //3: Pressure outlet
        //        //4: Wall lower

        //        //Inlet oxidizer
        //        if (Math.Abs(x - 0) < 1e-8)
        //            return 1;

        //        //upper Wall
        //        if (Math.Abs(y - 3) < 1e-8)
        //            return 2;

        //        //Outlet
        //        if (Math.Abs(x - 1) < 1e-8)
        //            return 3;

        //        //lower Wall
        //        if (Math.Abs(y + 0) < 1e-8)
        //            return 2;
        //        else throw new ArgumentOutOfRangeException();
        //    };

        //    bool periodic = true;

        //    C.GridFunc = delegate {
        //        var _xNodes = GenericBlas.Linspace(0, 1, reallyPseudo2D ? 4 : cells2 + 1);
        //        var _yNodes = GenericBlas.Linspace(0, 3, (cells2) * 3 + 1);
        //        //var _xNodes = GenericBlas.Linspace(0, 10, 2+ 1);
        //        //var _yNodes = GenericBlas.Linspace(-1, 1, 2+ 1);
        //        var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, periodicX: periodic);

        //        if (periodic) {
        //            grd.EdgeTagNames.Add(2, "Velocity_Inlet");
        //        } else {
        //            grd.EdgeTagNames.Add(1, "Velocity_Inlet");
        //            grd.EdgeTagNames.Add(2, "wall");
        //            grd.EdgeTagNames.Add(3, "Pressure_Outlet");
        //        }
        //        grd.DefineEdgeTags(GridEdgeTagFunc);
        //        return grd;
        //    };

        //    C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
        //    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.5);
        //    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
        //    C.InitialValues_Evaluators.Add(VariableNames.MixtureFraction, X => 0.0);

        //    C.InitialValues_Evaluators.Add("Phi", X => -1);
        //    // boundary conditions
        //    // ===================

        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(0) + "#A", (X, t) => 0.0);
        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);
        //    C.AddBoundaryValue("Velocity_Inlet", VariableNames.MixtureFraction + "#A", (X, t) => 0.0);

        //    return C;
        //}

    
    }
}
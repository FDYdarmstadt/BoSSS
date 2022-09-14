using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Queries;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using static System.Math;


namespace BoSSS.Application.XNSEC {

    static public partial class FullNSEControlExamples {


        static public XNSEC_Control NUnitSteadyCouetteFlowWithTemperatureGradient() {
            var C = SteadyCouetteFlowWithTemperatureGradient(DGp: 2, nCells:16 );
            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.None;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            //C.rhoOne = true;
            C.ImmediatePlotPeriod = 1;

            return C;
        }


        /// <summary>
        /// Convergence study
        /// </summary>
        static public XNSEC_Control[] ConvergenceStudyCouetteTempGradient(string _dbPath = @"C:\Databases\Couette_TempDiffConvStudy_PowerLawMu") {
            List<XNSEC_Control> R = new List<XNSEC_Control>();

            foreach(int DGp in new int[] { 1, 2, 3, 4 }) {
                foreach(int nCells in new int[] { 7 }) {
                    var C = SteadyCouetteFlowWithTemperatureGradient(DGp, nCells);
                    C.Paramstudy_ContinueOnError = true;
                    C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("dgDegree", DGp));
                    C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("refinement", nCells));
                    C.ProjectName = "SteadyCouetteTempGradient. hp-study. p:" + DGp.ToString() + " h:" + (1.0 / (nCells)).ToString();
                    C.savetodb = true;
            
                    C.DbPath = _dbPath;
                    R.Add(C);
                }
            }
            return R.ToArray();
        }


        /// <summary>
        /// SteadyCouetteFlowWithTemperatureGradient
        /// </summary>
        static public XNSEC_Control SteadyCouetteFlowWithTemperatureGradient(int DGp = 2, int nCells = 5, string dbpath = null, bool powerLaw = true) {
            XNSEC_Control C = new XNSEC_Control();
            //Console.WriteLine("SteadyCouetteFlowWithTemperatureGradient");
            C.NumberOfChemicalSpecies = 1;
            C.SetDGdegree(DGp);
            // Solver configuration
            // ==============
            C.DbPath = dbpath; // "D:\\bosss_db"; 
            C.ProjectName = "CouetteFlowTempGrad_DG" + DGp + "K" + nCells;
            C.ProjectDescription = "Steady Low Mach couette flow with temperature gradient";
            C.savetodb = dbpath == null ? false : true;
            C.AnalyticsolutionSwitch = true;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.verbose = true;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.PhysicalParameters.IncludeConvection = false;
            if (powerLaw) { 
            C.MatParamsMode = MaterialParamsMode.PowerLaw;
            } else {
                C.MatParamsMode = MaterialParamsMode.Constant;
            }
            C.ChemicalReactionActive = false;
            C.EnableMassFractions = false;


            //double h = Math.Pow(2, -resolutions + 1); // cell length
            //double cells = 1 / h;
            //int cells2 = (int)cells;
            // Grid declaration
            // ===============

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet O_2
                //2: Wall upper
                //3: Pressure outlet
                //4: Wall lower

                //left Inlet
                if(Math.Abs(x - 0) < 1e-8)
                    return 2;

                //bottom Wall
                if(Math.Abs(y - 0) < 1e-8)
                    return 1;

                // right outlet
                if(Math.Abs(x - 1) < 1e-8)
                    return 3;

                //top Wall
                if(Math.Abs(y - 1) < 1e-8)
                    return 4;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                //var _xNodes = GenericBlas.Linspace(0, 1, nCells + 1);
                //var _yNodes = GenericBlas.Linspace(0, 1, nCells + 1);
                var _xNodes = GenericBlas.Linspace(0, 1, nCells + 1);
                var _yNodes = GenericBlas.Linspace(0, 1, (nCells) + 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                //grd.EdgeTagNames.Add(1, "wall_bottom");
                grd.EdgeTagNames.Add(1, "velocity_inlet_bottom");
                grd.EdgeTagNames.Add(2, "velocity_inlet_left");
                grd.EdgeTagNames.Add(3, "velocity_inlet_right");
                grd.EdgeTagNames.Add(4, "velocity_inlet_top"); // moving wall

                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // Parameters
            // ==============
            C.Reynolds = 10.0;
            C.Prandtl = 0.71;
            C.Froude = Math.Sqrt(2 * C.Prandtl * (1.6 - 0.4) / (1.6 + 0.4));
            C.PressureReferencePoint = new double[] { 0.5*0, 0.5 *0};

            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;

            C.GravityDirection = new double[] { 0.0, 1.0, 0.0 };

            double Th = 1.6; // Adimensional temperature of the hot wall (moving wall)
            double Tc = 0.4; // Adimensional temperature of the cold wall

            double p0 = 1.0; // The thermodynamic pressure, constant for open systems
            double c_1 = p0 / (C.Froude * C.Froude * (Th - Tc)) * Math.Log((Th - Tc) * C.PressureReferencePoint[1] + Tc) * 1; // A constant dependent on the reference point for the pressure

            #region ViscosityVariable_PowerLaw

            double TC53 = Pow(Tc, 5.0 / 3.0);
            double TH53 = Pow(Th, 5.0 / 3.0);

            double THTC = TH53 - TC53;
            double NUM = Pow(TC53 / THTC, 3.0 / 5.0);

            double DENOM = Pow(TC53 / THTC, 3.0 / 5.0) - Pow(TH53 / THTC, 3.0 / 5.0);
            double DENOM2 = Pow(TH53 / THTC, 3.0 / 5.0) - Pow(TC53 / THTC, 3.0 / 5.0);

            double C1 = NUM / DENOM;
            double C2 = 1.0 / (DENOM2);
            double C4 = Pow(Tc, 5.0 / 3.0);
            double C5 = 3.0 / 5.0 * (TC53 - TH53);
            Func<double[], double> func = X => 5.0 / 2.0 * p0 / (C.Froude * C.Froude * THTC) * Pow((X[1] * THTC + TC53), 2.0 / 5.0);

            double C3 = (-1) * func(C.PressureReferencePoint); // A constant dependent on the reference point for the pressure

            #endregion ViscosityVariable_PowerLaw

            // Analytic solutions
            // ===================
            if(C.MatParamsMode == MaterialParamsMode.Constant) {
                C.AnalyticVelocityX = X => X[1];
                C.AnalyticVelocityY = X => 0.0;
                C.AnalyticPressure = X => (-1.0) * p0 / (C.Froude * C.Froude * (Th - Tc)) * Log((Th - Tc) * X[1] + Tc) + c_1;
                C.AnalyticTemperature = X => (Th - Tc) * X[1] + Tc;
            } else if(C.MatParamsMode == MaterialParamsMode.PowerLaw) {
                C.AnalyticVelocityX = X => C1 + C2 * Pow((X[1] + TC53 / THTC), 3.0 / 5.0);
                C.AnalyticVelocityY = X => 0.0;
                C.AnalyticPressure = X => -5.0 / 2.0 * p0 / (C.Froude * C.Froude * THTC) * Pow((X[1] * THTC + TC53), 2.0 / 5.0) + C3;
                C.AnalyticTemperature = X => Pow(C4 - 5.0 / 3.0 * C5 * X[1], 3.0 / 5.0);
            } else if(C.MatParamsMode == MaterialParamsMode.Sutherland) { // NOT WORKING???
                C.AnalyticVelocityX = X => -0.33333333333333333332 + 1.2523108062960316903 * (X[1] + 0.11013981986815616002).Pow(3.0 / 5.0);
                C.AnalyticVelocityY = X => 0.0;
                C.AnalyticPressure = X => -1.4882576522041879973 * (1.9716158024185796207 * X[1] + 0.21715340932759252572).Pow(2.0 / 5.0) + 1.5508248735452649660;
                C.AnalyticTemperature = X => (1.9716158024185796207 * X[1] + 0.21715340932759252572).Pow(3.0 / 5.0);
            } else { throw new NotImplementedException("Not implemented matparammode"); }

            // Analytical / manufactured solutions

            C.InitialValues_Evaluators.Add(VariableNames.VelocityX + "_an", C.AnalyticVelocityX);
            C.InitialValues_Evaluators.Add(VariableNames.VelocityY + "_an", C.AnalyticVelocityY);
            C.InitialValues_Evaluators.Add(VariableNames.Pressure + "_an", C.AnalyticPressure);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature + "_an", C.AnalyticTemperature);

            // initial values
            // ==============
            double mult1 = 0.5;
            double mult2 = 0.5;
            C.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => (X[1]) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => ((-1) * p0 / (C.Froude * C.Froude * (Th - Tc)) * Log((Th - Tc) * X[1] + Tc) + c_1) * mult1);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => ((Th - Tc) * X[1] + Tc) * mult2);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);

            // boundary conditions
            // ===================

            C.AddBoundaryValue("velocity_inlet_top", VariableNames.VelocityX, C.AnalyticVelocityX);
            C.AddBoundaryValue("velocity_inlet_top", VariableNames.VelocityY, C.AnalyticVelocityY);
            C.AddBoundaryValue("velocity_inlet_top", VariableNames.Temperature, C.AnalyticTemperature);
            C.AddBoundaryValue("velocity_inlet_top", VariableNames.MassFraction0, X => 1.0);
            C.AddBoundaryValue("velocity_inlet_left", VariableNames.VelocityX, C.AnalyticVelocityX);
            C.AddBoundaryValue("velocity_inlet_left", VariableNames.VelocityY, C.AnalyticVelocityY);
            C.AddBoundaryValue("velocity_inlet_left", VariableNames.Temperature, C.AnalyticTemperature);
            C.AddBoundaryValue("velocity_inlet_left", VariableNames.MassFraction0, X => 1.0);
            C.AddBoundaryValue("velocity_inlet_right", VariableNames.VelocityX, C.AnalyticVelocityX);
            C.AddBoundaryValue("velocity_inlet_right", VariableNames.VelocityY, C.AnalyticVelocityY);
            C.AddBoundaryValue("velocity_inlet_right", VariableNames.Temperature, C.AnalyticTemperature);
            C.AddBoundaryValue("velocity_inlet_right", VariableNames.MassFraction0, X => 1.0);

            C.AddBoundaryValue("velocity_inlet_bottom", VariableNames.VelocityX, C.AnalyticVelocityX);
            C.AddBoundaryValue("velocity_inlet_bottom", VariableNames.VelocityY, C.AnalyticVelocityY);
            C.AddBoundaryValue("velocity_inlet_bottom", VariableNames.Temperature, C.AnalyticTemperature);
            C.AddBoundaryValue("velocity_inlet_bottom", VariableNames.MassFraction0, X => 1.0);

            //C.AddBoundaryValue("wall_bottom", VariableNames.Temperature, C.AnalyticTemperature);

            int queryQuadOrder = 10;
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


        ///// <summary>
        ///// SteadyCouetteFlowWithTemperatureGradient
        ///// </summary>
        //static public XNSEC_Control SteadyCouetteFlowWithTemperatureGradientTEST(int DGp = 2, int resolutions = 4, string dbpath = @"C:\Databases\BoSSS_DB") {
        //    XNSEC_Control C = new XNSEC_Control();
        //    Console.WriteLine("SteadyCouetteFlowWithTemperatureGradient");
        //    C.NumberOfChemicalSpecies = 1;
        //    C.SetDGdegree(DGp);
        //    // Solver configuration
        //    // ==============
        //    C.DbPath = dbpath; // "D:\\bosss_db";
        //    C.ProjectName = "CouetteFlowTempGrad_DG" + DGp + "K" + resolutions;
        //    C.ProjectDescription = "Steady Low Mach couette flow with temperature gradient";
        //    C.savetodb = true;// dbpath == null ? false : true;
        //    C.AnalyticsolutionSwitch = true;
        //    C.physicsMode = PhysicsMode.Combustion;
        //    //C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
        //    //C.dtFixed = 1.0;
        //    //C.NoOfTimesteps = 100;
        //    //C.Endtime = 10000;
        //    C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

        //    C.NonLinearSolver.verbose = true;
        //    C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
        //    C.NonLinearSolver.ConvergenceCriterion = 1e-6;
        //    C.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();
        //    C.MatParamsMode = MaterialParamsMode.Constant;

        //    double h = Math.Pow(2, -resolutions + 1); // cell length
        //    double cells = 1 / h;
        //    int cells2 = (int)cells;

        //    // Grid declaration
        //    // ===============
        //    int mult = 10;
        //    Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
        //        double x = X[0];
        //        double y = X[1];

        //        //Edge tags
        //        //1: Velocity inlet O_2
        //        //2: Wall upper
        //        //3: Pressure outlet
        //        //4: Wall lower

        //        //left Inlet
        //        if(Math.Abs(x - 0) < 1e-8)
        //            return 2;

        //        //bottom Wall
        //        if(Math.Abs(y - 0) < 1e-8)
        //            return 1;

        //        // right outlet
        //        if(Math.Abs(x - 1 * mult) < 1e-8)
        //            return 3;

        //        //top Wall
        //        if(Math.Abs(y - 1) < 1e-8)
        //            return 4;
        //        else throw new ArgumentOutOfRangeException();
        //    };

        //    C.GridFunc = delegate {
        //        //var _xNodes = GenericBlas.Linspace(0, 1, nCells + 1);
        //        //var _yNodes = GenericBlas.Linspace(0, 1, nCells + 1);
        //        var _xNodes = GenericBlas.Linspace(0, 1 * mult, cells2 * mult + 1);
        //        var _yNodes = GenericBlas.Linspace(0, 1, (cells2) + 1);
        //        var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
        //        grd.EdgeTagNames.Add(1, "wall_bottom");
        //        grd.EdgeTagNames.Add(2, "velocity_inlet_left");
        //        grd.EdgeTagNames.Add(3, "Pressure_Outlet");
        //        grd.EdgeTagNames.Add(4, "velocity_inlet_top"); // moving wall

        //        grd.DefineEdgeTags(GridEdgeTagFunc);
        //        return grd;
        //    };

        //    // Parameters
        //    // ==============
        //    C.Reynolds = 10.0;
        //    C.Prandtl = 0.71;
        //    C.Froude = Math.Sqrt(2 * C.Prandtl * (1.6 - 0.4) / (1.6 + 0.4));
        //    C.PressureReferencePoint = new double[] { 0.5, 0.5 };

        //    C.Schmidt = 1.0;
        //    C.PenaltyViscMomentum = 1.0;
        //    C.PenaltyHeatConduction = 1.0;

        //    C.GravityDirection = new double[] { 0.0, 1.0, 0.0 };

        //    double Th = 1.6; // Adimensional temperature of the hot wall (moving wall)
        //    double Tc = 0.4; // Adimensional temperature of the cold wall

        //    double p0 = 1.0; // The thermodynamic pressure, constant for open systems
        //    double c_1 = p0 / (C.Froude * C.Froude * (Th - Tc)) * Math.Log((Th - Tc) * C.PressureReferencePoint[1] + Tc) * 1; // A constant dependent on the reference point for the pressure

        //    #region ViscosityVariable_PowerLaw

        //    double TC53 = Pow(Tc, 5.0 / 3.0);
        //    double TH53 = Pow(Th, 5.0 / 3.0);

        //    double THTC = TH53 - TC53;
        //    double NUM = Pow(TC53 / THTC, 3.0 / 5.0);

        //    double DENOM = Pow(TC53 / THTC, 3.0 / 5.0) - Pow(TH53 / THTC, 3.0 / 5.0);
        //    double DENOM2 = Pow(TH53 / THTC, 3.0 / 5.0) - Pow(TC53 / THTC, 3.0 / 5.0);

        //    double C1 = NUM / DENOM;
        //    double C2 = 1 / (DENOM2);
        //    double C4 = Pow(Tc, 5.0 / 3.0);
        //    double C5 = 3.0 / 5.0 * (TC53 - TH53);
        //    Func<double[], double> func = X => 5 / 2 * p0 / (C.Froude * C.Froude * THTC) * Pow((X[1] * THTC + TC53), 2.0 / 5.0);

        //    double C3 = (-1) * func(C.PressureReferencePoint); // A constant dependent on the reference point for the pressure

        //    #endregion ViscosityVariable_PowerLaw

        //    // Analytic solutions
        //    // ===================
        //    Func<double[], double> AnalyticVelocityX;
        //    Func<double[], double> AnalyticVelocityY;
        //    Func<double[], double> AnalyticPressure;
        //    Func<double[], double> AnalyticTemperature;
        //    if(C.MatParamsMode == MaterialParamsMode.Constant) {
        //        AnalyticVelocityX = X => X[1];
        //        AnalyticVelocityY = X => 0.0;
        //        AnalyticPressure = X => (-1) * p0 / (C.Froude * C.Froude * (Th - Tc)) * Log((Th - Tc) * X[1] + Tc) + c_1;
        //        AnalyticTemperature = X => (Th - Tc) * X[1] + Tc;
        //    } else if(C.MatParamsMode == MaterialParamsMode.PowerLaw) {
        //        AnalyticVelocityX = X => C1 + C2 * Pow((X[1] + TC53 / THTC), 3.0 / 5.0);
        //        AnalyticVelocityY = X => 0.0;
        //        AnalyticPressure = X => -5 / 2 * p0 / (C.Froude * C.Froude * THTC) * Pow((X[1] * THTC + TC53), 2.0 / 5.0) + C3;
        //        AnalyticTemperature = X => Pow(C4 - 5.0 / 3.0 * C5 * X[1], 3.0 / 5.0);
        //    } else if(C.MatParamsMode == MaterialParamsMode.Sutherland) { // NOT WORKING???
        //        AnalyticVelocityX = X => -0.33333333333333333332 + 1.2523108062960316903 * (X[1] + 0.11013981986815616002).Pow(3.0 / 5.0);
        //        AnalyticVelocityY = X => 0.0;
        //        AnalyticPressure = X => -1.4882576522041879973 * (1.9716158024185796207 * X[1] + 0.21715340932759252572).Pow(2.0 / 5.0) + 1.5508248735452649660;
        //        AnalyticTemperature = X => (1.9716158024185796207 * X[1] + 0.21715340932759252572).Pow(3.0 / 5.0);
        //    } else {
        //        throw new NotImplementedException("Not implemented matparammode");
        //    }

        //    // Analytical / manufactured solutions

        //    C.InitialValues_Evaluators.Add(VariableNames.VelocityX + "_an", AnalyticVelocityX);
        //    C.InitialValues_Evaluators.Add(VariableNames.VelocityY + "_an", AnalyticVelocityY);
        //    C.InitialValues_Evaluators.Add(VariableNames.Pressure + "_an", AnalyticPressure);
        //    C.InitialValues_Evaluators.Add(VariableNames.Temperature + "_an", AnalyticTemperature);

        //    // initial values
        //    // ==============
        //    double mult1 = 0.1;
        //    double mult2 = 0.1;
        //    C.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => (X[1]) * mult1);
        //    C.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
        //    C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => ((-1) * p0 / (C.Froude * C.Froude * (Th - Tc)) * Log((Th - Tc) * X[1] + Tc) + c_1) * mult1);
        //    C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => ((Th - Tc) * X[1] + Tc) * mult2);
        //    C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);

        //    // boundary conditions
        //    // ===================

        //    C.AddBoundaryValue("velocity_inlet_top", VariableNames.VelocityX, AnalyticVelocityX);
        //    C.AddBoundaryValue("velocity_inlet_top", VariableNames.VelocityY, AnalyticVelocityY);
        //    C.AddBoundaryValue("velocity_inlet_top", VariableNames.Temperature, AnalyticTemperature);
        //    C.AddBoundaryValue("velocity_inlet_top", VariableNames.MassFraction0, X => 1.0);

        //    C.AddBoundaryValue("velocity_inlet_left", VariableNames.VelocityX, AnalyticVelocityX);
        //    C.AddBoundaryValue("velocity_inlet_left", VariableNames.VelocityY, AnalyticVelocityY);
        //    C.AddBoundaryValue("velocity_inlet_left", VariableNames.Temperature, AnalyticTemperature);
        //    C.AddBoundaryValue("velocity_inlet_left", VariableNames.MassFraction0, X => 1.0);
        //    C.AddBoundaryValue("Pressure_Outlet");

        //    C.AddBoundaryValue("wall_bottom", VariableNames.Temperature, AnalyticTemperature);

        //    int queryQuadOrder = 1;
        //    C.Queries.Add(
        //        "SolL2err_u",
        //        QueryLibrary.L2Error(VariableNames.VelocityX, AnalyticVelocityX, queryQuadOrder));
        //    C.Queries.Add(
        //        "SolL2err_v",
        //        QueryLibrary.L2Error(VariableNames.VelocityY, AnalyticVelocityY, queryQuadOrder));
        //    C.Queries.Add(
        //        "SolL2err_p",
        //        QueryLibrary.L2Error(VariableNames.Pressure, AnalyticPressure, queryQuadOrder));
        //    C.Queries.Add(
        //        "SolL2err_T",
        //        QueryLibrary.L2Error(VariableNames.Temperature, AnalyticTemperature, queryQuadOrder));

        //    return C;
        //}

    }
}
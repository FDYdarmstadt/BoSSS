using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using static System.Math;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    /// </summary>
    static public partial class FullNSEControlExamples {
        static public XNSEC_Control NaturalConvectionSquareCavityTest_Homotopy() {
            var C = NaturalConvectionSquareCavity(2, 20, 1e4, -1, null);
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined;
            C.savetodb = false;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.MaxSolverIterations = 1500;
            C.myThermalWallType = SIPDiffusionTemperature.ThermalWallType.fixedTemperature;
            C.TRef = 600;
            C.UseSelfMadeTemporalOperator = false;
            C.timeDerivativeEnergyp0_OK = false;
            C.timeDerivativeConti_OK = false;
            //C.HomotopyApproach = XNSEC_Control.HomotopyType.Automatic;

            C.LinearSolver = Solution.Control.LinearSolverCode.exp_Kcycle_schwarz.GetConfig();
            //C.LinearSolver.NoOfMultigridLevels = 10;
            C.HomotopyVariable = XNSEC_Control.HomotopyVariableEnum.Reynolds;
            C.homotopieAimedValue = Math.Sqrt(1e6);
            C.StartingHomotopyValue = Math.Sqrt(1e4);
            C.HomotopyApproach = XNSEC_Control.HomotopyType.Automatic;
            C.AdaptiveMeshRefinement = false;
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnResiduals(_maxRefinementLevelval: 2,_tresholdValue:0.1));
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnResidual_GRADIENT(3, _tresh: 5));
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnPerssonSensor(VariableNames.Temperature, 2));



            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.MaxSolverIterations = 1000;
            

            //C.NonLinearSolver.Globalization = Solution.AdvancedSolvers.Newton.GlobalizationOption.LineSearch;
            C.NonLinearSolver.Globalization = Solution.AdvancedSolvers.Newton.GlobalizationOption.Dogleg;
            return C;
        }

        static public XNSEC_Control NaturalConvectionSquareCavityTest() {
            var C = NaturalConvectionSquareCavity(2, 10, 1e3, -1, null);
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined;
            C.savetodb = false;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.TRef = 600;
            C.UseSelfMadeTemporalOperator = false;
            C.timeDerivativeEnergyp0_OK = false;
            C.timeDerivativeConti_OK = false;
            C.PhysicalParameters.IncludeConvection = true;
            C.ImmediatePlotPeriod = 1;
            return C;
        } 

        static public XNSEC_Control NaturalConvectionSquareCavityTestUnstdy() {
            var C = NaturalConvectionSquareCavity(1, 3, 1e1, 1, @"C:\Databases\BoSSS_DB");
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined;
            C.TRef = 600;
            C.savetodb = true;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.NonLinearSolver.MaxSolverIterations = 40;
            //C.PlotNewtonIterations = true;
            C.EnableMassFractions = false;
            C.timeDerivativeConti_OK = true;
            //C.timeDerivativeEnergyp0_OK = true;
            C.UseSelfMadeTemporalOperator =  false;
            C.timeDerivativeEnergyp0_OK = false; //////////////////////////////////////////////////////////

            return C;
        }
        static public XNSEC_Control CavityNaturalConvection_AMR_eachNewtonIteration() {
            var C = NaturalConvectionSquareCavity(2, 4, 1e3, -1, @"C:\Databases\BoSSS_DB");
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined;
            C.savetodb = true;
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;
            C.TRef = 600;
            C.UseSelfMadeTemporalOperator = false;
            C.timeDerivativeEnergyp0_OK = false;
            C.timeDerivativeConti_OK = false;
            C.EnableMassFractions = false;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NoOfTimesteps = 10;
            C.NonLinearSolver.MaxSolverIterations = 1; // Do only one newton iteration before refining
            C.NonLinearSolver.MinSolverIterations = 1; // Do only one newton iteration before refining

            C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnPerssonSensor(VariableNames.Temperature,4));
            C.activeAMRlevelIndicators.Add(new AMR_BasedOnVariableLimits(VariableNames.Temperature, new double[] { 0.4, 1.6 },3));
            C.AMR_startUpSweeps = 2;

            return C;
        }
        /// <summary>
        /// ConvergenceStudy NaturalConvectionSquareCavity
        /// </summary>
        static public XNSEC_Control[] ConvergenceStudyNaturalConvectionSquareCavity() {
            List<XNSEC_Control> R = new List<XNSEC_Control>();
            
            foreach(int DGp in new int[] { /*1, 2,*/ 3/*, 4 */}) {
                foreach(int SizeFactor in new int[] { 2, 3, 4, 5 }) {
                    //foreach(int DGp in new int[] { 1, 2, 3, 4 }) {
                    //    foreach(int SizeFactor in new int[] { 1, 2 }) {
                    string dbPath = @"C:\Databases\NaturalConvection_XNSEC_ConvStudy_OnlyHeatDiffusion_ConstantHeatDiff";
                    var C = NaturalConvectionSquareCavity(DGp, (int)Math.Pow(2, SizeFactor + 1), 1e3, -1, dbPath);

                    C.rhoOne = false;
                    C.GravityDirection = new double[] { 0, 0, 0 };
                    C.MatParamsMode = MaterialParamsMode.Constant;

                    C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("dgDegree", DGp));
                    C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("refinements", SizeFactor));

                    C.ProjectName = "SquareCavity. hp-study. p:" + DGp.ToString() + " h:" + (1.0 / (Math.Pow(2, SizeFactor))).ToString();
                    C.savetodb = true;

                    R.Add(C);
                }
            }









            return R.ToArray();
        }
        /// <summary>
        /// Control file for calculation of the heated cavity for increasing Rayleigh number using homotopy
        /// </summary>
        static public XNSEC_Control NaturalConvectionSquareCavityRayleighHomotopy() {
            //string dbPath = @"C:\Databases\HomotopyStudyNatConv";
            var C = NaturalConvectionSquareCavity(2, 20, 1e7, -1, null);
            C.homotopieAimedValue = Math.Sqrt(1e7);
            C.savetodb = true;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.DbPath = @"C:\Databases\NatConvectionP0Study";
            IEnumerable<double> RayleighHomotopyArray = new double[] { 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 1e7 };

            C.HomotopyArray = RayleighHomotopyArray.Select(x => Math.Sqrt(x)).ToArray();
            C.HomotopyApproach = XNSEC_Control.HomotopyType.Manual;
            C.NoOfTimesteps = RayleighHomotopyArray.ToArray().Length;
            return C;
        }


        static public XNSEC_Control NaturalConvectionSquareCavity_ForSolvingConvercenceProblem(int dgp = 2, int ncells = 4, double Ra = 1e3, double dt = -1, string dbpath = @"\\hpccluster\hpccluster-scratch\gutierrez\rhoOne_1muConst0", int VariableTransportParameters = 0, int VariableRho = 1, int p0MassDetermined = 0, double T_left = 1.6, double T_right = 0.4) {
            var C = NaturalConvectionSquareCavity(dgp, ncells, Ra, dt, dbpath, T_left, T_right);
            if(VariableTransportParameters == 0) {
                C.MatParamsMode = MaterialParamsMode.Constant;
            } else if(VariableTransportParameters == 1) {
                C.MatParamsMode = MaterialParamsMode.Sutherland;
            } else {
                throw new Exception("wrong ");
            }

            C.EnableMassFractions = false;
            if(VariableRho == 0) {
                C.rhoOne = true;
                C.GravityDirection = new double[] { 0, 0, 0 };
            } else if(VariableRho == 1) {
                C.rhoOne = false;
            } else {
                throw new Exception("wrong ");
            }

            if(p0MassDetermined == 0) {
                C.ThermodynamicPressureMode = ThermodynamicPressureMode.Constant;
            } else if(p0MassDetermined == 1) {
                C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined;
            } else {
                throw new Exception("wrong ");
            }
            return C;
        }

        /// <summary>
        /// Test case for NaturalConvection in a Square Cavity due temperature differences and density changes...
        /// </summary>
        static public XNSEC_Control NaturalConvectionSquareCavity(int DGp = 2, int ncells = 30, double Ra = 1e3, double dt = -1, string dbpath = @"C:\Databases\BoSSS_DB", double T_left = 1.6, double T_right = 0.4, double penalty = 1.0) {

            XNSEC_Control C = new XNSEC_Control();

            //Console.WriteLine("Control File: NaturalConvectionSquareCavity. \t Ra = {0} ", Ra);
            //Console.WriteLine("DG degree:{0} \t Number of cells: {1} ", DGp, ncells);
            C.DbPath = dbpath;
            //dbpath = @"C:\BoSSS_DB";

            //C.DbPath = dbpath;
            C.savetodb = dbpath == null ? false : true;
            C.EnableMassFractions = false;
            C.NumberOfChemicalSpecies = 1;
            C.ChemicalReactionActive = false;
            C.MatParamsMode = MaterialParamsMode.Sutherland;
            C.physicsMode = PhysicsMode.Combustion;
            //C.BDFOrder = 1;

            ////////////
            // Solver configuration
            // ==============
            C.SetDGdegree(DGp);

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.ConvergenceCriterion = 1e-11;
            C.NonLinearSolver.verbose = true;
            C.NonLinearSolver.MaxSolverIterations = 50;

            C.PenaltyViscMomentum = 1.0 * penalty;
            C.PenaltyHeatConduction = 1.0 * penalty;
            C.PhysicalParameters.IncludeConvection = true;

            C.UseSelfMadeTemporalOperator = false;
            C.timeDerivativeEnergyp0_OK = false;
            C.timeDerivativeConti_OK = false;

            if(dt < 0) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                C.dtFixed = dt;
                C.Endtime = double.MaxValue;
                C.NoOfTimesteps = int.MaxValue;
            }

            // Homotopy
            //C.SetAddaptiveMeshRefinement(1, 5);

            //if (numOfSteps == 1) {
            //    C.useHomotopie = false;
            //} else if (numOfSteps > 1) {
            //    C.useHomotopie = true;
            //    C.homotopieVariable = LowMachCombustionNSEControl.HomotopieVariableNames.Reynolds;
            //    C.homotopieAimedValue = Math.Sqrt(Ra);
            //    C.m_HomotopyArray = new double[] { 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7 };
            //} else {
            //    throw new NotImplementedException("ups");
            //}

            bool restartOK = false;

            double xyRatio = 1;
            double stretchfactorX = 0.95*0;
            double stretchfactorY = 0.95*0;

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Adiabatic no slip wall
                //2: Temperature fixed no slip wall

                //right cold wall
                if(Math.Abs(x - 0.5) < 1e-8)
                    return 3;

                //bottom adiabatic Wall
                if(Math.Abs(y - 0.5 * xyRatio) < 1e-8)
                    return 1;

                // left hot wall
                if(Math.Abs(x + 0.5) < 1e-8)
                    return 2;

                //top adiabatic Wall
                if(Math.Abs(y + 0.5 * xyRatio) < 1e-8)
                    return 1;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.SinLinSpacing(-0.5, 0.5, stretchfactorX, ncells + 1);
                var _yNodes = GenericBlas.SinLinSpacing(-0.5 * xyRatio, 0.5 * xyRatio, stretchfactorY, ncells + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes/*, periodicY: true*/);
                grd.EdgeTagNames.Add(1, "NoSlipNeumann");
                grd.EdgeTagNames.Add(2, "wall_tempfixed_left");
                grd.EdgeTagNames.Add(3, "wall_tempfixed_right");
                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // ==============
            C.EdgeTagsNusselt = new string[] { "wall_tempfixed_left", "wall_tempfixed_right", "NoSlipNeumann" };

            C.Rayleigh = Ra;
            C.Reynolds = Sqrt(Ra);
            C.Prandtl = 0.71;
            C.Froude = Math.Sqrt(2 * C.Prandtl * (1.6 - 0.4) / (1.6 + 0.4));
            C.HeatCapacityRatio = 1.4;

            C.T_ref_Sutherland = 600;
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined; // Because its a closed system, i.e. p0 = p0(time)
            // C._MyThermodynamicPressureMode = LowMachCombustionNSEControl.MyThermodynamicPressureMode.MassDetermined ; // Because its a closed system, i.e. p0 = p0(time)

         
            double Th = T_left;// 1.6; // Adimensional temperature of the hot wall
            double Tc = T_right;// 0.4; // Adimensional temperature of the cold wall

            if(!restartOK) {
                // initial values
                // ==============
                C.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => (-1) * X[1] / (C.Froude * C.Froude));

                C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);

                //C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => (Tc - Th) / 1 * X[0] + Th);

                //C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => X[0] * X[0] + X[1] * X[1] + 1);
                C.InitialValues_Evaluators.Add(VariableNames.ThermodynamicPressure, X => 1);
            } else {
                C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(C.ID), -1);
            }

            // boundary conditions
            // ===================

            C.AddBoundaryValue("NoSlipNeumann", VariableNames.VelocityX, X => 0.0);
            C.AddBoundaryValue("NoSlipNeumann", VariableNames.VelocityY, X => 0.0);

            C.AddBoundaryValue("wall_tempfixed_left", VariableNames.Temperature, X => Th);
            C.AddBoundaryValue("wall_tempfixed_right", VariableNames.Temperature, X => Tc);

            C.AddBoundaryValue("wall_tempfixed_left", VariableNames.MassFraction0, X => 1.0);
            C.AddBoundaryValue("wall_tempfixed_right", VariableNames.MassFraction0, X => 1.0);
            return C;
        }

        /// <summary>
        ///  Thermodynamic pressure is calculated for a given temperature profile.
        ///  TODO
        /// </summary>
        static public XNSEC_Control ThermodynamicPressureTest(int DGp = 2, int ncells = 8, double Ra = 1e3, double dt = -1, string dbpath = @"C:\Databases\BoSSS_DB") {
            XNSEC_Control C = new XNSEC_Control();

            C.DbPath = dbpath;
            C.savetodb = dbpath == null ? false : true;
            C.NumberOfChemicalSpecies = 1;
            C.ChemicalReactionActive = false;
            C.SkipSolveAndEvaluateResidual = true;

            ////////////
            // Solver configuration
            // ==============
            C.SetDGdegree(DGp);
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];
                if(Math.Abs(x - 0.5) < 1e-8)
                    return 1;
                if(Math.Abs(y - 0.5) < 1e-8)
                    return 1;
                if(Math.Abs(x + 0.5) < 1e-8)
                    return 1;
                if(Math.Abs(y + 0.5) < 1e-8)
                    return 1;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-0.5, 0.5, ncells + 1);
                var _yNodes = GenericBlas.Linspace(-0.5, 0.5, ncells + 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes/*, periodicY: true*/);
                grd.EdgeTagNames.Add(1, "wall");
                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined; // Because its a closed system, i.e. p0 = p0(time)

            C.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => (-1) * X[1] / (C.Froude * C.Froude));
            C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
            //C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => X[0] * X[0] + X[1] * X[1] + 1);
            C.InitialValues_Evaluators.Add(VariableNames.ThermodynamicPressure, X => 1);

            // boundary conditions
            // ===================

            C.AddBoundaryValue("wall", VariableNames.Temperature, X => 1.0);
            C.AddBoundaryValue("wall", VariableNames.MassFraction0, X => 1.0);

            return C;
        }

        static public XNSEC_Control HeatedBoxTest() {
            var C = HeatedBox(1, 4, 1.5, 0.005, @"C:\Databases\BoSSS_DB");
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined;
            //C.GravityDirection = new double[] { 0, 1, 0 };
            C.NonLinearSolver.ConvergenceCriterion = 1e-6;

            C.timeDerivativeConti_OK = false;
            C.timeDerivativeEnergyp0_OK = true;

            C.rhoOne = false;
            C.PhysicalParameters.IncludeConvection = true;
            C.PlotNewtonIterations = false;
            C.EnableMassFractions = false;

            //C.NonLinearSolver.MaxSolverIterations = 1;

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            return C;
        }

        /// <summary>
        ///
        /// </summary>
        static public XNSEC_Control HeatedBox(int DGp = 1, int ncells = 4, double wallTemperature = 2, double dt = -1, string dbpath = @"C:\Databases\BoSSS_DB") {
            XNSEC_Control C = new XNSEC_Control();

            C.DbPath = dbpath;
            C.savetodb = dbpath == null ? false : true;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMR_BasedOnFieldGradient() { maxRefinementLevel = 4, FieldName = VariableNames.Temperature });

            C.NumberOfChemicalSpecies = 1;
            C.ChemicalReactionActive = false;
            C.MatParamsMode = MaterialParamsMode.Sutherland;

            // ====================
            // Solver configuration
            // ====================
            C.SetDGdegree(DGp);

            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.verbose = true;

            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;

            if(dt < 0) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                C.NoOfTimesteps = 100;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                C.dtFixed = dt;
                C.Endtime = double.MaxValue;
                C.NoOfTimesteps = int.MaxValue;
            }

            bool restartOK = false;

            double xyRatio = 1;
            double stretchfactorX = 0.9;
            double stretchfactorY = 0.9;

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //right cold wall
                if(Math.Abs(x - 0.5) < 1e-8)
                    return 1;

                //bottom adiabatic Wall
                if(Math.Abs(y - 0.5 * xyRatio) < 1e-8)
                    return 1;

                // left hot wall
                if(Math.Abs(x + 0.5) < 1e-8)
                    return 1;

                //top adiabatic Wall
                if(Math.Abs(y + 0.5 * xyRatio) < 1e-8)
                    return 1;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.SinLinSpacing(-0.5, 0.5, stretchfactorX, ncells + 1);
                var _yNodes = GenericBlas.SinLinSpacing(-0.5 * xyRatio, 0.5 * xyRatio, stretchfactorY, ncells + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes/*, periodicY: true*/);
                grd.EdgeTagNames.Add(1, "WallHeated");
                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // ==============
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };
            C.Reynolds = 1;
            C.Prandtl = 0.71;

            C.ThermodynamicPressureMode = ThermodynamicPressureMode.Constant; // Because its a closed system, i.e. p0 = p0(time)

            // initial values
            // ==============
            C.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.ThermodynamicPressure, X => 1);

            // boundary conditions
            // ===================

            C.AddBoundaryValue("WallHeated", VariableNames.Temperature, X => wallTemperature);

            return C;
        }

        /// <summary>
        /// Test case for NaturalConvection in a Tall Cavity due temperature differences and density changes...
        /// </summary>
        static public XNSEC_Control NaturalConvectionTallCavity(int DGp = 5, int ncells = 26, int constantNewtonIterations = 1, double Ra = 1e6, double dt = 0.05, int numOfSteps = 1, int SolverSwitch = 1, int FDJOK = 1, double nonLinConvCrit = 1e-5, string dbPath = null) {
            XNSEC_Control C = new XNSEC_Control();

            Console.WriteLine("Control File: NaturalConvectionSquareCavity. \t Ra = {0} ", Ra);
            Console.WriteLine("DG degree:{0} \t Number of cells: {1} ", DGp, ncells);

            // Solver configuration
            // ==============
            C.SetDGdegree(DGp);
            C.NonLinearSolver.constantNewtonIterations = constantNewtonIterations; // Number of iterations in the newton algorithm where the Jacobi matrix stays constant
                                                                                   // C.NonLinearSolver.MaxSolverIterations = 40;
            C.NonLinearSolver.ConvergenceCriterion = nonLinConvCrit;
            //C.NonLinearSolver.ConvergenceCriterion = 1e-5;
            C.NonLinearSolver.verbose = true;

            //C.NonLinearSolver.PrecondSolver.verbose =true;
            //C.NonLinearSolver.PrecondSolver.SolverCode = LinearSolverCode.;

            switch(SolverSwitch) {
                case 0:
                C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
                break;

                case 1:
                C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
                break;

                case 2:
                C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
                C.NonLinearSolver.UnderRelax = 0.2;
                break;

                default:
                throw new Exception("notimplemented");
            }

            //C.NonLinearSolver.PrecondSolver.SolverCode = LinearSolverCode.automatic;

            if(FDJOK == 1) {
                C.UseFDJ = true;
            } else {
                C.UseFDJ = false;
            }

            if(dt < 0) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
                C.dtFixed = dt;
                C.Endtime = 400.0;
                C.NoOfTimesteps = int.MaxValue;
            }

            C.savetodb = true;
            //dbPath = @"\\hpccluster\hpccluster-scratch\gutierrez\NatConvStudyTEST";
            C.DbPath = dbPath;// "D:\\bosss_db_NatConvection";
           

            C.ID = "70b79f17-1347-405d-b7b9-c584c648d9f0";
            //C.ProjectName = "NaturalConvectionSquareCavity";
            // C.ProjectDescription = "";

            bool restartOK = false;

            int factor = 8;
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Adiabatic no slip wall
                //2: Temperature fixed no slip wall

                //left Inlet
                if(Math.Abs(x - 0.5) < 1e-8)
                    return 3;

                //bottom Wall
                if(Math.Abs(y - 0.5 * factor) < 1e-8)
                    return 1;

                // right outlet
                if(Math.Abs(x + 0.5) < 1e-8)
                    return 2;

                //top Wall
                if(Math.Abs(y + 0.5 * factor) < 1e-8)
                    return 1;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.SinLinSpacing(-0.5, 0.5, 0.8, ncells + 1);
                var _yNodes = GenericBlas.SinLinSpacing(-0.5 * factor, 0.5 * factor, 0.8, ncells * factor + 1);
                //var _xNodes = GenericBlas.Linspace(-0.5, 0.5, ncells + 1);
                //var _yNodes = GenericBlas.Linspace(-0.5, 0.5, ncells + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                grd.EdgeTagNames.Add(1, "NoSlipNeumann");
                grd.EdgeTagNames.Add(2, "wall_tempfixed_left");
                grd.EdgeTagNames.Add(3, "wall_tempfixed_right");
                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            C.InitialMass = 8.0;
            C.EdgeTagsNusselt = new string[] { "wall_tempfixed_left", "wall_tempfixed_right", "NoSlipNeumann" };
            C.Rayleigh = Ra;
            C.Prandtl = 0.71;
            C.Reynolds = Sqrt(C.Rayleigh / C.Prandtl);
            C.T_ref_Sutherland = 600; ///////////////////////////////// check!!!!!!!!!!!!
            C.ThermodynamicPressureMode = ThermodynamicPressureMode.MassDetermined; // Because its a closed system, i.e. p0 = p0(time)
            C.MatParamsMode = MaterialParamsMode.Constant;

            double Th = 1.1; // Adimensional temperature of the hot wall
            double Tc = 0.9; // Adimensional temperature of the cold wall

            double Tinf = (Th + Tc) / 2;
            double eps = (Th - Tc) / (2 * Tinf);
            C.Froude = Math.Sqrt(2 * eps);

            if(!restartOK) {
                // initial values
                // ==============
                C.InitialValues_Evaluators.Add(VariableNames.VelocityX, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.VelocityY, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => (-1) * X[1] / (C.Froude * C.Froude));

                C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.ThermodynamicPressure, X => 1);
            } else {
                C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid(C.ID), -1);
            }

            // boundary conditions
            // ===================

            C.AddBoundaryValue("NoSlipNeumann", VariableNames.VelocityX, X => 0.0);
            C.AddBoundaryValue("NoSlipNeumann", VariableNames.VelocityY, X => 0.0);

            C.AddBoundaryValue("wall_tempfixed_left", VariableNames.Temperature, X => Th);
            C.AddBoundaryValue("wall_tempfixed_right", VariableNames.Temperature, X => Tc);

            return C;
        }
    }
}
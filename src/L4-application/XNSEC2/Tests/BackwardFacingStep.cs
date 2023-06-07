using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    /// </summary>
    public static partial class FullNSEControlExamples {

        /// <summary>
        /// NUnit test
        /// </summary>
        public static XNSEC_Control BackwardFacingStep_NUnit() {
            var C = BackwardFacingStep(DGp: 2, nCellsMult: 3, checkConsistency: false);
            C.NoOfMultigridLevels = 1;
            C.savetodb = false;
            C.DbPath = null;
            //C.savetodb = true;
            //C.DbPath = @"C:\Databases\BoSSS_DB";
            C.ChemicalReactionActive = false;
            //C.ImmediatePlotPeriod = 1;
            return C;
        }

        /// <summary>
        /// Channel flow
        /// </summary>
        public static XNSEC_Control BackwardFacingStep(int DGp = 2, int nCellsMult = 4, bool checkConsistency = false) {
            XNSEC_Control C = new XNSEC_Control();
            Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");
            Console.WriteLine("BackwardFacingStep");
            Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");

            // Solver configuration
            // ==============

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            C.ImmediatePlotPeriod = 1;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.DbPath = @"C:\Databases\BoSSS_DB";
            C.savetodb = true;
            C.ProjectName = "BackwardFacingStep";

            //C.physicsMode = PhysicsMode.LowMach;
            C.physicsMode = PhysicsMode.Combustion;
            C.MatParamsMode = MaterialParamsMode.Constant;
            C.rhoOne = false;
            C.AnalyticsolutionSwitch = false;
            C.EnableMassFractions = false;
            C.EnableTemperature = true;
            C.PhysicalParameters.IncludeConvection = true;

            C.NumberOfChemicalSpecies = 1;
            C.ChemicalReactionActive = false;
            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };
            C.ImmediatePlotPeriod = 1;
            C.HeatRelease = 0.0;
            // Parameters
            // ==============

            C.Reynolds = 100;
            C.Prandtl = 0.7132;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;


            C.HomotopyApproach = XNSEC_Control.HomotopyType.Automatic;

            //C.HomotopyVariable = XNSEC_Control.HomotopyVariableEnum.Reynolds;
            //C.homotopieAimedValue = Math.Sqrt(C.Reynolds);
            //C.StartingHomotopyValue = 10; // Suficiently easy to find solution
            //C.NonLinearSolver.HomotopyStepLongFail = 40;

            //C.ImmediatePlotPeriod = 1;
            // Grid declaration
            // ===============



            int Resolution = nCellsMult;
 
            // static double Hdim = 30; // mm channel height
            // static double H = 30 / Hdim;// nondim channel height
            // static double S = 15 / Hdim;  //nondim step height
            double Hdim = 10.1; // mm channel height
            double hdim = 5.2;
            double Sdim = Hdim - hdim;  //dim step height

            double LREF = 2 * hdim; // Hydraulic referemce
            double H = Hdim / LREF;// nondim channel height
            double S = Sdim / LREF;  //nondim step height
            double h = hdim / LREF;
            double L0 = S; // nondim upstream wall
            double L = 10 * S; // downstream wall
                                      //double L = 50 * S; // downstream wall

            int xNodesMultiplier = 15;
            double hc = Math.Pow(2, -Resolution + 1); // cell length
            double cells = 1 / hc;
            int Res = (int)cells;

            double expansionRatio = 1 + S / h;
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet
                //2: Adiabatic wall
                //3: Pressure outlet
                //4: Wall lower

                //Inlet oxidizer
                if (Math.Abs(x + L0) < 1e-8)
                    return 1;

                //Outlet
                if (Math.Abs(x - L) < 1e-8)
                    return 3;

                //Bottom Wall
                if (Math.Abs(y - 0) < 1e-8 && x > 1e-8)
                    return 4;

                return 2;
            };



            C.GridFunc = delegate {
                // Grid2D grd;

                double[] CutOut1Point1 = new double[2] { 0.0, 0.0 };
                double[] CutOut1Point2 = new double[2] { -L0, S };

                var CutOut1 = new BoSSS.Platform.Utils.Geom.BoundingBox(2);
                CutOut1.AddPoint(CutOut1Point1);
                CutOut1.AddPoint(CutOut1Point2);


                // left part
                var _leftNodes = GenericBlas.SinLinSpacing(-2 * L0, 0, 0.98, (2 * Res) + 1);
                var leftNodes = _leftNodes.GetSubVector(_leftNodes.Length / 2, _leftNodes.Length / 2 + 1);
                //right part
                var _rightNodes = GenericBlas.SinLinSpacing(0, 2 * L, 0.9, (xNodesMultiplier * Res) + 1);
                var rightNodes = _rightNodes.GetSubVector(0, _rightNodes.Length / 2 + 1);

                List<double> listX = new List<double>();
                listX.AddRange(leftNodes.Take(leftNodes.Count() - 1).ToList());
                listX.AddRange(rightNodes.Take(rightNodes.Count() + 1).ToList());
                var xNodes = listX.ToArray();


                //Bottom part
                var _bottomNodes = GenericBlas.SinLinSpacing(0, S, 0.8, (3 * Res) + 1);
                //upper part
                var _upperNodes = GenericBlas.SinLinSpacing(0, 2 * h, 0.9, 2 * (3 * Res) + 1);
                var upperNodes = _upperNodes.GetSubVector(0, _upperNodes.Length / 2 + 1);
                for (int i = 0; i < upperNodes.Count(); i++) {
                    upperNodes[i] = upperNodes[i] + S;
                }

                List<double> listY = new List<double>();
                listY.AddRange(_bottomNodes.Take(_bottomNodes.Count() - 1).ToList());
                listY.AddRange(upperNodes.Take(upperNodes.Count() + 1).ToList());


                var yNodes = listY.ToArray();
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, CutOuts: CutOut1);
                grd.EdgeTagNames.Add(1, "Velocity_Inlet");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");
                grd.EdgeTagNames.Add(2, "NoSlipNeumann");
                grd.EdgeTagNames.Add(4, "Wall");


                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

    

            // initial values
            // ==============
            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 3.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);

            C.InitialValues_Evaluators.Add("Phi", X => -1);
            // boundary conditions
            // ===================
            double umean = 1;
            double TRef = 273 + 10;
            double Tin = (273 + 10) / TRef;
            double Twall = (273 + 40) / TRef;


            C.AddBoundaryValue("Wall", VariableNames.Temperature + "#A", (X, t) => Twall);
            C.AddBoundaryValue("Wall", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);


            C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(0) + "#A", (X, t) => -6 * umean * (X[1] - S) * (X[1] - H) / (h * h));
            C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet", VariableNames.Temperature + "#A", (X, t) => Tin);
            C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);

            C.AddBoundaryValue("Pressure_Outlet");

            return C;
        }
    }
}
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Queries;
using ilPSP;
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
        static public XNSEC_Control ChannelFlowTest_NUnit(bool immersedBoundary) {
            var C = ChannelFlowTest(DGp: 2, nCellsMult:3 , checkConsistency: false, immersedBoundary);
            C.NoOfMultigridLevels = 1;
            C.savetodb = false;
            C.DbPath = null;
            //C.savetodb = true;
            //C.DbPath = @"C:\Databases\BoSSS_DB";
            C.ChemicalReactionActive = false;
            return C;
        }

        /// <summary>
        /// Channel flow
        /// </summary>
        static public XNSEC_Control ChannelFlowTest(int DGp = 1, int nCellsMult = 3, bool checkConsistency = false, bool immersedBoundary = false) {
            XNSEC_Control C = new XNSEC_Control();
            Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");
            Console.WriteLine("ChannelFlowTest");
            Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");

            // Solver configuration
            // ==============

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.DbPath = null;
            C.ProjectName = "ChannelFlowTest";

            //C.physicsMode = PhysicsMode.LowMach;
            C.physicsMode = PhysicsMode.Combustion;

            C.rhoOne = true;
            C.AnalyticsolutionSwitch = true;
            C.PhysicalParameters.IncludeConvection = false;

            C.NumberOfChemicalSpecies = 1;
            C.ChemicalReactionActive = false;
            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

            C.HeatRelease = 0.0;
            // Parameters
            // ==============

            C.Reynolds = 100.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;

            C.UseImmersedBoundary = immersedBoundary;
            C.ImmediatePlotPeriod = 1;
            double L = 1.0;
            
            if (immersedBoundary) {
                L = 0.75;
                Func<double[], double, double> PhiFunc2 = (X, t) => (X[1] - L);
                C.InitialValues_Evaluators_TimeDep.Add("Phi2", PhiFunc2);
                C.ThermalParameters.T_sat = 2; // boundary temperature
            }


            //C.ImmediatePlotPeriod = 1;
            // Grid declaration
            // ===============
            double h = Math.Pow(2, -nCellsMult + 1); // cell length
            double cells = 1 / h;
            int cells2 = (int)cells;

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet
                //2: Wall upper
                //3: Pressure outlet
                //4: Wall lower

                //Inlet oxidizer
                if(Math.Abs(x - 0) < 1e-8)
                    return 1;

                //upper Wall
                if(Math.Abs(y - 1) < 1e-8)
                    return 2;

                //Outlet
                if(Math.Abs(x - 5) < 1e-8)
                    return 3;

                //lower Wall
                if(Math.Abs(y + 0) < 1e-8)
                    return 2;
                else throw new ArgumentOutOfRangeException();
            };

            bool periodic = false;

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 5, cells2 *5+ 1);
                var _yNodes = GenericBlas.Linspace(0, 1, (cells2) + 1);
                //var _xNodes = GenericBlas.Linspace(0, 10, 2+ 1);
                //var _yNodes = GenericBlas.Linspace(-1, 1, 2+ 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, periodicX: periodic);

                if(periodic) {
                    grd.EdgeTagNames.Add(2, "wall");
                } else {
                    grd.EdgeTagNames.Add(1, "Velocity_Inlet");
                    grd.EdgeTagNames.Add(2, "wall");
                    grd.EdgeTagNames.Add(3, "Pressure_Outlet");
                }
                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // Analytic solutions
            // ===================
            //C.AnalyticPressure = X => -2.0 / 100.0 * X[0] + 0.2;
            //C.AnalyticVelocityX = X => -4 * (X[1] / L) * ((X[1] / L) - 1.0);
            double beta = 1.0; //??
            double h1 = 0.0;
            double h2 = L;

            C.AnalyticVelocityX = X =>  (X[1] >= h1 && X[1] <= h2) ? -beta * C.Reynolds / 2 * (X[1] * X[1] - (h1 + h2) * X[1] + h1 * h2): 0.0;

            C.AnalyticVelocityY = X => 0.0; // OK
            C.AnalyticPressure = X => (X[1] >= h1 && X[1] <= h2) ? -beta * X[0]+5: 0.0;

            // Analytical / manufactured solutions
            // ==============
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0) + "_an", C.AnalyticVelocityX);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1) + "_an", C.AnalyticVelocityY);
            C.InitialValues_Evaluators.Add(VariableNames.Pressure + "_an", C.AnalyticPressure);

            // initial values
            // ==============
            if(checkConsistency == true) {
                C.InitialValues_Evaluators.Add(VariableNames.Pressure, C.AnalyticPressure);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), C.AnalyticVelocityX);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);
            } else {
                C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.5);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 0.5);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 0.4);
            }
            C.InitialValues_Evaluators.Add("Phi", X => -1);
            // boundary conditions
            // ===================

            C.AddBoundaryValue("wall", VariableNames.Temperature + "#A", (X, t) => 1.0);
            C.AddBoundaryValue("wall", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
            
            if(!periodic) {
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(0) + "#A", C.AnalyticVelocityX);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Temperature + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);

                C.AddBoundaryValue("Pressure_Outlet");
            }

            return C;
        }


        /// <summary>
        /// NUnit test
        /// </summary>
        static public XNSEC_Control TwoPhaseChannelFlowTest_NUnit() {
            var C = ChannelFlow_WithInterface(2, 7, 0);
            //var C = TwoPhaseChannelFlowTest();
            C.NoOfMultigridLevels = 1;
            //C.savetodb = false;
            //C.DbPath = null;
            C.savetodb = true;
            C.DbPath = @"C:\Databases\BoSSS_DB";
            C.ChemicalReactionActive = false;
            //C.ImmediatePlotPeriod = 1;
            return C;
        }

        /// <summary>
        /// Channel flow
        /// </summary>
        static public XNSEC_Control TwoPhaseChannelFlowTest(int DGp = 2, int nCellsMult = 3, bool checkConsistency = false) {
            XNSEC_Control C = new XNSEC_Control();
            Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");
            Console.WriteLine("Two Phase channel flowTest ");
            Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");

            // Solver configuration
            // ==============

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            //C.dtFixed = 0.1;
            //C.NoOfTimesteps = 100;
            //C.Endtime = 10000;

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.DbPath = null;
            C.ProjectName = "ChannelFlowTest";

            //C.physicsMode = PhysicsMode.LowMach;
            C.physicsMode = PhysicsMode.Combustion;

            C.rhoOne = true;
            C.AnalyticsolutionSwitch = true;
            C.PhysicalParameters.IncludeConvection = false;
            //C.MomentumConvection_OK = false; // Only Stokes
            C.NumberOfChemicalSpecies = 1;
            C.ChemicalReactionActive = false;
            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

            C.HeatRelease = 0.0;
            // Parameters
            // ==============

            C.Reynolds = 1.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 4.0;
            C.PenaltyHeatConduction = 4.0;

            //C.ImmediatePlotPeriod = 1;
            // Grid declaration
            // ===============
            double h = Math.Pow(2, -nCellsMult + 1); // cell length
            double cells = 1 / h;
            int cells2 = (int)cells;

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet
                //2: Wall upper
                //3: Pressure outlet
                //4: Wall lower

                //Inlet oxidizer
                if(Math.Abs(x - 0) < 1e-8)
                    return 1;

                //upper Wall
                if(Math.Abs(y - 1) < 1e-8)
                    return 2;

                //Outlet
                if(Math.Abs(x - 10) < 1e-8)
                    return 1;

                //lower Wall
                if(Math.Abs(y + 1) < 1e-8)
                    return 2;
                else throw new ArgumentOutOfRangeException();
            };

            bool periodic = false;

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 10, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, cells2 + 1);
                //var _xNodes = GenericBlas.Linspace(0, 10, 2+ 1);
                //var _yNodes = GenericBlas.Linspace(-1, 1, 2+ 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, periodicX: periodic);

                if(periodic) {
                    grd.EdgeTagNames.Add(2, "wall");
                } else {
                    grd.EdgeTagNames.Add(1, "Velocity_Inlet");
                    grd.EdgeTagNames.Add(2, "wall");
                    //grd.EdgeTagNames.Add(3, "Pressure_Outlet");
                }
                grd.DefineEdgeTags(GridEdgeTagFunc);

                //var gDat = new GridData(grd);
                //var em1 = gDat.GetBoundaryEdges();
                //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);
                return grd;
            };

            // Analytic solutions
            // ===================
            C.AnalyticPressure = X => -2.0 / 100.0 * X[0] + 0.2;
            C.AnalyticVelocityX = X => 1.0 - 1 * Math.Pow(X[1], 2);
            C.AnalyticVelocityY = X => 0.0;

            // Analytical / manufactured solutions
            // ==============
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0) + "_an", C.AnalyticVelocityX);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1) + "_an", C.AnalyticVelocityY);
            C.InitialValues_Evaluators.Add(VariableNames.Pressure + "_an", C.AnalyticPressure);

            // initial values
            // ==============
            if(checkConsistency == true) {
                C.InitialValues_Evaluators.Add(VariableNames.Pressure, C.AnalyticPressure);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), C.AnalyticVelocityX);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 1.0);
            } else {
                C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.5);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 0.5);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 0.4);

                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0) + "#A", (X) => 1.0 - Math.Pow(X[1], 2));
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1) + "#A", (X) => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#A", (X) => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "#A", (X) => 1.0);

                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0) + "#B", (X) => (1.0 - 1 * Math.Pow(X[1], 2)) * 1);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1) + "#B", (X) => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#B", (X) => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "#B", (X) => 1.0);
            }
            /// <summary>
            /// the zero-level-set is identical to the x-axis
            /// </summary>
            /// <param name="time"></param>
            /// <returns></returns>
            Func<double[], double> GetPhi() {
                return delegate (double[] X) {
                    double x = X[0];
                    double y = X[1];
                    return y;
                };
            }
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.Option_LevelSetEvolution = LevelSetEvolution.Prescribed;

            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            //C.NonLinearSolver.UnderRelax = 0.2;

            C.InitialValues_Evaluators.Add("Phi", GetPhi());
            // boundary conditions
            // ===================

            C.AddBoundaryValue("wall", VariableNames.Temperature + "#A", (X, t) => 1.0);
            C.AddBoundaryValue("wall", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
            C.AddBoundaryValue("wall", VariableNames.Temperature + "#B", (X, t) => 1.0);
            C.AddBoundaryValue("wall", VariableNames.MassFraction0 + "#B", (X, t) => 1.0);

            if(!periodic) {
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(0) + "#A", (X, t) => (1.0 - Math.Pow(X[1], 2)) * 1);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Temperature + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);

                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(0) + "#B", (X, t) => (1.0 - 1 * Math.Pow(X[1], 2)) *1);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(1) + "#B", (X, t) => 0.0);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.Temperature + "#B", (X, t) => 1.0);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction0 + "#B", (X, t) => 1.0);

                //C.AddBoundaryValue("Pressure_Outlet");
            }

            return C;
        }

        /// <summary>
        /// control object for various testing
        /// </summary>
        /// <returns></returns>
        public static XNSEC_Control ChannelFlow_WithInterface(int p = 1, int kelem = 5, int wallBC = 0) {
            XNSEC_Control C = new XNSEC_Control();

            string _DbPath = null;

            int D = 2;
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            //C.ImmediatePlotPeriod = 1;
            // basic database options
            // ======================

            #region db

            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSEC/Channel";
            C.ProjectDescription = "Channel flow  testing";

            //C.ContinueOnIoError = false;
            //C.LogValues = XNSE_Control.LoggingValues.ChannelFlow;

            #endregion db

            // DG degrees
            // ==========
            C.SetDGdegree(p);

            // Physical Parameters
            // ===================

            #region physics

            C.PhysicalParameters.rho_A = 1;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 1;
            C.PhysicalParameters.mu_B = 10;
            double sigma = 0.0;
            C.PhysicalParameters.Sigma = sigma;

            //C.PhysicalParameters.beta_S = 0.05;
            //C.PhysicalParameters.theta_e = Math.PI / 2.0;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;
            C.rhoOne = true;
            C.ChemicalReactionActive = false;
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

            #endregion physics

            // grid generation
            // ===============

            #region grid

            double L = 5;
            double H = 2;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, L, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, H, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);
                //var grd = Grid2D.UnstructuredTriangleGrid(Xnodes, Ynodes);

                switch(wallBC) {
                    case 0:
                    goto default;
                    case 1:
                    grd.EdgeTagNames.Add(1, "velocity_inlet_lower");
                    grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                    break;

                    case 2:
                    grd.EdgeTagNames.Add(1, "navierslip_linear_lower");
                    grd.EdgeTagNames.Add(2, "navierslip_linear_upper");
                    break;

                    default:
                    grd.EdgeTagNames.Add(1, "wall_lower");
                    grd.EdgeTagNames.Add(2, "wall_upper");
                    break;
                }
                grd.EdgeTagNames.Add(3, "velocity_inlet_left");
                //grd.EdgeTagNames.Add(3, "pressure_outlet_left");
                grd.EdgeTagNames.Add(4, "pressure_outlet_right");

                //grd.EdgeTagNames.Add(3, "freeslip_left");
                //grd.EdgeTagNames.Add(4, "freeslip_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if(Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if(Math.Abs(X[1] - H) <= 1.0e-8)
                        et = 2;
                    if(Math.Abs(X[0]) <= 1.0e-8)
                        et = 3;
                    if(Math.Abs(X[0] - L) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion grid

            // Initial Values
            // ==============

            #region init

            Func<double[], double> PhiFunc = (X => X[0] - 2.0); // + (H/20)*Math.Cos(8 * Math.PI * X[0] / L));
            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            double U = 0.0;

            if(wallBC == 0) {
                U = 0.125; // Vel max
                //C.InitialValues_Evaluators.Add("VelocityX#A", X => (-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U);
                C.InitialValues_Evaluators.Add("VelocityX#B", X => (-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U);

                C.InitialValues_Evaluators.Add("VelocityX#A", X => (-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U);
            }

            #endregion init

            // exact solution
            // ==============

            #region exact

            //C.Phi = ((X, t) => PhiFunc(X));

            //C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //C.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => 1 - X[1] * X[1], (X, t) => 0 });
            //C.ExactSolutionVelocity.Add("B", new Func<double[], double, double>[] { (X, t) => 1 - X[1] * X[1], (X, t) => 0 });

            //C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            //C.ExactSolutionPressure.Add("A", (X, t) => 8 - 2 * X[0]);
            //C.ExactSolutionPressure.Add("B", (X, t) => 8 - 2 * X[0]);

            #endregion exact

            // boundary conditions
            // ===================

            #region BC

            switch(wallBC) {
                case 0:
                goto default;
                case 1:
                C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#A", X => U);
                C.AddBoundaryValue("velocity_inlet_lower", "VelocityX#B", X => U);
                C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#A", X => U);
                C.AddBoundaryValue("velocity_inlet_upper", "VelocityX#B", X => U);
                break;

                case 2:
                C.AddBoundaryValue("navierslip_linear_lower");
                C.AddBoundaryValue("navierslip_linear_upper");
                break;

                default:
                C.AddBoundaryValue("wall_lower");
                C.AddBoundaryValue("wall_upper");
                break;
            }

            if(wallBC == 0) {
                C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => ((-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0*Math.PI*(t/T)));
                C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", (X, t) => ((-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U)); // * Math.Sin(2.0 * Math.PI * (t / T)));
            }

            if(wallBC == 1) {
                C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", X => U);
                C.AddBoundaryValue("velocity_inlet_left", "VelocityX#B", X => U);
            }

            //C.AddBoundaryValue("velocity_inlet_left", "KineticEnergy#A", X => 1.0 * ((-4.0 * U / H.Pow2()) * (X[1] - H / 2.0).Pow2() + U).Pow2() / 2.0);
            ////C.AddBoundaryValue("velocity_inlet_left", "KineticEnergy#B", X => U.Pow2() / 2);
            ////C.AddBoundaryValue("pressure_outlet_left");
            C.AddBoundaryValue("pressure_outlet_right");

            #endregion BC

            // misc. solver options
            // ====================

            #region solver

            //C.ComputeEnergyProperties = true;
            //C.solveKineticEnergyEquation = true;
            ////C.CheckJumpConditions = true;
            //C.kinEViscousDiscretization = Solution.EnergyCommon.KineticEnergyViscousSourceTerms.laplaceKinE;
            //C.kinEPressureDiscretization = Solution.EnergyCommon.KineticEnergyPressureSourceTerms.divergence;
            //C.withDissipativePressure = false;

            //C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = 0.0;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LevelSet_ConvergenceCriterion = 1e-6;

            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;

            //C.Phi = (X,t) => ((X[0] - (center[0]+U*t)).Pow2() + (X[1] - center[1]).Pow2()).Sqrt() - radius;

            //C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
            //C.FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;
            //C.useFiltLevSetGradientForEvolution = true;
            //C.ReInitPeriod = 1;
            //C.ReInitOnRestart = true;

            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;

            //C.Option_LevelSetEvolution = LevelSetEvolution.Prescribed;
            //C.AdvancedDiscretizationOptions.FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;

            //C.Option_LevelSetEvolution = LevelSetEvolution.ExtensionVelocity;
            //C.EllipticExtVelAlgoControl.solverFactory = () => new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
            //C.EllipticExtVelAlgoControl.IsotropicViscosity = 1e-3;
            //C.fullReInit;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.XNSECommon.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.SkipSolveAndEvaluateResidual = true;

            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
            C.AMR_startUpSweeps = 1;
            //C.RefineStrategy = XNSE_Control.RefinementStrategy.constantInterface;
            //C.BaseRefinementLevel = 2;
            //C.RefinementLevel = 2;

            C.InitSignedDistance = false;
            C.adaptiveReInit = false;

            #endregion solver

            // Timestepping
            // ============

            #region time

            //C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            //C.Timestepper_BDFinit = TimeStepperInit.SingleInit;
            //C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.PhysicalParameters.IncludeConvection = false;
            //C.MomentumConvection_OK = false; // Only Stokes
            //double dt = 0.0138; // 5e-2; // 5e-2;
            //C.dtMax = dt;
            //C.dtMin = dt;
            //C.Endtime = 1000;
            //C.NoOfTimesteps = 100; // 500;
            C.saveperiod = 10;

            #endregion time

            return C;
        }

      

        /// <summary>
        /// Channel flow
        /// </summary>
        static public XNSEC_Control ChannelFlow_MixingTest(int DGp = 2, int nCellsMult =3) {
            XNSEC_Control C = new XNSEC_Control();
            Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");
            Console.WriteLine("ChannelFlowTestMixing");
            Console.WriteLine("////////////////////////////////////////////////////////////////////////////////////");

            //C.AdaptiveMeshRefinement = true;
            //C.RefinementLevel = 6;

            #region AMR

            C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnFieldGradient() { maxRefinementLevel = 4, FieldName = VariableNames.Temperature });
            C.activeAMRlevelIndicators.Add(new AMR_BasedOnPerssonSensor(VariableNames.Temperature, 4) );

            #endregion
            C.timeDerivativeEnergyp0_OK = false;
            C.UseSelfMadeTemporalOperator = true;
            //C.NonLinearSolver.MaxSolverIterations = 1;
            //C.NonLinearSolver.MinSolverIterations = 1;
            // Solver configuration
            // ==============
            //C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.dtFixed = 0.5;
            C.NoOfTimesteps = 100;
            C.Endtime = 10000;

            //C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.DbPath = @"C:\Databases\BoSSS_DB";
            C.ProjectName = "ChannelFlow_MixingTest";

            //C.physicsMode = PhysicsMode.LowMach;
            C.physicsMode = PhysicsMode.Combustion;

            C.rhoOne = false;

            C.NumberOfChemicalSpecies = 2;
            C.ChemicalReactionActive = false;

            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

            C.HeatRelease = 0.0;
            // Parameters
            // ==============

            C.Reynolds = 10.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 4.0;
            C.PenaltyHeatConduction = 4.0;

            //C.ImmediatePlotPeriod = 1;
            // Grid declaration
            // ===============
            double h = Math.Pow(2, -nCellsMult + 1); // cell length
            double cells = 1 / h;
            int cells2 = (int)cells;

            int lengthMult = 1;
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet
                //2: Wall upper
                //3: Pressure outlet
                //4: Wall lower

                //Inlet oxidizer
                if(Math.Abs(x - 0) < 1e-8)
                    return 1;

                //upper Wall
                if(Math.Abs(y - 1) < 1e-8)
                    return 2;

                //Outlet
                if(Math.Abs(x - 10 * lengthMult) < 1e-8)
                    return 4;

                //lower Wall
                if(Math.Abs(y + 1) < 1e-8)
                    return 3;
                else throw new ArgumentOutOfRangeException();
            };

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 10 * lengthMult, cells2 * lengthMult * 2 + 1);
                var _yNodes = GenericBlas.Linspace(-1, 1, (cells2 / 2) + 1);
                //var _xNodes = GenericBlas.Linspace(0, 10, 2+ 1);
                //var _yNodes = GenericBlas.Linspace(-1, 1, 2+ 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet");
                grd.EdgeTagNames.Add(2, "wall_top_hot");
                grd.EdgeTagNames.Add(3, "wall_bot_cold");
                grd.EdgeTagNames.Add(4, "Pressure_Outlet");

                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            C.MolarMasses = new double[] { 4.0, 1.0, 1.0, 1.0, 1.0 };

            C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => 1.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 0.0);
            C.InitialValues_Evaluators.Add(VariableNames.MassFraction1, X => 0.0);
            C.InitialValues_Evaluators.Add("Phi", X => -1);
            // boundary conditions
            // ===================
            C.AddBoundaryValue("wall_top_hot", VariableNames.Temperature + "#A", (X, t) => 2.0);
            C.AddBoundaryValue("wall_top_hot", VariableNames.MassFraction0 + "#A", (X, t) => 100.0); // shouldn't matter?

            C.AddBoundaryValue("wall_bot_cold", VariableNames.Temperature + "#A", (X, t) => 0.5);
            C.AddBoundaryValue("wall_bot_cold", VariableNames.MassFraction0 + "#A", (X, t) => 100.0); // shouldn't matter?

            C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(0) + "#A", (X, t) => 1.0 - Math.Pow(X[1], 2));
            C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet", VariableNames.Temperature + "#A", (X, t) => 1.0);
            C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction0 + "#A", (X, t) => (X[1] > 0) ? 1.0 : 0.0);
            C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction1 + "#A", (X, t) => (X[1] > 0) ? 0.0 : 1.0);

            C.AddBoundaryValue("Pressure_Outlet");

            //C.rhoOne = true;

            return C;
        }
    }
}
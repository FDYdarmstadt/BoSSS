using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using BoSSS.Foundation.IO;
using System.Linq;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    ///
    /// </summary>
    public static partial class FullNSEControlExamples {


        /// <summary>
        ///configaration for the simulation using the mixture fraction approach, where an infinite reaction rate is assumed.
        ///Used to find adequate starting solution for the full problem.
        /// </summary>
        /// <param name="DGp"></param>
        /// <param name="nCells"></param>
        /// <param name="dbPath"></param>
        /// <returns></returns>
        static public XNSEC_Control FS_XDG_pseudo2dCombustion(int DGp = 2, int nCells = 4, string dbPath = null) {
            var C = XDG_pseudo2dCombustion(DGp, nCells, dbPath, true); // true for mixture fraction calculation
            C.physicsMode = PhysicsMode.MixtureFraction;
            C.ProjectName = "XDGPseudo1DCombustion";
            string name = C.ProjectName + "P" + DGp + "K" + nCells;
            C.SessionName = "FS_" + name;

            //C.SkipSolveAndEvaluateResidual = true;

            Dictionary<string, Tuple<double, double>> Bounds = new Dictionary<string, Tuple<double, double>>();
            double eps = 0.01;
            //Mixture fraction should by definition be between 0 and 1
            Bounds.Add(VariableNames.MixtureFraction, new Tuple<double, double>(0 - eps, 1 + eps));
            C.VariableBounds = Bounds;

            C.AdaptiveMeshRefinement = true; // AMR at each time step
            if (C.AdaptiveMeshRefinement) {
                C.activeAMRlevelIndicators.Add(new AMR_onFlameSheet(C.zSt, 4));
                C.AMR_startUpSweeps = 4;
            }
                
            C.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Steady; // steady calculation, time steps used for AMR
            C.NoOfTimesteps = 2;
            C.savetodb = true;
            return C;
        }

        /// <summary>
        /// Calculation of the fully coupled problem using the mixture fraction approach.
        /// In order to "start" the combustion, a proper initial condition has to be used.
        /// The Burke-Schumann solution is used for this.
        /// </summary>
        /// <param name="DGp"></param>
        /// <param name="nCells"></param>
        /// <param name="dbPath"></param>
        /// <returns></returns>
        static public XNSEC_Control Full_XDG_pseudo2dCombustion(int DGp = 2, int nCells = 5, string dbPath = null) {
  
            var C = XDG_pseudo2dCombustion(DGp, nCells, dbPath, false); // true for MF calculation, false for combustion
            C.physicsMode = PhysicsMode.Combustion;
            C.ProjectName = "XDGPseudo1DCombustion";
            string jobName = C.ProjectName + "P" + DGp + "K" + nCells;
            C.SessionName = "Full_" + jobName;

            // limiting of variable values
            Dictionary<string, Tuple<double, double>> Bounds = new Dictionary<string, Tuple<double, double>>();
            double eps = 0.05;
            Bounds.Add(VariableNames.Temperature, new Tuple<double, double>(1.0 - eps, 10)); // Min temp should be the inlet temperature.
            Bounds.Add(VariableNames.MassFraction0, new Tuple<double, double>(0.0 - eps, 1.0 + eps)); // Between 0 and 1 per definition
            Bounds.Add(VariableNames.MassFraction1, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            Bounds.Add(VariableNames.MassFraction2, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            Bounds.Add(VariableNames.MassFraction3, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            C.VariableBounds = Bounds;

            C.NoOfTimesteps =  1; // The steady solution will be calculated again and do AMR
            C.myThermalWallType = SIPDiffusionTemperature.ThermalWallType.Adiabatic;// to do

            //C.UseMixtureFractionsForCombustionInitialization = true;

            // Select the database
            DatabaseInfo dbi = DatabaseInfo.Open(C.DbPath);

            string RestartSessionName = ("FS_" + jobName);
            var sess = dbi.Sessions.Where(s => s.Name.Equals(RestartSessionName)).ToList();   //find the session where the restart should be done from
            if (sess.Count == 0) {
                Console.WriteLine("===========================");
                Console.WriteLine("No session found for restart. The pure mixing case will be calculated.");
                Console.WriteLine("===========================");
                C.ChemicalReactionActive = false;// true for combustion, false for MF
            } else {
                C.SetRestart(sess[0]);
                Console.WriteLine("===========================");
                Console.WriteLine("A session was found for restart in the database. The case with combustion will be calculated");
                if (sess.Count > 1) {
                    Console.WriteLine("Warning: multiple jobs with the same name defined for restart. Using the most recent one");
                }
                C.ChemicalReactionActive = true;// true for combustion, false for MF
                Console.WriteLine("===========================");
            }
            return C;
        }


        /// <summary>
        /// XDG- Mixture Fraction pseudo 2d combustion (no immersed boundary)
        /// </summary>
        /// <param name="DGp"></param>
        /// <param name="nCellsMult"></param>
        /// <param name="dbpath"></param>
        /// <param name="MF"></param>
        /// <returns></returns>
        public static XNSEC_Control XDG_pseudo2dCombustion(int DGp = 2, int nCellsMult = 6, string dbpath = null, bool MF = false) {
            // --control cs:BoSSS.Application.XNSEC.FullNSEControlExamples.XDG_pseudo2dCombustion()
            XNSEC_Control C;
            if (MF) { //It it is MF calculation, use the adequate control file XNSEC_MixtureFractions.cs
                C = new XNSEC_MF_Control();
            } else {
                C = new XNSEC_Control(); // Use control file XNSEC.cs
            }

            // ==============
            // Solver configuration
            // ==============
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady; // steady calculation
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton; // nonlinear solver used
            C.NonLinearSolver.verbose = true;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.savetodb = dbpath == null ? false : true;
            C.DbPath = dbpath;

            //switch for different ters in equations
            C.EnableTemperature = true; // if true, calculation of temperature equation
            C.EnableMassFractions = true; // if true, calculation of mass fraction equation
            C.PhysicalParameters.IncludeConvection = true; // true for convective term in momentum equation
            C.ThermalParameters.IncludeConvection = true; // true for convective term in heat equation

            //Miscellaneous physical switches
            bool m_RecoilPressure = false;
            C.IncludeRecoilPressure = m_RecoilPressure;
            C.ChemicalReactionActive = false;
            C.NumberOfChemicalSpecies = C.EnableMassFractions ? 2 : 1;
            C.rhoOne = false; // if true, constant density will be considered in each phase (mainly for debugging)
 

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 4, levelSet = 0 });// refinement near interface
            //C.activeAMRlevelIndicators.Add(new BoSSS.Application.XNSEC.AMR_onFlameSheet(C.zSt, 2));
            C.SkipSolveAndEvaluateResidual = false;
            C.AgglomerationThreshold = 0.1;
            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };

            C.ImmediatePlotPeriod = 1;
            // ==============
            // Parameters
            // ==============
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;
            C.Reynolds = 1.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.Damk = 01e8*0 +1e4;
            C.ReactionRateConstants = new double[] { C.Damk, 15, 1, 1 };
            C.HeatRelease = 4.0; // used in last term of heat equation

            C.smoothingFactor = 10 * 0;
            double y_interface = 7.14;
            double[] FuelInletMassFractions = new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 };
            double[] OxidizerInletMassFractions = new double[] { 0.0, 1.0, 0.0, 0.0, 0.0 };
            C.StoichiometricCoefficients = new double[] { -1, -1, 1, 1, 0 };
            C.YFuelInlet = FuelInletMassFractions[0];
            C.YOxInlet = OxidizerInletMassFractions[1];
            C.s =1;
            C.phi = C.s * C.YFuelInlet / C.YOxInlet;
            C.zSt = 1.0 / (1.0 + C.phi);

            double prescribedMass = 1e-2; // used when mass flux is fixed
            C.prescribedMassflux_Evaluator = (X, t) => (prescribedMass); // deactivate for changing muss flux
            Console.WriteLine("The flamesheet is located at points with Z = "+ C.zSt);

            C.PlotNewtonIterations = false;
            C.ThermalParameters.T_sat = 1.0; // boundary temperature at the interface of the droplet

            C.NonLinearSolver.MaxSolverIterations = 20;
            C.NonLinearSolver.ConvergenceCriterion = 1e-7;
            C.ThermalParameters.hVap = 1;
            C.PhysicalParameters.rho_A = 2.0;
            C.PhysicalParameters.rho_B = 1;
            C.PhysicalParameters.mu_A = 0.7;
            C.PhysicalParameters.mu_B = 0.2;

            C.ThermalParameters.rho_A = 2.0; // used two different densities but physically it's the same
            C.ThermalParameters.rho_B = 1;
            C.ThermalParameters.k_A = 0.7;
            C.ThermalParameters.k_B = 0.2;

            // ===============
            // Grid declaration
            // ===============
            double h = Math.Pow(2, -nCellsMult + 1); // cell length
            double cells = 1 / h;
            int cells2 = (int)cells;
            //int cells2 = (int)nCellsMult;
            //bool reallyPseudo2D = true;
            double LLL = 30; // length of computational domain
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];
                //upper Wall
                if (Math.Abs(y - LLL) < 1e-8)
                    return 2;
                //lower Wall
                if (Math.Abs(y + 0) < 1e-8)
                    return 1;
                else throw new ArgumentOutOfRangeException();
            };

            // switch for boundary conditions
            bool m_BotPressureOutlet = false;
            bool m_TopPressureOutlet = true;

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(0, 5, 4); // 4 cells at least are needed
                var _yNodes = GenericBlas.Linspace(0, LLL, cells2 + 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, periodicX: true);

                if (m_BotPressureOutlet) {
                    grd.EdgeTagNames.Add(1, "ScalarDirichlet_PressureOutlet_bot");
                } else {
                    grd.EdgeTagNames.Add(1, "velocity_inlet_bot");
                }

                if (m_TopPressureOutlet) {
                    grd.EdgeTagNames.Add(2, "ScalarDirichlet_PressureOutlet_top");
                } else {
                    grd.EdgeTagNames.Add(2, "velocity_inlet_top");
                }

                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            //values for pressure
            Func<double[], double, double> GetPress(string species) {
                if (species == "A") {
                    return ((_2D)((x, y) => C.IncludeRecoilPressure ? -prescribedMass * prescribedMass * (1 / C.PhysicalParameters.rho_A - 1 / C.PhysicalParameters.rho_B) : 0.0)).Convert_xy2X().Convert_X2Xt();
                } else if (species == "B") {
                    return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }


            //values for velocity
            Func<double[], double, double> GetU(string species, int d) {
                if (d == 0) {
                    return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                } else if (d == 1) {
                    if (species == "A") {
                        return ((_2D)((x, y) => prescribedMass / C.PhysicalParameters.rho_A)).Convert_xy2X().Convert_X2Xt();
                    } else if (species == "B") {
                        return ((_2D)((x, y) => prescribedMass / C.PhysicalParameters.rho_B)).Convert_xy2X().Convert_X2Xt();
                    } else {
                        throw new ArgumentOutOfRangeException();
                    }
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }

            //values for phi
            Func<double[], double, double> GetPhi() {
                return ((_2D)(delegate (double x, double y) {
                    return (y - y_interface);
                })).Convert_xy2X().Convert_X2Xt();
            }

            // values for temperature
            Func<double[], double, double> GetTemperature(string species) {

                if (species == "A") {
                    return ((_2D)((x, y) => 1)).Convert_xy2X().Convert_X2Xt();
                } else if (species == "B") {
                    Func<double[], double, double> T = delegate (double[] X, double t) {              
                        double T = 5;
                        return T;
                    };

                    return T;
                } else {
                    throw new ArgumentOutOfRangeException();
                }

            };

            // =================================
            // initial values and exact solution
            // =================================
            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionTemperature = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionMassFractions = new Dictionary<string, Func<double[], double, double>[]>();
            int D = 2;

            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionPressure.Add(spc, GetPress(spc));
                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => GetU(spc, d)));
                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, GetU(spc, d).Convert_Xt2X(0.0));
                }

                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, GetPress(spc).Convert_Xt2X(0.0));
                C.InitialValues_Evaluators.Add(VariableNames.MixtureFraction + "#" + spc, X => 0.5);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, GetTemperature(spc).Convert_Xt2X(0.0));
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "#" + spc, X => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction1 + "#" + spc, X => 0.0);
            }


            //C.Phi = GetPhi();
            C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCG, GetPhi());

            // ===================
            // boundary conditions
            // =================== 

            if (!m_BotPressureOutlet) { // bc for inner boundary
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.Velocity_d(0) + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.Velocity_d(1) + "#A", (X, t) => prescribedMass / C.PhysicalParameters.rho_A);
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.MixtureFraction + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.Temperature + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.MassFraction1 + "#A", (X, t) => 0.0);
            } else {
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.Pressure + "#A", (X, t) => m_RecoilPressure ? -prescribedMass * prescribedMass * (1 / C.PhysicalParameters.rho_A - 1 / C.PhysicalParameters.rho_B) : 0.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.MixtureFraction + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.Temperature + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.MassFraction1 + "#A", (X, t) => 0.0);
            }
            if (!m_TopPressureOutlet) { // bc for outer boundary
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.Velocity_d(0) + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.Velocity_d(1) + "#A", (X, t) => prescribedMass / C.PhysicalParameters.rho_B);
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.MixtureFraction + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.Temperature + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.MassFraction0 + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.MassFraction1 + "#A", (X, t) => 1.0);
            } else {
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.Pressure + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.MixtureFraction + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.Temperature + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.MassFraction0 + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.MassFraction1 + "#A", (X, t) => 1.0);
            }
            return C;
        }

         

    }








}
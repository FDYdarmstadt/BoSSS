using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

//using BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases;
//using log4net.Filter;

namespace BoSSS.Application.XNSEC {

    public static partial class FullNSEControlExamples {

        /// <summary>
        /// Calculation of the fully coupled problem.
        /// In order to "start" the combustion, a proper initial condition has to be used.
        /// The Burke-Schumann solution is used for this.
        /// </summary>
        /// <param name="DGp"></param>
        /// <param name="nCells"></param>
        /// <param name="Refinements"></param>
        /// <param name="velMultiplier"></param>
        /// <param name="NewtonConvergenceCriterion"></param>
        /// <param name="dbPath"></param>
        /// <returns></returns>
        public static XNSEC_Control Full_DropletFlame(int DGp = 3, int nCells = 10, string dbPath = @"C:\Databases\CounterFlowFlame_StrainSweep_DilutedFuel") {
            var C = DropletCombustionHardcoded(DGp, nCells, false, dbPath);
            C.ProjectName = "DropletFlame";
            string jobName = C.ProjectName + "P" + DGp + "K" + nCells;
            C.SessionName = "Full_" + jobName;

            C.physicsMode = PhysicsMode.Combustion;
            C.VariableOneStepParameters = false;
            C.UseSelfMadeTemporalOperator = false;
            C.HeatCapacityMode = MaterialLaw_MultipleSpecies.CpCalculationMode.constant;

            //// limiting of variable values
            //Dictionary<string, Tuple<double, double>> Bounds = new Dictionary<string, Tuple<double, double>>();
            //double eps = 0.05;
            //Bounds.Add(VariableNames.Temperature, new Tuple<double, double>(1.0 - eps, 10)); // Min temp should be the inlet temperature.
            //Bounds.Add(VariableNames.MassFraction0, new Tuple<double, double>(0.0 - eps, 1.0 + eps)); // Between 0 and 1 per definition
            //Bounds.Add(VariableNames.MassFraction1, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            //Bounds.Add(VariableNames.MassFraction2, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            //Bounds.Add(VariableNames.MassFraction3, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            //C.VariableBounds = Bounds;

            C.AMR_startUpSweeps = 1;
            C.NoOfTimesteps = 5; // The steady solution will be calculated again and do AMR

            C.UseMixtureFractionsForCombustionInitialization = true;

            C.NoOfTimesteps = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MinSolverIterations = 10;

            // Select the database
            DatabaseInfo dbi = DatabaseInfo.Open(C.DbPath);

            string RestartSessionName = ("FS_" + jobName);
            var sess = dbi.Sessions.Where(s => s.Name.Equals(RestartSessionName)).ToList();   //find the session where the restart should be done from
            if (sess.Count == 0) {
                Console.WriteLine("===========================");
                Console.WriteLine("No session found for restart. The pure mixing case will be calculated.");
                Console.WriteLine("===========================");
                C.ChemicalReactionActive = false;
            } else {
                C.SetRestart(sess[0]);
                Console.WriteLine("===========================");
                Console.WriteLine("A session was found for restart in the database. The case with combustion will be calculated");
                if (sess.Count > 1) {
                    Console.WriteLine("Warning: multiple jobs with the same name defined for restart. Using the most recent one");
                }
                Console.WriteLine("===========================");
            }
            return C;
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="DGp"></param>
        /// <param name="nCells"></param>
        /// <param name="numOfRefinements"></param>
        /// <param name="NewtonConvergenceCriterion"></param>
        /// <returns></returns>
        public static XNSEC_MF_Control FS_DropletFlame(int DGp = 2, int nCells = 16, string dbPath = @"C:\Databases\TestCDFlame") {
            var C = (XNSEC_MF_Control)DropletCombustionHardcoded(DGp, nCells, true, dbPath);
            C.physicsMode = PhysicsMode.MixtureFraction;
            C.ProjectName = "DropletFlame";
            string name = C.ProjectName + "P" + DGp + "K" + nCells;
            C.SessionName = "FS_" + name;

            C.UseSelfMadeTemporalOperator = false;

            Dictionary<string, Tuple<double, double>> Bounds = new Dictionary<string, Tuple<double, double>>();
            double eps = 0.05;
            Bounds.Add(VariableNames.MixtureFraction, new Tuple<double, double>(0 - eps, 1 + eps));
            C.VariableBounds = Bounds;

            C.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Steady;
            C.NoOfTimesteps = 5;
            C.savetodb = true;
            return C;
        }

        /// <summary>
        /// XDG- Mixture Fraction pseudo 2d combustion (no immersed boundary)
        /// </summary>
        /// <param name="DGp"></param>
        /// <param name="nCellsMult"></param>
        /// <param name="MF"></param>
        /// <returns></returns>
        public static XNSEC_Control XDG_DropletCombustion(int DGp = 1, int nCellsMult = 5, bool MF = true) {
            XNSEC_Control C = new XNSEC_Control();

            // Solver configuration
            // ==============
            C.ImmediatePlotPeriod = 1;

            C.physicsMode = PhysicsMode.Combustion;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.DbPath = null;
            bool m_RecoilPressure = true;
            C.IncludeRecoilPressure = m_RecoilPressure;
            C.savetodb = true;
            C.DbPath = @"C:\Databases\BoSSS_DB\";
            C.ChemicalReactionActive = false;
            C.EnableTemperature = false;
            C.EnableMassFractions = false;
            C.NumberOfChemicalSpecies = C.EnableMassFractions ? 2 : 1;

            C.rhoOne = true;
            C.PhysicalParameters.IncludeConvection = true;
            C.ThermalParameters.IncludeConvection = true;

            C.SkipSolveAndEvaluateResidual = false;
            //C.AgglomerationThreshold = 0.3;
            C.ChemicalReactionActive = false;
            C.SetDGdegree(DGp);
            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 };
            C.HeatRelease = 0.0;
            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.NonLinearSolver.MaxSolverIterations = 10;
            // Parameters
            // ==============
            //C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            C.Reynolds = 1.0;
            C.Prandtl = 1.0;
            C.Schmidt = 1.0;
            C.PenaltyViscMomentum = 1.0;
            C.PenaltyHeatConduction = 1.0;

            C.ImmediatePlotPeriod = 1;
            C.AdaptiveMeshRefinement = false;
            C.AMR_startUpSweeps = 2;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1, levelSet = 0 });
            //C.activeAMRlevelIndicators.Add(new BoSSS.Application.XNSEC.AMR_onFlameSheet(C.zSt, 2));

            double prescribedMass = 1.0;
            C.prescribedMassflux_Evaluator = (X, t) => (prescribedMass);

            C.ThermalParameters.T_sat = 1.0; // boundary temperature

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            C.ThermalParameters.hVap = 1;
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.rho_B = 0.25;
            C.PhysicalParameters.mu_A = 10;
            C.PhysicalParameters.mu_B = 1;
            C.ThermalParameters.rho_A = 1.0;
            C.ThermalParameters.rho_B = 0.25;
            // Grid declaration
            // ===============
            double h = Math.Pow(2, -nCellsMult + 1); // cell length
            double cells = 1 / h;
            int cells2 = (int)cells;
            //C.UseMixtureFractionsForCombustionInitialization = true;
            double LL = 50;
            double R = 10; //droplet radius

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                if (Math.Abs(x - LL) < 1e-8 | Math.Abs(y - LL) < 1e-8 || Math.Abs(x + LL) < 1e-8 || Math.Abs(y + LL) < 1e-8) {
                    return 2;
                } else
                    return 1;
            };

            bool m_InsideBC_PressureOutlet = false;
            bool m_OutsideBC_PressureOutlet = true;

            C.GridFunc = delegate {
                double co = 1; // bounding box, has to be smaller than droplet radius

                int n = 10;
                int nbb = 1;

                IEnumerable<double> xs_R__ = GenericBlas.SinLinSpacing(co, co + 2 * (LL - co), 0.8, 2 * n + 1).ToList().Take(n + 1);
                double[] xs_R = xs_R__.Skip(1).ToArray(); // this

                double[] xs_M2 = GenericBlas.Linspace(-co, co, 2 * nbb + 1).ToList().Skip(1).ToArray(); // this

                var xs_L_ = xs_R__.Reverse().ToList();
                xs_L_.ScaleV(-1);
                var xs_L = xs_L_.ToArray();

                double[] nodes = xs_L.Concat(xs_M2).Concat(xs_R).ToArray();     // All X Nodes

                double[] CutOut1Point1 = new double[2] { -co, -co };
                double[] CutOut1Point2 = new double[2] { co, co };

                var CutOut1 = new BoSSS.Platform.Utils.Geom.BoundingBox(2);
                CutOut1.AddPoint(CutOut1Point1);
                CutOut1.AddPoint(CutOut1Point2);
                var grd = Grid2D.Cartesian2DGrid(nodes, nodes, CutOuts: CutOut1);

                if (m_InsideBC_PressureOutlet) {
                    grd.EdgeTagNames.Add(1, "ScalarDirichlet_PressureOutlet_bot");
                } else {
                    grd.EdgeTagNames.Add(1, "velocity_inlet_bot");
                }

                if (m_OutsideBC_PressureOutlet) {
                    grd.EdgeTagNames.Add(2, "ScalarDirichlet_PressureOutlet_top");
                } else {
                    grd.EdgeTagNames.Add(2, "velocity_inlet_top");
                }

                grd.DefineEdgeTags(GridEdgeTagFunc);

                //var gDat = new GridData(grd);
                //var em1 = gDat.GetBoundaryEdges();
                //em1.SaveToTextFile("alledges2.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                return grd;
            };

            Func<double[], double, double> GetPress(string species) {
                if (species == "A") {
                    return ((_2D)((x, y) => m_RecoilPressure ? -prescribedMass * prescribedMass * (1 / C.PhysicalParameters.rho_A - 1 / C.PhysicalParameters.rho_B) : 0.0)).Convert_xy2X().Convert_X2Xt();
                } else if (species == "B") {
                    return ((_2D)((x, y) => 0.0)).Convert_xy2X().Convert_X2Xt();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }

            Func<double[], double, double> GetU(string species, int d) {
                if (d == 0) {
                    if (species == "A") {
                        return ((_2D)(
                             delegate (double x, double y) {
                                 return x / (Math.Sqrt(x * x + y * y)) *prescribedMass/ C.PhysicalParameters.rho_A;
                             })
                            ).Convert_xy2X().Convert_X2Xt();
                    } else if (species == "B") {
                        return ((_2D)(
                  delegate (double x, double y) {
                      return x / (Math.Sqrt(x * x + y * y)) * prescribedMass / C.PhysicalParameters.rho_B;
                  })
                 ).Convert_xy2X().Convert_X2Xt();
                    } else {
                        throw new ArgumentOutOfRangeException();
                    }
                } else if (d == 1) {
                    if (species == "A") {
                        return ((_2D)(
                 delegate (double x, double y) {
                     return y / (Math.Sqrt(x * x + y * y)) * prescribedMass / C.PhysicalParameters.rho_A;
                 })
                ).Convert_xy2X().Convert_X2Xt();
                    } else if (species == "B") {
                        return ((_2D)(
                 delegate (double x, double y) {
                     return y / (Math.Sqrt(x * x + y * y)) * prescribedMass / C.PhysicalParameters.rho_B;
                 })
                ).Convert_xy2X().Convert_X2Xt();
                    } else {
                        throw new ArgumentOutOfRangeException();
                    }
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            }

            Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()).Sqrt() - R);         // signed distance
            //Func<double[], double> PhiFunc = (X => ((X[0] - 0.0).Pow2() + (X[1] - 0.0).Pow2()) - r.Pow2());         // quadratic

            Func<double[], double, double> GetPhi() {
                return ((_2D)(delegate (double x, double y) {
                    return ((x - 0.0).Pow2() + (y - 0.0).Pow2()).Sqrt() - R;
                })).Convert_xy2X().Convert_X2Xt();
            }

            Func<double[], double, double> GetTemperature(string species) {
                if (species == "A") {
                    return ((_2D)((x, y) => 1)).Convert_xy2X().Convert_X2Xt();
                } else if (species == "B") {
                    return ((_2D)((x, y) => 1)).Convert_xy2X().Convert_X2Xt();
                } else {
                    throw new ArgumentOutOfRangeException();
                }
            };

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
                //C.ExactSolutionMassFractions.Add(spc, NoChemSpc.ForLoop(q => tst.GetMassFractions(spc, q)));

                //C.ExactSolutionTemperature.Add(spc, tst.GetTemperature(spc));
                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, GetU(spc, d).Convert_Xt2X(0.0));
                }

                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, GetPress(spc).Convert_Xt2X(0.0));
                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, GetTemperature(spc).Convert_Xt2X(0.0));
                //C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, X => 1.0) ;
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "#" + spc, X => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction1 + "#" + spc, X => 0.5);
            }

            C.Phi = GetPhi();
            C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCG, GetPhi());

            // boundary conditions
            // ===================

            if (!m_InsideBC_PressureOutlet) {
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.Velocity_d(0) + "#A", GetU("A",0).Convert_Xt2X(0.0));
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.Velocity_d(1) + "#A", GetU("A",1).Convert_Xt2X(0.0));
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.Temperature + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("velocity_inlet_bot", VariableNames.MassFraction1 + "#A", (X, t) => 0.0);
            } else {
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.Pressure + "#A", (X, t) => m_RecoilPressure ? -prescribedMass * prescribedMass * (1 / C.PhysicalParameters.rho_A - 1 / C.PhysicalParameters.rho_B) : 0.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.Temperature + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_bot", VariableNames.MassFraction1 + "#A", (X, t) => 0.0);
            }
            if (!m_OutsideBC_PressureOutlet) {
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.Velocity_d(0) + "#A", GetU("A", 0).Convert_Xt2X(0.0));
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.Velocity_d(1) + "#A", GetU("A", 1).Convert_Xt2X(0.0));
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.Temperature + "#A", (X, t) => 4.0);
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.MassFraction0 + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("velocity_inlet_top", VariableNames.MassFraction1 + "#A", (X, t) => 0.23);
            } else {
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.Pressure + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.Temperature + "#A", (X, t) => 4.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.MassFraction0 + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("ScalarDirichlet_PressureOutlet_top", VariableNames.MassFraction1 + "#A", (X, t) => 0.5);
            }
            return C;
        }

        public static XNSEC_Control DropletCombustionHardcoded(int DGp = 2, int nCells = 16 * 2, bool MF = false, string _dbPath = @"C:\Databases\BoSSS_DB_COMBUSTION") {
            string dbPath = _dbPath;

            // Problem Definition
            //===================
            double TemperatureInFuel = 300; // Tsat?
            double TemperatureInOxidizer = 300;
            double R = 0.5; // droplet length
            double AtmPressure = 101325; // Pa
            var CC = new ChemicalConstants();

            double[] MWs = new double[] { CC.MW_CH4, CC.MW_O2, CC.MW_CO2, CC.MW_H2O, CC.MW_N2 };
            double[] FuelInletMassFractions = new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 };
            double[] OxidizerInletMassFractions = new double[] { 0.0, 0.23, 0.0, 0.0, 0.77 };

            double mwFuel = CC.getAvgMW(MWs, FuelInletMassFractions);
            double mwAir = CC.getAvgMW(MWs, OxidizerInletMassFractions);

            double uInFuel = 0;
            double uInAir = 0;

            bool restartOK = false;

            // Reference values
            //===================
            // Basic units to be used: Kg, m, s, mol, pa,
            //double LRef = (Dth / uRef); // Length based on heat diffusion.
            double TRef = TemperatureInOxidizer;
            double pRef = AtmPressure; // Pa
            double uRef = Math.Max(uInFuel, uInAir); //
            double LRef = R; // droplet length
            bool useAdimensional = true;
            var C = MF ? new XNSEC_MF_Control() : new XNSEC_Control();

            //C = new XNSEC_Control(DGp, pRef, uRef, TRef, useAdimensional, LRef, FuelInletMassFractions, OxidizerInletMassFractions);

            C.rhoOne = true;
            C.physicsMode = PhysicsMode.Combustion;
            C.PhysicalParameters.IncludeConvection = false;
            C.Reynolds = 10;
            C.Damk = 1e10;
            C.HeatRelease = 5;
            C.Prandtl = 1;
            C.Lewis = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };
            C.MatParamsMode = MaterialParamsMode.Constant;

            C.SetTimeSteppingOptions(-1, 100);
            C.UseImmersedBoundary = true;

            C.AdaptiveMeshRefinement = true;
            C.AMR_startUpSweeps = 4;
            C.zSt = 0.4;

            C.activeAMRlevelIndicators.Add(new BoSSS.Application.XNSEC.AMR_onFlameSheet(C.zSt, 6));
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 4, levelSet = 1 });

            Func<double[], double, double> PhiFunc2 = (X, t) => -(Math.Sqrt(X[0] * X[0] + X[1] * X[1]) - R);
            C.InitialValues_Evaluators_TimeDep.Add("Phi2", PhiFunc2);
            C.InitialValues_Evaluators.Add("Phi", X => -1);

            C.ThermalParameters.T_sat = 1.3; // boundary temperature
            C.SetSaveOptions(dbPath, savePeriod: 1);

            // Grid declaration
            // ===============
            double h = Math.Pow(2, -4 + 1); // cell length
            double cells = 1 / h;
            int cells2 = (int)cells;
            //C.UseMixtureFractionsForCombustionInitialization = true;
            double LL = 20;

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet
                //2: Wall upper
                //3: Pressure outlet
                //4: Wall lower

                if (Math.Abs(x - LL) < 1e-8)
                    return 1;

                if (Math.Abs(y - LL) < 1e-8)
                    return 1;

                if (Math.Abs(x + LL) < 1e-8)
                    return 1;

                if (Math.Abs(y + LL) < 1e-8)
                    return 1;
                else throw new ArgumentOutOfRangeException();
            };

            bool periodic = true;

            C.GridFunc = delegate {
                var _xNodes = GenericBlas.Linspace(-LL, LL, cells2 + 1);
                var _yNodes = GenericBlas.Linspace(-LL, LL, (cells2) + 1);

                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, periodicX: false);

                grd.EdgeTagNames.Add(1, "Velocity_Inlet");

                grd.DefineEdgeTags(GridEdgeTagFunc);
                return grd;
            };

            // ==============
            if (!restartOK) {
                C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => TemperatureInOxidizer / C.TRef);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction1, X => C.YOxInlet);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction2, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction3, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction4, X => 1.0 - C.YOxInlet);
                C.InitialValues_Evaluators.Add(VariableNames.MixtureFraction, X => 1.0);

                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0);
            } else {
                C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid(C.ID), new TimestepNumber(-1));
                C.AMR_startUpSweeps = 0;
            }

            // boundary conditions
            // ===================

            C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(0) + "#A", (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet", VariableNames.Velocity_d(1) + "#A", (X, t) => 0.0);

            if (MF) {
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.MixtureFraction + "#A", (X, t) => 0.0);
            } else {
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction0 + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction1 + "#A", (X, t) => C.YOxInlet);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction2 + "#A", (X, t) => 0.0);
                C.AddBoundaryValue("Velocity_Inlet", VariableNames.MassFraction3 + "#A", (X, t) => 0.0);
            }

            return C;
        }
    }
}
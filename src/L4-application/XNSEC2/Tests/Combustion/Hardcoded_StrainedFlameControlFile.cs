using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

//using BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases;
//using log4net.Filter;

namespace BoSSS.Application.XNSEC {

    static public partial class FullNSEControlExamples {

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
        static public XNSEC_Control Full_CounterDiffusionFlame(int DGp = 3, int nCells = 10, double velMultiplier = 1, string dbPath = @"C:\Databases\CounterFlowFlame_StrainSweep_DilutedFuel") {
            double desiredVelMultiplier = velMultiplier;

            var C = LaminarStrainedFlameHardcoded(DGp, nCells, 1.0, dbPath);
            C.ProjectName = "CounterDifFlame";
            string jobName = C.ProjectName + "P" + DGp + "K" + nCells;
            C.SessionName = "Full_" + jobName;
            
            //C.cpOne = true;
            C.physicsMode = PhysicsMode.Combustion;
            C.VariableOneStepParameters = false;
            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.None;
            C.UseSelfMadeTemporalOperator = false;
            C.HeatCapacityMode = MaterialLaw_MultipleSpecies.CpCalculationMode.onlyN2;
            //C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new AMR_onFlameSheet(C.zSt, 2));
            C.activeAMRlevelIndicators.Add(new AMR_BasedOnVariableLimits("Temperature",new double[] { -100,5},3));

            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnFieldGradient(3, 0.85, VariableNames.Temperature));
            //C.activeAMRlevelIndicators.Add(new AMR_onReactiveZones(C.MolarMasses, 3, 0.8));
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnPerssonSensor(VariableNames.Temperature, 3));
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnPerssonSensor(VariableNames.MassFraction0, 3));
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnPerssonSensor(VariableNames.MassFraction2, 3));
            //C.activeAMRlevelIndicators.Add(new BoSSS.Application.XNSEC.AMR_onProblematicPoints(troubledPoints, C.AMR_startUpSweeps));

            // limiting of variable values
            Dictionary<string, Tuple<double, double>> Bounds = new Dictionary<string, Tuple<double, double>>();
            double eps = 0.05;
            Bounds.Add(VariableNames.Temperature, new Tuple<double, double>(1.0 - eps, 10)); // Min temp should be the inlet temperature.
            Bounds.Add(VariableNames.MassFraction0, new Tuple<double, double>(0.0 - eps, 1.0 + eps)); // Between 0 and 1 per definition
            Bounds.Add(VariableNames.MassFraction1, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            Bounds.Add(VariableNames.MassFraction2, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            Bounds.Add(VariableNames.MassFraction3, new Tuple<double, double>(0.0 - eps, 1.0 + eps));
            C.VariableBounds = Bounds;

            C.AMR_startUpSweeps = 1;
            C.NoOfTimesteps = 5; // The steady solution will be calculated again and do AMR
            C.myThermalWallType = SIPDiffusionTemperature.ThermalWallType.Adiabatic;

            C.UseMixtureFractionsForCombustionInitialization = true;
            C.LinearSolver.SolverCode = Solution.Control.LinearSolverCode.exp_Kcycle_schwarz;
            C.LinearSolver.NoOfMultigridLevels = 10;
            C.LinearSolver.verbose = true;

            C.NoOfTimesteps = 10;
            C.NonLinearSolver.MaxSolverIterations = 10;
            C.NonLinearSolver.MinSolverIterations = 10;
            bool useHomotopy = false;
            if(useHomotopy) {
                C.HomotopyApproach = XNSEC_Control.HomotopyType.Automatic;
                C.HomotopyVariable = XNSEC_Control.HomotopyVariableEnum.VelocityInletMultiplier;
                C.homotopieAimedValue = desiredVelMultiplier;
            }

            // Select the database
            DatabaseInfo dbi = DatabaseInfo.Open(C.DbPath);

            string RestartSessionName = ("FS_" + jobName);
            var sess = dbi.Sessions.Where(s => s.Name.Equals(RestartSessionName)).ToList();   //find the session where the restart should be done from
            if(sess.Count == 0) {
                Console.WriteLine("===========================");
                Console.WriteLine("No session found for restart. The pure mixing case will be calculated.");
                Console.WriteLine("===========================");
                C.ChemicalReactionActive = false;
            } else {
                C.SetRestart(sess[0]);
                Console.WriteLine("===========================");
                Console.WriteLine("A session was found for restart in the database. The case with combustion will be calculated");
                if(sess.Count > 1) {
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
        static public XNSEC_Control FS_CounterDiffusionFlame(int DGp = 2, int nCells = 16, double velMultiplier = 1, string dbPath = @"C:\Databases\TestCDFlame") {
            var C = LaminarStrainedFlameHardcoded(DGp, nCells, velMultiplier, dbPath);
            C.physicsMode = PhysicsMode.MixtureFraction;
            C.ProjectName = "CounterDifFlame";
            string name = C.ProjectName + "P" + DGp + "K" + nCells /*+ "Ne" + NewtonConvergenceCriterion*/;
            C.SessionName = "FS_" + name;

            C.UseSelfMadeTemporalOperator = false;
            C.AdaptiveMeshRefinement = true;



            bool useHomotopy = false;
            if(useHomotopy) {
                C.HomotopyApproach = XNSEC_Control.HomotopyType.Automatic;
                C.HomotopyVariable = XNSEC_Control.HomotopyVariableEnum.VelocityInletMultiplier;
                C.homotopieAimedValue = velMultiplier;
            }



            Dictionary<string, Tuple<double, double>> Bounds = new Dictionary<string, Tuple<double, double>>();
            double eps = 0.05;
            Bounds.Add(VariableNames.MixtureFraction, new Tuple<double, double>(0- eps, 1 + eps));      
            C.VariableBounds = Bounds;

            C.NonLinearSolver.MaxSolverIterations = 500;
            C.AMR_startUpSweeps = 2;
            double radius_inlet = 0.5;
            var troubledPoints = new System.Collections.Generic.List<double[]>();
            troubledPoints.Add(new double[] { 0, +radius_inlet });
            troubledPoints.Add(new double[] { 0, -radius_inlet });
            troubledPoints.Add(new double[] { 1, +radius_inlet });
            troubledPoints.Add(new double[] { 1, -radius_inlet });
            //C.activeAMRlevelIndicators.Add(new AMR_onFlameSheet(C.zSt, 4));
            //C.activeAMRlevelIndicators.Add(new BoSSS.Application.XNSEC.AMR_onProblematicPoints(troubledPoints, C.AMR_startUpSweeps));
            C.activeAMRlevelIndicators.Add(new BoSSS.Application.XNSEC.AMR_RefineAroundProblematicPoints(troubledPoints, C.AMR_startUpSweeps, 0.01));
            //C.ImmediatePlotPeriod = 1;

            //C.activeAMRlevelIndicators.Add(new BoSSS.Application.XNSEC.AMR_AroundCenterline(C.AMR_startUpSweeps));
            C.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Steady;
            C.NoOfTimesteps = 5;
            C.savetodb = true;
            return C;
        }

        /// <summary>
        /// CounterDifussionFlameTest2 Flame
        /// </summary>
        static public XNSEC_Control LaminarStrainedFlameHardcoded(int DGp = 2, int nCells = 16 * 2, double velMultiplier = 1, string _dbPath = @"C:\Databases\BoSSS_DB_COMBUSTION") {
            //Console.WriteLine(" ===========================================");
            //Console.WriteLine("LaminarStrainedFlame");
            //Console.WriteLine("===========================================");

            // Parameters
            // ==============
            //var Conc1 = new double[] { 0.2, 0.8 }; // inlet concentrations...

            string dbPath = _dbPath;

            double L = 20.0 / 1000; // separation between the two inlets
            // Problem Definition
            //===================
            double TemperatureInFuel = 300;
            double TemperatureInOxidizer = 300;

            double massFuelIn = 0.02400 * velMultiplier; //kg/m2s
            double massAirIn = 0.02400 * 3 * velMultiplier; //kg/m2s
            double AtmPressure = 101325; // Pa
            var CC = new ChemicalConstants();

            double[] MWs = new double[] { CC.MW_CH4, CC.MW_O2, CC.MW_CO2, CC.MW_H2O, CC.MW_N2 };
            double[] FuelInletMassFractions = new double[] { 0.2, 0.0, 0.0, 0.0, 0.8 };
            double[] OxidizerInletMassFractions = new double[] { 0.0, 0.23, 0.0, 0.0, 0.77 };

            double mwFuel = CC.getAvgMW(MWs, FuelInletMassFractions);
            double mwAir = CC.getAvgMW(MWs, OxidizerInletMassFractions);
            double densityAirIn = AtmPressure * mwAir / (CC.R_gas * TemperatureInOxidizer * 1000); // Kg/m3. ok
            double densityFuelIn = AtmPressure * mwFuel / (CC.R_gas * TemperatureInFuel * 1000); // Kg/m3. ok
            double uInFuel = massFuelIn / densityFuelIn;
            double uInAir = massAirIn / densityAirIn;
            //double uInAir = Math.Sqrt(densityFuelIn / densityAirIn) * uInFuel;
            bool restartOK = false;

            // Reference values
            //===================
            // Basic units to be used: Kg, m, s, mol, pa,
            //double LRef = (Dth / uRef); // Length based on heat diffusion.
            double TRef = TemperatureInOxidizer;// Reference temperature  is the inlet temperature, (K)
            double pRef = AtmPressure; // Pa
            double uRef = Math.Max(uInFuel, uInAir); //
            double LRef = L;
            bool useAdimensional = true;

            XNSEC_Control C = new XNSEC_Control(DGp, pRef, uRef, TRef, useAdimensional, LRef, FuelInletMassFractions, OxidizerInletMassFractions);
            C.physicsMode = PhysicsMode.Combustion;
            C.PhysicalParameters.IncludeConvection = true;

            C.MatParamsMode = MaterialParamsMode.Sutherland;
            C.NonLinearSolver.MaxSolverIterations = 30;

            C.SetTimeSteppingOptions(-1, 100);

            C.Lewis = new double[] { 0.97, 1.11, 1.39, 0.83, 1.0 }; ////////////////////////////////
            //C.Lewis = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 }; ////////////////////////////////

            C.SetSaveOptions(dbPath, savePeriod: 1);

            double separation = L / (C.LRef);
            double xleft = 0;
            double xright = separation;
            double ybot = -separation * 3;
            double ytop = separation * 3;
            double R_dim = L / 2; ///////////////////////////////////////////  L / 4;

            double radius_inlet = R_dim / C.LRef;

            C.usesimplifiedGeometry = false;

            //C.ImmediatePlotPeriod = 1;
            Console.WriteLine("Velocity fuel: " + uInFuel);
            Console.WriteLine("Velocity air: " + uInAir);
            Console.WriteLine("Strain: {0}, 1/s", (Math.Abs(uInFuel) + Math.Abs(uInAir)) / L);

            //Grid declaration
            // ===============

            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];
                //Edge tags
                //1: Velocity inlet O_2
                //2: Velocity inlet CH_4
                //3: Pressure outlet
                if(Math.Abs(x - xleft) < 1e-8) { // Left boundary
                    if(Math.Abs(y) - radius_inlet < 1e-8) { // Fuel Inlet
                        return 1;
                    } else {
                        return 4;//4
                    }
                }
                if(Math.Abs(x - xright) < 1e-8) { // right boundary
                    if(Math.Abs(y) - radius_inlet < 1e-8) { // oxy Inlet
                        return 2;//2
                    } else {
                        return 4;//4
                    }
                } else {
                    return 3; //3 Pressure outlet
                }
            };

            C.GridFunc = delegate {
                double[] _xNodes;
                double[] _yNodes;

                if(C.usesimplifiedGeometry) {
                    _xNodes = GenericBlas.Linspace(xleft, xright, nCells + 1);
                    _yNodes = GenericBlas.Linspace(-radius_inlet, radius_inlet, (int)nCells / 4 + 1);
                } else {
                    // X-NODES
                    _xNodes = GenericBlas.Linspace(xleft, xright, (int)(nCells * 2.5 + 1));
                    // Y-NODES
                    bool equidisant = false;
                    if(equidisant) {
                        var _yNodes1 = GenericBlas.Linspace(ybot, -radius_inlet, (int)(nCells * 4) / 2 + 1).ToList();
                        var _yNodes2 = GenericBlas.Linspace(-radius_inlet, radius_inlet, (int)(nCells * 1) + 1).ToList();
                        var _yNodes3 = GenericBlas.Linspace(radius_inlet, ytop, (int)(nCells * 4) / 2 + 1).ToList();
                        _yNodes2.RemoveAt(_yNodes2.Count() - 1); //removes last element
                        _yNodes2.RemoveAt(0); //removes first element
                        _yNodes = ArrayTools.Cat(_yNodes1, _yNodes2, _yNodes3);
                    } else {
                        double sf3 = 0.80;
                        var _yNodes2 = GenericBlas.Linspace(-radius_inlet, radius_inlet, nCells + 1).ToList(); // center
                        List<double> yNodesTop = (GenericBlas.SinLinSpacing(radius_inlet, (ytop - radius_inlet) * 2 + radius_inlet, sf3, nCells * 4).ToList()); // Nodes corresponding to the oxidizer inlet, right part
                        var myYnodesTop = yNodesTop.GetSubVector(0, yNodesTop.Count / 2 + 1); // Take only "bottom side" of node array
                        var myYNodesBot = myYnodesTop.CloneAs();
                        myYNodesBot.ScaleV(-1.0);
                        Array.Reverse(myYNodesBot);
                        List<double> list2 = new List<double>();
                        list2.AddRange(myYNodesBot.Take(myYNodesBot.Length - 1).ToList());
                        list2.AddRange(_yNodes2.Take(_yNodes2.Count() - 1).ToList());
                        list2.AddRange(myYnodesTop.ToList());
                        _yNodes = list2.ToArray();
                    }
                }

                Console.WriteLine("Number of cells in the X direction: {0}", _xNodes.Length - 1);
                Console.WriteLine("Number of cells in the Y direction: {0}", _yNodes.Length - 1);
                var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes);
                grd.EdgeTagNames.Add(1, "Velocity_Inlet_CH4");
                grd.EdgeTagNames.Add(2, "Velocity_Inlet_O2");
                grd.EdgeTagNames.Add(3, "Pressure_Outlet");

                if(!C.usesimplifiedGeometry) {
                    grd.EdgeTagNames.Add(4, "Wall");
                }
                grd.DefineEdgeTags(GridEdgeTagFunc);

                //var gDat = new GridData(grd);
                //var em1 = gDat.GetBoundaryEdges();
                //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                return grd;
            };

            //Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
            //    double x = X[0];
            //    double y = X[1];
            //    //Edge tags
            //    //1: Velocity inlet O_2
            //    //2: Velocity inlet CH_4
            //    //3: Pressure outlet
            //    if(Math.Abs(x - separation) < 1e-8 && Math.Abs(y) - radius_inlet < 1e-8) {
            //        return 1;
            //    }
            //    if(Math.Abs(x - 2 * separation) < 1e-8 && Math.Abs(y) - radius_inlet < 1e-8) {
            //        return 2;
            //    }

            //    if(Math.Abs(y - radius_inlet) < 1e-8 && x > 0 && x < separation) {
            //        return 4;
            //    }
            //    if(Math.Abs(y - radius_inlet) < 1e-8 && x > 2 * separation && x < 3 * separation) {
            //        return 4;
            //    }
            //    if(Math.Abs(y + radius_inlet) < 1e-8 && x > 0 && x < separation) {
            //        return 4;
            //    }
            //    if(Math.Abs(y + radius_inlet) < 1e-8 && x > 2 * separation && x < 3 * separation) {
            //        return 4;
            //    }

            //    return 3;
            //};

            //C.GridFunc = delegate {
            //    double[] _xNodes;
            //    double[] _yNodes;

            //    // X-NODES
            //    _xNodes = GenericBlas.Linspace(xleft, xright * 3, (int)(nCells * 1.5 + 1));
            //    // Left and right part of x nodes
            //    double[] XNODES1 = GenericBlas.Linspace(0, separation, (int)(nCells / 3 + 1));
            //    double[] XNODES2 = GenericBlas.Linspace(separation * 1, separation * 2, (int)(nCells + 1));
            //    double[] XNODES3 = GenericBlas.Linspace(separation * 2, separation * 3, (int)(nCells / 3 + 1));

            //    List<double> xNodesList = new List<double>();
            //    xNodesList.AddRange(XNODES1.Take(XNODES1.Length - 1));
            //    xNodesList.AddRange(XNODES2.Take(XNODES2.Length - 1));
            //    xNodesList.AddRange(XNODES3);
            //    _xNodes = xNodesList.ToArray();

            //    // Y-NODES
            //    double sf3 = 0.80;
            //    var _yNodes2 = GenericBlas.Linspace(-radius_inlet, radius_inlet, nCells + 1).ToList(); // center
            //    List<double> yNodesTop = (GenericBlas.SinLinSpacing(radius_inlet, (ytop - radius_inlet) * 2 + radius_inlet, sf3, nCells * 1).ToList()); // Nodes corresponding to the oxidizer inlet, right part
            //    var myYnodesTop = yNodesTop.GetSubVector(0, yNodesTop.Count / 2 + 1); // Take only "bottom side" of node array
            //    var myYNodesBot = myYnodesTop.CloneAs();
            //    myYNodesBot.ScaleV(-1.0);
            //    Array.Reverse(myYNodesBot);
            //    List<double> list2 = new List<double>();
            //    list2.AddRange(myYNodesBot.Take(myYNodesBot.Length - 1).ToList());
            //    list2.AddRange(_yNodes2.Take(_yNodes2.Count() - 1).ToList());
            //    list2.AddRange(myYnodesTop.ToList());
            //    _yNodes = list2.ToArray();

            //    double[] CutOut1Point1 = new double[2] { separation, separation / 2 };
            //    double[] CutOut1Point2 = new double[2] { 0, -separation / 2 };

            //    double[] CutOut2Point1 = new double[2] { 2 * separation, -separation / 2 };
            //    double[] CutOut2Point2 = new double[2] { 3 * separation, separation / 2 };

            //    var CutOut1 = new BoSSS.Platform.Utils.Geom.BoundingBox(2);
            //    CutOut1.AddPoint(CutOut1Point1);
            //    CutOut1.AddPoint(CutOut1Point2);

            //    var CutOut2 = new BoSSS.Platform.Utils.Geom.BoundingBox(2);
            //    CutOut2.AddPoint(CutOut2Point1);
            //    CutOut2.AddPoint(CutOut2Point2);
            //    var CutOuts = new BoSSS.Platform.Utils.Geom.BoundingBox[] { CutOut1, CutOut2 };
            //    var grd = Grid2D.Cartesian2DGrid(_xNodes, _yNodes, CutOuts: CutOuts);

            //    grd.EdgeTagNames.Add(1, "Velocity_Inlet_CH4");
            //    grd.EdgeTagNames.Add(2, "Velocity_Inlet_O2");
            //    grd.EdgeTagNames.Add(3, "Pressure_Outlet");
            //    grd.EdgeTagNames.Add(4, "Wall");
            //    grd.DefineEdgeTags(GridEdgeTagFunc);

            //    //var gDat = new GridData(grd);
            //    //var em1 = gDat.GetBoundaryEdges();
            //    //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

            //    return grd;
            //};
            // initial values
            // ==============
            if(!restartOK) {
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
            // Species 0: Fuel
            // Species 1: Oxidizer
            // Species 2: Product1-
            // Species 3: Product2

            Func<double[], double, double> parabolaTotalFuel = (X, t) => (1.0 - Math.Pow(X[1] / radius_inlet, 2)) * uInFuel / uRef * 0 + 1 * uInFuel / C.uRef;
            Func<double[], double, double> parabolaTotalAir = (X, t) => -(1.0 - Math.Pow(X[1] / radius_inlet, 2)) * uInAir / uRef * 0 - 1 * uInAir / C.uRef;

            double sigma = 60; //regularization factor            
                Func<double[], double, double> RegularizedPlugFlowFuel= delegate(double[] X, double t) {
                    double velFuel = uInFuel / C.uRef;
                    double wall = 0;
                    double res;

                    if(X[1] > 0) {
                        double H = 0.5 * (1.0 + Math.Tanh(sigma * (X[1] - radius_inlet)));
                        res = velFuel * (1 - H) + wall * H;
                    } 
                    else {
                        double H = 0.5 * (1.0 + Math.Tanh(sigma * (X[1] - (-radius_inlet))));
                        res = wall * (1 - H) + velFuel * H;
                    }

                    return res;
                }

;
            Func<double[], double, double> RegularizedPlugFlowAir = (X, t) => -(1.0 - Math.Pow(X[1] / radius_inlet, 2)) * uInAir / uRef * 0 - 1 * uInAir / C.uRef;


            C.AddBoundaryValue("Pressure_Outlet");
            if(!C.usesimplifiedGeometry) {
                //C.AddBoundaryValue("NoSlipNeumann", VariableNames.VelocityX, X => 0.0);
                //C.AddBoundaryValue("NoSlipNeumann", VariableNames.VelocityY, X => 0.0);
                C.AddBoundaryValue("Wall", VariableNames.VelocityX, X => 0.0);
                C.AddBoundaryValue("Wall", VariableNames.VelocityY, X => 0.0);
                C.AddBoundaryValue("Wall", VariableNames.Temperature, X => 1.0);
            }

            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Velocity_d(0), RegularizedPlugFlowFuel);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Velocity_d(1), (X, t) => 0.0);

            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Velocity_d(0), RegularizedPlugFlowAir);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Velocity_d(1), (X, t) => 0.0);

            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Temperature, (X, t) => TemperatureInFuel / C.TRef);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction0, (X, t) => FuelInletMassFractions[0]); // CH4
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction1, (X, t) => FuelInletMassFractions[1]); // O2
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction2, (X, t) => FuelInletMassFractions[2]); // CO2
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction3, (X, t) => FuelInletMassFractions[3]); // H2O
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MixtureFraction, (X, t) => 1.0); //

            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Temperature, (X, t) => TemperatureInOxidizer / C.TRef);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction0, (X, t) => OxidizerInletMassFractions[0]); // CH4
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction1, (X, t) => OxidizerInletMassFractions[1]); // O2
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction2, (X, t) => OxidizerInletMassFractions[2]); // CO2
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction3, (X, t) => OxidizerInletMassFractions[3]); // H2O
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MixtureFraction, (X, t) => 0.0); //
            return C;
        }
    }
}
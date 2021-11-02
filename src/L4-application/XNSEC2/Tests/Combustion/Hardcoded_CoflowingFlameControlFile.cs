using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
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
        static public XNSEC_Control Full_CoFlowDiffusionFlame(int DGp = 2, int nCells = 7, int Refinements = 1, double velMultiplier = 0.4, string dbPath = @"C:\Databases\CoflowingFlame") {
            var C = LaminarCoFlowFlameHardcoded(DGp, nCells, velMultiplier, dbPath);
            C.ProjectName = "CoFlowDiffusionFlame";
            string jobName = C.ProjectName + "P" + DGp + "K" + nCells;
            C.SessionName = "Full_" + jobName;
            //C.SetAdaptiveMeshRefinement(amrLevel: Refinements, pseudoTimeStepsNo: Refinements);
            C.physicsMode = PhysicsMode.Combustion;
            C.HeatCapacityMode = MaterialLaw_MultipleSpecies.CpCalculationMode.mixture;
            C.VariableOneStepParameters = false; /////////////////////////////////////////////////////////////////////////////////
            // Select the database
            DatabaseInfo dbi = DatabaseInfo.Open(C.DbPath);
            C.Timestepper_LevelSetHandling = Solution.XdgTimestepping.LevelSetHandling.None;
            C.AdaptiveMeshRefinement = false;
            //C.activeAMRlevelIndicators.Add(new AMR_onFlameSheet(C.zSt, 2));
            //C.activeAMRlevelIndicators.Add(new AMR_BasedOnPerssonSensor(VariableNames.Temperature, 3));
            C.activeAMRlevelIndicators.Add(new AMR_onProblematicPoints(C.troubledPoints, 3));
            C.activeAMRlevelIndicators.Add(new AMR_BasedOnVariableLimits(VariableNames.Temperature, new double[] { 0.98, C.AdiabaticTemperature },3));
            C.activeAMRlevelIndicators.Add(new AMR_BasedOnVariableLimits(VariableNames.MassFraction0, new double[] { 0.0, 1.0 },3));
            C.activeAMRlevelIndicators.Add(new AMR_BasedOnVariableLimits(VariableNames.MassFraction1, new double[] { 0.0, 1.0 },3));

            C.AMR_startUpSweeps = 1;

            C.UseMixtureFractionsForCombustionInitialization = true;

            //C.NonLinearSolver.MinSolverIterations = 1;
            //C.NonLinearSolver.MaxSolverIterations = 1;

            //C.NoOfTimesteps = 5;


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
        static public XNSEC_Control FS_CoFlowDiffusionFlame(int DGp = 2, int nCells = 7, double velMultiplier = 0.4, string dbPath = @"C:\Databases\CoflowingFlame") {
            var C = LaminarCoFlowFlameHardcoded(DGp, nCells, velMultiplier, dbPath);
            C.physicsMode = PhysicsMode.MixtureFraction;
            C.ProjectName = "CoFlowDiffusionFlame";
            string name = C.ProjectName + "P" + DGp + "K" + nCells /*+ "Ne" + NewtonConvergenceCriterion*/;
            C.SessionName = "FS_" + name;
            C.VariableOneStepParameters = false;
            C.AdaptiveMeshRefinement = false;

            return C;
        }

        /// <summary>
        /// CounterDifussionFlameTest2 Flame
        /// </summary>
        static public XNSEC_Control LaminarCoFlowFlameHardcoded(int DGp = 3, int nCells = 16, double velMultiplier = 1, string _dbPath = @"C:\Databases\BoSSS_DB") {
            Console.WriteLine("===========================================");
            Console.WriteLine("CoFlow Diffusion Flame");
            Console.WriteLine("===========================================");

            // Parameters
            // ==============
            //var Conc1 = new double[] { 0.2, 0.8 }; // inlet concentrations...

            string dbPath = _dbPath;

            // Problem Definition
            //===================
            double TemperatureInFuel = 300;
            double TemperatureInOxidizer = 300;

            double massFuelIn = 0.2400 * velMultiplier; //kg/m2s
            double massAirIn = 0.2400 * 3 * velMultiplier; //kg/m2s
            double AtmPressure = 101325; // Pa
            var CC = new ChemicalConstants();

            double[] MWs = new double[] { CC.MW_CH4, CC.MW_O2, CC.MW_CO2, CC.MW_H2O, CC.MW_N2 };
            double[] FuelInletMassFractions = new double[] { 0.2, 0.0, 0.0, 0.0, 0.8 };
            double[] OxidizerInletMassFractions = new double[] { 0.0, 0.23, 0.0, 0.0, 0.77 };
            //double[] FuelInletMassFractions = new double[] { 1.0, 0.0, 0.0, 0.0, 0.0 };
            //double[] OxidizerInletMassFractions = new double[] { 1.0*0, 1.0, 0.0, 0.0, 0.0 };

            double factor = 1.0;
            double r = 0.635 / 100 * factor; // Radius inner cylinder, m
            double R = 2.54 / 100 * factor; // Radius outter cylinder, m

            double LRef = r;
            double mwFuel = CC.getAvgMW(MWs, FuelInletMassFractions);
            double mwAir = CC.getAvgMW(MWs, OxidizerInletMassFractions);
            double densityAirIn = AtmPressure * mwAir / (CC.R_gas * TemperatureInOxidizer * 1000); // Kg/m3. ok
            double densityFuelIn = AtmPressure * mwFuel / (CC.R_gas * TemperatureInFuel * 1000); // Kg/m3. ok
            double uInFuel = (4.5 / 100) * 0.1; // massFuelIn / densityFuelIn; // Avg value, m/s
            //double uInAir = (9.88 / 100) * 0.5; // 0.6162;//  massAirIn / densityAirIn; // Avg value m/s
            double uInAir = uInFuel; // Burke schuhmann flame

            //double uInAir = Math.Sqrt(densityFuelIn / densityAirIn )* uInFuel ;

            // Reference values
            //===================
            // Basic units to be used: Kg, m, s, mol, pa,

            double TRef = TemperatureInOxidizer;// Reference temperature  is the inlet temperature, (K)
            double pRef = AtmPressure; // Pa
            double uRef = Math.Max(uInFuel, uInAir); //

            XNSEC_Control C = new XNSEC_Control(DGp, pRef, uRef, TRef, true, LRef, FuelInletMassFractions, OxidizerInletMassFractions);
            C.physicsMode = PhysicsMode.Combustion;
            C.rhoOne = false;
            C.MatParamsMode = MaterialParamsMode.Constant;
            C.NonLinearSolver.MaxSolverIterations = 30;
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            int numOfRefinements = 1;
            bool restartOK = false;

            C.GravityDirection = new double[] { 0.0, 0.0, 0.0 }; //No gravity.
            C.SetSaveOptions(dbPath, savePeriod: 1);
            C.PhysicalParameters.IncludeConvection = true;

            // Geometry
            // All lenghts are non.dimensionalized with Lref = fuel inlet radius (0.2 cm)
            double zlength = 20.0 / 100 * factor;// meters
            double xmin = -R / LRef;
            double xmax = +R / LRef;
            double ymin = 0;
            double ymax = zlength / LRef;
            double rAd = r / LRef;
            double R_Outer = 7.5 / 100;
            double R_Outer_Ad = R_Outer / LRef;

            double leftmidpoint = (xmin - rAd) / 2;
            double rightmidpoint = (xmax + rAd) / 2;

            C.troubledPoints = new double[][] {
            new double[]{ rAd, 0.0 },
            new double[]{-rAd, 0.0 }
            };
            //double[][] troubledPoints = new double[][] {
            // new double[]{ rAd, 0.0 -eps },
            // new double[]{ rAd, 0.0  +eps},
            // new double[]{-rAd, 0.0 -eps},
            // new double[]{-rAd, 0.0 +eps }
            // };

            // Grid declaration
            // ===============
            Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
                double x = X[0];
                double y = X[1];

                //Edge tags
                //1: Velocity inlet O_2
                //2: Velocity inlet CH_4
                //3: Pressure outlet

                //Inlet fuel
                if(Math.Abs(y - ymin) < 1e-8 && Math.Abs(x - 0.0) < rAd + 1e-8)
                    return 2;
                //Inlet oxidizer
                if((Math.Abs(y - ymin) < 1e-8 && ((Math.Abs(x - rightmidpoint) < ((xmax - rAd) * 0.5 + rAd) + 1e-8) || Math.Abs(x - leftmidpoint) < Math.Abs((xmin - rAd) * 0.5 + rAd) + 1e-8)))
                    return 1;

                // Pressure outlet

                if((Math.Abs(y - ymax) < 1e-8))
                    return 3;

                // Pressure outlet

                if((Math.Abs(x - xmin) < 1e-8) || (Math.Abs(x - xmax) < 1e-8))
                    return 4;
                else throw new ArgumentOutOfRangeException();
            };


            //Func<double[], byte> GridEdgeTagFunc = delegate (double[] X) {
            //    double x = X[0];
            //    double y = X[1];

            //    //Edge tags
            //    //1: Velocity inlet O_2
            //    //2: Velocity inlet CH_4
            //    //3: Pressure outlet
            //    //4: Wall

            //    if(Math.Abs(y - ymin) < 1e-8) {
            //        if(Math.Abs(x - 0.0) < rAd + 1e-8) {
            //            return 1;
            //            //} else if(((Math.Abs(x - rightmidpoint) < ((xmax - rAd) * 0.5 + rAd) + 1e-8) || Math.Abs(x - leftmidpoint) < Math.Abs((xmin - rAd) * 0.5 + rAd) + 1e-8)) {
            //        } else if(x > xmax + 1e-8 || x < xmin- 1e-8) {
            //            return 4;
            //        } else {

            //            return 2; // air outlet
            //        }
            //    }

            //    // Pressure outlet

            //    if((Math.Abs(y - ymax) < 1e-8))
            //        return 3;

            //    // Pressure outlet

            //    if((Math.Abs(x - R_Outer_Ad) < 1e-8) || (Math.Abs(x+ R_Outer_Ad) < 1e-8))
            //        return 3;
            //    else throw new ArgumentOutOfRangeException();

            //};




            double stretchfactorY = 0.98 * 1;

            // initial values
            // ==============
            if(!restartOK) {
                C.GridFunc = delegate {
                    //List<double> list = new List<double>();
                    //list.AddRange(myXnodesLeft);
                    //list.AddRange(myXnodesRight);
                    //double[] _xNodes = list.ToArray();
                    double sf1 = 0.97 * 1 * 0;
                    double sf2 = 0.95 * 1 * 0;
                    double sf3 = 0.97 * 1 * 0;
                    int n1 = (int)2.5 * nCells;
                    int n2 = (int)1.0 * nCells;
                    int n3 = (int)2.5 * nCells;
                    //var xNodes1 = GenericBlas.SinLinSpacing(xmin, -rAd, sf1, n1 + 1);
                    //var xNodes2 = GenericBlas.SinLinSpacing(-rAd, rAd, sf2, n2 + 1);
                    //var xNodes3 = GenericBlas.SinLinSpacing(rAd, xmax, sf3, n3 + 1);

                    List<double> xNodes2 = (GenericBlas.SinLinSpacing(-rAd, rAd, sf2, n2 + 1)).ToList(); // nodes corresponding to the fuel inlet
                    List<double> xNodes3 = (GenericBlas.SinLinSpacing(rAd, (xmax - rAd) * 2 + rAd, sf3, n1 * 2 + 1).ToList()); // Nodes corresponding to the oxidizer inlet, right part
                    var myXnodes3 = xNodes3.GetSubVector(0, xNodes3.Count / 2 + 1); // Take only "left side" of node array
                    var myxNodes1 = myXnodes3.CloneAs();
                    myxNodes1.ScaleV(-1.0);
                    Array.Reverse(myxNodes1);

                    List<double> list2 = new List<double>();
                    list2.AddRange(myxNodes1.Take(n1 + 0).ToList());
                    list2.AddRange(xNodes2.Take(n2 + 0).ToList());
                    list2.AddRange(myXnodes3.Take(n3 + 1).ToList());
                    double[] _xNodes = list2.ToArray();
                    //Debug.Assert(_xNodes.Contains(xmin));
                    //Debug.Assert(_xNodes.Contains(-rAd));
                    //Debug.Assert(_xNodes.Contains(rAd));
                    //Debug.Assert(_xNodes.Contains(xmax));
                    var _yNodes = GenericBlas.SinLinSpacing(ymin, ymax * 2, stretchfactorY, (2 * nCells) * 4 + 1);
                    var myYnodes = _yNodes.GetSubVector(0, _yNodes.Length / 2 + 1); // I just want a fine mesh in the bottom part of the grid.
                                                                                    //var _xNodes = GenericBlas.Linspace(xmin, xmax, 2 * meshScaling + 1);
                                                                                    //var myYnodes = GenericBlas.Linspace(ymin, ymax, 2 * meshScaling + 1);
                    var grd = Grid2D.Cartesian2DGrid(_xNodes, myYnodes, periodicX: false);
                    grd.EdgeTagNames.Add(1, "Velocity_Inlet_O2");
                    grd.EdgeTagNames.Add(2, "Velocity_Inlet_CH4");
                    grd.EdgeTagNames.Add(3, "Pressure_Outlet");
                    //grd.EdgeTagNames.Add(4, "wall");
                    //grd.EdgeTagNames.Add(4, "NoSlipNeumann");
                    grd.EdgeTagNames.Add(4, "Velocity_Inlet_outer"); // We want a constant velocity field.

                    grd.DefineEdgeTags(GridEdgeTagFunc);

                    //var gDat = new GridData(grd);
                    //var em1 = gDat.GetBoundaryEdges();
                    //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                    return grd;
                };

                //C.GridFunc = delegate {

                //    double sf1 = 0.97 * 1 * 0;
                //    double sf2 = 0.95 * 1 * 0;
                //    double sf3 = 0.97 * 1 * 0;
                //    int n1 = (int)2.5 * nCells;
                //    int n2 = (int)1.0 * nCells;
                //    int n3 = (int)2.5 * nCells;
                //    //var xNodes1 = GenericBlas.SinLinSpacing(xmin, -rAd, sf1, n1 + 1);
                //    //var xNodes2 = GenericBlas.SinLinSpacing(-rAd, rAd, sf2, n2 + 1);
                //    //var xNodes3 = GenericBlas.SinLinSpacing(rAd, xmax, sf3, n3 + 1);

                //    List<double> myXNodes1 = (GenericBlas.SinLinSpacing(0, rAd * 2, sf3, n2 * 2 + 1).ToList()); // Nodes corresponding to the oxidizer inlet, right part
                //    var xNodes1 = myXNodes1.GetSubVector(0, myXNodes1.Count / 2 + 1); // Take only "left side" of node array

                //    List<double> myXNodes2 = (GenericBlas.SinLinSpacing(rAd, (xmax - rAd) * 2 + rAd, sf3, n2 * 2 + 1).ToList()); // Nodes corresponding to the oxidizer inlet, right part
                //    var xNodes2 = myXNodes2.GetSubVector(0, myXNodes2.Count / 2 + 1); // Take only "left side" of node array

                //    var xNodes3 = GenericBlas.Linspace(xmax, R_Outer_Ad, n1);

                //    List<double> RightSide = new List<double>();
                //    RightSide.AddRange(xNodes1.SkipLast(1).ToList());
                //    RightSide.AddRange(xNodes2.SkipLast(1).ToList());
                //    RightSide.AddRange(xNodes3.SkipLast(0).ToList());
                //    var RightSideArr = RightSide.ToArray();

                //    var LeftSideArr = RightSideArr.CloneAs();
                //    LeftSideArr.ScaleV(-1.0);
                //    Array.Reverse(LeftSideArr);

                //    List<double> _xNodes = new List<double>();
                //    _xNodes.AddRange(LeftSideArr.SkipLast(1));
                //    _xNodes.AddRange(RightSide);

                //    var _yNodes = GenericBlas.SinLinSpacing(ymin, ymax * 2, stretchfactorY, (2 * nCells) * 4 + 1);
                //    var myYnodes = _yNodes.GetSubVector(0, _yNodes.Length / 2 + 1); // I just want a fine mesh in the bottom part of the grid.
                //                                                                    //var _xNodes = GenericBlas.Linspace(xmin, xmax, 2 * meshScaling + 1);
                //                                                                    //var myYnodes = GenericBlas.Linspace(ymin, ymax, 2 * meshScaling + 1);
                //    var grd = Grid2D.Cartesian2DGrid(_xNodes.ToArray(), myYnodes, periodicX: false);

                //    grd.EdgeTagNames.Add(1, "Velocity_Inlet_O2");
                //    grd.EdgeTagNames.Add(2, "Velocity_Inlet_CH4");
                //    grd.EdgeTagNames.Add(3, "Pressure_Outlet");
                //    grd.EdgeTagNames.Add(4, "wall");

                //    grd.DefineEdgeTags(GridEdgeTagFunc);

                //    //var gDat = new GridData(grd);
                //    //var em1 = gDat.GetBoundaryEdges();
                //    //em1.SaveToTextFile("alledges.csv", false, (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => (double)gDat.iGeomEdges.EdgeTags[GeomItemIndex]);

                //    return grd;
                //};



                // initial values
                //// ==============
                //Func<double[],  double> TempInit = (X) => (Math.Abs((X[0] - 0.2  < 0) && (X[1] == 0)) ? 5 : 1);

                C.InitialValues_Evaluators.Add(VariableNames.Pressure, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(0), X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(1), X => 0);
                C.InitialValues_Evaluators.Add(VariableNames.Temperature, X => TRef / TRef/*TempInit*/);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction1, X => 1.0); // initially everywhere just O2
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction2, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction3, X => 0.0);
                C.InitialValues_Evaluators.Add(VariableNames.MassFraction4, X => 1.0);
                C.InitialValues_Evaluators.Add(VariableNames.MixtureFraction, X => 0.0);
                #region
                Func<double[], double> BSFlameSolut = delegate (double[] _x) {
                    double x = _x[0] * C.LRef;
                    double z = _x[1] * C.LRef;
                    double D = C.DRef;
                    double v = uInAir;
                    double a = r;
                    double b = R;

                    double c = a / b;

                    double psi = x / b;

                    double eta = z * D / (v * b * b);

                    int N = 100; // number of terms considered for approximation of infinite sum
                    double zBs = c;
                    for(int i = 1; i < N; i++) {
                        zBs += 2 * ((1.0 / (i * Math.PI)) * Math.Sin(i * Math.PI * c) * Math.Cos(i * Math.PI * psi) * Math.Exp(-1 * i * i * Math.PI * Math.PI * eta));
                    }

                    return zBs;
                };

                C.InitialValues_Evaluators.Add("Z_an", BSFlameSolut);
                #endregion
            } else {
                C.RestartInfo = new Tuple<Guid, TimestepNumber>(new Guid(C.ID), new TimestepNumber(-1));
            }

            // boundary conditions
            // ===================
            // Species 0: Fuel
            // Species 1: Oxidizer
            // Species 2: Product1
            // Species 3: Product2

            // Analytical velocity profile for flow in annulus of a cylinder
            double k = r / R;
            double k1 = (1 - k * k * k * k) / (1 - k * k);
            double k2 = (1 - k * k) / Math.Log(1 / k);
            double kappaTerm = k1 - k2;//((1 -Math.Pow(k,4)) / (1 - k * k) - (1 - k - k) / (Math.Log(1 / k)));
            double A = 2 * uInAir / kappaTerm;

            Func<double[], double, double> parabolaVelocityAir = (X, t) => (X[0] > 0 ?
                        (A * (1 - Math.Pow(X[0] / xmax, 2) - (1 - k * k) / Math.Log(1 / k) * Math.Log(xmax / X[0]))) / uRef : // Right Parabola
                        (A * (1 - Math.Pow(X[0] / xmin, 2) - (1 - k * k) / Math.Log(1 / k) * Math.Log(xmin / X[0]))) / uRef   // Left Parabola
           );

            Func<double[], double, double> parabolaVelocityFuel = (X, t) => (1.0 - Math.Pow(X[0] / rAd, 2)) * uInFuel / uRef;

            Func<double[], double, double> parabolaTotal = (X, t) => (1.0 - (1 / (xmax * xmax)) * Math.Pow(X[0], 2));

            Func<double[], double, double> TempFunc = (X, t) => t < 0.04 ? 3 : TRef / TRef;

            double constantvelocity = 1;
            C.rhoOne = true; // the whole velociy field should be constant...

            C.AddBoundaryValue("Pressure_Outlet");
            ////C.AddBoundaryValue("Wall", VariableNames.MixtureFraction, X => 0.0); // Kept at reference temperature
            //C.AddBoundaryValue("NoSlipNeumann");

            //C.AddBoundaryValue("NoSlipNeumann", VariableNames.VelocityX, X => 0.0);
            //C.AddBoundaryValue("NoSlipNeumann", VariableNames.VelocityY, X => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Velocity_d(0), (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Velocity_d(1), (X, t) => uInFuel / uRef);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.Temperature, (X, t) => 1.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction0, (X, t) => C.YFuelInlet);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction1, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction2, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction3, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MassFraction4, (X, t) => 1.0 - C.YFuelInlet);
            C.AddBoundaryValue("Velocity_Inlet_CH4", VariableNames.MixtureFraction, (X, t) => 1.0);

            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Velocity_d(0), (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Velocity_d(1), (X, t) => uInAir / uRef);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.Temperature, (X, t) => 1.0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction0, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction1, (X, t) => C.YOxInlet);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction2, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction3, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MassFraction4, (X, t) => 1.0 - C.YOxInlet);
            C.AddBoundaryValue("Velocity_Inlet_O2", VariableNames.MixtureFraction, (X, t) => 0.0);

            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.Velocity_d(0), (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.Velocity_d(1), (X, t) => uInAir / uRef);
            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.Temperature, (X, t) => 1.0);
            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.MassFraction0, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.MassFraction1, (X, t) => C.YOxInlet);
            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.MassFraction2, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.MassFraction3, (X, t) => 0.0);
            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.MassFraction4, (X, t) => 1.0 - C.YOxInlet);
            C.AddBoundaryValue("Velocity_Inlet_outer", VariableNames.MixtureFraction, (X, t) => 0.0);

            return C;
        }
    }
}
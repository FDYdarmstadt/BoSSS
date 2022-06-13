using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.TestCases;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using System.IO;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.LsTest {
    public static partial class LevelSetUnitTests {

        public static void LevelSetZalesakDisc(
            [Values(2, 3, 4)] int LSdegree,
            [Values(0, 1, 2)] int AMRlevel,
            [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
            [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling,
            [Values(false, true)] bool reversed) {
            Solution.Application.DeleteOldPlotFiles();

            int gridResolution;
            gridResolution = 1;
            //switch (LSdegree) {
            //    case 2: gridResolution = 3; break;
            //    case 3: gridResolution = 2; break;
            //    //case 4: gridResolution = 1; break;
            //    default:
            //        gridResolution = 1; break;
            //}

            if (LSdegree == 4 && AMRlevel > 0 && levelSetEvolution == LevelSetEvolution.StokesExtension)
                return;

            var Tst = new LevelSetZalesakDiscTest(2, LSdegree, reversed);
            var C = LSTstObj2CtrlObj(Tst, int.MaxValue, levelSetEvolution, levelSetHandling, gridResolution, AMRlevel);
            C.PostprocessingModules.Add(new SolverWithLevelSetUpdaterLogging());

            string IO = $"LSAdvectionTest2D-deg{LSdegree}-amrLvl{AMRlevel}-lsEvo{levelSetEvolution}-rev{reversed}-grdRes{gridResolution}";
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 3;

            LevelSetUnitTests.LevelSetTest(Tst, C);

            Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!
        }
        public static SolverWithLevelSetUpdaterTestControl SwirlingFlowTemporalConvergence(int degree, int gridRes, int tempRes, LevelSetEvolution lsEvo, string ProjectName, string dbPath){
            var Tst = new LevelSetSwirlingFlowTest(2, degree, false, 0.1);
            var C = LSTstObj2CtrlObj(Tst, int.MaxValue, lsEvo, LevelSetHandling.LieSplitting, gridRes, 0, tempRes);

            C.SessionName = "SwirlingFlow_T_p" + degree + "_H" + gridRes + "_t" + tempRes + "_Evo" + lsEvo.ToString();
            C.ProjectName = ProjectName;
            C.savetodb = dbPath != null;
            C.DbPath = dbPath;
            C.saveperiod = 50; // in principal we only need the last timestep, save now and then for potential restarts

            var db = DatabaseInfo.Open(dbPath);
            var grd = C.GridFunc();
            db.Controller.DBDriver.SaveGridIfUnique(ref grd, out bool found, db);
            if (found) {
                Console.WriteLine("Found equivalent grid in database, grid will not be saved");
            }
            C.SetGrid(grd);
            C.GridFunc = null;

            C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Res", gridRes));
            C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Degree", degree));
            C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Evo", lsEvo.ToString()));
            C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("dt", C.dtFixed));

            return C;
        }
        public static SolverWithLevelSetUpdaterTestControl SwirlingFlowSpatialConvergence(int degree, int gridRes, LevelSetEvolution lsEvo, string ProjectName, string dbPath) {
            var Tst = new LevelSetSwirlingFlowTest(2, degree, false, 8);
            var C = LSTstObj2CtrlObj(Tst, int.MaxValue, lsEvo, LevelSetHandling.LieSplitting, gridRes, 0);

            C.SessionName = "SwirlingFlow_H_p" + degree + "_H" + gridRes + "_Evo" + lsEvo.ToString();
            C.ProjectName = ProjectName;
            C.savetodb = dbPath != null;
            C.DbPath = dbPath;
            C.saveperiod = 50; // in principal we only need the last timestep, save now and then for potential restarts
            C.dtFixed = 1.0 / (100 * degree * degree * 4); // gridres max is 4, and tempres is 4 for all simulations (all using the same timestep)

            var db = DatabaseInfo.Open(dbPath);
            var grd = C.GridFunc();
            db.Controller.DBDriver.SaveGridIfUnique(ref grd, out bool found, db);
            if (found) {
                Console.WriteLine("Found equivalent grid in database, grid will not be saved");
            }
            C.SetGrid(grd);
            C.GridFunc = null;

            C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Res", gridRes));
            C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Degree", degree));
            C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("Evo", lsEvo.ToString()));
            C.Paramstudy_CaseIdentification.Add(new Tuple<string, object>("dt", C.dtFixed));

            return C;
        }

        public static void LevelSetSwirlingFlow(
            [Values(2, 3, 4)] int LSdegree,
            [Values(0, 1, 2)] int AMRlevel,
            [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
            [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling,
            [Values(false, true)] bool reversed) {
            Solution.Application.DeleteOldPlotFiles();

            int gridResolution = 1;
            //switch (LSdegree) {
            //    case 2: gridResolution = 3; break;
            //    case 3: gridResolution = 2; break;
            //    //case 4: gridResolution = 1; break;
            //    default:
            //        gridResolution = 1; break;
            //}

            if (LSdegree == 4 && AMRlevel > 0 && levelSetEvolution == LevelSetEvolution.StokesExtension)
                return;

            var Tst = new LevelSetSwirlingFlowTest(2, LSdegree, reversed);
            var C = LSTstObj2CtrlObj(Tst, int.MaxValue, levelSetEvolution, levelSetHandling, gridResolution, AMRlevel);
            //C.PostprocessingModules.Add(new SolverWithLevelSetUpdaterLogging());
            //C.PostprocessingModules.Add(new SolverWithLevelSetUpdaterQualityLogging());
            //C.PostprocessingModules.Add(new LevelSetTestLogging() { LogPeriod = C.NoOfTimesteps, interface_measure = 2 * Math.PI * 0.15, exact_levelset = (X,t) => Math.Sqrt(Math.Pow(X[0] - 0.5,2) + Math.Pow(X[1]-0.75,2)) - 0.15});

            var db = BoSSS.Foundation.IO.DatabaseInfo.CreateOrOpen("./TempLsTestDB");
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid("b9d6231e-44be-488d-8de2-d6fe7df6fc99"), -1); 
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid("18389f0e-e32a-4ae5-b440-400b851e7dd1"), -1);
            //C.RestartInfo = new Tuple<Guid, BoSSS.Foundation.IO.TimestepNumber>(new Guid("55407481-b2d1-48ff-ba24-b865caed0e84"), 3078);            
            //C.GridFunc = null;

            C.SetDatabase(db);
            C.savetodb = true;
            C.saveperiod = 10;
            //C.rollingSaves = true;

            //string IO = $"LSAdvectionTest2D-deg{LSdegree}-amrLvl{AMRlevel}-lsEvo{levelSetEvolution}-rev{reversed}-grdRes{gridResolution}";
            //C.ImmediatePlotPeriod = 1;
            //C.SuperSampling = 3;

            LevelSetUnitTests.LevelSetTest(Tst, C);

            Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!
        }        

        public static void LevelSetSwirlingFlowConvergenceTest(
            [Values(2, 3, 4)] int LSdegree,
            [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
            [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling,
            [Values(false, true)] bool reversed
            ) {
            Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!

            int[] GridResolutionS = new int[] { 1, 2, 3, 4 };

            var Tst = new LevelSetSwirlingFlowTest(2, LSdegree, reversed);
            var CS = new SolverWithLevelSetUpdaterTestControl[GridResolutionS.Length];
            for (int i = 0; i < GridResolutionS.Length; i++) {
                var C = LSTstObj2CtrlObj(Tst, int.MaxValue, levelSetEvolution, levelSetHandling, GridResolutionS[i], 0);
                C.PostprocessingModules.Add(new SolverWithLevelSetUpdaterLogging());
                CS[i] = C;
            }

            var RegressionBounds = new (string Name, double Slope, double intercept, double interceptTol)[6];
            RegressionBounds[0] = (VariableNames.LevelSetCGidx(0), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[1] = (VariableNames.LevelSetDGidx(0), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[2] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Volume#A"), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[3] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Volume#B"), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[4] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Surface-Error"), LSdegree, 0.740, 0.1);
            RegressionBounds[5] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Contour-Error"), 1.0, 0.740, 0.1);

            LevelSetConvergenceTest(Tst, CS, true, RegressionBounds);
        }

        public static void LevelSetZalasakDiscConvergenceTest(
            [Values(2, 3, 4)] int LSdegree,
            [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
            [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling,
            [Values(false, true)] bool reversed
            ) {
            Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!

            int[] GridResolutionS = new int[] { 1, 2, 3, 4 };

            var Tst = new LevelSetZalesakDiscTest(2, LSdegree, reversed);
            var CS = new SolverWithLevelSetUpdaterTestControl[GridResolutionS.Length];
            for (int i = 0; i < GridResolutionS.Length; i++) {
                var C = LSTstObj2CtrlObj(Tst, int.MaxValue, levelSetEvolution, levelSetHandling, GridResolutionS[i], 0);
                C.PostprocessingModules.Add(new SolverWithLevelSetUpdaterLogging());
                CS[i] = C;
            }

            var RegressionBounds = new (string Name, double Slope, double intercept, double interceptTol)[6];
            RegressionBounds[0] = (VariableNames.LevelSetCGidx(0), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[1] = (VariableNames.LevelSetDGidx(0), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[2] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Volume#A"), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[3] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Volume#B"), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[4] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Surface-Error"), LSdegree, 0.740, 0.1);
            RegressionBounds[5] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Contour-Error"), 1.0, 0.740, 0.1);

            LevelSetConvergenceTest(Tst, CS, true, RegressionBounds);
        }

        /// <summary>
        /// Either smooth or sdf phi, compare the approximation error
        /// </summary>
        /// <param name="LSdegree"></param>
        /// <param name="levelSetEvolution"></param>
        /// <param name="levelSetHandling"></param>
        /// <param name="reversed"></param>
        public static void LevelSetCircleProjectionConvergenceTest(
            [Values(2, 3, 4)] int LSdegree,
            [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
            [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling,
            [Values(false, true)] bool reversed,
            [Values(false, true)] bool signeddistance
            ) {
            Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!

            int[] GridResolutionS = new int[] { 1, 2, 3, 4 };

            var Tst = new LevelSetCircleProjectionTest(2, LSdegree, reversed, signeddistance);
            var CS = new SolverWithLevelSetUpdaterTestControl[GridResolutionS.Length];
            for (int i = 0; i < GridResolutionS.Length; i++) {
                var C = LSTstObj2CtrlObj(Tst, int.MaxValue, levelSetEvolution, levelSetHandling, GridResolutionS[i], 0);
                CS[i] = C;
            }

            var RegressionBounds = new (string Name, double Slope, double intercept, double interceptTol)[6];
            RegressionBounds[0] = (VariableNames.LevelSetCGidx(0), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[1] = (VariableNames.LevelSetDGidx(0), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[2] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Volume#A"), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[3] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Volume#B"), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[4] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Surface-Error"), LSdegree, 0.740, 0.1);
            RegressionBounds[5] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Contour-Error"), 1.0, 0.740, 0.1);

            LevelSetConvergenceTest(Tst, CS, true, RegressionBounds);
        }

        /// <summary>
        /// Given as SDF, Bramble Hilbert gives us a prediction on the aproximation error convergence
        /// </summary>
        /// <param name="LSdegree"></param>
        /// <param name="levelSetEvolution"></param>
        /// <param name="levelSetHandling"></param>
        /// <param name="reversed"></param>
        public static void LevelSetCubeProjectionConvergenceTest(
            [Values(2, 3, 4)] int LSdegree,
            [Values(LevelSetEvolution.FastMarching, LevelSetEvolution.StokesExtension)] LevelSetEvolution levelSetEvolution,
            [Values(LevelSetHandling.LieSplitting)] LevelSetHandling levelSetHandling,
            [Values(false, true)] bool reversed
            ) {
            Solution.Application.DeleteOldPlotFiles(); // delete plot files if we don't throw an exception!

            int[] GridResolutionS = new int[] { 1, 2, 3 };

            var Tst = new LevelSetCubeProjectionTest(3, LSdegree, reversed);
            var CS = new SolverWithLevelSetUpdaterTestControl[GridResolutionS.Length];
            for (int i = 0; i < GridResolutionS.Length; i++) {
                var C = LSTstObj2CtrlObj(Tst, int.MaxValue, levelSetEvolution, levelSetHandling, GridResolutionS[i], 0);
                CS[i] = C;
            }

            var RegressionBounds = new (string Name, double Slope, double intercept, double interceptTol)[6];
            RegressionBounds[0] = (VariableNames.LevelSetCGidx(0), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[1] = (VariableNames.LevelSetDGidx(0), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[2] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Volume#A"), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[3] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Volume#B"), LSdegree + 1, 0.740, 0.1);
            RegressionBounds[4] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Surface-Error"), LSdegree, 0.740, 0.1);
            RegressionBounds[5] = (VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(0), "Contour-Error"), 1.0, 0.740, 0.1);

            LevelSetConvergenceTest(Tst, CS, true, RegressionBounds);
        }

        /// <summary>
        /// Basic test case for advecting the level-set field
        /// in a constant velocity field
        /// </summary>
        public class LevelSetZalesakDiscTest : ILevelSetTest {

            protected double T = 3; // how many cycles, later 3
            protected double sign;

            /// <summary>
            /// ctor
            /// </summary>
            public LevelSetZalesakDiscTest(int spatDim, int LevelSetDegree, bool reversed) {
                this.SpatialDimension = spatDim;
                this.LevelsetPolynomialDegree = LevelSetDegree;
                sign = (reversed) ? -1 : 1;
            }

            public double dt {
                get {
                    return -1.0;    // will be set in LevelSetTest() according to level set cfl 
                }
            }

            /// <summary>
            /// computes the timestep size according to the level-set CFL condition
            /// </summary>
            /// <param name="Resolution"></param>
            /// <param name="LSdegree"></param>
            /// <returns></returns>
            public double ComputeTimestep(int Resolution, int LSdegree, int AMRlevel, int temporalResolution) {
                int gridCells1D = (25 * Resolution) * (AMRlevel + 1);
                double h = 1.0 * 1.0 / (double)gridCells1D;
                double dt = h / (Math.Sqrt(2) * Math.PI); // this is grid width divided by maximum velocity
                dt /= (double)(LSdegree * LSdegree);

                int timesteps = temporalResolution * Math.Max((int)Math.Ceiling(T / dt), 1); // make sure that the singular point in time is exactly on one timestep

                return (T / (double)timesteps);
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public double getEndTime() {
                return T;
            }

            /// <summary>
            /// creates a square in 2D
            /// </summary>
            /// <param name="Resolution"></param>
            /// <returns></returns>
            public virtual GridCommons CreateGrid(int Resolution) {
                if (Resolution < 1)
                    throw new ArgumentException();

                var xNodes = GenericBlas.Linspace(0, 1, 25 * Resolution + 1);
                var yNodes = GenericBlas.Linspace(0, 1, 25 * Resolution + 1);

                GridCommons grd;
                switch (SpatialDimension) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                        break;
                    case 3:                       
                    default:
                        throw new ArgumentOutOfRangeException();
                }

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(delegate (double[] X) {
                    return 1;
                });

                return grd;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public virtual IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
                var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

                config.Add("wall", new AppControl.BoundaryValueCollection());

                return config;
            }

            /// <summary>
            /// Level-Set movement.
            /// </summary>
            public Func<double[], double, double>[] GetPhi() {
                return new Func<double[], double, double>[] { delegate (double[] X, double time) {

                double[] cutoutx = new double[] { -0.025, 0.025 };
                double cutouty = -0.1125;
                double Radius = 0.15;
                double phi = Math.PI;

                ZalesaksDisk zalesaksDisk = new ZalesaksDisk(cutoutx, cutouty, Radius);
                BoSSS.Platform.Utils.Geom.Rotation2D rotation = new Platform.Utils.Geom.Rotation2D(phi);
                Func<double[], double, double> PhiFunc;

                switch (SpatialDimension) {
                    case 2: {
                        PhiFunc = new Func<double[], double, double>((X, t) => zalesaksDisk.SignedDistanceLevelSet(rotation.Transform(new double[] { X[0] - 0.5, X[1] - 0.75 }, phi + 2 * Math.PI * t)));
                        break; 
                    }
                    case 3:
                    default:
                    throw new ArgumentOutOfRangeException();
                }
                return sign * PhiFunc(X,time);

            }};
            }

            /// <summary>
            /// velocity: solid body rotation
            /// </summary>
            public Func<double[], double, double>[][] GetU() {
                var Ret = new Func<double[], double, double>[this.NoOfLevelsets][];
                switch (SpatialDimension) {
                    case 2:
                        Ret[0] = new Func<double[], double, double>[] { (X, t) => -2 * Math.PI * (X[1] - 0.5), (X, t) => 2 * Math.PI * (X[0] - 0.5) };
                        break;
                    case 3:                       
                    default:
                        throw new ArgumentOutOfRangeException();
                }
                return Ret;
            }

            public int LevelsetPolynomialDegree {
                get;
                private set;
            }

            public double[,] AcceptableError {
                get {
                    return new double[,] { { 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3 } };
                }
            }

            /// <summary>
            /// returns spatial dimension
            /// </summary>
            public int SpatialDimension {
                get;
                private set;
            }

            public int NoOfLevelsets => 1;
        }

        /// <summary>
        /// Basic test case for advecting the level-set field
        /// in a constant velocity field
        /// </summary>
        public class LevelSetSwirlingFlowTest : ILevelSetTest {

            protected double T; // how much distortion, later 8
            protected double sign;

            /// <summary>
            /// ctor
            /// </summary>
            public LevelSetSwirlingFlowTest(int spatDim, int LevelSetDegree, bool reversed, double T = 8.0) {
                this.SpatialDimension = spatDim;
                this.LevelsetPolynomialDegree = LevelSetDegree;
                this.T = T;
                sign = (reversed) ? -1 : 1;
            }

            public double dt {
                get {
                    return -1.0;    // will be set in LevelSetTest() according to level set cfl 
                }
            }

            /// <summary>
            /// computes the timestep size according to the level-set CFL condition
            /// </summary>
            /// <param name="Resolution"></param>
            /// <param name="LSdegree"></param>
            /// <returns></returns>
            public double ComputeTimestep(int Resolution, int LSdegree, int AMRlevel, int temporalResolution) {
                int gridCells1D = (25 * Resolution) * (AMRlevel + 1);
                double h = 1.0 * 1.0 / (double)gridCells1D;
                double dt = h / (1.0) * 1.0; // this is grid width divided by maximum velocity, times safety factor of 1.0 (1.0 no saftey)
                dt /= (double)(LSdegree * LSdegree);

                int timesteps = temporalResolution * Math.Max((int)Math.Ceiling(T / dt), 1); // make sure that the singular point in time is exactly on one timestep

                return (T / (double)timesteps);
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public double getEndTime() {
                return T;
            }

            /// <summary>
            /// creates a square in 2D
            /// </summary>
            /// <param name="Resolution"></param>
            /// <returns></returns>
            public virtual GridCommons CreateGrid(int Resolution) {
                if (Resolution < 1)
                    throw new ArgumentException();

                var xNodes = GenericBlas.Linspace(0, 1, 25 * Resolution + 1);
                var yNodes = GenericBlas.Linspace(0, 1, 25 * Resolution + 1);

                GridCommons grd;
                switch (SpatialDimension) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                        break;
                    case 3:
                    default:
                        throw new ArgumentOutOfRangeException();
                }

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(delegate (double[] X) {
                    return 1;
                });
                
                return grd;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public virtual IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
                var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

                config.Add("wall", new AppControl.BoundaryValueCollection());

                return config;
            }

            /// <summary>
            /// Level-Set movement.
            /// </summary>
            public Func<double[], double, double>[] GetPhi() {
                return new Func<double[], double, double>[] { delegate (double[] X, double time) {

                double Radius = 0.15;

                Func<double[], double, double> PhiFunc;

                switch (SpatialDimension) {
                    case 2: {
                        PhiFunc = new Func<double[], double, double>((X, t) => Math.Sqrt(Math.Pow(X[0] - 0.5, 2) + Math.Pow(X[1] - 0.75, 2)) - Radius);
                        break;
                    }
                    case 3:
                    default:
                    throw new ArgumentOutOfRangeException();
                }
                return sign * PhiFunc(X,T);

            }};
            }

            /// <summary>
            /// velocity: solid body rotation
            /// </summary>
            public Func<double[], double, double>[][] GetU() {
                var Ret = new Func<double[], double, double>[this.NoOfLevelsets][];
                switch (SpatialDimension) {
                    case 2:
                        Ret[0] = new Func<double[], double, double>[] { (X, t) => -2 * Math.Sin(Math.PI * X[1]) * Math.Cos(Math.PI * X[1]) * Math.Pow(Math.Sin(Math.PI * X[0]) , 2) * Math.Cos(Math.PI * t / T), (X, t) => 2 * Math.Sin(Math.PI * X[0]) * Math.Cos(Math.PI * X[0]) * Math.Pow(Math.Sin(Math.PI * X[1]), 2) * Math.Cos(Math.PI * t / T) };
                        break;
                    case 3:
                    default:
                        throw new ArgumentOutOfRangeException();
                }
                return Ret;
            }

            public int LevelsetPolynomialDegree {
                get;
                private set;
            }

            public double[,] AcceptableError {
                get {
                    return new double[,] { { 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3 } };
                }
            }

            /// <summary>
            /// returns spatial dimension
            /// </summary>
            public int SpatialDimension {
                get;
                private set;
            }

            public int NoOfLevelsets => 1;
        }        

        /// <summary>
        /// No advection, just the projection of a circle with some high order even degree polynomial (or signed distance) function
        /// </summary>
        public class LevelSetCircleProjectionTest : ILevelSetTest {

            protected double T = 0.0;
            protected double sign;
            protected bool signed;


            /// <summary>
            /// ctor
            /// </summary>
            public LevelSetCircleProjectionTest(int spatDim, int LevelSetDegree, bool reversed, bool signed) {
                this.SpatialDimension = spatDim;
                this.LevelsetPolynomialDegree = LevelSetDegree;
                sign = (reversed) ? -1 : 1;
                this.signed = signed;
            }

            public double dt {
                get {
                    return -1.0;    // will be set in LevelSetTest() according to level set cfl 
                }
            }

            /// <summary>
            /// computes the timestep size according to the level-set CFL condition
            /// </summary>
            /// <param name="Resolution"></param>
            /// <param name="LSdegree"></param>
            /// <returns></returns>
            public double ComputeTimestep(int Resolution, int LSdegree, int AMRlevel, int temporalResolution) {
                int gridCells1D = (9 * Resolution) * (AMRlevel + 1);
                double h = 1.0 * 1.0 / (double)gridCells1D;
                double dt = h / 1.0; // this is grid width divided by maximum velocity
                dt /= (double)(LSdegree * LSdegree);

                int timesteps = temporalResolution * Math.Max((int)Math.Ceiling(T / dt), 1); // make sure that the singular point in time is exactly on one timestep

                return (T / (double)timesteps);
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public double getEndTime() {
                return T;
            }

            /// <summary>
            /// creates a square in 2D
            /// </summary>
            /// <param name="Resolution"></param>
            /// <returns></returns>
            public virtual GridCommons CreateGrid(int Resolution) {
                if (Resolution < 1)
                    throw new ArgumentException();

                var xNodes = GenericBlas.Linspace(0, 1, 9 * Resolution + 1);
                var yNodes = GenericBlas.Linspace(0, 1, 9 * Resolution + 1);

                GridCommons grd;
                switch (SpatialDimension) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                        break;
                    case 3:
                    default:
                        throw new ArgumentOutOfRangeException();
                }

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(delegate (double[] X) {
                    return 1;
                });

                return grd;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public virtual IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
                var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

                config.Add("wall", new AppControl.BoundaryValueCollection());

                return config;
            }

            int implicitFunctionDegree = 6;
            /// <summary>
            /// Level-Set movement.
            /// </summary>
            public Func<double[], double, double>[] GetPhi() {
                return new Func<double[], double, double>[] { delegate (double[] X, double time) {

                double Radius = 0.25;

                Func<double[], double, double> PhiFunc;

                switch (SpatialDimension) {
                    case 2: {
                                if(signed)
                                    PhiFunc = new Func<double[], double, double>((X, t) => Math.Sqrt(Math.Pow(X[0] - 0.5, 2) + Math.Pow(X[1] - 0.5, 2)) - Radius);
                                else
                                    PhiFunc = new Func<double[], double, double>((X, t) => Math.Pow(Math.Sqrt(Math.Pow(X[0] - 0.5, 2) + Math.Pow(X[1] - 0.5, 2)), implicitFunctionDegree) - Math.Pow(Radius, implicitFunctionDegree));
                        break;
                    }
                    case 3:
                    default:
                    throw new ArgumentOutOfRangeException();
                }
                return sign * PhiFunc(X,T);

            }};
            }

            /// <summary>
            /// velocity: solid body rotation
            /// </summary>
            public Func<double[], double, double>[][] GetU() {
                var Ret = new Func<double[], double, double>[this.NoOfLevelsets][];
                switch (SpatialDimension) {
                    case 2:
                        Ret[0] = new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 };
                        break;
                    case 3:
                    default:
                        throw new ArgumentOutOfRangeException();
                }
                return Ret;
            }

            public int LevelsetPolynomialDegree {
                get;
                private set;
            }

            public double[,] AcceptableError {
                get {
                    return new double[,] { { 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3 } };
                }
            }

            /// <summary>
            /// returns spatial dimension
            /// </summary>
            public int SpatialDimension {
                get;
                private set;
            }

            public int NoOfLevelsets => 1;
        }

        /// <summary>
        /// No advection, just the projection of a cube with a signed distance function
        /// </summary>
        public class LevelSetCubeProjectionTest : ILevelSetTest {

            protected double T = 0.0;
            protected double sign;


            /// <summary>
            /// ctor
            /// </summary>
            public LevelSetCubeProjectionTest(int spatDim, int LevelSetDegree, bool reversed) {
                this.SpatialDimension = spatDim;
                this.LevelsetPolynomialDegree = LevelSetDegree;
                sign = (reversed) ? -1 : 1;
            }

            public double dt {
                get {
                    return -1.0;    // will be set in LevelSetTest() according to level set cfl 
                }
            }

            /// <summary>
            /// computes the timestep size according to the level-set CFL condition
            /// </summary>
            /// <param name="Resolution"></param>
            /// <param name="LSdegree"></param>
            /// <returns></returns>
            public double ComputeTimestep(int Resolution, int LSdegree, int AMRlevel, int temporalResolution) {
                int gridCells1D = (9 * Resolution) * (AMRlevel + 1);
                double h = 1.0 * 1.0 / (double)gridCells1D;
                double dt = h / 1.0; // this is grid width divided by maximum velocity
                dt /= (double)(LSdegree * LSdegree);

                int timesteps = temporalResolution * Math.Max((int)Math.Ceiling(T / dt), 1); // make sure that the singular point in time is exactly on one timestep

                return (T / (double)timesteps);
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public double getEndTime() {
                return T;
            }

            /// <summary>
            /// creates a square in 2D
            /// </summary>
            /// <param name="Resolution"></param>
            /// <returns></returns>
            public virtual GridCommons CreateGrid(int Resolution) {
                if (Resolution < 1)
                    throw new ArgumentException();

                var xNodes = GenericBlas.Linspace(0, 1, 9 * Resolution + 1);
                var yNodes = GenericBlas.Linspace(0, 1, 9 * Resolution + 1);

                GridCommons grd;
                switch (SpatialDimension) {
                    case 2:
                        grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                        break;
                    case 3:
                        var zNodes = GenericBlas.Linspace(0, 1, 9 * Resolution + 1);
                        grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
                        break;
                    default:
                        throw new ArgumentOutOfRangeException();
                }

                grd.EdgeTagNames.Add(1, "wall");

                grd.DefineEdgeTags(delegate (double[] X) {
                    return 1;
                });

                return grd;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public virtual IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
                var config = new Dictionary<string, AppControl.BoundaryValueCollection>();

                config.Add("wall", new AppControl.BoundaryValueCollection());

                return config;
            }

            /// <summary>
            /// Level-Set movement.
            /// </summary>
            public Func<double[], double, double>[] GetPhi() {
                return new Func<double[], double, double>[] { delegate (double[] X, double time) {

                double Length = 0.25;

                Func<double[], double, double> PhiFunc;
                double x = Math.Abs(X[0] - 0.5) - Length;
                double y = Math.Abs(X[1] - 0.5) - Length;

                switch (SpatialDimension) {
                    case 2: {                                
                        PhiFunc = new Func<double[], double, double>((X, t) => Math.Sqrt(Math.Pow(Math.Max(x,0.0),2.0) + Math.Pow(Math.Max(y,0.0),2.0)) + Math.Min(Math.Max(x,y), 0.0));
                        break;
                    }
                    case 3:
                        double z = Math.Abs(X[2] - 0.5) - Length;
                        PhiFunc = new Func<double[], double, double>((X, t) => Math.Sqrt(Math.Pow(Math.Max(x,0.0),2.0) + Math.Pow(Math.Max(y,0.0),2.0) + Math.Pow(Math.Max(z,0.0),2.0)) + Math.Min(Math.Max(z,Math.Max(x,y)), 0.0));
                        break;
                    default:
                    throw new ArgumentOutOfRangeException();
                }
                return sign * PhiFunc(X,T);

            }};
            }

            /// <summary>
            /// velocity: solid body rotation
            /// </summary>
            public Func<double[], double, double>[][] GetU() {
                var Ret = new Func<double[], double, double>[this.NoOfLevelsets][];
                switch (SpatialDimension) {
                    case 2:
                        Ret[0] = new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0 };
                        break;
                    case 3:
                        Ret[0] = new Func<double[], double, double>[] { (X, t) => 0.0, (X, t) => 0.0, (X, t) => 0.0 };
                        break;
                    default:
                        throw new ArgumentOutOfRangeException();
                }
                return Ret;
            }

            public int LevelsetPolynomialDegree {
                get;
                private set;
            }

            public double[,] AcceptableError {
                get {
                    return new double[,] { { 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3 } };
                }
            }

            /// <summary>
            /// returns spatial dimension
            /// </summary>
            public int SpatialDimension {
                get;
                private set;
            }

            public int NoOfLevelsets => 1;
        }
    }
}

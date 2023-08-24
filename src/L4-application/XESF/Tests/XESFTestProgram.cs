using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using XESF.Fluxes;
using ilPSP.LinSolvers;
using ApplicationWithIDT;
using ApplicationWithIDT.OptiLevelSets;
using MathNet.Numerics.Interpolation;
using System.Linq;
using NUnit.Framework;

namespace XESF.Tests {
    public static class XESFTestProgram {
        #region NUnit stuff
        //[OneTimeSetUp]
        //static public void Init() {
        //    BoSSS.Solution.Application.InitMPI();
        //}

        //[OneTimeTearDown]
        //static public void Cleanup() {
        //}
        #endregion

        #region SupersonicWedgeFlow using one LS on a Triangle Mesh (here the wedge is "cut out")
        public static void XDG_SWF_OneLS_Cutout() {
            using(var p = new XESFMain()) {
                var C = XESFHardCodedControl.XDGWedgeFlow_OneLs_Base_GridFromDB(
                    dbPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\WedgeFlow\WedgeFlow_db",
                    lsDegree: 1,
                    shocksetup: ApplicationWithIDT.GetLevelSet.FromFunction,
                    initialValue: ApplicationWithIDT.GetInitialValue.FromFunctionPerSpecies,
                    MaxIterations: 30,
                    dgDegree: 0,
                    initialAngle_shockLS: 42,
                    PlotInterval: 1,
                    interfaceFluxLS1: ConvectiveInterfaceFluxes.GodunovInterface,
                    
                    optiLevelSetType: OptiLevelSetType.GlobalLevelSet,
                    FluxVersion: Fluxes.FluxVersion.Optimized,
                    bulkFlux: ConvectiveBulkFluxes.OptimizedHLLC,
                    iVFromShockRelations: false,
                    meshPath: @"..\..\..\Meshes\WFLastTry5.msh",
                    agg: 0.2,
                    globalization: ApplicationWithIDT.GlobalizationStrategy.LineSearch
                    );



                p.Init(C);
                p.RunSolverMode();
            }
        }
        #endregion
        #region SupersonicWedgeFlow using two LS on a Cartesian Mesh
        public static void XDG_SWF_TwoLs() {
            using(var p = new XESFMain()) {
                var C = XESFHardCodedControl.XDGWedgeFlow_TwoLs_Base(
                    optiLSDegree: 1,
                    //dbPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\WedgeFlow\WedgeFlow_db",
                    lsDegree: 1,
                    shocksetup: ApplicationWithIDT.GetLevelSet.FromParams,
                    optiLevelSetType: OptiLevelSetType.GlobalLevelSet,
                    initialValue: ApplicationWithIDT.GetInitialValue.FromFunctionPerSpecies,
                    MaxIterations: 200,
                    dgDegree: 0,
                    numOfCellsX: 15,
                    numOfCellsY: 10,
                    initialAngle_shockLS: 32,
                    PlotInterval: -1,
                    interfaceFluxLS2: ConvectiveInterfaceFluxes.GodunovInterface,
                    bulkFlux: ConvectiveBulkFluxes.OptimizedHLLC,
                    FluxVersion: Fluxes.FluxVersion.Optimized,
                    agg: 0.4,
                    globalization: ApplicationWithIDT.GlobalizationStrategy.LineSearch
                    );



                p.Init(C);
                p.RunSolverMode();
                var tol1 = 1e-02;
                var tol2 = 2.5;
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < tol2 && p.ResidualVector.MPI_L2Norm() < tol1), $"the L2 Error is greater than {tol1} (Residual {p.ResidualVector.MPI_L2Norm()}, Enriched Residual {p.obj_f_vec.MPI_L2Norm()}");

            }
        }
        public static void XDG_SWF_TwoLs_ExactSol() {
            using(var p = new XESFMain()) {
                var C = XESFHardCodedControl.XDGWedgeFlow_TwoLs_Base(
                    dbPath: null,
                    optiLSDegree: 1,
                    lsDegree: 1,
                    MaxIterations: 50,
                    dgDegree: 0,
                    numOfCellsX: 15,
                    numOfCellsY: 10,
                    //initialAngle_shockLS: 39.3139318, //exact angle
                    initialAngle_shockLS: 32, //exact angle
                    PlotInterval: 1,
                    iVFromShockRelations:true,
                    interfaceFluxLS1: ConvectiveInterfaceFluxes.RoeWall,
                    interfaceFluxLS2: ConvectiveInterfaceFluxes.RoeInterface,
                    bulkFlux: ConvectiveBulkFluxes.Roe,
                    FluxVersion: Fluxes.FluxVersion.NonOptimized,
                    agg: 0.2,
                    globalization: ApplicationWithIDT.GlobalizationStrategy.LineSearch
                    );



                p.Init(C);
                p.RunSolverMode();
                var tol1 = 1e-02;
                var tol2 = 2.5;
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < tol2 && p.ResidualVector.MPI_L2Norm() < tol1), $"the L2 Error is greater than {tol1} (Residual {p.ResidualVector.MPI_L2Norm()}, Enriched Residual {p.obj_f_vec.MPI_L2Norm()}");

            }
        }
        #endregion
        #region SupersonicWedgeFlow using one LS on a Cartesian rotated mesh
        public static void XDG_SWF_OneLs_Cart() {
            using(var p = new XESFMain()) {
                var C = XESFHardCodedControl.XDGWedgeFlow_OneLs_Rotation(
                    //dbPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\WedgeFlow\WedgeFlow_db",
                    MaxIterations:100,
                    dgDegree:0,
                    numOfCellsX:10,
                    numOfCellsY:15,
                    lsDegree:1
                    );

                p.Init(C);
                p.RunSolverMode();
                var tol = 1e-07;
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < tol && p.ResidualVector.MPI_L2Norm() < tol), $"the L2 Error is greater than {tol} (Residual {p.ResidualVector.MPI_L2Norm()}, Enriched Residual {p.obj_f_vec.MPI_L2Norm()}");

            }
        }
        /// <summary>
        /// Runs the Supersonic Wedge Flow Configuration with one Level Set (and rotated mesh) starting from the exact solution and should return therefor residuals of 1e-8 immediatly
        /// </summary>
        internal static void XDG_SWF_OneLs_Cart_Exact() {
            var C = XESFHardCodedControl.XDGWedgeFlow_OneLs_Rotation(
                    dbPath: null,
                    MaxIterations: 100,
                    dgDegree: 0,
                    numOfCellsX: 10,
                    numOfCellsY: 10,
                    initialAngle_shockLS: 39.3139318, //exact angle
                    //initialAngle_shockLS: 32, //exact angle
                    PlotInterval: 1,
                    interfaceFluxLS1: ConvectiveInterfaceFluxes.RoeInterface,
                    bulkFlux: ConvectiveBulkFluxes.Roe
                    );
            var p = new XESFMain();
            p.Init(C);
            p.RunSolverMode();
            var tol = 1e-07;
            Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < tol && p.ResidualVector.MPI_L2Norm() < tol), $"the L2 Error is greater than {tol} (Residual {p.ResidualVector.MPI_L2Norm()}, Enriched Residual {p.obj_f_vec.MPI_L2Norm()}");
        }
        #endregion
        /// <summary>
        ///  BowShock with Trivial first guess, does not work so far
        /// </summary>
        public static void XDGBS_Local(int xCells, int yCells, int deg, int sDeg) {
            using(var p = new XESFMain()) {
                var C = XESF.XESFHardCodedControl.XDGBS_Local(numX:xCells,numY:yCells,DegE:deg,DegS:sDeg,plotInterval:1,iflux:1,cflux:1);
                p.Init(C);
                p.RunSolverMode();
            }
        }


            /// <summary>
            /// BowShock, where the initial conditions and LevelSet are loaded from a Database, mostly an AV run
            /// </summary>
            public static void XDGBowShockFromDB(int xCells, int yCells, int deg, int sDeg) {
            using(var p = new XESFMain()) {
                var C = XESFHardCodedControl.XDGBowShock_TwoLs_LSFromDB(
                    optiLSDegree: 3,
                    lsTwoDegree: 3,
                    lsOneDegree: 4,

                    //dbPath: null,
                    //dbPath: @"C:\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_db",
                    dbPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_db",

                    //22x10 BowShock *************
                    ///Uni PC
                    //shockLevelSet_Db: @"C:\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_db",
                    //shockLevelSet_SessionId: @"eacec144-46bb-449d-a7cc-7ba8f6a4e3db",
                    //Home PC
                    //shockLevelSet_Db: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\XESF\BowShock\BowShock_db",
                    //IVtsNumber: 0,
                    //pointPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-mag\XDGShock\Tests\BowShockPoints.txt",
                    //initialValue: GetInitialValue.FromDBXDG,
                    //getLevelSet: GetLevelSet.DirectyFromTimestep,
                    //shockLevelSet_SessionId: @"e939131f-e3a1-4199-a49c-028b0c59dc6b", // <--- after p0 timestepping
                    //shockLevelSet_SessionId: @"2a419621-cf42-4217-996f-045a59bde8c2", // <---- best solution run so far agg=0.2
                    //shockLevelSet_SessionId: @"012530cd-de44-4938-afc9-ba3d67933d2d", // <---- pretty much the same as above but does some p=4 iterations, also agg=0.2
                    //shockLevelSet_SessionId: @"9f07b706-4c79-4290-b297-79ec83ce1fe0", // <---- bad run agg=0.1, stagnates after changing to p=1
                    //shockLevelSet_SessionId: @"5292909a-4416-4f83-bbb8-d28ae20382d3", // <- after p0 Timestepping with NewSpline
                    //shockLevelSet_SessionId: @"eacec144-46bb-449d-a7cc-7ba8f6a4e3db",


                    //MArkus AV Run *************
                    ///Uni PC
                    //shockLevelSet_Db: @"C:\experimental\internal\src\private-mag\XDGShock\Tests\bosss_db_levelSets.zip",
                    //shockLevelSet_SessionId: @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2",
                    //pointPath: @"C:\experimental\internal\src\private-mag\XDGShock\Tests\BowShockPoints.txt",
                    //initialValue: GetInitialValue.FromDBSinglePhase,
                    ///Home PC
                    shockLevelSet_Db: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\Databases\bosss_db_levelSets\bosss_db_levelSets",
                    shockLevelSet_SessionId: @"9c45ebf9-f3e0-4d1d-bf91-776bf46e4fc2",
                    pointPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-mag\XDGShock\Tests\BowShockPoints.txt",
                    initialValue: GetInitialValue.FromDBSinglePhase,

                    MaxIterations: 200,
                    dgDegreeStart:sDeg,
                    dgDegreeEnd:deg,
                    MinPIter:new int[] {20,20,20,20,20,20},
                    agg: 0.4,
                    numOfCellsX: xCells,
                    numOfCellsY: yCells,
                    applyReInit: true,
                    ReInitTols:new double[] { -1.5, -1.5, -1.5, -1.5, -1.5 },
                    //solverRunType: SolverRunType.Standard,
                    solverRunType: SolverRunType.Staggerd,
                    staggeredTS: new int[] { 20, 35, 45, 55 ,65},
                    //numOfCellsX: 5,
                    //numOfCellsY: 10,
                    PlotInterval: 1,
                    bulkFlux: ConvectiveBulkFluxes.OptimizedHLLC,
                    interfaceFluxLS1: ConvectiveInterfaceFluxes.OptimizedHLLCWall_Separate_For_Each_Var,
                    interfaceFluxLS2: ConvectiveInterfaceFluxes.GodunovInterface,
                    FluxVersion: XESF.Fluxes.FluxVersion.Optimized
                ) ;
                p.Init(C);
                p.RunSolverMode();
                SaveIsoContourToTextFile(p, $"BowShock_{C.SolDegree}_10x22.txt");
            }
        }

        /// <summary>
        /// Helper function to store interpolating points of an SplineLevelSet
        /// </summary>
        /// <param name="p"></param>
        /// <param name="filename"></param>
        public static void SaveIsoContourToTextFile(XESFMain p, string filename) {
            if(p.LevelSetOpti is SplineOptiLevelSet spliny) {
                spliny.GetSpline();
                if(spliny.Spline is CubicSpline cSpliny) {
                    var yMax = ((GridData) p.Grid.iGridData).Vertices.Coordinates.ExtractSubArrayShallow(-1, 1).To1DArray().Max();
                    var yMin = ((GridData)p.Grid.iGridData).Vertices.Coordinates.ExtractSubArrayShallow(-1, 1).To1DArray().Min();
                    var yPoints = GenericBlas.Linspace(yMin, yMax, 200);
                    var xPoints = new double[100];
                    var cPoints = new double[100];
                    var allPoints = MultidimensionalArray.Create(100, 3);
                    for(int i = 0; i < 100; i++) {
                        allPoints[i, 0] = yPoints[i];

                        xPoints[i] = cSpliny.Interpolate(yPoints[i]);
                        allPoints[i, 1] = xPoints[i];

                        cPoints[i] = cSpliny.Differentiate(yPoints[i]);
                        allPoints[i, 2] = cPoints[i];
                    }
                    allPoints.SaveToTextFile(filename);

                } else {
                    throw new NotSupportedException("not supported - but one could implemet a Newton Rootfinding along lines of type x -> (x,y_i) for that");
                }
            }
        }
        
    }
}

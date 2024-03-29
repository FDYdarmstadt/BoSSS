using ApplicationWithIDT;
using ApplicationWithIDT.OptiLevelSets;
using ilPSP.Utils;
using NUnit.Framework;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using System;
using BUIDT.Fluxes;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using System.Linq;

namespace BUIDT.Tests
{
    [TestFixture]
    public static class BUIDTTestProgram
    {
        #region NUnit stuff
        [OneTimeSetUp]
        static public void Init() {
            BoSSS.Solution.Application.InitMPI();
        }

        [OneTimeTearDown]
        static public void Cleanup() {
        }
        #endregion

        [Test]
        /// <summary>
        /// Test for the StraightShock presented at Eccomas22
        /// </summary>
        public static void StraightShockCurvedStart_Eccomas22()
        {
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            using (var p = new BUIDTMain())
            {
                var C = BUIDTHardCodedControl.StraightShockCurvedStart_Eccomas22(
                    //dbPath: @"C:\Users\jakob\Documents\Uni\Promotion\Programmieren\BoSSS\experimental\internal\src\private-seb\Notebooks\BUIDT\AcceleratingShock\BUIDT_db",
                    dbPath: null,
                    MaxIterations: 50,
                    dgDegree: 0,
                    numOfCellsX: 10,
                    numOfCellsY: 10,
                    agg:0.2,
                    ImmediatePlotPeriod: -1
                    );
                p.Init(C);
                p.RunSolverMode();
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 1e-09 && p.ResidualVector.MPI_L2Norm() < 1e-09), System.String.Format("the L2 Error is greater than 1e-10 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(),p.obj_f_vec.MPI_L2Norm()));

            }
        }

        [Test]
        /// <summary>
        /// Test for the AcceleratingShock
        /// </summary>
        public static void AcceleratingShock() {
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            using (var p = new BUIDTMain()) {
                var C = BUIDTHardCodedControl.AccShock(
                    dbPath:null,
                    MaxIterations: 100,
                    dgDegree: 3,
                    numOfCellsX: 10,
                    numOfCellsY: 10,
                    OptiNumOfCellsX: 10,
                    OptiNumOfCellsY: 10,
                    linearization: Linearization.FD,
                    agg: 0.4,
                    ImmediatePlotPeriod: -1,
                    optiLevelSetType: OptiLevelSetType.SplineLevelSet,
                    getLevelSet: GetLevelSet.FromFunction,
                    solverRunType: SolverRunType.PContinuation
                    ) ;
                p.Init(C);
                p.RunSolverMode();
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 5e-01 && p.ResidualVector.MPI_L2Norm() < 5e-02), System.String.Format("the L2 Error is greater than 1e-02,1e-01 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(), p.obj_f_vec.MPI_L2Norm()));

            }
        }
        /// <summary>
        /// Utility function to save the Explicit surface described by a spline LevelSet
        /// </summary>
        /// <param name="p"></param>
        /// <param name="filename"></param>
        public static void SaveIsoContourToTextFile(BUIDTMain p, string filename) {
            if(p.LevelSetOpti is SplineOptiLevelSet spliny) {
                spliny.GetSpline();
                if(spliny.Spline is CubicSpline cSpliny) {

                    var yPoints = GenericBlas.Linspace(0, 1.2, 100);
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

                }
            }
        }
    }

}

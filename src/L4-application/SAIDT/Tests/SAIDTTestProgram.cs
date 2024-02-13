using ApplicationWithIDT;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SAIDT.Tests {
    [TestFixture]
    public static class SAIDTTestProgram {
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
        //This test checks if the results presented on Eccomas 2022 can still be reproduced
        public static void CurvedShock_Eccomas22() {
            using(var p = new SAIDTMain()) {
                var C = SAIDTHardCodedControl.CurvedShock_Eccomas22(
                        dbPath: null,
                        MaxIterations: 50,
                        dgDegree: 0,
                        numOfCellsX: 10,
                        numOfCellsY: 10,
                        agg: 0.4,
                        optiLevelSetType: OptiLevelSetType.SplineLevelSet
                        );
                C.MinPIter = new int[] {40};
                p.Init(C);
                p.RunSolverMode();
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 1e-09 && p.ResidualVector.MPI_L2Norm() < 1e-09), System.String.Format("the L2 Error is greater than 1e-09 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(), p.obj_f_vec.MPI_L2Norm()));
            }
        }

        [Test]
        //Some Example with good first guess, but p=2
        public static void StraightShock_p2() {
            using(var p = new SAIDTMain()) {
                var C = SAIDTHardCodedControl.StraightShock(
                    dbPath: null,
                    MaxIterations: 50,
                    dgDegree: 2,
                    numOfCellsX: 5,
                    numOfCellsY: 5,
                    OptiNumOfCellsX: 1,
                    OptiNumOfCellsY: 1,
                    agg: 0.1,
                    ImmediatePlotPeriod: -1,
                    optiLevelSetType: OptiLevelSetType.GlobalLevelSet,
                    LSDegree: 1,
                    withReInit: true
                    );
                //SAIDTMain.DeleteOldPlotFiles();
                p.Init(C);
                p.RunSolverMode();
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 1e-09 && p.ResidualVector.MPI_L2Norm() < 1e-09), System.String.Format("the L2 Error is greater than 1e-09 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(), p.obj_f_vec.MPI_L2Norm()));
            }
        }        //Some Example with good first guess, but p=3
        [Test]
        //Some straight Example with good first guess and SinglePhaseFieldLevelSet
        // very bad convergence to the actual solution
        public static void StraightShock_p0_SInglePhaseFieldLS() {
            using(var p = new SAIDTMain()) {
                var C = SAIDTHardCodedControl.StraightShock(
                    dbPath: null,
                    MaxIterations: 50,
                    dgDegree: 0,
                    numOfCellsX: 5,
                    numOfCellsY: 5,
                    OptiNumOfCellsX: 5,
                    OptiNumOfCellsY: 5,
                    agg: 0.1,
                    ImmediatePlotPeriod: -1,
                    optiLevelSetType: OptiLevelSetType.SinglePhaseField,
                    LSDegree: 1,
                    withReInit: true
                    );
                C.Gamma_Min = 1e-2;
                //SAIDTMain.DeleteOldPlotFiles();
                p.Init(C);
                p.RunSolverMode();
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 1e-01 && p.ResidualVector.MPI_L2Norm() < 1e-01), System.String.Format("the L2 Error is greater than 1e-09 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(), p.obj_f_vec.MPI_L2Norm()));
            }
        }
        [Test]
        //Some straight Example with good first guess and SinglePhaseFieldLevelSet
        public static void StraightShock_p0_SplineLevelSet() {
            using(var p = new SAIDTMain()) {
                var C = SAIDTHardCodedControl.StraightShock(
                    MaxIterations: 150,
                    dgDegree: 0,
                    numOfCellsX: 5,
                    numOfCellsY: 5,
                    OptiNumOfCellsX: 5,
                    OptiNumOfCellsY: 5,
                    agg: 0.1,
                    ImmediatePlotPeriod: -1,
                    optiLevelSetType: OptiLevelSetType.SplineLevelSet,
                    LSDegree: 3,
                    withReInit: true,
                    isFarConfig: true
                    ) ;
                C.MinPIter = new int[] { 100, 100, 100 }; //we need to ensure enough iterations
                C.Gamma_Start = 1e-1;
                C.Gamma_Min = 1e-2;
                p.Init(C);
                p.RunSolverMode();
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 1e-10 && p.ResidualVector.MPI_L2Norm() < 1e-13), System.String.Format("the L2 Error is greater than 1e-09 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(), p.obj_f_vec.MPI_L2Norm()));
            }
        }

    }

}

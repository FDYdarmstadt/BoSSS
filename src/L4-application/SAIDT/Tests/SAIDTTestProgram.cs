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
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            using(var p = new SAIDTMain()) {
                var C = SAIDTHardCodedControl.CurvedShock_Eccomas22(
                        dbPath: null,
                        MaxIterations: 50,
                        dgDegree: 0,
                        numOfCellsX: 10,
                        numOfCellsY: 10,
                        ImmediatePlotPeriod:-1,
                        agg: 0.4,
                        optiLevelSetType: OptiLevelSetType.SplineLevelSet
                        );
                C.minimalSQPIterations = new int[] {40};
                p.Init(C);
                p.RunSolverMode();
                p.Gammas.SaveToTextFile("Gammas.txt");
                p.Alphas.SaveToTextFile("Alphas.txt");
                p.Residuals.SaveToTextFile("Residuals");
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 3e-08 && p.ResidualVector.MPI_L2Norm() < 3e-08), System.String.Format("the L2 Error is greater than 1e-09 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(), p.obj_f_vec.MPI_L2Norm()));
            }
        }     

        [Test]
        //Some straight Example with good first guess and SinglePhaseFieldLevelSet
        public static void StraightShock_p0_SplineLevelSet() {
            BoSSS.Solution.Application.InitMPI();

            BoSSS.Solution.Application.DeleteOldPlotFiles();
            using (var p = new SAIDTMain()) {
                var C = SAIDTHardCodedControl.StraightShock(
                    MaxIterations: 150,
                    dgDegree: 0,
                    numOfCellsX:10,
                    numOfCellsY: 10,
                    agg: 0.2,
                    ImmediatePlotPeriod: 1,
                    optiLevelSetType: OptiLevelSetType.SplineLevelSet,
                    LSDegree: 3,
                    withReInit: true,
                    isFarConfig: true
                    ) ;
                C.minimalSQPIterations = new int[] { 100, 100, 100 }; //we need to ensure enough iterations
                C.Gamma_Start = 1;
                p.Init(C);
                p.RunSolverMode();
                Assert.IsTrue((p.obj_f_vec.MPI_L2Norm() < 1e-4 && p.ResidualVector.MPI_L2Norm() < 1e-4), System.String.Format("the L2 Error is greater than 1e-09 (Residual {0}, Enriched Residual {1}", p.ResidualVector.MPI_L2Norm(), p.obj_f_vec.MPI_L2Norm()));
            }
        }

    }

}

/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.Rheology
{
    [TestFixture]
    static class RheologyTestProgram
    {

        [TestFixtureSetUp]
        public static void Init()
        {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);


        }

        [TestFixtureTearDown]
        public static void Cleanup()
        {
            //Console.Out.Dispose();
            //csMPI.Raw.mpiFinalize();
        }

        //TESTS_CHANNEL__________________________________________________________________________________________________
        //Test 1: Insert exact solution and only compute residual.
        [Test]
        public static void ChannelTestComputeRes(
#if DEBUG
            [Values(3)] int GridRes,
            [Values(2)] int deg,
            [Values(0.0, 1.0)] double beta
#else
            [Values(3)] int GridRes,
            [Values(2)] int deg,
            [Values(0.0, 0.5, 1.0)] double beta
#endif
            )
        {

            var C = RheologyChannelTest.ChannelComputeRes(GridRes, deg, beta);
            double AcceptableRes = 5E-8;
            double AcceptableResError = 5E-8;
            bool checkPressure = true;
            GenericTest(C, AcceptableRes, AcceptableResError, checkPressure);
        }

        //Test 2: Compute Stokes system
        [Test]
        public static void ChannelTestStokes(
#if DEBUG
            [Values(3, 4)] int GridRes,
            [Values(2)] int deg,
            [Values(0.0, 1.0)] double beta
#else
            [Values(3, 4, 5)] int GridRes,
            [Values(2, 3)] int deg,
            [Values(0.0, 0.5, 1.0)] double beta
#endif
            )
        {

            var C = RheologyChannelTest.ChannelStokes(GridRes, deg, beta);

            double AcceptableRes = 5E-8;
            double AcceptableResError = 5E-8;
            bool checkPressure = true;
            GenericTest(C, AcceptableRes, AcceptableResError, checkPressure);
        }

     
        //TESTS_CHANNEL_END_______________________________________________________________________________________________________________

        //TESTS_LOCAL_DG__________________________________________________________________________________________________________________
        
        //Test 1: Compute Stokes system
        [Test]
        public static void LocalDGTestStokes(
#if DEBUG
            [Values(4)] int GridRes,
            [Values(2)] int deg,
            [Values(0.0, 1.0)] double beta
#else
            [Values(3, 4, 5)] int GridRes,
            [Values(1, 2, 3)] int deg,
            [Values(0.0, 0.5, 1.0)] double beta
#endif
            )
        {

            var C = RheologyLocalDGTest.LocalDGStokes(GridRes, deg, beta);

            double AcceptableRes = 5E-8;
            double AcceptableResError = 5;
            bool checkPressure = true;
            GenericTest(C, AcceptableRes, AcceptableResError, checkPressure);
        }

        //Test 2: Convergence order test
        [Test]
        static public void LocalDGConvergence([Range(1, 3)] int order)
        {

            RheologyControl fine, coarse;

            coarse = RheologyLocalDGTest.LocalDGStokes(3, order, 0); //(GridRes, Deg, Beta)

            fine = RheologyLocalDGTest.LocalDGStokes(6, order, 0); //(GridRes, Deg, Beta)

            TestConvergence(order, coarse, fine);

        }
        //TESTS_LOCAL_DG_END____________________________________________________________________________________________________________________________

        //TESTS_CONSISTENCY_CONSTITUTIVE________________________________________________________________________________________________________________
        //Test 1: Insert exact polynomial solution and only compute residual.
        [Test]
        public static void ConsistencyConstitutiveTestComputeRes(
#if DEBUG
            [Values(3)] int GridRes,
            [Values(2)] int deg,
            [Values(0)] double beta
#else
            [Values(3)] int GridRes,
            [Values(2)] int deg,
            [Values(0)] double beta
#endif
            )
        {

            var C = RheologyConsistencyTest.ConsistencyConstitutiveComputeRes(GridRes, deg, beta);
            double AcceptableRes = 1e-12;
            double AcceptableResError = 1e-12;
            bool checkPressure = false;
            GenericTest(C, AcceptableRes, AcceptableResError, checkPressure);
        }

        //Test 2: Insert exact polynomial solution and test if nonlinear solver (NEWTON) stays stable.
        [Test]
        public static void ConsistencyConstitutiveTestStability(
#if DEBUG
            [Values(3)] int GridRes,
            [Values(2)] int deg,
            [Values(0)] double beta
#else
            [Values(3)] int GridRes,
            [Values(2)] int deg,
            [Values(0, 0.5, 1.0)] double beta
#endif
            )
        {

            var C = RheologyConsistencyTest.ConsistencyConstitutiveStability(GridRes, deg, beta);
            double AcceptableRes = 10;
            double AcceptableResError = 0.125;
            bool checkPressure = false;
            GenericTest(C, AcceptableRes, AcceptableResError, checkPressure);
        }

        //TESTS_CONSISTENCY_CONSTITUTIVE_END________________________________________________________________________________________________________

        private static void GenericTest(RheologyControl C, double AcceptableRes, double AcceptableResError, bool checkPressure)
        {
            using (var solver = new Rheology())
            {
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 3;
                Solution.Application.DeleteOldPlotFiles();

                solver.Init(C);
                solver.RunSolverMode();

                // matrix analysis
                // ===============

                BlockMsrMatrix OpMatrix;
                double[] OpAffine;
                solver.AssembleMatrix(out OpMatrix, out OpAffine, solver.CurrentSolution.Mapping.ToArray(), true);
                if (OpMatrix.NoOfRows < 1000)
                {
                    //int rnk = (int)Math.Round(OpMatrix.rank());
                    //Assert.IsTrue(OpMatrix.NoOfRows == rnk, "Matrix is degenerated: Rank out of range.");

                    //double Cond = OpMatrix.cond();
                    //Assert.LessOrEqual(Cond, 1.0e8, "Matrix is degenerated: Condition no out of range.");
                }


                // check residuals and errors
                // ==========================

                double ResThresh = AcceptableRes;
                double[] ResNorms = new double[solver.CurrentSolution.Mapping.Fields.Count()];
                for (int i = 0; i < ResNorms.Length; i++)
                {

                        ResNorms[i] = solver.CurrentResidual.Mapping.Fields[i].L2Norm();
                        Console.WriteLine("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Mapping.Fields[i].Identification, ResNorms[i]);
                }

                if (checkPressure == true)
                {
                    for (int i = 0; i < ResNorms.Length; i++)
                        Assert.LessOrEqual(ResNorms[i], ResThresh, "L2 norm of " + solver.CurrentResidual.Mapping.Fields[i].Identification + " too high.");
                }
                else
                {
                    for (int i = 0; i < 2; i++)
                        Assert.LessOrEqual(ResNorms[i], ResThresh, "L2 norm of " + solver.CurrentResidual.Mapping.Fields[i].Identification + " too high.");
                    for (int i = 3; i < ResNorms.Length; i++)
                        Assert.LessOrEqual(ResNorms[i], ResThresh, "L2 norm of " + solver.CurrentResidual.Mapping.Fields[i].Identification + " too high.");
                }

                double ResThresh2 = AcceptableResError;

                double L2VelX = (double)solver.QueryHandler.QueryResults["L2err_VelocityX"];
                Assert.LessOrEqual(L2VelX, ResThresh2, "L2 Error of VelocityX seems wrong");
                Console.WriteLine("L2err_VelocityX is" + L2VelX);

                double L2VelY = (double)solver.QueryHandler.QueryResults["L2err_VelocityY"];
                Assert.LessOrEqual(L2VelX, ResThresh2, "L2 Error of VelocityY seems wrong");
                Console.WriteLine("L2err_VelocityY is" + L2VelY);

                double L2Pres = (double)solver.QueryHandler.QueryResults["L2err_Pressure"];
                Assert.LessOrEqual(L2VelX, ResThresh2, "L2 Error of Pressure seems wrong");
                Console.WriteLine("L2err_Pressure is" + L2Pres);

                if (C.beta != 1)
                {

                    double L2StressXX = (double)solver.QueryHandler.QueryResults["L2err_StressXX"];
                    Assert.LessOrEqual(L2VelX, ResThresh2, "L2 Error of StressXX seems wrong");

                    double L2StressXY = (double)solver.QueryHandler.QueryResults["L2err_StressXY"];
                    Assert.LessOrEqual(L2VelX, ResThresh2, "L2 Error of StressXY seems wrong");

                    double L2StressYY = (double)solver.QueryHandler.QueryResults["L2err_StressYY"];
                    Assert.LessOrEqual(L2VelX, ResThresh2, "L2 Error of StressYY seems wrong");

                }
            }
        }

        private static void TestConvergence(int order, RheologyControl coarseRun, RheologyControl fineRun)
        {
            var solver = new Rheology();
            solver.Init(coarseRun);
            solver.RunSolverMode();

            double L2VelXCoarse = (double)solver.QueryHandler.QueryResults["L2err_VelocityX"];
            double L2VelYCoarse = (double)solver.QueryHandler.QueryResults["L2err_VelocityY"];
            double L2PresCoarse = (double)solver.QueryHandler.QueryResults["L2err_Pressure"];
            double L2StressXXCoarse = (double)solver.QueryHandler.QueryResults["L2err_StressXX"];
            double L2StressXYCoarse = (double)solver.QueryHandler.QueryResults["L2err_StressXY"];
            double L2StressYYCoarse = (double)solver.QueryHandler.QueryResults["L2err_StressYY"];

            solver = new Rheology();
            solver.Init(fineRun);
            solver.RunSolverMode();

            double L2VelXFine = (double)solver.QueryHandler.QueryResults["L2err_VelocityX"];
            double L2VelYFine = (double)solver.QueryHandler.QueryResults["L2err_VelocityY"];
            double L2PresFine = (double)solver.QueryHandler.QueryResults["L2err_Pressure"];
            double L2StressXXFine = (double)solver.QueryHandler.QueryResults["L2err_StressXX"];
            double L2StressXYFine = (double)solver.QueryHandler.QueryResults["L2err_StressXY"];
            double L2StressYYFine = (double)solver.QueryHandler.QueryResults["L2err_StressYY"];

            //Debugger.Launch();

            double rateL2VelX = (Math.Log(L2VelXCoarse) - Math.Log(L2VelXFine)) / Math.Log(2);
            double rateL2VelY = (Math.Log(L2VelYCoarse) - Math.Log(L2VelYFine)) / Math.Log(2);
            double rateL2Pres = (Math.Log(L2PresCoarse) - Math.Log(L2PresFine)) / Math.Log(2);
            double rateL2StressXX = (Math.Log(L2StressXXCoarse) - Math.Log(L2StressXXFine)) / Math.Log(2);
            double rateL2StressXY = (Math.Log(L2StressXYCoarse) - Math.Log(L2StressXYFine)) / Math.Log(2);
            double rateL2StressYY = (Math.Log(L2StressYYCoarse) - Math.Log(L2StressYYFine)) / Math.Log(2);

            Console.WriteLine("Convergence rate Vecolity X L2 error: " + rateL2VelX);
            Console.WriteLine("Convergence rate Vecolity Y L2 error: " + rateL2VelY);
            Console.WriteLine("Convergence rate Pressure L2 error: " + rateL2Pres);
            Console.WriteLine("Convergence rate StressXX L2 error: " + rateL2StressXX);
            Console.WriteLine("Convergence rate StressXY L2 error: " + rateL2StressXY);
            Console.WriteLine("Convergence rate StressYY L2 error: " + rateL2StressYY);

            double tolerance = 0.4;

            Assert.That(rateL2VelX >= order + 1 - tolerance,
            String.Format("Convergence rate of L2 error of Velocity X lower than expected (was {0} but should be {1} +- {2})", rateL2VelX, order + 1, tolerance));
            Assert.That(rateL2VelY >= order + 1 - tolerance,
            String.Format("Convergence rate of L2 error of Velocity Y lower than expected (was {0} but should be {1} +- {2})", rateL2VelY, order + 1, tolerance));

            Assert.That(rateL2Pres >= order - tolerance,
            String.Format("Convergence rate of L2 error of Pressure lower than expected (was {0} but should be {1} +- {2})", rateL2Pres, order, tolerance));

            Assert.That(rateL2StressXX >= order + 1 - tolerance,
            String.Format("Convergence rate of L2 error of Stress XX lower than expected (was {0} but should be {1} +- {2})", rateL2StressXX, order + 1, tolerance));
            Assert.That(rateL2StressXY >= order + 1 - tolerance,
            String.Format("Convergence rate of L2 error of Stress XY lower than expected (was {0} but should be {1} +- {2})", rateL2StressXY, order + 1, tolerance));
            Assert.That(rateL2StressYY >= order + 1 - tolerance,
            String.Format("Convergence rate of L2 error of Stress YY lower than expected (was {0} but should be {1} +- {2})", rateL2StressYY, order + 1, tolerance));

        }
    }
}
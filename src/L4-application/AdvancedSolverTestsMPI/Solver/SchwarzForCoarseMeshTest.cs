using AdvancedSolverTests;
using BoSSS.Solution.AdvancedSolvers;
using AdvancedSolverTests.SubBlocking;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AdvancedSolverTestsMPI.Solver {

    [TestFixture]
    public static class SchwarzForCoarseMeshTest {


        // --test=AdvancedSolverTestsMPI.Solver.SchwarzForCoarseMeshTest.TestInit
        [Test]
        public static void TestInit(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(4)] int Res,
            [Values(-2, -1, 0, 1, 2)] int BlockVariation,
            [Values(0, 1, 2)] int SwzOverlap
            ) { //
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if (MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            if (MPISize + BlockVariation <= 0) {
                Console.WriteLine($"Skipping test for {MPISize + BlockVariation} blocks...");
                return;
            }


            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("ExternalIndexTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- Get global index by mapping
            using (var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MatrixShape.laplace, Res)) {
                MultigridOperator MGOp = O.MGOp;

                var s = new SchwarzForCoarseMesh();
                s.config.NoOfBlocks = MPISize + BlockVariation;
                s.config.Overlap = SwzOverlap;
                s.InitWithTest(MGOp, true);
            }
        }
    }
}

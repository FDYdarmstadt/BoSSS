using AdvancedSolverTests;
using BoSSS.Solution.AdvancedSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using MPI.Wrappers;
using ilPSP.LinSolvers;
using System.Diagnostics;
using NUnit.Framework;
using AdvancedSolverTests.SubBlocking;

namespace AdvancedSolverTests.Solver {

    [TestFixture]
    public class addSchwarzTest  {

        internal class TestBlockingStrat : Schwarz.BlockingStrategy {

            public TestBlockingStrat() {
            }

            private int IdentifyBlock(int iCell, MultigridOperator op) {
                //int NoOfNodes = cell.TransformationParams.NoOfRows;
                //int D = op.GridData.SpatialDimension;
                //var center = new double[D];
                //for (int k = 0; k < NoOfNodes; k++) {
                //    for (int d = 0; d < D; d++) {
                //        center[d] += cell.TransformationParams[k, d];
                //    }
                //}
                var center = op.Mapping.AggGrid.ParentGrid.iLogicalCells.GetCenter(iCell);
                if (center.x < 0 && center.y < 0)
                    return 0;
                if (center.x > 0 && center.y < 0)
                    return 1;
                if (center.x > 0 && center.y > 0)
                    return 2;
                if (center.x < 0 && center.y > 0)
                    return 3;
                return -1;
            }

            public override IEnumerable<List<int>> GetBlocking(MultigridOperator op) {
                if (op.Mapping.MpiSize > 1) {
                    //int MPIrank = op.Mapping.MpiRank;
                    //return new List<int>[] { SblockList[MPIrank] };
                    NoOfBlocks = 1;
                    return new List<int>[] { op.Mapping.LocalNoOfBlocks.ForLoop(i => i).ToList() };
                } else {
                    List<int>[] SblockList = NoOfBlocks.ForLoop(s => new List<int>());
                    long N_cells = op.Mapping.TotalNoOfBlocks;
                    for (int iCell = 0; iCell < N_cells; iCell++) {
                        int iBlock = IdentifyBlock(iCell, op);
                        SblockList[iBlock].Add(iCell);
                    }
                    return SblockList;
                }
            }

            private int NoOfBlocks = 4;

            public override int GetNoOfBlocks(MultigridOperator op) {
                return NoOfBlocks;
            }
        }

        internal class TestSchwarz : Schwarz {

            public BlockMsrMatrix[] GetBlockMatrices() {
                return base.blockSolvers.Select(bs => bs.SubsysMatrix).ToArray();
            }

            public BlockMask[] GetBlockMasks() {
                return base.blockSolvers.Select(bs => bs.RestrictionMask).ToArray();
            }

            protected override ISubsystemSolver InitBlockSolver(MultigridOperator op, int iPart, GridAndDegRestriction BlockSolver) {
                var ds = new DirectSolver() {
                    
                };
                ds.config.UseDoublePrecision = true;
                ds.Init(BlockSolver.OperatorRestriction);
                return ds;
            }
        }

        internal static BlockMsrMatrix GetLocalMatrix(BlockMsrMatrix Matrix) {
            int rank;
            MPI.Wrappers.csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out rank);

            // extraction of local Matrix block, ILU will be executed process local
            var LocBlocki0s = new List<long>();
            var LocBlockLen = new List<int>();
            long IdxOffset = Matrix._RowPartitioning.i0;
            for (int i = 0; i < Matrix._RowPartitioning.LocalNoOfBlocks; i++) {
                long iBlock = i + Matrix._RowPartitioning.FirstBlock;
                long i0 = Matrix._RowPartitioning.GetBlockI0(iBlock) - IdxOffset;
                int Len = Matrix._RowPartitioning.GetBlockLen(iBlock);
                LocBlocki0s.Add(i0);
                LocBlockLen.Add(Len);
            }
            long[] RowISrc = Matrix._RowPartitioning.LocalLength.ForLoop(i => i + IdxOffset);

            var part = new BlockPartitioning(Matrix._RowPartitioning.LocalLength, LocBlocki0s, LocBlockLen, csMPI.Raw._COMM.SELF);
            BlockMsrMatrix localMatrix = new BlockMsrMatrix(part);

            Matrix.WriteSubMatrixTo(localMatrix, RowISrc, default(long[]), RowISrc, default(long[]));
            return localMatrix;
        }

        internal static TestSchwarz InitSchwarzTestSolver(MultigridOperator op, bool overlap) {
            
            var solver = new TestSchwarz() {
                //m_BlockingStrategy = ,
                m_BlockingStrategy = new TestBlockingStrat(),
            };
            solver.config.Overlap = overlap ? 1 : 0;
            solver.config.EnableOverlapScaling = overlap;
            solver.Init(op);

            return solver;
        }

        public static MultigridOperator CreateTestMGOperator() {
            var MgoPair = Utils.CreateTestMGOperator(XDGusage.all, 3, MatrixShape.diagonal_spec, 40);
            MgoPair.MGOp.OperatorMatrix.Clear();
            int L = MgoPair.MGOp.Mapping.LocalLength;
            long i0 = MgoPair.MGOp.Mapping.i0;
            var V = CreateRndVector(L, MgoPair.MGOp.Mapping.MpiRank);
            for (int i = 0; i < L; i++) {
                MgoPair.MGOp.OperatorMatrix.SetDiagonalElement(i+i0, V[i]);
            }
            return MgoPair.MGOp;
        }

        private static double[] CreateRndVector(int Length, int Seed) {
            var rnd = new Random(Seed);
            return Length.ForLoop(s => rnd.NextDouble());
        }

        public static void CreateRndSolAndRHS(MultigridOperator op, out double[] X, out double[] B) {
            int locLength = op.Mapping.LocalLength;
            X = CreateRndVector(locLength,op.Mapping.MpiRank);
            B = new double[X.Length];
            op.OperatorMatrix.SpMV(1.0, X, 0, B);
        }

        /// <summary>
        /// Only available for non overlapping blocks.
        /// To get overlapping blocks is complex and is covered by sub blocking tests anyway.
        /// </summary>
        [Test]
        public static void TestBlocksAreCorrect_NonOverlapping() {
            Test_Blocks_AreCorrect(false);
        }

        [Test]
        public static void TestSolutionIsCorrect_NonOverlapping() {
            Test_Solution_IsCorrect(false);
        }

        [Test]
        public static void TestSolutionIsCorrect_Overlapping() {
            Test_Solution_IsCorrect(true);
        }

        private static void Test_Blocks_AreCorrect(bool overlapOn) {
            var op = CreateTestMGOperator();
            var solver = InitSchwarzTestSolver(op, overlap: overlapOn);
            var Sblocks = solver.GetBlockMatrices();
            var localblock = GetLocalMatrix(op.OperatorMatrix);
            Debug.Assert(Sblocks.Length==1);
            Sblocks[0].Acc(-1.0, localblock);
            double infNorm = Sblocks[0].InfNorm();
            Assert.IsTrue(infNorm == 0);
        }

        private static void Test_Solution_IsCorrect(bool overlapOn) {
            // Arrange
            var op = CreateTestMGOperator();
            var solver = InitSchwarzTestSolver(op, overlap: overlapOn);
            var Sblocks = solver.GetBlockMatrices();
            var localblock = GetLocalMatrix(op.OperatorMatrix);
            double[] X, B;
            CreateRndSolAndRHS(op, out X, out B);
            double[] Xcheck = X;
            X = new double[X.Length];

            // Act
            Console.WriteLine("NormL2:" + Xcheck.MPI_L2Norm());
            solver.Solve(X, B);
            Xcheck.AccV(-1.0, X);
            double norm = Xcheck.MPI_L2Norm();
            Console.WriteLine("NormL2:" + norm);
            Assert.IsTrue(norm < 1E-14);
        }



        public static bool RunTest() {

            System.Threading.Thread.Sleep(10000);

            TestSolutionIsCorrect_Overlapping();

            return true;
        }
    }
}

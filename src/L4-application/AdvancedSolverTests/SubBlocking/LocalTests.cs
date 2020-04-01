using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Foundation;
using System.Threading;
using System.Globalization;
using BoSSS.Foundation.XDG;
using AdvancedSolverTests;

namespace AdvancedSolverTests.SubBlocking
{
    [TestFixture]
    static public class LocalTests
    {

        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);
            Thread.CurrentThread.CurrentCulture = CultureInfo.CurrentCulture;
        }

        /// <summary>
        /// MPI shutdown.
        /// </summary>
        [TestFixtureTearDown]
        public static void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        [Test]
        public static void LocalIndexTest(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder) {

            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("LocalLIndexTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- Get global index by mapping
            MultigridOperator MGOp = Utils.CreateTestMGOperator(UseXdg, DGOrder);
            var map = MGOp.Mapping;
            int[] fields = map.NoOfVariables.ForLoop(i => i);
            int[] GlobalIdxMap_loc = map.GetSubvectorIndices(fields);

            //Arrange --- Prepare stuff for mask
            var selector = new SubBlockSelector(map);
            var stw = new Stopwatch();
            stw.Reset();

            //Act --- do the masking to get index lists
            stw.Start();
            var mask = new TestMask(selector,null);
            stw.Stop();
            int[] GlobalIdxMask_loc = mask.Global_IList_LocalCells;

            //Assert --- Idx lists are of same length
            Assert.IsTrue(GlobalIdxMap_loc.Length == GlobalIdxMask_loc.Length);

            //Assert --- Compare map and mask indices
            for (int iLoc = 0; iLoc < GlobalIdxMask_loc.Length; iLoc++) {
                Assert.True(GlobalIdxMap_loc[iLoc] == GlobalIdxMask_loc[iLoc]);
            }
        }

  

        [Test]
        public static void MapConsistencyTest(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder
            ) {

            //if no selection is chosen mapping should be the same as of origin
            //assume: indexing is right 'cause of other tests
            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("MapConsistencyTest({0},{1})", UseXdg, DGOrder);

            //Arrange
            MultigridOperator MGOp = Utils.CreateTestMGOperator(UseXdg, DGOrder);
            var sbs = new SubBlockSelector(MGOp.Mapping);
            var mask = new BlockMask(sbs);
            var stw = new Stopwatch();
            stw.Restart();

            //Act --- Create Mapping from mask
            stw.Start();
            var submatrix = mask.GetSubBlockMatrix(MGOp.OperatorMatrix);
            stw.Stop();
            var rowpart = submatrix._RowPartitioning;
            var colpart = submatrix._ColPartitioning;

            //Assert --- Equal Partition of mask and origin
            Assert.AreEqual(rowpart, colpart);
            Assert.IsTrue(rowpart.IsLocallyEqual(MGOp.Mapping));
        }



        /// <summary>
        /// Test for extracting masked diagonal cell blocks and respective vector matrix multiplication
        /// </summary>
        [Test]
        public static void SubBlockExtractionWithCoupling(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal, MatrixShape.diagonal_var, MatrixShape.diagonal_spec, MatrixShape.diagonal_var_spec)] MatrixShape MShape
            ) {

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("ExtractDiagonalBlocks({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- get multigridoperator
            MultigridOperator MGOp = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape);
            BlockMsrMatrix M = MGOp.OperatorMatrix;
            MultigridMapping map = MGOp.Mapping;

            //Arrange --- setup masking
            SubBlockSelector SBS = new SubBlockSelector(map);
            BlockMask mask = new BlockMask(SBS, null);
            bool[] coupling = Utils.SetCoupling(MShape);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Arrange --- setup auxiliary matrix
            //this will show us if more is extracted, than it should ...
            var Mprep = new BlockMsrMatrix(map);
            Mprep.Acc(1.0, M);

            //Act --- diagonal subblock extraction
            stw.Start();
            var blocks = mask.GetSubBlocks(Mprep, coupling[0], coupling[1], coupling[2]);
            stw.Stop();

            //Assert --- all diagonal blocks are extracted
            Assert.IsTrue(blocks.Length == map.LocalNoOfBlocks);


            for (int i = 0; i < map.LocalNoOfBlocks; i++) {

                //Arrange --- get ith diagonal block of M: M_i
                int iBlock = i + map.AggGrid.CellPartitioning.i0;
                int L = map.GetBlockLen(iBlock);
                int i0 = map.GetBlockI0(iBlock);
                var Mblock = MultidimensionalArray.Create(L, L);
                M.ReadBlock(i0, i0, Mblock);

                //Act --- M_i-Mones_i
                Mblock.Acc(-1.0, blocks[i]);

                //Assert --- are extracted blocks and 
                Assert.IsTrue(Mblock.InfNorm() == 0.0, String.Format("infNorm of block {0} neq 0!", i));
            }



            //BlockMsrMatrix all1;
            //all1.SetAll(1);
            //Generate broken diagonal matrix, die zur Maske passt: M
            //M+all1=M_prep
            //Wende Extraction auf M_prep an, Man sollte nun M bekommen
            //Test: M_prep-extract(M_prep)=all1
            //Test-crit: Result.SumEntries=DOF^2 oder Result.Max()==Result.Min()==1
            //oder (besser)
            //Test: M-extract(M_prep)=zeros
            //Test-crit: Result.InfNorm()==0

            //Der Test kann für ExtractSubMatrix mit ignore coupling wiederholt werden
            //eventuell: Testmatrix finden mit brauchbaren Nebendiagonalen für einen Fall

            //Was wird getestet: funktioniert ignorecoupling richtig?
        }

        [Test]
        public static void SubMatrixExtractionWithCoupling(
        [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
        [Values(2)] int DGOrder,
        [Values(MatrixShape.full, MatrixShape.full_var, MatrixShape.full_spec, MatrixShape.full_var_spec)] MatrixShape MShape
        ) {
            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("ExtractSubMatrixAndIgnoreCoupling({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- get multigridoperator
            MultigridOperator MGOp = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape);
            BlockMsrMatrix M = MGOp.OperatorMatrix;
            MultigridMapping map = MGOp.Mapping;

            //Arrange --- setup masking
            SubBlockSelector SBS = new SubBlockSelector(map);
            BlockMask mask = new BlockMask(SBS, null);
            bool[] coup = Utils.SetCoupling(MShape);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- establish submatrix
            stw.Start();
            var Mext = mask.GetSubBlockMatrix(M, false, coup[1], coup[2]);
            stw.Stop();
            Mext.Acc(-1.0, M);

            //Assert --- Mext conains only diagonal blocks of M
            Assert.IsTrue(Mext.InfNorm() == 0);
        }

        [Test]
        public static void FastSubMatrixExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.laplace)] MatrixShape MShape
            ) {

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("ExtractSubMatrixAndIgnoreCoupling({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- get multigridoperator
            MultigridOperator MGOp = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape);
            BlockMsrMatrix M = MGOp.OperatorMatrix;
            MultigridMapping map = MGOp.Mapping;

            //Arrange --- setup masking
            SubBlockSelector SBS = new SubBlockSelector(map);
            BlockMask mask = new BlockMask(SBS, null);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- establish submatrix
            stw.Start();
            var Mext = mask.GetSubBlockMatrix(M);
            stw.Stop();
            Mext.Acc(-1.0, M);

            //Assert --- Mext conains only diagonal blocks of M
            Assert.IsTrue(Mext.InfNorm() == 0);
        }

        [Test]
        public static void CellBlockVectorOperations(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal, MatrixShape.diagonal_var, MatrixShape.diagonal_spec, MatrixShape.diagonal_var_spec)] MatrixShape MShape
            ) {
            //matrix Erzeugung wie in ExtractDiagonalCellBlocks...
            //Auf der HierarchieEbene, auf der Kopplung ausgesetzt wird kann Auswahl vorgenommen werden
            //bei var: 0 / 1, bei DG: <=1 / >1, bei spec: A / B, bei Cells: odd / even
            //accumulierte Teilergebnisse sind dann == fullM*fullX

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("CellBlockVectorOperations({0},{1},{2})", UseXdg, DGOrder, MShape);

            var mop = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape);
            var map = mop.Mapping;
            double[] Vec = Utils.GetRandomVector(mop.Mapping.LocalLength);

            //Arrange --- setup masking
            SubBlockSelector SBS = new SubBlockSelector(map);
            
            BlockMask mask = new BlockMask(SBS, null);
            bool[] coupling = Utils.SetCoupling(MShape);
            //var blocks = mask.GetSubBlocks(M, coupling[0], coupling[1], coupling[2]);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Assert --- all diagonal blocks are extracted
            //Assert.IsTrue(blocks.Length == map.LocalNoOfBlocks);

            double[] Vec_col = new double[map.LocalLength];

            for (int i = 0; i < map.LocalNoOfBlocks; i++) {
                stw.Start();
                double[] Vec_i = mask.GetVectorCellwise(Vec,i);
                mask.AccVecCellwiseToFull(Vec_i,i, Vec_col);
                stw.Stop();
            }
            Vec_col.AccV(-1.0, Vec);

            //Assert --- are extracted blocks and 
            Assert.IsTrue(Vec_col.L2Norm() == 0.0, String.Format("L2Norm neq 0!"));
        }

        

        [Test]
        public static void SplitVectorOperations(
        [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
        [Values(2)] int DGOrder,
        [Values(MatrixShape.full, MatrixShape.full_spec, MatrixShape.full_var)] MatrixShape MShape
        ) {

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("SplitVectorOperations({0},{1},{2})", UseXdg, DGOrder, MShape);

            //matrix Erzeugung wie in ExtractDiagonalCellBlocks...
            //Auf der HierarchieEbene, auf der Kopplung ausgesetzt wird kann Auswahl vorgenommen werden
            //bei var: 0 / 1, bei DG: <=1 / >1, bei spec: A / B, bei Cells: odd / even
            //accumulierte Teilergebnisse sind dann == fullM*fullX
            var mop = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape);
            var map = mop.Mapping;
            double[] Vec = Utils.GetRandomVector(mop.Mapping.LocalLength);

            //Arrange --- setup masking
            SubBlockSelector sbsA = new SubBlockSelector(map);
            sbsA.SetDefaultSplitSelection(MShape, true);
            BlockMask maskA = new BlockMask(sbsA, null);

            SubBlockSelector sbsB = new SubBlockSelector(map);
            sbsB.SetDefaultSplitSelection(MShape, false);
            BlockMask maskB = new BlockMask(sbsB, null);

            double[] VecAB = new double[Vec.Length];

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- 
            stw.Start();
            var VecA = maskA.GetSubBlockVec(Vec);
            var VecB = maskB.GetSubBlockVec(Vec);

            maskA.AccVecToFull(VecA, VecAB);
            maskB.AccVecToFull(VecB, VecAB);
            stw.Stop();

            Debug.Assert(Vec.L2Norm() != 0);
            double fac = ((MShape == MatrixShape.full_var || MShape == MatrixShape.diagonal_var)  && UseXdg == XDGusage.none) ? -2.0 : -1.0;
            VecAB.AccV(fac, Vec);

            //Assert --- are extracted blocks and 
            Assert.IsTrue(VecAB.L2Norm() == 0.0, String.Format("L2Norm neq 0!"));
        }

        /// <summary>
        /// Test for extracting masked diagonal cell blocks and respective vector matrix multiplication
        /// </summary>
        [Test]
        public static void CellwiseSubSelection(
            [Values(SelectionType.all_combined, SelectionType.degrees, SelectionType.species, SelectionType.variables)] SelectionType SType
            ) {

            Utils.TestInit((int)SType);
            Console.WriteLine("SubSelection({0})", SType);

            //Arrange --- extracts entries of matrix according to hardcoded selection
            int sampleCellA = 0;
            int sampleCellB = 2;
            int DGdegree = 2;
            int GridResolution = 4;
            var mgo = Utils.CreateTestMGOperator(XDGusage.all, DGdegree, MatrixShape.full_var_spec, GridResolution);
            BlockMsrMatrix compA = Utils.GetCellCompMatrix(SType, mgo, sampleCellA);
            BlockMsrMatrix compB = Utils.GetCellCompMatrix(SType, mgo, sampleCellB);

            //Arrange --- setup masking, which correspond to hardcoded
            SubBlockSelector sbsA = new SubBlockSelector(mgo.Mapping);
            sbsA.GetDefaultSelection(SType, sampleCellA); // single spec
            BlockMask maskA = new BlockMask(sbsA, null);
            SubBlockSelector sbsB = new SubBlockSelector(mgo.Mapping);
            sbsB.GetDefaultSelection(SType, sampleCellB); // double spec
            BlockMask maskB = new BlockMask(sbsB, null);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- subblock extraction
            stw.Start();
            var blocksA = maskA.GetSubBlocks(mgo.OperatorMatrix, false, false, false);
            var blocksB = maskB.GetSubBlocks(mgo.OperatorMatrix, false, false, false);
            stw.Stop();

            //Assert --- 
            Assert.IsTrue(blocksA.Length==1);
            Assert.IsTrue(blocksB.Length==1);

            //Assert --- compare masking of single spec cell
            Debug.Assert(compA.InfNorm() != 0.0);
            compA.AccBlock(0, 0, -1.0, blocksA[0]);
            Assert.IsTrue(compA.InfNorm() == 0.0);

            //Assert --- compare masking of double spec cell
            Debug.Assert(compB.InfNorm() != 0.0);
            compB.AccBlock(0, 0, -1.0, blocksB[0]);
            Assert.IsTrue(compB.InfNorm() == 0.0);
        }

        [Test]
        public static void SubSelection(
        [Values(SelectionType.all_combined, SelectionType.degrees, SelectionType.species, SelectionType.variables)] SelectionType SType
        ) {

            Utils.TestInit((int)SType);
            Console.WriteLine("SubSelection({0})", SType);

            //Arrange --- extracts entries of matrix according to hardcoded selection
            int sampleCellA = 0;
            int sampleCellB = 2;
            int DGdegree = 2;
            int GridResolution = 4;
            var mgo = Utils.CreateTestMGOperator(XDGusage.all, DGdegree, MatrixShape.full_var_spec, GridResolution);
            BlockMsrMatrix compA = Utils.GetCellCompMatrix(SType, mgo, sampleCellA);
            BlockMsrMatrix compB = Utils.GetCellCompMatrix(SType, mgo, sampleCellB);

            //Arrange --- setup masking, which correspond to hardcoded
            SubBlockSelector sbsA = new SubBlockSelector(mgo.Mapping);
            sbsA.GetDefaultSelection(SType, sampleCellA); // single spec
            BlockMask maskA = new BlockMask(sbsA, null);
            SubBlockSelector sbsB = new SubBlockSelector(mgo.Mapping);
            sbsB.GetDefaultSelection(SType, sampleCellB); // double spec
            BlockMask maskB = new BlockMask(sbsB, null);

            //Arrange --- stop the watch
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- get subblocks
            stw.Start();
            BlockMsrMatrix subA = maskA.GetSubBlockMatrix(mgo.OperatorMatrix);
            BlockMsrMatrix subB = maskB.GetSubBlockMatrix(mgo.OperatorMatrix);
            stw.Stop();

            //Assert --- compare masking of single spec cell
            Debug.Assert(compA.InfNorm() != 0.0);
            subA.Acc(-1.0, compA);
            Assert.IsTrue(subA.InfNorm() == 0.0);

            //Assert --- compare masking of double spec cell
            Debug.Assert(compB.InfNorm()!=0.0);
            subB.Acc(-1.0, compB);
            Assert.IsTrue(subB.InfNorm() == 0.0);
        }

    }
}

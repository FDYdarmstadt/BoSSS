using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using ilPSP.LinSolvers;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;

namespace AdvancedSolverTests.SubBlocking
{
    [TestFixture]
    class ExternalTests
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
        public static void GetExternalRowsTest(
        [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
        [Values(2)] int DGOrder,
        [Values(4)] int Res) {

            //Matlabaufruf --> gesamte Matrix nach Matlab schreiben
            //Teilmatritzen gemäß Globalid extrahieren
            //Mit ExternalRows vergleichen
            //Die große Frage: funktioniert der batchmode connector parallel? Beim rausschreiben beachten

            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("GetExternalRowsTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- setup mgo and mask
            MultigridOperator mgo = Utils.CreateTestMGOperator(UseXdg, DGOrder, MatrixShape.laplace, Res);
            MultigridMapping map = mgo.Mapping;
            BlockMsrMatrix M = mgo.OperatorMatrix;
            
            //Delete this plz ...
            //M.SaveToTextFileSparse("M");
            //int[] A = Utils.GimmeAllBlocksWithSpec(map, 9);
            //int[] B = Utils.GimmeAllBlocksWithSpec(map, 18);
            //if (map.MpiRank == 0) {
            //    A.SaveToTextFileDebug("ACells");
            //    B.SaveToTextFileDebug("BCells");
            //}
            var selector = new SubBlockSelector(map);
            var dummy = new BlockMsrMatrix(map); // we are only interested in getting indices, so a dummy is sufficient
            var mask = new TestMask(selector, dummy);

            //Arrange --- get stuff to put into matlab
            int[] GlobalIdx_ext = Utils.GetAllExtCellIdc(map);
            double[] GlobIdx = GlobalIdx_ext.Length.ForLoop(i => (double)GlobalIdx_ext[i]+1.0);

            //Arrange --- get external rows by mask
            BlockMsrMatrix extrows = BlockMask.GetAllExternalRows(mgo.Mapping, mgo.OperatorMatrix);

            //Assert --- idc and rows of extrows have to be the same
            Assert.IsTrue(GlobIdx.Length == extrows._RowPartitioning.LocalLength);

            //Arrange --- get external rows by matlab
            var infNorm = MultidimensionalArray.Create(1, 1);
            using (BatchmodeConnector matlab = new BatchmodeConnector()) {
                //note: BatchmodeCon maybe working on proc0 but savetotxt file, etc. (I/O) is full mpi parallel
                //so concider this as full mpi-parallel
                matlab.PutSparseMatrix(M, "M");
                matlab.PutSparseMatrix(extrows,"M_test");
                matlab.PutVector(GlobIdx, "Idx");                
                matlab.Cmd(String.Format("M_ext = M(Idx, :);"));
                matlab.Cmd("n=norm(M_test-M_ext,inf)");
                matlab.GetMatrix(infNorm,"n");
                matlab.Execute();
            }

            //Assert --- test if we actually got the right Matrix corresponding to Index
            Assert.IsTrue(infNorm[0,0]==0.0);
        }

        [Test]
        public static void FastSubMatrixExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.laplace)] MatrixShape MShape,
            [Values(4)] int Res
            ) {

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("FastSubMatrixExtraction({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange ---
            MultigridOperator mgo = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res);
            MultigridMapping map = mgo.Mapping;
            BlockMsrMatrix M = mgo.OperatorMatrix;

            var sbs = new SubBlockSelector(map);
            int[] extcells = sbs.AllExternalCellsSelection();
            var M_ext = BlockMask.GetAllExternalRows(map, M);
            var mask = new BlockMask(sbs, M_ext);

            //Arrange --- get index list of all external cells
            int[] idc = Utils.GetAllExtCellIdc(map);
            double[] GlobIdx = idc.Count().ForLoop(i => (double)idc[i] + 1.0);

            //Arrange --- stopwatch
            var stw = new Stopwatch();
            stw.Reset();

            //Act --- Extract SubMatrix
            stw.Start();
            BlockMsrMatrix subM = mask.GetSubBlockMatrix(M);
            stw.Stop();

            //Arrange --- Extract Blocks in Matlab and substract
            var infNorm = MultidimensionalArray.Create(4, 1);
            int rank = map.MpiRank;
            using (BatchmodeConnector matlab = new BatchmodeConnector()) {
                matlab.PutSparseMatrix(M, "M");
                // note: M_sub lives on Comm_Self, therefore we have to distinguish between procs ...
                matlab.PutSparseMatrixRankExclusive(subM, "M_sub");
                matlab.PutVectorRankExclusive(GlobIdx, "Idx");
                matlab.Cmd("M_0 = M(Idx_0, Idx_0);");
                matlab.Cmd("M_1 = M(Idx_1, Idx_1);");
                matlab.Cmd("M_2 = M(Idx_2, Idx_2);");
                matlab.Cmd("M_3 = M(Idx_3, Idx_3);");
                matlab.Cmd("n=[0; 0; 0; 0];");
                matlab.Cmd("n(1,1)=norm(M_0-M_sub_0,inf);");
                matlab.Cmd("n(2,1)=norm(M_1-M_sub_1,inf);");
                matlab.Cmd("n(3,1)=norm(M_2-M_sub_2,inf);");
                matlab.Cmd("n(4,1)=norm(M_3-M_sub_3,inf);");
                matlab.GetMatrix(infNorm, "n");
                matlab.Execute();
            }

            //Assert --- mask blocks and extracted blocks are the same
            Assert.IsTrue(infNorm[rank, 0] == 0.0);
        }


        [Test]
        public static void SubMatrixExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.full_var_spec, MatrixShape.full_var, MatrixShape.full_spec, MatrixShape.full)] MatrixShape MShape,
            [Values(4)] int Res
            ) {

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("SubMatrixExtraction({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange ---
            MultigridOperator mgo = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res);
            MultigridMapping map = mgo.Mapping;
            BlockMsrMatrix M = mgo.OperatorMatrix;
            var sbs = new SubBlockSelector(map);
            int[] extcells = sbs.AllExternalCellsSelection();
            var M_ext=BlockMask.GetAllExternalRows(map,M);
            var mask = new BlockMask(sbs, M_ext);
            bool[] coup = Utils.SetCoupling(MShape);


            //Arrange --- get index list of all external cells
            int[] idc = Utils.GetAllExtCellIdc(map);
            double[] GlobIdx = idc.Count().ForLoop(i => (double)idc[i]+1.0);

            //Arrange --- stopwatch
            var stw = new Stopwatch();
            stw.Reset();

            //Act --- Extract SubMatrix
            stw.Start();
            BlockMsrMatrix subM = mask.GetSubBlockMatrix(M, false , coup[1], coup[2]);
            stw.Stop();

            //Arrange --- Extract Blocks in Matlab and substract
            var infNorm = MultidimensionalArray.Create(4, 1);
            int rank = map.MpiRank;
            using (BatchmodeConnector matlab = new BatchmodeConnector()) {
                matlab.PutSparseMatrix(M, "M");
                // note: M_sub lives on Comm_Self, therefore we have to distinguish between procs ...
                matlab.PutSparseMatrixRankExclusive(subM, "M_sub"); 
                matlab.PutVectorRankExclusive(GlobIdx, "Idx");
                matlab.Cmd("M_0 = M(Idx_0, Idx_0);");
                matlab.Cmd("M_1 = M(Idx_1, Idx_1);");
                matlab.Cmd("M_2 = M(Idx_2, Idx_2);");
                matlab.Cmd("M_3 = M(Idx_3, Idx_3);");
                matlab.Cmd("n=[0; 0; 0; 0];");
                matlab.Cmd("n(1,1)=norm(M_0-M_sub_0,inf);");
                matlab.Cmd("n(2,1)=norm(M_1-M_sub_1,inf);");
                matlab.Cmd("n(3,1)=norm(M_2-M_sub_2,inf);");
                matlab.Cmd("n(4,1)=norm(M_3-M_sub_3,inf);");
                matlab.GetMatrix(infNorm, "n");
                matlab.Execute();
            }

            //Assert --- mask blocks and extracted blocks are the same
            Assert.IsTrue(infNorm[rank,0] == 0.0);
        }

        [Test]
        public static void SubBlockExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal_var_spec, MatrixShape.diagonal_spec, MatrixShape.diagonal_var, MatrixShape.diagonal)] MatrixShape MShape,
            [Values(4)] int Res
            ) {

            Utils.TestInit((int) UseXdg, DGOrder, (int) MShape);
            Console.WriteLine("SubMatrixIgnoreCoupling({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- create test matrix and MG mapping
            MultigridOperator mgo = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res);
            MultigridMapping map = mgo.Mapping;
            BlockMsrMatrix M = mgo.OperatorMatrix;

            //Arrange --- masking of all external cells
            var sbs = new SubBlockSelector(map);
            sbs.AllExternalCellsSelection();
            var M_ext = BlockMask.GetAllExternalRows(map, M);
            var mask = new BlockMask(sbs, M_ext);
            //bool[] coup = Utils.SetCoupling(MShape);

            //Arrange --- get index dictonary of all external cell indices
            Dictionary<int,int[]> Didc = Utils.GetDictOfAllExtCellIdc(map);

            //Arrange --- stopwatch
            var stw = new Stopwatch();
            stw.Reset();

            //Act --- Extract subblocks
            stw.Start();
            //var eblocks = mask.GetSubBlocks(M,coup[0],coup[1],coup[2]);
            var eblocks = mask.GetSubBlocks(M, true, false, false);
            stw.Stop();

            //Assert --- same number of blocks?
            Assert.IsTrue(eblocks.Length==M_ext._RowPartitioning.LocalNoOfBlocks);

            bool test = eblocks.Length.MPIEquals();
            Debug.Assert(test);
            for (int iBlock = 0; iBlock < eblocks.Length; iBlock++) {
                
               
                var infNorm = MultidimensionalArray.Create(4, 1);
                int rank = map.MpiRank;
                int ExtBlockIdx = iBlock + map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                Didc.TryGetValue(ExtBlockIdx, out int[] idc);

                using (BatchmodeConnector matlab = new BatchmodeConnector()) {

                    double[] GlobIdx = idc.Count().ForLoop(i => (double)idc[i] + 1.0);

                    matlab.PutSparseMatrix(M, "M");
                    // note: M_sub lives on Comm_Self, therefore we have to distinguish between procs ...
                    matlab.PutMatrixRankExclusive(eblocks[iBlock], "M_sub");
                    matlab.PutVectorRankExclusive(GlobIdx, "Idx");
                    matlab.Cmd("M_0 = full(M(Idx_0, Idx_0));");
                    matlab.Cmd("M_1 = full(M(Idx_1, Idx_1));");
                    matlab.Cmd("M_2 = full(M(Idx_2, Idx_2));");
                    matlab.Cmd("M_3 = full(M(Idx_3, Idx_3));");
                    matlab.Cmd("n=[0; 0; 0; 0];");
                    matlab.Cmd("n(1,1)=norm(M_0-M_sub_0,inf);");
                    matlab.Cmd("n(2,1)=norm(M_1-M_sub_1,inf);");
                    matlab.Cmd("n(3,1)=norm(M_2-M_sub_2,inf);");
                    matlab.Cmd("n(4,1)=norm(M_3-M_sub_3,inf);");
                    matlab.GetMatrix(infNorm, "n");
                    matlab.Execute();
                }
                Assert.IsTrue(infNorm[rank, 0] < 1e-13);
            }
        }


        public static void CellblocksIgnoreCoupling() {

        }

        public static void SubSelection() { 
        
        }

        public static void VectorCellwiseOperation() {
        
        }

        public static void VectorSplitOperation() {
        
        }

        [Test]
        public static void ExternalIndexTest(
        [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
        [Values(2)] int DGOrder,
        [Values(4)] int Res
        ) {

            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("ExternalIndexTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- Get global index by mapping
            MultigridOperator MGOp = Utils.CreateTestMGOperator(UseXdg, DGOrder,MatrixShape.laplace, Res);
            var map = MGOp.Mapping;
            int[] GlobalIdxMap_ext = Utils.GetAllExtCellIdc(map);

            //Arrange --- Prepare stuff for mask
            var selector = new SubBlockSelector(map);
            var dummy = new BlockMsrMatrix(map); // we are only interested in getting indices, so a dummy is sufficient
            var stw = new Stopwatch();
            stw.Reset();

            //Act --- do the masking to get index lists
            stw.Start();
            var mask = new TestMask(selector, dummy);
            stw.Stop();
            int[] GlobalIdxMask_ext = mask.Global_IList_ExternalCells;

            //Assert --- Idx lists are of same length
            Assert.IsTrue(GlobalIdxMap_ext.Length == GlobalIdxMask_ext.Length);

            //Assert --- Compare map and mask indices
            for (int iLoc = 0; iLoc < GlobalIdxMask_ext.Length; iLoc++) {
                Assert.IsTrue(GlobalIdxMask_ext[iLoc] == GlobalIdxMap_ext[iLoc]);
            }
        }

    }
}

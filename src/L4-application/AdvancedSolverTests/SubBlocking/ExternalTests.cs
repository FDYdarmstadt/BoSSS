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
        [Values(2)] int DGOrder) {

            //Matlabaufruf --> gesamte Matrix nach Matlab schreiben
            //Teilmatritzen gemäß Globalid extrahieren
            //Mit ExternalRows vergleichen
            //Die große Frage: funktioniert der batchmode connector parallel? Beim rausschreiben beachten

            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("GetExternalRowsTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- setup mgo and mask
            MultigridOperator mgo = Utils.CreateTestMGOperator(UseXdg, DGOrder, MatrixShape.full_var_spec, 4);
            MultigridMapping map = mgo.Mapping;
            BlockMsrMatrix M = mgo.OperatorMatrix;
            var selector = new SubBlockSelector(map);
            var dummy = new BlockMsrMatrix(map); // we are only interested in getting indices, so a dummy is sufficient
            var mask = new TestMask(selector, dummy);

            //Arrange --- get stuff to put into matlab
            int[] GlobalIdxMask_ext = mask.Global_IList_ExternalCells;
            double[] GlobIdx = GlobalIdxMask_ext.Length.ForLoop(i => (double)GlobalIdxMask_ext[i]+1.0);

            //Arrange --- get external rows by mask
            BlockMsrMatrix extrows = BlockMask.GetAllExternalRows(mgo.Mapping, mgo.OperatorMatrix);

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
        public static void SubMatrixExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.full_var_spec)] MatrixShape MShape
            ) {

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("ExternalVecOperation({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange ---
            MultigridOperator mgo = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, 2);
            MultigridMapping map = mgo.Mapping;
            BlockMsrMatrix M = mgo.OperatorMatrix;
            var sbs = new SubBlockSelector(map);
            int[] extcells = sbs.AllExternalCellsSelection();
            var M_ext=BlockMask.GetAllExternalRows(map,M);
            var mask = new BlockMask(sbs, M_ext);

            //Arrange --- get index list of all external cells
            List<int> idc = new List<int>();
            for (int i = 0; i < extcells.Length; i++) {
                int iCell = extcells[i];
                idc.AddRange(map.GetIndcOfExtCell(iCell));
            }
            double[] GlobIdx = idc.Count().ForLoop(i => (double)idc[i]+1.0);

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
            Assert.IsTrue(infNorm[rank,0] == 0.0);
        }


    


        [Test]
        public static void ExternalIndexTest(
        [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
        [Values(2)] int DGOrder) {

            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("ExternalIndexTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- Get global index by mapping
            MultigridOperator MGOp = Utils.CreateTestMGOperator(UseXdg, DGOrder);
            var map = MGOp.Mapping;
            int[] fields = map.NoOfVariables.ForLoop(i => i);
            int[] GlobalIdxMap_ext = map.GetSubvectorIndices_Ext(fields); // maybe not sorted, so do it ...
            Array.Sort(GlobalIdxMap_ext);

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

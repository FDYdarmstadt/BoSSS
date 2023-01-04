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

namespace AdvancedSolverTests.SubBlocking {
    [TestFixture]
    static public class ExternalTests {

        [Test]
        public static void GetExternalRowsTest(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(4)] int Res) {
            

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("GetExternalRowsTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- setup mgo and mask
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MatrixShape.laplace, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;

                //Arrange --- get stuff to put into matlab
                long[] GlobalIdx_ext = Utils.GetAllExtCellIdc(map);
                double[] GlobIdx = GlobalIdx_ext.Length.ForLoop(i => (double)GlobalIdx_ext[i] + 1.0);

                //Arrange --- get external rows by mask
                BlockMsrMatrix extrows = BlockMask.GetAllExternalRows(mgo.Mapping, mgo.OperatorMatrix);

                //Assert --- idc and rows of extrows have to be the same
                Assert.IsTrue(GlobIdx.Length == extrows._RowPartitioning.LocalLength);

                //Arrange --- get external rows by matlab
                var infNorm = MultidimensionalArray.Create(1, 1);
                using(BatchmodeConnector matlab = new BatchmodeConnector()) {
                    //note: BatchmodeCon maybe working on proc0 but savetotxt file, etc. (I/O) is full mpi parallel
                    //so concider this as full mpi-parallel
                    matlab.PutSparseMatrix(M, "M");
                    matlab.PutSparseMatrix(extrows, "M_test");
                    matlab.PutVector(GlobIdx, "Idx");
                    matlab.Cmd(String.Format("M_ext = M(Idx, :);"));
                    matlab.Cmd("n=norm(M_test-M_ext,inf)");
                    matlab.GetMatrix(infNorm, "n");
                    matlab.Execute();
                }

                //Assert --- test if we actually got the right Matrix corresponding to Index
                Assert.IsTrue(infNorm[0, 0] == 0.0);
            }
        }

        [Test]
        public static void FastSubMatrixExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.laplace)] MatrixShape MShape,
            [Values(4)] int Res
            ) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("FastSubMatrixExtraction({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange ---
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;

                var sbs = new SubBlockSelector(map);
                int[] extcells = sbs.AllExternalCellsSelection();
                var M_ext = BlockMask.GetAllExternalRows(map, M);
                var mask = new BlockMask(sbs, M_ext);

                //Arrange --- get index list of all external cells
                long[] idc = Utils.GetAllExtCellIdc(map);
                double[] GlobIdx = idc.Count().ForLoop(i => (double)idc[i] + 1.0);

                //Arrange --- stopwatch
                var stw = new Stopwatch();
                stw.Reset();

                //Act --- Extract SubMatrix
                stw.Start();
                BlockMsrMatrix subM = mask.GetSubBlockMatrix_MpiSelf(M);
                stw.Stop();

                //Arrange --- Extract Blocks in Matlab and substract
                var infNorm = MultidimensionalArray.Create(4, 1);
                int rank = map.MpiRank;
                using(BatchmodeConnector matlab = new BatchmodeConnector()) {
                    matlab.PutSparseMatrix(M, "M");
                    // note: M_sub lives on Comm_Self, therefore we have to distinguish between procs ...
                    matlab.PutSparseMatrixPerMPIrank(subM, "M_sub");
                    matlab.PutVectorPerMPIrank(GlobIdx, "Idx");
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
        }

        [Test]
        public static void SubMatrixExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.full_var_spec, MatrixShape.full_var, MatrixShape.full_spec, MatrixShape.full)] MatrixShape MShape,
            [Values(4)] int Res
            ) {
            // --test=AdvancedSolverTests.SubBlocking.ExternalTests.SubMatrixExtraction

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("SubMatrixExtraction({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange ---
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;
                var sbs = new SubBlockSelector(map);
                int[] extcells = sbs.AllExternalCellsSelection();
                var M_ext = BlockMask.GetAllExternalRows(map, M);
                var mask = new BlockMask(sbs, M_ext);
                bool[] coup = Utils.SetCoupling(MShape);

                //Arrange --- get index list of all external cells
                long[] idc = Utils.GetAllExtCellIdc(map);
                double[] GlobIdx = idc.Count().ForLoop(i => (double)idc[i] + 1.0);

                //Arrange --- stopwatch
                var stw = new Stopwatch();
                stw.Reset();

                //Act --- Extract SubMatrix
                stw.Start();
                BlockMsrMatrix subM = mask.GetSubBlockMatrix(M, false, coup[0], false);
                // ignoring of species and including cell coupling can not be tested with this setup
                stw.Stop();

                //Arrange --- Extract Blocks in Matlab and substract
                var infNorm = MultidimensionalArray.Create(4, 1);
                int rank = map.MpiRank;
                using(BatchmodeConnector matlab = new BatchmodeConnector()) {
                    matlab.PutSparseMatrix(M, "M");
                    // note: M_sub lives on Comm_Self, therefore we have to distinguish between procs ...
                    matlab.PutSparseMatrixPerMPIrank(subM, "M_sub");
                    matlab.PutVectorPerMPIrank(GlobIdx, "Idx");
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
        }

        [Test]
        public static void SubBlockExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal_var_spec, MatrixShape.diagonal_spec, MatrixShape.diagonal_var, MatrixShape.diagonal)] MatrixShape MShape,
            [Values(4)] int Res
            ) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("SubMatrixIgnoreCoupling({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- create test matrix and MG mapping
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;

                //Arrange --- masking of all external cells
                var sbs = new SubBlockSelector(map);
                sbs.AllExternalCellsSelection();
                var M_ext = BlockMask.GetAllExternalRows(map, M);
                var mask = new BlockMask(sbs, M_ext);
                bool[] coup = Utils.SetCoupling(MShape);

                //Arrange --- get index dictonary of all external cell indices
                Dictionary<int, long[]> Didc = Utils.GetDictOfAllExtCellIdc(map);

                //Arrange --- stopwatch
                var stw = new Stopwatch();
                stw.Reset();

                //Act --- Extract subblocks
                stw.Start();
                var eblocks = mask.GetDiagonalBlocks(M, coup[0], coup[1]);
                stw.Stop();

                //Assert --- same number of blocks?
                Assert.IsTrue(eblocks.Length == M_ext._RowPartitioning.LocalNoOfBlocks);

                bool test = eblocks.Length.MPIEquals();
                Debug.Assert(test);
                for(int iBlock = 0; iBlock < eblocks.Length; iBlock++) {


                    var infNorm = MultidimensionalArray.Create(4, 1);
                    int rank = map.MpiRank;
                    int ExtBlockIdx = iBlock + map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                    Didc.TryGetValue(ExtBlockIdx, out long[] idc);

                    using(BatchmodeConnector matlab = new BatchmodeConnector()) {

                        double[] GlobIdx = idc.Count().ForLoop(i => (double)idc[i] + 1.0);
                        Assert.IsTrue(GlobIdx.Length == eblocks[iBlock].Lengths[0]);
                        MsrMatrix M_sub = eblocks[iBlock].ConvertToMsr();

                        matlab.PutSparseMatrix(M, "M");
                        // note: M_sub lives on Comm_Self, therefore we have to distinguish between procs ...
                        matlab.PutSparseMatrixPerMPIrank(M_sub, "M_sub");
                        matlab.PutVectorPerMPIrank(GlobIdx, "Idx");
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
                    Assert.IsTrue(infNorm[rank, 0] == 0.0); //
                }
            }
        }

        [Test]
        public static void VectorCellwiseOperation(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal_var_spec, MatrixShape.diagonal_spec, MatrixShape.diagonal_var, MatrixShape.diagonal)] MatrixShape MShape,
            [Values(4)] int Res
            ) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("SubMatrixIgnoreCoupling({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- create test matrix, MG mapping
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;


                //Arrange --- masking and subblock extraction of external cells
                var sbs = new SubBlockSelector(map);
                sbs.AllExternalCellsSelection();
                var M_ext = BlockMask.GetAllExternalRows(map, M);
                var mask = new BlockMask(sbs, M_ext);
                var eblocks = mask.GetDiagonalBlocks(M, false, false);
                //Dictionary<int, int[]> Didc = Utils.GetDictOfAllExtCellIdc(map);

                //Arrange --- generate rnd vector and distribute it
                double[] vec = new double[map.LocalLength];
                vec = Utils.GetRandomVector(map.LocalLength);
                var vec_ex = new MPIexchange<double[]>(map, vec);
                vec_ex.TransceiveStartImReturn();
                vec_ex.TransceiveFinish(0.0);
                Debug.Assert(vec_ex.Vector_Ext.L2Norm() != 0);

                //Arrange --- stopwatch
                var stw = new Stopwatch();
                stw.Reset();

                //Arrange --- get extended (loc+external cells) vector
                double[] Vec_ext = new double[vec.Length + vec_ex.Vector_Ext.Length];
                mask.AccSubVec(vec_ex.Vector_Ext, Vec_ext);

                bool test = eblocks.Length.MPIEquals();
                Debug.Assert(test);
                //Act --- calculate blockwise result: M_i*vec_i=Res_i
                double[] Res_ext = new double[Vec_ext.Length];
                stw.Start();
                for(int i = 0; i < eblocks.Length; i++) {
                    //int iBlock = i + map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                    double[] vec_i = mask.GetSubVecOfCell(Vec_ext, i);
                    double[] Res_i = new double[vec_i.Length];
                    eblocks[i].MatVecMul(1.0, vec_i, 0.0, Res_i);
                    mask.AccSubVecOfCell(Res_i, i, Res_ext);
                    if(map.MpiRank == 0) {
                        eblocks[i].ConvertToMsr().SaveToTextFileSparseDebug(String.Format("block_{0}_{1}", i, map.MpiRank));
                        vec_i.SaveToTextFileDebug(String.Format("vec_{0}_{1}", i, map.MpiRank));
                        Res_i.SaveToTextFileDebug(String.Format("Res_{0}_{1}", i, map.MpiRank));
                    }
                }
                stw.Stop();

                //Act --- project Res_i onto Res_g and Res_g=M_ext*vec_ext-Res_g
                double[] Res_g = mask.GetSubVec(Res_ext);
                var qM_ext = M_ext.ConvertToQuadraticBMsr(mask.GlobalIndices_External.ToArray(), false);
                qM_ext.SpMV(1.0, vec_ex.Vector_Ext, -1.0, Res_g);

                if(map.MpiRank == 0) {
                    vec_ex.Vector_Ext.SaveToTextFileDebug("vec_g");
                    Res_g.SaveToTextFileDebug("Res_g");
                    M_ext.SaveToTextFileSparseDebug("M_ext");
                    qM_ext.SaveToTextFileSparseDebug("qM_ext");
                }

                //Assert --- |Res_g| should be at least near to zero
                Console.WriteLine("VectorCellwiseOperation: Residual is: " + Res_g.L2Norm());
                Assert.Less(Res_g.L2Norm(), 1e-13); // not exactly zero, can depend on processor architecture
            }
        }

        [Test]
        public static void SubSelection(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.full_var_spec, MatrixShape.full_spec, MatrixShape.full_var, MatrixShape.full)] MatrixShape MShape,
            [Values(4)] int Res) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }


            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape, Res);
            Console.WriteLine("SubSelection({0},{1},{2},{3})", UseXdg, DGOrder, MShape, Res);

            //Arrange --- create test matrix, MG mapping
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;

                //Arrange --- get mask
                int[] cells = Utils.GetCellsOfOverlappingTestBlock(map);
                Array.Sort(cells);
                var sbs = new SubBlockSelector(map);
                sbs.CellSelector(cells, false);
                BlockMsrMatrix M_ext = BlockMask.GetAllExternalRows(map, M);
                var mask = new BlockMask(sbs, M_ext);

                //Arrange --- get GlobalIdxList
                long[] idc = Utils.GetIdcOfSubBlock(map, cells);
                bool[] coup = Utils.SetCoupling(MShape);

                var M_sub = mask.GetSubBlockMatrix(M, false, coup[0], coup[1]);

                var infNorm = MultidimensionalArray.Create(4, 1);
                int rank = map.MpiRank;
                using(BatchmodeConnector matlab = new BatchmodeConnector()) {

                    double[] GlobIdx = idc.Count().ForLoop(i => (double)idc[i] + 1.0);
                    Assert.IsTrue(GlobIdx.Length == M_sub.NoOfRows);

                    matlab.PutSparseMatrix(M, "M");
                    // note: M_sub lives on Comm_Self, therefore we have to distinguish between procs ...
                    matlab.PutSparseMatrixPerMPIrank(M_sub, "M_sub");
                    matlab.PutVectorPerMPIrank(GlobIdx, "Idx");
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
                Console.WriteLine("SubSelection: Inf-Norm is: " + infNorm[rank, 0]);
                Assert.IsTrue(infNorm[rank, 0] == 0.0);
            }
        }

        /// <summary>
        /// this test is for the version, which has local and external vector as imput
        /// </summary>
        /// <param name="UseXdg"></param>
        /// <param name="DGOrder"></param>
        /// <param name="MShape"></param>
        /// <param name="Res"></param>
        [Test]
        public static void VectorSplitOperation_Op1(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.full_var_spec, MatrixShape.full_spec, MatrixShape.full)] MatrixShape MShape,
            [Values(4)] int Res) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape, Res);
            Console.WriteLine("VectorSplitOperation({0},{1},{2},{3})", UseXdg, DGOrder, MShape, Res);

            //Arrange --- create test matrix, MG mapping
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;
                BlockMsrMatrix M_ext = BlockMask.GetAllExternalRows(map, M);
                double[] Vec = Utils.GetRandomVector(M_ext.RowPartitioning.LocalLength);

                //Arrange --- setup masking
                SubBlockSelector sbsA = new SubBlockSelector(map);
                sbsA.SetDefaultSplitSelection(MShape, true, false);
                BlockMask maskA = new BlockMask(sbsA, M_ext);

                SubBlockSelector sbsB = new SubBlockSelector(map);
                sbsB.SetDefaultSplitSelection(MShape, false, false);
                BlockMask maskB = new BlockMask(sbsB, M_ext);

                double[] VecAB = new double[Vec.Length];

                //Arrange --- some time measurement
                Stopwatch stw = new Stopwatch();
                stw.Reset();

                //Act --- 
                stw.Start();
                var VecA = maskA.GetSubVec(Vec, new double[0]);
                var VecB = maskB.GetSubVec(Vec, new double[0]);

                maskA.AccSubVec(VecA, VecAB, new double[0]);
                maskB.AccSubVec(VecB, VecAB, new double[0]);
                stw.Stop();

                Debug.Assert(Vec.L2Norm() != 0);
                double fac = ((MShape == MatrixShape.full_var || MShape == MatrixShape.diagonal_var) && UseXdg == XDGusage.none) ? -2.0 : -1.0;
                VecAB.AccV(fac, Vec);

                //Assert --- are extracted blocks and 
                Console.WriteLine("VectorSplitOperation_Op1: L2 norm is " + VecAB.L2Norm());
                Assert.IsTrue(VecAB.L2Norm() == 0.0, String.Format("L2Norm neq 0!"));
            }
        }

        /// <summary>
        /// This test is for the version, which has an extended vector as imput.
        /// Extended vector contains local and external parts
        /// </summary>
        /// <param name="UseXdg"></param>
        /// <param name="DGOrder"></param>
        /// <param name="MShape"></param>
        /// <param name="Res"></param>
        [Test]
        public static void VectorSplitOperation_Op2([Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.full_var_spec, MatrixShape.full_spec, MatrixShape.full)] MatrixShape MShape,
            [Values(16)] int Res) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape, Res);
            Console.WriteLine("VectorSplitOperation({0},{1},{2},{3})", UseXdg, DGOrder, MShape, Res);

            //Arrange --- create test matrix, MG mapping
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;
                BlockMsrMatrix M_ext = BlockMask.GetAllExternalRows(map, M);

                //Arrange --- setup masking
                SubBlockSelector sbsA = new SubBlockSelector(map);
                sbsA.SetDefaultSplitSelection(MShape, true, false);
                BlockMask maskA = new BlockMask(sbsA, M_ext);
                SubBlockSelector sbsB = new SubBlockSelector(map);
                sbsB.SetDefaultSplitSelection(MShape, false, false);
                BlockMask maskB = new BlockMask(sbsB, M_ext);
                SubBlockSelector sbsFull = new SubBlockSelector(map);
                BlockMask maskFull = new BlockMask(sbsFull, M_ext);

                //Arrange --- the test vectors
                double[] locVec = Utils.GetRandomVector(M.RowPartitioning.LocalLength);
                MPIexchange<double[]> Exchange;
                Exchange = new MPIexchange<double[]>(map, locVec);
                Exchange.TransceiveStartImReturn();
                Exchange.TransceiveFinish(0.0);
                var extVec = Exchange.Vector_Ext;
                var fullVec = maskFull.GetSubVec(extVec, new double[locVec.Length]);
                var test = new double[fullVec.Length];
                Debug.Assert(fullVec.L2Norm() != 0);

                //Act --- check operations
                var subA = maskA.GetSubVec(fullVec);
                var subB = maskB.GetSubVec(fullVec);

                maskA.AccSubVec(subA, test);
                maskB.AccSubVec(subB, test);

                double fac = ((MShape == MatrixShape.full_var || MShape == MatrixShape.diagonal_var) && UseXdg == XDGusage.none) ? -2.0 : -1.0;
                test.AccV(fac, fullVec);

                //Assert --- are extracted blocks and 
                Console.WriteLine("VectorSplitOperation_Op2: L2 norm is " + test.L2Norm());
                Assert.IsTrue(test.L2Norm() == 0.0, String.Format("L2Norm neq 0!"));
            }
        }

        /// <summary>
        /// Delete this, if one of the operations is dismissed
        /// </summary>
        /// <param name="UseXdg"></param>
        /// <param name="DGOrder"></param>
        /// <param name="MShape"></param>
        /// <param name="Res"></param>
        [Test]
        public static void CompareVectorOps([Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.full_var_spec, MatrixShape.full_spec, MatrixShape.full)] MatrixShape MShape,
            [Values(16)] int Res) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder, (int)MShape, Res);
            Console.WriteLine("VectorSplitOperation({0},{1},{2},{3})", UseXdg, DGOrder, MShape, Res);

            //Arrange --- create test matrix, MG mapping
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MShape, Res)) {
                MultigridOperator mgo = O.MGOp;
                MultigridMapping map = mgo.Mapping;
                BlockMsrMatrix M = mgo.OperatorMatrix;
                BlockMsrMatrix M_ext = BlockMask.GetAllExternalRows(map, M);

                //Arrange --- setup masking
                SubBlockSelector sbsA = new SubBlockSelector(map);
                sbsA.SetDefaultSplitSelection(MShape, true, false);
                BlockMask maskA = new BlockMask(sbsA, M_ext);
                SubBlockSelector sbsFull = new SubBlockSelector(map);
                BlockMask maskFull = new BlockMask(sbsFull, M_ext);

                //Arrange --- the test vectors
                double[] locVec = Utils.GetRandomVector(M.RowPartitioning.LocalLength);
                MPIexchange<double[]> Exchange;
                Exchange = new MPIexchange<double[]>(map, locVec);
                Exchange.TransceiveStartImReturn();
                Exchange.TransceiveFinish(0.0);
                var extVec = Exchange.Vector_Ext;
                var exttmp = new double[extVec.Length];
                var fullVec = maskFull.GetSubVec(extVec, new double[locVec.Length]);
                var test_two = new double[fullVec.Length];
                Debug.Assert(fullVec.L2Norm() != 0);

                //Act --- Compare operations
                var sub_two = maskA.GetSubVec(fullVec);
                var sub_one = maskA.GetSubVec(extVec, new double[0]);
                sub_one.AccV(-1.0, sub_two);

                maskA.AccSubVec(sub_two, test_two);
                maskA.AccSubVec(sub_two, exttmp, new double[0]);
                var test_one = maskFull.GetSubVec(exttmp, new double[locVec.Length]);
                test_one.AccV(-1.0, test_two);

                //Assert --- check equality
                Assert.IsTrue(sub_one.L2Norm() == 0.0, String.Format("operations differ in getsubvec"));
                Assert.IsTrue(test_one.L2Norm() == 0.0, String.Format("operations differ in accsubvec"));
            }
        }

        [Test]
        public static void ExternalIndexTest(
        [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
        [Values(2)] int DGOrder,
        [Values(4)] int Res
        ) {

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int MPISize);
            if(MPISize <= 1) {
                Console.WriteLine("Terminating: this test is supposed to be run in parallel, but currently running on only 1 processor.");
                return;
            }

            Utils.TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("ExternalIndexTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- Get global index by mapping
            using(var O = Utils.CreateTestMGOperator(UseXdg, DGOrder, MatrixShape.laplace, Res)) {
                MultigridOperator MGOp = O.MGOp;
                var map = MGOp.Mapping;
                long[] GlobalIdxMap_ext = Utils.GetAllExtCellIdc(map);

                //Arrange --- Prepare stuff for mask
                var selector = new SubBlockSelector(map);
                var dummy = new BlockMsrMatrix(map); // we are only interested in getting indices, so a dummy is sufficient
                var stw = new Stopwatch();
                stw.Reset();

                //Act --- do the masking to get index lists
                stw.Start();
                var mask = new BlockMask(selector, dummy);
                stw.Stop();
                long[] GlobalIdxMask_ext = mask.GlobalIndices_External.ToArray();

                //Assert --- Idx lists are of same length
                Assert.IsTrue(GlobalIdxMap_ext.Length == GlobalIdxMask_ext.Length);

                //Assert --- Compare map and mask indices
                for(int iLoc = 0; iLoc < GlobalIdxMask_ext.Length; iLoc++) {
                    Assert.IsTrue(GlobalIdxMask_ext[iLoc] == GlobalIdxMap_ext[iLoc]);
                }
            }
        }
    }
}

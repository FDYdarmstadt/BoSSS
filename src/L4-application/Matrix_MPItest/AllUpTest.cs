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

using BoSSS.Platform;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.Matrix_MPItest {

    [TestFixture]
    static public class AllUpTest {

        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp() {
            bool MpiInit;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out MpiInit);

            if (System.Environment.MachineName.ToLowerInvariant().EndsWith("rennmaschin")
                //|| System.Environment.MachineName.ToLowerInvariant().Contains("jenkins")
                ) {
                // This is Florians Laptop;
                // he is to poor to afford MATLAB, so he uses OCTAVE
                BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;
                BatchmodeConnector.MatlabExecuteable = "C:\\cygwin64\\bin\\bash.exe";
            }

        }

        /// <summary>
        /// MPI shutdown.
        /// </summary>
        [TestFixtureTearDown]
        public static void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

         public static void Main(string[] args) {
            SetUp();

            /*
            int counter = 0;
            foreach (var x in new XDGusage[] { XDGusage.all, XDGusage.mixed1, XDGusage.mixed2, XDGusage.none }) {
                foreach (int p in new int[] { 2 }) {
                    foreach (bool b in new bool[] { false }) {
                        foreach (bool c in new bool[] { false }) {

                            SubMatrixTest(x, p, b, c);
                            MultiplyTest(x, p, b, c);
                            SpMVTest(x, p, b, c);
                            counter++;

                        }
                    }
                }
            }
            //*/


            BoSSS.Application.Matrix_MPItest.AllUpTest.MultiplyTest(XDGusage.none, 3, true, false);
            //MultiplyTest(XDGusage.none, 1, false, false);
            //SubMatrixTest(XDGusage.none, 2, false, false);
            //MultiplyTest(XDGusage.none, 2, false, false);
            //SpMVTest(XDGusage.none, 2, false, false);

            TestFixtureTearDown();

            Console.WriteLine("TOTAL Time spend in matrix operations: " + TotTime_MatrixOp.TotalSeconds + " sec.");

        }

        
        /// <summary>
        /// Test for matrix/matrix multiplication.
        /// </summary>
        [Test]
        public static void MultiplyTest(
            [Values(XDGusage.none, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(false)] bool compressL1,
            [Values(false, true)] bool compressL2) { 

            unsafe
            {
                int[] Params = new int[8], ParamsGlob = new int[8];
                fixed(int* pParams = Params, pParamsGlob = ParamsGlob)
                {
                    pParams[0] = (int)UseXdg;
                    pParams[1] = DGOrder;
                    pParams[2] = compressL1 ? 1 : 0;
                    pParams[3] = compressL2 ? 1 : 0;
                    pParams[4] = -pParams[0];
                    pParams[5] = -pParams[1]; 
                    pParams[6] = -pParams[2];
                    pParams[7] = -pParams[3];

                    csMPI.Raw.Allreduce((IntPtr)pParams, (IntPtr)pParamsGlob, 8, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                }

                int[] ParamsMin = ParamsGlob.GetSubVector(0, 4);
                int[] ParamsMax = ParamsGlob.GetSubVector(4, 4);
                for(int i = 0; i < 4; i++) {
                    if (Params[i] != ParamsMin[i])
                        throw new ApplicationException();
                    if (Params[i] != -ParamsMax[i])
                        throw new ApplicationException();
                }

                Console.WriteLine("MultiplyTest({0},{1},{2},{3})", UseXdg, DGOrder, compressL1, compressL2);
            }

            using (var solver = new Matrix_MPItestMain() { m_UseXdg = UseXdg, m_DGorder = DGOrder }) {

                // create the test data
                // ====================

                solver.Init(null);
                solver.RunSolverMode();

                Stopwatch stw = new Stopwatch();
                stw.Reset();
                stw.Start();

                BlockMsrMatrix M = solver.OperatorMatrix;
                
                int[] Ilist1 = solver.ProblemMapping.GetSubvectorIndices(false, 0);
                int[] Ilist2 = solver.ProblemMapping.GetSubvectorIndices(false, 1);

                foreach(int i in Ilist1) {
                    Assert.IsTrue(solver.ProblemMapping.IsInLocalRange(i));
                }
                foreach(int i in Ilist2) {
                    Assert.IsTrue(solver.ProblemMapping.IsInLocalRange(i));
                }

                var Blk1 = solver.ProblemMapping.GetSubBlocking(Ilist1, csMPI.Raw._COMM.WORLD, compressL1 ? -1 : 0);
                var Blk2 = solver.ProblemMapping.GetSubBlocking(Ilist2, csMPI.Raw._COMM.WORLD, compressL2 ? -1 : 0);



                int[] Tlist1 = compressL1 ? default(int[]) : Blk1.GetOccupiedIndicesList();
                int[] Tlist2 = compressL2 ? default(int[]) : Blk2.GetOccupiedIndicesList();
                if (Tlist1 != null) {
                    Assert.AreEqual(Tlist1.Length, Ilist1.Length);
                    foreach (int i in Tlist1) {
                        Assert.IsTrue(Blk1.IsInLocalRange(i));
                    }
                }
                if (Tlist2 != null) {
                    Assert.AreEqual(Tlist2.Length, Ilist2.Length);
                    foreach (int i in Tlist2) {
                        Assert.IsTrue(Blk2.IsInLocalRange(i));
                    }
                }
                BlockMsrMatrix M11 = new BlockMsrMatrix(Blk1, Blk1);
                BlockMsrMatrix M12 = new BlockMsrMatrix(Blk1, Blk2);
                BlockMsrMatrix M21 = new BlockMsrMatrix(Blk2, Blk1);
                BlockMsrMatrix M22 = new BlockMsrMatrix(Blk2, Blk2);

                M.AccSubMatrixTo(1.0, M11, Ilist1, Tlist1, Ilist1, Tlist1);
                M.AccSubMatrixTo(1.0, M12, Ilist1, Tlist1, Ilist2, Tlist2);
                M.AccSubMatrixTo(1.0, M21, Ilist2, Tlist2, Ilist1, Tlist1);
                M.AccSubMatrixTo(1.0, M22, Ilist2, Tlist2, Ilist2, Tlist2);

                /*
                MultidimensionalArray CheckRes2 = MultidimensionalArray.Create(1, 4);
                using (var MatlabRef = new BatchmodeConnector()) {

                    MatlabRef.PutVector(Ilist1.Select(i => (double)i + 1.0).ToArray(), "Ilist1");
                    MatlabRef.PutVector(Ilist2.Select(i => (double)i + 1.0).ToArray(), "Ilist2");
                    MatlabRef.PutVector(Tlist1 == null ? Ilist1.Length.ForLoop(i => (double)i + 1.0 + Blk1.i0) : Tlist1.Select(i => (double)i + 1.0).ToArray(), "Tlist1");
                    MatlabRef.PutVector(Tlist2 == null ? Ilist2.Length.ForLoop(i => (double)i + 1.0 + Blk2.i0) : Tlist2.Select(i => (double)i + 1.0).ToArray(), "Tlist2");

                    MatlabRef.PutSparseMatrix(solver.AltOperatorMatrix, "M");

                    
                    MatlabRef.Cmd("L1 = {0};", Blk1.TotalLength);
                    MatlabRef.Cmd("L2 = {0};", Blk2.TotalLength);
                    //MatlabRef.Cmd("refM11 = sparse(L1, L1);");
                    //MatlabRef.Cmd("refM12 = sparse(L1, L2);");
                    MatlabRef.Cmd("refM21 = sparse(L2, L1);");
                    //MatlabRef.Cmd("refM22 = sparse(L2, L2);");

                    //MatlabRef.Cmd("refM11(Tlist1, Tlist1) = M(Ilist1, Ilist1);");
                    //MatlabRef.Cmd("refM12(Tlist1, Tlist2) = M(Ilist1, Ilist2);");
                    MatlabRef.Cmd("refM21(Tlist2, Tlist1) = M(Ilist2, Ilist1);");
                    //MatlabRef.Cmd("refM22(Tlist2, Tlist2) = M(Ilist2, Ilist2);");

                    //MatlabRef.Cmd("err11 = norm(refM11 - M11, inf);");
                    //MatlabRef.Cmd("err12 = norm(refM12 - M12, inf);");
                    //MatlabRef.Cmd("err21 = norm(refM21 - M21, inf);");
                    //MatlabRef.Cmd("err22 = norm(refM22 - M22, inf);");

                    MatlabRef.Cmd("CheckRes = [refM21(1339, 1321), 0.0, 1.567, 0 ];");
                    MatlabRef.GetMatrix(CheckRes2, "CheckRes");

                    MatlabRef.Execute();
                }
                */

                // test multipliation (later verified by matlab)
                BlockMsrMatrix M11xM12 = new BlockMsrMatrix(M11._RowPartitioning, M12._ColPartitioning);
                M11xM12.Acc(1.0, M12);
                BlockMsrMatrix.Multiply(M11xM12, M11, M12);

                BlockMsrMatrix M22xM21 = new BlockMsrMatrix(M22._RowPartitioning, M21._ColPartitioning);
                BlockMsrMatrix.Multiply(M22xM21, M22, M21);
                double ProdNorm = M22xM21.InfNorm();
                


                stw.Stop();
                
                //M.SaveToTextFileSparse(@"C:\tmp\M.txt");
                //M11.SaveToTextFileSparse(@"C:\tmp\M11.txt");
                //M12.SaveToTextFileSparse(@"C:\tmp\M12.txt");
                //M21.SaveToTextFileSparse(@"C:\tmp\M21.txt");
                //M22.SaveToTextFileSparse(@"C:\tmp\M22.txt");
                //M22xM21.SaveToTextFileSparse(@"C:\tmp\M22xM21.txt");

                
                using(var MatlabRef = new BatchmodeConnector()) {
                    MultidimensionalArray CheckRes = MultidimensionalArray.Create(1, 4);

                    MatlabRef.PutSparseMatrix(M11, "M11");
                    MatlabRef.PutSparseMatrix(M12, "M12");
                    MatlabRef.PutSparseMatrix(M21, "M21");
                    MatlabRef.PutSparseMatrix(M22, "M22");
                    MatlabRef.PutSparseMatrix(M11xM12, "M11xM12");
                    MatlabRef.PutSparseMatrix(M22xM21, "M22xM21");

                    MatlabRef.Cmd("refM11xM12 = M12 + M11*M12;");
                    MatlabRef.Cmd("refM22xM21 = M22*M21;");

                    MatlabRef.Cmd("err1112 = norm(refM11xM12 - M11xM12, inf);");
                    MatlabRef.Cmd("err2221 = norm(refM22xM21 - M22xM21, inf);");

                    MatlabRef.Cmd("CheckRes = [err1112, err2221, 0, 0];");
                    MatlabRef.GetMatrix(CheckRes, "CheckRes");

                    MatlabRef.Execute();

                    Console.WriteLine("Matlab check M11*M12: " + CheckRes[0, 0]);
                    Console.WriteLine("Matlab check M22*M21: " + CheckRes[0, 1]);

                    Assert.IsTrue(CheckRes[0, 0] == 0.0);
                    Assert.IsTrue(CheckRes[0, 1] < 1.0e-10*ProdNorm);
                    //Assert.IsTrue(CheckRes[0, 2] == 0.0);
                    //Assert.IsTrue(CheckRes[0, 3] == 0.0);
                }

               
                Console.WriteLine("Time spend in matrix operations: " + stw.Elapsed.TotalSeconds + " sec.");

                TotTime_MatrixOp += stw.Elapsed;
            }
        }

        /*
        static void Oneify(IMutableMatrixEx M) {
            int i0 = M.RowPartitioning.i0;
            int iE = M.RowPartitioning.iE;

            int[] ci = null;
            for(int i = i0; i < iE; i++) {
                int L = M.GetOccupiedColumnIndices(i, ref ci);

                for(int l = 0; l < L; l++) {
                    M[i, ci[l]] = 1.0;
                }
            }
        }
        */

        /// <summary>
        /// Test for submatrix extraction.
        /// </summary>
        [Test]
        public static void SubMatrixTest(
            [Values(XDGusage.none, XDGusage.mixed1, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(1, 3)] int DGOrder,
            [Values(false, true)] bool compressL1,
            [Values(false, true)] bool compressL2) { 

            unsafe
            {
                int[] Params = new int[8], ParamsGlob = new int[8];
                fixed(int* pParams = Params, pParamsGlob = ParamsGlob)
                {
                    pParams[0] = (int)UseXdg;
                    pParams[1] = DGOrder;
                    pParams[2] = compressL1 ? 1 : 0;
                    pParams[3] = compressL2 ? 1 : 0;
                    pParams[4] = -pParams[0];
                    pParams[5] = -pParams[1]; 
                    pParams[6] = -pParams[2];
                    pParams[7] = -pParams[3];

                    csMPI.Raw.Allreduce((IntPtr)pParams, (IntPtr)pParamsGlob, 8, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                }

                int[] ParamsMin = ParamsGlob.GetSubVector(0, 4);
                int[] ParamsMax = ParamsGlob.GetSubVector(4, 4);
                for(int i = 0; i < 4; i++) {
                    if (Params[i] != ParamsMin[i])
                        throw new ApplicationException();
                    if (Params[i] != -ParamsMax[i])
                        throw new ApplicationException();
                }

                Console.WriteLine("SubMatrixTest({0},{1},{2},{3})", UseXdg, DGOrder, compressL1, compressL2);
            }

            using (var solver = new Matrix_MPItestMain() { m_UseXdg = UseXdg, m_DGorder = DGOrder }) {

                // create the test data
                // ====================


                solver.Init(null);
                solver.RunSolverMode();

                Stopwatch stw = new Stopwatch();
                stw.Reset();
                stw.Start();

                BlockMsrMatrix M = solver.OperatorMatrix;
                
                int[] Ilist1 = solver.ProblemMapping.GetSubvectorIndices(false, 0);
                int[] Ilist2 = solver.ProblemMapping.GetSubvectorIndices(false, 1);

                foreach(int i in Ilist1) {
                    Assert.IsTrue(solver.ProblemMapping.IsInLocalRange(i));
                }
                foreach(int i in Ilist2) {
                    Assert.IsTrue(solver.ProblemMapping.IsInLocalRange(i));
                }

                var Blk1 = solver.ProblemMapping.GetSubBlocking(Ilist1, csMPI.Raw._COMM.WORLD, compressL1 ? -1 : 0);
                var Blk2 = solver.ProblemMapping.GetSubBlocking(Ilist2, csMPI.Raw._COMM.WORLD, compressL2 ? -1 : 0);



                int[] Tlist1 = compressL1 ? default(int[]) : Blk1.GetOccupiedIndicesList();
                int[] Tlist2 = compressL2 ? default(int[]) : Blk2.GetOccupiedIndicesList();
                if (Tlist1 != null) {
                    Assert.AreEqual(Tlist1.Length, Ilist1.Length);
                    foreach (int i in Tlist1) {
                        Assert.IsTrue(Blk1.IsInLocalRange(i));
                    }
                }
                if (Tlist2 != null) {
                    Assert.AreEqual(Tlist2.Length, Ilist2.Length);
                    foreach (int i in Tlist2) {
                        Assert.IsTrue(Blk2.IsInLocalRange(i));
                    }
                }
                BlockMsrMatrix M11 = new BlockMsrMatrix(Blk1, Blk1);
                BlockMsrMatrix M12 = new BlockMsrMatrix(Blk1, Blk2);
                BlockMsrMatrix M21 = new BlockMsrMatrix(Blk2, Blk1);
                BlockMsrMatrix M22 = new BlockMsrMatrix(Blk2, Blk2);

                M.AccSubMatrixTo(1.0, M11, Ilist1, Tlist1, Ilist1, Tlist1);
                M.AccSubMatrixTo(1.0, M12, Ilist1, Tlist1, Ilist2, Tlist2);
                M.AccSubMatrixTo(1.0, M21, Ilist2, Tlist2, Ilist1, Tlist1);
                M.AccSubMatrixTo(1.0, M22, Ilist2, Tlist2, Ilist2, Tlist2);

                BlockMsrMatrix restored_M = new BlockMsrMatrix(M._RowPartitioning, M._ColPartitioning);
                int[] Idx1 = compressL1 ? Blk1.LocalLength.ForLoop(i => i + Blk1.i0) : Tlist1;
                int[] Idx2 = compressL2 ? Blk2.LocalLength.ForLoop(i => i + Blk2.i0) : Tlist2;
                M11.AccSubMatrixTo(1.0, restored_M, Idx1, Ilist1, Idx1, Ilist1);
                M12.AccSubMatrixTo(1.0, restored_M, Idx1, Ilist1, Idx2, Ilist2);
                M21.AccSubMatrixTo(1.0, restored_M, Idx2, Ilist2, Idx1, Ilist1);
                M22.AccSubMatrixTo(1.0, restored_M, Idx2, Ilist2, Idx2, Ilist2);

                // test transpose-operator
                var M_TT = M.Transpose().Transpose();
                var M11_TT = M11.Transpose().Transpose();
                var M12_TT = M12.Transpose().Transpose();
                var M21_TT = M21.Transpose().Transpose();
                var M22_TT = M22.Transpose().Transpose();
                M_TT.Acc(-1.0, M);
                M11_TT.Acc(-1.0, M11);
                M12_TT.Acc(-1.0, M12);
                M21_TT.Acc(-1.0, M21);
                M22_TT.Acc(-1.0, M22);
                double M_TT_norm = M_TT.InfNorm();
                double M11_TT_norm = M11_TT.InfNorm();
                double M12_TT_norm = M12_TT.InfNorm();
                double M21_TT_norm = M21_TT.InfNorm();
                double M22_TT_norm = M22_TT.InfNorm();
                Assert.IsTrue(M_TT_norm == 0.0, "Transpose^2 is not identity.");
                Assert.IsTrue(M11_TT_norm == 0.0, "Transpose^2 is not identity.");
                Assert.IsTrue(M12_TT_norm == 0.0, "Transpose^2 is not identity.");
                Assert.IsTrue(M21_TT_norm == 0.0, "Transpose^2 is not identity.");
                Assert.IsTrue(M22_TT_norm == 0.0, "Transpose^2 is not identity.");

                //M.SaveToTextFileSparse(@"C:\tmp\M.txt");
                //M11.SaveToTextFileSparse(@"C:\tmp\M11.txt");
                //M12.SaveToTextFileSparse(@"C:\tmp\M12.txt");
                //M21.SaveToTextFileSparse(@"C:\tmp\M21.txt");
                //M22.SaveToTextFileSparse(@"C:\tmp\M22.txt");
                //restored_M.SaveToTextFileSparse(@"C:\tmp\Mr.txt");

                stw.Stop();
                
                using(var MatlabRef = new BatchmodeConnector()) {
                    MatlabRef.PutVector(Ilist1.Select(i => (double)i + 1.0).ToArray(), "Ilist1");
                    MatlabRef.PutVector(Ilist2.Select(i => (double)i + 1.0).ToArray(), "Ilist2");
                    MatlabRef.PutVector(Tlist1 == null ? Ilist1.Length.ForLoop(i => (double)i + 1.0 + Blk1.i0) : Tlist1.Select(i => (double)i + 1.0).ToArray(), "Tlist1");
                    MatlabRef.PutVector(Tlist2 == null ? Ilist2.Length.ForLoop(i => (double)i + 1.0 + Blk2.i0) : Tlist2.Select(i => (double)i + 1.0).ToArray(), "Tlist2");

                    MultidimensionalArray CheckRes = MultidimensionalArray.Create(1, 4);

                    MatlabRef.PutSparseMatrix(M, "M");
                    MatlabRef.PutSparseMatrix(M11, "M11");
                    MatlabRef.PutSparseMatrix(M12, "M12");
                    MatlabRef.PutSparseMatrix(M21, "M21");
                    MatlabRef.PutSparseMatrix(M22, "M22");

                    MatlabRef.Cmd("L1 = {0};", Blk1.TotalLength);
                    MatlabRef.Cmd("L2 = {0};", Blk2.TotalLength);
                    MatlabRef.Cmd("refM11 = sparse(L1, L1);");
                    MatlabRef.Cmd("refM12 = sparse(L1, L2);");
                    MatlabRef.Cmd("refM21 = sparse(L2, L1);");
                    MatlabRef.Cmd("refM22 = sparse(L2, L2);");

                    MatlabRef.Cmd("refM11(Tlist1, Tlist1) = M(Ilist1, Ilist1);");
                    MatlabRef.Cmd("refM12(Tlist1, Tlist2) = M(Ilist1, Ilist2);");
                    MatlabRef.Cmd("refM21(Tlist2, Tlist1) = M(Ilist2, Ilist1);");
                    MatlabRef.Cmd("refM22(Tlist2, Tlist2) = M(Ilist2, Ilist2);");

                    MatlabRef.Cmd("err11 = norm(refM11 - M11, inf);");
                    MatlabRef.Cmd("err12 = norm(refM12 - M12, inf);");
                    MatlabRef.Cmd("err21 = norm(refM21 - M21, inf);");
                    MatlabRef.Cmd("err22 = norm(refM22 - M22, inf);");

                    MatlabRef.Cmd("CheckRes = [err11, err12, err21, err22];");
                    MatlabRef.GetMatrix(CheckRes, "CheckRes");

                    MatlabRef.Execute();

                    Console.WriteLine("Matlab check 11: " + CheckRes[0, 0]);
                    Console.WriteLine("Matlab check 12: " + CheckRes[0, 1]);
                    Console.WriteLine("Matlab check 21: " + CheckRes[0, 2]);
                    Console.WriteLine("Matlab check 22: " + CheckRes[0, 3]);

                    Assert.IsTrue(CheckRes[0, 0] == 0.0);
                    Assert.IsTrue(CheckRes[0, 1] == 0.0);
                    Assert.IsTrue(CheckRes[0, 2] == 0.0);
                    Assert.IsTrue(CheckRes[0, 3] == 0.0);
                }

                stw.Start();

                restored_M.Acc(-1.0, M);
                double err = restored_M.InfNorm();
                Console.WriteLine("Submatrix operations error: " + err);
                Assert.IsTrue(err == 0.0);

                restored_M.Clear();
                restored_M.Acc(1.0, M);
                IMutuableMatrixEx_Extensions.Acc(restored_M, -1.0, M);
                double err2 = restored_M.InfNorm();
                Console.WriteLine("Submatrix operations error: " + err2);
                Assert.IsTrue(err2 == 0.0);

                stw.Stop();

                Console.WriteLine("Time spend in matrix operations: " + stw.Elapsed.TotalSeconds + " sec.");

                TotTime_MatrixOp += stw.Elapsed;
            }
        }

        static TimeSpan TotTime_MatrixOp = new TimeSpan(0);


        /// <summary>
        /// Test for matrix/vector multiplication.
        /// </summary>
        [Test]
        public static void SpMVTest(
            [Values(XDGusage.none, XDGusage.mixed1, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(1, 3)] int DGOrder,
            [Values(false, true)] bool compressL1,
            [Values(false, true)] bool compressL2) {

            unsafe
            {
                int[] Params = new int[8], ParamsGlob = new int[8];
                fixed (int* pParams = Params, pParamsGlob = ParamsGlob)
                {
                    pParams[0] = (int)UseXdg;
                    pParams[1] = DGOrder;
                    pParams[2] = compressL1 ? 1 : 0;
                    pParams[3] = compressL2 ? 1 : 0;
                    pParams[4] = -pParams[0];
                    pParams[5] = -pParams[1];
                    pParams[6] = -pParams[2];
                    pParams[7] = -pParams[3];

                    csMPI.Raw.Allreduce((IntPtr)pParams, (IntPtr)pParamsGlob, 8, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                }

                int[] ParamsMin = ParamsGlob.GetSubVector(0, 4);
                int[] ParamsMax = ParamsGlob.GetSubVector(4, 4);
                for (int i = 0; i < 4; i++) {
                    if (Params[i] != ParamsMin[i])
                        throw new ApplicationException();
                    if (Params[i] != -ParamsMax[i])
                        throw new ApplicationException();
                }

                Console.WriteLine("SpMVTest({0},{1},{2},{3})", UseXdg, DGOrder, compressL1, compressL2);
            }

            using (var solver = new Matrix_MPItestMain() { m_UseXdg = UseXdg, m_DGorder = DGOrder }) {

                // create the test data
                // ====================

                solver.Init(null);
                solver.RunSolverMode();

                Stopwatch stw = new Stopwatch();
                stw.Reset();

                BlockMsrMatrix M = solver.OperatorMatrix;
                double[] B = new double[M.RowPartitioning.LocalLength];
                double[] X = new double[M.ColPartition.LocalLength];

                Random R = new Random();
                for (int i = 0; i < X.Length; i++)
                    X[i] = R.NextDouble();
                for (int i = 0; i < B.Length; i++)
                    B[i] = R.NextDouble();

                double[] Bb4 = B.CloneAs();

                double RefNorm = B.L2NormPow2().MPISum().Sqrt()*1e-10;

                stw.Start();

                M.SpMV(1.6, X, 0.5, B);
                
                stw.Stop();

                //M.SaveToTextFileSparse(@"C:\tmp\M.txt");
                //M11.SaveToTextFileSparse(@"C:\tmp\M11.txt");
                //M12.SaveToTextFileSparse(@"C:\tmp\M12.txt");
                //M21.SaveToTextFileSparse(@"C:\tmp\M21.txt");
                //M22.SaveToTextFileSparse(@"C:\tmp\M22.txt");
                //M22xM21.SaveToTextFileSparse(@"C:\tmp\M22xM21.txt");


                using (var MatlabRef = new BatchmodeConnector()) {
                    MultidimensionalArray CheckRes = MultidimensionalArray.Create(1, 1);

                    MatlabRef.PutSparseMatrix(M, "M");
                    MatlabRef.PutVector(Bb4, "Bref");
                    MatlabRef.PutVector(B, "B");
                    MatlabRef.PutVector(X, "X");


                    MatlabRef.Cmd("Bref = Bref*0.5 + M*X*1.6;");
                    MatlabRef.Cmd("errB = norm(B - Bref, 2);");

                    MatlabRef.Cmd("CheckRes = [errB];");
                    MatlabRef.GetMatrix(CheckRes, "CheckRes");
   


                    MatlabRef.Execute();

                    Console.WriteLine("Matlab check SpMV: " + CheckRes[0, 0]);


                    Assert.LessOrEqual(CheckRes[0, 0], RefNorm, "Error in SpMV");
                }


                Console.WriteLine("Time spend in matrix operations: " + stw.Elapsed.TotalSeconds + " sec.");

                TotTime_MatrixOp += stw.Elapsed;
            }
        }


       
    }
}

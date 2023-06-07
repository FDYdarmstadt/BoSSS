using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ilPSP {

    /// <summary>
    /// Various tests for Linear Algebra functionality in <see cref="IMatrixExtensions"/>
    /// </summary>
    [TestFixture]
    public static class MatrixOp_Tests {


        [Test]
        static public void GEMMTest_Frame([Values(1, 0.5)] double alpha, [Values(1, 2, 4)] int mult, [Values(false, true)] bool transpA, [Values(false, true)] bool transpB) {
            //            ilPSP.MatrixOp_Tests.GEMMTest

            int N = 5 * mult, M = 6 * mult, L = 13 * mult;

            var Abig = MultidimensionalArray.Create(N*4, L*3); Abig.Storage.FillRandom();
            var Bbig = MultidimensionalArray.Create(L*5, M*6); Bbig.Storage.FillRandom();

            var Cbig = MultidimensionalArray.Create(L*mult, L*mult); Cbig.Storage.FillRandom();

            for (int pass = 0; pass < 1; pass++) {

                int offsetA = N;
                int rowsA = !transpA ? N : L;
                int colsA = !transpA ? L : N;
                var A = Abig.ExtractSubArrayShallow(new int[] { offsetA, offsetA }, new int[] { offsetA + rowsA - 1, offsetA + colsA - 1 });

                int offsetB = M;
                int rowsB = !transpB ? L : M;
                int colsB = !transpB ? M : L;
                var B = Bbig.ExtractSubArrayShallow(new int[] { offsetB, offsetB }, new int[] { offsetB + rowsB - 1, offsetB + colsB - 1 });


                int offsetC = 2;
                int rowsC = !transpA ? A.NoOfRows : A.NoOfCols;
                int colsC = !transpB ? B.NoOfCols : B.NoOfRows;
                var C = Cbig.ExtractSubArrayShallow(new int[] { offsetC, offsetC }, new int[] { offsetC + rowsC - 1, offsetC + colsC - 1 });

                var D = MultidimensionalArray.Create(rowsC, colsC); D.Acc(1.0, C);
                var E = MultidimensionalArray.Create(rowsC, colsC); E.Acc(1.0, C);

                C.GEMM(alpha, A, B, 0.5, transpA, transpB);
                
                D.RefGEMM1(alpha, A, B, 0.5, transpA, transpB);
                RefGEMM2(E, alpha, A, B, 0.5, transpA, transpB);

                D.Acc(-1, C);
                E.Acc(-1, C);

                Assert.LessOrEqual(D.InfNorm(), 1.0e-12, "Error in optimized GEMM (1)");
                Assert.LessOrEqual(E.InfNorm(), 1.0e-12, "Error in optimized GEMM (2)");

                Console.WriteLine($"{(transpA ? 't' : '0')}{(transpB ? 't' : '0')} error: {D.InfNorm()} error2: {E.InfNorm()}");
            }
        }



        [Test]
        static public void GEMMTest([Values(1, 0.5)] double alpha, [Values(0.0, 0.5, 1.0)] double beta, [Values(1, 2, 4)] int mult, [Values(false, true)] bool transpA, [Values(false, true)] bool transpB) {
            //            ilPSP.MatrixOp_Tests.GEMMTest
            for (int pass = 0; pass < 2; pass++) {
                Stopwatch stw = new Stopwatch();

                int N = 9 * mult, M = 6 * mult, L = 16 * mult;
                var A = MultidimensionalArray.Create(N, L); A.Storage.FillRandom(); A = transpA ? A.TransposeTo() : A;

                var B = MultidimensionalArray.Create(L, M); B.Storage.FillRandom(); B = transpB ? B.TransposeTo() : B;

                var C = MultidimensionalArray.Create(!transpA ? A.NoOfRows : A.NoOfCols, !transpB ? B.NoOfCols : B.NoOfRows); C.Storage.FillRandom();
                var D = C.CloneAs();
                var E = C.CloneAs();

                stw.Reset();
                stw.Start();
                C.GEMM(alpha, A, B, beta, transpA, transpB);
                stw.Stop();
                double timeINTR = stw.Elapsed.TotalSeconds;

                stw.Reset();
                stw.Start();
                D.RefGEMM1(alpha, A, B, beta, transpA, transpB);
                stw.Stop();
                double timeBLAS = stw.Elapsed.TotalSeconds;

                RefGEMM2(E, alpha, A, B, beta, transpA, transpB);

                D.Acc(-1, C);
                E.Acc(-1, C);

                Assert.LessOrEqual(D.InfNorm(), 1.0e-12, "Error in optimized GEMM (1)");
                Assert.LessOrEqual(E.InfNorm(), 1.0e-12, "Error in optimized GEMM (2)");

                if (pass > 0)
                    //Console.WriteLine($"mult = {mult}, C-Size = {C.NoOfCols*C.NoOfRows} INTR = {timeINTR} , BLAS = {timeBLAS}, Factor = {timeINTR/timeBLAS}");
                    Console.WriteLine($"{(transpA ? 't' : '0')}{(transpB ? 't' : '0')} error: {D.InfNorm()} error2: {E.InfNorm()}");
            }
        }

        /// <summary>
        /// General matrix/matrix multiplication, 
        /// based on Level 3 Blas routine `dgemm`:
        /// 
        /// <paramref name="C"/> = <paramref name="alpha"/>*<paramref name="A"/>*<paramref name="B"/> + <paramref name="beta"/>*<paramref name="C"/>;
        /// 
        /// Matrix A or B should be used as transpose by setting <paramref name="transA"/> and <paramref name="transB"/>
        /// </summary>
        static void RefGEMM2(this MultidimensionalArray C, double alpha, MultidimensionalArray A, MultidimensionalArray B, double beta, bool transA = false, bool transB = false)        {
            if (!transA && !transB) {
                C.Multiply(alpha, A, B, beta, ref GEMMnn_Prog);
            } else if (transA && !transB) {
                C.Multiply(alpha, A, B, beta, ref GEMMtn_Prog);
            } else if (!transA && transB) {
                C.Multiply(alpha, A, B, beta, ref GEMMnt_Prog);
            } else if (transA && transB) {
                C.Multiply(alpha, A, B, beta, ref GEMMtt_Prog);
            }

        }

        static MultidimensionalArray.MultiplyProgram GEMMnn_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ik", "kj"); // A*B
        static MultidimensionalArray.MultiplyProgram GEMMtn_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ki", "kj"); // A^T * B
        static MultidimensionalArray.MultiplyProgram GEMMnt_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ik", "jk"); // A   * B^T
        static MultidimensionalArray.MultiplyProgram GEMMtt_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ki", "jk"); // A^T * B^T

        /// <summary>
        /// General matrix/matrix multiplication, 
        /// based on Level 3 Blas routine `dgemm`:
        /// 
        /// <paramref name="C"/> = <paramref name="alpha"/>*<paramref name="A"/>*<paramref name="B"/> + <paramref name="beta"/>*<paramref name="C"/>;
        /// 
        /// Matrix A or B should be used as transpose by setting <paramref name="transA"/> and <paramref name="transB"/>
        /// </summary>
        static public void RefGEMM1<Matrix1, Matrix2, Matrix3>(this Matrix1 C, double alpha, Matrix2 A, Matrix3 B, double beta, bool transA = false, bool transB = false)
            where Matrix1 : IMatrix
            where Matrix2 : IMatrix
            where Matrix3 : IMatrix //
        {
            if (!transA && !transB) {
                if (A.NoOfCols != B.NoOfRows)
                    throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
                if (A.NoOfRows != C.NoOfRows)
                    throw new ArgumentException("A.NoOfRows != C.NoOfRows", "A,C");
                if (B.NoOfCols != C.NoOfCols)
                    throw new ArgumentException("B.NoOfCols != C.NoOfCols", "B,C");
            } else if (transA && !transB) {
                if (A.NoOfRows != B.NoOfRows)
                    throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
                if (A.NoOfCols != C.NoOfRows)
                    throw new ArgumentException("A.NoOfRows != C.NoOfRows", "A,C");
                if (B.NoOfCols != C.NoOfCols)
                    throw new ArgumentException("B.NoOfCols != C.NoOfCols", "B,C");
            } else if (!transA && transB) {
                if (A.NoOfCols != B.NoOfCols)
                    throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
                if (A.NoOfRows != C.NoOfRows)
                    throw new ArgumentException("A.NoOfRows != C.NoOfRows", "A,C");
                if (B.NoOfRows != C.NoOfCols)
                    throw new ArgumentException("B.NoOfCols != C.NoOfCols", "B,C");
            } else if (transA && transB) {
                if (A.NoOfRows != B.NoOfCols)
                    throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
                if (A.NoOfCols != C.NoOfRows)
                    throw new ArgumentException("A.NoOfRows != C.NoOfRows", "A,C");
                if (B.NoOfRows != C.NoOfCols)
                    throw new ArgumentException("B.NoOfCols != C.NoOfCols", "B,C");
            }

            unsafe {
                int TRANSA = transA ? 't' : 'n';
                int TRANSB = transB ? 't' : 'n';

                int M = transA ? A.NoOfCols : A.NoOfRows;
                int N = transB ? B.NoOfRows : B.NoOfCols;
                int K = transA ? A.NoOfRows : A.NoOfCols;

                int LDA = transA ? Math.Max(1, K) : Math.Max(1, M);
                int LDB = transB ? Math.Max(1, N) : Math.Max(1, K);
                int LDC = Math.Max(1, M);

                int i0, i1, i2;
                double[] __A = TempBuffer.GetTempBuffer(out i0, M * K);
                double[] __B = TempBuffer.GetTempBuffer(out i1, N * K);
                double[] __C = TempBuffer.GetTempBuffer(out i2, M * N);
                fixed (double* _A = __A, _B = __B, _C = __C) {
                    IMatrixExtensions.CopyToUnsafeBuffer(A, _A, true);
                    IMatrixExtensions.CopyToUnsafeBuffer(B, _B, true);
                    IMatrixExtensions.CopyToUnsafeBuffer(C, _C, true);

                    BLAS.dgemm(TRANSA, TRANSB, M, N, K, alpha, _A, LDA, _B, LDB, beta, _C, LDC);

                    IMatrixExtensions.CopyFromUnsafeBuffer(C, _C, true);
                }
                TempBuffer.FreeTempBuffer(i0);
                TempBuffer.FreeTempBuffer(i1);
                TempBuffer.FreeTempBuffer(i2);
            }
        }
    }
}

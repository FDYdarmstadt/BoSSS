using ilPSP;
using ilPSP.Utils;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AdvancedSolverTests.Solver {
    public class mklILU {
        public class testilu : ilPSP.LinSolvers.ILU.ILUSolver {
            public double[] _ForwardSubstitution(double[] rhs) {
                return ForwardSubstitution(rhs);
            }

            public double[] _BackwardSubstitution(double[] rhs) {
                return BackwardSubstitution(rhs);
            }
        } 


        private static (MultidimensionalArray L_matlab, MultidimensionalArray U_matlab, MultidimensionalArray A) GenMatrix() {

            int basis = 5;
            int size = (int)Math.Floor(Math.Pow(basis, 2));
            var L_matlab = MultidimensionalArray.Create(size, size);
            var U_matlab = MultidimensionalArray.Create(size, size);
            var A = MultidimensionalArray.Create(size, size);

            using (BatchmodeConnector matlab = new BatchmodeConnector()) {
                //note: BatchmodeCon maybe working on proc0 but savetotxt file, etc. (I/O) is full mpi parallel
                //so concider this as full mpi-parallel

                matlab.Cmd($"A = gallery('neumann', {size}) + speye({size});");
                matlab.Cmd($"options.type = 'nofill';");
                matlab.Cmd($"[L, U] = ilu(A, options);");
                matlab.Cmd($"Lfull = full(L);");
                matlab.Cmd($"Ufull = full(U);");
                matlab.Cmd($"Afull = full(A);");
                matlab.GetMatrix(A, "Afull");
                matlab.GetMatrix(L_matlab, "Lfull");
                matlab.GetMatrix(U_matlab, "Ufull");
                matlab.Execute();
            }
            Console.WriteLine(L_matlab.InfNorm());
            return (L_matlab, U_matlab, A);
        }

        private static MsrMatrix ExecuteILU(MultidimensionalArray A) {
            var ilu = new testilu();
            var A_BMsr = new MsrMatrix(A.Lengths[0],A.Lengths[1]);
            A_BMsr.AccDenseMatrix(1.0, A);
            ilu.DefineMatrix(A_BMsr);
            var LU=ilu.GetILUFactorization;
            return LU;
        }

        [Test]
        public static void CompareFactorization() {
            var matrices=GenMatrix();
            var mklILU = ExecuteILU(matrices.A);
            Console.WriteLine($"dim mkl ILU: {mklILU.RowPartitioning.LocalLength}, dim matlab L: {matrices.L_matlab.Lengths[0]}, dim matlab U: {matrices.U_matlab.Lengths[0]}");
            //matrices.L_matlab.SaveToTextFile("Lmatlab");
            //matrices.U_matlab.SaveToTextFile("Umatlab");
            //ILU.SaveToTextFileSparse("ILU");
            mklILU.AccBlock(0,0,-1.0, matrices.L_matlab);
            mklILU.AccBlock(0,0,-1.0, matrices.U_matlab);
            Assert.IsTrue(mklILU.InfNorm()-1<1E-14);
        }

        [Test]
        public static void CompareSubstitution() {
            // matlab part
            var matrices = GenMatrix();
            int L = matrices.A.GetLength(0);
            var rnd = new Random(0);
            var rndB = L.ForLoop(i => rnd.NextDouble());
            var Ycheck = ForwardSubstitution(matrices.L_matlab, rndB);
            var Xcheck = BackwardSubstitution(matrices.U_matlab, Ycheck);

            // mkl part
            var ilu = new testilu();
            var A_BMsr = new MsrMatrix(matrices.A.GetLength(0), matrices.A.GetLength(1));
            A_BMsr.AccDenseMatrix(1.0, matrices.A);
            ilu.DefineMatrix(A_BMsr);
            var Y = ilu._ForwardSubstitution(rndB);
            var X = ilu._BackwardSubstitution(Y);

            Ycheck.AccV(-1.0, Y);
            Xcheck.AccV(-1.0, X);
            Assert.IsTrue(Ycheck.L2Norm() < 1E-14);
            Assert.IsTrue(Xcheck.L2Norm() < 1E-14);
        }

        private static double[] ForwardSubstitution(MultidimensionalArray L, double[] RHS) {
            var X = new double[RHS.Length];
            X.SetV(RHS, 1.0);
            int N = L.GetLength(0);
            for(int iRow = 0; iRow < N; iRow++) {
                double tmp = 0.0;
                for(int j = 0; j < iRow; j++) {
                    tmp += X[j]*L[iRow,j] ;
                }
                X[iRow] = (X[iRow] - tmp) / L[iRow, iRow];
            }
            return X;
        }


        private static double[] BackwardSubstitution(MultidimensionalArray U, double[] RHS) {
            var X = new double[RHS.Length];
            X.SetV(RHS, 1.0);
            int N = U.GetLength(0);
            for (int iRow = N-1; iRow >= 0; iRow--) {
                double tmp = 0.0;
                for (int j = iRow + 1; j < N; j++) {
                    tmp += X[j] * U[iRow, j];
                }
                X[iRow] = (X[iRow] - tmp) / U[iRow, iRow];
            }
            return X;
        }

    }
}

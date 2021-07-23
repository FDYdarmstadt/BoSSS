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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Kraypis;

namespace BoSSS.Foundation.ConstrainedDGprojection {

    public class myCG : IDisposable {
        public void Init(BlockMsrMatrix M) {
            m_Matrix = M;
            //var M_test = M.CloneAs();
            //M_test.Acc(-1.0, M.Transpose());
            //Console.WriteLine("Symm-test: " + M_test.InfNorm());
            //Console.WriteLine("Inf-Norm: " + M.InfNorm());
            //m_Matrix.SaveToTextFileSparse("M");
            PrecondInit();
        }

        private MsrMatrix ILU_M;
        public void GetLocalMatrix() {
            int rank;
            MPI.Wrappers.csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out rank);

            // extraction of local Matrix block, ILU will be executed process local
            var LocBlocki0s = new List<long>();
            var LocBlockLen = new List<int>();
            long IdxOffset = m_Matrix._RowPartitioning.i0;
            for (int i = 0; i < m_Matrix._RowPartitioning.LocalNoOfBlocks; i++) {
                long iBlock = i + m_Matrix._RowPartitioning.FirstBlock;
                long i0 = m_Matrix._RowPartitioning.GetBlockI0(iBlock) - IdxOffset;
                int Len = m_Matrix._RowPartitioning.GetBlockLen(iBlock);
                LocBlocki0s.Add(i0);
                LocBlockLen.Add(Len);
            }
            long[] RowISrc = m_Matrix._RowPartitioning.LocalLength.ForLoop(i => i + IdxOffset);

            var part = new BlockPartitioning(m_Matrix._RowPartitioning.LocalLength, LocBlocki0s, LocBlockLen, csMPI.Raw._COMM.SELF);
            BlockMsrMatrix localMatrix = new BlockMsrMatrix(part);

            m_Matrix.WriteSubMatrixTo(localMatrix, RowISrc, default(long[]), RowISrc, default(long[]));
            ILU_M = localMatrix.ToMsrMatrix();

            //ILU_M.Multicondest(printout:true);
            //Console.WriteLine("cond num:"+m_Matrix.condest());
        }

        ISparseSolver dirsolver;

        private void DirectInit() {
            dirsolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver() {
                CacheFactorization = true,
                UseDoublePrecision = false
            };
            dirsolver.DefineMatrix(ILU_M);
        }

        private void DirectSolve<P, Q>(P X, Q B)
            where P : IList<double>
            where Q : IList<double> {
            if (ILU_M.RowPartitioning.LocalLength == 0 && ILU_M.RowPartitioning.MPI_Comm == csMPI.Raw._COMM.SELF)
                return; // there is nothing to do and we can safely skip the solve step
            dirsolver.Solve(X, B);
        }

        private void ILUDecomposition() {


            if (this.ILU_M.RowPartitioning.MpiSize != 1)
                throw new NotSupportedException();
            long n = ILU_M.RowPartitioning.LocalLength;

            // Zeros on diagonal elements because of saddle point structure
            for (int bla = 0; bla < n; bla++) {
                if (ILU_M.GetDiagonalElement(bla) == 0)
                    ILU_M.SetDiagonalElement(bla, 1);
            }

            Func<double, bool> ZeroPattern = delegate (double e) {
                return e == 0;
                //return false;
            };

            // ++++++++++++++++++++
            // ILU(0) decomposition
            // ++++++++++++++++++++

            for (int k = 0; k < n - 1; k++) {
                for (int i = k + 1; i < n; i++) {
                    if (ZeroPattern(ILU_M[i, k])) continue;
                    ILU_M[i, k] = ILU_M[i, k] / ILU_M[k, k];
                    for (int j = k + 1; j < n; j++) {
                        if (ZeroPattern(ILU_M[i, j])) continue;
                        ILU_M[i, j] = ILU_M[i, j] - ILU_M[i, k] * ILU_M[k, j];
                    }
                }
            }

            //var part = ILU_M._RowPartitioning;
            //var L = new BlockMsrMatrix(part);
            //var U = new BlockMsrMatrix(part);

            //for (int i = 0; i < n; i++) {
            //    for (int j = 0; j < i; j++) {
            //        L[i, j] = ILU_M[i, j];
            //        U[j, i] = ILU_M[j, i];
            //    }
            //    L[i, i] = 1;
            //    U[i, i] = ILU_M[i, i];
            //}
        }

        private void ICholDecomposition() {

            if (this.ILU_M.RowPartitioning.MpiSize != 1)
                throw new NotSupportedException();
            long n = ILU_M.RowPartitioning.LocalLength;

            // Zeros on diagonal elements because of saddle point structure
            for (int bla = 0; bla < n; bla++) {
                if (ILU_M.GetDiagonalElement(bla) == 0)
                    ILU_M.SetDiagonalElement(bla, 1);
            }

            Func<double, bool> ZeroPattern = delegate (double e) {
                return e == 0;
                //return false;
            };

            // ++++++++++++++++++++
            // ICHOL(0) decomposition
            // ++++++++++++++++++++

            for (int k = 0; k < n; k++) {
                ILU_M[k, k] = Math.Sqrt(ILU_M[k, k]);
                for (int i = k + 1; i < n; i++) {
                    if (ZeroPattern(ILU_M[i, k])) continue;
                    ILU_M[i, k] = ILU_M[i, k] / ILU_M[k, k];
                }
                for (int j = k + 1; j < n; j++) {
                    for (int i = j; i < n; i++) {
                        if (ZeroPattern(ILU_M[i, j])) continue;
                        ILU_M[i, j] = ILU_M[i, j] - ILU_M[i, k] * ILU_M[j, k];
                    }
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    ILU_M[i, j] = ILU_M[j, i];
                    if (Double.IsNaN(ILU_M[i, j]) || Double.IsInfinity(ILU_M[i, j]))
                        throw new Exception("Ich habs doch gewusst");
                }
            }
        }

        public void DecompSolve<P, Q>(P X, Q B)
            where P : IList<double>
            where Q : IList<double> {

            long n = ILU_M.RowPartitioning.LocalLength;
            double buffer = 0;

            // find solution of Ly = b
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                buffer = 0;
                for (int k = 0; k < i; k++)
                    buffer += ILU_M[i, k] * y[k];
                //y[i] = (1/ ILU_M[i, i])*(B[i] - buffer);
                y[i] = (B[i] - buffer);
            }
            // find solution of Ux = y
            for (long i = n - 1; i >= 0; i--) {
                buffer = 0;
                for (long k = i + 1; k < n; k++)
                    buffer += ILU_M[i, k] * X[(int)k];
                X[(int)i] = (1 / ILU_M[i, i]) * (y[i] - buffer);
            }
        }

        public void BlockJacInit() {
            BlockMsrMatrix M = m_Matrix;

            int L = M.RowPartitioning.LocalLength;

            Diag = new BlockMsrMatrix(M._RowPartitioning, M._ColPartitioning);
            invDiag = new BlockMsrMatrix(M._RowPartitioning, M._ColPartitioning);
            int Jloc = M._RowPartitioning.LocalNoOfBlocks;
            long j0 = M._RowPartitioning.FirstBlock;
            MultidimensionalArray temp = null;
            for (int j = 0; j < Jloc; j++) {
                long jBlock = j + j0;
                int Nblk = M._RowPartitioning.GetBlockLen(jBlock);
                long i0 = M._RowPartitioning.GetBlockI0(jBlock);

                if (temp == null || temp.NoOfCols != Nblk)
                    temp = MultidimensionalArray.Create(Nblk, Nblk);

                M.ReadBlock(i0, i0, temp);
                Diag.AccBlock(i0, i0, 1.0, temp, 0.0);
                temp.InvertInPlace();
                invDiag.AccBlock(i0, i0, 1.0, temp, 0.0);
            }
        }
        private double omega = 0.5;

        public void BlockJacSolve<U, V>(U xl, V bl)
            where U : IList<double>
            where V : IList<double> {
            int L = xl.Count;
            double[] ql = new double[L];

            for (int iIter = 0; iIter < 5; iIter++) {
                if (L > 0) ql.SetV(bl);

                m_Matrix.SpMV(-1.0, xl, 1.0, ql);
                //double ResNorm = ql.L2NormPow2().MPISum().Sqrt();

                Diag.SpMV(1.0, xl, 1.0, ql);
                invDiag.SpMV(omega, ql, 1.0 - omega, xl);

            }
        }

        BlockMsrMatrix m_Matrix;
        BlockMsrMatrix Diag;
        BlockMsrMatrix invDiag;


        private void HYPRE_ILU_Init() {
            dirsolver = new ilPSP.LinSolvers.HYPRE.Euclid() {
                Level = 0,
                Comm = csMPI.Raw._COMM.SELF
            };
            if (ILU_M != null && ILU_M.RowPartitioning.LocalLength > 0) dirsolver.DefineMatrix(ILU_M);
        }


        private void PrecondInit() {
            //BlockJacInit();
            GetLocalMatrix();
            //MKL_ILU_Init();
            //DirectInit();
            //ILUDecomposition();
            //ICholDecomposition();
            //MKL_ILU_Init();
            HYPRE_ILU_Init();
        }

        private void PrecondSolve<U, V>(U xl, V bl) where U : IList<double>
            where V : IList<double> {
            //BlockJacSolve(xl,bl);
            //DecompSolve(xl, bl);
            DirectSolve(xl, bl);
        }

        public void Solve<Vec1, Vec2>(Vec1 _x, Vec2 _R)
            where Vec1 : IList<double>
            where Vec2 : IList<double> //
        {

            double[] x, R;
            if (_x is double[]) {
                x = _x as double[];
            } else {
                x = _x.ToArray();
            }
            if (_R is double[]) {
                R = _R as double[];
            } else {
                R = _R.ToArray();
            }

            int L = x.Length;


            double[] P = new double[L];
            double[] V = new double[L];
            double[] Z = new double[L];

            int[] Lengths = L.MPIGather(0);
            if (ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                //for (int i = 0; i < Lengths.Length; i++)
                //    Console.WriteLine("L from proc " + i + " : " + Lengths[i]);
            }
            // compute P0, R0
            // ==============
            if (x.Length != 0)
                GenericBlas.dswap(L, x, 1, P, 1);

            m_Matrix.SpMV(-1.0, P, 1.0, R);

            if (x.Length != 0)
                GenericBlas.dswap(L, x, 1, P, 1);

            //if (x.Length != 0) P.SetV(R);
            PrecondSolve(Z, R);
            if (x.Length != 0) P.SetV(Z);

            double alpha_loc = x.Length != 0 ? R.InnerProd(P) : 0;
            double alpha = alpha_loc.MPISum();
            double alpha_0 = alpha;
            Console.WriteLine(alpha);
            double ResNorm;

            var ResReal = _R.ToArray();
            var Xdummy = new double[R.Length];

            ResNorm = Math.Sqrt(alpha);
            double ResNorm0 = ResNorm;



            // iterate
            // =======
            for (int n = 1; true; n++) {

                if (n % 1 == 0) {

                    var theResidual = new double[R.Length];
                    if (x.Length != 0) {
                        theResidual.SetV(ResReal);
                        Xdummy.SetV(x);
                    }
                    m_Matrix.SpMV(-1.0, Xdummy, 1.0, theResidual);
                    double bla = (x.Length != 0 ? theResidual.L2NormPow2() : 0).MPISum().Sqrt();
                    Console.WriteLine("Res real at n" + n + ":" + bla);
                }

                Console.WriteLine("ResNorm at n" + n + ":" + ResNorm);
                if (ResNorm / ResNorm0 + ResNorm < 1E-6 || ResNorm < 1E-6 || n >= 100) {
                    if (n > 1000) Console.WriteLine("maximum number of iterations reached. Solution maybe not been converged.");
                    break;
                }

                if (Math.Abs(alpha) <= double.Epsilon) {
                    // numerical breakdown
                    break;
                }


                m_Matrix.SpMV(1.0, P, 0, V);
                double VxP_loc = x.Length != 0 ? V.InnerProd(P) : 0;
                double VxP = VxP_loc.MPISum();
                if (double.IsNaN(VxP) || double.IsInfinity(VxP))
                    throw new ArithmeticException();
                double lambda = alpha / VxP;
                if (double.IsNaN(lambda) || double.IsInfinity(lambda))
                    throw new ArithmeticException();

                if (x.Length != 0) {
                    x.AccV(lambda, P);
                    R.AccV(-lambda, V);
                }

                //if (x.Length != 0) Z.SetV(R);
                if (x.Length != 0) Z.Clear();
                PrecondSolve(Z, R);


                double alpha_neu_loc = x.Length != 0 ? R.InnerProd(Z) : 0;
                double alpha_neu = alpha_neu_loc.MPISum();

                // compute residual norm
                ResNorm = (x.Length != 0 ? R.L2NormPow2() : 0).MPISum().Sqrt();

                if (x.Length != 0) {
                    P.ScaleV(alpha_neu / alpha);
                    P.AccV(1.0, Z);

                    alpha = alpha_neu;
                }

            }

            if (!object.ReferenceEquals(_x, x))
                _x.SetV(x);
            if (!object.ReferenceEquals(_R, R))
                _R.SetV(R);



            return;

        }

        public void Dispose() {
            m_Matrix.Clear();
        }
    }

}

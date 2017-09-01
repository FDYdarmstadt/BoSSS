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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using BoSSS.Foundation;

namespace BoSSS.Solution.Multigrid
{
    /// <summary>
    /// Standard preconditioned GMRES.
    /// </summary>
    public class SoftGMRES : ISolverSmootherTemplate, ISolverWithCallback
    {


        public void Init(MultigridOperator op)
        {

            var Mtx = op.OperatorMatrix;
            var MgMap = op.Mapping;
            this.m_mgop = op;

            if (!Mtx.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!Mtx.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");
            this.Matrix = Mtx;
            if (Precond != null)
            {
                Precond.Init(op);
            }
        }

        MultigridOperator m_mgop;

        public ISolverSmootherTemplate Precond;

        BlockMsrMatrix Matrix;

        public double m_Tolerance = 1.0e-10;
        public int m_MaxIterations = 10000;

        public int NoOfIterations = 0;


        public Action<int, double[], double[], MultigridOperator> IterationCallback
        {
            get;
            set;
        }


        /// <summary>
        /// number of iterations between restarts
        /// </summary>
        public int MaxKrylovDim = 80;

        /// <summary>
        /// Compute the Givens rotation matrix parameters for a and b.
        /// </summary>
        static void rotmat(out double c, out double s, double a, double b)
        {
            double temp;
            if (b == 0.0)
            {
                c = 1.0;
                s = 0.0;
            }
            else if (Math.Abs(b) > Math.Abs(a))
            {
                temp = a / b;
                s = 1.0 / Math.Sqrt(1.0 + temp * temp);
                c = temp * s;
            }
            else
            {
                temp = b / a;
                c = 1.0 / Math.Sqrt(1.0 + temp * temp);
                s = temp * c;
            }
        }

        public bool m_Converged = false;

        public void Solve<V1, V2>(V1 _X, V2 _B)
            where V1 : IList<double>
            where V2 : IList<double>
        {
            double[] X, B;
            if (_X is double[])
            {
                X = _X as double[];
            }
            else
            {
                X = _X.ToArray();
            }
            if (_B is double[])
            {
                B = _B as double[];
            }
            else
            {
                B = _B.ToArray();
            }


            double bnrm2 = B.L2NormPow2().MPISum().Sqrt();
            if (bnrm2 == 0.0)
            {
                bnrm2 = 1.0;
            }

            int Nloc = Matrix.RowPartitioning.LocalLength;
            int Ntot = Matrix.RowPartitioning.TotalLength;

            double[] r = new double[Nloc];
            double[] z = new double[Nloc];

            //r = M \ ( b-A*x );, where M is the precond
            z.SetV(B);
            Matrix.SpMV(-1.0, X, 1.0, z);
            if (IterationCallback != null)
            {
                IterationCallback(0, X.CloneAs(), z.CloneAs(), this.m_mgop);
            }
            if (this.Precond != null)
            {
                r.Clear();
                this.Precond.Solve(r, z);
            }
            else
            {
                r.SetV(z);
            }

            // Inserted for real residual
            double error2 = z.L2NormPow2().MPISum().Sqrt();

            double error = (r.L2NormPow2().MPISum().Sqrt()) / bnrm2;
            if (error < this.m_Tolerance)
            {
                if (!object.ReferenceEquals(_X, X))
                    _X.SetV(X);
                B.SetV(z);
                this.m_Converged = true;
                return;
            }

            if (MaxKrylovDim <= 0)
                throw new NotSupportedException("unsupported restart length.");

            int m = MaxKrylovDim;
            double[][] V = (m + 1).ForLoop(i => new double[Nloc]); //   V(1:n,1:m+1) = zeros(n,m+1);
            MultidimensionalArray H = MultidimensionalArray.Create(m + 1, m); //   H(1:m+1,1:m) = zeros(m+1,m);

            double[] cs = new double[m];
            double[] sn = new double[m];
            //double[] s = new double[m+1];
            double[] e1 = new double[Nloc];
            //if(Matrix.RowPartitioning.Rank == 0)
            e1[0] = 1.0;

            double[] s = new double[Nloc], w = new double[Nloc], y;
            double temp;
            int iter;
            for (iter = 1; iter <= m_MaxIterations; iter++)
            { // GMRES iterations
                // r = M \ ( b-A*x );
                z.SetV(B);
                Matrix.SpMV(-1.0, X, 1.0, z);

                error2 = z.L2NormPow2().MPISum().Sqrt();

                if (this.Precond != null)
                {
                    r.Clear();
                    this.Precond.Solve(r, z);
                }
                else
                {
                    r.SetV(z);
                }

                // V(:,1) = r / norm( r );
                double norm_r = r.L2NormPow2().MPISum().Sqrt();
                V[0].SetV(r, alpha: (1.0 / norm_r));

                //s = norm( r )*e1;
                s.SetV(e1, alpha: norm_r);

                int i;

                # region Gram-Schmidt (construct orthonormal  basis using Gram-Schmidt)
                for (i = 1; i <= m; i++)
                {
                    this.NoOfIterations++;

                    #region Arnoldi procdure

                    //w = M \ (A*V(:,i));                         
                    Matrix.SpMV(1.0, V[i - 1], 0.0, z);
                    if (this.Precond != null)
                    {
                        w.Clear();
                        this.Precond.Solve(w, z);
                    }
                    else
                    {
                        w.SetV(z);
                    }

                    for (int k = 1; k <= i; k++)
                    {
                        H[k - 1, i - 1] = GenericBlas.InnerProd(w, V[k - 1]).MPISum();
                        //w = w - H(k,i)*V(:,k);
                        w.AccV(-H[k - 1, i - 1], V[k - 1]);
                    }

                    double norm_w = w.L2NormPow2().MPISum().Sqrt();
                    H[i + 1 - 1, i - 1] = norm_w; // the +1-1 actually makes me sure I haven't forgotten to subtract -1 when porting the code
                    //V(:,i+1) = w / H(i+1,i);
                    V[i + 1 - 1].SetV(w, alpha: (1.0 / norm_w));

                    #endregion

                    #region Givens rotation

                    for (int k = 1; k <= i - 1; k++)
                    {
                        // apply Givens rotation, H is Hessenbergmatrix
                        temp = cs[k - 1] * H[k - 1, i - 1] + sn[k - 1] * H[k + 1 - 1, i - 1];
                        H[k + 1 - 1, i - 1] = -sn[k - 1] * H[k - 1, i - 1] + cs[k - 1] * H[k + 1 - 1, i - 1];
                        H[k - 1, i - 1] = temp;
                    }
                    //	 [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
                    rotmat(out cs[i - 1], out sn[i - 1], H[i - 1, i - 1], H[i + 1 - 1, i - 1]);
                    temp = cs[i - 1] * s[i - 1]; //                       % approximate residual norm
                    H[i - 1, i - 1] = cs[i - 1] * H[i - 1, i - 1] + sn[i - 1] * H[i + 1 - 1, i - 1];
                    H[i + 1 - 1, i - 1] = 0.0;

                    #endregion

                    // update the residual vector (s == beta in many pseudocodes)
                    s[i + 1 - 1] = -sn[i - 1] * s[i - 1];
                    s[i - 1] = temp;
                    error = Math.Abs(s[i + 1 - 1]) / bnrm2;
                    //{
                    //    int rootRank = Matrix.RowPartitioning.FindProcess(i + 1 - 1);
                    //    if (Matrix.RowPartitioning.Rank == rootRank) {


                    //    } else {
                    //        error = double.NaN;
                    //    }
                    //    unsafe {
                    //        csMPI.Raw.Bcast((IntPtr)(&error), 1, csMPI.Raw._DATATYPE.DOUBLE, rootRank, Matrix.RowPartitioning.MPI_Comm);
                    //    }
                    //}



                    if (error <= m_Tolerance)
                    {
                        // update approximation and exit

                        //y = H(1:i,1:i) \ s(1:i);    
                        y = new double[i];
                        H.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { i - 1, i - 1 })
                            .Solve(y, s.GetSubVector(0, i));



                        // x = x + V(:,1:i)*y;
                        for (int ii = 0; ii < i; ii++)
                        {
                            X.AccV(y[ii], V[ii]);
                        }
                        this.m_Converged = true;
                        break;
                    }
                }
                #endregion
                //Debugger.Launch();


                if (error <= this.m_Tolerance)
                {
                    this.m_Converged = true;
                    break;
                }


                // y = H(1:m,1:m) \ s(1:m);
                y = new double[m];
                H.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { m - 1, m - 1 })
                    .Solve(y, s.GetSubVector(0, m));
                // update approximation: x = x + V(:,1:m)*y;  
                for (int ii = 0; ii < m; ii++)
                {
                    X.AccV(y[ii], V[ii]);
                }

                // compute residual: r = M \ ( b-A*x )     
                z.SetV(B);
                Matrix.SpMV(-1.0, X, 1.0, z);
                if (IterationCallback != null)
                {
                    error2 = z.L2NormPow2().MPISum().Sqrt();
                    IterationCallback(iter, X.CloneAs(), z.CloneAs(), this.m_mgop);
                    //IterationCallback(this.NoOfIterations, X.CloneAs(), z.CloneAs(), this.m_mgop);
                }
                if (this.Precond != null)
                {
                    r.Clear();
                    this.Precond.Solve(r, z);
                }
                else
                {
                    r.SetV(z);
                }


                norm_r = r.L2NormPow2().MPISum().Sqrt();
                s[i + 1 - 1] = norm_r;
                error = s[i + 1 - 1] / bnrm2;        // % check convergence
                if (error2 <= m_Tolerance)
                    break;
            }

            if (IterationCallback != null)
            {
                z.SetV(B);
                Matrix.SpMV(-1.0, X, 1.0, z);
                IterationCallback(iter, X.CloneAs(), z.CloneAs(), this.m_mgop);
            }


            if (!object.ReferenceEquals(_X, X))
                _X.SetV(X);
            B.SetV(z);
        }


        public int IterationsInNested
        {
            get
            {
                if (this.Precond != null)
                    return this.Precond.IterationsInNested + this.Precond.ThisLevelIterations;
                else
                    return 0;
            }
        }

        public int ThisLevelIterations
        {
            get
            {
                return this.NoOfIterations;
            }
        }

        public bool Converged
        {
            get
            {
                return m_Converged;
            }
        }

        public void ResetStat()
        {
            this.m_Converged = false;
            this.NoOfIterations = 0;
            if (this.Precond != null)
                this.Precond.ResetStat();
        }
    }
}

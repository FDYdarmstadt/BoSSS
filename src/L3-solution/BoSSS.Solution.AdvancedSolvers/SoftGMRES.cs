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
using System.IO;
using ilPSP.Tracing;

namespace BoSSS.Solution.AdvancedSolvers {
    /// <summary>
    /// Standard preconditioned GMRES.
    /// </summary>
    public class SoftGMRES : ISubsystemSolver, ISolverWithCallback, IProgrammableTermination {

        Func<int, double, double, (bool bNotTerminate, bool bSuccess)> m_TerminationCriterion;

        /// <summary>
        /// ~
        /// </summary>
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get {
                return m_TerminationCriterion;
            }
            set {
                m_TerminationCriterion = value;
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        public SoftGMRES() {
            m_TerminationCriterion = (iIter, r0, ri) => {
                return (iIter <= 500 && ri > 1e-7, ri <= 1e-7);
            };
        }


        /// <summary>
        /// ~
        /// </summary>
        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        /// <summary>
        /// ~
        /// </summary>
        public void Init(MultigridOperator op) {
            InitImpl(op);
        }

        void InitImpl(IOperatorMappingPair op) {
            using(new FuncTrace()) {
                if(object.ReferenceEquals(op, this.m_mgop))
                    return; // already initialized
                else
                    this.Dispose(); // must re-initialize

                var Mtx = op.OperatorMatrix;
                var MgMap = op.DgMapping;
                this.m_mgop = op;

                if(!Mtx.RowPartitioning.EqualsPartition(MgMap))
                    throw new ArgumentException("Row partitioning mismatch.");
                if(!Mtx.ColPartition.EqualsPartition(MgMap))
                    throw new ArgumentException("Column partitioning mismatch.");
                if(Precond != null) {
                    if(Precond is ISubsystemSolver sssol) {
                        sssol.Init(m_mgop);
                    } else if(m_mgop is MultigridOperator mgOp) {
                        Precond.Init(mgOp);
                    } else {
                        throw new NotSupportedException($"Unable to initialize preconditioner if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                    }
                }
            }
        }


        IOperatorMappingPair m_mgop;

        public ISolverSmootherTemplate Precond;


        /*
        ISparseMatrix m_Matrix;

        ISparseMatrix Matrix {
            get {
                if(m_Matrix == null) {
                    m_Matrix = m_mgop.OperatorMatrix; // activate for BlockMatrixSpMV

                    //var hypreMtx = new ilPSP.LinSolvers.HYPRE.IJMatrix(m_mgop.OperatorMatrix); // HYPRE
                    //m_Matrix = hypreMtx;

                    //var monkeyMtx = new ilPSP.LinSolvers.monkey.CPU.RefMatrix(m_mgop.OperatorMatrix.ToMsrMatrix()); // Monkey
                    //m_Matrix = monkeyMtx;

                }
                return m_Matrix; 
            }
        }
        */

        BlockMsrMatrix Matrix => m_mgop.OperatorMatrix;

        public string m_SessionPath;

        //public double m_Tolerance = 1.0e-10;
        //public int m_MaxIterations = 10000;

        private int NoOfIterations = 0;

        /// <summary>
        /// ~
        /// </summary>
        public Action<int, double[], double[], MultigridOperator> IterationCallback {
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
        static void rotmat(out double c, out double s, double a, double b) {
            double temp;
            if(b == 0.0) {
                c = 1.0;
                s = 0.0;
            } else if(Math.Abs(b) > Math.Abs(a)) {
                temp = a / b;
                s = 1.0 / Math.Sqrt(1.0 + temp * temp);
                c = temp * s;
            } else {
                temp = b / a;
                c = 1.0 / Math.Sqrt(1.0 + temp * temp);
                s = temp * c;
            }
        }

        public bool m_Converged = false;

        /// <summary>
        /// ~
        /// </summary>
        public void Solve<V1, V2>(V1 _X, V2 _B)
            where V1 : IList<double>
            where V2 : IList<double> {

            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = false;
                tr.Info("IterationCallback set? " + (this.IterationCallback != null));
                tr.Info("Preconditioner is " + (this.Precond?.ToString() ?? "null"));
                double[] X, B;
                if(_X is double[]) {
                    X = _X as double[];
                } else {
                    X = _X.ToArray();
                }
                if(_B is double[]) {
                    B = _B as double[];
                } else {
                    B = _B.ToArray();
                }


                double bnrm2 = B.MPI_L2Norm(Matrix.MPI_Comm);
                if(bnrm2 == 0.0) {
                    bnrm2 = 1.0;
                }

                int Nloc = Matrix.RowPartitioning.LocalLength;
                long Ntot = Matrix.RowPartitioning.TotalLength;

                double[] r = new double[Nloc];
                double[] z = new double[Nloc];

                //r = M \ ( b-A*x );, where M is the precond
                z.SetV(B);
                Matrix.SpMV(-1.0, X, 1.0, z);
                IterationCallback?.Invoke(0, X.CloneAs(), z.CloneAs(), this.m_mgop as MultigridOperator);

                if(this.Precond != null) {
                    r.Clear();
                    this.Precond.Solve(r, z);
                } else {
                    r.SetV(z);
                }

                // Inserted for real residual
                double error2 = z.MPI_L2Norm(Matrix.MPI_Comm);
                double iter0_error2 = error2;

                double error = (r.L2NormPow2().MPISum(Matrix.MPI_Comm).Sqrt()) / bnrm2;
                var term0 = TerminationCriterion(0, error2, error2);
                if(!term0.bNotTerminate) {

                    if(!object.ReferenceEquals(_X, X))
                        _X.SetV(X);
                    B.SetV(z);
                    tr.Info($"convergence reached before iterations, success? {term0.bSuccess}: iter0_err = {iter0_error2}");
                    //Console.WriteLine("convergence reached (1);");
                    this.m_Converged = term0.bSuccess;
                    return;
                }

                if(MaxKrylovDim <= 0)
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
                int totIterCounter = 0;
                for(iter = 1; true; iter++) { // GMRES iterations
                                              // r = M \ ( b-A*x );
                    z.SetV(B);
                    Matrix.SpMV(-1.0, X, 1.0, z);

                    error2 = z.MPI_L2Norm(Matrix.MPI_Comm);

                    if(this.Precond != null) {
                        r.Clear();
                        this.Precond.Solve(r, z);
                    } else {
                        r.SetV(z);
                    }

                    // V(:,1) = r / norm( r );
                    double norm_r = r.MPI_L2Norm(Matrix.MPI_Comm);
                    V[0].SetV(r, alpha: (1.0 / norm_r));

                    //s = norm( r )*e1;
                    s.SetV(e1, alpha: norm_r);

                    int i;

                    #region Gram-Schmidt (construct orthonormal basis using Gram-Schmidt)
                    for(i = 1; i <= m; i++) {
                        this.NoOfIterations++;

                        #region Arnoldi procdure

                        //w = M \ (A*V(:,i));                         
                        Matrix.SpMV(1.0, V[i - 1], 0.0, z);
                        if(this.Precond != null) {
                            w.Clear();
                            this.Precond.Solve(w, z);
                        } else {
                            w.SetV(z);
                        }

                        for(int k = 1; k <= i; k++) {
                            H[k - 1, i - 1] = GenericBlas.InnerProd(w, V[k - 1]).MPISum(Matrix.MPI_Comm);
                            //w = w - H(k,i)*V(:,k);
                            w.AccV(-H[k - 1, i - 1], V[k - 1]);
                        }

                        double norm_w = w.L2NormPow2().MPISum(Matrix.MPI_Comm).Sqrt();
                        H[i + 1 - 1, i - 1] = norm_w; // the +1-1 actually makes me sure I haven't forgotten to subtract -1 when porting the code
                                                      //V(:,i+1) = w / H(i+1,i);
                        V[i + 1 - 1].SetV(w, alpha: (1.0 / norm_w));
                        #endregion

                        #region Givens rotation

                        for(int k = 1; k <= i - 1; k++) {
                            // apply Givens rotation, H is Hessenberg-Matrix
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
                        error2 = Math.Abs(s[i + 1 - 1]);

                        // For Residual tracking, do not delete
                        /*
                        y = new double[i];
                        H.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { i - 1, i - 1 }).Solve(y, s.GetSubVector(0, i));
                        double[] Xtmp = X.CloneAs();
                        for (int ii = 0; ii < i; ii++) {
                            Xtmp.AccV(y[ii], V[ii]);
                        }
                        double[] ztmp = z.CloneAs();
                        ztmp.SetV(B);
                        Matrix.SpMV(-1.0, Xtmp, 1.0, ztmp);
                        IterationCallback?.Invoke(iter, Xtmp.CloneAs(), ztmp.CloneAs(), this.m_mgop);
                        */

                        var term1 = TerminationCriterion(totIterCounter, iter0_error2, error2);
                        if(!term1.bNotTerminate) {
                            // update approximation and exit
                            //y = H(1:i,1:i) \ s(1:i);    
                            y = new double[i];
                            H.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { i - 1, i - 1 })
                                .Solve(y, s.GetSubVector(0, i));

                            // x = x + V(:,1:i)*y;
                            for(int ii = 0; ii < i; ii++) {
                                X.AccV(y[ii], V[ii]);
                            }
                            this.m_Converged = term1.bSuccess;
                            tr.Info($"convergence reached (2), success? {term1.bSuccess}: iter {totIterCounter}:   iter0_err = {iter0_error2}, err = {error2}");
                            break;
                        }

                        totIterCounter++; // every preconditioner call is count as one iteration
                    }
                    #endregion

                    var term2 = TerminationCriterion(totIterCounter, iter0_error2, error2);
                    if(!term2.bNotTerminate) {
                        this.m_Converged = term2.bSuccess;
                        tr.Info($"convergence reached (3), success? {term2.bSuccess}: iter {totIterCounter}:   iter0_err = {iter0_error2}, err = {error2}");
                        break;
                    }

                    // y = H(1:m,1:m) \ s(1:m);
                    y = new double[m];
                    H.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { m - 1, m - 1 })
                        .Solve(y, s.GetSubVector(0, m));
                    // update approximation: x = x + V(:,1:m)*y;  
                    for(int ii = 0; ii < m; ii++) {
                        X.AccV(y[ii], V[ii]);
                    }

                    // compute residual: r = M \ ( b-A*x )     
                    z.SetV(B);
                    Matrix.SpMV(-1.0, X, 1.0, z);
                    error2 = z.MPI_L2Norm(Matrix.MPI_Comm);
                    IterationCallback?.Invoke(totIterCounter, X.CloneAs(), z.CloneAs(), this.m_mgop as MultigridOperator);


                    if(this.Precond != null) {
                        r.Clear();
                        this.Precond.Solve(r, z);
                    } else {
                        r.SetV(z);
                    }

                    norm_r = r.MPI_L2Norm(Matrix.MPI_Comm);
                    s[i + 1 - 1] = norm_r;
                    error = s[i + 1 - 1] / bnrm2;        // % check convergence
                                                         //  if (error2 <= m_Tolerance) Check for error not error2
                                                         //break;
                }

                if(IterationCallback != null) {
                    z.SetV(B);
                    Matrix.SpMV(-1.0, X, 1.0, z);
                    IterationCallback(totIterCounter, X.CloneAs(), z.CloneAs(), this.m_mgop as MultigridOperator);
                }


                if(!object.ReferenceEquals(_X, X))
                    _X.SetV(X);
                //B.SetV(z);

                // Disposing should not be done here, 
                // otherwise the Precond must be initialized very often.
                //if (this.Precond is IDisposable)
                //    (this.Precond as IDisposable).Dispose();
            }
        }

        /// <summary>
        /// ~
        /// </summary>
        public int IterationsInNested {
            get {
                if(this.Precond != null)
                    return this.Precond.IterationsInNested + this.Precond.ThisLevelIterations;
                else
                    return 0;
            }
        }

        /// <summary>
        /// ~
        /// </summary>
        public int ThisLevelIterations {
            get {
                return this.NoOfIterations;
            }
        }

        /// <summary>
        /// ~
        /// </summary>
        public bool Converged {
            get {
                return m_Converged;
            }
        }

        public void ResetStat() {
            this.m_Converged = false;
            this.NoOfIterations = 0;
            if(this.Precond != null)
                this.Precond.ResetStat();
        }


        public object Clone() {
            SoftGMRES Clone = new SoftGMRES();
            Clone.IterationCallback = this.IterationCallback;
            Clone.MaxKrylovDim = this.MaxKrylovDim;
            Clone.TerminationCriterion = this.TerminationCriterion;
            if(this.Precond != null)
                Clone.Precond = this.Precond.CloneAs();
            return Clone;
        }

        public void Dispose() {
            if(this.Precond != null)
                this.Precond.Dispose();
        }

        public long UsedMemory() {
            long Memory = 0;
            if(Precond != null)
                Memory += Precond.UsedMemory();
            return Memory;
        }
    }
}
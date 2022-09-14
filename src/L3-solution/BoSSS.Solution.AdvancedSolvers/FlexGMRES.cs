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
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Tracing;

namespace BoSSS.Solution.AdvancedSolvers {
    
    /// <summary>
    /// Flexible GMRES (FGMRES), i.e. GMRES with flexible preconditioner,
    /// accordting to
    /// @article{saad_flexible_1993,
	///     title = {A Flexible Inner-Outer Preconditioned {GMRES} Algorithm},
    ///     volume = {14},
    ///     issn = {1064-8275, 1095-7197},
    ///     url = {http://epubs.siam.org/doi/abs/10.1137/0914028},
    ///     doi = {10.1137/0914028},
    ///     number = {2},
    ///     journal = {{SIAM} Journal on Scientific Computing},
    ///     author = {Saad, Youcef},
    ///     year = {1993},
    ///     pages = {461--469}
    /// }
    /// </summary>
    public class FlexGMRES : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination {

        MultigridOperator m_mgop;

        public void Init(MultigridOperator op) {
            using(var tr = new FuncTrace()) {
                if(object.ReferenceEquals(op, this.m_mgop))
                    return; // already initialized
                else
                    this.Dispose();

                var M = op.OperatorMatrix;
                var MgMap = op.DgMapping;
                this.m_mgop = op;

                if(!M.RowPartitioning.EqualsPartition(MgMap))
                    throw new ArgumentException("Row partitioning mismatch.");
                if(!M.ColPartition.EqualsPartition(MgMap))
                    throw new ArgumentException("Column partitioning mismatch.");

                foreach(var pc in PrecondS) {
                    pc.Init(m_mgop);
                }
            }
        }

        /// <summary>
        /// preconditioners used 
        /// </summary>
        public ISolverSmootherTemplate[] PrecondS;

        /// <summary>
        /// ~
        /// </summary>
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get;
            set;
        }

        /// <summary>
        /// ctor
        /// </summary>
        public FlexGMRES() {
            TerminationCriterion = (iIter, r0, ri) => (iIter <= 1, true);
        }


          
        /// <summary>
        /// Number of solution vectors in the internal Krylov-Space
        /// </summary>
        public int MaxKrylovDim = 50;

        /// <summary>
        /// ~
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            int L = B.Count;

            if(X.Count != this.m_mgop.OperatorMatrix.ColPartition.LocalLength)
                throw new ArgumentException("Wrong length of unknowns vector.", "X");
            if(B.Count != this.m_mgop.OperatorMatrix.RowPartitioning.LocalLength)
                throw new ArgumentException("Wrong length of RHS vector.", "B");
            
            MultidimensionalArray H = MultidimensionalArray.Create(MaxKrylovDim + 1, MaxKrylovDim);

            double[] X0 = X.ToArray();
            double[] R0 = new double[B.Count];
            
            List<double[]> Vau = new List<double[]>();
            List<double[]> Zed = new List<double[]>();

            int iPrecond = 0; // counter which iterates through preconditioners

            int iIter = 0;
            double iter0_l2Residual = 0;
            double iter_l2Residual = 0;
            while(true) {

                // init for Arnoldi
                // -----------------

                R0.SetV(B);
                this.m_mgop.OperatorMatrix.SpMV(-1.0, X0, 1.0, R0);
                iter0_l2Residual = R0.MPI_L2Norm();
                if (iIter == 0) {
                    iter_l2Residual = iter0_l2Residual;
                }

                // callback
                this.IterationCallback?.Invoke(iIter, X0.CloneAs(), R0.CloneAs(), this.m_mgop);

                // termination condition
                var term = TerminationCriterion(iIter, iter0_l2Residual, iter_l2Residual);
                if (!term.bNotTerminate) {
                    this.m_Converged = term.bSuccess;
                    break;
                }
                                    
                
                double beta = R0.L2NormPow2().MPISum().Sqrt();
                Vau.Add(R0.CloneAs());
                Vau[0].ScaleV(1.0 / beta);

                // inner iteration (Arnoldi)
                // --------------------------

                for(int j = 0; j < MaxKrylovDim; j++) {
                    Debug.Assert(Vau.Count == Zed.Count + 1);

                    double[] Zj = new double[L];
                    this.PrecondS[iPrecond].Solve(Zj, Vau[j]);
                    iPrecond++;
                    if(iPrecond >= this.PrecondS.Length)
                        iPrecond = 0;
                    iIter++; // we count preconditioner calls as iterations

                    Zed.Add(Zj);

                    double[] W = new double[L];
                    this.m_mgop.OperatorMatrix.SpMV(1.0, Zj, 0.0, W);

                    for(int i = 0; i <= j; i++) {
                        double hij = GenericBlas.InnerProd(W, Vau[i]).MPISum();
                        H[i, j] = hij;
                        W.AccV(-hij, Vau[i]);
                    }

                    double wNorm = W.L2NormPow2().MPISum().Sqrt();
                    H[j + 1, j] = wNorm;
                    W.ScaleV(1.0 / wNorm); // W is now Vau[j + 1] !!

                    Vau.Add(W);

                    Debug.Assert(Vau.Count == Zed.Count + 1);

                    //if (j >= 1) {
                    //    var y = new double[j];
                    //    H.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { j - 1, j - 1 }).Solve(y, W.GetSubVector(0, j));
                    //    double[] Xtmp = X.ToArray();
                    //    for (int ii = 0; ii < j; ii++) {
                    //        Xtmp.AccV(y[ii], Vau[ii]);
                    //    }
                    //    double[] ztmp = new double[Xtmp.Length];
                    //    ztmp.SetV(B);
                    //    this.m_mgop.OperatorMatrix.SpMV(-1.0, Xtmp, 1.0, ztmp);
                    //    IterationCallback?.Invoke(j, Xtmp, ztmp, this.m_mgop);
                    //}
                }

                // compute minimized-Residual solution over 'Zed'-Vectors
                // -------------------------------------------------------
                double[] alpha;
                {
                    alpha = new double[MaxKrylovDim];
                    //var alpha2 = new double[MaxKrylovDim];
                    

                    Debug.Assert(Vau.Count == H.NoOfRows);
                    Debug.Assert(Zed.Count == H.NoOfCols);
                    
                    // solve small minimization problem 
                    double[] minimiRHS = new double[MaxKrylovDim + 1];
                    minimiRHS[0] = beta;
                    H.LeastSquareSolve(alpha, minimiRHS);  // using LAPACK DGELSY
                    //H.CloneAs().LeastSquareSolve(alpha2, minimiRHS.CloneAs());

                    //var Q = H;  // 
                    //var Qt = H.Transpose();
                    //var lhs = Qt.GEMM(Q);
                    //var rhs = new double[Qt.NoOfRows];
                    //Qt.gemv(1.0, minimiRHS, 0.0, rhs);
                    //lhs.Solve(alpha, rhs);

                    alpha.ClearEntries();
                    alpha[0] = beta;
                }

                for(int i = 0; i < MaxKrylovDim; i++)
                    X0.AccV(alpha[i], Zed[i]);
                // now, X0 = 'old X0' + sum(Z[i]*alpha[i]),
                // i.e. X0 is the new X for the next iteration


                // restart
                // -------
                Vau.Clear();
                Zed.Clear();
                H.Clear();

            }
        }

        

        bool m_Converged = false;
        int m_ThisLevelIterations = 0;

        public int IterationsInNested {
            get {
                if(this.PrecondS != null)
                    return this.PrecondS.Sum(pc => pc.IterationsInNested + pc.ThisLevelIterations);
                else
                    return 0;
            }
        }

        public int ThisLevelIterations {
            get { return this.m_ThisLevelIterations; }
        }

        public bool Converged {
            get { return this.m_Converged; }
        }

        public void ResetStat() {
            m_Converged = false;
            m_ThisLevelIterations = 0;
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        public object Clone() {
            var clone = new FlexGMRES();
            clone.TerminationCriterion = this.TerminationCriterion;
            clone.MaxKrylovDim = this.MaxKrylovDim;
            var tmp = new List<ISolverSmootherTemplate>();
            foreach (ISolverSmootherTemplate precond in this.PrecondS) {
                tmp.Add(precond.CloneAs());
            }
            clone.PrecondS = tmp.ToArray();
            return clone;
        }

        public void Dispose() {
            throw new NotImplementedException();
        }

        public long UsedMemory() {
            throw new NotImplementedException();
        }
    }
}

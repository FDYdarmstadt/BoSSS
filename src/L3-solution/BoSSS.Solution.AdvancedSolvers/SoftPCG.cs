﻿/* =======================================================================
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
using ilPSP.Tracing;

namespace BoSSS.Solution.AdvancedSolvers {
    /// <summary>
    /// Pre-conditioned (<see cref="Precond"/>) conjugate gradient algorithm.
    /// </summary>
    public class SoftPCG : ISolverSmootherTemplate, ISolverWithCallback {

        /*
        static public ISolverSmootherTemplate InitMultigridChain(MultigridOperator MgOp,
            Action<int, SoftPCG> ParamsSeter,
            Func<ISolverSmootherTemplate> CoarsestSolverFactory) {
            //
            if(MgOp.CoarserLevel == null) {
                return CoarsestSolverFactory();
            } else {
                var MgTop = new SoftPCG();
                ParamsSeter(MgOp.LevelIndex, MgTop);

                var R = new GenericRestriction();
                MgTop.Precond = R;
                R.CoarserLevelSolver = InitMultigridChain(MgOp.CoarserLevel, ParamsSeter, CoarsestSolverFactory);
                                
                return MgTop;
            }
        }
        */
        
        //public double m_Tolerance = 1.0e-10;
        //public int m_MaxIterations = 10000;
        //public int m_MinIterations = 5;

        /// <summary>
        /// ~
        /// </summary>
        public Func<int, double, double, bool> TerminationCriterion {
            get;
            set;
        }

        /// <summary>
        /// ctor
        /// </summary>
        public SoftPCG() {
            TerminationCriterion = (iIter, r0, ri) => iIter <= 1;
        }


        /// <summary>
        /// ~
        /// </summary>
        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        /// <summary>
        /// ~
        /// </summary>
        public ISolverSmootherTemplate Precond {
            get;
            set;
        }

        public int NoOfIterations = 0;

        /// <summary>
        /// implementation of the CG algorithm
        /// </summary>
        public void Solve<Vec1, Vec2>(Vec1 _x, Vec2 _R)
            where Vec1 : IList<double>
            where Vec2 : IList<double> //
        {
            using (new FuncTrace()) {
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

                //double[] R = rhs; // rhs is only needed once, so we can use it to store residuals
                double[] V = new double[L];
                double[] Z = new double[L];


                // compute P0, R0
                // ==============
                GenericBlas.dswap(L, x, 1, P, 1);
                m_Matrix.SpMV(-1.0, P, 1.0, R);
                IterationCallback?.Invoke(NoOfIterations, P.CloneAs(), R.CloneAs(), this.m_MgOp as MultigridOperator);

                GenericBlas.dswap(L, x, 1, P, 1);
                if (Precond != null) {
                    Precond.Solve(Z, R);
                    P.SetV(Z);
                } else {
                    P.SetV(R);
                }

                double alpha = R.InnerProd(P).MPISum();
                double alpha_0 = alpha;
                double ResNorm;

                ResNorm = Math.Sqrt(alpha);
                double ResNorm0 = ResNorm;
                if (TerminationCriterion(0, ResNorm0, ResNorm)) {
                    this.m_Converged = true;
                    return;
                }

                // iterate
                // =======
                NoOfIterations++; // one iteration has already been performed (P0, R0)
                for (int n = 1; true; n++) {

                    if (TerminationCriterion(n,ResNorm0, ResNorm)) {
                        this.m_Converged = true;
                        break;
                    }
                    NoOfIterations++;

                    if (Math.Abs(alpha) <= double.Epsilon) {
                        // numerical breakdown
                        break;
                    }


                    m_Matrix.SpMV(1.0, P, 0, V);
                    double VxP = V.InnerProd(P).MPISum();
                    //Console.WriteLine("VxP: {0}", VxP);
                    if (double.IsNaN(VxP) || double.IsInfinity(VxP))
                        throw new ArithmeticException();
                    double lambda = alpha / VxP;
                    if (double.IsNaN(lambda) || double.IsInfinity(lambda))
                        throw new ArithmeticException();


                    x.AccV(lambda, P);

                    R.AccV(-lambda, V);

                    if (IterationCallback != null) {
                        IterationCallback(NoOfIterations, x.CloneAs(), R.CloneAs(), this.m_MgOp as MultigridOperator);
                    }

                    if (Precond != null) {
                        Z.Clear();
                        Precond.Solve(Z, R);
                    } else {
                        Z.SetV(R);
                    }

                    double alpha_neu = R.InnerProd(Z).MPISum();

                    // compute residual norm
                    ResNorm = R.L2NormPow2().MPISum().Sqrt();

                    P.ScaleV(alpha_neu / alpha);
                    P.AccV(1.0, Z);

                    alpha = alpha_neu;
                }

                if (!object.ReferenceEquals(_x, x))
                    _x.SetV(x);
                if (!object.ReferenceEquals(_R, R))
                    _R.SetV(R);


                return;
            }
        }

        IOperatorMappingPair m_MgOp;
        BlockMsrMatrix m_Matrix => m_MgOp.OperatorMatrix;


        /// <summary>
        /// ~
        /// </summary>
        public void Init(MultigridOperator op) {
            using(new FuncTrace()) {
                using(new FuncTrace()) {
                    if(object.ReferenceEquals(op, this.m_MgOp))
                        return; // already initialized
                    else
                        this.Dispose(); // must re-initialize

                    var M = op.OperatorMatrix;
                    var MgMap = op.DgMapping;
                    this.m_MgOp = op;

                    if(!M.RowPartitioning.EqualsPartition(MgMap))
                        throw new ArgumentException("Row partitioning mismatch.");
                    if(!M.ColPartition.EqualsPartition(MgMap))
                        throw new ArgumentException("Column partitioning mismatch.");

                    
                    if(Precond != null)
                        Precond.Init(op);
                }
            }
        }
        public int IterationsInNested {
            get {
                if(this.Precond != null)
                    return this.Precond.IterationsInNested + this.Precond.ThisLevelIterations;
                else
                    return 0;
            }
        }

        public int ThisLevelIterations {
            get {
                return this.NoOfIterations;
            }
        }

        public bool Converged {
            get {
                return m_Converged;
            }
        }

        bool m_Converged = false;

        public void ResetStat() {
            this.m_Converged = false;
            this.NoOfIterations = 0;
            this?.Precond.ResetStat();
        }

        /// <summary>
        /// ToDo: Cloning of Preconditioner is not supported.
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            var clone = new SoftPCG();
            clone.IterationCallback = this.IterationCallback;
            clone.NoOfIterations = this.NoOfIterations;
            clone.TerminationCriterion = this.TerminationCriterion;
            clone.Precond = null;
            return clone;
        }

        public void Dispose() {
            Precond?.Dispose();
            m_MgOp = null;  // setting this to null ensures that Init(...) will actually initialize the solvers
        }

        public long UsedMemory() {
            throw new NotImplementedException();
        }
    }
}

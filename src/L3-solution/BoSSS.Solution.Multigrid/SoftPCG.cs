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

namespace BoSSS.Solution.Multigrid {
    public class SoftPCG : ISolverSmootherTemplate, ISolverWithCallback {

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


        BlockMsrMatrix m_Matrix;
        ilPSP.LinSolvers.monkey.MatrixBase m_OptMatrix;

        public double m_Tolerance = 1.0e-10;
        public int m_MaxIterations = 10000;
        public int m_MinIterations = 5;

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        public ISolverSmootherTemplate Precond {
            get;
            set;
        }

        int NoOfIterations = 0;

        /// <summary>
        /// implementation of the CG algorithm
        /// </summary>
        /// <param name="_x"></param>
        /// <param name="_R"></param>
        /// <param name="stats"></param>
        public void Solve<Vec1, Vec2>(Vec1 _x, Vec2 _R)
            where Vec1 : IList<double>
            where Vec2 : IList<double> {
            double[] x, R;
            if(_x is double[]) {
                x = _x as double[];
            } else {
                x = _x.ToArray();
            }
            if(_R is double[]) {
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
            m_OptMatrix.SpMV(-1.0, P, 1.0, R);
            if(IterationCallback != null) {
                IterationCallback(NoOfIterations, P.CloneAs(), R.CloneAs(), this.m_MgOp);
            }
            GenericBlas.dswap(L, x, 1, P, 1);
            if(Precond != null) {
                Precond.Solve(Z, R);
                P.SetV(Z);
            } else {
                P.SetV(R);
            }

            double alpha = R.InnerProd(P).MPISum();
            double alpha_0 = alpha;
            double ResNorm;

            ResNorm = Math.Sqrt(alpha);

            // iterate
            // =======
            NoOfIterations++; // one iteration has allready been performed (P0, R0)
            for(int n = 1; n < this.m_MaxIterations; n++) {

                if(ResNorm <= m_Tolerance && NoOfIterations >= m_MinIterations) {
                    this.m_Converged = true;
                    break;
                }
                NoOfIterations++;

                if(Math.Abs(alpha) <= double.Epsilon) {
                    // numerical breakdown
                    break;
                }

                
                m_Matrix.SpMV(1.0, P, 0, V);
                double VxP = V.InnerProd(P).MPISum();
                //Console.WriteLine("VxP: {0}", VxP);
                if(double.IsNaN(VxP) || double.IsInfinity(VxP))
                    throw new ArithmeticException();
                double lambda = alpha / VxP;
                if(double.IsNaN(lambda) || double.IsInfinity(lambda))
                    throw new ArithmeticException();


                x.AccV(lambda, P);

                R.AccV(-lambda, V);

                if(IterationCallback != null) {
                    IterationCallback(NoOfIterations, x.CloneAs(), R.CloneAs(), this.m_MgOp);
                }

                if(Precond != null) {
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

            if(!object.ReferenceEquals(_x, x))
                _x.SetV(x);
            if(!object.ReferenceEquals(_R, R))
                _R.SetV(R);


            return;
        }

        MultigridOperator m_MgOp;

        public void Init(MultigridOperator op) {
            var M = op.OperatorMatrix;
            var MgMap = op.Mapping;
            this.m_MgOp = op;
                        
            if(!M.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if(!M.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");
            this.m_Matrix = M;
            m_OptMatrix = new ilPSP.LinSolvers.monkey.CPU.RefMatrix(M.ToMsrMatrix());
            if(Precond != null)
                Precond.Init(op);
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
            if(this.Precond != null)
                this.Precond.ResetStat();
        }
    }


}

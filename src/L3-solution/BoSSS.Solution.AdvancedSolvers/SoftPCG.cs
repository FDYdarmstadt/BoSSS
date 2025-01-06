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
using ilPSP.Tracing;
using System.Runtime.Serialization;
using ilPSP.LinSolvers.PARDISO;

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

                var Res = new GenericRestriction();
                MgTop.Precond = Res;
                Res.CoarserLevelSolver = InitMultigridChain(MgOp.CoarserLevel, ParamsSeter, CoarsestSolverFactory);
                                
                return MgTop;
            }
        }
        */

		//public double m_Tolerance = 1.0e-10;
		[DataMember]
		public int MaxIterations = 10000;
		//public int m_MinIterations = 5;

		[DataMember]
		public double ConvergenceCriterion = 1e-10;

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
        public SoftPCG(bool calculateWithMatrix = true) {
            this.CalculateWithMatrix = calculateWithMatrix;
            TerminationCriterion = (int iter, double R0_l2, double R_l2) => R_l2 < R0_l2 * ConvergenceCriterion + ConvergenceCriterion;
        }


		/// <summary>
		/// ~
		/// </summary>
		public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

		public Action<double[], double[]> InnerIterBefore;
		public Action<double[], double[]> InnerIterAfter;


		/// <summary>
		/// ~
		/// </summary>
		public PARDISOSolver Precond {
            get;
            set;
        }

		/// <summary>
		/// When the matrix is not explicilty available, an inner iteration can be defined to calculate Matrix * Search Direction (m_matrix x P)
		/// </summary>
		public ISparseSolverExt InnerCycle {
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
                double[] x, Res;
                if (_x is double[]) {
                    x = _x as double[];
                } else {
                    x = _x.ToArray();
                }
                if (_R is double[]) {
                    Res = _R as double[];
                } else {
                    Res = _R.ToArray();
                }

                int L = x.Length;

                // Search direction
                double[] P = new double[L];

                //double[] Res = rhs; // rhs is only needed once, so we can use it to store residuals
                double[] V = new double[L];
                double[] PrecondRes = new double[L]; //Preconditioned Res


                // compute P0, R0
                // ==============
                GenericBlas.dswap(L, x, 1, P, 1);
                
                //Initial res
                if (CalculateWithMatrix)
                    m_Matrix.SpMV(-1.0, P, 1.0, Res);

                IterationCallback?.Invoke(NoOfIterations, P.CloneAs(), Res.CloneAs(), this.m_MgOp as MultigridOperator);

                GenericBlas.dswap(L, x, 1, P, 1);

                ApplyPreconditioner(PrecondRes, Res);
                P.SetV(PrecondRes);

                double alpha = Res.InnerProd(P).MPISum();
                double alpha_0 = alpha;
                double ResNorm;

                ResNorm = Math.Sqrt(alpha);
                double ResNorm0 = ResNorm;

                bool ShouldContinue = CheckIteration(0, ResNorm0, ResNorm,alpha); // one iteration has already been performed (P0, R0)

				// iterate
				// =======
				for (int n = 1; ShouldContinue; n++) {
                    double VxP = CalculateVxP(V, P);
					double lambda = alpha / VxP;
                    CheckIfAcceptable(lambda);
					x.AccV(lambda, P);
                    Res.AccV(-lambda, V);

					IterationCallback?.Invoke(NoOfIterations, P.CloneAs(), Res.CloneAs(), this.m_MgOp as MultigridOperator);

					ApplyPreconditioner(PrecondRes, Res);
                    double alpha_neu = Res.InnerProd(PrecondRes).MPISum();

                    // compute residual norm
                    ResNorm = Res.L2NormPow2().MPISum().Sqrt();

                    P.ScaleV(alpha_neu / alpha);
                    P.AccV(1.0, PrecondRes);

                    alpha = alpha_neu;
					ShouldContinue = CheckIteration(n, ResNorm0, ResNorm, alpha);
				}
                Console.WriteLine("ResNorm " + ResNorm);
				if (!object.ReferenceEquals(_x, x))
                    _x.SetV(x);
                if (!object.ReferenceEquals(_R, Res))
                    _R.SetV(Res);


                return;
            }
        }

		double CalculateVxP(double[] V, double[] P) {
			if (m_Matrix != null) {
				m_Matrix.SpMV(1.0, P, 0, V);
			} else {
                double[] Sol = new double[InnerCycle.GetMatrix().RowPartitioning.LocalLength];
                double[] Y = new double[InnerCycle.GetMatrix().RowPartitioning.LocalLength];
				InnerIterBefore(P,Y);
                InnerCycle.Solve(Sol, Y);
				InnerIterAfter(Sol,V);
			}

			double VxP = V.InnerProd(P).MPISum();
			CheckIfAcceptable(VxP);
			return VxP;
		}

        bool CheckIteration(int n, double ResNorm0, double ResNorm, double alpha) {
			NoOfIterations++;
			bool ret = true;

			if (TerminationCriterion(n, ResNorm0, ResNorm)) {
				this.m_Converged = true;
                ret = false;
			}

            if (n >= MaxIterations) { 
                ret = false;
			}

			if (Math.Abs(alpha) <= double.Epsilon) {
				// numerical breakdown
				ret = false;
			}
            return ret;
		}

        void CheckIfAcceptable(double X) {
			if (double.IsNaN(X) || double.IsInfinity(X))
				throw new ArithmeticException();
		}

        void ApplyPreconditioner(double[] Z, double[] Res) {
			if (Precond != null) {
				Z.Clear();
				Precond.Solve(Z, Res);
			} else {
				Z.SetV(Res);
			}
		}

		IOperatorMappingPair m_MgOp;

        bool CalculateWithMatrix;

		BlockMsrMatrix m_Matrix {
            get {
				return CalculateWithMatrix ? m_MgOp.OperatorMatrix : null;
			}
        }


		//BlockMsrMatrix m_Matrix => m_MgOp.OperatorMatrix;


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

                    
                    //if(Precond != null)
                    //    Precond.Init(op);
                }
            }
        }
        public int IterationsInNested {
            get {
                if (this.Precond != null)
                    return 0;// this.Precond.IterationsInNested + this.Precond.ThisLevelIterations;
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
            //this?.Precond.ResetStat();
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

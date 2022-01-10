using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers
{
    public class kcycle : ISolverSmootherTemplate, ISolverWithCallback
    {
        private MultigridOperator m_MgOperator;
        public ISolverSmootherTemplate CoarserLevelSolver;
        public ISolverSmootherTemplate PreSmoother;
        public ISolverSmootherTemplate PostSmoother;
        public ISolverSmootherTemplate ThisLevelKrylovMethod;
        BlockMsrMatrix OpMatrix;

        public void Init(MultigridOperator op) {
            using (var tr = new FuncTrace()) {
                this.m_MgOperator = op;
                var Mtx = op.OperatorMatrix;
                var MgMap = op.Mapping;

                if (!Mtx.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!Mtx.ColPartition.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Column partitioning mismatch.");


                // set operator
                // ============
                this.OpMatrix = Mtx;


                // initiate coarser level
                // ======================
                if (this.CoarserLevelSolver == null) {
                    throw new NotSupportedException("Missing coarse level solver.");
                    //Console.WriteLine("OrthonormalizationMultigrid: running without coarse solver.");
                } else {
                    if (op.CoarserLevel != null) {
                        this.CoarserLevelSolver.Init(op.CoarserLevel);
                    }
                }

                // initiate krylov method
                // ======================
                ThisLevelKrylovMethod = new FlexGMRES() {
                    PrecondS = new ISolverSmootherTemplate[] { this.CoarserLevelSolver },
                    MaxKrylovDim = 1,
                    TerminationCriterion = (int iter, double r0, double r) => iter <= 1,
                };
                ThisLevelKrylovMethod.Init(op.CoarserLevel);

                // init smoother
                // =============
                if (PreSmoother != null)
                    PreSmoother.Init(op);
                if (PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                    PostSmoother.Init(op);
            }
        }

        int MGDepth {
            get { return m_MGDepth; }
            set { m_MGDepth = Math.Max(value, m_MGDepth); }
        }
        int m_MGDepth = 0;

        /// <summary>
        /// computes the residual on this level.
        /// </summary>
        public double Residual<V1,V2,V3>(V1 rl, V2 xl, V3 bl) 
            where V1 : IList<double>
            where V2 : IList<double>
            where V3 : IList<double>
        {
            rl.SetV( bl);
            OpMatrix.SpMV(-1.0, xl, 1.0, rl);
            
            return rl.MPI_L2Norm();
        }

        /// <summary>
        /// the multigrid iterations for a linear problem
        /// </summary>
        /// <param name="xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="bl">the right-hand-side of the problem</param>
        public void Solve<U, V>(U xl, V bl)
            where U : IList<double>
            where V : IList<double> {
            using (new FuncTrace()) {

                int N = xl.Count;
                int NN = m_MgOperator.CoarserLevel.Mapping.LocalLength; // RestrictionOperator.RowPartitioning.LocalLength;
                double[] rl = new double[N];
                double[] rlp1 = new double[NN];
                var x1 = new double[N];
                var x2 = new double[N];
                var x3 = new double[N];

                double iter0ResidualNorm = Residual(rl, xl, bl);
                double iterNorm = iter0ResidualNorm;

                for (int iIter = 0; true; iIter++) {
                    if (!TerminationCriterion(iIter, iter0ResidualNorm, iterNorm))
                        return;

                    if (PreSmoother != null) {
                        PreSmoother.Solve(x1, rl); // Vorglättung
                        iterNorm = Residual(rl, x1, bl); // Residual on this level
                    }

                    this.m_MgOperator.CoarserLevel.Restrict(rl, rlp1);

                    // Berechnung der Grobgitterkorrektur
                    var vlp1 = new double[NN];

                    if (m_MGDepth-2<=m_MgOperator.LevelIndex)
                        this.CoarserLevelSolver.Solve(vlp1, rlp1);
                    else
                        ThisLevelKrylovMethod.Solve(vlp1, rlp1);
                    
                    // Prolongation der Grobgitterkorrektur
                    this.m_MgOperator.CoarserLevel.Prolongate(1.0, x2, 1.0, vlp1);

                    // check residual after coarse grid correction
                    Residual(rl, x2, bl); // Residual on this level

                    // Nachglättung
                    if (PostSmoother != null) {
                        PostSmoother.Solve(x3, rl);
                        iterNorm = Residual(rl, x3, bl);
                    }

                    xl.SumV(0.0, 1.0, x1, 1.0, x2, 1.0, x3);

                    // update residual
                    this.IterationCallback?.Invoke(iIter + 1, xl.ToArray(), rl.CloneAs(), this.m_MgOperator);
                    this.ThisLevelIterations++;
                }
            }
        }

        public bool Converged {
            get;
            private set;
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }
        public Func<int, double, double, bool> TerminationCriterion {
            get;
            set;
        }

        public int IterationsInNested {
            get {
                int iter = 0;

                if (PreSmoother != null)
                    iter += this.PreSmoother.IterationsInNested + this.PreSmoother.ThisLevelIterations;

                if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                    iter += this.PostSmoother.IterationsInNested + this.PostSmoother.ThisLevelIterations;

                iter += this.CoarserLevelSolver.IterationsInNested + this.CoarserLevelSolver.ThisLevelIterations;

                return iter;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public int ThisLevelIterations {
            get;
            private set;
        }

        public void ResetStat() {
            this.Converged = false;
            this.ThisLevelIterations = 0;
            if (this.PreSmoother != null)
                this.PreSmoother.ResetStat();
            if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                this.PostSmoother.ResetStat();
            if (this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.ResetStat();
        }
        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

        public void Dispose() {
            throw new NotImplementedException();
        }

        public long UsedMemory() {
            throw new NotImplementedException();
        }
    }
}

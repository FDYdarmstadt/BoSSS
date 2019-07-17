using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Executes two or more solvers in sequence;
    /// This could be e.g. some predictor/corrector pair, or some coarse solver/fine smoother pair.
    /// </summary>
    public class SolverSquence : ISolverSmootherTemplate, ISolverWithCallback {
        
        MultigridOperator m_MgOperator;

        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(MultigridOperator op) {
            
            // init & check
            // ============

            this.m_MgOperator = op;
            //OpMatrix = op.OperatorMatrix;
            OpMatrix = new ilPSP.LinSolvers.monkey.CPU.RefMatrix(op.OperatorMatrix.ToMsrMatrix());
            var MgMap = op.Mapping;

            if(!OpMatrix.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if(!OpMatrix.ColPartition.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Column partitioning mismatch.");

            // init solver chain
            // =================
            if (SolverChain == null || SolverChain.Length <= 0)
                throw new NotSupportedException("Illegal configuration.");
            for(int i = 0; i < SolverChain.Length; i++) {
                SolverChain[i].Init(op);
            }
        }

        ilPSP.LinSolvers.monkey.CPU.RefMatrix OpMatrix;
        //BlockMsrMatrix OpMatrix;

        /// <summary>
        /// A list of solvers which is used sequentially, after each other.
        /// </summary>
        public ISolverSmootherTemplate[] SolverChain;

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }


        

        /// <summary>
        /// computes the residual on this level.
        /// </summary>
        public void Residual<V1, V2, V3>(V1 rl, V2 xl, V3 bl)
            where V1 : IList<double>
            where V2 : IList<double>
            where V3 : IList<double> {
            OpMatrix.SpMV(-1.0, xl, 0.0, rl);
            rl.AccV(1.0, bl);
        }


        int m_Gamma = 1;


        public int Gamma {
            get {
                return m_Gamma;
            }
            set {
                if (value < 1)
                    throw new ArgumentException();
                m_Gamma = value;
            }
        }

        /// <summary>
        /// How often we cycle over the <see cref="SolverChain"/>.
        /// </summary>
        public int Fixedterations = 1;


        //bool m_converged = false;
        int NoOfIterations = 0;


        /// <summary>
        /// the multigrid iterations for a linear problem
        /// </summary>
        /// <param name="xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="bl">the right-hand-side of the problem</param>
        public void Solve<U, V>(U xl, V bl)
            where U : IList<double>
            where V : IList<double> {
            //

            int N = xl.Count;
            double[] ri = new double[N];
            double[] C = new double[N];

            if (this.IterationCallback != null) {
                var _bl = bl.ToArray();
                this.OpMatrix.SpMV(-1.0, xl, 1.0, _bl);
                this.IterationCallback(0, xl.ToArray(), _bl, this.m_MgOperator);
            }

            if (Fixedterations <= 0)
                throw new NotSupportedException("illegal configuration");

            for (int iIter = 0; iIter < Fixedterations; iIter++) {

                //var DGBasis = m_MgOperator.BaseGridProblemMapping.BasisS[0];
                //XDGField ResB4Jacobi = new XDGField((XDGBasis)DGBasis, "Resi_B4Jacobi");
                //m_MgOperator.TransformRhsFrom(ResB4Jacobi.CoordinatesAsVector, bl);

                //if (PreSmoother != null)
                //    PreSmoother.Solve(xl, bl); // Vorglättung
                //Residual(rl, xl, bl); // Residual on this level

                
                //this.NoOfIterations++;

                //this.m_MgOperator.CoarserLevel.Restrict(rl, rlp1);

                //// Berechnung der Grobgitterkorrektur
                //double[] vlp1 = new double[NN];

                //for (int j = 0; j < m_Gamma; j++) {
                //    this.CoarserLevelSolver.Solve(vlp1, rlp1);
                //}

                //// Prolongation der Grobgitterkorrektur
                //this.m_MgOperator.CoarserLevel.Prolongate(1.0, xl, 1.0, vlp1);


                //// Nachglättung
                //if (PostSmoother != null)
                //    PostSmoother.Solve(xl, bl);

                for(int i = 0; i < this.SolverChain.Length; i++) {

                    // compute residual 
                    Residual(ri, xl, bl);

                    // compute correction
                    C.Clear();
                    this.SolverChain[i].Solve(C, ri);

                    // update x
                    xl.AccV(1.0, C);
                }



                if (this.IterationCallback != null) {
                    var _bl = bl.ToArray();
                    this.OpMatrix.SpMV(-1.0, xl, 1.0, _bl);
                    this.IterationCallback(iIter + 1, xl.ToArray(), _bl, this.m_MgOperator);
                }
            }

            //m_converged = false;
        }

        /// <summary>
        /// Sum of iterations in entire <see cref="SolverChain"/>.
        /// </summary>
        public int IterationsInNested {
            get {
                return this.SolverChain.Sum(solver => solver.IterationsInNested + solver.ThisLevelIterations);
                    //((this.PreSmoother != null) ? (this.PreSmoother.IterationsInNested + this.PreSmoother.ThisLevelIterations) : 0)
                    //+ ((this.PostSmoother != null) ? (this.PostSmoother.IterationsInNested + this.PostSmoother.ThisLevelIterations) : 0)
                    //+ (this.CoarserLevelSolver.IterationsInNested + this.CoarserLevelSolver.ThisLevelIterations);
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public int ThisLevelIterations {
            get {
                return this.NoOfIterations;
            }
        }

        /// <summary>
        /// Only true, if all solvers in <see cref="SolverChain"/> are converged.
        /// </summary>
        public bool Converged {
            get {
                return this.SolverChain.All(solver => solver.Converged);
            }
        }

        public void ResetStat() {
            //this.m_converged = false;
            this.NoOfIterations = 0;
            foreach (var s in this.SolverChain)
                s.ResetStat();
        }

        public ISolverSmootherTemplate Clone() {
            var clone = new SolverSquence();
            List<ISolverSmootherTemplate> clonelist =new List<ISolverSmootherTemplate>();
            foreach (ISolverSmootherTemplate solver in this.SolverChain)
                clonelist.Add(solver.Clone());
            clone.SolverChain=clonelist.ToArray();
            return clone;
        }

    }

}

using BoSSS.Foundation;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers {

    /*

    /// <summary>
    /// A multi-grid method, based on orthonormalization. Based on the current solution,
    /// it always tries to use the coarsest (thus cheapest) solver possible,
    /// which can give a significant reduction of the residual
    /// </summary>
    public class DynamicMultigrid  : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination {

        MultigridOperator m_mgop;

        /// <summary>
        /// Control over the lowest multi-grid level: if, at a certain grid level, the number of DOF's is 
        /// lower than this threshold,
        /// a direct solver is used
        /// </summary>
        public int LowestLevelDOFthreshold = 100;

        /// <summary>
        /// Multigrid Operators for each grid level
        /// </summary>
        MultigridOperator[] Op4Level;

        /// <summary>
        /// Preconditioner's for each level
        /// </summary>
        ISolverSmootherTemplate[] PrecondS;

        /// <summary>
        /// Operator matrix for each level
        /// </summary>
        BlockMsrMatrix[] Mtx4Level; 

        /// <summary>
        /// Number of actually used multi-grid levels.
        /// </summary>
        public int NoOfUsedLevels {
            get {
                return Op4Level.Length;
            }
        }

        /// <summary>
        /// Initialization of each multi-grid level
        /// </summary>
        public void Init(MultigridOperator op) {
            using (new FuncTrace()) {
                // checks
                // ======

                var M = op.OperatorMatrix;
                var MgMap = op.Mapping;
                this.m_mgop = op;
                int MpiSize = M._RowPartitioning.MpiSize;

                if (!M.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!M.ColPartition.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Column partitioning mismatch.");

                // determine which levels to use
                // =============================
                {
                    var tempOp4Level = new List<MultigridOperator>();
                    var Op = op;
                    do {
                        tempOp4Level.Add(Op);
                        Op = Op.CoarserLevel;
                        if (Op == null)
                            break;
                    } while (tempOp4Level.Last().Mapping.TotalLength > LowestLevelDOFthreshold);
                    Op4Level = tempOp4Level.ToArray();

                    Mtx4Level = new BlockMsrMatrix[Op4Level.Length];
                    for(int i = 0; i < Mtx4Level.Length; i++) {
                        Mtx4Level[i] = Op4Level[i].OperatorMatrix;
                    }
                }

                // define preconditioner's
                // =======================
                {

                    PrecondS = new ISolverSmootherTemplate[Op4Level.Length];

                    // all levels except the coarsest one 
                    for (int i = 0; i < PrecondS.Length - 1; i++) {
                        int NoOfBlocks = (int)Math.Ceiling((double)(Op4Level[i].Mapping.LocalLength) / (double)LowestLevelDOFthreshold);

                        // we want at least two blocks - otherwise we could use the direct solver directly.
                        if (MpiSize > 1) {
                            // more than one MPI core -> at least one block per core -> globally at least 2 cores
                            NoOfBlocks = Math.Max(1, NoOfBlocks);
                        } else {
                            // singe MPI core -> at least two blocks
                            NoOfBlocks = Math.Max(2, NoOfBlocks);
                        }

                        PrecondS[i] = new Schwarz() {
                            FixedNoOfIterations = 1,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsOnCurrentProcess = NoOfBlocks
                            },
                            Overlap = 1
                        };
                    }

                    // coarsest level
                    PrecondS[PrecondS.Length - 1] = new DirectSolver() {
                        WhichSolver = DirectSolver._whichSolver.PARDISO,
                        TestSolution = false
                    };

                     

                    // init each level
                    for (int i = 0; i < PrecondS.Length; i++) {
                        PrecondS[i].Init(Op4Level[i]);
                    }
                }
            }
        }


        int FindLevel(int L) {
            for(int iLv = 0; iLv < NoOfUsedLevels; iLv++) {
                if(L == Op4Level[iLv].Mapping.LocalLength) {
                    return iLv;
                }
            }
            return -1;
        }


        SinglePhaseField ProlongateToDg(double[] V, string name) {
            double[] Curr = ProlongateToTop(V);

            var gdat = m_mgop.BaseGridProblemMapping.GridDat;
            var basis = m_mgop.BaseGridProblemMapping.BasisS[0];
            var dgCurrentSol = new SinglePhaseField(basis, name);
            Op4Level[0].TransformSolFrom(dgCurrentSol.CoordinateVector, Curr);
            return dgCurrentSol;
        }

        private double[] ProlongateToTop(double[] V) {
            int iLv = FindLevel(V.Length);


            double[] Curr = V;
            for (; iLv > 0; iLv--) {
                double[] Next = new double[Op4Level[iLv - 1].Mapping.LocalLength];
                Op4Level[iLv].Prolongate(1.0, Next, 0.0, Curr);
                Curr = Next;
            }

            return Curr;
        }

        int counter = 0;
        void PlotVectors(IEnumerable<double[]> VV, string[] names) {

            DGField[] all = new DGField[names.Length];
            for(int i = 0; i < VV.Count(); i++) {
                all[i] = ProlongateToDg(VV.ElementAt(i), names[i]);
            }

            Tecplot.Tecplot.PlotFields(all, "OrthoScheme-" + counter, counter, 2);

        }


        /// <summary>
        /// ~
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (new FuncTrace()) {
                // init 
                // ====

                int[] L = this.Op4Level.Select(op => op.Mapping.LocalLength).ToArray();
                int NoLevels = L.Length;
                if (X.Count != L[0])
                    throw new ArgumentException();
                if (B.Count != L[0])
                    throw new ArgumentException();


                double[][] CurrSol = new double[NoLevels][];
                double[][] CurrRestRes = new double[NoLevels][];  // restricted residual
               
                double[][] Rhs = new double[NoLevels][];
                for (int iLv = 0; iLv < NoLevels; iLv++) {
                    CurrSol[iLv] = new double[L[iLv]];
                    CurrRestRes[iLv] = new double[L[iLv]];
                    Rhs[iLv] = new double[L[iLv]];

                    Console.WriteLine("    Level " + iLv);

                    if (iLv == 0) {
                        CurrSol[iLv].SetV(X);
                        Rhs[iLv].SetV(B);

                        CurrRestRes[iLv].SetV(Rhs[iLv]);
                        Mtx4Level[iLv].SpMV(-1.0, CurrSol[iLv], 1.0, CurrRestRes[iLv]);


                    } else {
                        Op4Level[iLv].Restrict(CurrSol[iLv - 1], CurrSol[iLv]);
                        Op4Level[iLv].Restrict(Rhs[iLv - 1], Rhs[iLv]);
                        Op4Level[iLv].Restrict(CurrRestRes[iLv - 1], CurrRestRes[iLv]);

                        
                        //Console.WriteLine("       dist b level " + iLv + " residual and rest. residual " + ResDist);
                    }
                    //CurrRes[iLv] = LevelRes;

                    //double l2_RestRes_iLv = CurrRestRes[iLv].L2NormPow2().MPISum().Sqrt();
                    //double l2_Res_iLv = LevelRes.L2NormPow2().MPISum().Sqrt();
                    //Console.WriteLine("          " + l2_RestRes_iLv + "    " + CurrSol[iLv].L2Norm());
                }

                double[][] LevelRestRes = new double[NoLevels][];  // residual on level
                for (int iLv = 0; iLv < NoLevels - 1; iLv++) {
                    LevelRestRes[iLv] = CurrRestRes[iLv].CloneAs();
                    Op4Level[iLv + 1].Prolongate(-1.0, LevelRestRes[iLv], 1.0, CurrRestRes[iLv + 1]);
                }
                LevelRestRes[NoLevels - 1] = CurrRestRes[NoLevels - 1].CloneAs();

                for (int iLv = 0; iLv < NoLevels; iLv++) {
                    Console.WriteLine("       level " + iLv + ":     " + CurrRestRes[iLv].L2Norm() + "       " + LevelRestRes[iLv].L2Norm());
                }

                

                PlotVectors(new[] { CurrRestRes[0], CurrRestRes[1], LevelRestRes[0] }, new[] { "Res0", "Res1", "Res0f" });


                

                List<double[]> SolHistory = new List<double[]>();
                List<double[]> MxxHistory = new List<double[]>();


               

                throw new NotImplementedException("Work in progress - todo");
            }
        }

       
        /// <summary>
        /// Number of solution vectors in the internal Krylov-Space
        /// </summary>
        public int MaxKrylovDim = 80;
       

        /// <summary>
        /// .
        /// </summary>
        public int IterationsInNested {
            get {
                if (this.PrecondS != null)
                    return this.PrecondS.Sum(pc => pc.IterationsInNested + pc.ThisLevelIterations);
                else
                    return 0;
            }
        }

        /// <summary>
        /// ~
        /// </summary>
        public int ThisLevelIterations {
            get;
            private set;
        }

        /// <summary>
        /// ~
        /// </summary>
        public bool Converged {
            get;
            private set;
        }

        /// <summary>
        /// ~
        /// </summary>
        public void ResetStat() {
            Converged = false;
            ThisLevelIterations = 0;
        }

        /// <summary>
        /// %
        /// </summary>
        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        /// <summary>
        /// ~
        /// </summary>
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get;
            set;
        }

        /// <summary>
        /// ~
        /// </summary>
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

    */
}

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

    /// <summary>
    /// A multi-grid method, based on orthonormalization. Based on the current solution,
    /// it always tries to use the coarsest (thus cheapest) solver possible,
    /// which can give a significant reduction of the residual
    /// </summary>
    public class DynamicMultigrid  : ISolverSmootherTemplate, ISolverWithCallback {

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
                            m_MaxIterations = 1,
                            CoarseSolver = null,
                            m_BlockingStrategy = new Schwarz.METISBlockingStrategy() {
                                NoOfPartsPerProcess = NoOfBlocks
                            },
                            Overlap = 1
                        };
                    }

                    // coarsest level
                    PrecondS[PrecondS.Length - 1] = new SparseSolver() {
                        WhichSolver = SparseSolver._whichSolver.PARDISO,
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
                //double[][] CurrRes = new double[NoLevels][];  // residual on level
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

                /*
                int CsLv = NoLevels - 1; // index of coarsest level
                double[] Corr_CsLv = new double[L[CsLv]];
                PrecondS[CsLv].Solve(Corr_CsLv, CurrRestRes[CsLv]);

                //double[] Corr2_CsLv = new double[L[CsLv]];
                //PrecondS[CsLv].Solve(Corr2_CsLv, CurrRes[CsLv]);

                double[] Corr_0Lv = Prolongate(Corr_CsLv);
                //double[] Corr2_0Lv = Prolongate(Corr2_CsLv);
                Corr_0Lv.AccV(1.0, CurrSol[0]);
                //Corr2_0Lv.AccV(1.0, CurrSol[0]);
                */


                List<double[]> SolHistory = new List<double[]>();
                List<double[]> MxxHistory = new List<double[]>();



                for (int iIter = 0; iIter < MaxIter; iIter++) {


                }


                /*
                var Mtx = this.m_mgop.OperatorMatrix;


                // residual of initial guess
                // =========================

                // history of solutions and residuals (max vector length 'MaxKrylovDim')
                List<double[]> SolHistory = new List<double[]>();
                List<double[]> MxxHistory = new List<double[]>();

                double[] Correction = new double[L];
                double[] Mxx = new double[L];
                double[] CurrentSol = new double[L];
                double[] CurrentRes = new double[L];

                CurrentSol.SetV(X, 1.0);
                CurrentRes.SetV(B, 1.0);
                Mtx.SpMV(-1.0, CurrentSol, 1.0, CurrentRes);
                int KrylovDim = 0;

                double[] Residual0 = CurrentRes.CloneAs();
                double[] Solution0 = CurrentSol.CloneAs();

                List<double> _R = new List<double>();

                // diagnostic output
                if (this.IterationCallback != null)
                    this.IterationCallback(0, CurrentSol.CloneAs(), CurrentRes.CloneAs(), this.m_mgop);

                // iterations...
                // =============
                double[] PreviousRes = new double[L];

                //MultidimensionalArray raw_Mxx = MultidimensionalArray.Create(L, MaxIter + 1);
                //MultidimensionalArray ortho_Mxx = MultidimensionalArray.Create(L, MaxIter + 1);

                MultidimensionalArray MassMatrix = MultidimensionalArray.Create(MaxKrylovDim, MaxKrylovDim);

                int PCcounter = 0;
                double[] prevAlpha = null;
                for (int iIter = 0; iIter < MaxIter; iIter++) {
                    Debug.Assert(SolHistory.Count == MxxHistory.Count);
                    Debug.Assert(SolHistory.Count == KrylovDim);

                    // select preconditioner
                    var Precond = PrecondS[PCcounter];
                    PCcounter++;
                    if (PCcounter >= PrecondS.Length) {
                        PCcounter = 0;
                        m_ThisLevelIterations++; // because we abuse the Orthonormalization to do some multi-grid stuff, 
                        //                          we only count every full cycle of preconditiones.
                    }

                    // solve the residual equation: M*Correction = prev. Residual
                    PreviousRes.SetV(CurrentRes);
                    Correction.ClearEntries();
                    Precond.Solve(Correction, PreviousRes);

                    // compute M*Correction
                    Mtx.SpMV(1.0, Correction, 0.0, Mxx);

                    // orthonormalize the Mxx -- vector with respect to the previous ones.
                    Debug.Assert(KrylovDim == MxxHistory.Count);
                    Debug.Assert(KrylovDim == SolHistory.Count);

                    //raw_Mxx.SetColumn(KrylovDim, Mxx);

                    for (int i = 0; i < KrylovDim; i++) {
                        Debug.Assert(!object.ReferenceEquals(Mxx, MxxHistory[i]));
                        double beta = GenericBlas.InnerProd(Mxx, MxxHistory[i]).MPISum();
                        Mxx.AccV(-beta, MxxHistory[i]);
                        Correction.AccV(-beta, SolHistory[i]);
                    }
                    {
                        double gamma = 1.0 / GenericBlas.L2NormPow2(Mxx).MPISum().Sqrt();
                        Mxx.ScaleV(gamma);
                        Correction.ScaleV(gamma);
                    }

                    // the following lines should produce the identity matrix
                    for (int i = 0; i < KrylovDim; i++) {
                        MassMatrix[i, KrylovDim] = GenericBlas.InnerProd(Mxx, MxxHistory[i]).MPISum();
                    }
                    MassMatrix[KrylovDim, KrylovDim] = GenericBlas.L2NormPow2(Mxx).MPISum();
                    //


                    //ortho_Mxx.SetColumn(KrylovDim, Mxx);

                    MxxHistory.Add(Mxx.CloneAs());
                    SolHistory.Add(Correction.CloneAs());
                    KrylovDim++;



                    bool updateEveryIteration = false;


                    // RHS of the minimization problem (LHS is identity matrix)
                    if (!updateEveryIteration) {
                        _R.Add(GenericBlas.InnerProd(MxxHistory.Last(), Residual0).MPISum());
                    } else {
                        _R.Clear();
                        for (int i = 0; i < KrylovDim; i++) {
                            _R.Add(GenericBlas.InnerProd(MxxHistory[i], Residual0).MPISum());
                        }
                    }

                    // compute accelerated solution
                    //double[] alpha = _R.ToArray(); // factors for re-combining solutions
                    double[] alpha;
                    {
                        double[] minimi_rhs = _R.ToArray();
                        var minimi_lhs = MassMatrix.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { KrylovDim - 1, KrylovDim - 1 }).CloneAs();
                        alpha = new double[KrylovDim];
                        minimi_lhs.Solve(alpha, minimi_rhs);
                    }
                    if (prevAlpha != null) {
                        var del = alpha.GetSubVector(0, prevAlpha.Length);
                        del.AccV(-1.0, prevAlpha);
                    }
                    prevAlpha = alpha;

                    Console.WriteLine("Correction factor: " + alpha.Last() + ", solution " + SolHistory.Last().L2Norm() + " Resi " + Mxx.L2Norm());


                    Debug.Assert(alpha.Length == SolHistory.Count);
                    Debug.Assert(alpha.Length == MxxHistory.Count);
                    Debug.Assert(alpha.Length == KrylovDim);
                    CurrentSol.SetV(Solution0, 1.0);
                    for (int i = 0; i < KrylovDim; i++)
                        CurrentSol.AccV(alpha[i], SolHistory[i]);

                    // compute new Residual
                    CurrentRes.SetV(B);
                    Mtx.SpMV(-1.0, CurrentSol, 1.0, CurrentRes);
                    double crL2 = CurrentRes.L2Norm();

                    // diagnostic output
                    if (this.IterationCallback != null)
                        this.IterationCallback(iIter + 1, CurrentSol.CloneAs(), CurrentRes.CloneAs(), this.m_mgop);

                    //{
                    //    var gdat = m_mgop.BaseGridProblemMapping.GridDat;
                    //    var basis = m_mgop.BaseGridProblemMapping.BasisS[0];

                    //    var dgCurrentSol = new Foundation.SinglePhaseField(basis, "Solution");
                    //    var dgResidual = new Foundation.SinglePhaseField(basis, "Residual");
                    //    var dgCorrection = new Foundation.SinglePhaseField(basis, "Correction");

                    //    m_mgop.TransformRhsFrom(dgResidual.CoordinateVector, CurrentRes);
                    //    m_mgop.TransformSolFrom(dgCurrentSol.CoordinateVector, CurrentSol);
                    //    m_mgop.TransformSolFrom(dgCorrection.CoordinateVector, SolHistory.Last());
                    //    dgCorrection.Scale(alpha.Last());

                    //    Tecplot.Tecplot.PlotFields(new Foundation.DGField[] { dgCurrentSol, dgResidual, dgCorrection},  "OrthoScheme-" + iIter, iIter, 2);

                    //}


                    if (crL2 < Tolerance) {
                        //Console.WriteLine("    Kcy converged:");
                        //for (int iii = 0; iii < KrylovDim; iii++) {
                        //    Console.WriteLine("       fac #" + iii + "  :  " + alpha[iii]);
                        //}
                        if (PCcounter > 0)
                            m_ThisLevelIterations += 1;

                        m_Converged = true;
                        break;
                    }

                    if (updateEveryIteration) {
                        Solution0.SetV(CurrentSol);
                        Residual0.SetV(CurrentRes);
                    }

                    if (KrylovDim >= MaxKrylovDim) {
                        if (this.Restarted) {
                            // restarted version of the algorithm
                            // ++++++++++++++++++++++++++++++++++

                            MxxHistory.Clear();
                            SolHistory.Clear();
                            _R.Clear();
                            KrylovDim = 0;
                            Residual0.SetV(CurrentRes);
                            Solution0.SetV(CurrentSol);
                        } else {
                            // throw-away version of the algorithm
                            // +++++++++++++++++++++++++++++++++++

                            int i_leastSig = alpha.IndexOfMin(x => x.Abs());
                            MxxHistory.RemoveAt(i_leastSig);
                            SolHistory.RemoveAt(i_leastSig);
                            KrylovDim--;

                            for (int i = i_leastSig; i < KrylovDim; i++) {
                                for (int j = 0; j <= KrylovDim; j++) {
                                    MassMatrix[i, j] = MassMatrix[i + 1, j];
                                }
                            }
                            for (int i = i_leastSig; i < KrylovDim; i++) {
                                for (int j = 0; j <= KrylovDim; j++) {
                                    MassMatrix[j, i] = MassMatrix[j, i + 1];
                                }
                            }

                            Residual0.SetV(CurrentRes);
                            Solution0.SetV(CurrentSol);

                            _R.Clear();
                            foreach (double[] mxx in MxxHistory) {
                                _R.Add(GenericBlas.InnerProd(mxx, Residual0).MPISum());
                            }
                        }
                    }
                }


                X.SetV(CurrentSol, 1.0);
                //raw_Mxx.SaveToTextFile("C:\\temp\\raw_Mxx.txt");
                //ortho_Mxx.SaveToTextFile("C:\\temp\\ortho_Mxx.txt");

                */
            }
        }

       
        



        /// <summary>
        /// Maximum number of iterations.
        /// </summary>
        public int MaxIter = 100;

        public double Tolerance = 1E-10;

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
        public ISolverSmootherTemplate Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }
    }
}

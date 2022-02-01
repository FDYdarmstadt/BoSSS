using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// p-Multigrid on a single grid level
    /// </summary>
    public class LevelPmg : ISolverSmootherTemplate, ISolverWithCallback {

        public bool UseDiagonalPmg = true;

        public bool EqualOrder = false;

        /// <summary>
        /// ctor
        /// </summary>
        public LevelPmg() {
            UseHiOrderSmoothing = true;
        }

        /// <summary>
        /// always 0, because there is no nested solver
        /// </summary>
        public int IterationsInNested {
            get {
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
        /// If true, cell-local solvers will be used to approximate a solution to high-order modes
        /// </summary>
        public bool UseHiOrderSmoothing {
            get;
            set;
        }

        public object Clone() {
            throw new NotImplementedException();
        }

        private MultigridOperator m_op;

        /// <summary>
        /// - 1st index: cell
        /// </summary>
        MultidimensionalArray[] HighOrderBlocks_LU;
        int[][] HighOrderBlocks_LUpivots;

        public int OrderOfCoarseSystem {
            get { return m_LowOrder; }
            set { m_LowOrder = value; }
        }


        /// <summary>
        /// DG degree at low order blocks. This degree is the border, which divides into low order and high order blocks
        /// </summary>
        private int m_LowOrder = 1;

        /// <summary>
        /// If true blocks/cells containing more than one species are completely assigned to low order block solver.
        /// This hopefully is better than the default approach
        /// </summary>
        public bool FullSolveOfCutcells {
            get;
            set;
        }

        private bool AnyHighOrderTerms {
            get {
                Debug.Assert(m_op != null, "there is no matrix given yet!");
                return m_op.Mapping.DgDegree.Any(p => p > m_LowOrder);
            }
        }

        /// <summary>
        /// Krankplätze müssen verdichtet werden
        /// -Kranführer Ronny, ProSieben Reportage
        /// </summary>
        public void Init(MultigridOperator op) {

            //            //System.Threading.Thread.Sleep(10000);
            //            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            m_op = op;

            if (m_LowOrder > m_op.Mapping.DgDegree.Max())
                throw new ArgumentOutOfRangeException("CoarseLowOrder is higher than maximal DG degree");

#if TEST
            var debugerSW = new StreamWriter(String.Concat("debug_of_", ilPSP.Environment.MPIEnv.MPI_Rank));
            Console.WriteLine("variable TEST is defined");
            //debugerSW.WriteLine("proc {0} reporting Num of Blocks {1}", ilPSP.Environment.MPIEnv.MPI_Rank, HighOrderBlocks_LUpivots.Length);
#endif

            int D = this.m_op.GridData.SpatialDimension;


            var DGlowSelect = new SubBlockSelector(op.Mapping);
            Func<int, int, int, int, bool> lowFilter = (int iCell, int iVar, int iSpec, int pDeg) => pDeg <= (iVar != D && !EqualOrder ? OrderOfCoarseSystem : OrderOfCoarseSystem - 1); // containd the pressure hack
            DGlowSelect.ModeSelector(lowFilter);

            if (FullSolveOfCutcells)
                ModifyLowSelector(DGlowSelect, op);

            lMask = new BlockMask(DGlowSelect);
            m_lowMaskLen = lMask.GetNoOfMaskedRows;

            if (UseHiOrderSmoothing && AnyHighOrderTerms) {
                var DGhighSelect = new SubBlockSelector(op.Mapping);
                Func<int, int, int, int, bool> highFilter = (int iCell, int iVar, int iSpec, int pDeg) => pDeg > (iVar != D && !EqualOrder ? OrderOfCoarseSystem : OrderOfCoarseSystem - 1);
                DGhighSelect.ModeSelector(highFilter);

                if (FullSolveOfCutcells)
                    ModifyHighSelector(DGhighSelect, op);

                hMask = new BlockMask(DGhighSelect);
                m_highMaskLen = hMask.GetNoOfMaskedRows;

                BlockMsrMatrix P01HiMatrix = null;

                if (UseDiagonalPmg) {
                    HighOrderBlocks_LU = hMask.GetDiagonalBlocks(op.OperatorMatrix, false, false);
                    int NoOfBlocks = HighOrderBlocks_LU.Length;
                    HighOrderBlocks_LUpivots = new int[NoOfBlocks][];

                    for (int jLoc = 0; jLoc < NoOfBlocks; jLoc++) {
                        int len = HighOrderBlocks_LU[jLoc].NoOfRows;
                        HighOrderBlocks_LUpivots[jLoc] = new int[len];
                        HighOrderBlocks_LU[jLoc].FactorizeLU(HighOrderBlocks_LUpivots[jLoc]);
                    }
                } else {
                    P01HiMatrix = hMask.GetSubBlockMatrix(op.OperatorMatrix, csMPI.Raw._COMM.SELF);

                    hiSolver = new PARDISOSolver() {
                        CacheFactorization = true,
                        UseDoublePrecision = true, // keep it true, experiments showed, that this leads to fewer iterations
                        Parallelism = Parallelism.OMP
                    };
                    hiSolver.DefineMatrix(P01HiMatrix);
                }
            }

            var P01SubMatrix = lMask.GetSubBlockMatrix(op.OperatorMatrix, csMPI.Raw._COMM.WORLD);

            lowSolver = new PARDISOSolver() {
                CacheFactorization = true,
                UseDoublePrecision = false, // no difference towards =true observed for XDGPoisson
                Parallelism = Parallelism.OMP
            };
            lowSolver.DefineMatrix(P01SubMatrix);


            Debug.Assert(UseDiagonalPmg && lowSolver != null);
            Debug.Assert(UseDiagonalPmg || (!UseDiagonalPmg && hiSolver != null));
            Debug.Assert(m_lowMaskLen > 0);
            //Debug.Assert(AnyHighOrderTerms && m_highMaskLen > 0);
#if TEST
            P01SubMatrix.SaveToTextFileSparseDebug("lowM");
            P01SubMatrix.SaveToTextFileSparse("lowM_full");
            if (!UseDiagonalPmg) {
                //P01HiMatrix.SaveToTextFileSparseDebug("hiM");
                //P01HiMatrix.SaveToTextFileSparse("hiM_full");
            }
            m_op.OperatorMatrix.SaveToTextFileSparseDebug("M");
            m_op.OperatorMatrix.SaveToTextFileSparse("M_full");
            debugerSW.Flush();
            debugerSW.Close();

            long[] bla = m_op.BaseGridProblemMapping.GridDat.CurrentGlobalIdPermutation.Values;
            bla.SaveToTextFileDebug("permutation_");

            List<int> BlockI0 = new List<int>();
            List<int> Block_N = new List<int>();
            foreach (long Block in bla) {
                BlockI0.Add(m_op.Mapping.GetBlockI0((int)Block));
                Block_N.Add(m_op.Mapping.GetBlockLen((int)Block));
            }
            BlockI0.SaveToTextFileDebug("BlockI0");
            Block_N.SaveToTextFileDebug("Block_N");
#endif
        }

        private void ModifyLowSelector(SubBlockSelector sbs, MultigridOperator op) {
            AssignXdgBlocksModification(sbs, op, true);
        }

        private void ModifyHighSelector(SubBlockSelector sbs, MultigridOperator op) {
            AssignXdgBlocksModification(sbs, op, false);
        }

        private void AssignXdgBlocksModification(SubBlockSelector sbs, MultigridOperator op, bool IsLowSelector) {
            var Filter = sbs.ModeFilter;
            Func<int, int, int, int, bool> Modification = delegate (int iCell, int iVar, int iSpec, int pDeg) {
                int NoOfSpec = op.Mapping.AggBasis[0].GetNoOfSpecies(iCell);
                if (NoOfSpec >= 2)
                    return IsLowSelector;
                else
                    return Filter(iCell, iVar, iSpec, pDeg);
            };
            sbs.ModeSelector(Modification);
        }

        /// <summary>
        /// Solver of low order system.
        /// The low order system is defined by <see cref="OrderOfCoarseSystem"/>
        /// </summary>
        private ISparseSolver lowSolver;

        /// <summary>
        /// experimental, used if <see cref="UseDiagonalPmg"/> is not set.
        /// Then low order and high order blocks are both solved by direct solver.
        /// </summary>
        private ISparseSolver hiSolver;

        int m_Iter = 0;

        private int m_lowMaskLen = 0;
        private int m_highMaskLen = 0;

        private BlockMask hMask;
        private BlockMask lMask;

        /// <summary>
        /// ~
        /// </summary>
        public void ResetStat() {
            Converged = false;
            ThisLevelIterations = 0;
        }

        double[] Res_f;
        double[] Cor_c;

        public bool SkipLowOrderSolve = false;

        /// <summary>
        /// ~
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> // 
        {
            using(var tr = new FuncTrace()) {
                int Lf = m_op.Mapping.LocalLength; // DOF's in entire system
                int Lc = m_lowMaskLen; //             DOF's in low-order system

                if(Res_f == null || Res_f.Length != Lf) {
                    Res_f = new double[Lf];
                }
                if(Cor_c == null || Cor_c.Length != Lc) {
                    Cor_c = new double[Lc];
                }
                var Cor_f = new double[Lf];
                Cor_f.SetV(X);
                var Mtx = m_op.OperatorMatrix;

                // compute fine residual
                Res_f.SetV(B);
                Mtx.SpMV(-1.0, Cor_f, 1.0, Res_f);

                if(!SkipLowOrderSolve) {
                    // project to low-p/coarse
                    double[] Res_c = lMask.GetSubVec(Res_f);

                    // low-p solve
                    lowSolver.Solve(Cor_c, Res_c);

                    // accumulate low-p correction
                    lMask.AccSubVec(Cor_c, Cor_f);

                    // compute residual of low-order solution
                    Res_f.SetV(B);
                    Mtx.SpMV(-1.0, Cor_f, 1.0, Res_f);

                    /*
                    DGField[] afterLow = this.m_op.ProlongateSolToDg(Cor_f, "Correction");
                    Tecplot.Tecplot.PlotFields(afterLow, "LevelPMG-0", 0.0, 2);
                    */
                }


                if(UseHiOrderSmoothing && AnyHighOrderTerms) {
                    // solver high-order 

                    Console.WriteLine("UseDiagonalPmg: " + UseDiagonalPmg);
                    if(UseDiagonalPmg) {
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // solve the high-order blocks diagonally, i.e. use a DENSE direct solver for each cell
                        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        var Map = m_op.Mapping;
                        int NoVars = Map.AggBasis.Length;
                        long j0 = Map.FirstBlock;
                        int J = HighOrderBlocks_LU.Length;
                        int[] degs = m_op.Degrees;
                        var BS = Map.AggBasis;

                        long Mapi0 = Map.i0;
                        double[] x_hi = null;
                        for(int j = 0; j < J; j++) { // loop over cells

                            if(HighOrderBlocks_LU[j] != null) {
                                int NpTotHi = HighOrderBlocks_LU[j].NoOfRows;
                                x_hi = new double[NpTotHi];

                                double[] b_f = hMask.GetSubVecOfCell(Res_f, j);
                                Debug.Assert(b_f.Length == NpTotHi);
                                HighOrderBlocks_LU[j].BacksubsLU(HighOrderBlocks_LUpivots[j], x_hi, b_f);
                                hMask.AccSubVecOfCell(x_hi, j, X);
                            }

                        }
                    } else {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // solver the high-order system at once, using a SPARSE direct solver for all high-order modes
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        if(m_highMaskLen > 0) {
                            int Hc = m_highMaskLen;
                            // project to low-p/coarse
                            double[] hi_Res_c = hMask.GetSubVec(Res_f);
                            Debug.Assert(hi_Res_c.Length == m_highMaskLen);
                            double[] hi_Cor_c = new double[Hc];
                            hiSolver.Solve(hi_Cor_c, hi_Res_c);
                            hMask.AccSubVec(hi_Cor_c, X);
                        }
                    }

                    //compute residual for Callback
                    Res_f.SetV(B);
                    Mtx.SpMV(-1.0, Cor_f, 1.0, Res_f);
                }

                // combine fine mode correction with low order solution
                X.ScaleV(0.5);
                X.AccV(1.0, Cor_f);
                
                
                /*
                DGField[] afterHigh = this.m_op.ProlongateSolToDg(X, "Correction");
                Tecplot.Tecplot.PlotFields(afterHigh, "LevelPMG-1", 0.0, 2);


                double[] exact = new double[Cor_f.Length];
                this.m_op.OperatorMatrix.Solve_Direct(exact, B);

                DGField[] exactDG = this.m_op.ProlongateSolToDg(exact, "Correction");
                Tecplot.Tecplot.PlotFields(exactDG, "LevelPMG-2", 0.0, 2);

                var st = DateTime.Now;
                int lastEll = -1;
                while(true) {
                    int ell = (int)((DateTime.Now - st).TotalSeconds)/5;
                    if(ell > lastEll)
                        Console.WriteLine("infinity loop...");
                    lastEll = ell;
                }
                */

                //IterationCallback?.Invoke(m_Iter, X.ToArray(), Res_f, m_op);

                m_Iter++;
            }
        }

        /// <summary>
        /// Called upon each iteration
        /// </summary>
        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        /// <summary>
        /// 
        /// </summary>
        public void Dispose() {
            if(lowSolver != null) {
                lowSolver.Dispose();
                lowSolver = null;
            }
            if(hiSolver != null) {
                hiSolver.Dispose();
                hiSolver = null;
            }


        }

        /// <summary>
        /// 
        /// </summary>
        public long UsedMemory() {
            long r = 0;

            foreach(var mda in this.HighOrderBlocks_LU) {
                if(mda != null) {
                    r += mda.Length * sizeof(double);
                }
            }

            foreach(var ia in this.HighOrderBlocks_LUpivots) {
                if(ia != null) {
                    r += ia.Length * sizeof(int);
                }
            }


            return r;
        }
    }
}

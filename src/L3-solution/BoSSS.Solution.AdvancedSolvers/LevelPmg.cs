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
    public class LevelPmg : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination, ISubsystemSolver {

        /// <summary>
        /// Configuration
        /// </summary>
        [Serializable]
        public class Config : ISolverFactory {
            public string Name => "p-two-grid";

            public string Shortname => "pTG";


            /// <summary>
            /// If true, the high order system is solved cell-by-cell (i.e. a Block-Jacobi)
            /// </summary>
            public bool UseDiagonalPmg = true;

      
            /// <summary>
            /// If true, cell-local solvers will be used to approximate a solution to high-order modes
            /// </summary>
            public bool UseHiOrderSmoothing = true;

            /// <summary>
            /// DG degree at low order blocks. This degree is the border, which divides into low order and high order blocks
            /// </summary>
            public int OrderOfCoarseSystem = 1;

            /// <summary>
            /// If true blocks/cells containing more than one species are completely assigned to low order block solver.
            /// This hopefully is better than the default approach
            /// </summary>
            public bool FullSolveOfCutcells = true;


            /// <summary>
            /// 
            /// </summary>
            public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                var instance = new LevelPmg();
                instance.m_Config = this;
                instance.Init(level);
                return instance;
            }

            /// <summary>
            /// 
            /// </summary>
            public bool Equals(ISolverFactory other) {
                return EqualsImpl(other);
            }

            /// <summary>
            /// 
            /// </summary>
            public override bool Equals(object obj) {
                return EqualsImpl(obj);
            }

            private bool EqualsImpl(object o) {
                var other = o as Config;

                return (this.UseDiagonalPmg == other.UseDiagonalPmg)
                    && (this.FullSolveOfCutcells == other.FullSolveOfCutcells)
                    && (this.OrderOfCoarseSystem == other.OrderOfCoarseSystem)
                    && (this.UseHiOrderSmoothing == other.UseHiOrderSmoothing);
            }

            public override int GetHashCode() {
                return this.OrderOfCoarseSystem;
            }

        }

        Config m_Config = new Config();

        /// <summary>
        /// configuration
        /// </summary>
        public Config config {
            get {
                return m_Config;
            }
        }



        /// <summary>
        /// ctor
        /// </summary>
        public LevelPmg() {
            TerminationCriterion = (int iter, double r0, double r) => (iter < 1, true);
        }

        /// <summary>
        /// User-Programmable termination criterion: 
        /// - 1st argument: iteration index
        /// - 2nd argument: l2-Norm of residual of initial solution 
        /// - 3rd argument: l2-Norm of residual of solution in current iteration
        /// - return value, 1st item: true to continue, false to terminate
        /// - return value, 2nd item: true for successful convergence (e.g. convergence criterion reached), false for failure (e.g. maximum number of iterations reached)
        /// </summary>
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get;
            set;
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

       

        public object Clone() {
            throw new NotImplementedException();
        }

        private IOperatorMappingPair m_op;

        /// <summary>
        /// - 1st index: cell
        /// </summary>
        MultidimensionalArray[] HighOrderBlocks_LU;
        int[][] HighOrderBlocks_LUpivots;


        MultidimensionalArray[] HighOrderBlocks;


        private bool AnyHighOrderTerms {
            get {
                Debug.Assert(m_op != null, "there is no matrix given yet!");
                return m_op.DgMapping.DgDegree.Any(p => p > config.OrderOfCoarseSystem);
            }
        }


        int[] GetBestFitLowOrder(int pLow) {
            if (m_op is MultigridOperator _op) {
                var _degs = _op.DGpolynomialDegreeHierarchy;

                int pBestDist = int.MaxValue;
                int iBest = -1;
                for (int i = 0; i < _degs.Length; i++) {
                    int pMax = _degs[i].Max();

                    if (pMax == pLow)
                        return _degs[i];

                    int pdist = Math.Abs(pLow - pMax);
                    if (pdist <= pBestDist) {
                        pBestDist = pdist;
                        iBest = i;
                    }
                }

                return _degs[iBest];
            } else {
                var _degs = new List<int[]>();
                _degs.Add(m_op.DgMapping.DgDegree);

                int pBestDist = Math.Abs(_degs[0].Max() - pLow);
                int iBest = 0;
                for (int i = 1; i <= _degs[0].Max(); i++) {
                    int[] _degs_i = _degs[i - 1].Select(k => k - 1).ToArray();
                    int pMax = _degs_i.Max();
                    if (_degs_i.Min() < 0)
                        break;
                    _degs.Add(_degs_i);

                    if (pMax == pLow)
                        return _degs[i];

                    int pdist = Math.Abs(pLow - pMax);
                    if (pdist <= pBestDist) {
                        pBestDist = pdist;
                        iBest = i;
                    }

                }

                return _degs[iBest];

            }
        }

        /// <summary>
        /// 
        /// </summary>
        void InitImpl(IOperatorMappingPair op) {
            if (object.ReferenceEquals(op, m_op))
                return;
            if (m_op != null)
                this.Dispose();

            //            //System.Threading.Thread.Sleep(10000);
            //            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            m_op = op;

            if (config.OrderOfCoarseSystem > m_op.DgMapping.DgDegree.Max())
                throw new ArgumentOutOfRangeException("CoarseLowOrder is higher than maximal DG degree");



            //int D = this.m_op.GridData.SpatialDimension;

            int[] lowDegs = GetBestFitLowOrder(config.OrderOfCoarseSystem);
            bool LowSelector(int iCell, int iVar, int iSpec, int pDeg) {
                return pDeg <= lowDegs[iVar];
            }



            var DGlowSelect = new SubBlockSelector(op.DgMapping);
            //Func<int, int, int, int, bool> lowFilter = delegate (int iCell, int iVar, int iSpec, int pDeg) {
            //    return pDeg <= (iVar != D && !config.EqualOrder ? config.OrderOfCoarseSystem : config.OrderOfCoarseSystem - 1); // containd the pressure hack
            //};
            DGlowSelect.SetModeSelector(LowSelector);

            if (config.FullSolveOfCutcells)
                ModifyLowSelector(DGlowSelect, op);

            lMask = new BlockMask(DGlowSelect);
            int m_lowMaskLen = lMask.NoOfMaskedRows;

            if (config.UseHiOrderSmoothing && AnyHighOrderTerms) {
                var DGhighSelect = new SubBlockSelector(op.DgMapping);
                Func<int, int, int, int, bool> highFilter = (int iCell, int iVar, int iSpec, int pDeg) => !LowSelector(iCell, iVar, iSpec, pDeg);
                //Func<int, int, int, int, bool> highFilter = (int iCell, int iVar, int iSpec, int pDeg) => pDeg > (iVar != D && !config.EqualOrder ? config.OrderOfCoarseSystem : config.OrderOfCoarseSystem - 1);
                //Func<int, int, int, int, bool> highFilter = (int iCell, int iVar, int iSpec, int pDeg) => pDeg >= 0;
                DGhighSelect.SetModeSelector(highFilter);

                if (config.FullSolveOfCutcells)
                    ModifyHighSelector(DGhighSelect, op);

                hMask = new BlockMask(DGhighSelect);
                int m_highMaskLen = hMask.NoOfMaskedRows;

                BlockMsrMatrix P01HiMatrix = null;

                if (config.UseDiagonalPmg) {
                    HighOrderBlocks_LU = hMask.GetDiagonalBlocks(op.OperatorMatrix, false, false);
                    HighOrderBlocks = HighOrderBlocks_LU.Select(b => b.CloneAs()).ToArray();
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
            if (P01SubMatrix.MPI_Comm != op.OperatorMatrix.MPI_Comm)
                throw new ApplicationException("MPI comm messup");

            lowSolver = new PARDISOSolver() {
                CacheFactorization = true,
                UseDoublePrecision = false, // no difference towards =true observed for XDGPoisson

                // fk, 24jun22, i have confirmed, at least on my laptop, that serial is faster when only used as some local block solver
                Parallelism = (op.OperatorMatrix.MPI_Comm == csMPI.Raw._COMM.SELF) ? Parallelism.SEQ : Parallelism.OMP
            };
            lowSolver.DefineMatrix(P01SubMatrix);


            Debug.Assert(config.UseDiagonalPmg && lowSolver != null);
            Debug.Assert(config.UseDiagonalPmg || (!config.UseDiagonalPmg && hiSolver != null));
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

        private void ModifyLowSelector(SubBlockSelector sbs, IOperatorMappingPair op) {
            AssignXdgBlocksModification(sbs, op, true);
        }

        private void ModifyHighSelector(SubBlockSelector sbs, IOperatorMappingPair op) {
            AssignXdgBlocksModification(sbs, op, false);
        }

        private void AssignXdgBlocksModification(SubBlockSelector sbs, IOperatorMappingPair op, bool IsLowSelector) {
            var Filter = sbs.ModeFilter;
            //var Mask = (op.BaseGridProblemMapping.BasisS[0] as XDGBasis).Tracker.Regions.GetNearFieldMask(0);
            //var bMask = Mask.GetBitMask();
            //Console.WriteLine($"Fine solution in {Mask.NoOfItemsLocally} of {Mask.GridData.iLogicalCells.NoOfLocalUpdatedCells} cells.");

            Func<int, int, int, int, bool> Modification = delegate (int iCell, int iVar, int iSpec, int pDeg) {
                //if(bMask[iCell])
                //    return true;

                int NoOfSpec = op.DgMapping.GetNoOfSpecies(iCell);
                if (NoOfSpec >= 2)
                    return IsLowSelector;
                else
                    return Filter(iCell, iVar, iSpec, pDeg);
            };
            sbs.SetModeSelector(Modification);
        }

        /// <summary>
        /// Solver of low order system.
        /// The low order system is defined by <see cref="Config.OrderOfCoarseSystem"/>
        /// </summary>
        private ISparseSolver lowSolver;

        /// <summary>
        /// experimental, used if <see cref="Config.UseDiagonalPmg"/> is not set.
        /// Then low order and high order blocks are both solved by direct solver.
        /// </summary>
        private ISparseSolver hiSolver;

        int m_Iter = 0;
                

        private BlockMask hMask;
        private BlockMask lMask;

        /// <summary>
        /// ~
        /// </summary>
        public void ResetStat() {
            Converged = false;
            ThisLevelIterations = 0;
        }


        /// <summary>
        /// Computes the coarse-grid correction
        /// </summary>
        /// <param name="x_out">output: coarse level solution, prolongated to fine level</param>
        /// <param name="in_rhs">input: RHS on fine level</param>
        void CoarseSolve(double[] x_out, double[] in_rhs) {
            // project to low-p/coarse
            double[] rhs_c = lMask.GetSubVec(in_rhs);

            // low-p solve
            double[] x_c = new double[rhs_c.Length];
            lowSolver.Solve(x_c, rhs_c);

            // accumulate low-p correction
            lMask.AccSubVec(x_c, x_out);

            //// test: if we use an exact solution, we should terminate in one iteration!
            //double[] xtest = new double[in_rhs.Length];
            //m_op.OperatorMatrix.Solve_Direct(xtest, in_rhs);
            //x_out.AccV(1.0, xtest);

        }

        /// <summary>
        /// smoothing/solving on high level
        /// </summary>
        void FineSolve(double[] x_in_out, double[] in_rhs) {
            using(var tr = new FuncTrace()) {


                // compute residual
                double[] Res_f = in_rhs.CloneAs();
                this.m_op.OperatorMatrix.SpMV(-1.0, x_in_out, 1.0, Res_f);


                // test: if we use an exact solution, we should terminate in one iteration!
                /*
                double[] xtest = new double[in_rhs.Length];
                m_op.OperatorMatrix.Solve_Direct(xtest, Res_f);
                x_in_out.AccV(1.0, xtest);
                Res_f.ClearEntries();
                */
    

                if(config.UseHiOrderSmoothing && AnyHighOrderTerms) {
                    // solver high-order 

                    tr.Info("UseDiagonalPmg: " + config.UseDiagonalPmg);
                    if(config.UseDiagonalPmg) {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // solve the high-order blocks diagonally, i.e. use a DENSE direct solver LOCALLY IN EACH CELL
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        int J = HighOrderBlocks_LU.Length;

                        for(int j = 0; j < J; j++) { // loop over cells

                            if(HighOrderBlocks_LU[j] != null) {
                                int NpTotHi = HighOrderBlocks_LU[j].NoOfRows;
                                var x_hi = new double[NpTotHi];

                                double[] b_f = hMask.GetSubVecOfCell(Res_f, j);
                                Debug.Assert(b_f.Length == NpTotHi);
                                HighOrderBlocks_LU[j].BacksubsLU(HighOrderBlocks_LUpivots[j], x_hi, b_f);
                                hMask.AccSubVecOfCell(x_hi, j, x_in_out);
                            }
                        }
                    } else {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        // solver the high-order system at once, using a SPARSE direct solver for all high-order modes
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                        int Hc = (hMask?.NoOfMaskedRows) ?? 0;

                        if(Hc > 0) {
                            // project to low-p/coarse
                            double[] hi_Res_c = hMask.GetSubVec(Res_f);
                            Debug.Assert(hi_Res_c.Length == Hc);
                            double[] hi_Cor_c = new double[Hc];
                            hiSolver.Solve(hi_Cor_c, hi_Res_c);
                            hMask.AccSubVec(hi_Cor_c, x_in_out);
                        }
                    }
                }
            }
        }


        /// <summary>
        /// Algorithm 3 in:
        /// 
        /// p-Multigrid matrix-free discontinuous Galerkin solution strategies for the under-resolved simulation of incompressible turbulent flows
        /// M. Franciolini, L. Botti, A. Colombo, A. Crivellini
        /// </summary>
        double[] MGfull(int l, double[] gl, double[] wl) {
            if(l >= 1) {
                // 
                throw new Exception("should not happen");
            } else {
                double[] wl_hat = new double[gl.Length];
                CoarseSolve(wl_hat, gl.CloneAs());

                var dl = gl.CloneAs();
                m_op.OperatorMatrix.SpMV(-1.0, wl_hat, 1.0, dl);

                var el = MGv(l, dl, new double[gl.Length]);

                // wl_dash = wl_hat + el
                var wl_dash = el.CloneAs();
                wl_dash.AccV(1, wl_hat);

                return wl_dash;
            }
        }

        /// <summary>
        /// Algorithm 2 in:
        /// 
        /// p-Multigrid matrix-free discontinuous Galerkin solution strategies for the under-resolved simulation of incompressible turbulent flows
        /// M. Franciolini, L. Botti, A. Colombo, A. Crivellini
        /// </summary>
        double[] MGv(int l, double[] gl, double[] wl) {
            if(l >= 1) {
                throw new Exception("should not happen");
            } else {
                var wl_dash = wl.CloneAs();
                FineSolve(wl_dash, gl);

                var dl = gl.CloneAs();
                m_op.OperatorMatrix.SpMV(-1.0, wl_dash, 1.0, dl);

                double[] el = new double[wl.Length];
                CoarseSolve(el, dl);
                double[] wl_hat = wl_dash.CloneAs();
                wl_hat.AccV(1.0, el);


                FineSolve(wl_hat, gl);
                wl_dash = wl_hat;
                return wl_dash;
            }
        }


        /// <summary>
        /// ~
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> // 
        {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                int Lf = m_op.DgMapping.LocalLength; // DOF's in entire system
                //int Lc = this.lMask.NoOfMaskedRows; // DOF's in low-order system


                var Mtx = m_op.OperatorMatrix;

                // compute fine residual: Res_f = B - Mtx*X
                double[] Res_f = new double[Lf];
                Res_f.SetV(B);
                Mtx.SpMV(-1.0, X, 1.0, Res_f);
                
                // solve for coarse correction
                var Cor_f = new double[Lf];
                CoarseSolve(Cor_f, Res_f);

                // accumulate smoothing 
                FineSolve(Cor_f, Res_f);
                
                // solution update
                X.AccV(1.0, Cor_f);
                m_Iter++;
                
              

                if(IterationCallback != null) {
                    Res_f.SetV(B);
                    Mtx.SpMV(-1.0, X, 1.0, Res_f);
                    IterationCallback(m_Iter, X.ToArray(), Res_f, m_op as MultigridOperator);
                }


                //var Res_f_0 = Res_f.CloneAs();

               

                /*
                if(!SkipLowOrderSolve) {
                    // project to low-p/coarse
                    double[] Res_c = lMask.GetSubVec(Res_f);

                    // low-p solve
                    lowSolver.Solve(Cor_c, Res_c);

                    // accumulate low-p correction
                    lMask.AccSubVec(Cor_c, Cor_f);

                    // compute residual of low-order solution (on fine Level)
                    Res_f.SetV(B);
                    Mtx.SpMV(-1.0, Cor_f, 1.0, Res_f); 
                }
                */


                

                m_Iter++;
            }
        }

        /*
        double[] GetVariableDOFs(double[] X, int SelVar) {
            var DGSelect = new SubBlockSelector(m_op.Mapping);
            DGSelect.VariableSelector(SelVar);
                      

            var lMask = new BlockMask(DGSelect);

            var X_SelVar = lMask.GetSubVec(X);

            double[] Ret = new double[X.Length];
            lMask.AccSubVec(X_SelVar, Ret);

            return Ret;
        }
        */
        



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
            HighOrderBlocks = null;
            HighOrderBlocks_LU = null;
            HighOrderBlocks_LUpivots = null;
            lMask = null;
            hMask = null;
        }

        /// <summary>
        /// 
        /// </summary>
        public long UsedMemory() {
            long r = 0;

            if(this.HighOrderBlocks_LUpivots != null) {
                foreach(var mda in this.HighOrderBlocks_LU) {
                    if(mda != null) {
                        r += mda.Length * sizeof(double);
                    }
                }
            }

            if(this.HighOrderBlocks_LU != null) {
                foreach(var ia in this.HighOrderBlocks_LUpivots) {
                    if(ia != null) {
                        r += ia.Length * sizeof(int);
                    }
                }
            }


            var lowPard = this.lowSolver as PARDISOSolver;
            if(lowPard != null) {
                r += lowPard.UsedMemory();
            }

            var hiPard = this.hiSolver as PARDISOSolver;
            if(hiPard != null) {
                r += hiPard.UsedMemory();
            }


            return r;
        }

        public void Init(MultigridOperator op) {
            InitImpl(op);
        }

        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }
    }
}

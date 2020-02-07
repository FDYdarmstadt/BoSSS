using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
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

        /// <summary>
        /// ctor
        /// </summary>
        public LevelPmg() {
            UseHiOrderSmoothing = false;
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

        MultigridOperator m_op;

        /// <summary>
        /// - 1st index: cell
        /// </summary>
        MultidimensionalArray[] HighOrderBlocks_LU;
        int[][] HighOrderBlocks_LUpivots;

        public int CoarseLowOrder {
            get { return m_CoarseLowOrder; }
            set { m_CoarseLowOrder = value; }
        }


        /// <summary>
        /// DG degree at low order blocks. This degree is the border, which divides into low order and high order blocks
        /// </summary>
        private int m_CoarseLowOrder = 1;



        /// <summary>
        /// Krankplätze müssen verdichtet werden
        /// -Kranführer Ronny, ProSieben Reportage
        /// </summary>
        public void Init(MultigridOperator op) {

            //            //System.Threading.Thread.Sleep(10000);
            //            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            m_op = op;

            if (m_CoarseLowOrder > m_op.Mapping.DgDegree.Max())
                throw new ArgumentOutOfRangeException("CoarseLowOrder is higher than maximal DG degree");

#if TEST
            var debugerSW = new StreamWriter(String.Concat("debug_of_", ilPSP.Environment.MPIEnv.MPI_Rank));
            Console.WriteLine("variable TEST is defined");
            //debugerSW.WriteLine("proc {0} reporting Num of Blocks {1}", ilPSP.Environment.MPIEnv.MPI_Rank, HighOrderBlocks_LUpivots.Length);
#endif

            int D = this.m_op.GridData.SpatialDimension;
           

            var DGlowSelect = new SubBlockSelector(op.Mapping);
            DGlowSelect.ModeSelector((int iCell, int iVar, int iSpec, int pDeg) => pDeg <= (iVar != D ? m_CoarseLowOrder : m_CoarseLowOrder - 1)); // dirty hack for mixed order stokes
            lMask = new BlockMask(DGlowSelect,false);
            m_lowMaskLen = lMask.GetNoOfMaskedRows;

            if (UseHiOrderSmoothing) {
                var DGhighSelect = new SubBlockSelector(op.Mapping);
                DGhighSelect.ModeSelector((int iCell, int iVar, int iSpec, int pDeg) => pDeg > (iVar != D ? m_CoarseLowOrder : m_CoarseLowOrder - 1));
                hMask = new BlockMask(DGhighSelect,false);
                m_highMaskLen = hMask.GetNoOfMaskedRows;

                BlockMsrMatrix P01HiMatrix = null;

                if (UseDiagonalPmg) {
                    HighOrderBlocks_LU = hMask.GetSubBlocks(op.OperatorMatrix, true, false, false);
                    int NoOfBlocks = HighOrderBlocks_LU.Length;
                    HighOrderBlocks_LUpivots = new int[NoOfBlocks][];

                    for (int jLoc = 0; jLoc < NoOfBlocks; jLoc++) {
                        int len = HighOrderBlocks_LU[jLoc].NoOfRows;
                        HighOrderBlocks_LUpivots[jLoc] = new int[len];
                        HighOrderBlocks_LU[jLoc].FactorizeLU(HighOrderBlocks_LUpivots[jLoc]);
                    }
                } else {
                    P01HiMatrix = hMask.GetSubBlockMatrix(op.OperatorMatrix, true, false, true);

                    hiSolver = new PARDISOSolver() {
                        CacheFactorization = true,
                        UseDoublePrecision = true, // keep it true, experiments showed, that this leads to fewer iterations
                        SolverVersion = Parallelism.OMP
                    };
                    hiSolver.DefineMatrix(P01HiMatrix);
                }
            }

            var P01SubMatrix = lMask.GetSubBlockMatrix(op.OperatorMatrix);

            intSolver = new PARDISOSolver() {
                CacheFactorization = true,
                UseDoublePrecision = false, // no difference towards =true observed for XDGPoisson
            };
            intSolver.DefineMatrix(P01SubMatrix);
            
#if TEST
            P01SubMatrix.SaveToTextFileSparseDebug("lowM");
            P01SubMatrix.SaveToTextFileSparse("lowM_full");
            if (!UseDiagonalPmg) {
                P01HiMatrix.SaveToTextFileSparseDebug("hiM");
                P01HiMatrix.SaveToTextFileSparse("hiM_full");
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

        ISparseSolver intSolver;
        ISparseSolver hiSolver;

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

        /// <summary>
        /// ~
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> // 
        {
            int Lf = m_op.Mapping.LocalLength;
            int Lc = m_lowMaskLen;

            if (Res_f == null || Res_f.Length != Lf) {
                Res_f = new double[Lf];
            }
            if (Cor_c == null || Cor_c.Length != Lc) {
                Cor_c = new double[Lc];
            }
            var Cor_f = new double[Lf];
            Cor_f.SetV(X);
            var Mtx = m_op.OperatorMatrix;

            // compute fine residual
            Res_f.SetV(B);
            Mtx.SpMV(-1.0, Cor_f, 1.0, Res_f);

            // project to low-p/coarse
            double[] Res_c = lMask.GetSubBlockVec(Res_f);

            // low-p solve
            intSolver.Solve(Cor_c, Res_c);

            // accumulate low-p correction
            lMask.AccVecToFull(Cor_c, Cor_f);

            if (UseHiOrderSmoothing) {
                // solver high-order 
                // compute residual of low-order solution
                Res_f.SetV(B);
                Mtx.SpMV(-1.0, Cor_f, 1.0, Res_f);
                if (UseDiagonalPmg) {
                    var Map = m_op.Mapping;
                    int NoVars = Map.AggBasis.Length;
                    int j0 = Map.FirstBlock;
                    int J = Map.LocalNoOfBlocks;
                    int[] degs = m_op.Degrees;
                    var BS = Map.AggBasis;

                    int Mapi0 = Map.i0;
                    double[] x_hi = null;
                    for (int j = 0; j < J; j++) {

                        if (HighOrderBlocks_LU[j] != null) {
                            int NpTotHi = HighOrderBlocks_LU[j].NoOfRows;
                            x_hi = new double[NpTotHi];

                            double[] b_f = hMask.GetVectorCellwise(Res_f, j);
                            Debug.Assert(b_f.Length == NpTotHi);
                            HighOrderBlocks_LU[j].BacksubsLU(HighOrderBlocks_LUpivots[j], x_hi, b_f);
                            hMask.AccVecCellwiseToFull(x_hi, j, X);
                        }

                    }
                } else {
                    if (m_highMaskLen > 0) {
                        int Hc = m_highMaskLen;
                        // project to low-p/coarse
                        double[] hi_Res_c = hMask.GetSubBlockVec(Res_f);
                        Debug.Assert(hi_Res_c.Length == m_highMaskLen);
                        double[] hi_Cor_c = new double[Hc];
                        hiSolver.Solve(hi_Cor_c, hi_Res_c);
                        hMask.AccVecToFull(hi_Cor_c, X);
                    }
                }
            }

            //compute residual for Callback
            Res_f.SetV(B);
            Mtx.SpMV(-1.0, Cor_f, 1.0, Res_f);

            X.AccV(1.0,Cor_f);

            if (IterationCallback != null) {
                IterationCallback(m_Iter, X.ToArray(), Res_f, m_op);
            }
            m_Iter++;
        }
        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }
    }
}

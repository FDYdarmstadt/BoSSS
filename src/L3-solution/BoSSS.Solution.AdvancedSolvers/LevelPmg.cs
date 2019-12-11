//#define TEST

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
    /// Utility class for visualization of intermediate results.
    /// </summary>
    internal class MGViz {

        public MGViz(MultigridOperator op) {
            m_op = op;
        }

        MultigridOperator m_op;

        public int FindLevel(int L) {
            int iLv = 0;
            for (var Op4Level = m_op.FinestLevel; Op4Level != null; Op4Level = Op4Level.CoarserLevel) {
                if (L == Op4Level.Mapping.LocalLength) {
                    Debug.Assert(Op4Level.LevelIndex == iLv);
                    return iLv;
                }
                iLv++;
            }
            return -1;
        }

        public DGField ProlongateToDg(double[] V, string name) {
            double[] Curr = ProlongateToTop(V);

            var gdat = m_op.BaseGridProblemMapping.GridDat;
            var basis = m_op.BaseGridProblemMapping.BasisS[0];
            DGField dgCurrentSol;
            if (basis is XDGBasis)
                dgCurrentSol = new XDGField((XDGBasis)basis, name);
            else 
                dgCurrentSol = new SinglePhaseField(basis, name);
            m_op.FinestLevel.TransformSolFrom(dgCurrentSol.CoordinateVector, Curr);
            return dgCurrentSol;
        }

        public double[] ProlongateToTop(double[] V) {
            int iLv = FindLevel(V.Length);

            MultigridOperator op_iLv = m_op.FinestLevel;
            for (int i = 0; i < iLv; i++)
                op_iLv = op_iLv.CoarserLevel;
            Debug.Assert(op_iLv.LevelIndex == iLv);
            Debug.Assert(V.Length == op_iLv.Mapping.LocalLength);

            double[] Curr = V;
            for (var Op4Level = op_iLv; Op4Level.FinerLevel != null; Op4Level = Op4Level.FinerLevel) {
                double[] Next = new double[Op4Level.FinerLevel.Mapping.LocalLength];
                Op4Level.Prolongate(1.0, Next, 0.0, Curr);
                Curr = Next;
            }

            return Curr;
        }

        int counter = 0;
        public void PlotVectors(IEnumerable<double[]> VV, string[] names) {

            DGField[] all = new DGField[names.Length];
            for (int i = 0; i < VV.Count(); i++) {
                all[i] = ProlongateToDg(VV.ElementAt(i), names[i]);
            }

            Tecplot.Tecplot.PlotFields(all, "MGviz-" + counter, counter, 2);
            counter++;
        }
    }

    public class MultiIndexBlocking {
        public struct MultiIndex {
            int iBlock;
            int iVariable;
            int iSpecies;
        }
        public MultiIndexBlocking(MultigridOperator mop) {
            var bla = new MultiIndex();
            
        }
    }

    /// <summary>
    /// p-Multigrid on a single grid level
    /// </summary>
    public class LevelPmg : ISolverSmootherTemplate, ISolverWithCallback {

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
        int[][] HighOrderBlocks_indices;

        //MultidimensionalArray[][] HighLoOrderBlocks;

        /// <summary>
        /// DG degree at low order blocks. This degree is the border, which divides into low order and high order blocks
        /// </summary>
        public int CoarseLowOrder = 1;

        /// <summary>
        /// Krankboden muss verdichtet werden
        /// -Kranführer Ronny, ProSieben Reportage
        /// </summary>
        public void Init(MultigridOperator op) {

            //System.Threading.Thread.Sleep(10000);
            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
            m_op = op;

            var Map = op.Mapping;
            int NoVars = Map.AggBasis.Length;
            int j0 = Map.FirstBlock;
            int J = Map.LocalNoOfBlocks;
            int[] degs = m_op.Degrees;

            //for (int ideg = 0; ideg < degs.Length; ideg++)
                //if (degs[ideg] == 0)
                //    throw new ArgumentException(String.Format("DGdegree for Variable {0} is 0, p-multigrid not possible", ideg));

            var BS = Map.AggBasis;

            if (UseHiOrderSmoothing) {
                HighOrderBlocks_LU = new MultidimensionalArray[J];
                HighOrderBlocks_LUpivots = new int[J][];
                HighOrderBlocks_indices = new int[J][];
            }


            var GsubIdx = new List<int>();
            var LsubIdx = new List<int>();
            var HiIdxdummy = new List<int>();
            var GhiIdx = new List<int>();

            var lowLocalBlocks_i0 = new List<int>();
            var lowLocalBlocks__N = new List<int>();
            var hiLocalBlocks_i0 = new List<int>();
            var hiLocalBlocks__N = new List<int>();

            int cnt = 0;
            int NpHiTot = 0;

#if TEST
            var debugerSW = new StreamWriter(String.Concat("debug_of_", ilPSP.Environment.MPIEnv.MPI_Rank));
            debugerSW.WriteLine("proc {0} reporting Num of Blocks {1}",ilPSP.Environment.MPIEnv.MPI_Rank, HighOrderBlocks_LUpivots.Length);   
#endif

            for (int jLoc = 0; jLoc < J; jLoc++) {
                lowLocalBlocks_i0.Add(cnt);
                lowLocalBlocks__N.Add(0);

                var LhiIdx = new List<int>();
                var IdxHighBlockOffset = new int[NoVars][];
                var IdxHighOffset = new int[NoVars][];

                int NpHiLoc = 0;
                for (int iVar = 0; iVar < NoVars; iVar++) {

                    int pReq;

                    if (degs[iVar] <= 1)
                        pReq = 0;
                    else
                        pReq = CoarseLowOrder;

                    int NoOfSpc = BS[iVar].GetNoOfSpecies(jLoc);
                    int Np1 = BS[iVar].GetLength(jLoc, pReq);
                    int Np = BS[iVar].GetLength(jLoc, degs[iVar]);
                    
                    int NpBase = Np / NoOfSpc;
                    int NpBaseLow = Np1 / NoOfSpc;
                    IdxHighBlockOffset[iVar] = new int[NoOfSpc+1];
                    IdxHighOffset[iVar] = new int[NoOfSpc];

                    //this is a hack, blocks with levelset are considered as low order blocks
                    //if (NoOfSpc > 1) {
                    //    NpBaseLow = NpBase;
                    //    lowLocalBlocks__N[jLoc] += Np;
                    //} else {
                        lowLocalBlocks__N[jLoc] += Np1;
                    //}

                for (int iSpc = 0; iSpc < NoOfSpc; iSpc++)
                    {
                        int n = 0;
                        int cellOffset = NpBase * iSpc;

                        if (NpBase == NpBaseLow) {
                            IdxHighOffset[iVar][iSpc] = Map.GlobalUniqueIndex(iVar, jLoc, 0); //there are no high order blocks
                        } else {
                            IdxHighOffset[iVar][iSpc] = Map.GlobalUniqueIndex(iVar, jLoc, cellOffset + NpBaseLow);
                        }

                        for (; n < NpBaseLow; n++)
                        {
                            int Lidx = Map.LocalUniqueIndex(iVar, jLoc, n + cellOffset);
                            LsubIdx.Add(Lidx); //local block mapping Coarse Matrix (low order entries) to original matrix

                            int Gidx = Map.GlobalUniqueIndex(iVar, jLoc, n + cellOffset);
                            GsubIdx.Add(Gidx); //global mapping Coarse Matrix (low order entries) to original matrix
#if TEST
                            if(NoOfSpc==2)
                                debugerSW.WriteLine("{0}", Lidx);
#endif
                        }

                        for (; n < NpBase; n++)
                        {
                            int Lidx = Map.LocalUniqueIndex(iVar, jLoc, n + cellOffset);
                            LhiIdx.Add(Lidx); //local block mapping of high order entries to original matrix


                            int Gidx = Map.GlobalUniqueIndex(iVar, jLoc, n + cellOffset);
                            GhiIdx.Add(Gidx);
#if TEST
                            if (NoOfSpc == 2)
                                debugerSW.WriteLine("{0}", Lidx);
#endif
                        }

                        IdxHighBlockOffset[iVar][iSpc] = NpHiLoc;
                        NpHiLoc += (NpBase- NpBaseLow);
                        
                    }
                    IdxHighBlockOffset[iVar][NoOfSpc] = NpHiLoc; //probably not so nice: necessary to calculate block length

                    Debug.Assert(GsubIdx.Last() == Map.GlobalUniqueIndex(iVar, jLoc, (NoOfSpc-1) * NpBase + NpBaseLow - 1));
                    Debug.Assert(LhiIdx.Count==0 || LhiIdx.Last()==Map.LocalUniqueIndex(iVar, jLoc, NoOfSpc * NpBase-1));
                }


                // Save high order blocks for later smoothing
                if (NpHiLoc > 0) {
                    if (true) {
                        HighOrderBlocks_LU[jLoc] = MultidimensionalArray.Create(NpHiLoc, NpHiLoc);
                        HighOrderBlocks_indices[jLoc] = LhiIdx.ToArray();

                        for (int iVar = 0; iVar < NoVars; iVar++) {
                            for (int jVar = 0; jVar < NoVars; jVar++) {

                                for (int iSpc = 0; iSpc < BS[jVar].GetNoOfSpecies(jLoc); iSpc++) {

                                    int i0_hi = IdxHighOffset[iVar][iSpc];
                                    int j0_hi = IdxHighOffset[jVar][iSpc];
                                    int Row_i0 = IdxHighBlockOffset[iVar][iSpc];
                                    int Col_i0 = IdxHighBlockOffset[jVar][iSpc];
                                    int Row_ie = IdxHighBlockOffset[iVar][iSpc + 1]-1;
                                    int col_ie = IdxHighBlockOffset[jVar][iSpc + 1]-1;

                                    m_op.OperatorMatrix.ReadBlock(i0_hi, j0_hi,
                                                HighOrderBlocks_LU[jLoc].ExtractSubArrayShallow(new int[] { Row_i0, Col_i0 }, new int[] { Row_ie, col_ie }));

                                }
                            }
                        }

                        HighOrderBlocks_LUpivots[jLoc] = new int[NpHiLoc];
                        HighOrderBlocks_LU[jLoc].FactorizeLU(HighOrderBlocks_LUpivots[jLoc]);
                    } else {
                        hiLocalBlocks_i0.Add(NpHiTot);
                        hiLocalBlocks__N.Add(NpHiLoc);
                        Debug.Assert(NpHiLoc == IdxHighBlockOffset[NoVars - 1][BS[NoVars - 1].GetNoOfSpecies(jLoc)]);
                        HiIdxdummy.AddRange(LhiIdx);
                    }
                }

                NpHiTot += NpHiLoc;
                cnt += lowLocalBlocks__N[jLoc];
            }

            if (false) {
                BlockPartitioning localhiBlocking = new BlockPartitioning(GhiIdx.Count, hiLocalBlocks_i0, hiLocalBlocks__N, Map.MPI_Comm, i0isLocal: true);
                var P01HiMatrix = new BlockMsrMatrix(localhiBlocking);
                op.OperatorMatrix.AccSubMatrixTo(1.0, P01HiMatrix, GhiIdx, default(int[]), GhiIdx, default(int[]));
                hiSolver = new PARDISOSolver() {
                    CacheFactorization = true,
                    UseDoublePrecision = true // keep it true, experiments showed, that this leads to fewer iterations
                };
                hiSolver.DefineMatrix(P01HiMatrix);
                m_HiIdx = HiIdxdummy.ToArray();
            }

            if (NpHiTot == 0)
                Console.WriteLine("ATTENTION: No high order blocks, executing direct solve for whole matrix!!!");

            m_LsubIdx = LsubIdx.ToArray();

            BlockPartitioning localBlocking = new BlockPartitioning(GsubIdx.Count, lowLocalBlocks_i0, lowLocalBlocks__N, Map.MPI_Comm, i0isLocal:true);
            var P01SubMatrix = new BlockMsrMatrix(localBlocking);
            op.OperatorMatrix.AccSubMatrixTo(1.0, P01SubMatrix, GsubIdx, default(int[]), GsubIdx, default(int[]));
            
            
            intSolver = new PARDISOSolver() {
                CacheFactorization = true,
                UseDoublePrecision = false // no difference towards =true observed for XDGPoisson
            };
          
            intSolver.DefineMatrix(P01SubMatrix);
            
#if TEST
            P01SubMatrix.SaveToTextFileSparseDebug("lowM");
            P01SubMatrix.SaveToTextFileSparse("lowM_full");

            LsubIdx.SaveToTextFileDebug("LsubIdx");
            GsubIdx.SaveToTextFileDebug("GsubIdx");
            GhiIdx.SaveToTextFileDebug("GhiIdx");
            m_op.OperatorMatrix.SaveToTextFileSparseDebug("M");
            m_op.OperatorMatrix.SaveToTextFileSparse("M_full");
            debugerSW.Flush();
            debugerSW.Close();
            Console.WriteLine("Length HiIdxdummy: " + HiIdxdummy.Count());
            Console.WriteLine("Length GhiIdx:" + GhiIdx.Count());


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

            //debugerSW.WriteLine("Dim of lowMatrix: {0}",GsubIdx.Count);
#endif
        }

        //BlockMsrMatrix P01SubMatrix;

        ISparseSolver intSolver;
        ISparseSolver hiSolver;

        int m_Iter = 0;

        int[] m_LsubIdx;

        int[] m_HiIdx;

        /// <summary>
        /// ~
        /// </summary>
        public void ResetStat() {
            Converged = false;
            ThisLevelIterations = 0;
        }


        double[] Res_f;

        double[] Res_c;
        double[] Cor_c;

        /// <summary>
        /// ~
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> // 
        {
            int Lf = m_op.Mapping.LocalLength;
            int Lc = m_LsubIdx.Length;

            if (Res_f == null || Res_f.Length != Lf) {
                Res_f = new double[Lf];
            }
            if (Res_c == null || Res_c.Length != Lc) {
                Res_c = new double[Lc];
            } else {
                Res_c.ClearEntries();
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
            Res_c.AccV(1.0, Res_f, default(int[]), m_LsubIdx);

            // low-p solve
            intSolver.Solve(Cor_c, Res_c);

            // accumulate low-p correction
            Cor_f.AccV(1.0, Cor_c, m_LsubIdx, default(int[]));

            // solver high-order 
                // compute residual of low-order solution
                Res_f.SetV(B);
                Mtx.SpMV(-1.0, Cor_f, 1.0, Res_f);
            if (true) {
                var Map = m_op.Mapping;
                int NoVars = Map.AggBasis.Length;
                int j0 = Map.FirstBlock;
                int J = Map.LocalNoOfBlocks;
                int[] degs = m_op.Degrees;
                var BS = Map.AggBasis;

                int Mapi0 = Map.i0;
                double[] b_f = null, x_hi = null;
                for (int j = 0; j < J; j++) {

                    if (HighOrderBlocks_LU[j] != null) {
                        int NpTotHi = HighOrderBlocks_LU[j].NoOfRows;
                        if (b_f == null || b_f.Length != NpTotHi) {
                            b_f = new double[NpTotHi];
                            x_hi = new double[NpTotHi];
                        }

                        ArrayTools.GetSubVector<int[], int[], double>(Res_f, b_f, HighOrderBlocks_indices[j]);
                        HighOrderBlocks_LU[j].BacksubsLU(HighOrderBlocks_LUpivots[j], x_hi, b_f);
                        X.AccV(1.0, x_hi, HighOrderBlocks_indices[j], default(int[]));
                    }

                }
            } else {
                int Hc = m_HiIdx.Length;
                if (m_HiIdx.Length != 0) {
                    // project to low-p/coarse
                    double[] hi_Res_c = new double[Hc];
                    hi_Res_c.AccV(1.0, Res_f, default(int[]), m_HiIdx);
                    double[] hi_Cor_c = new double[Hc];
                    hiSolver.Solve(hi_Cor_c, hi_Res_c);
                    X.AccV(1.0, hi_Cor_c, m_HiIdx, default(int[]));
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

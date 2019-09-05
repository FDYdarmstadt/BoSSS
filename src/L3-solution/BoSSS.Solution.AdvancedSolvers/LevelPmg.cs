using BoSSS.Foundation;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
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

        public SinglePhaseField ProlongateToDg(double[] V, string name) {
            double[] Curr = ProlongateToTop(V);

            var gdat = m_op.BaseGridProblemMapping.GridDat;
            var basis = m_op.BaseGridProblemMapping.BasisS[0];
            var dgCurrentSol = new SinglePhaseField(basis, name);
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



    /// <summary>
    /// p-Multigrid on a single grid level
    /// </summary>
    public class LevelPmg : ISolverSmootherTemplate {

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


        public ISolverSmootherTemplate Clone() {
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
        /// ~
        /// </summary>
        public void Init(MultigridOperator op) {
            //var Mtx = op.OperatorMatrix;
            m_op = op;

            var Map = op.Mapping;
            int NoVars = Map.AggBasis.Length;
            int j0 = Map.FirstBlock;
            int J = Map.LocalNoOfBlocks;
            int[] degs = m_op.Degrees;

            var BS = Map.AggBasis;

            if (UseHiOrderSmoothing) {
                HighOrderBlocks_LU = new MultidimensionalArray[J];
                HighOrderBlocks_LUpivots = new int[J][];
                HighOrderBlocks_indices = new int[J][];
            }


            var GsubIdx = new List<int>();
            var LsubIdx = new List<int>();
            var lowLocalBlocks_i0 = new List<int>();
            var lowLocalBlocks__N = new List<int>();
            int cnt = 0;
            for (int jLoc = 0; jLoc < J; jLoc++) {
                lowLocalBlocks_i0.Add(cnt);
                lowLocalBlocks__N.Add(0);

                //if (UseHiOrderSmoothing) {
                //    HighOrderBlocks_LU[jLoc] = new MultidimensionalArray[NoVars];
                //    HighLoOrderBlocks[jLoc] = new MultidimensionalArray[NoVars];
                //    HighOrderBlocks_LUpivots[jLoc] = new int[NoVars][];
                //}

                int NpHiTot = 0;
                for (int iVar = 0; iVar < NoVars; iVar++) {
                    int pReq;
                    if (iVar == 2) // quick hack for Anne Kikker
                        pReq = 0;
                    else
                        pReq = 1;

                    int Np1 = BS[iVar].GetLength(jLoc, pReq);
                    int Np = BS[iVar].GetLength(jLoc, degs[iVar]);
                    lowLocalBlocks__N[jLoc] += Np1;
                    NpHiTot += (Np - Np1);

                    for (int n = 0; n < Np1; n++) {

                        int Lidx = Map.LocalUniqueIndex(iVar, jLoc, n);
                        LsubIdx.Add(Lidx);

                        int Gidx = Map.GlobalUniqueIndex(iVar, jLoc, n);
                        GsubIdx.Add(Gidx);

                    }
                }

                if (UseHiOrderSmoothing && NpHiTot > 0) {
                    HighOrderBlocks_LU[jLoc] = MultidimensionalArray.Create(NpHiTot, NpHiTot);
                    int[] idxs = new int[NpHiTot];
                    HighOrderBlocks_indices[jLoc] = idxs;

                    int RowOffset = 0;
                    for (int iVar = 0; iVar < NoVars; iVar++) {
                        int pReq;
                        if (iVar == 2) // quick hack for Anne Kikker
                            pReq = 0;
                        else
                            pReq = 1;

                        int Np1 = BS[iVar].GetLength(jLoc, pReq);
                        int Np = BS[iVar].GetLength(jLoc, degs[iVar]);
                        int Nphi = Np - Np1;

                        int i0_hi = Map.GlobalUniqueIndex(iVar, jLoc, Np1);

                        for (int i = 0; i < Nphi; i++)
                            idxs[i + RowOffset] = i0_hi + i - Map.i0;

                        int ColOffset = 0;
                        for (int jVar = 0; jVar < NoVars; jVar++) {
                            int _pReq;
                            if (jVar == 2) // quick hack for Anne Kikker
                                _pReq = 0;
                            else
                                _pReq = 1;

                            int _Np1 = BS[jVar].GetLength(jLoc, _pReq);
                            int _Np = BS[jVar].GetLength(jLoc, degs[jVar]);
                            int j0_hi = Map.GlobalUniqueIndex(jVar, jLoc, _Np1);
                            int _Nphi = _Np - _Np1;
                            
                            m_op.OperatorMatrix.ReadBlock(i0_hi, j0_hi,
                                HighOrderBlocks_LU[jLoc].ExtractSubArrayShallow(new int[] { RowOffset, ColOffset }, new int[] { RowOffset + Nphi - 1, ColOffset + _Nphi - 1}));
                            
                            ColOffset += _Nphi;
                        }
                        Debug.Assert(ColOffset == NpHiTot);
                        RowOffset += (Np - Np1);
                    }
                    Debug.Assert(RowOffset == NpHiTot);

                    HighOrderBlocks_LUpivots[jLoc] = new int[NpHiTot];
                    HighOrderBlocks_LU[jLoc].FactorizeLU(HighOrderBlocks_LUpivots[jLoc]);
                    //HighOrderBlocks_LU[jLoc][iVar].Invert();
                }



                cnt += lowLocalBlocks__N[jLoc];
            }
            m_LsubIdx = LsubIdx.ToArray();


            BlockPartitioning localBlocking = new BlockPartitioning(GsubIdx.Count, lowLocalBlocks_i0, lowLocalBlocks__N, Map.MPI_Comm, i0isLocal:true);
            var P01SubMatrix = new BlockMsrMatrix(localBlocking);
            op.OperatorMatrix.AccSubMatrixTo(1.0, P01SubMatrix, GsubIdx, default(int[]), GsubIdx, default(int[]));

            intSolver = new PARDISOSolver() {
                CacheFactorization = true,
                UseDoublePrecision = false
            };
            intSolver.DefineMatrix(P01SubMatrix);
        }

        //BlockMsrMatrix P01SubMatrix;

        ISparseSolver intSolver;

        int[] m_LsubIdx;

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

            var Mtx = m_op.OperatorMatrix;

            // compute fine residual
            Res_f.SetV(B);
            Mtx.SpMV(-1.0, X, 1.0, Res_f);

            // project to low-p/coarse
            Res_c.AccV(1.0, Res_f, default(int[]), m_LsubIdx);

            // low-p solve
            intSolver.Solve(Cor_c, Res_c);

            // accumulate low-p correction
            X.AccV(1.0, Cor_c, m_LsubIdx, default(int[]));


            double[] Xbkup = X.ToArray();


            // solver high-order 
            if (UseHiOrderSmoothing) {
                // compute residual of low-order solution
                Res_f.SetV(B);
                Mtx.SpMV(-1.0, X, 1.0, Res_f);



                var Map = m_op.Mapping;
                int NoVars = Map.AggBasis.Length;
                int j0 = Map.FirstBlock;
                int J = Map.LocalNoOfBlocks;
                int[] degs = m_op.Degrees;
                var BS = Map.AggBasis;

                int Mapi0 = Map.i0;
                double[] b_f = null, x_hi = null;
                for (int j = 0; j < J; j++) {

                    if(HighOrderBlocks_LU[j] != null) {
                        int NpTotHi = HighOrderBlocks_LU[j].NoOfRows;
                        if (b_f == null || b_f.Length != NpTotHi) {
                            b_f = new double[NpTotHi];
                            x_hi = new double[NpTotHi];
                        }

                        ArrayTools.GetSubVector<int[],int[],double>(Res_f, b_f, HighOrderBlocks_indices[j]);
                        HighOrderBlocks_LU[j].BacksubsLU(HighOrderBlocks_LUpivots[j], x_hi, b_f);
                        X.AccV(1.0, x_hi, HighOrderBlocks_indices[j], default(int[]));
                    }


                    /*
                    for (int iVar = 0; iVar < NoVars; iVar++) {
                        int pReq = 1;

                        int Np1 = BS[iVar].GetLength(j, pReq);
                        int Np = BS[iVar].GetLength(j, degs[iVar]);
                        int i0_lo = Map.GlobalUniqueIndex(iVar, j, 0);
                        int i0_hi = Map.GlobalUniqueIndex(iVar, j, Np1);

                        double[] xLo = X.GetSubVector(i0_lo - Mapi0, Np1);
                        double[] bHi = Res_f.GetSubVector(i0_hi - Mapi0, Np - Np1);

                        double[] xhi = new double[Np - Np1];

                        //HighOrderBlocks_LU[j][iVar].Solve(xhi, bHi);
                        HighOrderBlocks_LU[j][iVar].BacksubsLU(HighOrderBlocks_LUpivots[j][iVar], xhi, bHi);
                        //HighOrderBlocks_LU[j][iVar].gemv(1.0, bHi, 0.0, xhi);

                        X.AccV(1.0, xhi, offset_acc: i0_hi - Mapi0);
                    }
                    */
                }
            }

            /*
            double[] Xf = X.ToArray();
            double[] exSol = new double[Xf.Length];
            m_op.OperatorMatrix.Solve_Direct(exSol, B);

            PlotVectors(new[] { Xbkup, Xf, exSol }, new[] { "lowP", "smooth", "exsol" });
            */

        }
    }
}

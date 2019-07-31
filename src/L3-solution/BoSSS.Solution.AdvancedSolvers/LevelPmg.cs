using BoSSS.Foundation;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
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

            MultigridOperator op_iLv = m_op;
            for (int i = 0; i < iLv; i++)
                op_iLv = op_iLv.CoarserLevel;


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

        public bool UseHiOrderSmoothing {
            get;
            set;
        }


        public ISolverSmootherTemplate Clone() {
            throw new NotImplementedException();
        }

        MultigridOperator m_op;

        MultidimensionalArray[][] HighOrderBlocks_LU;
        int[][][] HighOrderBlocks_LUpivots;
        MultidimensionalArray[][] HighLoOrderBlocks;

        


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
                HighOrderBlocks_LU = new MultidimensionalArray[J][];
                HighOrderBlocks_LUpivots = new int[J][][];
                HighLoOrderBlocks = new MultidimensionalArray[J][];
            }


            var GsubIdx = new List<int>();
            var LsubIdx = new List<int>();
            var lowLocalBlocks_i0 = new List<int>();
            var lowLocalBlocks__N = new List<int>();
            int cnt = 0;
            for (int jLoc = 0; jLoc < J; jLoc++) {
                lowLocalBlocks_i0.Add(cnt);
                lowLocalBlocks__N.Add(0);

                if (UseHiOrderSmoothing) {
                    HighOrderBlocks_LU[jLoc] = new MultidimensionalArray[NoVars];
                    HighLoOrderBlocks[jLoc] = new MultidimensionalArray[NoVars];
                    HighOrderBlocks_LUpivots[jLoc] = new int[NoVars][];
                }

                for (int iVar = 0; iVar < NoVars; iVar++) {
                    int pReq;
                    if (iVar == 2) // quick hack for Anne Kikker
                        pReq = 0;
                    else
                        pReq = 1;

                    int Np1 = BS[iVar].GetLength(jLoc, pReq);
                    int Np = BS[iVar].GetLength(jLoc, degs[iVar]);
                    lowLocalBlocks__N[jLoc] += Np1;

                    if (UseHiOrderSmoothing) {
                        HighOrderBlocks_LU[jLoc][iVar] = MultidimensionalArray.Create(Np - Np1, Np - Np1);
                        int i0_lo = Map.GlobalUniqueIndex(iVar, jLoc, 0);
                        int i0_hi = Map.GlobalUniqueIndex(iVar, jLoc, Np1);
                        m_op.OperatorMatrix.ReadBlock(i0_hi, i0_hi, HighOrderBlocks_LU[jLoc][iVar]);
                        HighOrderBlocks_LUpivots[jLoc][iVar] = new int[Np - Np1];

                        HighOrderBlocks_LU[jLoc][iVar].FactorizeLU(HighOrderBlocks_LUpivots[jLoc][iVar]);
                        //HighOrderBlocks_LU[jLoc][iVar].Invert();


                        HighLoOrderBlocks[jLoc][iVar] = MultidimensionalArray.Create(Np - Np1, Np1);
                        m_op.OperatorMatrix.ReadBlock(i0_lo, i0_hi, HighLoOrderBlocks[jLoc][iVar]);
                    }

                    for (int n = 0; n < Np1; n++) {

                        int Lidx = Map.LocalUniqueIndex(iVar, jLoc, n);
                        LsubIdx.Add(Lidx);

                        int Gidx = Map.GlobalUniqueIndex(iVar, jLoc, n);
                        GsubIdx.Add(Gidx);

                    }
                }

                cnt += lowLocalBlocks__N[jLoc];
            }
            m_LsubIdx = LsubIdx.ToArray();


            BlockPartitioning localBlocking = new BlockPartitioning(GsubIdx.Count, lowLocalBlocks_i0, lowLocalBlocks__N, Map.MPI_Comm);
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
                Res_f.SetV(B);
                Mtx.SpMV(-1.0, X, 1.0, Res_f);



                var Map = m_op.Mapping;
                int NoVars = Map.AggBasis.Length;
                int j0 = Map.FirstBlock;
                int J = Map.LocalNoOfBlocks;
                int[] degs = m_op.Degrees;
                var BS = Map.AggBasis;

                int Mapi0 = Map.i0;
                for (int j = 0; j < J; j++) {
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

                        X.AccV(1.0, xhi, offset_acc: i0_hi);
                    }
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

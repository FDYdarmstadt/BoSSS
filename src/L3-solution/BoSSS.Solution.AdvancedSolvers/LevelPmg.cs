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
    /// p-Multigrid on a single grid level
    /// </summary>
    public class LevelPmg : ISolverSmootherTemplate {

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

        public ISolverSmootherTemplate Clone() {
            throw new NotImplementedException();
        }

        MultigridOperator m_op;

        public void Init(MultigridOperator op) {
            //var Mtx = op.OperatorMatrix;
            m_op = op;

            var Map = op.Mapping;
            int NoVars = Map.AggBasis.Length;
            int j0 = Map.FirstBlock;
            int J = Map.LocalNoOfBlocks;

            var BS = Map.AggBasis;

            var GsubIdx = new List<int>();
            var LsubIdx = new List<int>();
            var lowLocalBlocks_i0 = new List<int>();
            var lowLocalBlocks__N = new List<int>();
            int cnt = 0;
            for (int jLoc = 0; jLoc < J; jLoc++) {
                lowLocalBlocks_i0.Add(cnt);
                lowLocalBlocks__N.Add(0);

                for (int iVar = 0; iVar < NoVars; iVar++) {
                    int Np = BS[iVar].GetLength(jLoc, 1);
                    lowLocalBlocks__N[jLoc] += Np;

                    for (int n = 0; n < Np; n++) {

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

        public void ResetStat() {
            Converged = false;
            ThisLevelIterations = 0;
        }


        double[] Res_f;

        double[] Res_c;
        double[] Cor_c;

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
        }
    }
}

using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// Incomplete LU decomposition on a cell-level
    /// 
    /// </summary>
    public class CellILU : ISolverSmootherTemplate {
        
        /// <summary>
        /// always 0, no sub-calls here
        /// </summary>
        public int IterationsInNested {
            get {
                return 0;
            }
        }

        int m_ThisLevelIterations;

        /// <summary>
        /// 
        /// </summary>
        public int ThisLevelIterations {
            get {
                return m_ThisLevelIterations;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public bool Converged {
            get;
            private set;
        }

        /// <summary>
        /// 
        /// </summary>
        public object Clone() {
            throw new NotImplementedException();
        }

        /// <summary>
        /// 
        /// </summary>
        public void Dispose() {
            m_BlockL = null;
            m_BlockU = null;
        }

        MultigridOperator m_op;

        /// <summary>
        /// 
        /// </summary>
        public void Init(MultigridOperator op) {
            m_op = op;
            UpdateILU();
        }

        /// <summary>
        /// 
        /// </summary>
        public void ResetStat() {
            m_ThisLevelIterations = 0;
            Converged = false;
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {


            int L = X.Count;
            double[] y = new double[L];

            m_BlockL.Solve_Direct(y, B);
            m_BlockU.Solve_Direct(X, y);

            
            m_ThisLevelIterations += 1;
            Converged = true;
        }

        BlockMsrMatrix m_BlockL;

        BlockMsrMatrix m_BlockU;

        /// <summary>
        /// 
        /// </summary>
        public long UsedMemory() {
            return (m_BlockL?.UsedMemory ?? 0L) + (m_BlockU?.UsedMemory ?? 0L);
        }


        void UpdateILU() {
            if(m_op.Mapping.MpiSize > 1)
                throw new NotSupportedException();

            m_BlockU = m_op.OperatorMatrix.CloneAs();
            m_BlockL = new BlockMsrMatrix(m_BlockU._RowPartitioning, m_BlockU._ColPartitioning);
            m_BlockL.AccEyeSp(1.0);
            var grd = m_op.Mapping.AggGrid;

            long j0 = m_BlockU._ColPartitioning.FirstBlock;
            long J = m_BlockU._ColPartitioning.LocalNoOfBlocks; 
            
            long i0 = m_BlockU._RowPartitioning.FirstBlock;
            long I = m_BlockU._RowPartitioning.LocalNoOfBlocks;

            if(j0 != i0)
                throw new NotSupportedException();
            if(J != I)
                throw new NotSupportedException();
            var part = m_BlockU._RowPartitioning;

            MultidimensionalArray invUii = null;
            for(long i = j0; i < J-1; i++) { // Iteration over matrix (rows and columns)

                long iRow = part.GetBlockI0(i);
                int sz_i = part.GetBlockLen(i);
                if(sz_i <= 0)
                    continue; // cell i is empty


                if(invUii == null || invUii.NoOfCols != sz_i)
                    invUii = MultidimensionalArray.Create(sz_i, sz_i);
                m_BlockU.ReadBlock(iRow, iRow, invUii);
                invUii.InvertInPlace();

                // loop over remaining rows ....
                // For k = (i+1):N, and if (k,i) in P
                var Neighs_of_i = grd.GetCellNeighboursViaEdges(checked((int)(i - j0)));
                foreach(var tt in Neighs_of_i) {
                    long k = tt.jCellLoc + i0;
                    if(k <= i) // k > i
                        continue;
                    Debug.Assert(k > i);


                    // L(Bk,Bi) = U(Bk,Bi)*invUii
                    long kRow = part.GetBlockI0(k);
                    int sz_k = part.GetBlockLen(k);
                    if(sz_k <= 0)
                        continue; // cell k & block-row k are empty

                    MultidimensionalArray U_ki = MultidimensionalArray.Create(sz_k, sz_i);
                    m_BlockU.ReadBlock(kRow, iRow, U_ki);
                    MultidimensionalArray L_ki = U_ki.GEMM(invUii);

                    m_BlockL.ClearBlock(kRow, iRow, sz_k, sz_i);
                    m_BlockL.AccBlock(kRow, iRow, 1.0, L_ki);
                    
                    // loop over remaining columns...
                    foreach(var ttt in Neighs_of_i) {
                        long j = ttt.jCellLoc + i0;
                        if(j <= i)
                            continue;
                        Debug.Assert(j > i);

                        long jRow = part.GetBlockI0(j);
                        int sz_j = part.GetBlockLen(j);
                        if(sz_j <= 0)
                            continue; // cell j & block-col j are empty

                        MultidimensionalArray U_ij = MultidimensionalArray.Create(sz_i, sz_j);
                        m_BlockU.ReadBlock(iRow, jRow, U_ij);

                        MultidimensionalArray L_ki__x__U_ij = L_ki.GEMM(U_ij);
                        m_BlockU.AccBlock(kRow, jRow, -1.0, L_ki__x__U_ij);
                    }
                }
            }
        }
    }
}

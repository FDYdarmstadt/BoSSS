using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
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
            CheckILU();
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

        /// <summary>
        /// Full LU, out-of-place version, seems to be ok.
        /// </summary>
        void UpdateILU_Full() {
            if(m_op.Mapping.MpiSize > 1)
                throw new NotImplementedException();

            var grd = m_op.Mapping.GridData;
            var Mtx = m_op.OperatorMatrix;
            IBlockPartitioning part = m_op.Mapping;
            Debug.Assert(part.EqualsPartition(Mtx._RowPartitioning), "mismatch in row partition");
            Debug.Assert(part.EqualsPartition(Mtx._ColPartitioning), "mismatch in column partition");
            long cell0 = part.FirstBlock;
            int J = part.LocalNoOfBlocks;


            m_BlockU = m_op.OperatorMatrix.CloneAs();
            m_BlockL = new BlockMsrMatrix(part, part);
            m_BlockL.AccEyeSp(1.0);
           

            bool[,] P1 = GetPattern();

            MultidimensionalArray invUii = null;
            for(long i = cell0; i < (cell0 + J - 1); i++) { // Iteration over matrix (rows and columns)

                long iIdx = part.GetBlockI0(i);
                int sz_i = part.GetBlockLen(i);
                if(sz_i <= 0)
                    continue; // cell i is empty


                if(invUii == null || invUii.NoOfCols != sz_i)
                    invUii = MultidimensionalArray.Create(sz_i, sz_i);
                m_BlockU.ReadBlock(iIdx, iIdx, invUii);
                invUii.InvertInPlace();

                // loop over remaining rows ....
                // For k = (i+1):N, and if (k,i) in P
                var Neighs_of_i = grd.GetCellNeighboursViaEdges(checked((int)(i - cell0))).Select(ttt => ttt.jCellLoc).ToArray();
                //foreach(int kLoc in Neighs_of_i) {
                //    long k = kLoc + i0;
                //    if(k <= i) // k > i
                //        continue;

                for(long k = i + 1; k < (cell0 + J); k++) {
                    Debug.Assert(k > i);


                    // L(Bk,Bi) = U(Bk,Bi)*invUii
                    long kIdx = part.GetBlockI0(k);
                    int sz_k = part.GetBlockLen(k);
                    if(sz_k <= 0)
                        continue; // cell k & block-row k are empty

                    MultidimensionalArray U_ki = MultidimensionalArray.Create(sz_k, sz_i);
                    m_BlockU.ReadBlock(kIdx, iIdx, U_ki);
                    MultidimensionalArray L_ki = U_ki.GEMM(invUii);

                    m_BlockL.ClearBlock(kIdx, iIdx, sz_k, sz_i);
                    m_BlockL.AccBlock(kIdx, iIdx, 1.0, L_ki);

                    // loop over remaining columns...
                    // for j = i:N % loop over remaining columns
                    var Neighs_of_k = grd.GetCellNeighboursViaEdges(checked((int)(k - cell0))).Select(ttt => ttt.jCellLoc);
                    var colPattern = ArrayTools.Cat(Neighs_of_k, checked((int)(k - cell0)));
                    //foreach(int jLoc in colPattern) {
                    //    long j = jLoc + i0;
                    //    if(j < i)
                    //        continue;
                    for(long j = i; j < (cell0 + J); j++) {
                        Debug.Assert(j >= i);

                        long jIdx = part.GetBlockI0(j);
                        int sz_j = part.GetBlockLen(j);
                        if(sz_j <= 0)
                            continue; // cell j & block-col j are empty

                        MultidimensionalArray U_ij = MultidimensionalArray.Create(sz_i, sz_j);
                        m_BlockU.ReadBlock(iIdx, jIdx, U_ij);

                        MultidimensionalArray L_ki__x__U_ij = L_ki.GEMM(U_ij);
                        m_BlockU.AccBlock(kIdx, jIdx, -1.0, L_ki__x__U_ij);


                        if(k > j && j == i) {
                            MultidimensionalArray U_kj = MultidimensionalArray.Create(sz_k, sz_j);
                            m_BlockU.ReadBlock(kIdx, jIdx, U_kj);
                            Debug.Assert(U_kj.L2Norm() <= 1.0e-4);
                            //Console.WriteLine($"U[{k},{j}] norm : {U_kj.L2Norm()}");
                        }
                    }
                }
            }

            m_BlockL.SaveToTextFileSparse("L.txt");
            m_BlockU.SaveToTextFileSparse("U.txt");
        }




        /// <summary>
        /// ILU, in-place version, after book of Saad (p 303)
        /// </summary>
        void UpdateILU() {
            if(m_op.Mapping.MpiSize > 1)
                throw new NotImplementedException();

            var grd = m_op.Mapping.GridData;
            var Mtx = m_op.OperatorMatrix;
            IBlockPartitioning part = m_op.Mapping;
            Debug.Assert(part.EqualsPartition(Mtx._RowPartitioning), "mismatch in row partition");
            Debug.Assert(part.EqualsPartition(Mtx._ColPartitioning), "mismatch in column partition");
            long cell0 = part.FirstBlock;
            int J = part.LocalNoOfBlocks;


            


            MultidimensionalArray GetBlock(BlockMsrMatrix M, long iBlk, long jBlk) {
                int _sz_i = part.GetBlockLen(iBlk);
                int _sz_j = part.GetBlockLen(jBlk);
                long _idx_i = part.GetBlockI0(iBlk);
                long _idx_j = part.GetBlockI0(jBlk);
                var ret = MultidimensionalArray.Create(_sz_i, _sz_j);
                M.ReadBlock(_idx_i, _idx_j, ret);
                return ret;
            }

            void SetBlock(BlockMsrMatrix M, MultidimensionalArray Blk, long iBlk, long jBlk) {
                int _sz_i = part.GetBlockLen(iBlk);
                int _sz_j = part.GetBlockLen(jBlk);
                if(_sz_i != Blk.NoOfRows)
                    throw new ArgumentException();
                if(_sz_j != Blk.NoOfCols)
                    throw new ArgumentException();
                long _idx_i = part.GetBlockI0(iBlk);
                long _idx_j = part.GetBlockI0(jBlk);

                M.AccBlock(_idx_i, _idx_j, 1.0, Blk, 0.0);
            }

            var A = Mtx.CloneAs();
           

            bool[,] P1 = GetPattern();

            for(long k = cell0; k < (cell0 + J - 1); k++) { // Iteration over matrix (rows and columns)
                int sz_k = part.GetBlockLen(k);
                if(sz_k <= 0)
                    continue; // cell i is empty

                var invAkk = GetBlock(A, k, k);
                invAkk.InvertInPlace();

                for(long i = k + 1; i < (cell0 + J); i++) {
                    Debug.Assert(i > k);
                    if(P1[i, k]) {


                        var Aik = GetBlock(A, i, k);
                        Aik = Aik.GEMM(invAkk);
                        SetBlock(A, Aik, i, k);

                        for(long j = k + 1; j < (cell0 + J); j++) {
                            if(P1[i, j]) {
                                var Aij = GetBlock(A, i, j);
                                var Akj = GetBlock(A, k, j);
                                Aij.GEMM(-1, Aik, Akj, 1.0);
                                SetBlock(A, Aij, i, j);
                            }
                        }

                        // L(Bk,Bi) = U(Bk,Bi)*invUii
                    }
                }
            }


            m_BlockU = new BlockMsrMatrix(part, part);
            m_BlockL = new BlockMsrMatrix(part, part);
            m_BlockL.AccEyeSp(1.0);

            for(long j = cell0; j < J + cell0; j++) {
                for(long i = cell0; i < J + cell0; i++) {
                    if(P1[i,j]) {
                        var Aji = GetBlock(A, j, i);
                        if(j > i) {
                            // lower tri
                            SetBlock(m_BlockL, Aji, j, i);
                        } else {
                            // upper tri
                            SetBlock(m_BlockU, Aji, j, i);
                        }


                    } else {
                        var Aji = GetBlock(A, j, i);
                        double ErrNorm = Aji.L2Norm();
                        if(ErrNorm > 0)
                            throw new ArithmeticException();
                    }
                }
            }

            m_BlockL.SaveToTextFileSparse("L.txt");
            m_BlockU.SaveToTextFileSparse("U.txt");
        }

        
        
        private bool[,] GetPattern() {
            if(m_op.Mapping.MpiSize > 1)
                throw new NotImplementedException();

            var grd = m_op.Mapping.GridData;
            var Mtx = m_op.OperatorMatrix;
            IBlockPartitioning part = m_op.Mapping;
            Debug.Assert(part.EqualsPartition(Mtx._RowPartitioning), "mismatch in row partition");
            Debug.Assert(part.EqualsPartition(Mtx._ColPartitioning), "mismatch in column partition");
            long cell0 = part.FirstBlock;
            int J = part.LocalNoOfBlocks;

            bool[,] P1 = new bool[J, J];
            bool[,] P2 = new bool[J, J];
            for(int j = 0; j < J; j++) {

                var Neighs_of_j = grd.GetCellNeighboursViaEdges(checked((int)(j - cell0))).Select(ttt => ttt.jCellLoc).ToArray();
                var colPattern = ArrayTools.Cat(Neighs_of_j, checked((int)(j - cell0)));
                foreach(int i in colPattern)
                    P2[j, i] = true;

                for(int i = 0; i < J; i++) {
                    long iB = part.GetBlockI0(i);
                    long jB = part.GetBlockI0(j);
                    int sz_i = part.GetBlockLen(i);
                    int sz_j = part.GetBlockLen(j);

                    var Mji = MultidimensionalArray.Create(sz_i, sz_j);
                    m_op.OperatorMatrix.ReadBlock(jB, iB, Mji);
                    if(Mji.L2Norm() > 0)
                        P1[j, i] = true;
                }
            }

            for(int j = 0; j < J; j++) {
                for(int i = 0; i < J; i++) {
                    Debug.Assert(P1[i, j] == P2[i, j], "cell neighborship differs from occupancy pattern");
                    Debug.Assert(P1[i, j] == P1[j, i], "missing structural symmetry");
                }
            }

            return P1;
        }

        void CheckILU() {
            BlockMsrMatrix LxU = BlockMsrMatrix.Multiply(m_BlockL, m_BlockU);
            BlockMsrMatrix Lc = m_BlockL.CloneAs();
            BlockMsrMatrix Uc = m_BlockU.CloneAs();


            bool[,] P1 = GetPattern();


            var grd = m_op.Mapping.GridData;

            BlockMsrMatrix ErrMtx = LxU.CloneAs();
            ErrMtx.Acc(-1, m_op.OperatorMatrix);
            double totErr = ErrMtx.InfNorm();
            Console.WriteLine("|ILU - Mtx| : " + totErr);

            var part = LxU._RowPartitioning;
            var i0 = part.FirstBlock;
            var iE = part.LocalNoOfBlocks + i0;

            double Errsum = 0;
            double UloErr = 0;
            double LupErr = 0;
            for(long j = i0; j < iE; j++) {
                long jRow = part.GetBlockI0(j);
                int sz_j = part.GetBlockLen(j);
                if(sz_j <= 0)
                    continue;

                //var Neighs_of_j = grd.GetCellNeighboursViaEdges(checked((int)(j - i0))).Select(tt => tt.jCellLoc);
                //var rowPattern = ArrayTools.Cat(Neighs_of_j, checked((int)(j - i0)));
                //foreach(int kLoc in rowPattern) {
                //    long k = kLoc + i0;
                for(long k = i0; k < iE; k++) {
                    
                    long kRow = part.GetBlockI0(k);
                    int sz_k = part.GetBlockLen(k);
                    if(sz_k <= 0)
                        continue;

                    if(P1[j, k]) {
                        var A = MultidimensionalArray.Create(sz_j, sz_k);
                        var B = MultidimensionalArray.Create(sz_j, sz_k);

                        LxU.ReadBlock(jRow, kRow, A);
                        m_op.OperatorMatrix.ReadBlock(jRow, kRow, B);

                        double dist_kj = A.L2Dist(B);
                        //Console.WriteLine($"dist[{k}, {j}] = {dist_kj}");

                        Errsum += dist_kj.Pow2();
                    }


                    if(j > k) {
                        // lower triangle: U must be 0
                        var Ujk = MultidimensionalArray.Create(sz_j, sz_k);
                        m_BlockU.ReadBlock(jRow, kRow, Ujk);
                        double shouldBe0 = Ujk.L2Norm(); 
                        if(shouldBe0 > 1.0e-4)
                            Console.WriteLine($"U[{j}, {k}] = {shouldBe0}");
                        UloErr += shouldBe0.Pow2();
                    } else {
                        // clear upper part + diagonal => the remainder of Uc must be 0
                        Uc.ClearBlock(jRow, kRow, sz_j, sz_k);
                    }

                    if(j < k) {
                        // upper triangle: L must be 0
                        var Ljk = MultidimensionalArray.Create(sz_j, sz_k);
                        m_BlockL.ReadBlock(jRow, kRow, Ljk);
                        LupErr += Ljk.L2Norm().Pow2();
                    } else {
                        // clear lower part + diagonal => the remainder of Lc must be 0
                        Lc.ClearBlock(jRow, kRow, sz_j, sz_k);
                    }




                }
            }
            Errsum = Math.Sqrt(Errsum);

            double Uerr = Uc.InfNorm();
            double Lerr = Lc.InfNorm();

            Console.WriteLine("Error in occupied: " + Errsum);
            Console.WriteLine("Low Error in U: " + Uerr);
            Console.WriteLine("Upp error in L: " + Lerr);

            System.Environment.Exit(-1);
        }
    
    }
}

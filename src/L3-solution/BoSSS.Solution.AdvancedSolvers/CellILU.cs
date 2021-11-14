using BoSSS.Foundation.Grid;
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
            m_ILUp_pattern = null;
        }

        MultigridOperator m_op;

        /// <summary>
        /// 
        /// </summary>
        public void Init(MultigridOperator op) {
            m_op = op;
            Console.WriteLine($"CellILU, MG level {op.LevelIndex}...");
            UpdateILU();
            //CheckILU();
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
            where V : IList<double>  //
        {


            int L = X.Count;

            //double[] _B;
            //if(B.GetType() == typeof(double[]))
            //    _B = B as double[];
            //else
            //    _B = B.ToArray();


            double[] y = new double[L];
            //double[] yref = new double[L];
            //m_BlockL.Solve_Direct(yref, B.ToArray().CloneAs());
            LoTriDiagonalSolve(m_BlockL, y, B, true);
            //double check1 = GenericBlas.L2Dist(yref, y);
            //Console.WriteLine("Check value (lo solve) is: " + check1);


            //double[] Xref = new double[y.Length];
            //m_BlockU.Solve_Direct(Xref, y.CloneAs());
            UpTriDiagonalSolve(m_BlockU, X, y, false);
            //double check2 = GenericBlas.L2Dist(Xref, X);
            //Console.WriteLine("Check value (hi solve) is: " + check2);

            
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

        /*
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
           

            bool[,] P1 = GetPattern(out _);

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
        */

        /// <summary>
        /// ILU, in-place version, after book of Saad (p 303)
        /// </summary>
        void UpdateILU() {
            using(var ft = new FuncTrace()) {
                if(m_op.Mapping.MpiSize > 1)
                    throw new NotImplementedException();

                var grd = m_op.Mapping.GridData;
                var Mtx = m_op.OperatorMatrix;
                IBlockPartitioning part = m_op.Mapping;
                Debug.Assert(part.EqualsPartition(Mtx._RowPartitioning), "mismatch in row partition");
                Debug.Assert(part.EqualsPartition(Mtx._ColPartitioning), "mismatch in column partition");
                long cell0 = part.FirstBlock;
                int J = part.LocalNoOfBlocks;


                var A = Mtx.CloneAs();


                var ILUp = GetPattern(out int Occupied);
                var ILUpT = m_ILUp_pattern.Transpose();
                double occupancy = (double)Occupied / (double)(J * J);
                ft.Info($"CellILU, lv {m_op.LevelIndex}, ILU-{ILU_level}, occupancy = {Math.Round(occupancy * 100)}%.");

                for(long k = cell0; k < (cell0 + J - 1); k++) { // Iteration over matrix (rows and columns)
                    int sz_k = part.GetBlockLen(k);
                    if(sz_k <= 0)
                        continue; // cell i is empty

                    var invAkk = A.GetBlock(k, k);
                    invAkk.InvertInPlace();

                    long[] occColumn_k = ILUpT.GetOccupiedColumnIndices(k);
                    foreach(long i in occColumn_k) { 
                    //for(long i = k + 1; i < (cell0 + J); i++) {
                    //    Debug.Assert(i > k);
                    //    if(P1[i, k]) {
                        if(i >= k + 1) { 
                            //if(!occColumn_k.Contains(i))
                            //    throw new Exception($"fuck: P1[{i},{k}] = {P1[i, k]} //  {ILUpT[i, k]}  {m_ILUp_pattern[i, k]} ");

                            var Aik = A.GetBlock(i, k);
                            Aik = Aik.GEMM(invAkk);
                            A.SetBlock(Aik, i, k);

                            long[] occRow_i = m_ILUp_pattern.GetOccupiedColumnIndices(i);
                            //for(long j = k + 1; j < (cell0 + J); j++) {
                            //    if(P1[i, j] != occRow_i.Contains(j))
                            //        throw new Exception("fuck 2");
                            //    if(P1[i, j]) {
                            foreach(long j in occRow_i) { 
                                if(j >= k + 1) { 
                                    var Aij = A.GetBlock(i, j);
                                    var Akj = A.GetBlock(k, j);
                                    Aij.GEMM(-1, Aik, Akj, 1.0);
                                    A.SetBlock(Aij, i, j);
                                }
                            }

                            // L(Bk,Bi) = U(Bk,Bi)*invUii
                        } else {
                            //if(occColumn_k.Contains(i))
                            //    throw new Exception($"fuck: P1[{i},{k}] = {P1[i, k]} //  {ILUpT[i, k]}  {m_ILUp_pattern[i, k]} ");
                        }
                    }
                }


                m_BlockU = new BlockMsrMatrix(part, part);
                m_BlockL = new BlockMsrMatrix(part, part);
                m_BlockL.AccEyeSp(1.0);

                for(long j = cell0; j < J + cell0; j++) {
                    //for(long i = cell0; i < J + cell0; i++) {
                    //    if(P1[i, j]) {
                    long[] occRow_j = ILUp.GetOccupiedColumnIndices(j);
                    foreach(long i in occRow_j) {
                        {
                            var Aji = A.GetBlock(j, i);
                            if(j > i) {
                                // lower tri
                                m_BlockL.SetBlock(Aji, j, i);
                            } else {
                                // upper tri
                                m_BlockU.SetBlock(Aji, j, i);
                            }
                        }
                        //} else {
                        //    var Aji = A.GetBlock(j, i);
                        //    double ErrNorm = Aji.L2Norm();
                        //    if(ErrNorm > 0)
                        //        throw new ArithmeticException();
                        //}
                    }
                }

                //m_BlockL.SaveToTextFileSparse("L.txt");
                //m_BlockU.SaveToTextFileSparse("U.txt");
            }
        }

        /*
        /// <summary>
        /// Very slow reference version, probably does not work in the presence of agglomeration
        /// </summary>
        private bool[,] GetPattern_ref() {
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
            for(int i = 0; i < J; i++) { // loop over block rows

                
                //var Neighs_of_j = grd.GetCellNeighboursViaEdges(checked((int)(j - cell0))).Select(ttt => ttt.jCellLoc).ToArray();
                //var colPattern = ArrayTools.Cat(Neighs_of_j, checked((int)(j - cell0)));
                //foreach(int i in colPattern)
                //    P2[j, i] = true;
                

                long[] colBlocks = Mtx.GetOccupiedRowBlockIndices(i);
                foreach(long j in colBlocks) { // loop over block columns
                    var Blk = Mtx.GetBlock(j, i);
                    if(Blk.L2Norm() > 0)
                        P2[i, j] = true;
                    else
                        Console.WriteLine("Mem occupied, but the block is 0");
                }


                for(int j = 0; j < J; j++) { // loop over block columns
                    long iB = part.GetBlockI0(j);
                    long jB = part.GetBlockI0(i);
                    int sz_i = part.GetBlockLen(j);
                    int sz_j = part.GetBlockLen(i);

                    var Mji = MultidimensionalArray.Create(sz_i, sz_j);
                    m_op.OperatorMatrix.ReadBlock(jB, iB, Mji);
                    if(Mji.L2Norm() > 0)
                        P1[i, j] = true;
                    
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
        */


        public int ILU_level = 1;

        /// <summary>
        /// Optimized version
        /// </summary>
        private MsrMatrix GetPattern(out int Occupied) {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = true;
                if(m_op.Mapping.MpiSize > 1)
                    throw new NotImplementedException();
                if(ILU_level < 0)
                    throw new ArgumentException("ILU-Level must be >= 0, got " + ILU_level);
                if(ILU_level > Math.Min((int)(byte.MaxValue), 10))
                    throw new ArgumentException("ILU-Level must be <= 10, got " + ILU_level);

                var grd = m_op.Mapping.GridData;
                var Mtx = m_op.OperatorMatrix;
                IBlockPartitioning part = m_op.Mapping;
                Debug.Assert(part.EqualsPartition(Mtx._RowPartitioning), "mismatch in row partition");
                Debug.Assert(part.EqualsPartition(Mtx._ColPartitioning), "mismatch in column partition");
                long cell0 = part.FirstBlock;
                int J = part.LocalNoOfBlocks;

                byte[,] lev = new byte[J, J];
                lev.SetAll(byte.MaxValue);

                // determinte ILU-0 pattern
                tr.Info("computing ILU(0) pattern...");
                var p = new Partitioning(J, csMPI.Raw._COMM.SELF);
                var ILU0_pattern = new MsrMatrix(p, p);
                for(long i = cell0; i < (J + cell0); i++) { // loop over block rows
                    //long iB = part.GetBlockI0(i);
                    //int sz_i = part.GetBlockLen(i);

                    long[] colBlocks = Mtx.GetOccupiedRowBlockIndices(i);
                    foreach(long j in colBlocks) { // loop over block columns
                        var Blk = Mtx.GetBlock(j, i);
                        if(Blk.L2Norm() > 0) {
                            lev[i, j] = 0;
                            ILU0_pattern[i, j] = 1;
                        } else {
                            tr.Info("Mem occupied, but the block is 0");
                        }
                    }
                }
                //ILU0pattern.SaveToTextFileSparse("ILU0.txt");
                Console.WriteLine("done.");
                MsrMatrix ILUp_pattern = ILU0_pattern;
                for(int iLevel = 1; iLevel <= ILU_level; iLevel++) {
                    tr.Info($"computing ILU({iLevel}) pattern...");

                    ILUp_pattern = MsrMatrix.Multiply(ILUp_pattern, ILU0_pattern);

                    /*
                    for(int j = 0; j < J; j++) {
                        for(int i = 0; i < J; i++) {
                            byte lev_ij = lev[i, j];

                            for(int k = 0; k < J; k++) {
                                lev_ij = checked((byte)Math.Min(lev_ij, (int)lev[i, k] + (int)lev[k, j] + 1));
                            }

                            lev[i, j] = lev_ij;
                        }
                    }
                    */
                    tr.Info("done.");
                }

                Occupied = 0;
                bool[,] ret = new bool[J, J];
                for(int j = 0; j < J; j++) {

                    long[] rowPattern = ILUp_pattern.GetOccupiedColumnIndices(j);
                    /*
                    var checkPattern = new List<long>();

                    for(int i = 0; i < J; i++) {
                        if(lev[i, j] != lev[j, i])
                            throw new ArithmeticException("missing structural symmetry of operator matrix");

                        if(lev[i, j] <= ILU_level) {
                            Occupied++;
                            ret[i, j] = true;
                            checkPattern.Add(i);
                        }
                    }

                    if(!rowPattern.SetEquals(checkPattern)) {
                        throw new Exception();
                    }
                    */
                }

                m_ILUp_pattern = ILUp_pattern;
                return m_ILUp_pattern;
            }
        }

        MsrMatrix m_ILUp_pattern;

        void CheckILU() {
            using(new FuncTrace()) {
                BlockMsrMatrix LxU = BlockMsrMatrix.Multiply(m_BlockL, m_BlockU);
                BlockMsrMatrix Lc = m_BlockL.CloneAs();
                BlockMsrMatrix Uc = m_BlockU.CloneAs();


                var P1 = GetPattern(out _);


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

                        if(P1[j, k] != 0) {
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

                //System.Environment.Exit(-1);
            }
        }

        void LoTriDiagonalSolve<U, V>(BlockMsrMatrix Mtx, U X, V B, bool diagEye)
            where U : IList<double>
            where V : IList<double>  //
        {
            using(new FuncTrace()) {
                var rowPart = Mtx._RowPartitioning;
                var colPart = Mtx._ColPartitioning;
                long i0 = rowPart.FirstBlock;
                long J = rowPart.LocalNoOfBlocks;
                long i0Idx = rowPart.i0;
                long j0Idx = colPart.i0;


                for(long i = i0; i < i0 + J; i++) {
                    int sz_i = rowPart.GetBlockLen(i);
                    if(sz_i <= 0)
                        continue; // empty cell;

                    if(colPart.GetBlockLen(i) <= 0)
                        throw new NotSupportedException("cannot do (lower) trigonal solve on non-quadratic zero diagonal block");

                    long iIdx = rowPart.GetBlockI0(i);
                    int iIdxLoc = checked((int)(iIdx - i0Idx));

                    double[] Bi = B.GetSubVector(iIdxLoc, sz_i);

                    long[] jS = Mtx.GetOccupiedRowBlockIndices(i);
                    foreach(long j in jS) {
                        if(j >= i) // ignore everything above the diagonal
                            continue;
                        int sz_j = colPart.GetBlockLen(j);
                        if(sz_j <= 0)
                            continue;

                        long jIdx = colPart.GetBlockI0(j);
                        int jIdxLoc = checked((int)(jIdx - j0Idx));

                        double[] Xj = X.GetSubVector(jIdxLoc, sz_j);
                        var Mtx_ij = Mtx.GetBlock(i, j);

                        Mtx_ij.GEMV(-1.0, Xj, 1.0, Bi);
                    }

                    if(diagEye) {
                        X.SetSubVector<double, U, double[]>(Bi, iIdxLoc, sz_i);
                    } else {
                        double[] Xi = new double[sz_i];
                        var Mtx_ii = Mtx.GetBlock(i, i);
                        Mtx_ii.Solve(Xi, Bi);
                        X.SetSubVector<double, U, double[]>(Xi, iIdxLoc, sz_i);
                    }
                }
            }
        }


        void UpTriDiagonalSolve<U, V>(BlockMsrMatrix Mtx, U X, V B, bool diagEye)
            where U : IList<double>
            where V : IList<double>  //
        {
            using(new FuncTrace()) {
                var rowPart = Mtx._RowPartitioning;
                var colPart = Mtx._ColPartitioning;
                long i0 = rowPart.FirstBlock;
                long J = rowPart.LocalNoOfBlocks;
                long i0Idx = rowPart.i0;
                long j0Idx = colPart.i0;


                for(long i = i0 + J - 1; i >= i0; i--) {
                    int sz_i = rowPart.GetBlockLen(i);
                    if(sz_i <= 0)
                        continue; // empty cell;

                    if(colPart.GetBlockLen(i) <= 0)
                        throw new NotSupportedException("cannot do (upper) trigonal solve on non-quadratic zero diagonal block");

                    long iIdx = rowPart.GetBlockI0(i);
                    int iIdxLoc = checked((int)(iIdx - i0Idx));

                    double[] Bi = B.GetSubVector(iIdxLoc, sz_i);

                    long[] jS = Mtx.GetOccupiedRowBlockIndices(i);
                    foreach(long j in jS) {
                        if(j <= i) // ignore everything below the diagonal
                            continue;
                        int sz_j = colPart.GetBlockLen(j);
                        if(sz_j <= 0)
                            continue;

                        long jIdx = colPart.GetBlockI0(j);
                        int jIdxLoc = checked((int)(jIdx - j0Idx));

                        double[] Xj = X.GetSubVector(jIdxLoc, sz_j);
                        var Mtx_ij = Mtx.GetBlock(i, j);

                        Mtx_ij.GEMV(-1.0, Xj, 1.0, Bi);
                    }

                    if(diagEye) {
                        X.SetSubVector<double, U, double[]>(Bi, iIdxLoc, sz_i);
                    } else {
                        double[] Xi = new double[sz_i];
                        var Mtx_ii = Mtx.GetBlock(i, i);
                        Mtx_ii.Solve(Xi, Bi);
                        X.SetSubVector<double, U, double[]>(Xi, iIdxLoc, sz_i);
                    }
                }
            }
        }
    }
}

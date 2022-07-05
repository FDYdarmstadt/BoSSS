using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.Connectors.Matlab;
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
    public class CellILU : ISubsystemSolver {
        
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

        IOperatorMappingPair m_op;

        /// <summary>
        /// 
        /// </summary>
        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        /// <summary>
        /// 
        /// </summary>
        public void Init(MultigridOperator op) {
            InitImpl(op);
        }

        public string id;

        /// <summary>
        /// 
        /// </summary>
        public void InitImpl(IOperatorMappingPair op) {
            using(new FuncTrace()) {
                if(object.ReferenceEquals(op, m_op))
                    return; // already initialized
                else
                    this.Dispose();

                m_op = op;

                (this.m_BlockL, this.m_BlockL_compressed, this.m_BlockL_OccupiedBlockColumnsPerRow,
                    this.m_BlockU, this.m_BlockU_OccupiedBlockColumnsPerRow,
                    this.m_BlockU_lublks, this.m_BlockU_luipiv) = UpdateILU(this.ILU_level, op.OperatorMatrix);
                //CheckILU();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public void ResetStat() {
            m_ThisLevelIterations = 0;
            Converged = false;
        }

        bool written = true;


        public static void Verify(string _iD) {
            var Blocks = VectorIO.LoadFromTextFile("part" + _iD + ".txt").Select(d => (int)d).ToArray();
            int J = Blocks.Length;
            long[] i0 = new long[J];
            for (int j = 1; j < J; j++) {
                i0[j] = i0[j - 1] + Blocks[j - 1];
            }
            var part = new BlockPartitioning(checked((int)(i0[J - 1] + Blocks[J - 1])), i0, Blocks, csMPI.Raw._COMM.SELF);


            var opMtx_tmp = MsrMatrix.LoadFromFile("M" + _iD + ".mtx", csMPI.Raw._COMM.SELF, part, part);
            var L_tmp = MsrMatrix.LoadFromFile("L" + _iD + ".mtx", csMPI.Raw._COMM.SELF, part, part);
            var U_tmp = MsrMatrix.LoadFromFile("U" + _iD + ".mtx", csMPI.Raw._COMM.SELF, part, part);

            var opMtx = new BlockMsrMatrix(part); opMtx.Acc(1.0, opMtx_tmp);
            var L = new BlockMsrMatrix(part); L.Acc(1.0, L_tmp);
            var U = new BlockMsrMatrix(part); U.Acc(1.0, U_tmp);


            for (int i = 0; i < 1; i++) {
                Console.WriteLine("pass " + i + " on " + _iD +  " ...");
                var (test_L1, _, _, test_U1, _, _, _) = UpdateILU(0, opMtx);
                var (test_L2, _, _, test_U2, _, _, _) = UpdateILU(0, opMtx);

                //double cond_Matlab = opMtx.condest();
                //double cond_Arnold = opMtx.condestArnoldi();

                double LErr1 = L.MatrixDist(test_L1);
                double UErr1 = U.MatrixDist(test_U1);
                double LErr2 = L.MatrixDist(test_L2);
                double UErr2 = U.MatrixDist(test_U2);
                double LErr12 = test_L1.MatrixDist(test_L2);
                double UErr12 = test_U1.MatrixDist(test_U2);

                var LDiff = test_L1.Minus(test_L2);
                var UDiff = test_U1.Minus(test_U2);


                Console.WriteLine($" {LErr1}, {UErr1}");
                Console.WriteLine($" {LErr2}, {UErr2}");
                Console.WriteLine($" {LErr12}, {UErr12}");
            }

            //Console.WriteLine($"ID: {_iD}: LErr = {LErr}, UErr = {UErr}; cond = {cond_Matlab} (Matlab), {cond_Arnold} (BoSSS Arnoldi)");

        }


        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double>  //
        {
            if (!written) {
                var part = m_op.OperatorMatrix._RowPartitioning;
                if (part.MpiSize != 1)
                    throw new NotSupportedException();
                if (part.LocalNoOfBlocks != part.TotalNoOfBlocks)
                    throw new NotSupportedException();

                int NoOfBlocks = part.LocalNoOfBlocks;
                int[] BL = NoOfBlocks.ForLoop(iblk => part.GetBlockLen(iblk));
                BL.SaveToTextFile("part" + id + ".txt", part.MPI_Comm);

                m_op.OperatorMatrix.SaveToTextFileSparse("M" + id + ".txt");
                m_BlockL.SaveToTextFileSparse("L" + id + ".txt");
                m_BlockU.SaveToTextFileSparse("U" + id + ".txt");

                m_op.OperatorMatrix.ToMsrMatrix().SaveToFile("M" + id + ".mtx");
                m_BlockL.ToMsrMatrix().SaveToFile("L" + id + ".mtx");
                m_BlockU.ToMsrMatrix().SaveToFile("U" + id + ".mtx");

                written = true;
                Console.Error.WriteLine("Written: " + id);
            }
            int L = X.Count;

            //double[] _B;
            //if(B.GetType() == typeof(double[]))
            //    _B = B as double[];
            //else
            //    _B = B.ToArray();


            double[] y = new double[L];
            //double[] yref = new double[L];
            //m_BlockL.Solve_Direct(yref, B.ToArray().CloneAs());
            LoTriDiagonalSolve(m_BlockL, m_BlockL_compressed, m_BlockL_OccupiedBlockColumnsPerRow, y, B, true);
            //double check1 = GenericBlas.L2Dist(yref, y);
            //Console.Error.WriteLine("Check value (lo solve) is: " + check1);


            //double[] Xref = new double[y.Length];
            //m_BlockU.Solve_Direct(Xref, y.CloneAs());
            UpTriDiagonalSolve(m_BlockU, m_BlockU_OccupiedBlockColumnsPerRow, X, y, m_BlockU_lublks, m_BlockU_luipiv);
            //double check2 = GenericBlas.L2Dist(Xref, X);
            //Console.Error.WriteLine("Check value (hi solve) is: " + check2);

            
            m_ThisLevelIterations += 1;
            Converged = true;
        }

        /// <summary>
        /// Lower triangular part of ILU-decomposition; Diagonal blocks are Identity matrices
        /// </summary>
        BlockMsrMatrix m_BlockL;

        /// <summary>
        /// 
        /// </summary>
        MultidimensionalArray[] m_BlockL_compressed;

        /// <summary>
        /// occupied block-columns for each block-row for <see cref="m_BlockL"/>
        /// </summary>
        long[][] m_BlockL_OccupiedBlockColumnsPerRow;

        /// <summary>
        /// Upper triangular part of ILU-decomposition; Diagonal blocks are fully occupied
        /// </summary>
        BlockMsrMatrix m_BlockU;

        /// <summary>
        /// occupied block-columns for each block-row for <see cref="m_BlockL"/>
        /// </summary>
        long[][] m_BlockU_OccupiedBlockColumnsPerRow;


        /// <summary>
        /// LU-decomposition of the diagonal blocks of <see cref="m_BlockU"/>
        /// </summary>
        MultidimensionalArray[] m_BlockU_lublks;


        /// <summary>
        /// Pivot indices for LU-decomposition of the diagonal blocks of <see cref="m_BlockU"/>
        /// </summary>
        int[][] m_BlockU_luipiv;

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
        /// ILU, in-place version, after book of Saad (p 303);
        /// This only considers the MPI-local part of the matrix, i.e. the ILU decomposition is MPI-rank diagonal.
        /// </summary>
        static (BlockMsrMatrix L, MultidimensionalArray[] BlockL_compressed, long[][] L_OccupiedBlockColumnsPerRow, BlockMsrMatrix U, long[][] U_OccupiedBlockColumnsPerRow, MultidimensionalArray[] U_lublks, int[][] U_luipiv) 
            UpdateILU(int iluLevel, BlockMsrMatrix Mtx) {
            using(var ft = new FuncTrace()) {
                IBlockPartitioning part = Mtx._RowPartitioning;
                Debug.Assert(part.EqualsPartition(Mtx._RowPartitioning), "mismatch in row partition");
                Debug.Assert(part.EqualsPartition(Mtx._ColPartitioning), "mismatch in column partition");
                long cell0 = part.FirstBlock;
                int J = part.LocalNoOfBlocks;


                var A = Mtx.CloneAs();

                // dbg_launch();
                var ILUp_pattern = GetPattern(iluLevel, Mtx, out int Occupied);
                csMPI.Raw.Barrier(Mtx.MPI_Comm);
                csMPI.Raw.Barrier(Mtx.MPI_Comm);
                csMPI.Raw.Barrier(Mtx.MPI_Comm);
                var ILUpT = ILUp_pattern.Transpose();
                double occupancy = (double)Occupied / (double)(J * J);
                //ft.Info($"CellILU, lv {m_op.LevelIndex}, ILU-{ILU_level}, occupancy = {Math.Round(occupancy * 100)}%.");

                MultidimensionalArray LU_Akk_T = null, 
                    Aik = null, 
                    Aij = null, Akj = null;

                for (long k = cell0; k < (cell0 + J - 1); k++) { // Iteration over matrix (rows and columns)
                    int sz_k = part.GetBlockLen(k);
                    if(sz_k <= 0)
                        continue; // cell i is empty

                    //var invAkk = A.GetBlock(k, k);
                    
                    //LU_Akk_T = A.GetBlock(k, k);
                    int Szk = part.GetBlockLen(k);
                    long Idxk = part.GetBlockI0(k);
                    LU_Akk_T = LU_Akk_T.ReuseTemp(Szk, Szk);
                    A.ReadBlock(Idxk, Idxk, LU_Akk_T);
                    LU_Akk_T.TransposeInPlace();

                    int[] _ipiv = new int[LU_Akk_T.NoOfRows];
                    try {
                        //invAkk.InvertInPlace();

                        LU_Akk_T.FactorizeLU(_ipiv);
                    } catch(ArithmeticException ae) {
                        Console.Error.WriteLine(ae.GetType() + ": " + ae.Message);
                        continue;
                    }

                    long[] occColumn_k = ILUpT.GetOccupiedColumnIndices(k);
                    foreach(long i in occColumn_k) { 
                    //for(long i = k + 1; i < (cell0 + J); i++) {
                    //    Debug.Assert(i > k);
                    //    if(P1[i, k]) {
                        if(i >= k + 1) {
                            
                            //Aik = A.GetBlock(i, k);
                            int Szi = part.GetBlockLen(i);
                            long Idxi = part.GetBlockI0(i);
                            Aik = Aik.ReuseTemp(Szi, Szk);
                            A.ReadBlock(Idxi, Idxk, Aik);
                            
                            // Compute: Aik = Aik*inv(Akk) via LU decomposition
                            Aik.TransposeInPlace();
                            LU_Akk_T.BacksubsLU(_ipiv, Aik, Aik);
                            Aik.TransposeInPlace();
                            A.SetBlock(Aik, i, k); 

                            long[] occRow_i = ILUp_pattern.GetOccupiedColumnIndices(i);
                            
                            // L(Bk,Bi) = U(Bk,Bi)*invUii
                            foreach(long j in occRow_i) { 
                                if(j >= k + 1) {
                                    //Aij = A.GetBlock(i, j);
                                    //Akj = A.GetBlock(k, j);
                                    int Szj = part.GetBlockLen(j);
                                    long Idxj = part.GetBlockI0(j);
                                    Aij = Aij.ReuseTemp(Szi, Szj);
                                    Akj = Akj.ReuseTemp(Szk, Szj);
                                    A.ReadBlock(Idxi, Idxj, Aij);
                                    A.ReadBlock(Idxk, Idxj, Akj);

                                    Aij.GEMM(-1, Aik, Akj, 1.0);
                                    A.SetBlock(Aij, i, j);
                                }
                            }

                        } else {
                            //if(occColumn_k.Contains(i))
                            //    throw new Exception($"fuck: P1[{i},{k}] = {P1[i, k]} //  {ILUpT[i, k]}  {m_ILUp_pattern[i, k]} ");
                        }
                    }
                }


                var BlockU = new BlockMsrMatrix(part, part);
                var BlockL = new BlockMsrMatrix(part, part);
                MultidimensionalArray[] BlockL_compressed = new MultidimensionalArray[J];
                BlockL.AccEyeSp(1.0);

                var Udiag_lublks = new MultidimensionalArray[J];
                int[][] Udiag_lubl_ipiv = new int[J][];
                long[][] L_OccupiedBlockColumnsPerRow = new long[J][];
                long[][] U_OccupiedBlockColumnsPerRow = new long[J][];

                for (long j = cell0; j < J + cell0; j++) {
                    //for(long i = cell0; i < J + cell0; i++) {
                    //    if(P1[i, j]) {
                    long[] occRow_j = ILUp_pattern.GetOccupiedColumnIndices(j);
                    Debug.Assert(A.GetOccupiedRowBlockIndices(j).SetEquals(occRow_j));

        

                    int Szj = part.GetBlockLen(j);
                    if (Szj > 0) {
                        int lenL = 0, lenU = 0;
                        for (int zz = 0; zz < occRow_j.Length; zz++) {
                            long i = occRow_j[zz];
                            if (zz > 0 && occRow_j[zz - 1] >= occRow_j[zz])
                                throw new ApplicationException("expecting strictly increasing data");

                            int Szi = part.GetBlockLen(i);
                            if (j > i) {
                                // lower tri
                                lenL += Szi;
                            }
                            if (i > j) {
                                // upper tri
                                lenU += Szi;
                            }
                        }
                        MultidimensionalArray BlockL_compressed_j = null;
                        if (lenL > 0) {
                            BlockL_compressed_j = MultidimensionalArray.Create(Szj, lenL);
                            BlockL_compressed[j] = BlockL_compressed_j;
                        }
                        int oL = 0;
                        long j0 = part.GetBlockI0(j);
                        for (int zz = 0; zz < occRow_j.Length; zz++) {
                            long i = occRow_j[zz];
                            int Szi = part.GetBlockLen(i);
                            long i0 = part.GetBlockI0(i);
                            if (Szi > 0) {
                                if (j > i) {
                                    A.ReadBlock(j0, i0, BlockL_compressed_j.ExtractSubArrayShallow(new int[] { 0, oL }, new int[] { Szj - 1, oL + Szi - 1 }));
                                    oL += Szi;
                                }
                                if (i > j) {

                                }
                            }
                        }
                        Debug.Assert(oL == (BlockL_compressed_j?.NoOfCols ?? 0));
                    }




                    foreach (long i in occRow_j) {
                        {
                            var Aji = A.GetBlock(j, i);
                            if(j > i) {
                                // lower tri
                                BlockL.SetBlock(Aji, j, i);
                                
                            } else {
                                // upper tri
                                BlockU.SetBlock(Aji, j, i);
                            }

                            if(i == j) {
                                int[] ipiv = new int[Aji.NoOfRows];
                                Aji.FactorizeLU(ipiv);
                                Udiag_lublks[j] = Aji;
                                Udiag_lubl_ipiv[j] = ipiv;
                            }
                        }
                        //} else {
                        //    var Aji = A.GetBlock(j, i);
                        //    double ErrNorm = Aji.L2Norm();
                        //    if(ErrNorm > 0)
                        //        throw new ArithmeticException();
                        //}
                    }

                    L_OccupiedBlockColumnsPerRow[j - cell0] = BlockL.GetOccupiedRowBlockIndices(j);
                    U_OccupiedBlockColumnsPerRow[j - cell0] = BlockU.GetOccupiedRowBlockIndices(j);
                }

                //m_BlockL.SaveToTextFileSparse("L.txt");
                //m_BlockU.SaveToTextFileSparse("U.txt");

                return (BlockL, BlockL_compressed, L_OccupiedBlockColumnsPerRow, BlockU, U_OccupiedBlockColumnsPerRow, Udiag_lublks, Udiag_lubl_ipiv);
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
        static private MsrMatrix GetPattern(int ILU_level, BlockMsrMatrix Mtx, out int Occupied) {
            using(var tr = new FuncTrace()) {
                //tr.InfoToConsole = true;

                if(ILU_level < 0)
                    throw new ArgumentException("ILU-Level must be >= 0, got " + ILU_level);
                if(ILU_level > Math.Min((int)(byte.MaxValue), 10))
                    throw new ArgumentException("ILU-Level must be <= 10, got " + ILU_level);

                //var grd = m_op.Mapping.GridData;
                //var Mtx = m_op.OperatorMatrix;
                IBlockPartitioning part = Mtx._RowPartitioning;
                Debug.Assert(part.EqualsPartition(Mtx._RowPartitioning), "mismatch in row partition");
                Debug.Assert(part.EqualsPartition(Mtx._ColPartitioning), "mismatch in column partition");
                long cell0 = part.FirstBlock;
                int J = part.LocalNoOfBlocks;

                //byte[,] lev = new byte[J, J];
                //lev.SetAll(byte.MaxValue);

                // determine ILU-0 pattern
                tr.Info("computing ILU(0) pattern...");
                var p = new Partitioning(J, Mtx.MPI_Comm);
                var ILU0_pattern = new MsrMatrix(p, p);
                for(long i = cell0; i < (J + cell0); i++) { // loop over block rows
                    //long iB = part.GetBlockI0(i);
                    //int sz_i = part.GetBlockLen(i);

                    long[] colBlocks = Mtx.GetOccupiedRowBlockIndices(i, alsoExternal: false);
                    foreach(long j in colBlocks) { // loop over block columns
                        var Blk = Mtx.GetBlock(i, j);
                        if(!Blk.IsZeroValued()) {
                            //lev[i, j] = 0;
                            ILU0_pattern[i, j] = 1;
                        } else {
                            /*
                             * seems to happen sometimes in the XDG case
                             * 
                             * 
                            int RowType = Mtx._RowPartitioning.GetBlockType(i);
                            int ColType = Mtx._ColPartitioning.GetBlockType(j);
                            int[] rowI0s = Mtx._RowPartitioning.GetSubblk_i0(RowType);
                            int[] colI0s = Mtx._ColPartitioning.GetSubblk_i0(ColType);
                            int[] rowLns = Mtx._RowPartitioning.GetSubblkLen(RowType);
                            int[] colLns = Mtx._ColPartitioning.GetSubblkLen(ColType);
                            int II = Mtx._RowPartitioning.GetBlockLen(i);
                            int JJ = Mtx._RowPartitioning.GetBlockLen(j);

                            tr.Info("Mem occupied, but the block is 0");
                            */
                        }
                    }
                }
                //ILU0pattern.SaveToTextFileSparse("ILU0.txt");

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
                //bool[,] ret = new bool[J, J];
                /*
                for(long j = cell0; j < (J + cell0); j++) {

                    long[] rowPattern = ILUp_pattern.GetOccupiedColumnIndices(j);
                    
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
                    
                }
                */
                //m_ILUp_pattern = ILUp_pattern;
                return ILUp_pattern;
            }
        }

        
        void CheckILU() {
            using(new FuncTrace()) {
                BlockMsrMatrix LxU = BlockMsrMatrix.Multiply(m_BlockL, m_BlockU);
                BlockMsrMatrix Lc = m_BlockL.CloneAs();
                BlockMsrMatrix Uc = m_BlockU.CloneAs();


                var P1 = GetPattern(this.ILU_level, m_op.OperatorMatrix, out _);


                //var grd = m_op.Mapping.GridData;

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

        static void LoTriDiagonalSolve<U, V>(
            BlockMsrMatrix Mtx,
            MultidimensionalArray[] m_BlockL_compressed,
            long[][] OccupiedBlockColumnsPerRow, U X, V B, bool diagEye)
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


                double[] Xtemp = null;

                for(long i = i0; i < i0 + J; i++) {
                    int sz_i = rowPart.GetBlockLen(i);
                    if(sz_i <= 0)
                        continue; // empty cell;
                    
                    long iIdx = rowPart.GetBlockI0(i);
                    int iIdxLoc = checked((int)(iIdx - i0Idx));

                    long[] occBlk_i = OccupiedBlockColumnsPerRow[i];
                    int NB = occBlk_i.Length;

                    double[] Bi = B.GetSubVector(iIdxLoc, sz_i);

                    if (NB > 0) {
                        var Lcompr = m_BlockL_compressed[i - i0];
                        int K = Lcompr?.NoOfCols ?? 0;
                        if (Xtemp == null || Xtemp.Length != K)
                            Array.Resize(ref Xtemp, K);
                        int oj = 0;
                        for (int jx = 0; jx < NB; jx++) {
                            long j = occBlk_i[jx];
                            if (j != i) {
                                long jIdx = colPart.GetBlockI0(j);
                                int jIdxLoc = checked((int)(jIdx - j0Idx));
                                int Szj = colPart.GetBlockLen(j);

                                for (int k = 0; k < Szj; k++)
                                    Xtemp[oj + k] = X[jIdxLoc + k];
                                oj += Szj;
                            }
                        }

                        if(K > 0)
                            m_BlockL_compressed[i].GEMV(-1.0, Xtemp, 1.0, Bi);

                    }


                    /*
                    //long[] jS = Mtx.GetOccupiedRowBlockIndices(i);
                    long[] jS = OccupiedBlockColumnsPerRow[i - i0];
                    foreach (long j in jS) {
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
                    */

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


        
        static void UpTriDiagonalSolve<U, V>(BlockMsrMatrix Mtx, long[][] OccupiedBlockColumnsPerRow, U X, V B, MultidimensionalArray[] Mtx_lublks, int[][] Mtx_luipiv)
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

                bool diagEye = Mtx_lublks == null;


                for (long i = i0 + J - 1; i >= i0; i--) {
                    int sz_i = rowPart.GetBlockLen(i);
                    if(sz_i <= 0)
                        continue; // empty cell;

                    if(colPart.GetBlockLen(i) <= 0)
                        throw new NotSupportedException("cannot do (upper) trigonal solve on non-quadratic zero diagonal block");

                    long iIdx = rowPart.GetBlockI0(i);
                    int iIdxLoc = checked((int)(iIdx - i0Idx));

                    double[] Bi = B.GetSubVector(iIdxLoc, sz_i);

                    //long[] jS = Mtx.GetOccupiedRowBlockIndices(i);
                    long[] jS = OccupiedBlockColumnsPerRow[i - i0];
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
                        //var Mtx_ii = Mtx.GetBlock(i, i);
                        //Mtx_ii.Solve(Xi, Bi);
                        Mtx_lublks[i - i0].BacksubsLU(Mtx_luipiv[i - i0], Xi, Bi);
                        X.SetSubVector<double, U, double[]>(Xi, iIdxLoc, sz_i);
                    }
                }
            }
        }


    }
}

using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
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
            m_backSubs = null;
            m_op = null;
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

                //(this.m_BlockL, this.m_BlockL_compressed, this.m_BlockL_OccupiedBlockColumnsPerRow,
                //    this.m_BlockU, this.m_BlockU_OccupiedBlockColumnsPerRow,
                //    this.m_BlockU_lublks, this.m_BlockU_luipiv) = ComputeILU(this.ILU_level, op.OperatorMatrix);
                //CheckILU();


                //m_Perm = GetPermMatrix(op.OperatorMatrix, true);
                //m_PermT = m_Perm.Transpose();
                //var op_perm = BlockMsrMatrix.Multiply(m_PermT, BlockMsrMatrix.Multiply(op.OperatorMatrix, m_Perm));
                //var ILU = ComputeILU(this.ILU_level, op_perm);

                
            }
        }



        //BlockMsrMatrix m_Perm;
        //BlockMsrMatrix m_PermT;

        BackSubs m_backSubs;

        /// <summary>
        /// 
        /// </summary>
        public void ResetStat() {
            m_ThisLevelIterations = 0;
            Converged = false;
        }

        //bool written = true;

        
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
                var test1 = new BackSubs_Reference(ComputeILU(0, opMtx));
                var test2 = new BackSubs_Reference(ComputeILU(0, opMtx));
                var (test_L1, test_U1) = (test1.BlockL, test1.BlockU);
                var (test_L2, test_U2) = (test2.BlockL, test2.BlockU);

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


        //BackSubs m_backSubsRef;

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double>  //
        {
            if(m_backSubs == null) {
                var ILU = ComputeILU(this.ILU_level, m_op.OperatorMatrix);
                m_backSubs = new BackSubs_Optimized_SinglePrec(ILU);
                //m_backSubs = new BackSubs_Optimized(ILU);
                //m_backSubsRef = new BackSubs_Reference(ILU);
            }
 
            /*
            if (!written) {
                
                var part = m_op.OperatorMatrix._RowPartitioning;
                if (part.MpiSize != 1)
                    throw new NotSupportedException();
                if (part.LocalNoOfBlocks != part.TotalNoOfBlocks)
                    throw new NotSupportedException();

                int NoOfBlocks = part.LocalNoOfBlocks;
                int[] BL = NoOfBlocks.ForLoop(iblk => part.GetBlockLen(iblk));
                BL.SaveToTextFile("part" + id + ".txt", part.MPI_Comm);

                var ILU = ComputeILU(this.ILU_level, BlockMsrMatrix.Multiply(m_PermT, BlockMsrMatrix.Multiply(m_op.OperatorMatrix, m_Perm)));
                var test = new BackSubs_Reference(ILU);


                m_op.OperatorMatrix.SaveToTextFileSparse("M" + id + ".txt");
                test.BlockL.SaveToTextFileSparse("L" + id + ".txt");
                test.BlockU.SaveToTextFileSparse("U" + id + ".txt");
                m_Perm.SaveToTextFileSparse("P" + id + ".txt");

                m_op.OperatorMatrix.ToMsrMatrix().SaveToFile("M" + id + ".mtx");
                test.BlockL.ToMsrMatrix().SaveToFile("L" + id + ".mtx");
                test.BlockU.ToMsrMatrix().SaveToFile("U" + id + ".mtx");

                written = true;
                Console.Error.WriteLine("Written: " + id);
                
            }*/
            int L = X.Count;

            //double[] _B = new double[B.Count];
            //m_PermT.SpMV(1.0, B, 0.0, _B);


            double[] y = new double[L];
            
            //double[] yref = new double[L];
            //m_backSubsRef.LoTriDiagonalSolve(yref, B.ToArray().CloneAs());
            m_backSubs.LoTriDiagonalSolve(y, B);
            //var ERR = yref.Minus(y);
            //double check1 = GenericBlas.L2Dist(yref, y);
            //Console.Error.WriteLine("Check value (lo solve) is: " + check1);


            //double[] Xref = new double[y.Length];
            //m_backSubsRef.UpTriDiagonalSolve(Xref, y.CloneAs());
            m_backSubs.UpTriDiagonalSolve(X, y);
            //double check2 = GenericBlas.L2Dist(Xref, X);
            //Console.Error.WriteLine("Check value (hi solve) is: " + check2);

            
            m_ThisLevelIterations += 1;
            Converged = true;
        }

        
        /// <summary>
        /// 
        /// </summary>
        public long UsedMemory() {
            return m_backSubs?.UsedMemory() ?? 0;
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

        abstract class BackSubs {

            public BackSubs(BlockMsrMatrix LUdecomp) {

            }

            public abstract void UpTriDiagonalSolve<U, V>(U X, V B)
                where U : IList<double>
                where V : IList<double>;
        

            public abstract void LoTriDiagonalSolve<U, V>(U X, V B)
                where U : IList<double>
                where V : IList<double>;

            abstract public long UsedMemory();
        }


        class BackSubs_Reference : BackSubs {

            IBlockPartitioning part;

            public BackSubs_Reference(BlockMsrMatrix A) : base(A) {
                part = A._RowPartitioning;
                if (!A._ColPartitioning.EqualsPartition(part))
                    throw new ArgumentException();
                
                

                BlockU = new BlockMsrMatrix(part, part);
                BlockL = new BlockMsrMatrix(part, part);
                BlockL.AccEyeSp(1.0);

                long cell0 = part.FirstBlock;
                int J = part.LocalNoOfBlocks;
                for (long j = cell0; j < J + cell0; j++) {
                    //for(long i = cell0; i < J + cell0; i++) {
                    //    if(P1[i, j]) {
                    long[] occRow_j = A.GetOccupiedRowBlockIndices(j);

                   


                    foreach (long i in occRow_j) {
                        {
                            var Aji = A.GetBlock(j, i);
                            if (j > i) {
                                // lower tri
                                BlockL.SetBlock(Aji, j, i);

                            } else {
                                // upper tri
                                BlockU.SetBlock(Aji, j, i);
                            }
                        }
                        //} else {
                        //    var Aji = A.GetBlock(j, i);
                        //    double ErrNorm = Aji.L2Norm();
                        //    if(ErrNorm > 0)
                        //        throw new ArithmeticException();
                        //}
                    }

                    //L_OccupiedBlockColumnsPerRow[j - cell0] = BlockL.GetOccupiedRowBlockIndices(j);
                    //U_OccupiedBlockColumnsPerRow[j - cell0] = BlockU.GetOccupiedRowBlockIndices(j);
                }
            }

            /*
            /// <summary>
            /// Global Column indices of non-zero blocks in the L-factor
            /// - 1st index: local block-row/cell index
            /// - 2nd index: enumeration
            /// </summary>
            long[][] L_OccupiedBlockColumnsPerRow;

            /// <summary>
            /// Global Column indices of non-zero blocks in the U-factor
            /// - 1st index: local block-row/cell index
            /// - 2nd index: enumeration
            /// </summary>
            long[][] U_OccupiedBlockColumnsPerRow;
            */

            /// <summary>
            /// L-factor of the ILU decomposition
            /// </summary>
            internal BlockMsrMatrix BlockL;

            /// <summary>
            /// U-factor of the ILU decomposition
            /// </summary>
            internal BlockMsrMatrix BlockU;

            public override void UpTriDiagonalSolve<U, V>(U X, V B) {
                using (new FuncTrace()) {
                    var Mtx = BlockU;

                    var rowPart = Mtx._RowPartitioning;
                    var colPart = Mtx._ColPartitioning;
                    long i0 = rowPart.FirstBlock;
                    long J = rowPart.LocalNoOfBlocks;
                    long i0Idx = rowPart.i0;
                    long j0Idx = colPart.i0;

                   
                    for (long i = i0 + J - 1; i >= i0; i--) {
                        int sz_i = rowPart.GetBlockLen(i);
                        if (sz_i <= 0)
                            continue; // empty cell;

                        if (colPart.GetBlockLen(i) <= 0)
                            throw new NotSupportedException("cannot do (upper) trigonal solve on non-quadratic zero diagonal block");

                        long iIdx = rowPart.GetBlockI0(i);
                        int iIdxLoc = checked((int)(iIdx - i0Idx));

                        double[] Bi = B.GetSubVector(iIdxLoc, sz_i);

                        long[] jS = Mtx.GetOccupiedRowBlockIndices(i);
                        foreach (long j in jS) {
                            if (j <= i) // ignore everything below the diagonal
                                continue;
                            int sz_j = colPart.GetBlockLen(j);
                            if (sz_j <= 0)
                                continue;

                            long jIdx = colPart.GetBlockI0(j);
                            int jIdxLoc = checked((int)(jIdx - j0Idx));

                            double[] Xj = X.GetSubVector(jIdxLoc, sz_j);
                            var Mtx_ij = Mtx.GetBlock(i, j);

                            Mtx_ij.GEMV(-1.0, Xj, 1.0, Bi);
                        }

                        //if (diagEye) {
                        //    X.SetSubVector<double, U, double[]>(Bi, iIdxLoc, sz_i);
                        //} else 
                        {
                            double[] Xi = new double[sz_i];
                            var Mtx_ii = Mtx.GetBlock(i, i);
                            Mtx_ii.Solve(Xi, Bi);
                            //Mtx_lublks[i - i0].BacksubsLU(Mtx_luipiv[i - i0], Xi, Bi);
                            //X.SetSubVector<double, U, double[]>(Xi, iIdxLoc, sz_i);
                        }
                    }
                }
            }

            public override void LoTriDiagonalSolve<U, V>(U X, V B) {
                using (new FuncTrace()) {
                    var Mtx = BlockL;

                    var rowPart = Mtx._RowPartitioning;
                    var colPart = Mtx._ColPartitioning;
                    long i0 = rowPart.FirstBlock;
                    long J = rowPart.LocalNoOfBlocks;
                    long i0Idx = rowPart.i0;
                    long j0Idx = colPart.i0;


                    for (long i = i0; i < i0 + J; i++) {
                        int sz_i = rowPart.GetBlockLen(i);
                        if (sz_i <= 0)
                            continue; // empty cell;

                        if (colPart.GetBlockLen(i) <= 0)
                            throw new NotSupportedException("cannot do (lower) trigonal solve on non-quadratic zero diagonal block");

                        long iIdx = rowPart.GetBlockI0(i);
                        int iIdxLoc = checked((int)(iIdx - i0Idx));

                        double[] Bi = B.GetSubVector(iIdxLoc, sz_i);

                        long[] jS = Mtx.GetOccupiedRowBlockIndices(i);
                        foreach (long j in jS) {
                            if (j >= i) // ignore everything above the diagonal
                                continue;
                            int sz_j = colPart.GetBlockLen(j);
                            if (sz_j <= 0)
                                continue;

                            long jIdx = colPart.GetBlockI0(j);
                            int jIdxLoc = checked((int)(jIdx - j0Idx));

                            double[] Xj = X.GetSubVector(jIdxLoc, sz_j);
                            var Mtx_ij = Mtx.GetBlock(i, j);

                            Mtx_ij.GEMV(-1.0, Xj, 1.0, Bi);
                        }

                        //if (diagEye) {
                        X.SetSubVector<double, U, double[]>(Bi, iIdxLoc, sz_i);
                        //} else {
                        //    double[] Xi = new double[sz_i];
                        //    var Mtx_ii = Mtx.GetBlock(i, i);
                        //    Mtx_ii.Solve(Xi, Bi);
                        //    X.SetSubVector<double, U, double[]>(Xi, iIdxLoc, sz_i);
                        //}
                    }
                }
            }

            public override long UsedMemory() {
                return BlockU.UsedMemory + BlockL.UsedMemory;
            }
        }


        class BackSubs_Optimized : BackSubs {

            IBlockPartitioning part;

            public BackSubs_Optimized(BlockMsrMatrix A) : base(A) {
                part = A._RowPartitioning;
                if (!A._ColPartitioning.EqualsPartition(part))
                    throw new ArgumentException();

                int J = part.LocalNoOfBlocks;
                this.m_BlockU_compressed = new MultidimensionalArray[J];
                this.m_BlockL_compressed = new MultidimensionalArray[J];

                this.Udiag_lublks = new MultidimensionalArray[J];
                this.Udiag_inverse = new MultidimensionalArray[J];
                this.Udiag_luipiv = new int[J][];
                this.L_OccupiedBlockColumnsPerRow = new long[J][];
                this.U_OccupiedBlockColumnsPerRow = new long[J][];

                long cell0 = part.FirstBlock;
                for (long j = cell0; j < J + cell0; j++) {
                    //for(long i = cell0; i < J + cell0; i++) {
                    //    if(P1[i, j]) {
                    long[] occRow_j = A.GetOccupiedRowBlockIndices(j);

                    int Szj = part.GetBlockLen(j);
                    if (Szj > 0) {
                        int lenL = 0, lenU = 0, zzL = 0, zzU = 0;
                        for (int zz = 0; zz < occRow_j.Length; zz++) {
                            long i = occRow_j[zz];
                            if (zz > 0 && occRow_j[zz - 1] >= occRow_j[zz])
                                throw new ApplicationException("expecting strictly increasing data");

                            int Szi = part.GetBlockLen(i);
                            if (j > i) {
                                // lower tri
                                lenL += Szi;
                                zzL++;
                            }
                            if (i > j) {
                                // upper tri
                                lenU += Szi;
                                zzU++;
                            }
                        }
                        MultidimensionalArray BlockL_compressed_j = null;
                        MultidimensionalArray BlockU_compressed_j = null;
                        if (lenU > 0) {
                            BlockU_compressed_j = MultidimensionalArray.Create(Szj, lenU);
                            this.m_BlockU_compressed[j - cell0] = BlockU_compressed_j;
                        }
                        if (lenL > 0) {
                            BlockL_compressed_j = MultidimensionalArray.Create(Szj, lenL);
                            m_BlockL_compressed[j - cell0] = BlockL_compressed_j;
                        }
                        long[] BlockL_GetOccupiedRowBlockIndices = new long[zzL];
                        long[] BlockU_GetOccupiedRowBlockIndices = new long[zzU];

                        int oL = 0, oU = 0;
                        zzL = 0; zzU = 0;
                        long j0 = part.GetBlockI0(j);
                        for (int zz = 0; zz < occRow_j.Length; zz++) {
                            long i = occRow_j[zz];
                            int Szi = part.GetBlockLen(i);
                            long i0 = part.GetBlockI0(i);
                            if (Szi > 0) {
                                if (j > i) {
                                    A.ReadBlock(j0, i0, BlockL_compressed_j.ExtractSubArrayShallow(new int[] { 0, oL }, new int[] { Szj - 1, oL + Szi - 1 }));
                                    oL += Szi;
                                    BlockL_GetOccupiedRowBlockIndices[zzL] = i; zzL++;
                                }
                                if (i > j) {
                                    A.ReadBlock(j0, i0, BlockU_compressed_j.ExtractSubArrayShallow(new int[] { 0, oU }, new int[] { Szj - 1, oU + Szi - 1 }));
                                    oU += Szi;
                                    BlockU_GetOccupiedRowBlockIndices[zzU] = i; zzU++;
                                }
                            }
                        }
                        Debug.Assert(oL == (BlockL_compressed_j?.NoOfCols ?? 0));
                        Debug.Assert(oU == (BlockU_compressed_j?.NoOfCols ?? 0));


                        L_OccupiedBlockColumnsPerRow[j - cell0] = BlockL_GetOccupiedRowBlockIndices;
                        U_OccupiedBlockColumnsPerRow[j - cell0] = BlockU_GetOccupiedRowBlockIndices;


                        {

                            var Ajj = A.GetBlock(j, j);
                            this.Udiag_inverse[j - cell0] = Ajj.CloneAs();
                            this.Udiag_inverse[j - cell0].InvertInPlace();

                            int[] ipiv = new int[Ajj.NoOfRows];
                            Ajj.FactorizeLU(ipiv);
                            this.Udiag_lublks[j - cell0] = Ajj;
                            this.Udiag_luipiv[j - cell0] = ipiv;
                            

                        }
                    }

                }
            }

            /// <summary>
            /// Concatenation of all upper diagonal parts of the U-factor for each block row
            /// </summary>
            internal MultidimensionalArray[] m_BlockU_compressed;

            /// <summary>
            /// Concatenation of all lower diagonal parts of the L-factor for each block row
            /// </summary>
            internal MultidimensionalArray[] m_BlockL_compressed;

            /// <summary>
            /// Global Column indices of non-zero blocks in the L-factor
            /// - 1st index: local block-row/cell index
            /// - 2nd index: enumeration
            /// </summary>
            long[][] L_OccupiedBlockColumnsPerRow;

            /// <summary>
            /// Global Column indices of non-zero blocks in the U-factor
            /// - 1st index: local block-row/cell index
            /// - 2nd index: enumeration
            /// </summary>
            long[][] U_OccupiedBlockColumnsPerRow;

            /// <summary>
            /// LU-decompositions for each diagonal block in the U-factor (upper triangle)
            /// </summary>
            MultidimensionalArray[] Udiag_lublks;

            /// <summary>
            /// inverse for each diagonal block in the U-factor (upper triangle)
            /// </summary>
            MultidimensionalArray[] Udiag_inverse;

            /// <summary>
            /// pivot indices corresponding with <see cref="Udiag_lublks"/>
            /// </summary>
            int[][] Udiag_luipiv;


            public override long UsedMemory() {
                long ret = 0;
                int J = part.LocalNoOfBlocks;
                for(int j = 0; j < J; j++) {
                    ret += (m_BlockU_compressed[j]?.Length ?? 0) * sizeof(double);
                    ret += (m_BlockL_compressed[j]?.Length ?? 0) * sizeof(double);

                    ret += (L_OccupiedBlockColumnsPerRow[j]?.Length ?? 0) * sizeof(long);
                    ret += (U_OccupiedBlockColumnsPerRow[j]?.Length ?? 0) * sizeof(long);

                    ret += (Udiag_lublks[j]?.Length ?? 0) * sizeof(double);
                    ret += (Udiag_luipiv[j]?.Length ?? 0) * sizeof(int);
                }


                return ret;
            }


            public override void UpTriDiagonalSolve<U, V>(U X, V B) {
                using (new FuncTrace()) {
                    var rowPart = part;
                    var colPart = part;
                    long i0 = rowPart.FirstBlock;
                    long J = rowPart.LocalNoOfBlocks;
                    long i0Idx = rowPart.i0;
                    long j0Idx = colPart.i0;


                    double[] Xtemp = null;

                    for (long i = i0 + J - 1; i >= i0; i--) { // reverse loop over block-rows
                        int sz_i = rowPart.GetBlockLen(i);
                        if (sz_i <= 0)
                            continue; // empty cell;

                        long iIdx = rowPart.GetBlockI0(i); // global matrix row index
                        int iIdxLoc = checked((int)(iIdx - i0Idx)); // local row index

                        long[] occBlk_i = U_OccupiedBlockColumnsPerRow[i - i0];
                        int NB = occBlk_i.Length;

                        double[] Bi = B.GetSubVector(iIdxLoc, sz_i);

                        if (NB > 0) {
                            var Ucompr = m_BlockU_compressed[i - i0];
                            int K = Ucompr?.NoOfCols ?? 0;
                            if (K > 0) {
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

                                
                                m_BlockU_compressed[i - i0].GEMV(-1.0, Xtemp, 1.0, Bi);
                                

                                
                                //// testcode for single precision
                                //{
                                //    var MTX = m_BlockU_compressed[i - i0];
                                //    int N = MTX.NoOfRows;
                                //    int M = MTX.NoOfCols;
                                //    if (M != Xtemp.Length)
                                //        throw new ArgumentException();

                                //    for(int n = 0; n < N; n++) {
                                //        float acc = 0;
                                //        for(int m = 0; m < M; m++) {
                                //            float MTX_nm = (float)MTX[n, m];
                                //            float X_m = (float) Xtemp[m];
                                //            acc += MTX_nm * X_m;
                                //        }

                                //        Bi[n] += (-1.0) * acc;
                                //    }
                                //}

                            }
                        }


                        {
                            double[] Xi = new double[sz_i];
                            //var Mtx_ii = Mtx.GetBlock(i, i);
                            //Mtx_ii.Solve(Xi, Bi);
                            //this.Udiag_lublks[i - i0].BacksubsLU(this.Udiag_luipiv[i - i0], Xi, Bi);

                            if (Bi.CheckForNanOrInfV(ExceptionIfFound: false) > 0) {
                                Console.Error.WriteLine("ILU breakdown in U1 slove");
                                return;
                            }

                            this.Udiag_inverse[i - i0].GEMV(1.0, Bi, 1.0, Xi);

                            if (Xi.CheckForNanOrInfV(ExceptionIfFound: false) > 0) {
                                Console.Error.WriteLine("ILU breakdown in U2 slove");
                                return;
                            }

                            X.SetSubVector<double, U, double[]>(Xi, iIdxLoc, sz_i);
                        }

                    }
                }
            }

            public override void LoTriDiagonalSolve<U, V>(U X, V B) {
                using (new FuncTrace()) {
                    var rowPart = part;
                    var colPart = part;
                    long i0 = rowPart.FirstBlock;
                    long J = rowPart.LocalNoOfBlocks;
                    long i0Idx = rowPart.i0;
                    long j0Idx = colPart.i0;


                    double[] Xtemp = null;

                    for (long i = i0; i < i0 + J; i++) { // loop over block-rows
                        int sz_i = rowPart.GetBlockLen(i);
                        if (sz_i <= 0)
                            continue; // empty cell;

                        long iIdx = rowPart.GetBlockI0(i); // global matrix row index
                        int iIdxLoc = checked((int)(iIdx - i0Idx)); // local row index

                        long[] occBlk_i = L_OccupiedBlockColumnsPerRow[i - i0];
                        int NB = occBlk_i.Length;

                        double[] Bi = B.GetSubVector(iIdxLoc, sz_i);

                        if (NB > 0) {
                            var Lcompr = m_BlockL_compressed[i - i0];
                            int K = Lcompr?.NoOfCols ?? 0;
                            if (K > 0) {
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

                                m_BlockL_compressed[i - i0].GEMV(-1.0, Xtemp, 1.0, Bi);

                                
                                //// testcode for single precision
                                //{
                                //    var MTX = m_BlockL_compressed[i - i0];
                                //    int N = MTX.NoOfRows;
                                //    int M = MTX.NoOfCols;
                                //    if (M != Xtemp.Length)
                                //        throw new ArgumentException();

                                //    for (int n = 0; n < N; n++) {
                                //        float acc = 0;
                                //        for (int m = 0; m < M; m++) {
                                //            float MTX_nm = (float)MTX[n, m];
                                //            float X_m = (float)Xtemp[m];
                                //            acc += MTX_nm * X_m;
                                //        }

                                //        Bi[n] += (-1.0) * acc;
                                //    }
                                //}
                            }
                        }

                        if (Bi.CheckForNanOrInfV(ExceptionIfFound:false) > 0) {
                            Console.Error.WriteLine("ILU breakdown in L slove");
                            return;
                        }
                        X.SetSubVector<double, U, double[]>(Bi, iIdxLoc, sz_i);
                        
                    }
                }
            }
        }


        class BackSubs_Optimized_SinglePrec : BackSubs {

            IBlockPartitioning part;

            //BackSubs_Optimized m_backSubsRef;

            public BackSubs_Optimized_SinglePrec(BlockMsrMatrix A) : base(A) {
                part = A._RowPartitioning;
                if (!A._ColPartitioning.EqualsPartition(part))
                    throw new ArgumentException();
                //m_backSubsRef = new BackSubs_Optimized(A.CloneAs());

                int J = part.LocalNoOfBlocks;
                this.m_BlockU_compressed = new float[J][,];
                this.m_BlockL_compressed = new float[J][,];

                this.Udiag_inverse = new float[J][,];
                this.L_OccupiedBlockColumnsPerRow = new long[J][];
                this.U_OccupiedBlockColumnsPerRow = new long[J][];

                MultidimensionalArray tempAji = null;
                max_sz_i = 0;
                max_row = 0;

                long cell0 = part.FirstBlock;
                for (long j = cell0; j < J + cell0; j++) {
                    //for(long i = cell0; i < J + cell0; i++) {
                    //    if(P1[i, j]) {
                    long[] occRow_j = A.GetOccupiedRowBlockIndices(j);

                    int Szj = part.GetBlockLen(j);
                    max_sz_i = Math.Max(max_sz_i, Szj);
                    if (Szj > 0) {
                        int lenL = 0, lenU = 0, zzL = 0, zzU = 0;
                        for (int zz = 0; zz < occRow_j.Length; zz++) {
                            long i = occRow_j[zz];
                            if (zz > 0 && occRow_j[zz - 1] >= occRow_j[zz])
                                throw new ApplicationException("expecting strictly increasing data");

                            int Szi = part.GetBlockLen(i);
                            if (j > i) {
                                // lower tri
                                lenL += Szi;
                                zzL++;
                            }
                            if (i > j) {
                                // upper tri
                                lenU += Szi;
                                zzU++;
                            }
                        }
                        float[,] BlockL_compressed_j = null;
                        float[,] BlockU_compressed_j = null;
                        if (lenU > 0) {
                            BlockU_compressed_j = new float[Szj, lenU];
                            this.m_BlockU_compressed[j - cell0] = BlockU_compressed_j;
                        }
                        if (lenL > 0) {
                            BlockL_compressed_j = new float[Szj, lenL];
                            m_BlockL_compressed[j - cell0] = BlockL_compressed_j;
                        }
                        max_row = Math.Max(max_row, lenU);
                        max_row = Math.Max(max_row, lenL);
                        long[] BlockL_GetOccupiedRowBlockIndices = new long[zzL];
                        long[] BlockU_GetOccupiedRowBlockIndices = new long[zzU];

                        int oL = 0, oU = 0;
                        zzL = 0; zzU = 0;
                        long j0 = part.GetBlockI0(j);
                        for (int zz = 0; zz < occRow_j.Length; zz++) {
                            long i = occRow_j[zz];
                            int Szi = part.GetBlockLen(i);
                            long i0 = part.GetBlockI0(i);
                            if (Szi > 0) {
                                if (j > i) {
                                    tempAji = tempAji.ReuseTemp(Szj, Szi);
                                    //A.ReadBlock(j0, i0, BlockL_compressed_j.ExtractSubArrayShallow(new int[] { 0, oL }, new int[] { Szj - 1, oL + Szi - 1 }));
                                    A.ReadBlock(j0, i0, tempAji);
                                    for (int n = 0; n < Szj; n++)
                                        for (int m = 0; m < Szi; m++)
                                            BlockL_compressed_j[n, m + oL] = (float)tempAji[n, m];
                                    oL += Szi;
                                    BlockL_GetOccupiedRowBlockIndices[zzL] = i; zzL++;
                                }
                                if (i > j) {
                                    tempAji = tempAji.ReuseTemp(Szj, Szi);
                                    //A.ReadBlock(j0, i0, BlockU_compressed_j.ExtractSubArrayShallow(new int[] { 0, oU }, new int[] { Szj - 1, oU + Szi - 1 }));
                                    A.ReadBlock(j0, i0, tempAji);
                                    for (int n = 0; n < Szj; n++)
                                        for (int m = 0; m < Szi; m++)
                                            BlockU_compressed_j[n, m + oU] = (float)tempAji[n, m];
                                    oU += Szi;
                                    BlockU_GetOccupiedRowBlockIndices[zzU] = i; zzU++;
                                }
                            }
                        }
                        Debug.Assert(oL == (BlockL_compressed_j?.GetLength(1) ?? 0));
                        Debug.Assert(oU == (BlockU_compressed_j?.GetLength(1) ?? 0));


                        L_OccupiedBlockColumnsPerRow[j - cell0] = BlockL_GetOccupiedRowBlockIndices;
                        U_OccupiedBlockColumnsPerRow[j - cell0] = BlockU_GetOccupiedRowBlockIndices;

                        {

                            var invAjj = A.GetBlock(j, j);
                            invAjj.InvertInPlace();
                            Debug.Assert(invAjj.NoOfRows == Szj);
                            Debug.Assert(invAjj.NoOfCols == Szj);
                            float[,] invAjjf = new float[Szj, Szj];
                            this.Udiag_inverse[j - cell0] = invAjjf;
                            for (int n = 0; n < Szj; n++)
                                for (int m = 0; m < Szj; m++)
                                    invAjjf[n, m] = (float)invAjj[n, m];

                            //int[] ipiv = new int[Ajj.NoOfRows];
                            //Ajj.FactorizeLU(ipiv);
                            //this.Udiag_lublks[j - cell0] = Ajj;
                            //this.Udiag_luipiv[j - cell0] = ipiv;


                        }
                    }

                }
            }

            /// <summary>
            /// Concatenation of all upper diagonal parts of the U-factor for each block row
            /// </summary>
            float[][,] m_BlockU_compressed;

            /// <summary>
            /// Concatenation of all lower diagonal parts of the L-factor for each block row
            /// </summary>
            float[][,] m_BlockL_compressed;

            /// <summary>
            /// Global Column indices of non-zero blocks in the L-factor
            /// - 1st index: local block-row/cell index
            /// - 2nd index: enumeration
            /// </summary>
            long[][] L_OccupiedBlockColumnsPerRow;

            /// <summary>
            /// Global Column indices of non-zero blocks in the U-factor
            /// - 1st index: local block-row/cell index
            /// - 2nd index: enumeration
            /// </summary>
            long[][] U_OccupiedBlockColumnsPerRow;

            /// <summary>
            /// inverse for each diagonal block in the U-factor (upper triangle)
            /// </summary>
            float[][,] Udiag_inverse;



            public override long UsedMemory() {
                long ret = 0;
                int J = part.LocalNoOfBlocks;
                for (int j = 0; j < J; j++) {
                    ret += (m_BlockU_compressed[j]?.Length ?? 0) * sizeof(double);
                    ret += (m_BlockL_compressed[j]?.Length ?? 0) * sizeof(double);

                    ret += (L_OccupiedBlockColumnsPerRow[j]?.Length ?? 0) * sizeof(long);
                    ret += (U_OccupiedBlockColumnsPerRow[j]?.Length ?? 0) * sizeof(long);


                }


                return ret;
            }

            int max_sz_i;
            int max_row;

            public override void UpTriDiagonalSolve<U, V>(U X, V B) {
                using (new FuncTrace()) {
                    var rowPart = part;
                    var colPart = part;
                    long i0 = rowPart.FirstBlock;
                    long J = rowPart.LocalNoOfBlocks;
                    long i0Idx = rowPart.i0;
                    long j0Idx = colPart.i0;

                    unsafe {
                        float* Xi = stackalloc float[max_sz_i];
                        float* Bi = stackalloc float[max_sz_i];
                        float* Xtemp = stackalloc float[max_row];

                        for (long i = i0 + J - 1; i >= i0; i--) { // reverse loop over block-rows
                            int sz_i = rowPart.GetBlockLen(i);
                            if (sz_i <= 0)
                                continue; // empty cell;

                            long iIdx = rowPart.GetBlockI0(i); // global matrix row index
                            int iIdxLoc = checked((int)(iIdx - i0Idx)); // local row index

                            long[] occBlk_i = U_OccupiedBlockColumnsPerRow[i - i0];
                            int NB = occBlk_i.Length;

                            for (int n = 0; n < sz_i; n++)
                                Bi[n] = (float)B[iIdxLoc + n];

                            if (NB > 0) {
                                var Ucompr = m_BlockU_compressed[i - i0];
                                int K = Ucompr?.GetLength(1) ?? 0;
                                if (K > 0) {

                                    int oj = 0;
                                    for (int jx = 0; jx < NB; jx++) {
                                        long j = occBlk_i[jx];
                                        if (j != i) {
                                            long jIdx = colPart.GetBlockI0(j);
                                            int jIdxLoc = checked((int)(jIdx - j0Idx));
                                            int Szj = colPart.GetBlockLen(j);

                                            for (int k = 0; k < Szj; k++)
                                                Xtemp[oj + k] = (float)X[jIdxLoc + k];
                                            oj += Szj;
                                        }
                                    }


                                    fixed (float* pMtx = Ucompr) {
                                        BLAS.sgemv('T', K, sz_i, -1.0f, pMtx, K, Xtemp, 1, 1.0f, Bi, 1);
                                    }    
                                }
                            }



                            Debug.Assert(this.Udiag_inverse[i - i0].GetLength(0) == sz_i);
                            Debug.Assert(this.Udiag_inverse[i - i0].GetLength(1) == sz_i);
                            fixed (float* pInvUdiag = this.Udiag_inverse[i - i0]) {
                                BLAS.sgemv('T', sz_i, sz_i, 1.0f, pInvUdiag, sz_i, Bi, 1, 0.0f, Xi, 1);
                            }

                            for (int n = 0; n < sz_i; n++)
                                X[iIdxLoc + n] = Xi[n];


                        }
                    }
                }
            }

          
            public override void LoTriDiagonalSolve<U, V>(U X, V B) {
                using (new FuncTrace()) {
                    var rowPart = part;
                    var colPart = part;
                    long i0 = rowPart.FirstBlock;
                    long J = rowPart.LocalNoOfBlocks;
                    long i0Idx = rowPart.i0;
                    long j0Idx = colPart.i0;

                    unsafe {
                        float* Xtemp = stackalloc float[max_row];
                        float* Bi = stackalloc float[max_sz_i];

                        for (long i = i0; i < i0 + J; i++) { // loop over block-rows
                            int sz_i = rowPart.GetBlockLen(i);
                            if (sz_i <= 0)
                                continue; // empty cell;

                            long iIdx = rowPart.GetBlockI0(i); // global matrix row index
                            int iIdxLoc = checked((int)(iIdx - i0Idx)); // local row index

                            long[] occBlk_i = L_OccupiedBlockColumnsPerRow[i - i0];
                            int NB = occBlk_i.Length;

                            for (int n = 0; n < sz_i; n++)
                                Bi[n] = (float)B[iIdxLoc + n];

                            if (NB > 0) {
                                var Lcompr = m_BlockL_compressed[i - i0];
                                //var Lcompr_Ref = m_backSubsRef.m_BlockL_compressed[i - i0];
                                int K = Lcompr?.GetLength(1) ?? 0;
                                if (K > 0) {
                                    //if (Xtemp == null || Xtemp.Length != K)
                                    //    Array.Resize(ref Xtemp, K);
                                    int oj = 0;
                                    for (int jx = 0; jx < NB; jx++) {
                                        long j = occBlk_i[jx];
                                        if (j != i) {
                                            long jIdx = colPart.GetBlockI0(j);
                                            int jIdxLoc = checked((int)(jIdx - j0Idx));
                                            int Szj = colPart.GetBlockLen(j);

                                            for (int k = 0; k < Szj; k++)
                                                Xtemp[oj + k] = (float) X[jIdxLoc + k];
                                            oj += Szj;
                                        }
                                    }

                                    fixed (float* pMtx = Lcompr) {
                                        BLAS.sgemv('T', K, sz_i, -1.0f, pMtx, K, Xtemp, 1, 1.0f, Bi, 1);
                                    }
                                }
                            }
                            for (int n = 0; n < sz_i; n++)
                                X[iIdxLoc + n] = Bi[n];
                        }
                    }
                }
            }
        }

        /// <summary>
        /// ILU, in-place version, after book of Saad (p 303);
        /// This only considers the MPI-local part of the matrix, i.e. the ILU decomposition is MPI-rank diagonal.
        /// </summary>
        static BlockMsrMatrix ComputeILU(int iluLevel, BlockMsrMatrix Mtx) {
            using(var ft = new FuncTrace()) {
                IBlockPartitioning part = Mtx._RowPartitioning;
                Debug.Assert(part.EqualsPartition(Mtx._RowPartitioning), "mismatch in row partition");
                Debug.Assert(part.EqualsPartition(Mtx._ColPartitioning), "mismatch in column partition");
                long cell0 = part.FirstBlock;
                int J = part.LocalNoOfBlocks;


                var A = Mtx.CloneAs();

                // dbg_launch();
                var ILUp_pattern = GetPattern(iluLevel, Mtx);
                csMPI.Raw.Barrier(Mtx.MPI_Comm);
                var ILUpT = ILUp_pattern.Transpose();
                
                MultidimensionalArray LU_Akk_T = null, 
                    Aik = null, AikT = null,
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
                        LU_Akk_T.CheckForNanOrInf();
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
                            AikT = AikT.ReuseTemp(Szk, Szi);
                            Aik.TransposeTo(AikT);
                            LU_Akk_T.BacksubsLU(_ipiv, AikT, AikT);
                            AikT.TransposeTo(Aik);
                            Aik.CheckForNanOrInf();
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

                return A;

                
            }
        }

       

        public int ILU_level = 1;


        static int[] RandomPermute(int JC) {
            // random permute of aggregation cells
            int[] R = new int[JC];
            R.SetAll(-1);
            Random rnd = new Random();
            //Console.Write(" perm: ");
            for (int jc = 0; jc < JC; jc++) {
                int jDest = rnd.Next(JC);
                while (R[jDest] >= 0) {
                    jDest = (jDest + 1) % JC;
                }

                //Console.Write($" {jc}>{jDest}");
                R[jDest] = jc;
            }
            //Console.WriteLine();

            return R;
        }

        static private int[] CuthillMcKeyPerm(BlockMsrMatrix Mtx, bool reverse) {
            if (!Mtx._RowPartitioning.EqualsPartition(Mtx._ColPartitioning))
                throw new ApplicationException();

            var part = Mtx._RowPartitioning;
            MsrMatrix Adj = GetPattern(0, Mtx);

            int JC = part.LocalNoOfBlocks;

            BitArray added = new BitArray(JC);
            int[] AdjRest(int jc) {
                long i0Adj = Adj.RowPartitioning.i0;
                long[] ColIdx = null;
                Adj.GetOccupiedColumnIndices(jc + i0Adj, ref ColIdx);

                return ColIdx.Where(jneigh => jneigh != (jc + i0Adj)  && !added[checked((int)(jneigh - i0Adj))])
                    .Select(jneigh => checked((int)(jneigh - i0Adj)))
                    .ToArray();
            }

            int Degree(int jc) {
                long i0Adj = Adj.RowPartitioning.i0;
                return Adj.GetNoOfOffDiagonalNonZerosPerRow(jc + i0Adj);
            }

            List<int> R = new List<int>(new int[] { 0 }); added[0] = true;
            for (int i = 0; R.Count < JC; i++) {
                int Ri = R[i];

                int[] Ai = AdjRest(Ri);
                Debug.Assert(Ai.Where(a => added[a]).Count() == 0);
                if (Ai.Length > 0) {
                    int[] Degs = Ai.Select(a => Degree(a)).ToArray();

                    Array.Sort(Degs, Ai);

                    R.AddRange(Ai);
                    foreach (var k in Ai)
                        added[k] = true;
                }
            }

            if (R.Count != JC)
                throw new ApplicationException("Cuthill-McKey internal error.");

            /*
            int[][] ret = new int[JC][];
            for (int j = 0; j < JC; j++) {
                if (reverse)
                    ret[JC - j - 1] = AggCells[R[j]];
                else
                    ret[j] = AggCells[R[j]];
            }*/
            if (reverse) {
                int[] ret = new int[JC];
                for (int j = 0; j < JC; j++) {
                    ret[JC - j - 1] = R[j];
                }
                return ret;
            } else {
                return R.ToArray();
            }
        }


        static private BlockMsrMatrix GetPermMatrix(BlockMsrMatrix Mtx, bool reverse) {
            if (!Mtx._RowPartitioning.EqualsPartition(Mtx._ColPartitioning))
                throw new ApplicationException();

            var part = Mtx._RowPartitioning;
            int J = Mtx._RowPartitioning.LocalNoOfBlocks;

            int[] perm = CuthillMcKeyPerm(Mtx, reverse);
            //int[] perm = RandomPermute(J);
            //perm = J.ForLoop(j => J - j - 1);
            BitArray check = new BitArray(J);
            for(int j = 0; j < J; j++) {
                if (check[perm[j]])
                    throw new ArgumentException("double");
                check[perm[j]] = true;
            }
            for (int j = 0; j < J; j++) {
                if (!check[perm[j]])
                    throw new ArgumentException("missing");
            }

            BlockMsrMatrix P = new BlockMsrMatrix(part, part);


            //perm = J.ForLoop(i => i);

            for(int j = 0; j < J; j++) {
                int i = perm[j];
                //long row_idx = part.GetBlockI0(j0 + j);
                //long col_idx = part.GetBlockI0(j0 + i);

                int NoRow = part.GetBlockLen(j);
                int NoCol = part.GetBlockLen(i);

                var Eye = MultidimensionalArray.Create(NoRow, NoCol);
                Eye.AccEye(1.0);
                P.SetBlock(Eye, j + part.FirstBlock, i + part.FirstBlock);
            }

            double[] test = new double[part.LocalLength];
            long partI0 = part.i0;
            for (int j = 0; j < J; j++) {
                int NoRow = part.GetBlockLen(j);
                long BlkI0 = part.GetBlockI0(j);
                for(int i = 0; i < NoRow; i++) {
                    test[i + BlkI0 - partI0] = j;
                }
            }

            double[] testOut = new double[test.Length];
            P.SpMV(1.0, test, 0.0, testOut);
            for (int j = 0; j < J; j++) {
                int NoRow = part.GetBlockLen(j);
                long BlkI0 = part.GetBlockI0(j);
                for (int i = 0; i < NoRow; i++) {
                    if(testOut[i + BlkI0 - partI0] != perm[j]) {
                        throw new ArgumentException("perm matrix fubar");
                    }
                }
            }

            return P;
        }



        /// <summary>
        /// Obtain occupancy pattern of a <see cref="BlockMsrMatrix"/>;
        /// this is also an adjacency matrix.
        /// </summary>
        /// <returns>
        /// A matrix with one entry for each block of the input matrix <paramref name="Mtx"/>
        /// </returns>
        static private MsrMatrix GetPattern(int ILU_level, BlockMsrMatrix Mtx) {
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

        
        /// <summary>
        /// Verifies that, within the occupancy pattern of <paramref name="OperatorMatrix"/>,
        /// the product of <paramref name="m_BlockU"/>*<paramref name="m_BlockL"/> is equal to the <paramref name="OperatorMatrix"/>
        /// </summary>
        static void CheckILU(int ILU_level, BlockMsrMatrix OperatorMatrix, BlockMsrMatrix m_BlockU, BlockMsrMatrix m_BlockL) {
            using(new FuncTrace()) {
                BlockMsrMatrix LxU = BlockMsrMatrix.Multiply(m_BlockL, m_BlockU);
                BlockMsrMatrix Lc = m_BlockL.CloneAs();
                BlockMsrMatrix Uc = m_BlockU.CloneAs();


                var P1 = GetPattern(ILU_level, OperatorMatrix);


                //var grd = m_op.Mapping.GridData;

                BlockMsrMatrix ErrMtx = LxU.CloneAs();
                ErrMtx.Acc(-1, OperatorMatrix);
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
                            OperatorMatrix.ReadBlock(jRow, kRow, B);

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

     

        
   

    }
}

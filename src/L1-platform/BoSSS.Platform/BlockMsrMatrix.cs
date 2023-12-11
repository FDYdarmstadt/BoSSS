﻿/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
using System;
using System.Collections.Generic;
using System.Diagnostics;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System.Collections;
using System.Linq;
using System.Runtime.InteropServices;
using ilPSP.Connectors;

namespace ilPSP.LinSolvers {

    /// <summary>
    /// Extension functions for the <see cref="BlockMsrMatrix"/>.
    /// </summary>
    public static class BMext {

        /// <summary>
        /// Driver routine for <see cref="BlockMsrMatrix.AccSubMatrixTo{V1, V2, V3, V4}(double, IMutableMatrixEx, V1, V2, V3, V4)"/>
        /// </summary>
        static public BlockMsrMatrix GetSubMatrix<V1, V3>(this BlockMsrMatrix OrgMtx, V1 RowIndicesSource, V3 ColumnIndiceSource)
            where V1 : IList<long>
            where V3 : IList<long> //
        {
            return GetSubMatrix(OrgMtx, RowIndicesSource, ColumnIndiceSource, csMPI.Raw._COMM.WORLD);
        }


        /// <summary>
        /// Driver routine for <see cref="BlockMsrMatrix.AccSubMatrixTo{V1, V2, V3, V4}(double, IMutableMatrixEx, V1, V2, V3, V4)"/>
        /// </summary>
        static public BlockMsrMatrix GetSubMatrix<V1, V3>(this BlockMsrMatrix OrgMtx, V1 RowIndicesSource, V3 ColumnIndiceSource, MPI_Comm destComm)
            where V1 : IList<long>
            where V3 : IList<long> //
        {

            var Blk1 = OrgMtx._RowPartitioning.GetSubBlocking(RowIndicesSource, destComm, -1);
            var Blk2 = OrgMtx._ColPartitioning.GetSubBlocking(ColumnIndiceSource, destComm, -1);

            long[] Tlist1 = default(long[]);
            long[] Tlist2 = default(long[]);


            BlockMsrMatrix SUB = new BlockMsrMatrix(Blk1, Blk2);
            OrgMtx.AccSubMatrixTo(1.0, SUB, RowIndicesSource, Tlist1, ColumnIndiceSource, Tlist2);

            return SUB;
        }


        /// <summary>
        /// converts an arbitrary mutable matrix to an <see cref="BlockMsrMatrix"/>.
        /// </summary>
        static public BlockMsrMatrix ToBlockMsrMatrix(this IMutableMatrixEx M, IBlockPartitioning rowmap, IBlockPartitioning colmap) {
            using (new FuncTrace()) {

                if (!rowmap.EqualsPartition(M.RowPartitioning))
                    throw new ArgumentException();
                if (!colmap.EqualsPartition(M.ColPartition))
                    throw new ArgumentException();

                BlockMsrMatrix R = new BlockMsrMatrix(rowmap, colmap);

                long[] col = null;
                double[] val = null;
                int i0 = (int)R.RowPartitioning.i0, L = R.RowPartitioning.LocalLength;
                for (int i = 0; i < L; i++) {
                    int iRow = i0 + i;

                    int Lr = M.GetRow(iRow, ref col, ref val);
                    //R.SetRow(iRow, col, val, Lr); 
                    for( int l = 0; l < Lr; l++) {
                        R[iRow, col[l]] = val[l];
                    }
                }

                return R;
            }
        }

        /// <summary>
        /// Adds <paramref name="factor"/> to all diagonal entries of <paramref name="M"/>.
        /// </summary>
        static public void AccEyeSp(this BlockMsrMatrix M, double factor = 1.0) {
            using (new FuncTrace()) {
                if (M._RowPartitioning.LocalLength != M._ColPartitioning.LocalLength)
                    throw new ArgumentException("supported only for quadratic matrices");
                if (M._RowPartitioning.LocalNoOfBlocks != M._ColPartitioning.LocalNoOfBlocks)
                    throw new ArgumentException("supported only for quadratic matrices");
                               

                long[] RowIdx = M._RowPartitioning.GetOccupiedIndicesList();
                long[] ColIdx = M._ColPartitioning.GetOccupiedIndicesList();

                int ir = 0, ic = 0;
                while(ir < RowIdx.Length && ic < ColIdx.Length) {
                    if (RowIdx[ir] == ColIdx[ic]) {
                        M[RowIdx[ir], ColIdx[ic]] += factor;
                        ir++;
                        ic++;
                    } else if (RowIdx[ir] < ColIdx[ic]) {
                        ir++;
                    } else if (RowIdx[ir] > ColIdx[ic]) {
                        ic++;
                    } else {
                        throw new ApplicationException();
                    }
                }
            }
        }


        /// <summary>
        /// Extraction of a block from a <see cref="BlockMsrMatrix"/>;
        /// 
        /// This is a convenience routine, which
        /// causes memory allocation, and therefore impacts performance slightly negatively.
        /// For best performance, use <see cref="BlockMsrMatrix.ReadBlock"/> instead.
        /// </summary>
        static public MultidimensionalArray GetBlock(this BlockMsrMatrix M, long iBlk, long jBlk) {
            var rpart = M._RowPartitioning;
            var cpart = M._ColPartitioning;
            int _sz_i = rpart.GetBlockLen(iBlk);
            int _sz_j = cpart.GetBlockLen(jBlk);
            long _idx_i = rpart.GetBlockI0(iBlk);
            long _idx_j = cpart.GetBlockI0(jBlk);
            var ret = MultidimensionalArray.Create(_sz_i, _sz_j);
            M.ReadBlock(_idx_i, _idx_j, ret);
            return ret;
        }

        /// <summary>
        /// Setting a block from a <see cref="BlockMsrMatrix"/>;
        /// 
        /// This is a convenience routine, which
        /// causes memory allocation, and therefore impacts performance slightly negatively.
        /// For best performance, use <see cref="BlockMsrMatrix.AccBlock(long, long, double, MultidimensionalArray, double)"/> instead.
        /// </summary>
        public static void SetBlock(this BlockMsrMatrix M, MultidimensionalArray Blk, long iBlk, long jBlk) {
            var rpart = M._RowPartitioning;
            var cpart = M._ColPartitioning;
            int _sz_i = rpart.GetBlockLen(iBlk);
            int _sz_j = cpart.GetBlockLen(jBlk);
            if(_sz_i != Blk.NoOfRows)
                throw new ArgumentException();
            if(_sz_j != Blk.NoOfCols)
                throw new ArgumentException();
            long _idx_i = rpart.GetBlockI0(iBlk);
            long _idx_j = cpart.GetBlockI0(jBlk);

            M.AccBlock(_idx_i, _idx_j, 1.0, Blk, 0.0); // last entry 0.0 will also clear NAN, INF, etc.
        }



    }


    /// <summary>
    /// Class for sparse matrices - it is mutable and used during equation
    /// assembly. When the system of equations is complete, this object can be
    /// handed over to an <see cref="ISparseSolver"/> - implementation that
    /// solves the system. For performance reasons, this solver typically will
    /// convert this matrix into his internal format.<br/>
    /// This matrix is addressed by global row/column indices; 
    /// </summary>
    /// <remarks>
    /// MSR stands for 'M'utuble 'S'parse 'R'ow;
    /// </remarks>
    public class BlockMsrMatrix : IMutableMatrixEx {

        /*
        static void Debug_Assert(bool cond) {
            if(!cond) {
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
                string err = $"BlockMsrMatrix rank{rank} error";

                
                bool o = ilPSP.Environment.StdoutOnlyOnRank0;
                ilPSP.Environment.StdoutOnlyOnRank0 = false;
                Console.Error.WriteLine(err);
                ilPSP.Environment.StdoutOnlyOnRank0 = o;
                

                //throw new ApplicationException(err);

            }
        }

        static void Debug_Assert(bool cond, string message) {
            if(!cond) {
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
                string err = $"BlockMsrMatrix rank{rank} error: {message}";
               
                bool o = ilPSP.Environment.StdoutOnlyOnRank0;
                ilPSP.Environment.StdoutOnlyOnRank0 = false;
                Console.Error.WriteLine(err);
                ilPSP.Environment.StdoutOnlyOnRank0 = o;
                

                //throw new ApplicationException(err);
            }
        }
        */

        /// <summary>
        /// MPI Communicator on which this object lives on
        /// </summary>
        public MPI_Comm MPI_Comm {
            get {
                return m_RowPartitioning.MPI_Comm;
            }
        }

        /// <summary>
        /// <see cref="RowPartitioning"/>;
        /// </summary>
        IBlockPartitioning m_RowPartitioning;

        /// <summary>
        /// distribution of matrix rows over MPI processors
        /// </summary>
        public IBlockPartitioning _RowPartitioning {
            get {
                Debug.Assert(m_RowPartitioning.IsMutable == false);
                return m_RowPartitioning;
            }
        }

        /// <summary>
        /// <see cref="RowPartitioning"/>;
        /// </summary>
        IBlockPartitioning m_ColPartitioning;

        /// <summary>
        /// distribution of matrix rows over MPI processors
        /// </summary>
        public IBlockPartitioning _ColPartitioning {
            get {
                Debug.Assert(m_ColPartitioning.IsMutable == false);
                return m_ColPartitioning;
            }
        }

        /// <summary>
        /// distribution of matrix rows over MPI processors
        /// </summary>
        public IPartitioning RowPartitioning {
            get {
                Debug.Assert(m_RowPartitioning.IsMutable == false);
                return m_RowPartitioning;
            }
        }

        /// <summary>
        /// distribution of matrix rows over MPI processors
        /// </summary>
        public IPartitioning ColPartition {
            get {
                Debug.Assert(m_ColPartitioning.IsMutable == false);
                return m_ColPartitioning;
            }
        }

        /// <summary>
        /// Constructor for a quadratic matrix.
        /// </summary>
        /// <param name="Part">Used for <see cref="_RowPartitioning"/> and <see cref="_ColPartitioning"/>.</param>
        public BlockMsrMatrix(IBlockPartitioning Part) :
            this(Part, Part) {
        }

        /// <summary>
        /// Constructor for a quadratic or rectangular matrix.
        /// </summary>
        /// <param name="rowPart">See <see cref="_RowPartitioning"/>.</param>
        /// <param name="colPart">See <see cref="_ColPartitioning"/>.</param>
        public BlockMsrMatrix(IBlockPartitioning rowPart, IBlockPartitioning colPart) {
            this.m_RowPartitioning = rowPart.IsMutable ? rowPart.GetImmutableBlockPartitioning() : rowPart;
            if(Object.ReferenceEquals(rowPart, colPart)) {
                this.m_ColPartitioning = this.m_RowPartitioning;
            } else {
                this.m_ColPartitioning = colPart.IsMutable ? colPart.GetImmutableBlockPartitioning() : colPart;
            }
            if(m_RowPartitioning.IsMutable)
                throw new ApplicationException();
            if(m_ColPartitioning.IsMutable)
                throw new ApplicationException();

            this.m_BlockRows = new SortedDictionary<long, BlockEntry>[rowPart.LocalNoOfBlocks];
            this.m_ExternalBlock = new BitArray(rowPart.LocalNoOfBlocks);
        }

        /// <summary>
        /// private/empty constructor.
        /// </summary>
        BlockMsrMatrix() { }

        /// <summary>
        /// A utility function to aid the XDG moving interface time-stepper, 
        /// i.e. the application of a row- and column permutation which is local to certain blocks.
        /// The memory of the actual matrix is 'recycled' i.e. unchanged rows are copied shallow. 
        /// After this operation, this matrix becomes unusable.
        /// </summary>
        /// <param name="OldColIndices2New"></param>
        /// <param name="OldRowIndices2New"></param>
        /// <param name="newColPart">Column blocking for the returned matrix.</param>
        /// <param name="newRowPart">Row blocking for the returned matrix.</param>
        /// <returns></returns>
        /// <remarks>
        /// Since this is a very special-purpose function, which is only useful for the XDG moving interface time-stepper, it is questionable
        /// to put it here; on the other hand, we need to access tons of private members to be efficient, 
        /// therefore it is difficult to implement this method in some other place.
        /// </remarks>
        public BlockMsrMatrix RecyclePermute(IBlockPartitioning newRowPart, IBlockPartitioning newColPart,
            long[,] OldRowIndices2New, long[,] OldColIndices2New
            ) {
            using(new FuncTrace()) {
                if(newRowPart.LocalNoOfBlocks != this._RowPartitioning.LocalNoOfBlocks)
                    throw new ArgumentException();
                if(newRowPart.LocalLength != this._RowPartitioning.LocalLength)
                    throw new ArgumentException();
                if(newColPart.LocalNoOfBlocks != this._ColPartitioning.LocalNoOfBlocks)
                    throw new ArgumentException();
                if(newColPart.LocalLength != this._ColPartitioning.LocalLength)
                    throw new ArgumentException();
                if(OldRowIndices2New.GetLength(1) != 2)
                    throw new ArgumentException("Second index is only 0 or 1: 0 for the old, 1 for the new row index of the permutation operation.");
                if(OldColIndices2New.GetLength(1) != 2)
                    throw new ArgumentException("Second index is only 0 or 1: 0 for the old, 1 for the new column index of the permutation operation.");
                int LR = OldRowIndices2New.GetLength(0);
                int LC = OldColIndices2New.GetLength(0);
#if DEBUG_EXTENDED
                for (int l = 1; l < LC; l++) {
                    Debug.Assert(OldColIndices2New[l, 0] > OldColIndices2New[l - 1, 0]);
                }
                for (int l = 1; l < LR; l++) {
                    Debug.Assert(OldRowIndices2New[l, 0] > OldRowIndices2New[l - 1, 0]);
                }

#endif

                this.VerifyDataStructure("Recycle_this");

                BlockMsrMatrix Ret = new BlockMsrMatrix();
                Ret.m_RowPartitioning = newRowPart.IsMutable ? newRowPart.GetImmutableBlockPartitioning() : newRowPart;
                if(Object.ReferenceEquals(newRowPart, newColPart)) {
                    Ret.m_ColPartitioning = Ret.m_RowPartitioning;
                } else {
                    Ret.m_ColPartitioning = newColPart.IsMutable ? newColPart.GetImmutableBlockPartitioning() : newColPart;
                }
                if(Ret.m_RowPartitioning.IsMutable)
                    throw new ApplicationException();
                if(Ret.m_ColPartitioning.IsMutable)
                    throw new ApplicationException();

                if(newRowPart.LocalNoOfBlocks != _RowPartitioning.LocalNoOfBlocks)
                    throw new ArgumentException("Number of blocks must stay the same.");
                if(newColPart.LocalNoOfBlocks != _ColPartitioning.LocalNoOfBlocks)
                    throw new ArgumentException("Number of blocks must stay the same.");

                Ret.m_AssumeSymmetric = this.m_AssumeSymmetric;
                Ret.m_ExternalBlock = this.m_ExternalBlock;
                Ret.m_ExternalBlockIndicesByProcessor = this.m_ExternalBlockIndicesByProcessor;
                Ret.m_Membanks = this.m_Membanks;
                Ret.m_BlockRows = this.m_BlockRows;

#if DEBUG_EXTENDED
                if (!this.m_ExternalBlockIndicesByProcessor.Any(kv => kv.Value.Count > 0)) {
                    // no external block
                    Debug.Assert(this.ComPatternValid);
                }
#endif
                Ret.ComPatternValid = !this.m_ExternalBlockIndicesByProcessor.Any(kv => kv.Value.Count > 0);

                // Markers for all block rows and columns which we are going to permute
                // ====================================================================

                long FirstRowBlock = this._RowPartitioning.FirstBlock;
                long FirstColBlock = this._ColPartitioning.FirstBlock;

                BitArray RowBlockMarker = new BitArray(this._RowPartitioning.LocalNoOfBlocks);
                BitArray ColBlockMarker = new BitArray(this._ColPartitioning.LocalNoOfBlocks);
                HashSet<long> ExtColBlockMarker = new HashSet<long>();

                for(int l = 0; l < LR; l++) { // row indices loop
                    long iOld = OldRowIndices2New[l, 0];
                    long iNew = OldRowIndices2New[l, 1];
                    _RowPartitioning.TestIfInLocalRange(iOld);
                    newRowPart.TestIfInLocalRange(iNew);

                    long iBlock = _RowPartitioning.GetBlockIndex(iOld);
                    if(_RowPartitioning.GetBlockIndex(iOld) != newRowPart.GetBlockIndex(iNew))
                        throw new ArgumentException("Recycling only supports resorting within blocks.");

                    RowBlockMarker[(int)(iBlock - FirstRowBlock)] = true;
                }


                for(int l = 0; l < LC; l++) { // column indices loop
                    long jOld = OldColIndices2New[l, 0];
                    long jNew = OldColIndices2New[l, 1];
                    if(_ColPartitioning.IsInLocalRange(jOld)) {
                        newColPart.TestIfInLocalRange(jNew);
                        long iBlock = _ColPartitioning.GetBlockIndex(jOld);
                        if(_ColPartitioning.GetBlockIndex(jOld) != newColPart.GetBlockIndex(jNew))
                            throw new ArgumentException("Recycling only supports resorting within blocks.");

                        ColBlockMarker[(int)(iBlock - FirstColBlock)] = true;
                    } else {
                        if(newColPart.IsInLocalRange(jNew)) {
                            throw new ArgumentException("If old column index is external, new column index should also be.");
                        }

                        long jBlock, j0, jE;
                        GetExternalSubblockIndices(_ColPartitioning, jOld, out jBlock, out j0, out jE);
                        ExtColBlockMarker.Add(jBlock);
                    }
                }




                // loop over block rows
                // ====================

                int NoofRowBlocks = _RowPartitioning.LocalNoOfBlocks;
                Debug.Assert(m_BlockRows.Length == NoofRowBlocks);
                int lr = 0;
                for(int iBlockLoc = 0; iBlockLoc < NoofRowBlocks; iBlockLoc++) {
                    //this.VerifyDataStructure();
                    //Ret.VerifyDataStructure();
                    var Row = m_BlockRows[iBlockLoc];
                    if(Row != null) {
                        bool RowPermReq = RowBlockMarker[iBlockLoc]; // row permutation is required.

                        long i0 = _RowPartitioning.GetBlockI0(iBlockLoc + FirstRowBlock);
                        long iE = _RowPartitioning.GetBlockLen(iBlockLoc + FirstRowBlock) + i0;

                        int llr0, llre, LRinc;
                        if(RowPermReq) {
                            LRinc = Scan(lr, OldRowIndices2New, _RowPartitioning, false, iBlockLoc + FirstRowBlock, out llr0, out llre);
                            Debug.Assert(llre - llr0 >= 0);
                            Debug.Assert(LRinc > 0);
#if DEBUG_EXTENDED
                            for (int llr = llr0; llr <= llre; llr++) {
                                Debug.Assert(OldRowIndices2New[llr, 0] >= i0);
                                Debug.Assert(OldRowIndices2New[llr, 0] < iE);
                            }
#endif
                        } else {
                            LRinc = 0;
                            llre = int.MinValue + 10;
                            llr0 = int.MinValue;
                        }


                        long LastjBlock = -1;
                        int lc = 0;
                        var OrgRow = Row.ToArray();
                        foreach(var kv in OrgRow) {
                            long jBlock = kv.Key;
                            Debug.Assert(jBlock == kv.Value.jBlkCol);
                            Debug.Assert(jBlock > LastjBlock);
                            LastjBlock = jBlock;

                            long j0, jE;
                            bool colPermReq;
                            if(_ColPartitioning.IsLocalBlock(jBlock)) {
                                colPermReq = ColBlockMarker[(int)(jBlock - FirstColBlock)];
                                j0 = _ColPartitioning.GetBlockI0(jBlock);
                                jE = _ColPartitioning.GetBlockLen(jBlock) + j0;
                            } else {
                                colPermReq = ExtColBlockMarker.Contains(jBlock);
                                GetExternalSubblockIndices(_ColPartitioning, jBlock, out j0, out jE);
                                Debug.Assert(jE - j0 > 0);
                            }

                            int llc0, llce, LCinc;
                            if(colPermReq) {
                                LCinc = Scan(lc, OldColIndices2New, _ColPartitioning, true, jBlock, out llc0, out llce);
                                Debug.Assert(llce - llc0 >= 0);
                                Debug.Assert(LCinc > 0);
#if DEBUG_EXTENDED
                                for (int llc = llc0; llc <= llce; llc++) {
                                    Debug.Assert(OldColIndices2New[llc, 0] >= j0);
                                    Debug.Assert(OldColIndices2New[llc, 0] < jE);
                                }
#endif
                            } else {
                                LCinc = 0;
                                llce = int.MinValue + 10;
                                llc0 = int.MinValue;
                            }


                            if(RowPermReq || colPermReq) {
                                int I = (int)(iE - i0);
                                int J = (int)(jE - j0);

                                MultidimensionalArray tmpBlock = MultidimensionalArray.Create(I, J);
                                MultidimensionalArray prmBlock = MultidimensionalArray.Create(I, J);
                                this.ReadBlock(i0, j0, tmpBlock);
                                this.ClearBlock(i0, j0, I, J);
#if DEBUG_EXTENDED
                                {
                                    Debug.Assert(this.m_BlockRows[iBlockLoc].ContainsKey(jBlock) == false);
                                    Debug.Assert(Ret.m_BlockRows[iBlockLoc].ContainsKey(jBlock) == false);

                                    int RowBlockType = this._RowPartitioning.GetBlockType(iBlockLoc + FirstRowBlock);
                                    BlockEntry BE = kv.Value;
                                    Debug.Assert(BE.InMembnk.GetLength(0) == this._RowPartitioning.GetSubblk_i0(RowBlockType).Length);
                                    if (this._ColPartitioning.IsLocalBlock(jBlock)) {
                                        int ColBlockType = this._ColPartitioning.GetBlockType(iBlockLoc + FirstColBlock);
                                        Debug.Assert(BE.InMembnk.GetLength(1) == this._ColPartitioning.GetSubblk_i0(ColBlockType).Length);
                                    }
                                    Debug.Assert(BE.InMembnk.GetLength(0) == BE.MembnkIdx.GetLength(0));
                                    Debug.Assert(BE.InMembnk.GetLength(1) == BE.MembnkIdx.GetLength(1));
                                    for (int __i = 0; __i < BE.InMembnk.GetLength(0); __i++) {
                                        for (int __j = 0; __j < BE.InMembnk.GetLength(1); __j++) {
                                            Debug.Assert(kv.Value.InMembnk[__i, __j] < 0);
                                            Debug.Assert(kv.Value.MembnkIdx[__i, __j] < 0);
                                        }
                                    }
                                }
#endif                                

                                if(RowPermReq) {
                                    for(int llr = llr0; llr <= llre; llr++) {
                                        int iOrg = (int)(OldRowIndices2New[llr, 0] - i0);
                                        int iDst = (int)(OldRowIndices2New[llr, 1] - i0);
                                        Debug.Assert(iOrg >= 0);
                                        Debug.Assert(iDst >= 0);
                                        Debug.Assert(iOrg < I);
                                        Debug.Assert(iDst < I);
                                        for(int _j = 0; _j < J; _j++)
                                            prmBlock[iDst, _j] = tmpBlock[iOrg, _j];
                                    }

                                    if(colPermReq) {
                                        var kkkk = tmpBlock;
                                        tmpBlock = prmBlock;
                                        prmBlock = kkkk;
                                        prmBlock.Clear();
                                    }
                                }

                                if(colPermReq) {
                                    for(int llc = llc0; llc <= llce; llc++) {
                                        int jOrg = (int)(OldColIndices2New[llc, 0] - j0);
                                        int jDst = (int)(OldColIndices2New[llc, 1] - j0);
                                        Debug.Assert(jOrg >= 0);
                                        Debug.Assert(jDst >= 0);
                                        Debug.Assert(jOrg < J);
                                        Debug.Assert(jDst < J);
                                        for(int _i = 0; _i < I; _i++)
                                            prmBlock[_i, jDst] = tmpBlock[_i, jOrg];
                                    }
                                }

                                //this.VerifyDataStructure();
                                //Ret.VerifyDataStructure();
                                Ret.AccBlock(i0, j0, 1.0, prmBlock);
#if DEBUG_EXTENDED
                                {
                                    Debug.Assert(Ret.m_BlockRows[iBlockLoc].ContainsKey(jBlock) == true);

                                    BlockEntry BE = Ret.m_BlockRows[iBlockLoc][jBlock];
                                    int RowBlockType = Ret._RowPartitioning.GetBlockType(iBlockLoc + FirstRowBlock);
                                    Debug.Assert(BE.InMembnk.GetLength(0) == Ret._RowPartitioning.GetSubblk_i0(RowBlockType).Length);
                                    if (Ret._ColPartitioning.IsLocalBlock(jBlock)) {
                                        int ColBlockType = Ret._ColPartitioning.GetBlockType(iBlockLoc + FirstColBlock);
                                        Debug.Assert(BE.InMembnk.GetLength(1) == Ret._ColPartitioning.GetSubblk_i0(ColBlockType).Length);
                                    }
                                    Debug.Assert(BE.InMembnk.GetLength(0) == BE.MembnkIdx.GetLength(0));
                                    Debug.Assert(BE.InMembnk.GetLength(1) == BE.MembnkIdx.GetLength(1));
                                    for (int __i = 0; __i < BE.InMembnk.GetLength(0); __i++) {
                                        for (int __j = 0; __j < BE.InMembnk.GetLength(1); __j++) {
                                            Debug.Assert(kv.Value.InMembnk[__i, __j] < 0);
                                            Debug.Assert(kv.Value.MembnkIdx[__i, __j] < 0);
                                        }
                                    }
                                }
#endif
                            }

                            lc += LCinc;
                        }

                        lr += LRinc;
                    }
                }

                // invalidate this matrix and return
                // =================================

                m_ExternalBlock = null;
                m_ExternalBlockIndicesByProcessor = null;
                m_Membanks = null;
                m_BlockRows = null;
                m_RowPartitioning = null;
                m_ColPartitioning = null;

                Ret.VerifyDataStructure("Recycle_ret");

                return Ret;
            }
        }

        static int Scan(int l, long[,] old2new, IBlockPartitioning part, bool ExtOk, long iBlock, out int l0, out int le) {
            int Inc = 0;
            l0 = l;
            le = l - 1;
            int L = old2new.GetLength(0);

            if(l >= L)
                return Inc;

            // skip all entries of 'old2new' which are before 'iBlock'
            // -------------------------------------------------------
            long _iBlock, i;
            while(true) {
                i = old2new[l, 0];
                if(part.IsInLocalRange(i)) {
                    _iBlock = part.GetBlockIndex(i);
                } else {
                    Debug.Assert(ExtOk == true);
                    long d1, d2;
                    GetExternalSubblockIndices(part, i, out _iBlock, out d1, out d2);
                    Debug.Assert(d1 <= i);
                    Debug.Assert(i <= d2);
                }

                if(_iBlock < iBlock) {
                    l0++;
                    le++;
                    l++;
                    if(l >= L)
                        return Inc;
                } else {
                    break;
                }
            }

            // scan all entries of 'old2new' which in 'iBlock'
            // -----------------------------------------------
            Debug.Assert(l >= l0);
            Debug.Assert(le < l0);
            Debug.Assert(_iBlock >= iBlock);
            while(iBlock == _iBlock) {
                le++;
                Inc++;
                l++;

                if(l >= L)
                    return Inc;

                Debug.Assert(old2new[l, 0] > i, "Indices must be strictly ascending.");
                i = old2new[l, 0];
                if(part.IsInLocalRange(i)) {
                    _iBlock = part.GetBlockIndex(i);
                } else {
                    long Block_i0, Block_iE;
                    GetExternalSubblockIndices(part, i, out _iBlock, out Block_i0, out Block_iE);
                }
                Debug.Assert(_iBlock >= iBlock);
            }

            return Inc;
        }



        /// <summary>
        /// - array index: local block row
        /// - dictionary key: block column
        /// - dictionary value: a block entry, <see cref="BlockEntry.jBlkCol"/> equals the key.
        /// </summary>
        SortedDictionary<long, BlockEntry>[] m_BlockRows;

        /// <summary>
        /// Marker for rows with blocks/entries in the external range of the <see cref="_ColPartitioning"/>.
        /// - index: local block row.
        /// - value: true if the row contains an external block.
        /// </summary>
        BitArray m_ExternalBlock;

        /// <summary>
        /// Tracking of external blocks, i.e. blocks outside of the local range of the column partitioning.
        /// - key: MPI processor rank (within <see cref="MPI_Comm"/>)
        /// - value: set of block indices which are owned by the respective processor.
        /// </summary>
        Dictionary<int, HashSet<long>> m_ExternalBlockIndicesByProcessor = new Dictionary<int, HashSet<long>>();

        /// <summary>
        /// Flag which indicates that <see cref="ReceiveLists"/> resp. <see cref="SendLists"/> contains valid data.
        /// </summary>
        bool ComPatternValid = true;

        /// <summary>
        /// - key: MPI processor rank within the <see cref="MPI_Comm"/> communicator
        /// - value: index ranges which are received from respective processor.
        /// - 0th index into value:
        /// - 1st index into value: either 0 (start of index range) or 1 (end of index range)
        /// </summary>
        Dictionary<int, long[,]> ReceiveLists = new Dictionary<int, long[,]>();

        /// <summary>
        /// - key: MPI processor rank within the <see cref="MPI_Comm"/> communicator
        /// - value: index ranges which have to be sent to other processor
        /// </summary>
        Dictionary<int, long[,]> SendLists = new Dictionary<int, long[,]>();

        /// <summary>
        /// Update of <see cref="ReceiveLists"/> and <see cref="SendLists"/>
        /// </summary>
        void UpdateCommPattern(MPI_Comm comm) {
            MPICollectiveWatchDog.Watch(comm);

            if(comm != this.MPI_Comm)
                throw new NotSupportedException();
            // hint: MPI_Comm_group & MPI_Group_translate_ranks


            // check if the comm pattern is still valid on **all** Processors in the communicator
            // ----------------------------------------------------------------------------------

            //Console.WriteLine("P{0}: comm pattern valid: {1}.", m_RowPartitioning.MpiRank, ComPatternValid);
            int _ComPatternValid = ComPatternValid ? 11 : 0;
            _ComPatternValid = _ComPatternValid.MPIMin(this.MPI_Comm); // potentially super-expensive!!!!
            if(_ComPatternValid != 0)
                return;


            // clear old Lists
            // ---------------

            int MpiRank;
            int MpiSize;
            csMPI.Raw.Comm_Rank(comm, out MpiRank);
            csMPI.Raw.Comm_Rank(comm, out MpiSize);

            if(ReceiveLists != null)
                ReceiveLists.Clear();
            else
                ReceiveLists = new Dictionary<int, long[,]>();
            if(SendLists != null)
                SendLists.Clear();
            else
                SendLists = new Dictionary<int, long[,]>();

            int[] PeerProcs = m_ExternalBlockIndicesByProcessor.Keys.ToArray();

            // build receive lists
            // -------------------

            foreach(int ProcRank in PeerProcs) {
                Debug.Assert(ProcRank != MpiRank);
                long[,] IdxRange = GetExternalIndexRange(ProcRank);
                ReceiveLists.Add(ProcRank, IdxRange); //] = GetExternalIndexRange(ProcRank);
            }

            //foreach (var kv in ReceiveLists) {
            //    int Trnk = kv.Key;
            //    int NoOf = kv.Value.GetLength(0);
            //    Console.WriteLine("P{0} from P{1}: #itms {2}", m_RowPartitioning.MpiRank, Trnk, NoOf);
            //}
            //if (ReceiveLists.Count == 0)
            //    Console.WriteLine("P{0}: nothing to receive.", m_RowPartitioning.MpiRank);



            // MPI exchange
            // ------------

            IntBufferExchange(comm, ReceiveLists, out SendLists);


            // check
            // -----

#if DEBUG_EXTENDED
            foreach (long[,] slist in SendLists.Values) {
                int L = slist.GetLength(0);
                Debug.Assert(slist.GetLength(1) == 2);
                for (int l = 0; l < L; l++) {
                    long i0 = slist[l, 0];
                    long iE = slist[l, 1];
                    Debug.Assert(iE - i0 > 0);
                    Debug.Assert(m_ColPartitioning.IsInLocalRange(i0));
                    Debug.Assert(m_ColPartitioning.IsInLocalRange(iE - 1));
                }
            }
#endif
        }


        static void IntBufferExchange(MPI_Comm comm, Dictionary<int, long[,]> Sdata, out Dictionary<int, long[,]> Rdata) {
            int MpiRank;
            int MpiSize;
            csMPI.Raw.Comm_Rank(comm, out MpiRank);
            csMPI.Raw.Comm_Size(comm, out MpiSize);

            var PeerProcs = Sdata.Keys;


            unsafe {
                int[] NoOfItems_this = new int[MpiSize];
                int[,] NoOfItems = new int[MpiSize, MpiSize];

                foreach(var kv in Sdata) {
                    int ProcRank = kv.Key;
                    long[,] data = kv.Value;
                    if(data.GetLength(1) != 2)
                        throw new NotSupportedException();
                    Debug.Assert(ProcRank != MpiRank);
                    Debug.Assert(ProcRank >= 0 && ProcRank < MpiSize);
                    NoOfItems_this[ProcRank] = data.GetLength(0) * 2;
                }


                fixed(int* pNoOfItems = NoOfItems, pNoOfItems_this = NoOfItems_this) {
                    csMPI.Raw.Allgather((IntPtr)pNoOfItems_this, MpiSize, csMPI.Raw._DATATYPE.INT, (IntPtr)pNoOfItems, MpiSize, csMPI.Raw._DATATYPE.INT, comm);
                }
#if DEBUG_EXTENDED
                for (int i = 0; i < MpiSize; i++) {
                    Debug.Assert(NoOfItems[MpiRank, i] == NoOfItems_this[i]);
                    Debug.Assert(NoOfItems[i, i] == 0);
                }
#endif

                Rdata = new Dictionary<int, long[,]>();

                List<int> RcvRanks = new List<int>();
                List<long[,]> RcvBuffers = new List<long[,]>();
                List<Tuple<IntPtr, GCHandle>> RcvBuffersUnsafe = new List<Tuple<IntPtr, GCHandle>>();
                for(int i = 0; i < MpiSize; i++) {
                    int NoOfItm = NoOfItems[i, MpiRank];
                    if(NoOfItm > 0) {
                        RcvRanks.Add(i);
                        long[,] buf = new long[NoOfItm / 2, 2];
                        RcvBuffers.Add(buf);
                        GCHandle gch = GCHandle.Alloc(buf, GCHandleType.Pinned);
                        IntPtr pbuf = Marshal.UnsafeAddrOfPinnedArrayElement(buf, 0);
                        RcvBuffersUnsafe.Add(new Tuple<IntPtr, GCHandle>(pbuf, gch));
                        Rdata.Add(i, buf);
                    }
                }

                // all MPI requests for send *and* receive
                // - indices 0 to Sdata.Count -1 are used for send
                // - the other ones for receiving
                MPI_Request[] SendRcvReq = new MPI_Request[RcvBuffers.Count + Sdata.Count];

                for(int i = 0; i < RcvRanks.Count; i++) {
                    int SourceRank = RcvRanks[i];
                    csMPI.Raw.Irecv(RcvBuffersUnsafe[i].Item1, NoOfItems[SourceRank, MpiRank], csMPI.Raw._DATATYPE.LONG_LONG, SourceRank, 213678 + SourceRank * 13, comm, out SendRcvReq[Sdata.Count + i]);
                }

                GCHandle[] SndBuffersUnsafe = new GCHandle[Sdata.Count];
                int ii = 0;
                foreach(var kv in Sdata) {
                    int TargProc = kv.Key;
                    long[,] buf = kv.Value;
                    SndBuffersUnsafe[ii] = GCHandle.Alloc(buf, GCHandleType.Pinned);
                    IntPtr pbuf = Marshal.UnsafeAddrOfPinnedArrayElement(buf, 0);
                    csMPI.Raw.Issend(pbuf, buf.Length, csMPI.Raw._DATATYPE.LONG_LONG, TargProc, 213678 + MpiRank * 13, comm, out SendRcvReq[ii]);
                    ii++;
                }

                MPI_Status[] stats = new MPI_Status[SendRcvReq.Length];
                csMPI.Raw.Waitall(SendRcvReq.Length, SendRcvReq, stats);

                for(int i = RcvBuffersUnsafe.Count - 1; i >= 0; i--) {
                    RcvBuffersUnsafe[i].Item2.Free();
                }
                for(int i = SndBuffersUnsafe.Length - 1; i >= 0; i--) {
                    SndBuffersUnsafe[i].Free();
                }

            }
        }


        long[,] GetExternalIndexRange(int ProcRank) {
            //m_ExternalBlockIndicesByProcessor

            var _set = m_ExternalBlockIndicesByProcessor[ProcRank];
            int L = _set.Count;
            int[] set = new int[L];
            int l = 0;
            foreach(int jBlk in _set) {
                set[l] = jBlk;
                l++;
            }
            Array.Sort(set);

            long[,] Range = new long[L, 2];
            for(int i = 0; i < L; i++) {
                long i0, iE;
                GetExternalSubblockIndices(m_ColPartitioning, set[i], out i0, out iE);
                Debug.Assert(iE - i0 > 0);
                Range[i, 0] = i0;
                Range[i, 1] = iE;

                Debug.Assert(m_ColPartitioning.FindProcess(i0) == ProcRank);
                Debug.Assert(iE - i0 == 0 || m_ColPartitioning.FindProcess(iE - 1) == ProcRank);
            }

            return Range;
        }


        /// <summary>
        /// Pointers/indices for a matrix block into the memory bank.
        /// </summary>
        class BlockEntry : ICloneable {

            /// <summary>
            /// Block-column index.
            /// </summary>
            public long jBlkCol;

            /// <summary>
            /// Memory bank index.
            /// - 1st index: sub-block row 
            /// - 2nd index: sub-block column
            /// - content: index into <see cref="m_Membanks"/>.
            /// </summary>
            public int[,] MembnkIdx;

            /// <summary>
            /// Block index within memory bank
            /// - 1st index: sub-block row 
            /// - 2nd index: sub-block column
            /// - content: 1st index into <see cref="Membank.Mem"/>.
            /// </summary>
            public int[,] InMembnk;

            /// <summary>
            /// True, if all entries in <see cref="MembnkIdx"/> and <see cref="InMembnk"/> are negative.
            /// </summary>
            public bool IsEmpty {
                get {
                    Debug.Assert(MembnkIdx.GetLength(0) == InMembnk.GetLength(0));
                    Debug.Assert(MembnkIdx.GetLength(1) == InMembnk.GetLength(1));
                    int NR = MembnkIdx.GetLength(0);
                    int NC = InMembnk.GetLength(1);

                    for(int i = 0; i < NR; i++) {
                        for(int j = 0; j < NC; j++) {
                            Debug.Assert((MembnkIdx[i, j] < 0) == (InMembnk[i, j] < 0));
                            if(MembnkIdx[i, j] >= 0)
                                return false;
                        }
                    }
                    return true;
                }
            }

            /// <summary>
            /// Cloning.
            /// </summary>
            public object Clone() {
                var R = new BlockEntry();
                R.InMembnk = this.InMembnk.CloneAs();
                R.MembnkIdx = this.MembnkIdx.CloneAs();
                R.jBlkCol = this.jBlkCol;
                return R;
            }
        }


        /// <summary>
        /// get/set an entry; Setting a zero element (to a value unequal to zero)
        /// allocates a new entry; Setting an entry to 0.0 releases the
        /// corresponding memory;
        /// </summary>
        /// <param name="i">Global (over all MPI processes) row index.</param>
        /// <param name="j">
        /// global (over all MPI processes) column index
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// Note that setting entries to 0.0 does not releases memory, if memory for the entry was allocated before.
        /// In the block-matrix structure, this would cause to much overhead.
        /// </remarks>
        public double this[long i, long j] {
            set {
                int iSblk, jSblk; //   row/col index within sub-block, which correspond to (i,j)
                int ISblk, JSblk; //   sub-block size
                double[] Storage; //   sub-block memory: where the result should be accumulated
                int Offset, CI, CJ; // offset pointer and i,j cycles into 'Storage'

                // We do not de-allocate if an entry becomes zero,
                // therefore we would have to check the whole sub-block,
                // which could be expensive.
                // Maybe, we add some de-allocation routine for this.

                GetSetAlloc(value != 0.0, i, j, out iSblk, out jSblk, out ISblk, out JSblk, out Storage, out Offset, out CI, out CJ);
                if(Storage == null && value != 0.0)
                    throw new ArgumentException("Can not save non-zero entry in void-region.");

                Debug.Assert(Storage != null || value == 0.0);

                if(Storage != null) {
                    int StorageIdx = Offset + (iSblk) * CI + (jSblk) * CJ;
                    Storage[StorageIdx] = value;
                }
            }
            get {
                int iSblk, jSblk; //   row/col index within sub-block, which correspond to (i,j)
                int ISblk, JSblk; //   sub-block size
                double[] Storage; //   sub-block memory: where the result should be accumulated
                int Offset, CI, CJ; // offset pointer and i,j cycles into 'Storage'

                GetSetAlloc(false, i, j, out iSblk, out jSblk, out ISblk, out JSblk, out Storage, out Offset, out CI, out CJ);

                if(Storage == null) {
                    return 0.0;
                } else {
                    int StorageIdx = Offset + (iSblk) * CI + (jSblk) * CJ;
                    return Storage[StorageIdx];
                }

            }
        }

        void GetSetAlloc(bool bAlloc,
            long i, long j, //                         global row/column index 
            out int iSblk, out int jSblk, //           row/col index within sub-block, which correspond to i,j
            out int ISblk, out int JSblk, //           sub-block size
            out double[] Storage, //                   sub-block memory: where the result should be accumulated
            out int Offset, out int CI, out int CJ  // offset pointer and i,j cycles into 'Storage'
            ) {
            int MembnkIdx, InMembnk, rowSblkIdx, colSblkIdx;
            long iBlk, jBlk;
            GetSetAlloc(bAlloc,
                        i, j,
                        out iBlk, out jBlk,
                        out rowSblkIdx, out colSblkIdx,
                        out iSblk, out jSblk,
                        out ISblk, out JSblk,
                        out Storage,
                        out Offset, out CI, out CJ,
                        out MembnkIdx, out InMembnk
                        );
#if DEBUG_EXTENDED
            Debug.Assert((Storage != null) == (MembnkIdx >= 0 && InMembnk >= 0));
            if (Storage != null) {
                if (!_ColPartitioning.IsLocalBlock(jBlk)) {
                    Debug.Assert(m_ExternalBlock[checked((int)(iBlk - _RowPartitioning.FirstBlock))] == true);
                    Debug.Assert(_ColPartitioning.FindProcessForBlock(jBlk) == _ColPartitioning.FindProcess(j));
                    int OwnerProcRank = _ColPartitioning.FindProcessForBlock(jBlk);
                    Debug.Assert(m_ExternalBlockIndicesByProcessor.ContainsKey(OwnerProcRank));
                    Debug.Assert(m_ExternalBlockIndicesByProcessor[OwnerProcRank].Contains(jBlk));
                }
            }
#endif
        }

        void GetSetAlloc(bool bAlloc,
            long i, long j, //                         global row/column index 
            out long BlkRow, out long BlkCol, //       block row and block column index
            out int rowSblkIdx, out int colSblkIdx, // sub-block row and sub-block column index
            out int iSblk, out int jSblk, //           row/col index within sub-block, which correspond to i,j
            out int ISblk, out int JSblk, //           sub-block size
            out double[] Storage, //                   sub-block memory: where the result should be accumulated
            out int Offset, out int CI, out int CJ, // offset pointer and i,j cycles into 'Storage'
            out int MembnkIdx, out int InMembnk //     membank index and index within membank
            ) {
            // default values if not allocated, and no allocation requested (bAlloc == false).
            Storage = null;
            Offset = int.MinValue;
            CI = int.MinValue;
            CJ = int.MinValue;

            // =====================
            // translate row indices
            // =====================
            int iBlkLoc, iLoc, iBlkT, i0_Sblk, NoOfRowSblk, RemVoidLen;
            //int[] Sblk_i0, iSblkLen;
            TranslateIndex(i, m_RowPartitioning, false, out BlkRow, out iBlkLoc, out long i0, out iLoc, out iBlkT, out rowSblkIdx, out i0_Sblk, out ISblk, out NoOfRowSblk, out RemVoidLen);
            if(i0_Sblk < 0) {
                // in void row region: the remaining columns of the block are empty
                jSblk = 0;
                long __jColBlk, __Block0, __BlockEnd;
                if(m_ColPartitioning.IsInLocalRange(j)) {
                    __jColBlk = m_ColPartitioning.GetBlockIndex(j);
                    __Block0 = m_ColPartitioning.GetBlockI0(__jColBlk);
                    __BlockEnd = m_ColPartitioning.GetBlockLen(__jColBlk) + __Block0;
                } else {
                    GetExternalSubblockIndices(m_ColPartitioning, j, out __jColBlk, out __Block0, out __BlockEnd);
                }
                Debug.Assert(j >= __Block0 && j < __BlockEnd);
                JSblk = (int)(__BlockEnd - j);
                iSblk = 0;
                ISblk = RemVoidLen;
                MembnkIdx = -3;
                InMembnk = -4;
                BlkCol = -12323;
                colSblkIdx = -676;
                return;
                //}
            }
            iSblk = (int)(i - i0 - i0_Sblk);


            // =========
            // mem alloc
            // =========
            var BlockRows = this.m_BlockRows[iBlkLoc];
            if(BlockRows == null) {
                if(!bAlloc) {
                    jSblk = -23672;
                    JSblk = -2178326;
                    MembnkIdx = -3;
                    InMembnk = -4;
                    BlkCol = -12323;
                    colSblkIdx = -676;
                    return;
                } else {
                    this.m_BlockRows[iBlkLoc] = new SortedDictionary<long, BlockEntry>();
                    BlockRows = this.m_BlockRows[iBlkLoc];
                }
            }

            {

                // =====================
                // translate col indices
                // =====================

                int jBlkLoc, jLoc, jBlkT, j0_Sblk, NoOfColSblk, RemVoidCols;
                TranslateIndex(j, m_ColPartitioning, true, out BlkCol, out jBlkLoc, out long j0, out jLoc, out jBlkT, out colSblkIdx, out j0_Sblk, out JSblk, out NoOfColSblk, out RemVoidCols);
                if(j0_Sblk < 0) {
                    // void region
                    jSblk = 0;
                    JSblk = RemVoidCols;
                    MembnkIdx = -3;
                    InMembnk = -4;
                    return;
                }

                jSblk = (int)(j - j0 - j0_Sblk);
                Debug.Assert(jSblk >= 0);
                Debug.Assert(jSblk < JSblk);
                Debug.Assert(_ColPartitioning.IsInLocalRange(j) == _ColPartitioning.IsLocalBlock(BlkCol));
                Debug.Assert(_ColPartitioning.IsInLocalRange(j) == true || NoOfColSblk == 1);

                // ============
                // access block
                // ============
                BlockEntry BE;
                if(!BlockRows.TryGetValue(BlkCol, out BE)) {
                    if(!bAlloc) {
                        //jSblk = int.MinValue;
                        //JSblk = int.MinValue;
                        MembnkIdx = -3;
                        InMembnk = -4;
                        return;
                    } else {
                        // mem alloc necessary 
                        BE = new BlockEntry();
                        BE.jBlkCol = BlkCol;
                        BE.InMembnk = new int[NoOfRowSblk, NoOfColSblk];
                        BE.MembnkIdx = new int[NoOfRowSblk, NoOfColSblk];
                        ArrayTools.SetAll(BE.InMembnk, int.MinValue);
                        ArrayTools.SetAll(BE.MembnkIdx, int.MinValue);
                        BlockRows.Add(BlkCol, BE);
                    }
                }
                Debug.Assert(BE.InMembnk.GetLength(0) == NoOfRowSblk);
                Debug.Assert(BE.InMembnk.GetLength(1) == NoOfColSblk);
                Debug.Assert(BE.MembnkIdx.GetLength(0) == NoOfRowSblk);
                Debug.Assert(BE.MembnkIdx.GetLength(1) == NoOfColSblk);

                // ================
                // access sub-block
                // ================

                // 
                InMembnk = BE.InMembnk[rowSblkIdx, colSblkIdx];
                MembnkIdx = BE.MembnkIdx[rowSblkIdx, colSblkIdx];
                Membank B;
                if(MembnkIdx < 0) {
                    // a void entry/position
                    if(!bAlloc) {
                        //jSblk = int.MinValue;
                        //JSblk = int.MinValue;
                        return;
                    } else {
                        // mem alloc necessary 
                        AllocSblk(out B, out MembnkIdx, out InMembnk, ISblk, JSblk);
                        BE.InMembnk[rowSblkIdx, colSblkIdx] = InMembnk;
                        BE.MembnkIdx[rowSblkIdx, colSblkIdx] = MembnkIdx;


                        bool isExternal = !_ColPartitioning.IsLocalBlock(BlkCol);// < m_ColPartitioning.FirstBlock || BlkCol >= m_ColPartitioning.FirstBlock + m_ColPartitioning.LocalNoOfBlocks;
                        if(isExternal) {
                            // +++++++++++++++++++++++++++++++++++++
                            // allocated an *external* block 
                            // this *may* invalidate the communication pattern.
                            // +++++++++++++++++++++++++++++++++++++

                            // block range
                            int OwnerRank = m_ColPartitioning.FindProcess(j);
                            Debug.Assert(OwnerRank != m_ColPartitioning.MpiRank);
                            HashSet<long> BlockIndicesByProcessor;
                            if(!m_ExternalBlockIndicesByProcessor.TryGetValue(OwnerRank, out BlockIndicesByProcessor)) {
                                BlockIndicesByProcessor = new HashSet<long>();
                                m_ExternalBlockIndicesByProcessor.Add(OwnerRank, BlockIndicesByProcessor);
                            }
                            this.ComPatternValid &= !BlockIndicesByProcessor.Add(BlkCol);
                            this.m_ExternalBlock[iBlkLoc] |= true;
                        }

                    }
                } else {
                    B = m_Membanks[MembnkIdx];
                }

                // get access pattern
                bool isDense;

                B.GetFastBlockAccessInfo(out Storage, out Offset, out CI, out CJ, out isDense, InMembnk);

            }
        }


        void GetSetAlloc2(bool bAlloc,
            long iBlock, long jBlock, //              row/col block index
            int rowSblk, int colSblk, //              row/col index within sub-block, which correspond to i,j
            int ISblk, int JSblk, //                  sub-block size
            int NoOfRowSblk, int NoOfColSblk, //      number of sub-blocks
            out double[] Storage, //                  sub-block memory: where the result should be accumulated
            out int Offset, out int CI, out int CJ // offset pointer and i,j cycles into 'Storage'
            ) {

            Debug.Assert(jBlock >= 0);
            Debug.Assert(jBlock < this._ColPartitioning.TotalNoOfBlocks);
            int iBlockLoc = (int)(iBlock - this._RowPartitioning.FirstBlock);
            Debug.Assert(iBlockLoc >= 0);
            Debug.Assert(iBlockLoc < this._RowPartitioning.LocalNoOfBlocks);


#if DEBUG_EXTENDED
            {
                int rowBT = this._RowPartitioning.GetBlockType(iBlock);
                int[] i0_Sblk_row = this._RowPartitioning.GetSubblk_i0(rowBT);
                int[] LenSblk_row = this._RowPartitioning.GetSubblkLen(rowBT);
                Debug.Assert(i0_Sblk_row.Length == NoOfRowSblk);
                Debug.Assert(LenSblk_row.Length == NoOfRowSblk);
                Debug.Assert(ISblk == LenSblk_row[rowSblk]);

                if (this._ColPartitioning .IsLocalBlock(jBlock)) {
                    int colBT = this._ColPartitioning.GetBlockType(jBlock);
                    int[] i0_Sblk_col = this._ColPartitioning.GetSubblk_i0(colBT);
                    int[] LenSblk_col = this._ColPartitioning.GetSubblkLen(colBT);
                    Debug.Assert(i0_Sblk_col.Length == NoOfColSblk);
                    Debug.Assert(LenSblk_col.Length == NoOfColSblk);
                    Debug.Assert(JSblk == LenSblk_col[colSblk]);
                } else {
                    Debug.Assert(NoOfColSblk == 1);
                }
            }
#endif

            var BlockRow = this.m_BlockRows[iBlockLoc];
            if(BlockRow == null) {
                if(bAlloc) {
                    BlockRow = new SortedDictionary<long, BlockEntry>();
                    this.m_BlockRows[iBlockLoc] = BlockRow;
                } else {
                    Storage = null;
                    Offset = -12364;
                    CI = -12349;
                    CJ = -66;
                    return;
                }
            }


            BlockEntry BE;
            if(!BlockRow.TryGetValue(jBlock, out BE)) {
                if(!bAlloc) {
                    Storage = null;
                    Offset = -12364;
                    CI = -12349;
                    CJ = -666;
                    return;
                } else {
                    // mem alloc necessary 
                    BE = new BlockEntry();
                    BE.jBlkCol = jBlock;
                    BE.InMembnk = new int[NoOfRowSblk, NoOfColSblk];
                    BE.MembnkIdx = new int[NoOfRowSblk, NoOfColSblk];
                    ArrayTools.SetAll(BE.InMembnk, int.MinValue);
                    ArrayTools.SetAll(BE.MembnkIdx, int.MinValue);
                    BlockRow.Add(jBlock, BE);
                }
            }
            Debug.Assert(BE.InMembnk.GetLength(0) == NoOfRowSblk);
            Debug.Assert(BE.InMembnk.GetLength(1) == NoOfColSblk);
            Debug.Assert(BE.MembnkIdx.GetLength(0) == NoOfRowSblk);
            Debug.Assert(BE.MembnkIdx.GetLength(1) == NoOfColSblk);

            // ================
            // access sub-block
            // ================

            // 
            int InMembnk = BE.InMembnk[rowSblk, colSblk];
            int MembnkIdx = BE.MembnkIdx[rowSblk, colSblk];
            Membank B;
            if(MembnkIdx < 0) {
                // a void entry/position
                if(!bAlloc) {
                    Storage = null;
                    Offset = -12364;
                    CI = -12349;
                    CJ = -667;
                    return;
                } else {
                    // mem alloc necessary 
                    AllocSblk(out B, out MembnkIdx, out InMembnk, ISblk, JSblk);
                    BE.InMembnk[rowSblk, colSblk] = InMembnk;
                    BE.MembnkIdx[rowSblk, colSblk] = MembnkIdx;

                    bool isExternal = jBlock < m_ColPartitioning.FirstBlock || jBlock >= m_ColPartitioning.FirstBlock + m_ColPartitioning.LocalNoOfBlocks;
                    if(isExternal) {
                        // +++++++++++++++++++++++++++++++++++++
                        // allocated an *external* block 
                        // this *may* invalidate the communication pattern.
                        // +++++++++++++++++++++++++++++++++++++

                        // block range

                        long j0, jE;
                        BlockMsrMatrix.GetExternalSubblockIndices(this._ColPartitioning, jBlock, out j0, out jE);
                        int OwnerRank = m_ColPartitioning.FindProcess(j0);
                        Debug.Assert(OwnerRank != m_ColPartitioning.MpiRank);
                        HashSet<long> BlockIndicesByProcessor;
                        if(!m_ExternalBlockIndicesByProcessor.TryGetValue(OwnerRank, out BlockIndicesByProcessor)) {
                            BlockIndicesByProcessor = new HashSet<long>();
                            m_ExternalBlockIndicesByProcessor.Add(OwnerRank, BlockIndicesByProcessor);
                        }
                        this.ComPatternValid &= !BlockIndicesByProcessor.Add(jBlock);
                        this.m_ExternalBlock[iBlockLoc] |= true;
                    }

                }
            } else {
                B = m_Membanks[MembnkIdx];
            }

            // get access pattern
            bool isDense;
            B.GetFastBlockAccessInfo(out Storage, out Offset, out CI, out CJ, out isDense, InMembnk);

        }

        /// <summary>
        /// Accumulates <paramref name="Block"/>*<paramref name="alpha"/> to this matrix,
        /// at the row/column offset <paramref name="i0"/> resp. <paramref name="j0"/>.
        /// </summary>
        /// <param name="i0">Row offset.</param>
        /// <param name="j0">Column offset.</param>
        /// <param name="alpha">Scaling factor for the accumulation operation.</param>
        /// <param name="Block">Block to accumulate.</param>
        public void AccBlock(long i0, long j0, double alpha, MultidimensionalArray Block) {
            this.AccBlock(i0, j0, alpha, Block, 1.0);
        }


        /// <summary>
        /// Accumulates <paramref name="Block"/>*<paramref name="alpha"/> to this matrix,
        /// at the row/column offset <paramref name="i0"/> resp. <paramref name="j0"/>.
        /// </summary>
        /// <param name="i0">Row offset.</param>
        /// <param name="j0">Column offset.</param>
        /// <param name="alpha">Scaling factor for the accumulation operation.</param>
        /// <param name="Block">Block to accumulate.</param>
        /// <param name="beta">Scaling applied to this matrix before accumulation;
        /// Note: if set to 0.0, this will also clear NAN's, etc.
        /// </param>
        public void AccBlock(long i0, long j0, double alpha, MultidimensionalArray Block, double beta) {
            if(Block.Dimension != 2)
                throw new ArgumentException("Expecting a 2D array.");
            int I = Block.NoOfRows;
            int J = Block.NoOfCols;
            if(i0 < this.m_RowPartitioning.i0)
                throw new ArgumentException("Row index out of range.");
            if(i0 + I > this.m_RowPartitioning.iE)
                throw new ArgumentException("Row index out of range.");
            if(j0 < 0)
                throw new ArgumentException("Column index out of range.");
            if(j0 + J > this.ColPartition.TotalLength)
                throw new ArgumentException("Column index out of range.");
            if(I <= 0 || J <= 0)
                return;


#if DEBUG_EXTENDED
            BitArray bTouch = new BitArray(I * J);
#endif

            int j;
            int i = 0;
            while(i < I) {
                j = 0;

                int IWRT = -1;
                while(j < J) {

                    int iSblk, jSblk; //   row/col index within sub-block, which correspond to (i0 + i),(j0 + j)
                    int ISblk, JSblk; //   sub-block size
                    double[] Storage; //   sub-block memory: where the result should be accumulated
                    int Offset, CI, CJ; // offset pointer and i,j cycles into 'Storage'

                    GetSetAlloc(true, i + i0, j + j0, out iSblk, out jSblk, out ISblk, out JSblk, out Storage, out Offset, out CI, out CJ);

                    //this.ColPartition.MpiRank

                    int IWrt = Math.Min(I - i, ISblk - iSblk); // number of rows to write
                    int JWrt = Math.Min(J - j, JSblk - jSblk); // number of columns to write
                    Debug.Assert(IWRT < 0 || IWRT == IWrt);
                    IWRT = IWrt;

                    if(Storage != null) {
                        for(int iw = 0; iw < IWrt; iw++) {
                            for(int jw = 0; jw < JWrt; jw++) {
#if DEBUG_EXTENDED
                                Debug.Assert(bTouch[(i + iw) * J + j + jw] == false);
                                bTouch[(i + iw) * J + j + jw] = true;
#endif
                                int StorageIdx = Offset + (iSblk + iw) * CI + (jSblk + jw) * CJ;
                                if(beta == 0.0)
                                    Storage[StorageIdx] = alpha * Block[i + iw, j + jw];
                                else
                                    Storage[StorageIdx] = Storage[StorageIdx] * beta + alpha * Block[i + iw, j + jw];
                            }
                        }
                    } else {
                        for(int iw = 0; iw < IWrt; iw++) {
                            for(int jw = 0; jw < JWrt; jw++) {
#if DEBUG_EXTENDED
                                Debug.Assert(bTouch[(i + iw) * J + j + jw] == false);
                                bTouch[(i + iw) * J + j + jw] = true;
#endif
                                if(Block[i + iw, j + jw] != 0.0) {
                                    throw new ArgumentException("Can not save non-zero entry in void-region.");
                                }
                            }
                        }
                    }

                    j += JWrt;
                    Debug.Assert(JWrt > 0);
                }
                i += IWRT;
                Debug.Assert(IWRT > 0);
            }

#if DEBUG_EXTENDED
            for (i = 0; i < I; i++)
                for (j = 0; j < J; j++)
                    Debug.Assert(bTouch[i * J + j] == true);
#endif

        }


        /// <summary>
        /// Extracts a block from this matrix and stores it into <paramref name="Block"/>.
        /// </summary>
        /// <param name="i0">Row index to start reading form.</param>
        /// <param name="j0">Column index to start reading form.</param>
        /// <param name="Block">Output; allocated by caller.</param>
        public void ReadBlock(long i0, long j0, MultidimensionalArray Block) {
            if(Block.Dimension != 2)
                throw new ArgumentException("Expecting a 2D array.");
            int I = Block.NoOfRows;
            int J = Block.NoOfCols;
            if(i0 < this.m_RowPartitioning.i0)
                throw new ArgumentException("Row index out of range.");
            if(i0 + I > this.m_RowPartitioning.iE)
                throw new ArgumentException("Row index out of range.");
            if(j0 < 0)
                throw new ArgumentException("Column index out of range.");
            if(j0 + J > this.ColPartition.TotalLength)
                throw new ArgumentException("Column index out of range.");
            if(I <= 0 || J <= 0)
                return;

#if DEBUG_EXTENDED
            BitArray bTouch = new BitArray(I * J);
#endif

            int j;
            int i = 0;
            while(i < I) {
                j = 0;

                int IWRT = -1;
                while(j < J) {

                    int iSblk, jSblk; //   row/col index within sub-block, which correspond to (i0 + i),(j0 + j)
                    int ISblk, JSblk; //   sub-block size
                    double[] Storage; //   sub-block memory: where the result should be accumulated
                    int Offset, CI, CJ; // offset pointer and i,j cycles into 'Storage'

                    GetSetAlloc(false, i + i0, j + j0, out iSblk, out jSblk, out ISblk, out JSblk,
                        out Storage, out Offset, out CI, out CJ);

                    Debug.Assert(I - i > 0);
                    Debug.Assert(J - j > 0);
                    int IWrt = Math.Min(I - i, ISblk - iSblk); // number of rows to write
                    int JWrt = Math.Min(J - j, JSblk - jSblk); // number of columns to write
                    Debug.Assert(IWRT < 0 || IWRT == IWrt);
                    IWRT = IWrt;

                    if(Storage != null) {
                        Debug.Assert(IWrt > 0);
                        Debug.Assert(JWrt > 0);
                        for(int iw = 0; iw < IWrt; iw++) {
                            for(int jw = 0; jw < JWrt; jw++) {
#if DEBUG_EXTENDED
                                Debug.Assert(bTouch[(i + iw) * J + j + jw] == false);
                                bTouch[(i + iw) * J + j + jw] = true;
#endif
                                int StorageIdx = Offset + (iSblk + iw) * CI + (jSblk + jw) * CJ;
                                Block[i + iw, j + jw] = Storage[StorageIdx];
                            }
                        }
                    } else {
                        Debug.Assert(IWrt > 0);

                        if(JSblk <= 0) {
                            // nothing allocated in respective block-row

                            int __iBlkLoc, __iLoc, __BlkT, __Sblk_idx, __i0_Sblk, __ISblk, __NoOfSblk, __RemainingVoidWidth;
                            TranslateIndex(j + j0, m_ColPartitioning, true, out long __iBlk, out __iBlkLoc, out long __i0, out __iLoc, out __BlkT, out __Sblk_idx, out __i0_Sblk, out __ISblk, out __NoOfSblk, out __RemainingVoidWidth);
                            //    )

                            if(__RemainingVoidWidth > 0) {
                                JWrt = Math.Min(J - j, __RemainingVoidWidth);
                            } else {
                                JWrt = Math.Min(J - j, __ISblk + __i0_Sblk - __iLoc);
                            }

                            //JWrt = 1; // work column by column
                        }

                        Debug.Assert(JWrt > 0);
                        for(int iw = 0; iw < IWrt; iw++) {
                            for(int jw = 0; jw < JWrt; jw++) {
#if DEBUG_EXTENDED
                                Debug.Assert(bTouch[(i + iw) * J + j + jw] == false);
                                bTouch[(i + iw) * J + j + jw] = true;
#endif
                                Block[i + iw, j + jw] = 0.0;
                            }
                        }
                    }

                    Debug.Assert(JWrt >= 0);
                    j += JWrt;
                }
                Debug.Assert(IWRT > 0);
                i += IWRT;
            }

#if DEBUG_EXTENDED
            for (i = 0; i < I; i++)
                for (j = 0; j < J; j++)
                    Debug.Assert(bTouch[i * J + j] == true);
#endif
        }


        /// <summary>
        /// Clears a block in this matrix and frees the allocated memory.
        /// </summary>
        /// <param name="i0">Row index to start clearing.</param>
        /// <param name="j0">Column index to start clearing.</param>
        /// <param name="I">Number of rows to clear.</param>
        /// <param name="J">Number of columns to clear.</param>
        public void ClearBlock(long i0, long j0, int I, int J) {
            if(I < 0)
                throw new ArgumentOutOfRangeException();
            if(J < 0)
                throw new ArgumentOutOfRangeException();
            if(i0 < this.m_RowPartitioning.i0)
                throw new ArgumentException("Row index out of range.");
            if(i0 + I > this.m_RowPartitioning.iE)
                throw new ArgumentException("Row index out of range.");
            if(j0 < 0)
                throw new ArgumentException("Column index out of range.");
            if(j0 + J > this.ColPartition.TotalLength)
                throw new ArgumentException("Column index out of range.");
            if(I <= 0 || J <= 0)
                return;

#if DEBUG_EXTENDED
            BitArray bTouch = new BitArray(I * J);
#endif

            int j;
            int i = 0;
            while(i < I) {
                j = 0;

                int IWRT = -1;
                while(j < J) {

                    long BlockRow, BlockCol; //    block row and column
                    int rowSblkIdx, colSblkIdx; // sub-block row and column
                    int iSblk, jSblk; //           row/col index within sub-block, which correspond to (i0 + i),(j0 + j)
                    int ISblk, JSblk; //           sub-block size
                    double[] Storage; //           sub-block memory
                    int Offset, CI, CJ; //         offset pointer and i,j cycles into 'Storage'
                    int MembnkIdx, InMembnk; //    which bank, which index within bank

                    GetSetAlloc(false, i + i0, j + j0, out BlockRow, out BlockCol, out rowSblkIdx, out colSblkIdx, out iSblk, out jSblk, out ISblk, out JSblk, out Storage, out Offset, out CI, out CJ, out MembnkIdx, out InMembnk);
                    Debug.Assert((Storage != null) == (MembnkIdx >= 0 && InMembnk >= 0));

                    int IWrt = Math.Min(I - i, ISblk - iSblk); // number of rows to write
                    int JWrt = Math.Min(J - j, JSblk - jSblk); // number of columns to write
                    Debug.Assert(IWRT < 0 || IWRT == IWrt);
                    IWRT = IWrt;

                    if(Storage != null) {
                        // need to release memory

                        if(IWRT == ISblk && JWrt == JSblk) {
                            // release sub-block

                            Debug.Assert(iSblk == 0);
                            Debug.Assert(jSblk == 0);

                            FreeSblk(MembnkIdx, InMembnk);

#if DEBUG_EXTENDED
                            for (int iw = 0; iw < IWrt; iw++) {
                                for (int jw = 0; jw < JWrt; jw++) {
                                    Debug.Assert(bTouch[(i + iw) * J + j + jw] == false);
                                    bTouch[(i + iw) * J + j + jw] = true;
                                }
                            }
#endif
                            BlockEntry BE = m_BlockRows[BlockRow - m_RowPartitioning.FirstBlock][BlockCol];
                            Debug.Assert(BlockCol == BE.jBlkCol);
                            Debug.Assert(BE.MembnkIdx[rowSblkIdx, colSblkIdx] == MembnkIdx);
                            Debug.Assert(BE.InMembnk[rowSblkIdx, colSblkIdx] == InMembnk);
                            BE.MembnkIdx[rowSblkIdx, colSblkIdx] = -1;
                            BE.InMembnk[rowSblkIdx, colSblkIdx] = -1;

                            if(BE.IsEmpty) {
                                m_BlockRows[BlockRow - m_RowPartitioning.FirstBlock].Remove(BE.jBlkCol);

                                /*
                                bool isExternal = (BE.jBlkCol < m_ColPartitioning.FirstBlock) || (m_ColPartitioning.LocalNoOfBlocks + m_ColPartitioning.FirstBlock <= BE.jBlkCol);

                                if(isExternal) {
                                    int OwnerProc = m_ColPartitioning.FindProcessForBlock(BE.jBlkCol);
                                    Debug.Assert(OwnerProc != m_ColPartitioning.MpiRank);
                                    HashSet<int> ExtIdx = m_ExternalBlockIndicesByProcessor[OwnerProc];
                                    Debug.Assert(ExtIdx.Contains(BE.jBlkCol));
                                    ExtIdx.Remove(BE.jBlkCol);
                                }
                                */

                                // Reminder: the code above is wrong, for the following reason:
                                // 'm_ExternalBlockIndicesByProcessor[OwnerProc]' is a collection of external blocks owned by 'OwnerProc'
                                // for **all** rows -- we would have to check all other rows in order to safely remove 
                                // 'BE.jBlkCol' from 'm_ExternalBlockIndicesByProcessor[OwnerProc]'.
                                // On the other hand, having unused external block in the set may impact MPI communication, but 
                                // does not impact the numerical results.

                            }

                        } else {
                            // set entries to 0.0

                            for(int iw = 0; iw < IWrt; iw++) {
                                for(int jw = 0; jw < JWrt; jw++) {
#if DEBUG_EXTENDED
                                    Debug.Assert(bTouch[(i + iw) * J + j + jw] == false);
                                    bTouch[(i + iw) * J + j + jw] = true;
#endif
                                    int StorageIdx = Offset + (iSblk + iw) * CI + (jSblk + jw) * CJ;
                                    Storage[StorageIdx] = 0.0;
                                }
                            }
                        }
                    } else {
#if DEBUG_EXTENDED
                            for (int iw = 0; iw < IWrt; iw++) {
                                for (int jw = 0; jw < JWrt; jw++) {
                                    Debug.Assert(bTouch[(i + iw) * J + j + jw] == false);
                                    bTouch[(i + iw) * J + j + jw] = true;
                                }
                            }
#endif
                    }

                    Debug.Assert(JWrt > 0);
                    j += JWrt;
                }
                Debug.Assert(IWRT > 0);
                i += IWRT;
            }

#if DEBUG_EXTENDED
            for (i = 0; i < I; i++)
                for (j = 0; j < J; j++)
                    Debug.Assert(bTouch[i * J + j] == true);
#endif
        }

        /// <summary>
        /// Computes all kinds of indices from a global index.
        /// </summary>
        /// <param name="i">Row or column index.</param>
        /// <param name="part">Partitioning.</param>
        /// <param name="SupportExternal">
        /// If true, <paramref name="i"/> can be outside of the range of local indices of <paramref name="part"/>.
        /// </param>
        /// <param name="iBlk">On exit, global (over all MPI processes) block index.</param>
        /// <param name="iBlkLoc">On exit, local (on this MPI process) block index.</param>
        /// <param name="i0">Global start index of the block</param>
        /// <param name="iLoc">Index within block.</param>
        /// <param name="BlkT">Block type, see <see cref="IBlockPartitioning.GetBlockType"/></param>
        /// <param name="Sblk_idx">Sub block index.</param>
        /// <param name="i0_Sblk">First block index in sub block.</param>
        /// <param name="ISblk">Size/length of sub block</param>
        /// <param name="NoOfSblk">Number of sub blocks in block.</param>
        /// <param name="RemainingVoidWidth">If in void region (<paramref name="i0_Sblk"/> smaller than 0), the distance to the next sub-block or the end of the block.</param>
        static private void TranslateIndex(long i, IBlockPartitioning part, bool SupportExternal, out long iBlk, out int iBlkLoc, out long i0, out int iLoc, out int BlkT, out int Sblk_idx, out int i0_Sblk, out int ISblk, out int NoOfSblk, out int RemainingVoidWidth) {
            if(i >= part.i0 && i < part.iE) {

                iBlk = part.GetBlockIndex(i);
                iBlkLoc = (int)(iBlk - part.FirstBlock);
                Debug.Assert(iBlk >= part.FirstBlock && iBlk < part.LocalNoOfBlocks + part.FirstBlock);

                i0 = part.GetBlockI0(iBlk);
                iLoc = (int)(i - i0);
                int BlkLen = part.GetBlockLen(iBlk);
                Debug.Assert(iLoc >= 0);
                Debug.Assert(iLoc < BlkLen);

                BlkT = part.GetBlockType(iBlk);

                int[] Sblk_i0, SblkLen;
                Sblk_i0 = part.GetSubblk_i0(BlkT);
                SblkLen = part.GetSubblkLen(BlkT);
                Debug.Assert(Sblk_i0.Length == SblkLen.Length);
                NoOfSblk = Sblk_i0.Length;


                i0_Sblk = 0;
                ISblk = 0;
                int iE_Sblk = 0;
                int Prev_iE_Sblk = 0;
                for(Sblk_idx = 0; Sblk_idx < NoOfSblk; Sblk_idx++) {
                    i0_Sblk = Sblk_i0[Sblk_idx];
                    ISblk = SblkLen[Sblk_idx];
                    iE_Sblk = i0_Sblk + ISblk;
                    if(i0_Sblk <= iLoc && iLoc < iE_Sblk) {
                        RemainingVoidWidth = -1;
                        return;
                    }
                    /*
                    if (Sblk_idx < NoOfSblk - 1 && iLoc >= i0_Sblk && iLoc < Sblk_i0[Sblk_idx + 1]) {
                        // in void - region (between blocks)
                        i0_Sblk = int.MinValue;
                        Sblk_idx = -371827;
                        ISblk = -44266;
#if DEBUG_EXTENDED
                        for (int _Sblk_idx = 0; _Sblk_idx < NoOfSblk; _Sblk_idx++) {
                            int __i0_Sblk = Sblk_i0[_Sblk_idx];
                            int __ISblk = SblkLen[_Sblk_idx];
                            int __iE_Sblk = __i0_Sblk + __ISblk;
                            Debug.Assert((__i0_Sblk <= iLoc && iLoc < __iE_Sblk) == false); // not really in void-section
                        }
#endif
                        RemainingVoidWidth = Sblk_i0[Sblk_idx + 1] - iLoc;
                        return;
                    }
                    */
                    if(iLoc >= Prev_iE_Sblk && iLoc < i0_Sblk) {
                        // in void - region (between blocks)
                        RemainingVoidWidth = i0_Sblk - iLoc;
                        i0_Sblk = int.MinValue;
                        Sblk_idx = -371827;
                        ISblk = -44266;
#if DEBUG_EXTENDED
                        for (int _Sblk_idx = 0; _Sblk_idx < NoOfSblk; _Sblk_idx++) {
                            int __i0_Sblk = Sblk_i0[_Sblk_idx];
                            int __ISblk = SblkLen[_Sblk_idx];
                            int __iE_Sblk = __i0_Sblk + __ISblk;
                            Debug.Assert((__i0_Sblk <= iLoc && iLoc < __iE_Sblk) == false); // not really in void-section
                        }
#endif
                        return;
                    }
                    Prev_iE_Sblk = iE_Sblk;
                }

                // in void - region
                Debug.Assert(Sblk_idx == NoOfSblk);
                //Debug.Assert(iLoc >= Sblk_i0[NoOfSblk - 1] + SblkLen[NoOfSblk - 1]);
#if DEBUG_EXTENDED
                if (NoOfSblk > 0)
                    Debug.Assert(iLoc >= Sblk_i0[NoOfSblk - 1] + SblkLen[NoOfSblk - 1]);

                for (int _Sblk_idx = 0; _Sblk_idx < NoOfSblk; _Sblk_idx++) {
                    int __i0_Sblk = Sblk_i0[_Sblk_idx];
                    int __ISblk = SblkLen[_Sblk_idx];
                    int __iE_Sblk = __i0_Sblk + __ISblk;
                    Debug.Assert((__i0_Sblk <= iLoc && iLoc < __iE_Sblk) == false); // not really in void-section
                }
#endif
                i0_Sblk = int.MinValue;
                Sblk_idx = -371827;
                ISblk = -44266;
                RemainingVoidWidth = BlkLen - iLoc;
                return;
            } else {
                if(!SupportExternal) {
                    throw new ArgumentException("Index out of local range.");
                } else {
                    if(i < 0 || i >= part.TotalLength)
                        throw new ArgumentException("Index out of range.");
                }

                // column index outside local column range
                NoOfSblk = 1;


                long iE;
                GetExternalSubblockIndices(part, i, out iBlk, out i0, out iE);
                Debug.Assert(iBlk < part.FirstBlock || iBlk >= part.LocalNoOfBlocks + part.FirstBlock);


                BlkT = int.MinValue;
                i0_Sblk = 0;
                iBlkLoc = int.MinValue;
                iLoc = int.MinValue;
                ISblk = (int)(iE - i0);
                Sblk_idx = 0;
                RemainingVoidWidth = -1;
                return;
            }
        }

        /// <summary>
        /// Computes the blocking for external indices (outside of the local range, i.e. smaller than <see cref="IPartitioning.i0"/> or larger or equal to <see cref="IPartitioning.iE"/>).
        /// </summary>
        /// <param name="part"></param>
        /// <param name="iBlk">Input; an external block index.</param>
        /// <param name="i0">first index of the block</param>
        /// <param name="iE">first index of the next block</param>
        /// <returns>
        /// owner rank of the external sub-block
        /// </returns>
        static int GetExternalSubblockIndices(IBlockPartitioning part, long iBlk, out long i0, out long iE) {
            int Rank = part.FindProcessForBlock(iBlk);
            Debug.Assert(Rank != part.MpiRank);
            long NoBlkLoc = part.GetLocalNoOfBlocks(Rank);
            long L = part.GetLocalLength(Rank);
            long FistBlock = part.GetFirstBlock(Rank);
            int iBlkLoc = checked((int)(iBlk - FistBlock));
            Debug.Assert(iBlkLoc >= 0);
            Debug.Assert(iBlkLoc < NoBlkLoc);

            int i0Loc = (int)((iBlkLoc * L) / NoBlkLoc);
            int iELoc = (int)(((iBlkLoc + 1) * L) / NoBlkLoc);

            long i0_rank = part.GetI0Offest(Rank);
            i0 = i0Loc + i0_rank;
            iE = iELoc + i0_rank;

#if DEBUG_EXTENDED
            for (long i = i0; i < iE; i++) {
                long check_i0, check_iE, check_iBlk;
                GetExternalSubblockIndices(part, i, out check_iBlk, out check_i0, out check_iE);
                Debug.Assert(check_i0 == i0);
                Debug.Assert(check_iE == iE);
                Debug.Assert(check_iBlk == iBlk);
            }
#endif
            return Rank;
        }

        /// <summary>
        /// Computes the blocking for external indices (outside of the local range, i.e. smaller than <see cref="IPartitioning.i0"/> or larger or equal to <see cref="IPartitioning.iE"/>).
        /// </summary>
        /// <param name="part"></param>
        /// <param name="i">Input; an external index</param>
        /// <param name="iBlk">block index</param>
        /// <param name="i0">first index of the block</param>
        /// <param name="iE">last index of the block</param>
        static void GetExternalSubblockIndices(IBlockPartitioning part, long i, out long iBlk, out long i0, out long iE) {
            Debug.Assert(i < part.i0 || i >= part.iE);

            // we use a lot of integer operations like (a*b)/c here, which give a much more even block distribution than a*(b/c)
            // the downside is this can cause integer overflow

            int Rank = part.FindProcess(i);
            long i0_Rank = part.GetI0Offest(Rank);
            long iLoc = i - i0_Rank;
            long L = part.GetLocalLength(Rank);
            Debug.Assert(iLoc >= 0);
            Debug.Assert(iLoc < L);
            long NoBlkLoc = part.GetLocalNoOfBlocks(Rank);

            long iBlkLoc = (long)Math.Floor(((double)(iLoc * NoBlkLoc)) / ((double)L)); // get a first guess
            long i0Loc = iBlkLoc * L / NoBlkLoc;
            long iELoc = (iBlkLoc + 1) * L / NoBlkLoc;

            while(iLoc < i0Loc) {
                iBlkLoc--;
                i0Loc = (iBlkLoc * L) / NoBlkLoc;
                iELoc = ((iBlkLoc + 1) * L) / NoBlkLoc;
            }
            while(iLoc >= iELoc) {
                iBlkLoc++;
                i0Loc = (iBlkLoc * L) / NoBlkLoc;
                iELoc = ((iBlkLoc + 1) * L) / NoBlkLoc;
            }
            Debug.Assert(iBlkLoc >= 0);
            Debug.Assert(iBlkLoc < NoBlkLoc);
            Debug.Assert(iLoc >= i0Loc);
            Debug.Assert(iLoc < iELoc);

            iBlk = (int)(part.GetFirstBlock(Rank) + iBlkLoc);
            i0 = (int)(i0Loc + i0_Rank);
            iE = (int)(iELoc + i0_Rank);
        }


        /// <summary>
        /// Allocates memory for a sub block.
        /// </summary>
        /// <param name="B">On exit, the memory bank at which the sub block was allocated.</param>
        /// <param name="ISblk">Number of rows for sub block.</param>
        /// <param name="JSblk">Number of columns for sub block.</param>
        /// <param name="MembnkIdx">Memory bank index, see <see cref="BlockEntry.MembnkIdx"/>.</param>
        /// <param name="InMembnk">
        /// Block index within memory bank, 
        /// i.e. first index into <see cref="Membank.Mem"/>, see <see cref="BlockEntry.InMembnk"/>.
        /// </param>
        void AllocSblk(out Membank B, out int MembnkIdx, out int InMembnk, int ISblk, int JSblk) {
            int NoOfMembnk = m_Membanks.Count;
            for(int iMbnk = 0; iMbnk < NoOfMembnk; iMbnk++) { // loop over memory banks, trying to find a suitable one..
                B = m_Membanks[iMbnk];
                if(B.NextFree >= 0 && B.Mem.GetLength(1) == ISblk && B.Mem.GetLength(2) == JSblk) {
                    // found a memory bank which has still available space and 
                    // is of the correct size

                    MembnkIdx = iMbnk;
                    InMembnk = B.NextFree;

#if DEBUG_EXTENDED
                    // check that there are no invalid entries left
                    double[] RawMem;
                    int Offset, CI, CJ;
                    bool isDense;
                    B.GetFastBlockAccessInfo(out RawMem, out Offset, out CI, out CJ, out isDense, InMembnk);

                    for (int i = 0; i < ISblk; i++) {
                        for (int j = 0; j < JSblk; j++) {
                            Debug.Assert(RawMem[Offset + i * CI + j * CJ] == 0);
                        }
                    }
#endif
                    // find next free block (i.e. ensure that 'B.NextFree' is valid)
                    B.Occupied[InMembnk] = true;
                    int L = B.Occupied.Count;
                    Debug.Assert(B.Occupied.Count == B.Mem.GetLength(0));
                    B.NextFree = int.MinValue;
                    for(int i = InMembnk + 1; i < L; i++) {
                        if(B.Occupied[i] == false) {
                            B.NextFree = i;
                            break;
                        }
                    }

                    // work done
                    return;
                }
            }

            // found nothing so far - must allocate new memory bank.
            B = new Membank(this.m_RowPartitioning.LocalNoOfBlocks, ISblk, JSblk);
            InMembnk = 0;
            MembnkIdx = m_Membanks.Count;
            m_Membanks.Add(B);
            B.Occupied[0] = true;
            if(B.Occupied.Count > 1) {
                B.NextFree++;
            } else {
                B.NextFree = int.MinValue;
            }
        }

        /// <summary>
        /// Releases a sub-block.
        /// </summary>
        void FreeSblk(int MembnkIdx, int InMembnk) {
            // find mem bank
            Membank B = m_Membanks[MembnkIdx];

            // reset memory state
            double[] RawMem;
            int Offset, CI, CJ, ISblk = B.Mem.GetLength(1), JSblk = B.Mem.GetLength(2);
            bool isDense;
            B.GetFastBlockAccessInfo(out RawMem, out Offset, out CI, out CJ, out isDense, InMembnk);
            for(int i = 0; i < ISblk; i++) {
                for(int j = 0; j < JSblk; j++) {
                    RawMem[Offset + i * CI + j * CJ] = 0;
                }
            }

            // release block
            B.Occupied[InMembnk] = false;
            B.NextFree = Math.Min(B.NextFree, InMembnk);
        }

        /// <summary>
        /// see <see cref="IMutableMatrixEx.GetOccupiedColumnIndices"/>;
        /// the returned array is a non-shallow copy of the row;
        /// </summary>
        public int GetRow(long i, ref long[] ColumnIndices, ref double[] Values) {
            return GetRowInternal(i, ref ColumnIndices, ref Values, false);
        }

        /// <summary>
        /// Common implementation for <see cref="GetRow"/> and <see cref="GetOccupiedColumnIndices"/>.
        /// </summary>
        int GetRowInternal(long i, ref long[] ColumnIndices, ref double[] Values, bool NoValues) {
            Debug.Assert(m_RowPartitioning.IsInLocalRange(i));

            // =====================
            // translate row indices
            // =====================
            int iBlkLoc, iLoc, iBlkT, rowSblkIdx, i0_Sblk, ISblk, NoOfRowSblk, RemVoidRows;
            TranslateIndex(i, m_RowPartitioning, false, out long iBlk, out iBlkLoc, out long i0, out iLoc, out iBlkT, out rowSblkIdx, out i0_Sblk, out ISblk, out NoOfRowSblk, out RemVoidRows);
            if(i0_Sblk < 0) {
                // in void region
                Debug.Assert(RemVoidRows > 0);
                return 0;
            }

            int iSblk = (int)(i - i0 - i0_Sblk);
            Debug.Assert(iSblk >= 0);
            Debug.Assert(iSblk < m_RowPartitioning.GetSubblkLen(iBlkT)[rowSblkIdx]);

            // ==================================
            // compute number of occupied entries
            // ==================================
            int nnz = -10;

            var BlockRows = this.m_BlockRows[iBlkLoc];
            if(BlockRows == null) {
                nnz = 0;
            } else {
                GetRowInternal_pass(null, null, true, rowSblkIdx, iSblk, ref nnz, BlockRows);
            }

            // ======================
            // allocate return memory
            // ======================
            if(ColumnIndices == null || ColumnIndices.Length < nnz)
                ColumnIndices = new long[nnz];
            if((Values == null || Values.Length < nnz) && (NoValues == false))
                Values = new double[nnz];

            // ======================
            // set return values
            // ======================

            if(nnz > 0) {
                GetRowInternal_pass(ColumnIndices, Values, NoValues, rowSblkIdx, iSblk, ref nnz, BlockRows);
            }

            // ======================
            // return
            // ======================
            return nnz;
        }

        void GetRowInternal_pass(long[] ColumnIndices, double[] Values, bool NoValues, int rowSblkIdx, int iSblk, ref int nnz, SortedDictionary<long, BlockEntry> BlockRows) {
            int cnt = 0;
            foreach(var kv in BlockRows) { // loop over blocks in row...
                long jBlk = kv.Key;
                BlockEntry B = kv.Value;
                Debug.Assert(jBlk == B.jBlkCol);

                int[] ColSubblkLen, ColSubblk_i0;
                int NoOfSblkCol;
                long j0, jE;
                if(jBlk >= m_ColPartitioning.FirstBlock && jBlk < m_ColPartitioning.LocalNoOfBlocks + m_ColPartitioning.FirstBlock) {
                    int jBlkT = m_ColPartitioning.GetBlockType(jBlk);
                    ColSubblkLen = m_ColPartitioning.GetSubblkLen(jBlkT);
                    ColSubblk_i0 = m_ColPartitioning.GetSubblk_i0(jBlkT);
                    NoOfSblkCol = B.InMembnk.GetLength(1);
                    Debug.Assert(ColSubblkLen.Length == NoOfSblkCol);
                    Debug.Assert(ColSubblk_i0.Length == NoOfSblkCol);

                    j0 = m_ColPartitioning.GetBlockI0(jBlk);
                    jE = long.MinValue;
                } else {
                    GetExternalSubblockIndices(m_ColPartitioning, jBlk, out j0, out jE);
                    NoOfSblkCol = 1;
                    ColSubblkLen = null;
                    ColSubblk_i0 = null;

                }

                for(int colSblkIdx = 0; colSblkIdx < NoOfSblkCol; colSblkIdx++) { // loop over sub-blocks in block...
                    if(B.MembnkIdx[rowSblkIdx, colSblkIdx] >= 0) {
                        // occupied/allocated

                        Membank Mbnk = m_Membanks[B.MembnkIdx[rowSblkIdx, colSblkIdx]];
                        double[] Rawmem;
                        int Offset, CI, CJ;
                        bool isDense;
                        Mbnk.GetFastBlockAccessInfo(out Rawmem, out Offset, out CI, out CJ, out isDense, B.InMembnk[rowSblkIdx, colSblkIdx]);
                        Debug.Assert(iSblk >= 0);
                        Debug.Assert(iSblk < Mbnk.Mem.GetLength(1));

                        int RowOffset = Offset + CI * iSblk;

                        long JSb, j00;
                        if(jE < 0) {
                            j00 = ColSubblk_i0[colSblkIdx] + j0;
                            JSb = ColSubblkLen[colSblkIdx];
                        } else {
                            j00 = j0;
                            JSb = jE - j0;
                        }
                        Debug.Assert(Mbnk.Mem.GetLength(2) == JSb);

                        for(int jsb = 0; jsb < JSb; jsb++) { // loop over columns in sub-block...
                            long j = j00 + jsb;
                            Debug.Assert(jE >= 0 || j >= m_ColPartitioning.GetBlockI0(jBlk));
                            Debug.Assert(jE >= 0 || j < m_ColPartitioning.GetBlockLen(jBlk) + m_ColPartitioning.GetBlockI0(jBlk));

                            double vat = Rawmem[RowOffset + jsb * CJ];
                            if(vat != 0.0) {
                                if(ColumnIndices != null) {
                                    ColumnIndices[cnt] = j;
                                    if(NoValues == false)
                                        Values[cnt] = vat;
                                }
                                cnt++;
                            }
                        }
                    }
                }
            }

            if(nnz >= 0) {
                Debug.Assert(cnt == nnz);
            } else {
                nnz = cnt;
            }
        }

        /// <summary>
        /// All stores sub-blocks.
        /// </summary>         
        List<Membank> m_Membanks = new List<Membank>();

        /// <summary>
        /// Memory bank, stores sub-blocks of equal size,
        /// </summary>
        class Membank : ICloneable {
            /// <summary>
            /// ctor.
            /// </summary>
            public Membank(int L, int ISblk, int JSblk) {

                long MaxL = ((long)int.MaxValue) / (((long)L) * ISblk * JSblk);
                if (L > MaxL)
                    L = checked((int)MaxL);
                if (L < 1)
                    L = 1;

                Mem = MultidimensionalArray.Create(L, ISblk, JSblk);
                Occupied = new BitArray(L);
                NextFree = 0;
            }

            /// <summary>
            /// Storage to store sub-blocks.
            /// - 1st index: enumeration of sub-blocks
            /// - 2nd index: sub-block row index
            /// - 3rd index: sub-block column index.
            /// </summary>
            public MultidimensionalArray Mem;

            /// <summary>
            /// Number of rows in each sub-block
            /// </summary>
            public int SubblockNoRows {
                get {
                    return Mem.GetLength(1);
                }
            }

            /// <summary>
            /// Number of columns in each sub-block
            /// </summary>
            public int SubblockNoCols {
                get {
                    return Mem.GetLength(2);
                }
            }


            /// <summary>
            /// Marks which parts of <see cref="Mem"/> are occupied (true).
            /// Index correlates with first index of <see cref="Mem"/>,
            /// </summary>
            public BitArray Occupied;

            /// <summary>
            /// First free entry in <see cref="Mem"/>, correlates with first index of <see cref="Mem"/>.
            /// </summary>
            public int NextFree;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="RawMem"></param>
            /// <param name="Offset"></param>
            /// <param name="CI">Cycle for advancing the row index.</param>
            /// <param name="CJ">Cycle for advancing the column index.</param>
            /// <param name="isDense"></param>
            /// <param name="InMembnk">
            /// Memory bank index, i.e. first index into <see cref="Mem"/>.
            /// </param>
            public void GetFastBlockAccessInfo(out double[] RawMem, out int Offset, out int CI, out int CJ, out bool isDense, int InMembnk) {
                RawMem = Mem.Storage;
                Offset = Mem.Index(InMembnk, 0, 0);
                if(Mem.GetLength(1) > 1) {
                    CI = Mem.Index(InMembnk, 1, 0) - Offset;
                } else {
                    CI = Mem.GetLength(2);
                }

                if(Mem.GetLength(2) > 1) {
                    CJ = Mem.Index(InMembnk, 0, 1) - Offset;
                } else {
                    CJ = 1;
                }

                isDense = (CJ == 1) && (CI == Mem.GetLength(2));
            }

            private Membank() {
            }

            /// <summary>
            /// Cloning.
            /// </summary>
            public object Clone() {
                Membank R = new Membank();
                R.Mem = this.Mem.CloneAs();
                R.NextFree = this.NextFree;
                R.Occupied = this.Occupied.CloneAs();
                return R;
            }
        }

        /// <summary>
        /// A flag which indicates that a solver can assume this matrix to be symmetric;<br/>
        /// For performance reasons, we rely on the user to set this correctly (checking
        /// for symmetry is very expensive even on one pro
        /// </summary>
        public bool AssumeSymmetric {
            get {
                return m_AssumeSymmetric;
            }
            set {
                m_AssumeSymmetric = value;
            }
        }

        bool m_AssumeSymmetric;

        /// <summary>
        /// creates a non-shallow copy of this object
        /// </summary>
        public BlockMsrMatrix CloneAs() {
            using(new FuncTrace()) {
                BlockMsrMatrix Ret = new BlockMsrMatrix();

                Ret.m_RowPartitioning = this._RowPartitioning;
                Ret.m_ColPartitioning = this._ColPartitioning;

                Ret.m_AssumeSymmetric = this.m_AssumeSymmetric;
                Ret.m_ExternalBlock = this.m_ExternalBlock.CloneAs();
                Ret.m_ExternalBlockIndicesByProcessor = new Dictionary<int, HashSet<long>>();
                foreach(var kv in this.m_ExternalBlockIndicesByProcessor) {
                    Ret.m_ExternalBlockIndicesByProcessor.Add(kv.Key, new HashSet<long>(kv.Value));
                }

                Ret.m_BlockRows = new SortedDictionary<long, BlockEntry>[this.m_BlockRows.Length];
                for(int iRow = 0; iRow < Ret.m_BlockRows.Length; iRow++) {
                    var Row = this.m_BlockRows[iRow];
                    if(Row != null) {
                        var NewRow = new SortedDictionary<long, BlockEntry>();
                        foreach(var kv in Row) {
                            NewRow.Add(kv.Key, kv.Value.CloneAs());
                        }
                        Ret.m_BlockRows[iRow] = NewRow;
                    }
                }

                Ret.m_Membanks = new List<Membank>();
                foreach(var MB in this.m_Membanks) {
                    Ret.m_Membanks.Add(MB.CloneAs());
                }

                Ret.ComPatternValid = this.ComPatternValid;
                if(Ret.ComPatternValid) {
                    if(this.SendLists != null) {
                        Ret.SendLists = new Dictionary<int, long[,]>();
                        foreach(var kv in this.SendLists) {
                            Ret.SendLists.Add(kv.Key, kv.Value.CloneAs());
                        }
                    }

                    if(this.ReceiveLists != null) {
                        Ret.ReceiveLists = new Dictionary<int, long[,]>();
                        foreach(var kv in this.ReceiveLists) {
                            Ret.ReceiveLists.Add(kv.Key, kv.Value.CloneAs());
                        }
                    }
                }

                return Ret;
            }
        }

        /// <summary>
        /// creates a non-shallow copy of this object
        /// </summary>
        public object Clone() {
            return CloneAs();
        }

        /// <summary>
        /// Sets all entries to 0, all memory is released.
        /// </summary>
        public void Clear() {
            m_Membanks.Clear();
            m_ExternalBlockIndicesByProcessor.Clear();
            Array.Clear(m_BlockRows, 0, m_BlockRows.Length);
            SendLists = new Dictionary<int, long[,]>();
            ReceiveLists = new Dictionary<int, long[,]>();
            m_ExternalBlock.SetAll(false);
            ComPatternValid = true;
        }

        /// <summary>
        /// Number of Bytes used
        /// </summary>
        public long UsedMemory {
            get {
                GetMemoryInfo(out _, out long Used);
                return Used;
            }
        }

        /// <summary>
        /// Computes the amount of memory allocated resp. used by this matrix.
        /// </summary>
        /// <param name="Allocated">
        /// Allocated memory, in bytes
        /// </param>
        /// <param name="Used">
        /// Actually used memory (always smaller or equal to <paramref name="Allocated"/>) memory.
        /// </param>
        public void GetMemoryInfo(out long Allocated, out long Used) {
            Allocated = 0;
            foreach(var mbk in m_Membanks) {
                Allocated += mbk.Mem.Length * sizeof(double);
            }

            Used = 0;
            for(int iBlockLoc = 0; iBlockLoc < m_BlockRows.Length; iBlockLoc++) {
                var BlockRow = this.m_BlockRows[iBlockLoc];
                if(BlockRow != null) {
                    long iBlockGlb = iBlockLoc + this._RowPartitioning.FirstBlock;
                    long i0 = this._RowPartitioning.GetBlockI0(iBlockGlb);

                    foreach(var KV in BlockRow) {
                        long jBlockGlb = KV.Key;
                        BlockEntry Block = KV.Value;
                        Debug.Assert(jBlockGlb == Block.jBlkCol);

                        int BlockPointerMem = (Block.InMembnk.Length + Block.MembnkIdx.Length + 1) * sizeof(int);
                        Used += BlockPointerMem;
                        Allocated += BlockPointerMem;


                        if(_ColPartitioning.IsLocalBlock(jBlockGlb)) {
                            long j0 = this._ColPartitioning.GetBlockI0(jBlockGlb);

                            int NoOfSblk_rows = Block.InMembnk.GetLength(0);
                            int NoOfSblk_cols = Block.InMembnk.GetLength(1);
                            int RowBlockType = _RowPartitioning.GetBlockType(iBlockGlb);
                            int ColBlockType = _ColPartitioning.GetBlockType(jBlockGlb);

                            int[] RowSblkLen = _RowPartitioning.GetSubblkLen(RowBlockType);
                            Debug.Assert(RowSblkLen.Length == NoOfSblk_rows);

                            int[] ColSblkLen = _ColPartitioning.GetSubblkLen(ColBlockType);
                            Debug.Assert(ColSblkLen.Length == NoOfSblk_cols);

                            for(int sblkRow = 0; sblkRow < NoOfSblk_rows; sblkRow++) {
                                for(int sblkCol = 0; sblkCol < NoOfSblk_cols; sblkCol++) {
                                    if(Block.InMembnk[sblkCol, sblkRow] < 0 || Block.MembnkIdx[sblkCol, sblkRow] < 0)
                                        continue;

                                    int M = RowSblkLen[sblkRow];
                                    int N = ColSblkLen[sblkCol];

                                    Used += M * N * sizeof(double);
                                }
                            }
                        } else {
                            // todo 

                            GetExternalSubblockIndices(_ColPartitioning, jBlockGlb, out long j0, out long jE);

                            int NoOfSblk_rows = Block.InMembnk.GetLength(0);
                            int NoOfSblk_cols = Block.InMembnk.GetLength(1);

                            for(int sblkRow = 0; sblkRow < NoOfSblk_rows; sblkRow++) {
                                for(int sblkCol = 0; sblkCol < NoOfSblk_cols; sblkCol++) {
                                    int MembnkIdx = Block.MembnkIdx[sblkCol, sblkRow];
                                    int InMembnk = Block.InMembnk[sblkCol, sblkRow];

                                    if(MembnkIdx < 0 || InMembnk < 0)
                                        continue;

                                    int M = m_Membanks[MembnkIdx].SubblockNoRows;
                                    int N = m_Membanks[MembnkIdx].SubblockNoCols;

                                    Used += M * N * sizeof(double);
                                }
                            }

                        }


                    }
                }
            }

        }


        /// <summary>
        /// Ad-hoc performance instrumentation
        /// </summary>
        public static Stopwatch SPMV_tot = new Stopwatch();

        /// <summary>
        /// Ad-hoc performance instrumentation
        /// </summary>
        public static Stopwatch SPMV_inner = new Stopwatch();

        /// <summary>
        /// Ad-hoc performance instrumentation
        /// </summary>
        public static Stopwatch SpMV_local = new Stopwatch();

        /// <summary>
        /// Ad-hoc performance instrumentation
        /// </summary>
        public static Stopwatch SpMV_initSending = new Stopwatch();

        /// <summary>
        /// Ad-hoc performance instrumentation
        /// </summary>
        public static Stopwatch SpMV_receive = new Stopwatch();

        /// <summary>
        /// Ad-hoc performance instrumentation
        /// </summary>
        public static Stopwatch SpMV_external = new Stopwatch();

        /// <summary>
        /// Ad-hoc performance instrumentation
        /// </summary>
        public static Stopwatch SpMV_indextrans = new Stopwatch();

        /// <summary>
        /// Write the ad-hoc instrumentation to Console
        /// </summary>
        public static void PrintPerfStat() {
            if(BlockMsrMatrix.multiply != null)
                Console.WriteLine("  spmm total " + BlockMsrMatrix.multiply.Elapsed.TotalSeconds);
            if (BlockMsrMatrix.multiply_core != null)
                Console.WriteLine("  spmm core " + BlockMsrMatrix.multiply_core.Elapsed.TotalSeconds);
            
            Console.WriteLine("  spmv total     " + BlockMsrMatrix.SPMV_tot.Elapsed.TotalSeconds);
            Console.WriteLine("   spmv local    " + BlockMsrMatrix.SpMV_local.Elapsed.TotalSeconds);
            Console.WriteLine("    spmv inner   " + BlockMsrMatrix.SPMV_inner.Elapsed.TotalSeconds);
            Console.WriteLine("   spmv send     " + BlockMsrMatrix.SpMV_initSending.Elapsed.TotalSeconds);
            Console.WriteLine("   spmv receive  " + BlockMsrMatrix.SpMV_receive.Elapsed.TotalSeconds);
            Console.WriteLine("   spmv external " + BlockMsrMatrix.SpMV_external.Elapsed.TotalSeconds);
        }

        /// <summary>
        /// Sparse Matrix/Vector multiplication;
        /// the calculation 
        /// <paramref name="acc"/> = <paramref name="acc"/>*<paramref name="beta"/> + this*<paramref name="_a"/>*<paramref name="alpha"/>
        /// is performed;
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="_a">
        /// input; vector to be multiplied with this matrix from the right
        /// </param>
        /// <param name="beta"></param>
        /// <param name="acc">
        /// length of accumulator must be at least the update length 
        /// </param>
        public void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 _a, double beta, VectorType2 acc)
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> //
        {
            using(new FuncTrace()) {
                //SPMV_tot.Reset();
                //SPMV_inner.Reset();
                //SpMV_local.Reset();
                //SpMV_initSending.Reset();
                //SpMV_receive.Reset();
                //SpMV_external.Reset();
                //SpMV_indextrans.Reset();

                SPMV_tot.Start();

                if(_a.Count != this._ColPartitioning.LocalLength)
                    throw new ArgumentException("Mismatch in number of columns.");
                if(acc.Count != this._RowPartitioning.LocalLength)
                    throw new ArgumentException("Mismatch in number of rows.");
                if(object.ReferenceEquals(acc, _a) && (acc.Count > 0 || _a.Count > 0))
                    throw new ArgumentException("In-Place computation is not supported.", "acc");

                double[] a;
                if(_a is double[])
                    a = _a as double[];
                else
                    a = _a.ToArray(); // still faster than accessing via IList

                SpMV_initSending.Start();
                this.UpdateCommPattern(this.MPI_Comm);
                int MPIsize = this._RowPartitioning.MpiSize;
                int MPIrank = this._RowPartitioning.MpiRank;

                //Dictionary<int, double[]> SendBuf = new Dictionary<int, double[]>();

                int[] SendRanks = this.SendLists.Keys.ToArray();
                IntPtr[] SendBuffers = new IntPtr[SendRanks.Length];
                int[] SendBuffersLen = new int[SendRanks.Length];

                int[] RecvRanks = this.ReceiveLists.Keys.ToArray();
                IntPtr[] RecvBuffers = new IntPtr[RecvRanks.Length];
                int[] RecvBuffersLen = new int[RecvRanks.Length];

                MPI_Request[] Requests = new MPI_Request[SendRanks.Length + RecvRanks.Length];

                long iFirstRow = this._RowPartitioning.i0;
                long iFirstCol = this._ColPartitioning.i0;

                // prepare receiving
                // =================
                for(int i = 0; i < RecvRanks.Length; i++) {
                    int RecvRnk = RecvRanks[i];
                    long[,] RcvList = this.ReceiveLists[RecvRnk];

                    int Len = 0;
                    int L0 = RcvList.GetLength(0);
                    for(int l = 0; l < L0; l++) {
                        Len += (int)(RcvList[l, 1] - RcvList[l, 0]);
                    }
                    RecvBuffersLen[i] = Len;

                    RecvBuffers[i] = Marshal.AllocHGlobal(Len * sizeof(double));
                    csMPI.Raw.Irecv(RecvBuffers[i], Len, csMPI.Raw._DATATYPE.DOUBLE, RecvRnk, 55 * 19 * RecvRnk, this.MPI_Comm, out Requests[SendRanks.Length + i]);
                }

                // send data 
                // =========

                unsafe {
                    for(int i = 0; i < SendRanks.Length; i++) {
                        int SendRnk = SendRanks[i];
                        long[,] SndList = this.SendLists[SendRnk];

                        int Len = 0;
                        int L0 = SndList.GetLength(0);
                        for(int l = 0; l < L0; l++) {
                            Len += (int)(SndList[l, 1] - SndList[l, 0]);
                        }
                        SendBuffersLen[i] = Len;

                        SendBuffers[i] = Marshal.AllocHGlobal(Len * sizeof(double));
                        double* sb = (double*)(SendBuffers[i]);

                        for(int l = 0; l < L0; l++) { // loop over all chunks/blocks of the send list
                            long i0 = SndList[l, 0];
                            long iE = SndList[l, 1];

                            for(long iRow = i0; iRow < iE; iRow++) {
                                Debug.Assert(_ColPartitioning.IsInLocalRange(iRow));
                                int iRowLoc = (int)(iRow - iFirstCol);
                                Debug.Assert(iRowLoc >= 0 && iRowLoc < a.Length);
                                *sb = a[iRowLoc];
                                sb++;
                            }
                        }
                        Debug.Assert((sb - ((double*)(SendBuffers[i]))) == Len);

                        csMPI.Raw.Issend(SendBuffers[i], Len, csMPI.Raw._DATATYPE.DOUBLE, SendRnk, 55 * 19 * MPIrank, this.MPI_Comm, out Requests[i]);
                    }
                }
                SpMV_initSending.Stop();
                bool[] bTouch = new bool[acc.Count];

                // local multiplication
                // ====================

                SpMV_local.Start();

                int NoOfBlockRows = _RowPartitioning.LocalNoOfBlocks;
                Debug.Assert(NoOfBlockRows == m_BlockRows.Length);
                long FirstRowBlock = _RowPartitioning.FirstBlock;
                unsafe {
                    fixed(double* pa = a) {
                        double[] VecAccu = null;
                        for(int iBlockLoc = 0; iBlockLoc < NoOfBlockRows; iBlockLoc++) { // loop over block rows...
                            var BlockRow = m_BlockRows[iBlockLoc];

                            long iBlock = iBlockLoc + FirstRowBlock;
                            int RowBlockLength = _RowPartitioning.GetBlockLen(iBlock);
                            int locBlockRowOffset = (int)(_RowPartitioning.GetBlockI0(iBlock) - _RowPartitioning.i0);

                            if(BlockRow != null) {
                                int RowBlockType = _RowPartitioning.GetBlockType(iBlock);
                                int[] Row_i0Sblk = _RowPartitioning.GetSubblk_i0(RowBlockType);
                                int[] RowLenSblk = _RowPartitioning.GetSubblkLen(RowBlockType);
                                bool ContainsExternal = false;

                                if(VecAccu == null || VecAccu.Length < RowBlockLength) {
                                    VecAccu = new double[RowBlockLength];
                                } else {
                                    for(int i = 0; i < RowBlockLength; i++) {
                                        VecAccu[i] = 0.0;
                                    }
                                }
                                fixed(double* pVecAccu = VecAccu) {

                                    foreach(var kv in BlockRow) { // loop over block columns...
                                        BlockEntry BE = kv.Value;
                                        long jBlkCol = kv.Key;
                                        Debug.Assert(BE.jBlkCol == jBlkCol);
                                        if(_ColPartitioning.IsLocalBlock(jBlkCol)) {
                                            int locBlockColOffset = (int)(_ColPartitioning.GetBlockI0(jBlkCol) - _ColPartitioning.i0);

                                            int ColBlockType = _ColPartitioning.GetBlockType(jBlkCol);
                                            int[] Col_i0Sblk = _ColPartitioning.GetSubblk_i0(ColBlockType);
                                            int[] ColLenSblk = _ColPartitioning.GetSubblkLen(ColBlockType);

                                            Debug.Assert(BE.MembnkIdx.GetLength(0) == BE.InMembnk.GetLength(0));
                                            Debug.Assert(BE.MembnkIdx.GetLength(1) == BE.InMembnk.GetLength(1));
                                            int NoOfSblk_Rows = BE.MembnkIdx.GetLength(0);
                                            int NoOfSblk_Cols = BE.MembnkIdx.GetLength(1);
                                            Debug.Assert(Row_i0Sblk.Length == NoOfSblk_Rows);
                                            Debug.Assert(RowLenSblk.Length == NoOfSblk_Rows);
                                            Debug.Assert(Col_i0Sblk.Length == NoOfSblk_Cols);
                                            Debug.Assert(ColLenSblk.Length == NoOfSblk_Cols);
                                            for(int iSblkRow = 0; iSblkRow < NoOfSblk_Rows; iSblkRow++) { // loop over sub-block rows

                                                int _iRowLoc = locBlockRowOffset + Row_i0Sblk[iSblkRow]; // local row index
                                                int _iRowBlockLoc = Row_i0Sblk[iSblkRow]; // row index within block

                                                for(int jSblkCol = 0; jSblkCol < NoOfSblk_Cols; jSblkCol++) { // loop over sub-block columns
                                                    int MembnkIdx = BE.MembnkIdx[iSblkRow, jSblkCol];
                                                    int InMembnk = BE.InMembnk[iSblkRow, jSblkCol];
                                                    Debug.Assert((MembnkIdx >= 0) == (InMembnk >= 0));

                                                    if(InMembnk >= 0) {
                                                        int Offset, CI, CJ;
                                                        bool isDense;
                                                        m_Membanks[MembnkIdx].GetFastBlockAccessInfo(out double[] RawMem, out Offset, out CI, out CJ, out isDense, InMembnk);

                                                        fixed(double* pRawMem = RawMem) {

                                                            int _jColLoc = locBlockColOffset + Col_i0Sblk[jSblkCol]; // local storage index

                                                            int I = RowLenSblk[iSblkRow];
                                                            int J = ColLenSblk[jSblkCol];
                                                            Debug.Assert(m_Membanks[MembnkIdx].Mem.GetLength(1) == I);
                                                            Debug.Assert(m_Membanks[MembnkIdx].Mem.GetLength(2) == J);

                                                            int iRowLoc = _iRowLoc; // local row index
                                                            int iRowBlockLoc = _iRowBlockLoc; // row index within block

                                                            if(CJ != 1 || I * J < 40) { 
                                                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                                                // local multiplication branch:
                                                                // either the block is very small (less than 40 entries), or
                                                                // the fast rotating index has a skips some bytes
                                                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                                                for(int i = 0; i < I; i++) { // loop over sub-block rows...

                                                                    double Accu = 0;

                                                                    //int jColLoc = _jColLoc;
                                                                    double* pRow = pRawMem + Offset + CI * i;
                                                                    double* pCol = pa + _jColLoc;
                                                                    for(int j = 0; j < J; j++) { // loop over sub-block columns...
                                                                                                 // .
                                                                                                 //int iStorage = Offset + CI * i + CJ * j; // index into memory bank
                                                                                                 //Accu += RawMem[iStorage] * a[jColLoc];
                                                                                                 //jColLoc++;

                                                                        Accu += *pRow * *pCol;
                                                                        pCol++;
                                                                        pRow += CJ;
                                                                    }
                                                                    VecAccu[iRowBlockLoc] += Accu;
                                                                    iRowLoc++;
                                                                    iRowBlockLoc++;
                                                                }
                                                            } else {
                                                                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                                                                // default branch: try to use BLAS function DGEMV
                                                                // +++++++++++++++++++++++++++++++++++++++++++++++++++
                                                                SPMV_inner.Start();
                                                                BLAS.dgemv('t', J, I, 1.0, pRawMem + Offset, CI, pa + _jColLoc, 1, 1.0,
                                                                    pVecAccu + _iRowBlockLoc,
                                                                    1);
                                                                SPMV_inner.Stop();
                                                            }

                                                        }

                                                    }
                                                }
                                            }

                                        } else {
                                            ContainsExternal = true;
                                        }

                                    }

                                    for(int __i = 0; __i < RowBlockLength; __i++) {
                                        int __iRowLoc = __i + locBlockRowOffset;
                                        double ri = acc[__iRowLoc] * beta + alpha * VecAccu[__i];
                                        acc[__iRowLoc] = ri;
                                        Debug.Assert(acc[__iRowLoc] == ri, "ri = " + ri + "acc = " + acc[__iRowLoc]);
                                    }
                                }

                                Debug.Assert(m_ExternalBlock[iBlockLoc] == ContainsExternal);
                            } else {
                                for(int __i = 0; __i < RowBlockLength; __i++) {
                                    int __iRowLoc = __i + locBlockRowOffset;
                                    double ri = acc[__iRowLoc] * beta;
                                    acc[__iRowLoc] = ri;
                                    Debug.Assert(acc[__iRowLoc] == ri);
                                }
                            }
                        }
                    }
                }
                SpMV_local.Stop();

                // finish communication
                // ====================
                SpMV_receive.Start();
                MPI_Status[] Stats = new MPI_Status[Requests.Length];
                csMPI.Raw.Waitall(Requests.Length, Requests, Stats);
                SpMV_receive.Stop();
#if DEBUG_EXTENDED
                /*
                for (int index = 0; index < Requests.Length; index++) {
                    //MPI_Status Stat;
                    //int index;
                    //csMPI.Raw.Waitany(Requests.Length, Requests, out index, out Stat);
                    var Stat = Stats[index];

                    if (index < SendRanks.Length) {
                        // finished sending
                        Marshal.FreeHGlobal(SendBuffers[index]);
                        SendBuffers[index] = IntPtr.Zero;
                        Debug.Assert(Stat.MPI_SOURCE == SendRanks[index]);
                        Debug.Assert(Stat.MPI_TAG == 55 * 19 * MPIrank);
                        Debug.Assert(Stat.count == SendBuffersLen[index] * sizeof(double));
                    } else {
                        // finished receiving
                        Debug.Assert(Stat.MPI_SOURCE == RecvRanks[index - SendRanks.Length]);
                        Debug.Assert(Stat.MPI_TAG == 55 * 19 * RecvRanks[index - SendRanks.Length]);
                        Debug.Assert(Stat.count == RecvBuffersLen[index - SendRanks.Length] * sizeof(double));
                    }
                }
                */
#endif

                // external multiplication
                // ====================
                SpMV_external.Start();
                for(int iBlockLoc = 0; iBlockLoc < NoOfBlockRows; iBlockLoc++) { // loop over block rows...
                    if(m_ExternalBlock[iBlockLoc]) {

                        var BlockRow = m_BlockRows[iBlockLoc];
                        Debug.Assert(BlockRow != null);
                        long iBlock = iBlockLoc + FirstRowBlock;
                        int RowBlockType = _RowPartitioning.GetBlockType(iBlock);
                        int[] Row_i0Sblk = _RowPartitioning.GetSubblk_i0(RowBlockType);
                        int[] RowLenSblk = _RowPartitioning.GetSubblkLen(RowBlockType);
                        int locBlockRowOffset = (int)(_RowPartitioning.GetBlockI0(iBlock) - _RowPartitioning.i0);

                        foreach(var kv in BlockRow) { // loop over block columns...
                            BlockEntry BE = kv.Value;
                            long jBlkCol = kv.Key;
                            Debug.Assert(BE.jBlkCol == jBlkCol);

                            int OwnerRank = -1;
                            long[,] RcvList = null;
                            int RcvListBlockRow = 0;
                            IntPtr RecvBuffer = IntPtr.Zero;
                            int OffetInto_RecvBuffer = 0;

                            if(!_ColPartitioning.IsLocalBlock(jBlkCol)) {
                                //int OwnerProc = _ColPartitioning.FindProcessForBlock(jBlkCol);
                                //throw new NotImplementedException("para todo");

                                //int locBlockColOffset = _ColPartitioning.GetBlockI0(jBlkCol) - _ColPartitioning.i0;
                                long j0, jE;
                                int _OwnerRank = GetExternalSubblockIndices(_ColPartitioning, jBlkCol, out j0, out jE);
                                if(_OwnerRank != OwnerRank) {
                                    RcvList = ReceiveLists[_OwnerRank];
                                    RcvListBlockRow = 0;
                                    OffetInto_RecvBuffer = 0;
                                    OwnerRank = _OwnerRank;
                                    int bufferIdx = Array.IndexOf(RecvRanks, OwnerRank);
                                    RecvBuffer = RecvBuffers[bufferIdx];
                                }
                                while(RcvList[RcvListBlockRow, 0] < j0) {
                                    OffetInto_RecvBuffer += (int)(RcvList[RcvListBlockRow, 1] - RcvList[RcvListBlockRow, 0]);
                                    RcvListBlockRow++;
                                }
                                Debug.Assert(RcvList[RcvListBlockRow, 0] == j0);
                                Debug.Assert(RcvList[RcvListBlockRow, 1] == jE);

                                Debug.Assert(BE.MembnkIdx.GetLength(0) == BE.InMembnk.GetLength(0));
                                Debug.Assert(BE.MembnkIdx.GetLength(1) == BE.InMembnk.GetLength(1));
                                int NoOfSblk_Rows = BE.MembnkIdx.GetLength(0);
                                int NoOfSblk_Cols = BE.MembnkIdx.GetLength(1);
                                Debug.Assert(NoOfSblk_Cols == 1);
                                Debug.Assert(Row_i0Sblk.Length == NoOfSblk_Rows);
                                Debug.Assert(RowLenSblk.Length == NoOfSblk_Rows);

                                for(int iSblkRow = 0; iSblkRow < NoOfSblk_Rows; iSblkRow++) { // loop over sub-block rows
                                    int MembnkIdx = BE.MembnkIdx[iSblkRow, 0];
                                    int InMembnk = BE.InMembnk[iSblkRow, 0];
                                    Debug.Assert((MembnkIdx >= 0) == (InMembnk >= 0));

                                    if(InMembnk >= 0) {
                                        double[] RawMem;
                                        int Offset, CI, CJ;
                                        bool isDense;
                                        m_Membanks[MembnkIdx].GetFastBlockAccessInfo(out RawMem, out Offset, out CI, out CJ, out isDense, InMembnk);
                                        unsafe {
                                            double* dRecvBuffer = (double*)RecvBuffer;

                                            int I = RowLenSblk[iSblkRow];
                                            int J = (int)(jE - j0);
                                            Debug.Assert(I == m_Membanks[MembnkIdx].Mem.GetLength(1));
                                            Debug.Assert(J == m_Membanks[MembnkIdx].Mem.GetLength(2));
                                            for(int i = 0; i < I; i++) { // loop over sub-block rows...
                                                int iRowLoc = locBlockRowOffset + i + Row_i0Sblk[iSblkRow]; // local row index
                                                Debug.Assert(iRowLoc >= 0 && iRowLoc < _RowPartitioning.LocalLength);
                                                double Accu = 0;

                                                for(int j = 0; j < J; j++) { // loop over sub-block columns...
                                                    int iRcvBuff = OffetInto_RecvBuffer + j;
                                                    int iStorage = Offset + CI * i + CJ * j; // index into memory bank
#if DEBUG_EXTENDED
                                                    SpMV_indextrans.Start();
                                                    long jColGlob = j0 + j; // global column index
                                                    int ownerRank;
                                                    int iList = BruteForceExternalIndexTranslation(jColGlob, out ownerRank);
                                                    int __bufferIdx = Array.IndexOf(RecvRanks, OwnerRank);
                                                    SpMV_indextrans.Stop();
                                                    Debug.Assert(ownerRank == OwnerRank);
                                                    Debug.Assert(RecvBuffers[__bufferIdx] == RecvBuffer);
                                                    Debug.Assert(RcvList[RcvListBlockRow, 0] <= jColGlob);
                                                    Debug.Assert(RcvList[RcvListBlockRow, 1] > jColGlob);
                                                    Debug.Assert(iRcvBuff == iList);

#endif
                                                    //double* dRecvBuffer = (double*)RecvBuffers[bufferIdx];
                                                    //Accu += RawMem[iStorage] * dRecvBuffer[iList];

                                                    Accu += RawMem[iStorage] * dRecvBuffer[iRcvBuff];
                                                }

                                                acc[iRowLoc] += alpha * Accu;
                                            }

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                SpMV_external.Stop();

                // free temp buffers
                // =================
                {
                    foreach(IntPtr b in SendBuffers)
                        Marshal.FreeHGlobal(b);
                    foreach(IntPtr b in RecvBuffers)
                        Marshal.FreeHGlobal(b);
                }



                SPMV_tot.Stop();
            }
        }



        int BruteForceExternalIndexTranslation(long GlobIndex, out int OwnerProc) {
            OwnerProc = m_ColPartitioning.FindProcess(GlobIndex);
            long[,] Rcvlist = ReceiveLists[OwnerProc];

            int L = Rcvlist.GetLength(0);
            Debug.Assert(Rcvlist.GetLength(1) == 2);
            int ListIdx = 0;
            for(int l = 0; l < L; l++) {
                long i0 = Rcvlist[l, 0];
                long iE = Rcvlist[l, 1];
                for(long i = i0; i < iE; i++) {
                    if(i == GlobIndex)
                        return ListIdx;
                    ListIdx++;
                }
            }

            throw new ApplicationException();
        }

        /// <summary>
        /// Tests the integrity of the block matrix data structure.
        /// </summary>
        [Conditional("DEBUG")]
        private void VerifyDataStructure(string message) {
            for(int iBlockRowLoc = 0; iBlockRowLoc < this.m_BlockRows.Length; iBlockRowLoc++) {

                VerifyDataStructure_Row(iBlockRowLoc, message);
            }
        }

        [Conditional("DEBUG")]
        static void MyAssert(bool b, string message) {
            if(!b) {
                //Console.WriteLine();
                // dbg_launch();
                throw new ApplicationException("Data integrity of BlockMsrMatrix lost: " + message);
            }
        }

        [Conditional("DEBUG")]
        private void VerifyDataStructure_Row(int iBlockRowLoc, string matrixName) {
            long iBlockRow = iBlockRowLoc + m_RowPartitioning.FirstBlock;
            int RowBlockType = m_RowPartitioning.GetBlockType(iBlockRow);
            int[] _i0SubBlock = m_RowPartitioning.GetSubblk_i0(RowBlockType);
            int[] RowLenSubBlock = m_RowPartitioning.GetSubblkLen(RowBlockType);
            MyAssert(RowLenSubBlock.Length == _i0SubBlock.Length, matrixName + ": 1");

            if(m_BlockRows[iBlockRowLoc] != null) {
                foreach(var kv in m_BlockRows[iBlockRowLoc]) {
                    long jBlockCol = kv.Key;
                    BlockEntry BE = kv.Value;
                    MyAssert(BE.jBlkCol == jBlockCol, matrixName + ": 2");

                    MyAssert(BE.MembnkIdx.GetLength(0) == _i0SubBlock.Length, matrixName + ": 3");
                    MyAssert(BE.InMembnk.GetLength(0) == _i0SubBlock.Length, matrixName + ": 4");

                    int[] _j0SubBlock, ColLenSubBlock;
                    if(m_ColPartitioning.IsLocalBlock(jBlockCol)) {
                        int ColBlockType = m_ColPartitioning.GetBlockType(jBlockCol);
                        _j0SubBlock = m_ColPartitioning.GetSubblk_i0(ColBlockType);
                        ColLenSubBlock = m_ColPartitioning.GetSubblkLen(ColBlockType);
                    } else {
                        long j0, jE;
                        GetExternalSubblockIndices(m_ColPartitioning, jBlockCol, out j0, out jE);
                        _j0SubBlock = new int[] { 0 };
                        ColLenSubBlock = new int[] { (int)(jE - j0) };
                    }
                    MyAssert(ColLenSubBlock.Length == _j0SubBlock.Length, matrixName + ": 5");
                    MyAssert(BE.MembnkIdx.GetLength(1) == _j0SubBlock.Length, matrixName + ": 6");
                    MyAssert(BE.InMembnk.GetLength(1) == _j0SubBlock.Length, matrixName + ": 7");

                    for(int sblkRow = 0; sblkRow < _i0SubBlock.Length; sblkRow++) {
                        for(int sblkCol = 0; sblkCol < _j0SubBlock.Length; sblkCol++) {
                            MyAssert((BE.InMembnk[sblkRow, sblkCol] >= 0) == (BE.MembnkIdx[sblkRow, sblkCol] >= 0), matrixName + ": A");

                            if(BE.InMembnk[sblkRow, sblkCol] >= 0) {
                                var Mbnk = m_Membanks[BE.MembnkIdx[sblkRow, sblkCol]];
                                int inBank = BE.InMembnk[sblkRow, sblkCol];
                                MyAssert(Mbnk.Mem.GetLength(0) > inBank, matrixName + ": B");
                                MyAssert(Mbnk.Mem.GetLength(1) == RowLenSubBlock[sblkRow], matrixName + ": C");
                                MyAssert(Mbnk.Mem.GetLength(2) == ColLenSubBlock[sblkCol], matrixName + ": D");
                                MyAssert(Mbnk.Occupied[inBank] == true, matrixName + ": E");

                                double[] RawMem;
                                int Offset, CI, CJ;
                                bool isDense;
                                Mbnk.GetFastBlockAccessInfo(out RawMem, out Offset, out CI, out CJ, out isDense, inBank);

                                MyAssert(RawMem.Length > Offset + CI * (RowLenSubBlock[sblkRow] - 1) + CJ * (ColLenSubBlock[sblkCol] - 1), matrixName + ": F");
                            }
                        }
                    }
                }
            }
        }

        /// <summary> 
        /// Similar to <see cref="AccSubMatrixTo{V1,V2,V3,V4}"/>, 
        /// but the destination (<paramref name="target"/>) is cleared before the accumulation.
        /// </summary>
        public void WriteSubMatrixTo<V1, V2, V3, V4>(IMutableMatrixEx target,
            V1 RowIndicesSource, V2 RowIndicesTarget,
            V3 ColumnIndicesSource, V4 ColumnInidcesTarget)
            where V1 : IList<long>
            where V2 : IList<long>
            where V3 : IList<long>
            where V4 : IList<long> {

            target.Clear();
            this.AccSubMatrixTo(
                1.0,
                target,
                RowIndicesSource,
                RowIndicesTarget,
                ColumnIndicesSource,
                ColumnInidcesTarget,
                default(long[]), default(long[]));
        }

        /// <summary>
        /// See <see cref="AccSubMatrixTo{V1, V2, V3, V4, V5, V6}(double, IMutableMatrixEx, V1, V2, V3, V4, V5, V6)"/>, with less paramters.
        /// </summary>
        public void AccSubMatrixTo<V1, V2, V3, V4>(
            double alpha, IMutableMatrixEx Target,
            V1 RowIndicesSource, V2 RowIndicesTarget,
            V3 ColumnIndicesSource, V4 ColIndicesTarget)
            where V1 : IEnumerable<long>
            where V2 : IEnumerable<long>
            where V3 : IEnumerable<long>
            where V4 : IEnumerable<long>//
        {
            AccSubMatrixTo<V1, V2, V3, V4, long[], long[]>(alpha, Target, RowIndicesSource, RowIndicesTarget, ColumnIndicesSource, ColIndicesTarget, null, null);
        }

        /// <summary>
        /// Extracts a submatrix from this matrix and accumulates it to another matrix.
        /// </summary>
        /// <param name="alpha">
        /// scaling factor
        /// </param>
        /// <param name="RowIndicesSource">Row indices into this matrix, in the local range (<see cref="IPartitioning.IsInLocalRange"/>).</param>
        /// <param name="RowIndicesTarget">
        /// if null is specified, this array is assumed to be
        /// {0,1, ... , <paramref name="RowIndicesSource"/>.Length-1]};
        /// </param>
        /// <param name="ColumnIndicesSource">
        /// Column indices into this matrix, in the local range (<see cref="IPartitioning.IsInLocalRange"/>).
        /// </param>
        /// <param name="ColIndicesTarget">
        /// if null is specified, this array is assumed to be
        /// {0,1, ... , <paramref name="ColumnIndicesSource"/>.Length-1]};
        /// </param>
        /// <param name="Target">
        /// Output:
        /// The [<paramref name="RowIndicesTarget"/>[i],<paramref name="ColIndicesTarget"/>[j]]-th
        /// entry of the returned matrix is equal to the 
        /// [<paramref name="RowIndicesSource"/>[i],<paramref name="ColumnIndicesSource"/>[j]]-th
        /// entry of this matrix.
        /// </param>
        /// <param name="ExternalColIndicesTarget"></param>
        /// <param name="ExternalColumnIndicesSource">
        /// Additional/optional row indices into this matrix in the external range (<see cref="IPartitioning.IsInLocalRange"/> evaluates to false).
        /// These are not exchanged over MPI.
        /// </param>
        /// <remarks>
        /// Note on MPI parallelization: it is *not* necessary to exchange column indices, this can be done by the method internally.
        /// The parameters <paramref name="ExternalColumnIndicesSource"/> and <paramref name="ExternalColIndicesTarget"/> can be used if
        /// the <paramref name="ColIndicesTarget"/> is on another MPI communicator as this object (<see cref="MPI_Comm"/>).
        /// </remarks>
        public void AccSubMatrixTo<V1, V2, V3, V4, V5, V6>(
            double alpha, IMutableMatrixEx Target,
            V1 RowIndicesSource, V2 RowIndicesTarget,
            V3 ColumnIndicesSource, V4 ColIndicesTarget,
            V5 ExternalColumnIndicesSource, V6 ExternalColIndicesTarget)
            where V1 : IEnumerable<long>
            where V2 : IEnumerable<long>
            where V3 : IEnumerable<long>
            where V4 : IEnumerable<long>
            where V5 : IEnumerable<long>
            where V6 : IEnumerable<long> //
        {
            using(new FuncTrace("BlockMsrMatrix.AccSubMatrixTo")) {

                // check input arguments
                // =====================
#if DEBUG_EXTENDED
                int iPrev = -1;
                foreach (int i in RowIndicesSource) {
                    this.RowPartitioning.TestIfInLocalRange(i);
                    //if (i <= iPrev)
                    //    throw new ArgumentException("Indices must be in strictly ascending order.");
                    iPrev = i;
                }
                iPrev = -1;
                foreach (int i in ColumnIndicesSource) {
                    this.m_ColPartitioning.TestIfInLocalRange(i);
                    //if (i <= iPrev)
                    //    throw new ArgumentException("Indices must be in strictly ascending order.");
                    iPrev = i;
                }

                if (RowIndicesTarget != null && RowIndicesTarget.Count() != RowIndicesSource.Count()) {
                    throw new ArgumentException("Mismatch between number of source and target indices (rows).");
                }
                if (ColIndicesTarget != null && ColIndicesTarget.Count() != ColumnIndicesSource.Count()) {
                    throw new ArgumentException("Mismatch between number of source and target indices (columns).");
                }

                if((ExternalColumnIndicesSource == null) != (ExternalColIndicesTarget == null)) {
                    throw new ArgumentException("Mismatch between external column source and target indices.");
                }

                if(ExternalColumnIndicesSource != null) {
                    if(ExternalColumnIndicesSource.Count() != ExternalColIndicesTarget.Count()) {
                        throw new ArgumentException("Mismatch between number of external source and target indices.");
                    }

                    foreach(int i in ExternalColumnIndicesSource) {
                        if(this.m_ColPartitioning.IsInLocalRange(i) == true)
                            throw new ArgumentException("External column indices are expected to be external.");
                    }
                }
#endif
                // determine MPI communicator
                // ==========================
                if(Target.MPI_Comm != this.MPI_Comm && Target.MPI_Comm != csMPI.Raw._COMM.SELF) {
                    throw new ArgumentException("todo");
                }
                MPI_Comm comm = Target.MPI_Comm;
                MPICollectiveWatchDog.Watch(comm);

                Dictionary<int, long[,]> ExternalColIndices;
                if(comm != csMPI.Raw._COMM.SELF)
                    ExternalColIndices = DetermineExternColumnIndices(Target, ColumnIndicesSource, ColIndicesTarget, comm);
                else
                    ExternalColIndices = new Dictionary<int, long[,]>();

                // define permutation of column indices
                // ====================================

                long[] LocalColIndexPermutation; // mapping: [local column index into this matrix (index)] ----> [global index into target matrix (content)]
                {
                    int L = this.ColPartition.LocalLength;
                    long iCol0 = this.ColPartition.i0;
                    LocalColIndexPermutation = new long[L];
                    ArrayTools.SetAll(LocalColIndexPermutation, long.MinValue);

                    IEnumerator<long> EnuColTarget = ColIndicesTarget != null ? ColIndicesTarget.GetEnumerator() : null;
                    int Counter = -1;
                    foreach(long iColSource in ColumnIndicesSource) {
                        m_ColPartitioning.TestIfInLocalRange(iColSource);
                        Counter++;
                        long iColTarget;
                        if(EnuColTarget != null) {
                            bool hasNext = EnuColTarget.MoveNext();
                            if(!hasNext)
                                throw new ArgumentException("Source and target column indices enumeration must have the same length, if provided.");
                            iColTarget = EnuColTarget.Current;
                        } else {
                            iColTarget = Counter + Target.ColPartition.i0;
                        }

                        LocalColIndexPermutation[iColSource - iCol0] = iColTarget;
                    }
                }
                //throw new NotImplementedException("todo");

                Dictionary<long, long> ExternalColIndexPermutation; // mapping: [global/external column index into this matrix (key)] ----> [global index into target matrix (value)]
                {
                    ExternalColIndexPermutation = new Dictionary<long, long>();

                    foreach(var kv in ExternalColIndices) {
                        int rank = kv.Key;
                        long[,] IdxEs = kv.Value;
                        int L = IdxEs.GetLength(0);
                        Debug.Assert(IdxEs.GetLength(1) == 2);
                        for(int l = 0; l < L; l++) {
                            ExternalColIndexPermutation.Add(IdxEs[l, 0], IdxEs[l, 1]);
                        }
                    }

                    if(ExternalColumnIndicesSource != null) {
                        IEnumerator<long> enuSrc = ExternalColumnIndicesSource.GetEnumerator();
                        IEnumerator<long> enuTrg = ExternalColIndicesTarget.GetEnumerator();
                        while(enuSrc.MoveNext()) {
                            bool t = enuTrg.MoveNext();
                            Debug.Assert(t == true);

                            long iSrc = enuSrc.Current;
                            long iTrg = enuTrg.Current;
                            if(m_ColPartitioning.IsInLocalRange(iSrc) == true)
                                throw new ArgumentException("External column indices are expected to be external.");

                            ExternalColIndexPermutation.Add(iSrc, iTrg);
                        }
                        Debug.Assert(enuTrg.MoveNext() == false);

                    }
                }


                // get sub-matrix
                // ==============
                double[] CopyBufferArray = null;
                MultidimensionalArray CopyBuffer = null;
                long Colj0 = m_ColPartitioning.i0;
                long ColJE = m_ColPartitioning.iE;
                int TargetRowCounter = 0;
                using(IEnumerator<long> RowEnu = RowIndicesSource.GetEnumerator(), RowTargEnu = RowIndicesTarget != null ? RowIndicesTarget.GetEnumerator() : null) {
                    bool SkipRowEnu = false;
                    while(SkipRowEnu || RowEnu.MoveNext()) {
                        if(SkipRowEnu == false) {
                            TargetRowCounter++;
                            if(RowTargEnu != null) {
                                bool bt = RowTargEnu.MoveNext();
                                Debug.Assert(bt);
                            }
                        }
                        SkipRowEnu = false;

                        long i0ThisCopy = RowEnu.Current;
                        long iNext = i0ThisCopy;
                        int RowCopyLen = 1;


                        long iBlk; //     global (over all MPI processes) block index.
                        int iBlkLoc; //   local (on this MPI process) block index.
                        long i0; //       Global start index of the block
                        int iLoc; //      Index within block.
                        int BlkT; //      Block type, see <see cref="BlockPartitioning.GetBlockType(int)"/>, resp. <see cref="BlockPartitioning.Subblk_i0"/>, <see cref="BlockPartitioning.SubblkLen"/>.
                        int Sblk_idx; //  Sub block index.
                        int i0_Sblk; //   First block index in sub block.
                        int ISblk; //     Size/length of sub block.
                        int NoOfSblk; //  Number of sub-Blocks.
                        int RemVoidRs; // in void region, the number of remaining void rows
                        TranslateIndex(i0ThisCopy, m_RowPartitioning, false, out iBlk, out iBlkLoc, out i0, out iLoc, out BlkT, out Sblk_idx, out i0_Sblk, out ISblk, out NoOfSblk, out RemVoidRs);
                        if(i0_Sblk < 0) {
                            // void-region
                            continue;
                        }

                        int MaxRowBlockSize = m_RowPartitioning.GetBlockLen(iBlk) - iLoc;
                        Debug.Assert(MaxRowBlockSize > 0);

                        long i0Copy;
                        if(RowIndicesTarget == null) {
                            i0Copy = TargetRowCounter - 1 + Target.RowPartitioning.i0;
                        } else {
                            i0Copy = RowTargEnu.Current;
                        }

                        while(true) { // scan row numbers, but at max. to the end of the block
                            if(RowCopyLen >= MaxRowBlockSize)
                                break;
                            bool reMN = RowEnu.MoveNext();
                            if(RowTargEnu != null) {
                                bool rteMN = RowTargEnu.MoveNext();
                                Debug.Assert(rteMN == reMN);
                            }
                            if(!reMN)
                                break;
                            else
                                TargetRowCounter++;
                            if((RowEnu.Current - iNext != 1)
                                || (RowTargEnu != null && RowTargEnu.Current - i0Copy != RowCopyLen)) {
                                SkipRowEnu = true;
                                break;
                            }
                            RowCopyLen++;
                            iNext = RowEnu.Current;
                        }

                        Debug.Assert(RowCopyLen >= 1);
                        Debug.Assert(m_RowPartitioning.IsInLocalRange(i0ThisCopy));
                        Debug.Assert(m_RowPartitioning.IsInLocalRange(i0ThisCopy + RowCopyLen - 1));
                        Debug.Assert(Target.RowPartitioning.IsInLocalRange(i0Copy));
                        Debug.Assert(Target.RowPartitioning.IsInLocalRange(i0Copy + RowCopyLen - 1));

                        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        // ** now, we are going to extract the rows i0Copy to i0Copy+RowCopyLen-1 **
                        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

                        SortedDictionary<long, BlockEntry> BlockRow = m_BlockRows[iBlkLoc];
                        if(BlockRow != null) {
                            foreach(var kv in BlockRow) {
                                long jBlk = kv.Key;
                                BlockEntry Blk = kv.Value;

                                // determine block range
                                long j0, JE;
                                if(m_ColPartitioning.FirstBlock <= jBlk && jBlk < m_ColPartitioning.FirstBlock + m_ColPartitioning.LocalNoOfBlocks) {
                                    j0 = m_ColPartitioning.GetBlockI0(jBlk);
                                    int BlockNoOfCols = m_ColPartitioning.GetBlockLen(jBlk);
                                    JE = j0 + BlockNoOfCols;
                                } else {
                                    GetExternalSubblockIndices(m_ColPartitioning, jBlk, out j0, out JE);
                                }

                                bool LocalBlock = m_ColPartitioning.IsInLocalRange(j0);
                                Debug.Assert(LocalBlock == (Colj0 <= j0 && j0 < ColJE));
                                Debug.Assert(LocalBlock == (Colj0 < JE && JE <= ColJE));

                                long j = j0;
                                while(j < JE) { // scan the column block range

                                    long j0ThisCopy = long.MinValue, j0Copy = long.MinValue;
                                    int ColCopyLen = -1;
                                    long jPrev = -1;
                                    while(j < JE) {
                                        long jTarg;
                                        if(LocalBlock) {
                                            jTarg = LocalColIndexPermutation[j - Colj0];
                                        } else {
                                            //Debug.Assert(ExternalColIndexPermutation.ContainsKey(j));
                                            //jTarg = ExternalColIndexPermutation[j];
                                            if(ExternalColIndexPermutation.TryGetValue(j, out jTarg) == false)
                                                jTarg = -132343;
                                        }


                                        if(jTarg < 0) {
                                            if(j0ThisCopy < 0) {
                                                j++;
                                                continue;
                                            } else {
                                                // end of continuous region reached.
                                                break;
                                            }
                                        }
                                        if(j0ThisCopy < 0) {
                                            // start a continuous region
                                            j0ThisCopy = j;
                                            j0Copy = jTarg;
                                            ColCopyLen = 1;
                                            jPrev = jTarg;
                                            j++;
                                            continue;
                                        } else {
                                            // currently in a block
                                            if(jTarg >= 0 && jTarg - jPrev == 1) {
                                                // continue with the block
                                                jPrev = jTarg;
                                                ColCopyLen++;
                                                j++;
                                                continue;
                                            } else {
                                                // end of continuous region reached.
                                                break;
                                            }
                                        }
                                    }

                                    if(j0ThisCopy >= 0) {
                                        // copy columns j0Copy to j0Copy + jCopyLen - 1
                                        Debug.Assert(ColCopyLen > 0);
#if DEBUG_EXTENDED
                                        for (long _j = j0ThisCopy; _j < j0ThisCopy + ColCopyLen; _j++) {
                                            if (LocalBlock) {
                                                Debug.Assert(LocalColIndexPermutation[_j - Colj0] >= 0);
                                                if (_j < j0ThisCopy + ColCopyLen - 1) {
                                                    Debug.Assert(LocalColIndexPermutation[_j - Colj0 + 1] - LocalColIndexPermutation[_j - Colj0] == 1);
                                                }
                                            } else {
                                                Debug.Assert(ExternalColIndexPermutation.ContainsKey(_j) && ExternalColIndexPermutation[_j] >= 0);
                                                if (_j < j0ThisCopy + ColCopyLen - 1) {
                                                    Debug.Assert(ExternalColIndexPermutation[_j + 1] - ExternalColIndexPermutation[_j] == 1);
                                                }
                                            }
                                        }
#endif
                                        // --------------------------------------------------------------------
                                        // finally, we are ready to copy a block:
                                        // row index offset:    this matrix:    i0ThisCopy
                                        //                      Target matrix:  i0Copy
                                        // number of rows:                      RowCopyLen
                                        // col index offset:    this matrix:    j0ThisCopy
                                        //                      Target matrix:  j0Copy
                                        // number of columns:                   ColCopyLen
                                        // --------------------------------------------------------------------

                                        Debug.Assert(RowCopyLen > 0);
                                        Debug.Assert(ColCopyLen > 0);
                                        Debug.Assert(i0 + iLoc == i0ThisCopy);
                                        if(RowCopyLen == 1 && ColCopyLen == 1) {
                                            // scalar version
                                            Target[i0Copy, j0Copy] += alpha * this[i0ThisCopy, j0ThisCopy];
                                        } else {
                                            // vectorized version
                                            if(CopyBufferArray == null || CopyBufferArray.Length < ColCopyLen * RowCopyLen) {
                                                CopyBufferArray = new double[ColCopyLen * RowCopyLen];
                                                CopyBuffer = null;
                                            }
                                            if(CopyBuffer == null || CopyBuffer.NoOfRows != RowCopyLen || CopyBuffer.NoOfCols != ColCopyLen)
                                                CopyBuffer = MultidimensionalArray.CreateWrapper(CopyBufferArray, RowCopyLen, ColCopyLen);

                                            this.ReadBlock(i0ThisCopy, j0ThisCopy, CopyBuffer);
                                            Target.AccBlock(i0Copy, j0Copy, alpha, CopyBuffer);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        private Dictionary<int, long[,]> DetermineExternColumnIndices<V3, V4>(IMutableMatrixEx Target, V3 ColumnIndicesSource, V4 ColIndicesTarget, MPI_Comm comm)
            where V3 : IEnumerable<long>
            where V4 : IEnumerable<long> //
        {
            // check input arguments
            // =====================
#if DEBUG_EXTENDED
                int iPrev = -1;
                
                foreach (int i in ColumnIndicesSource) {
                    this.m_ColPartitioning.TestIfInLocalRange(i);
                    //if (i <= iPrev)
                    //    throw new ArgumentException("Indices must be in strictly ascending order.");
                    iPrev = i;
                }
                if (ColIndicesTarget != null && ColIndicesTarget.Count() != ColumnIndicesSource.Count()) {
                    throw new ArgumentException("Mismatch between number of source and target indices.");
                }
#endif


            // determine which column indices have to be sent to other processors
            // ==================================================================                

            Dictionary<int, long[,]> ColumnIndices2Send = new Dictionary<int, long[,]>(); //   key: target processor.
                                                                                          // value: mapping [this column index] --> [target column index]
                                                                                          //   first index: enumeration; second index 0 for this, 1 for target.
            {
                this.UpdateCommPattern(this.MPI_Comm);

                // convert send list into array:
                List<long>[] SendColSource = new List<long>[SendLists.Count];
                List<long>[] SendColTarget = new List<long>[SendLists.Count];
                int[] TargetProc = new int[SendLists.Count];
                long[][,] SendRange = new long[SendLists.Count][,];
                int Counter2 = -1;
                foreach(var kv in this.SendLists) {
                    Counter2++;
                    TargetProc[Counter2] = kv.Key;
                    SendRange[Counter2] = kv.Value;
                    SendColSource[Counter2] = new List<long>();
                    SendColTarget[Counter2] = new List<long>();
                }
                int[] RangePointer = new int[SendLists.Count];


                // loop over the column index list...
                using(IEnumerator<long> EnuColTarget = ColIndicesTarget != null ? ColIndicesTarget.GetEnumerator() : null) {
                    int Counter = -1;
                    foreach(long iColSource in ColumnIndicesSource) {
                        Counter++;
                        long iColTarget;
                        if(EnuColTarget != null) {
                            bool hasNext = EnuColTarget.MoveNext();
                            if(!hasNext)
                                throw new ArgumentException("Source and target column indices enumeration must have the same length, if provided.");
                            iColTarget = EnuColTarget.Current;
                        } else {
                            iColTarget = Counter + Target.ColPartition.i0;
                        }

                        m_ColPartitioning.TestIfInLocalRange(iColSource);
                        Target.ColPartition.TestIfInLocalRange(iColTarget);

                        // loop over all possible target processors...
                        bool AllFinish = true;
                        for(int i = 0; i < TargetProc.Length; i++) {
                            long[,] _SendRange = SendRange[i];
                            if(RangePointer[i] >= _SendRange.GetLength(0))
                                continue;

                            long Send_i0 = _SendRange[RangePointer[i], 0];
                            long Send_iE = _SendRange[RangePointer[i], 1];
                            Debug.Assert(m_ColPartitioning.IsInLocalRange(Send_i0));
                            Debug.Assert(Send_iE - Send_i0 >= 1);
                            Debug.Assert(m_ColPartitioning.IsInLocalRange(Send_iE - 1));
                            bool bFinished = false;
                            while(Send_iE <= iColSource) {
                                RangePointer[i]++;
                                if(RangePointer[i] >= _SendRange.GetLength(0)) {
                                    bFinished = true;
                                    break;
                                }
                                Send_i0 = _SendRange[RangePointer[i], 0];
                                Send_iE = _SendRange[RangePointer[i], 1];
                                Debug.Assert(m_ColPartitioning.IsInLocalRange(Send_i0));
                                Debug.Assert(Send_iE - Send_i0 >= 1);
                                Debug.Assert(m_ColPartitioning.IsInLocalRange(Send_iE - 1));
                            }

                            if(bFinished)
                                continue;
                            AllFinish = false;

                            Debug.Assert(iColSource < Send_iE);
                            if(Send_i0 <= iColSource) {
                                SendColSource[i].Add(iColSource);
                                SendColTarget[i].Add(iColTarget);
                            }
                        }

                        if(AllFinish)
                            break;
                    }
                }

                // finalize send list assembly
                for(int i = 0; i < TargetProc.Length; i++) {
                    Debug.Assert(SendColSource[i].Count == SendColTarget[i].Count);
                    if(SendColSource[i].Count > 0) {
                        long[,] data = new long[SendColSource[i].Count, 2];
                        data.SetColumn(SendColSource[i], 0);
                        data.SetColumn(SendColTarget[i], 1);
                        ColumnIndices2Send.Add(TargetProc[i], data);
                    }
                }
            }
            //foreach (var kv in ColumnIndices2Send) {
            //    int Trnk = kv.Key;
            //    int NoOf = kv.Value.GetLength(0);
            //    Console.WriteLine("P{0} to P{1}: #itms {2}", m_RowPartitioning.MpiRank, Trnk, NoOf);
            //}
            //if(ColumnIndices2Send.Count == 0)
            //    Console.WriteLine("P{0}: nothing to send, dork.", m_RowPartitioning.MpiRank);


            // MPI exchange of lists
            // =====================
            Dictionary<int, long[,]> ExternalColIndices; //   key: MPI processor rank
                                                         // value: mapping [this column index] --> [target column index]
                                                         //        first index: enumeration; second index 0 for this, 1 for target.
            IntBufferExchange(comm, ColumnIndices2Send, out ExternalColIndices);
            ColumnIndices2Send = null;
#if DEBUG_EXTENDED
            foreach (var kv in ExternalColIndices) {
                int rank = kv.Key;
                long[,] IdxEs = kv.Value;
                int L = IdxEs.GetLength(0);
                Debug.Assert(IdxEs.GetLength(1) == 2);
                for (int l = 0; l < L; l++) {
                    Debug.Assert(this.m_ColPartitioning.FindProcess(IdxEs[l, 0]) == rank);
                    Debug.Assert(Target.ColPartition.FindProcess(IdxEs[l, 1]) == rank);
                }
            }
#endif
            return ExternalColIndices;
        }

        /// <summary>
        /// Computes the transpose of this matrix.
        /// </summary>
        public BlockMsrMatrix Transpose() {
            using(new FuncTrace()) {
                //throw new NotImplementedException();
                BlockMsrMatrix Ret = new BlockMsrMatrix(this._ColPartitioning, this._RowPartitioning);

                Debug.Assert(this.m_BlockRows.Length == _RowPartitioning.LocalNoOfBlocks);
                int LocalNoOfBlocks = _RowPartitioning.LocalNoOfBlocks;
                long FirstBlock = _RowPartitioning.FirstBlock;

                // local transpose
                // ===============

                Dictionary<int, List<Tuple<long, long, MultidimensionalArray>>> ExchangeData = new Dictionary<int, List<Tuple<long, long, MultidimensionalArray>>>();

                MultidimensionalArray temp = null, tempT = null;
                for(int iBlockRowLoc = 0; iBlockRowLoc < LocalNoOfBlocks; iBlockRowLoc++) {
                    long iBlockRow = FirstBlock + iBlockRowLoc;
                    long i0 = _RowPartitioning.GetBlockI0(iBlockRow);
                    int I = _RowPartitioning.GetBlockLen(iBlockRow);
                    var Row = m_BlockRows[iBlockRowLoc];
                    if(Row != null && I > 0) {
                        foreach(var kv in Row) {
                            long jBlkCol = kv.Key;
                            BlockEntry BE = kv.Value;
                            Debug.Assert(BE.jBlkCol == jBlkCol);
                            int J;
                            long j0;
                            if(m_ColPartitioning.IsLocalBlock(jBlkCol)) {
                                J = _ColPartitioning.GetBlockLen(jBlkCol);
                                j0 = _ColPartitioning.GetBlockI0(jBlkCol);
                            } else {
                                GetExternalSubblockIndices(m_ColPartitioning, jBlkCol, out j0, out long _J);
                                J = checked((int)(_J - j0));
                            }
                            if(J > 0) {

                                if(temp == null || temp.GetLength(0) != I || temp.GetLength(1) != J) {
                                    temp = MultidimensionalArray.Create(I, J);
                                }
                                if(tempT == null || tempT.GetLength(0) != J || tempT.GetLength(1) != I) {
                                    tempT = MultidimensionalArray.Create(J, I);
                                }

                                this.ReadBlock(i0, j0, temp);
                                temp.TransposeTo(tempT);

                                if(_ColPartitioning.IsLocalBlock(jBlkCol)) {
                                    Ret.AccBlock(j0, i0, 1.0, tempT);
                                } else {
                                    int TargProc = _ColPartitioning.FindProcessForBlock(jBlkCol);
                                    List<Tuple<long, long, MultidimensionalArray>> ExchangeData_TargProc;
                                    if(!ExchangeData.TryGetValue(TargProc, out ExchangeData_TargProc)) {
                                        ExchangeData_TargProc = new List<Tuple<long, long, MultidimensionalArray>>();
                                        ExchangeData.Add(TargProc, ExchangeData_TargProc);
                                    }

                                    Debug.Assert(_ColPartitioning.FindProcess(j0) == TargProc);
                                    Debug.Assert(_ColPartitioning.FindProcess(j0 + J - 1) == TargProc);
                                    ExchangeData_TargProc.Add(new Tuple<long, long, MultidimensionalArray>(j0, i0, tempT));
                                    tempT = null;
                                }
                            }
                        }
                    }
                }

                // exchange
                // ========
                var RcvData = SerialisationMessenger.ExchangeData(ExchangeData, this.MPI_Comm);

                // received part
                // =============

                foreach(var vl in RcvData.Values) {
                    foreach(var t in vl) {
                        long i0 = t.Item1;
                        long j0 = t.Item2;
                        MultidimensionalArray data = t.Item3;
                        Ret.AccBlock(i0, j0, 1.0, data);
                    }
                }

                // return
                // ======
                return Ret;

            }
        }


        /// <summary>
        /// Computes the deviation of this matrix from symmetry
        /// </summary>
        /// <returns>
        /// The accumulated sum of the differences between corresponding
        /// off-diagonal entries.
        /// </returns>
        public double SymmetryDeviation() {
            var MT = this.Transpose();
            MT.Acc(-1.0, this);
            double sd = MT.InfNorm();
            return sd;
        }

        /// <summary>
        /// Performs the operation: this = this + <paramref name="Ascale"/>*<paramref name="A"/>;
        /// </summary>
        /// <param name="A">
        /// another matrix with same size and equal <see cref="RowPartitioning"/>;
        /// </param>
        /// <param name="Ascale">
        /// scaling
        /// </param>
        public void Acc(double Ascale, BlockMsrMatrix A) {
            if(m_RowPartitioning.LocalLength != A._RowPartitioning.LocalLength)
                throw new ArgumentException("Mismatch in number of rows.");
            if(m_ColPartitioning.TotalLength != A._ColPartitioning.TotalLength)
                throw new ArgumentException("Mismatch in number of columns.");

            long iBlk0 = A._RowPartitioning.FirstBlock;
            long iBlkE = A._RowPartitioning.LocalNoOfBlocks + iBlk0;
            for(long iBlk = iBlk0; iBlk < iBlkE; iBlk++) {
                var BlockRow = A.m_BlockRows[iBlk - iBlk0];
                if(BlockRow != null) {

                    long i0 = A._RowPartitioning.GetBlockI0(iBlk);
                    int BT = A._RowPartitioning.GetBlockType(iBlk);
                    int[] SblkOffset = A._RowPartitioning.GetSubblk_i0(BT);
#if DEBUG_EXTENDED
                    int[] SblkLen = A._RowPartitioning.GetSubblkLen(BT);
#endif

                    foreach(var kv in BlockRow) {
                        long jBlk = kv.Key;
                        BlockEntry B = kv.Value;

                        long j0;
                        int[] ColSblkOffset;
#if DEBUG_EXTENDED
                        int[] ColSblkLen;
#endif

                        if(jBlk >= A._ColPartitioning.FirstBlock && jBlk < A._ColPartitioning.LocalNoOfBlocks + A._ColPartitioning.FirstBlock) {
                            j0 = A._ColPartitioning.GetBlockI0(jBlk);
                            int CBT = A._ColPartitioning.GetBlockType(jBlk);
                            ColSblkOffset = A._ColPartitioning.GetSubblk_i0(CBT);
#if DEBUG_EXTENDED
                            ColSblkLen = A._ColPartitioning.GetSubblkLen(CBT);
#endif
                        } else {
                            long aa, bb;
                            BlockMsrMatrix.GetExternalSubblockIndices(A._ColPartitioning, jBlk, out aa, out bb);
                            ColSblkOffset = new int[] { 0 };
                            j0 = aa;
#if DEBUG_EXTENDED
                            ColSblkLen = new int[] { (int)(bb - aa) };
#endif
                        }

                        Debug.Assert(B.InMembnk.GetLength(0) == B.MembnkIdx.GetLength(0));
                        Debug.Assert(B.InMembnk.GetLength(1) == B.MembnkIdx.GetLength(1));

                        int NoOfRowSblk = B.InMembnk.GetLength(0);
                        int NoOfColSblk = B.InMembnk.GetLength(1);

                        for(int iSblk = 0; iSblk < NoOfRowSblk; iSblk++) {
                            for(int jSblk = 0; jSblk < NoOfColSblk; jSblk++) {
                                int InMembnk = B.InMembnk[iSblk, jSblk];
                                int MembnkIdx = B.MembnkIdx[iSblk, jSblk];
                                Debug.Assert((InMembnk >= 0) == (MembnkIdx >= 0));

                                if(InMembnk >= 0) {
                                    var Block = A.m_Membanks[MembnkIdx].Mem.ExtractSubArrayShallow(InMembnk, -1, -1);
#if DEBUG_EXTENDED
                                    Debug.Assert(Block.NoOfRows == SblkLen[iSblk]);
                                    Debug.Assert(Block.NoOfCols == ColSblkLen[jSblk]);

#endif
                                    this.AccBlock(i0 + SblkOffset[iSblk], j0 + ColSblkOffset[jSblk], Ascale, Block);
                                }
                            }
                        }

                    }
                }
            }
        }

        /// <summary>
        /// multiplies this matrix by the factor <paramref name="a"/>
        /// </summary>
        /// <param name="a"></param>
        public void Scale(double a) {
            if(a == 0.0) {
                this.Clear();
            } else {
                if(a == 1.0)
                    return;
                foreach(var Membnk in this.m_Membanks) {
                    Membnk.Mem.Scale(a);
                }
            }
        }

        /// <summary>
        /// matrix-matrix-multiplication, see <see cref="Multiply(BlockMsrMatrix, BlockMsrMatrix)"/>.
        /// </summary>
        static public BlockMsrMatrix operator *(BlockMsrMatrix left, BlockMsrMatrix right) {
            BlockMsrMatrix Res = new BlockMsrMatrix(left._RowPartitioning, right._ColPartitioning);
            BlockMsrMatrix.Multiply(Res, left, right);
            return Res;
        }

        /// <summary>
        /// Multiply matrix <paramref name="M"/> by a scalar <paramref name="a"/>.
        /// </summary>
        static public BlockMsrMatrix operator *(double a, BlockMsrMatrix M) {
            BlockMsrMatrix R = M.CloneAs();
            R.Scale(a);
            return R;
        }

        /// <summary>
        /// Multiply matrix <paramref name="M"/> by a scalar <paramref name="a"/>.
        /// </summary>
        static public BlockMsrMatrix operator *(BlockMsrMatrix M, double a) {
            BlockMsrMatrix R = M.CloneAs();
            R.Scale(a);
            return R;
        }

        /// <summary>
        /// summation of matrices, see <see cref="Acc"/>.
        /// </summary>
        static public BlockMsrMatrix operator +(BlockMsrMatrix left, BlockMsrMatrix right) {
            var ReturnMatrix = left.CloneAs();
            ReturnMatrix.Acc(1.0, right);
            return ReturnMatrix;
        }

        /// <summary>
        /// Difference of matrices, see <see cref="Acc"/>.
        /// </summary>
        static public BlockMsrMatrix operator -(BlockMsrMatrix left, BlockMsrMatrix right) {
            var ReturnMatrix = left.CloneAs();
            ReturnMatrix.Acc(-1.0, right);
            return ReturnMatrix;
        }

        /// <summary>
        /// total number of columns (over all MPI processors)
        /// </summary>
        public long NoOfCols {
            get {
                return m_ColPartitioning.TotalLength;
                //return m_Size; 
            }
        }

        /// <summary>
        /// total number of rows (over all MPI processors)
        /// </summary>
        public long NoOfRows {
            get {
                return m_RowPartitioning.TotalLength;
            }
        }

        void ClearBlockRow(long _iBlockRow) {
            long bi0 = this._RowPartitioning.FirstBlock;
            int NoB = this._RowPartitioning.LocalNoOfBlocks;
            int iBlockRow = (int)(_iBlockRow - bi0);
            if(iBlockRow < 0 || iBlockRow >= NoB)
                throw new IndexOutOfRangeException("Block row out of range!");

            var Row = this.m_BlockRows[iBlockRow];
            if(Row != null) {
                foreach(var kv in Row) {
                    BlockEntry Blk = kv.Value;
                    Debug.Assert(Blk.jBlkCol == kv.Key);

                    Debug.Assert(Blk.InMembnk.GetLength(0) == Blk.MembnkIdx.GetLength(0));
                    Debug.Assert(Blk.InMembnk.GetLength(1) == Blk.MembnkIdx.GetLength(1));
                    int IB = Blk.InMembnk.GetLength(0);
                    int JB = Blk.InMembnk.GetLength(1);

                    for(int ib = 0; ib < IB; ib++) {
                        for(int jb = 0; jb < JB; jb++) {
                            int inmembnk = Blk.InMembnk[ib, jb];
                            int membnkidx = Blk.MembnkIdx[ib, jb];
                            Debug.Assert((inmembnk < 0) == (membnkidx < 0));

                            if(membnkidx >= 0) {
                                Debug.Assert(this.m_Membanks[membnkidx].Mem.Dimension == 3);
                                Debug.Assert(this.m_Membanks[membnkidx].Mem.GetLength(0) == this.m_Membanks[membnkidx].Occupied.Length);
                                FreeSblk(membnkidx, inmembnk);
                            }

                        }
                    }

                    /*
                    bool isExternal = (Blk.jBlkCol < m_ColPartitioning.FirstBlock) || (m_ColPartitioning.LocalNoOfBlocks + m_ColPartitioning.FirstBlock <= Blk.jBlkCol);

                    if(isExternal) {
                        int OwnerProc = m_ColPartitioning.FindProcessForBlock(Blk.jBlkCol);
                        Debug.Assert(OwnerProc != m_ColPartitioning.MpiRank);
                        HashSet<int> ExtIdx = m_ExternalBlockIndicesByProcessor[OwnerProc];
                        Debug.Assert(ExtIdx.Contains(Blk.jBlkCol));
                        ExtIdx.Remove(Blk.jBlkCol);
                    }
                    */

                    // Reminder: the code above is wrong, for the following reason:
                    // 'm_ExternalBlockIndicesByProcessor[OwnerProc]' is a collection of external blocks owned by 'OwnerProc'
                    // for **all** rows -- we would have to check all other rows in order to safely remove 
                    // 'BE.jBlkCol' from 'm_ExternalBlockIndicesByProcessor[OwnerProc]'.
                    // On the other hand, having unused external block in the set may impact MPI communication, but 
                    // does not impact the numerical results.

                }


                this.m_BlockRows[iBlockRow].Clear();
                this.m_ExternalBlock[iBlockRow] = false;
            }
        }


        /*
        /// <summary>
        /// sets row number <paramref name="i"/>;
        /// All previous entries in this row are overwritten;
        /// </summary>
        /// <param name="i"> row index in global indices</param>
        /// <param name="row">
        /// </param>
        public void SetRow(int i, params MsrMatrix.MatrixEntry[] row) {
            throw new NotImplementedException();
        }
        

        /// <summary>
        /// Sets all entries in a row to 0
        /// </summary>
        /// <param name="i">row index in global indices</param>
        public void ClearRow(int i) {
            throw new NotImplementedException();
        }
        */

        /// <summary>
        /// Allocates memory for some entry, even if the entry is zero;
        /// </summary>
        /// <param name="i">global (over all MPI processes) row index</param>
        /// <param name="j">global (over all MPI processes) column index</param>
        public void AllocateEvenWhenZero(int i, int j) {
            double Mij = this[i, j];
            if(Mij == 0.0) {
                this[i, j] = 1.0;
                this[i, j] = 0.0;
            }
        }

        #region IMutuableMatrix Members

        /// <summary>
        /// see <see cref="IMutableMatrix.GetValues"/>;
        /// </summary>
        public double[] GetValues(long RowIndex, long[] ColumnIndices) {
            double[] ret = new double[ColumnIndices.Length];
            for(int i = 0; i < ret.Length; i++) {
                ret[i] = this[RowIndex, ColumnIndices[i]];
            }
            return ret;
        }

        /// <summary>
        /// see <see cref="IMutableMatrix.SetValues"/>;
        /// </summary>
        public void SetValues(long RowIndex, long[] ColumnIndices, double[] newValues) {
            if(newValues.Length != ColumnIndices.Length) {
                throw new ArgumentException("Mismatch in array length; column index array and value array must have the same length.");
            }
            int L = ColumnIndices.Length;
            for(int i = 0; i < L; i++) {
                this[RowIndex, ColumnIndices[i]] = newValues[i];
            }

        }

        /// <summary>
        /// always true, see <see cref="IMutableMatrix.OccupationMutable"/>;
        /// </summary>
        public bool OccupationMutable {
            get {
                return true;
            }
        }

        #endregion

        #region ISparseMatrix Members

        /// <summary>
        /// see <see cref="ISparseMatrix.GetDiagonalElement"/>;
        /// </summary>
        public double GetDiagonalElement(long row) {
            return this[row, row];
        }

        /// <summary>
        /// see <see cref="ISparseMatrix.SetDiagonalElement"/>;
        /// </summary>
        public void SetDiagonalElement(long row, double val) {
            this[row, row] = val;
        }

        #endregion

        #region IMutuableMatrixEx Members

        /// <summary>
        /// see <see cref="IMutableMatrixEx.GetOccupiedColumnIndices"/>
        /// </summary>
        public int GetOccupiedColumnIndices(long i, ref long[] ColumnIndices) {
            double[] dummy = null;
            return GetRowInternal(i, ref ColumnIndices, ref dummy, true);
        }

        #endregion

        /// <summary>
        /// Returns the block-column index of all non-zero blocks in block-row <paramref name="iBlk"/>
        /// </summary>
        /// <param name="alsoExternal">
        /// True: also return block indices in the external range 
        /// (<see cref="IBlockPartitioning.FirstBlock"/>, <see cref="IBlockPartitioning.LocalNoOfBlocks"/>), 
        /// i.e. coupling with other MPI processors
        /// </param>
        /// <param name="iBlk">row block index in the local range</param>
        public long[] GetOccupiedRowBlockIndices(long iBlk, bool alsoExternal = false) {
            var part = this._RowPartitioning;
            if(part.GetBlockLen(iBlk) <= 0)
                return new long[0];

            int iBlkLoc = (int)(iBlk - part.FirstBlock);

            if(alsoExternal == true)
                throw new NotSupportedException("Todo; I'm not sure how to handle this; ");

            var RowDict = this.m_BlockRows[iBlkLoc];
            if(RowDict == null)
                return new long[0];

            List<long> ret = new List<long>();
            foreach(var kv in RowDict) {
                long jBlk = kv.Key;
                if(alsoExternal == false) {
                    if(this._ColPartitioning.IsLocalBlock(jBlk) == false)
                        continue;
                }

                BlockEntry block = kv.Value;
                if(!block.IsEmpty) {
                    ret.Add(jBlk);
                }

            }

            return ret.ToArray();
        }


        #region MatrixMatrixMult

        struct Multiply_Helper_1 {
            public long i0;
            public long j0;
            public int MembankIdx;
            public int InMembank;
        }

        /// <summary>
        /// Utility class used by <see cref="Multiply(BlockMsrMatrix, BlockMsrMatrix, BlockMsrMatrix)"/>;
        /// basically, the implementation of the merge-algorithm which is crucial for the performance of the 
        /// sparse matrix-matrix multiplication.
        /// </summary>
        class TempBlockRow : IDisposable {
            public long[] m_j0S;
            public long[] m_JES;
            public int[] m_Offset;
            public double[] m_mem;
            public IntPtr m_pMem;
            GCHandle m_memLock;
            public int count = 0;
            public long i0;
            public long iE;

            /// <summary>
            /// The actually used number of bytes of this object in memory.
            /// </summary>
            public int GetSize() {
                return sizeof(long) * 2 + sizeof(int) // i0, ie, count 
                    + (sizeof(long) * 2 + sizeof(int)) * count // m_j0S, m_JES, m_Offset
                    + sizeof(double) * MemEntries;
            }

            int MemEntries {
                get {
                    if(count > 0)
                        return checked((int)(m_Offset[count - 1] + (iE - i0) * (m_JES[count - 1] - m_j0S[count - 1])));
                    else
                        return 0;
                }
            }

            /// <summary>
            /// Serializes this object into an unmanaged buffer, starting at address <paramref name="pBuffer"/>;
            /// the number of bytes is determined by <see cref="GetSize"/>.
            /// </summary>
            public IntPtr Serialize(IntPtr pBuffer) {
                unsafe {
                    long* pL = (long*)pBuffer;
                    *pL = i0;
                    pL++;
                    *pL = iE;
                    pL++;
                    int* pI = (int*)pL; pL = null;
                    *pI = count;
                    pI++;
                    pL = (long*)pI; pI = null;

                    Debug.Assert((count == 0 && m_j0S == null) || (m_j0S.Length >= count));
                    Debug.Assert((count == 0 && m_JES == null) || (m_JES.Length >= count));
                    Debug.Assert((count == 0 && m_Offset == null) || (m_Offset.Length >= count + 1));

                    fixed(long* p_j0S = m_j0S, p_jES = m_JES) {
                        fixed(int* p_Offset = m_Offset) {
                            long* _j0S = p_j0S, _jES = p_jES;
                            int* _Offset = p_Offset;
                            for(int l = count; l > 0; l--) {
                                *pL = *_j0S;
                                pL++;
                                _j0S++;
                            }
                            for(int l = count; l > 0; l--) {
                                *pL = *_jES;
                                pL++;
                                _jES++;
                            }
                            pI = (int*)pL; pL = null;
                            for(int l = count; l > 0; l--) {
                                *pI = *_Offset;
                                pI++;
                                _Offset++;
                            }
                        }
                    }

                    double* pD = (double*)pI;
                    double* _pM = (double*)m_pMem;

                    for(int l = MemEntries; l > 0; l--) {
                        *pD = *_pM;
                        pD++;
                        _pM++;
                    }

                    IntPtr r = (IntPtr)pD;
                    Debug.Assert(((byte*)r - (byte*)pBuffer) == this.GetSize());
                    return r;
                }
            }

            /// <summary>
            /// De-serializes a block row from an unmanaged buffer, starting at address <paramref name="pBuffer"/>;
            /// </summary>
            public IntPtr Deserialize(IntPtr pBuffer) {
                unsafe {
                    long* pL = (long*)pBuffer;
                    this.i0 = *pL;
                    pL++;
                    this.iE = *pL;
                    pL++;
                    int* pI = (int*)pL;
                    this.count = *pI; pL = null;
                    pI++;
                    pL = (long*)pI; pI = null;

                    int Capacity = this.m_j0S != null ? this.m_j0S.Length : 0;
                    if(count >= Capacity) {
                        Array.Resize(ref m_j0S, count + 5);
                        Array.Resize(ref m_JES, count + 5);
                        Array.Resize(ref m_Offset, count + 5);
                    }

                    fixed(long* p_j0S = m_j0S, p_jES = m_JES) {
                        fixed(int* p_Offset = m_Offset) {

                            long* _j0S = p_j0S, _jES = p_jES;
                            int* _Offset = p_Offset;
                            for(int l = count; l > 0; l--) {
                                *_j0S = *pL;
                                pL++;
                                _j0S++;
                            }
                            for(int l = count; l > 0; l--) {
                                *_jES = *pL;
                                pL++;
                                _jES++;
                            }
                            pI = (int*)pL; pL = null;
                            for(int l = count; l > 0; l--) {
                                *_Offset = *pI;
                                pI++;
                                _Offset++;
                            }
                        }
                    }

                    ResizeMem(this.MemEntries, 0);
                    Debug.Assert(m_mem.Length >= m_Offset[count]);


                    double* pD = (double*)pI;
                    double* _pM = (double*)m_pMem;

                    for(int l = MemEntries; l > 0; l--) {
                        *_pM = *pD;
                        pD++;
                        _pM++;
                    }

                    IntPtr r = (IntPtr)pD;
                    Debug.Assert(((byte*)r - (byte*)pBuffer) == this.GetSize());
                    return r;
                }
            }

            public void Dispose() {
                Debug.Assert((m_mem != null) == (m_pMem != IntPtr.Zero));

                if(m_mem != null) {
                    Debug.Assert(m_pMem == Marshal.UnsafeAddrOfPinnedArrayElement(m_mem, 0));
                    m_memLock.Free();
                    m_pMem = IntPtr.Zero;
                }
            }

            void ResizeMem(int minSz, int EnlargeSz) {
                if(m_mem != null && m_mem.Length >= minSz)
                    return;

                Debug.Assert((m_mem != null) == (m_pMem != IntPtr.Zero));
                if(m_mem != null) {
                    Debug.Assert(m_pMem == Marshal.UnsafeAddrOfPinnedArrayElement(m_mem, 0));
                    m_memLock.Free();
                }

                int CurSz = m_mem == null ? 0 : m_mem.Length;
                int NewSz = Math.Max(CurSz + EnlargeSz, minSz + 1);


                Array.Resize(ref m_mem, NewSz);
                m_memLock = GCHandle.Alloc(m_mem, GCHandleType.Pinned);
                m_pMem = Marshal.UnsafeAddrOfPinnedArrayElement(m_mem, 0);
            }

            public void AccessBlock(int jBlock, out int _Offset, out int _CI) {
                Debug.Assert(jBlock < count);
                _Offset = m_Offset[jBlock];
                _CI = (int)(m_JES[jBlock] - m_j0S[jBlock]);
            }

            void AllocBlock(long j0, long jE, out int _Offset, out int _CI) {
                Debug.Assert(jE > j0);
                Debug.Assert(this.iE > this.i0);
                Debug.Assert((m_j0S == null) == (m_JES == null));
                Debug.Assert((m_j0S == null) || (m_j0S.Length == m_JES.Length));

                int Capacity = m_j0S != null ? m_j0S.Length : 0;
                if(count >= Capacity - 1) {
                    Array.Resize(ref m_j0S, m_j0S == null ? 10 : m_j0S.Length * 2);
                    Array.Resize(ref m_JES, m_JES == null ? 10 : m_JES.Length * 2);
                    Array.Resize(ref m_Offset, m_Offset == null ? 10 : m_Offset.Length * 2);
                }

                int Offset;
                if(count > 0) {
                    Debug.Assert(j0 >= m_JES[count - 1]);
                    Offset = (int)(m_Offset[count - 1] + (m_JES[count - 1] - m_j0S[count - 1]) * (iE - i0));
                } else {
                    Offset = 0;
                }
                m_Offset[count] = Offset;

                m_j0S[count] = j0;
                m_JES[count] = jE;
                count++;

                int ReqMem = checked((int)((iE - i0) * (jE - j0)));
                bool bResize = (m_mem == null || m_mem.Length < Offset + ReqMem);
                if(bResize) {
                    ResizeMem(Offset + ReqMem, 5 * ReqMem);
                }
                Array.Clear(m_mem, Offset, ReqMem);

                _Offset = Offset;
                _CI = (int)(jE - j0);

#if DEBUG_EXTENDED
                for (int i = 0; i < (this.iE - this.i0); i++) {
                    for (int j = 0; j < (jE - j0); j++) {
                        Debug.Assert(m_mem[i * _CI + j + _Offset] == 0.0);
                        unsafe
                        {
                            Debug.Assert(((double*)(m_pMem))[i * _CI + j + _Offset] == 0.0);
                        }
                    }
                }
#endif

                //return Resize;
            }

            void EnlargeLastBlock(long new_JE, out int _Offset, out int _Ci) {
                Debug.Assert(this.count > 0);
                Debug.Assert(this.m_mem != null);
                int cc = this.count - 1;
                Debug.Assert(new_JE > this.m_JES[cc]);

                int Offset = m_Offset[cc];
                long j0 = m_j0S[cc];
                int ReqMem = checked((int)((iE - i0) * (new_JE - j0)));
                ResizeMem(Offset + ReqMem, 4 * ReqMem);

                long old_JE = m_JES[cc];
                int new_Ci = (int)(new_JE - j0);
                int old_Ci = (int)(old_JE - j0);
                _Ci = new_Ci;
                _Offset = Offset;

                for(int i = (int)((iE - i0) - 1); i >= 0; i--) {
                    for(int j = (int)((new_JE - j0) - 1); j >= 0; j--) {
                        int idx_new = Offset + new_Ci * i + j;
                        if(j >= old_Ci) {
                            m_mem[idx_new] = 0;
                        } else {
                            int idx_old = Offset + old_Ci * i + j;
                            m_mem[idx_new] = m_mem[idx_old];
                        }
                    }
                }

                m_JES[cc] = new_JE;
            }

            public void Clear(long new_i0, long new_iE) {
                count = 0;
                i0 = new_i0;
                iE = new_iE;
            }

            public void Save(BlockMsrMatrix BM, long iBlock) {
                long i0_Block = BM._RowPartitioning.GetBlockI0(iBlock);
                long iE_Block = BM._RowPartitioning.GetBlockLen(iBlock) + i0_Block;
                //if (i0_Block > i0)
                //    throw new NotSupportedException();
                //if (iE_Block < iE)
                //    throw new NotSupportedException();

                int RowBtüpe = BM._RowPartitioning.GetBlockType(iBlock);
                int[] RowSblk_i0 = BM._RowPartitioning.GetSubblk_i0(RowBtüpe);
                int[] RowSblkLen = BM._RowPartitioning.GetSubblkLen(RowBtüpe);
                Debug.Assert(RowSblk_i0.Length == RowSblkLen.Length);
                int NoOfRowSblk = RowSblk_i0.Length;

                double[] srcMem = this.m_mem;

                int[] dummy_i0 = new int[1];
                int[] dummy_iE = new int[1];

                for(int _c = 0; _c < this.count; _c++) { // loop over blocks in temporary row (row loop)
                    long _j0 = this.m_j0S[_c];
                    long _jE = this.m_JES[_c];
                    int _offset = this.m_Offset[_c];
                    int _CI = (int)(_jE - _j0);

                    long jPointer = _j0;
                    while(jPointer < _jE) {
                        int[] ColSblk_i0, ColSblkLen;
                        long jBlockCol, j0_Block, jE_Block;
                        if(BM.ColPartition.IsInLocalRange(jPointer)) {

                            jBlockCol = BM._ColPartitioning.GetBlockIndex(jPointer);
                            j0_Block = BM._ColPartitioning.GetBlockI0(jBlockCol);
                            jE_Block = BM._ColPartitioning.GetBlockLen(jBlockCol) + j0_Block;
                            int ColBtüpe = BM._ColPartitioning.GetBlockType(jBlockCol);
                            ColSblk_i0 = BM._ColPartitioning.GetSubblk_i0(ColBtüpe);
                            ColSblkLen = BM._ColPartitioning.GetSubblkLen(ColBtüpe);
                        } else {
                            BlockMsrMatrix.GetExternalSubblockIndices(BM._ColPartitioning, jPointer, out jBlockCol, out j0_Block, out jE_Block);
                            ColSblk_i0 = dummy_i0; // new int[] { 0 };
                            ColSblkLen = dummy_iE; // new int[] { jE_Block - j0_Block };
                            ColSblk_i0[0] = 0;
                            ColSblkLen[0] = (int)(jE_Block - j0_Block);
                        }
                        Debug.Assert(ColSblk_i0.Length == ColSblkLen.Length);
                        int NoOfColSblk = ColSblk_i0.Length;

                        for(int rowSblk = 0; rowSblk < NoOfRowSblk; rowSblk++) { // loop over row-sub-blocks...
                            long i0_Sblk = i0_Block + RowSblk_i0[rowSblk];
                            int NrwsSblk = RowSblkLen[rowSblk];
                            long iE_Sblk = i0_Sblk + NrwsSblk;
                            Debug.Assert(i0_Sblk >= i0_Block);
                            Debug.Assert(NrwsSblk >= 0);
                            Debug.Assert(iE_Sblk <= iE_Block);

                            if(iE_Sblk <= this.i0)
                                continue;
                            if(i0_Sblk >= this.iE)
                                continue;

                            long i0_clip = Math.Max(i0_Sblk, this.i0); // clip sub-block row range with block-row range
                            long iE_clip = Math.Min(iE_Sblk, this.iE); // 
                            int ICopy = (int)(iE_clip - i0_clip); //      number of rows to copy
                            if(ICopy <= 0)
                                continue;
                            int sblk_iOffset = (int)(i0_clip - i0_Sblk); // compute offset within sub-block ...
                            int this_iOffset = (int)(i0_clip - this.i0); // ...  and within block-row
                            Debug.Assert(sblk_iOffset >= 0);
                            Debug.Assert(this_iOffset >= 0);

                            for(int colSblk = 0; colSblk < NoOfColSblk; colSblk++) { // loop over column sub-blocks...
                                long j0_Sblk = j0_Block + ColSblk_i0[colSblk];
                                int NclsSblk = ColSblkLen[colSblk];
                                long jE_Sblk = j0_Sblk + NclsSblk;
                                Debug.Assert(j0_Sblk >= j0_Block);
                                Debug.Assert(NclsSblk >= 0);
                                Debug.Assert(jE_Sblk <= jE_Block);

                                if(jE_Sblk <= jPointer)
                                    continue;
                                if(j0_Sblk >= _jE)
                                    continue;

                                long j0_clip = Math.Max(j0_Sblk, jPointer);
                                long jE_clip = Math.Min(jE_Sblk, _jE);
                                int sblk_jOffset = (int)(j0_clip - j0_Sblk);
                                int this_jOffset = (int)(j0_clip - _j0);
                                Debug.Assert(sblk_jOffset >= 0);
                                Debug.Assert(this_jOffset >= 0);
                                int Jcopy = (int)(jE_clip - j0_clip);

                                double[] Storage;
                                int Offset, CI, CJ;
                                BM.GetSetAlloc2(true, iBlock, jBlockCol, rowSblk, colSblk, NrwsSblk, NclsSblk, NoOfRowSblk, NoOfColSblk, out Storage, out Offset, out CI, out CJ);

                                for(int i = 0; i < ICopy; i++) {
                                    for(int j = 0; j < Jcopy; j++) {
                                        int idxSrc = _offset + (i + this_iOffset) * _CI + (j + this_jOffset);
                                        int idxDst = Offset + (i + sblk_iOffset) * CI + (j + sblk_jOffset) * CJ;
                                        double V = srcMem[idxSrc];
                                        Storage[idxDst] += V;
                                    }
                                }
                            }
                        }

                        jPointer = jE_Block;
                    }
                }
            }

#if DEBUG_EXTENDED
            // facility to check that the block row is equal to a certain region of a matrix. 

            public long Check__i0_Block;
            public long Check__iE_Block;
            public long Check__jE_Block;

            void CheckEntries(BlockMsrMatrix BM) {
                Debug.Assert(Check__i0_Block >= this.i0);
                Debug.Assert(Check__iE_Block <= this.iE);
                
                int I = (int)(this.iE - this.i0);

                for (int bCol = 0; bCol < this.count; bCol++) { // loop over block-columns
                    int J = (int)(m_JES[bCol] - m_j0S[bCol]);

                    for (int i = 0; i < I; i++) { // loop over rows
                        long iM = i + i0; // row index into matrix

                        for (int j = 0; j < J; j++) { // loop over columns
                            long jM = j + m_j0S[bCol]; // column index into matrix
                            if (Check__jE_Block >= 0 && jM > Check__jE_Block)
                                continue;
                            
                            int idxB = m_Offset[bCol] + (i) * J + j; // index into the memory region of this object

                            if(iM >= Check__i0_Block && iM < Check__iE_Block)
                                Debug.Assert(BM[iM, jM] == m_mem[idxB]); 
                            else 
                                Debug.Assert(0.0 == m_mem[idxB]); 
                            unsafe
                            {
                                double* pMem = (double*)(this.m_pMem);
                                Debug.Assert(pMem[idxB] == m_mem[idxB]);
                            }
                        }
                    }
                }
            }


#endif
            public void Init(BlockMsrMatrix BM, long iBlock, bool doExternalCols) {

                long i0_Block = BM._RowPartitioning.GetBlockI0(iBlock);
                long iE_Block = BM._RowPartitioning.GetBlockLen(iBlock) + i0_Block;
                if(iE_Block <= i0 || iE <= i0_Block) {
                    count = 0;
                    return;
                }
#if DEBUG_EXTENDED
                //determine the clipping region
                Check__i0_Block = Math.Max(this.i0, i0_Block);
                Check__iE_Block = Math.Min(this.iE, iE_Block);
                Check__jE_Block = -1;
#endif

                long j0_Block = BM._ColPartitioning.FirstBlock;
                long jE_Block = BM._ColPartitioning.LocalNoOfBlocks + j0_Block;

                int BT = BM._RowPartitioning.GetBlockType(iBlock);
                int[] Sblk_i0 = BM._RowPartitioning.GetSubblk_i0(BT);
                int[] SblkLen = BM._RowPartitioning.GetSubblkLen(BT);
                long RowOffset = BM._RowPartitioning.GetBlockI0(iBlock);
                Debug.Assert(Sblk_i0.Length == SblkLen.Length);
                int NoOfSblk = Sblk_i0.Length;

                // clip the sub-blocks to the range [i0 .. iE[
                // -------------------------------------------
                int[] Sblk_i00 = new int[NoOfSblk];
                int[] SblkClen = new int[NoOfSblk];
                int[] SblkPntr = new int[NoOfSblk];
                int NoOfClipSblk = 0;

                int[] Dummy_i0 = new int[1];
                int[] Dummy_iE = new int[1];

                for(int SblkRow = 0; SblkRow < NoOfSblk; SblkRow++) { // loop over sub-blocks...
                    int rel_si0 = Sblk_i0[SblkRow]; // start of sub-block within block (relative)
                    int sLen = SblkLen[SblkRow]; //    sub-block length 
                    int rel_siE = rel_si0 + sLen; //   end of sub-block within block (relative)

                    long abs_si0 = RowOffset + rel_si0; // absolute start-index of sub-block
                    long abs_siE = RowOffset + rel_siE; // absolute end-index of sub-block

                    long Clip_abs_si0 = Math.Max(abs_si0, i0); // clip sub-block with (i0,iE)-range of block row
                    long Clip_abs_siE = Math.Min(abs_siE, iE);
                    int ClipLen = (int)(Clip_abs_siE - Clip_abs_si0); // length of clipped region
                    int Clipi00 = (int)(Clip_abs_si0 - RowOffset);  //   offset of clipped region within sub-block

                    if(ClipLen > 0) {
                        Debug.Assert(Clipi00 >= Sblk_i0[SblkRow]);
                        Debug.Assert(Clipi00 < Sblk_i0[SblkRow] + SblkLen[SblkRow]);
                        Debug.Assert(ClipLen <= SblkLen[SblkRow]);

                        Sblk_i00[NoOfClipSblk] = Clipi00;
                        SblkClen[NoOfClipSblk] = ClipLen;
                        SblkPntr[NoOfClipSblk] = SblkRow;
                        NoOfClipSblk++;
                    }
                }

                // loop over row columns
                // ---------------------
                if(NoOfClipSblk > 0) {
                    var Row = BM.m_BlockRows[iBlock - BM._RowPartitioning.FirstBlock];
                    if(Row != null) {
                        foreach(var kv in Row) {
                            long jBlock = kv.Key;
                            BlockEntry BE = kv.Value;
                            Debug.Assert(BE.jBlkCol == jBlock);

                            long ColBlockE;
                            int[] ColSblk_j0, ColSblkLen;
                            long ColOffset;
                            if(jBlock < j0_Block || jBlock >= jE_Block) {
                                if(doExternalCols) {
                                    long a, b;
                                    BlockMsrMatrix.GetExternalSubblockIndices(BM._ColPartitioning, jBlock, out a, out b);
                                    ColSblk_j0 = Dummy_i0;
                                    ColSblk_j0[0] = 0;
                                    ColSblkLen = Dummy_iE;
                                    ColSblkLen[0] = (int)(b - a);
                                    ColOffset = a;
                                    ColBlockE = b;
                                } else {
                                    continue;
                                }
                            } else {
                                int CBT = BM._ColPartitioning.GetBlockType(jBlock);
                                ColSblk_j0 = BM._ColPartitioning.GetSubblk_i0(CBT);
                                ColSblkLen = BM._ColPartitioning.GetSubblkLen(CBT);
                                ColOffset = BM._ColPartitioning.GetBlockI0(jBlock);
                                ColBlockE = ColOffset + BM._ColPartitioning.GetBlockLen(jBlock);
                            }
                            Debug.Assert(ColSblk_j0.Length == ColSblkLen.Length);
                            int NoOfColSblk = ColSblk_j0.Length;

                            int DestOffset, DestCI;
                            AllocBlock(ColOffset, ColBlockE, out DestOffset, out DestCI);
#if DEBUG_EXTENDED
                            
                            for (int i = 0; i < (this.iE - this.i0); i++) {
                                for(int j = 0; j < (ColBlockE - ColOffset); j++) {
                                    Debug.Assert(m_mem[i * DestCI + j + DestOffset] == 0.0);
                                    unsafe
                                    {
                                        Debug.Assert(((double*)(m_pMem))[i * DestCI + j + DestOffset] == 0.0);
                                    }
                                }
                            }
#endif

                            for(int SblkRow = 0; SblkRow < NoOfClipSblk; SblkRow++) {
                                for(int SblkCol = 0; SblkCol < NoOfColSblk; SblkCol++) {
                                    int _SblkRow = SblkPntr[SblkRow];

                                    int MembnkIdx = BE.MembnkIdx[_SblkRow, SblkCol];
                                    int InMembnk = BE.InMembnk[_SblkRow, SblkCol];
                                    Debug.Assert((InMembnk >= 0) == (MembnkIdx >= 0));
                                    if(MembnkIdx >= 0) {
                                        double[] Source;
                                        int Offset, CI, CJ;
                                        bool IsDense;
                                        BM.m_Membanks[MembnkIdx].GetFastBlockAccessInfo(out Source, out Offset, out CI, out CJ, out IsDense, InMembnk);

                                        int I = SblkClen[SblkRow];
                                        int di = Sblk_i00[SblkRow] - Sblk_i0[_SblkRow]; // offset within sub-block
                                        int J = ColSblkLen[SblkCol];
                                        int fj = ColSblk_j0[SblkCol];

                                        long ei = i0_Block + Sblk_i00[SblkRow]; // global start-index of the sub-block
                                        Debug.Assert(ei >= this.i0);
                                        Debug.Assert(ei < this.iE);
                                        Debug.Assert(ei + I <= this.iE);
                                        int fi = (int)(ei - this.i0); // relative start-index within this block row

                                        for(int i = 0; i < I; i++) {
#if DEBUG_EXTENDED
                                            long iM = i + Sblk_i00[SblkRow] + i0_Block; // row index into matrix
                                            Debug.Assert(i0_Block <= iM && iM < iE_Block);
#endif
                                            for(int j = 0; j < J; j++) {
#if DEBUG_EXTENDED
                                                long jM = j + ColSblk_j0[SblkCol] + ColOffset;
#endif
                                                int src_index = Offset + (i + di) * CI + j * CJ;
                                                double src_ij = Source[src_index];

#if DEBUG_EXTENDED
                                                Debug.Assert(BM[iM, jM] == src_ij);
#endif
                                                int dst_index = DestOffset + (i + fi) * DestCI + (j + fj);
                                                m_mem[dst_index] = src_ij;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }


#if DEBUG_EXTENDED
                Check__jE_Block = BM._ColPartitioning.TotalLength;
                CheckEntries(BM);
#endif

            }

          

            /// <summary>
            /// <paramref name="C"/> = <paramref name="_AscaleM_left"/>*<paramref name="A"/> + <paramref name="B"/>
            /// </summary>
            /// <param name="C">Output; the merge of <paramref name="A"/> and <paramref name="B"/></param>
            /// <param name="A">First Input.</param>
            /// <param name="B">Second Input.</param>
            /// <param name="_AscaleM_left"></param>
            /// <param name="AscaleM_I"></param>
            /// <param name="AscaleM_J"></param>
            unsafe static public void Merge(TempBlockRow C, TempBlockRow A, TempBlockRow B, double* _AscaleM_left, int AscaleM_I, int AscaleM_J) {
                Debug.Assert(!object.ReferenceEquals(C, A));
                Debug.Assert(!object.ReferenceEquals(C, B));
                Debug.Assert(!object.ReferenceEquals(A, B));
                Debug.Assert((A.m_mem == null && A.m_pMem == IntPtr.Zero) || (A.m_pMem == Marshal.UnsafeAddrOfPinnedArrayElement(A.m_mem, 0)));
                Debug.Assert((B.m_mem == null && B.m_pMem == IntPtr.Zero) || (B.m_pMem == Marshal.UnsafeAddrOfPinnedArrayElement(B.m_mem, 0)));
                Debug.Assert((C.m_mem == null && C.m_pMem == IntPtr.Zero) || (C.m_pMem == Marshal.UnsafeAddrOfPinnedArrayElement(C.m_mem, 0)));
#if DEBUG_EXTENDED
                Debug.Assert((_AscaleM_left != null) || ((A.iE - A.i0) == (B.iE - B.i0)));
                Debug.Assert((_AscaleM_left == null) || (AscaleM_I == (B.iE - B.i0) && AscaleM_J == (A.iE - A.i0)));
#endif
                C.Clear(B.i0, B.iE); // reset 'C'
                int NoRows = (int)(B.iE - B.i0);

                double* p_memA = (double*)A.m_pMem;
                double* p_memB = (double*)B.m_pMem;
                {
                    int counterA = 0;
                    int counterB = 0;

                    while(counterA < A.count || counterB < B.count) {
                        TempBlockRow NXT;
                        int counterNXT;
                        double* p_memNXT;
                        double* AscaleM_left;

                        // decide from which of the two rows we take the next block
                        // --------------------------------------------------------

                        if(counterA >= A.count) {
                            // finished A, but still entries in B
                            Debug.Assert(counterB < B.count);
                            NXT = B;
                            counterNXT = counterB;
                            p_memNXT = p_memB;
                            counterB++;
                            AscaleM_left = null;
                        } else if(counterB >= B.count) {
                            // finished B, but still entries in A
                            Debug.Assert(counterA < A.count);

                            NXT = A;
                            counterNXT = counterA;
                            p_memNXT = p_memA;
                            counterA++;
                            AscaleM_left = _AscaleM_left;
                        } else {
                            // pick the minimum of A and B
                            Debug.Assert(counterA < A.count && counterB < B.count);

                            if(A.m_j0S[counterA] <= B.m_j0S[counterB]) {
                                // pick A
                                NXT = A;
                                counterNXT = counterA;
                                p_memNXT = p_memA;
                                counterA++;
                                AscaleM_left = _AscaleM_left;
                            } else {
                                // pick B
                                NXT = B;
                                counterNXT = counterB;
                                p_memNXT = p_memB;
                                counterB++;
                                AscaleM_left = null;
                            }
                        }


                        // do the merge
                        // ------------

                        // access block:
                        int NXToffset, NXT_Ci;
                        NXT.AccessBlock(counterNXT, out NXToffset, out NXT_Ci);
                        int LNXT = NoRows * NXT_Ci;
                        double* p_memC, pS, pD;
                        int NoCols, C_Ci;

                        if(C.count == 0 || C.m_JES[C.count - 1] <= NXT.m_j0S[counterNXT]) {
                            // no overlap between the last entry in output C and next entry in A or B.

                            int Coffset;
                            C.AllocBlock(NXT.m_j0S[counterNXT], NXT.m_JES[counterNXT], out Coffset, out C_Ci);
                            Debug.Assert(C_Ci == NXT_Ci);

                            p_memC = (double*)C.m_pMem;
                            pS = p_memNXT + NXToffset;
                            pD = p_memC + Coffset;
                            NoCols = (int)(NXT.m_JES[counterNXT] - NXT.m_j0S[counterNXT]);
                            //Array.Copy(NXT.m_mem, NXToffset, C.m_mem, Coffset, LNXT);


                            Ker(pD, pS, C_Ci, NXT_Ci, 0.0, 1.0, AscaleM_left, NoRows, NoCols, AscaleM_I, AscaleM_J);


                        } else {
                            Debug.Assert(C.count >= 0);
                            Debug.Assert(NXT.m_j0S[counterNXT] < C.m_JES[C.count - 1]);

                            if(NXT.m_j0S[counterNXT] == C.m_j0S[counterNXT]
                                || NXT.m_JES[counterNXT] == C.m_JES[counterNXT]) {
                                // exact overlap of blocks

                                int Coffset;
                                C.AccessBlock(C.count - 1, out Coffset, out C_Ci);
                                Debug.Assert(C_Ci == NXT_Ci);

                                p_memC = (double*)C.m_pMem;
                                pS = p_memNXT + NXToffset;
                                pD = p_memC + Coffset;
                                NoCols = (int)(NXT.m_JES[counterNXT] - NXT.m_j0S[counterNXT]);
                                Ker(pD, pS, C_Ci, NXT_Ci, 1.0, 1.0, AscaleM_left, NoRows, NoCols, AscaleM_I, AscaleM_J);

                            } else {
                                // partial overlap

                                int Coffset;
                                if(C.m_JES[C.count - 1] < NXT.m_JES[counterNXT]) {
                                    // need to increase size of last block in C

                                    C.EnlargeLastBlock(NXT.m_JES[counterNXT], out Coffset, out C_Ci);
                                } else {
                                    C.AccessBlock(C.count - 1, out Coffset, out C_Ci);
                                }
                                Debug.Assert(NXT.m_j0S[counterNXT] >= C.m_j0S[C.count - 1]);
                                Debug.Assert(NXT.m_JES[counterNXT] <= C.m_JES[C.count - 1]);

                                p_memC = (double*)C.m_pMem;
                                NoCols = (int)(NXT.m_JES[counterNXT] - NXT.m_j0S[counterNXT]);
                                pS = p_memNXT + NXToffset;
                                pD = p_memC + Coffset + (NXT.m_j0S[counterNXT] - C.m_j0S[C.count - 1]);
                                Ker(pD, pS, C_Ci, NXT_Ci, 1.0, 1.0, AscaleM_left, NoRows, NoCols, AscaleM_I, AscaleM_J);

                            }
                        }
                    }

                    Debug.Assert((C.count > 0) == (A.count > 0 || B.count > 0));
                    Debug.Assert((A.count <= 0 && B.count <= 0) || (C.m_j0S != null));

                }
            }

            unsafe static void Ker(double* pC, double* pA, int C_Ci, int A_Ci, double Cscale, double Ascale, double* AscaleM_left, int NoRows, int NoCols, int AscaleM_I, int AscaleM_J) {
                if(AscaleM_left == null) {
                    if(Cscale != 0.0) {
                        for(int i = 0; i < NoRows; i++) {
                            for(int j = 0; j < NoCols; j++) {
                                pC[j] = pC[j] * Cscale + pA[j] * Ascale;
                            }
                            pA += A_Ci;
                            pC += C_Ci;
                        }
                    } else {
                        // by using a separate branch for Cscale == 0, we ensure that 
                        // e.g. NANs get cleared/overwritten (otherwise, NAN*0.0 = NAN).
                        for(int i = 0; i < NoRows; i++) {
                            for(int j = 0; j < NoCols; j++) {
                                pC[j] = pA[j] * Ascale;
                            }
                            pA += A_Ci;
                            pC += C_Ci;
                        }
                    }
                } else {
                    // C  <--  C*Cscale  +  AscaleM_left*pA*Ascale
                    // C                : NoRows * NoCols,     row skip: C_Ci
                    // AscaleM_left     : NoRows * AscaleM_J,  row skip: AscaleM_J
                    // pA               : AscaleM_J * NoCols,  row skip: A_Ci

                    /*
                    double[] pCcheck = new double[C_Ci * NoRows];

                    for (int i = 0; i < NoRows; i++) {
                        for (int j = 0; j < NoCols; j++) {
                            double a = 0; // accumulator

                            for (int k = 0; k < AscaleM_J; k++) {
                                a += AscaleM_left[AscaleM_J * i + k] * pA[A_Ci * k + j];
                            }

                            // --------------
                            if (Cscale != 0) {
                                //pC[i * C_Ci + j] *= Cscale;
                                //pC[i * C_Ci + j] += a * Ascale;
                                pCcheck[i * C_Ci + j] = pC[i * C_Ci + j] * Cscale + a * Ascale;

                            } else {
                                //pC[i * C_Ci + j] = a * Ascale;

                                pCcheck[i * C_Ci + j] = a * Ascale;
                            }
                        }
                    }
                    */

                    BLAS.dgemm('n', 'n', NoCols, NoRows, AscaleM_J, Ascale, pA, A_Ci, AscaleM_left, AscaleM_J, Cscale, pC, C_Ci);
                    /*
                    double eps = Math.Sqrt(BLAS.MachineEps);
                    for (int i = 0; i < NoRows; i++) {
                        for (int j = 0; j < NoCols; j++) {
                            if (!(Math.Abs(pCcheck[i * C_Ci + j] - pC[i * C_Ci + j]) <= Math.Max(eps, 1e-7 * Math.Abs(pCcheck[i * C_Ci + j]))))
                                throw new Exception("fuck blas.");
                        }
                    }
                    */
                }
            }

            public static TempBlockRow Init(BlockMsrMatrix BM, long i0, long iE, bool doExternalCols, TempBlockRow A1, TempBlockRow A2, TempBlockRow A3) {
                Debug.Assert(BM._RowPartitioning.i0 <= i0);
                Debug.Assert(iE <= BM._RowPartitioning.iE);
                Debug.Assert(iE > i0);
                unsafe {
                    TempBlockRow[] A = new TempBlockRow[] { A1, A2, A3 };
                    int Acnt = 0;
                    TempBlockRow OldAccu = null;

                    int initCount = 0, mergeCount = 0;
                    long _i = i0;
                    while(_i < iE) {
                        long iBlock = BM._RowPartitioning.GetBlockIndex(_i);

                        TempBlockRow Anext = A[Acnt];
                        Acnt++;
                        Acnt = Acnt % 3;

                        Anext.Clear(i0, iE);
                        Anext.Init(BM, iBlock, doExternalCols);
                        initCount++;

                        if(OldAccu == null) {
                            OldAccu = Anext;
                        } else {
                            TempBlockRow NewAccu = A[Acnt];
                            Acnt++;
                            Acnt = Acnt % 3;



                            NewAccu.Clear(i0, iE);

                            TempBlockRow.Merge(NewAccu, OldAccu, Anext, null, -1, -1);
                            mergeCount++;

                            OldAccu = NewAccu;
                            NewAccu = null;
                        }

#if DEBUG_EXTENDED
                        OldAccu.Check__jE_Block = BM._ColPartitioning.TotalLength;
                        OldAccu.Check__i0_Block = i0;
                        OldAccu.Check__iE_Block = Math.Min(OldAccu.iE, BM._RowPartitioning.GetBlockI0(iBlock) + BM._RowPartitioning.GetBlockLen(iBlock));
                        OldAccu.CheckEntries(BM);
#endif

                        int BlockLen = BM._RowPartitioning.GetBlockLen(iBlock);
                        _i = BM._RowPartitioning.GetBlockI0(iBlock) + BlockLen;
                    }

#if DEBUG_EXTENDED
                    Debug.Assert(OldAccu.iE == iE);
                    Debug.Assert(OldAccu.i0 == i0);
                    OldAccu.Check__jE_Block = BM._ColPartitioning.TotalLength;
                    OldAccu.Check__i0_Block = i0;
                    OldAccu.Check__iE_Block = iE;
                    OldAccu.CheckEntries(BM);
#endif

                    return OldAccu;
                }
            }


        }

        /// <summary>
        /// Sparse Matrix-Matrix multiplication.
        /// </summary>
        /// <param name="C">On exit, hopefully equal to <paramref name="A"/>*<paramref name="B"/>;</param>
        /// <param name="A">Left operand.</param>
        /// <param name="B">Right operand.</param>
        public static void Multiply(BlockMsrMatrix C, BlockMsrMatrix A, BlockMsrMatrix B) {
            if(multiply == null)
                multiply = new Stopwatch();
            if(multiply_core == null)
                multiply_core = new Stopwatch();
#if DEBUG_EXTENDED
            // Deactivated by Florian, 10feb2018
            // testcode was in here for more than a year, no complaints, should be fine
            //A.VerifyDataStructure("Multiply_A");
            //B.VerifyDataStructure("Multiply_B");
            //C.VerifyDataStructure("Multiply_C");

            //MsrMatrix _A = A.ToMsrMatrix();
            //MsrMatrix _B = B.ToMsrMatrix();
            //MsrMatrix _C_b4 = C.ToMsrMatrix();

            __Multiply(C, A, B);

            C.VerifyDataStructure("Multiply_Cout");

            //MsrMatrix _C_af = C.ToMsrMatrix();
            //MsrMatrix AB = MsrMatrix.Multiply(_A, _B);
            //double ABnorm = AB.InfNorm();
            //_C_af.Acc(-1.0, _C_b4);
            //_C_af.Acc(-1.0, AB);

            //double ErrAbs = _C_af.InfNorm();
            //double ErrRel = ABnorm > double.Epsilon ? ErrAbs / ABnorm : ErrAbs;
            ////Console.WriteLine("SpMM check: " + ErrRel);
            //if (ErrRel > 1.0e-8 || double.IsNaN(ErrRel) || double.IsInfinity(ErrRel))
            //    throw new ArithmeticException("Error in multiply");
#else
            multiply.Start();
            __Multiply(C, A, B);
            multiply.Stop();
#endif
        }

        public static Stopwatch multiply;
        public static Stopwatch multiply_core;


        /// <summary>
        /// Sparse Matrix-Matrix multiplication.
        /// </summary>
        /// <param name="C">On exit, hopefully equal to <paramref name="A"/>*<paramref name="B"/>;</param>
        /// <param name="A">Left operand.</param>
        /// <param name="B">Right operand.</param>
        private static void __Multiply(BlockMsrMatrix C, BlockMsrMatrix A, BlockMsrMatrix B) {
            using(new FuncTrace()) {
                if(C.RowPartitioning.LocalLength != A.RowPartitioning.LocalLength)
                    throw new ArgumentException("Number of rows mismatch between C and A.");
                if(C.RowPartitioning.i0 != A.RowPartitioning.i0)
                    throw new ArgumentException("Partitioning mismatch between C and A.");
                if(A.ColPartition.LocalLength != B.RowPartitioning.LocalLength)
                    throw new ArgumentException("Number of columns mismatch between A and B.");
                if(B.ColPartition.LocalLength != C._ColPartitioning.LocalLength)
                    throw new ArgumentException("Number of columns mismatch between C and B.");
                MPI_Comm comm = C.MPI_Comm;
                if(A.MPI_Comm != comm)
                    throw new ArgumentException("MPI communicator mismatch between C and A.");
                if(B.MPI_Comm != comm)
                    throw new ArgumentException("MPI communicator mismatch between C and B.");
                MPICollectiveWatchDog.Watch(comm);
                int MpiRank = C.RowPartitioning.MpiRank;
                int MpiSize = C.RowPartitioning.MpiSize;

                // set up communication
                // ---------------------

                A.UpdateCommPattern(comm);

                long A_iBlk0 = A._RowPartitioning.FirstBlock;
                int A_NoBlk = A._RowPartitioning.LocalNoOfBlocks;
                long A_IBlkE = A_iBlk0 + A_NoBlk;

                long A_coliBlk0 = A._ColPartitioning.FirstBlock;
                int A_colNoBlk = A._ColPartitioning.LocalNoOfBlocks;
                long A_colIBlkE = A_coliBlk0 + A_colNoBlk;



                // collect list of blocks which have to be sent to other processors
                // ----------------------------------------------------------------

                // key: target proc. / value: stuff to copy
                Dictionary<int, List<Multiply_Helper_1>> MpiSendCollection = new Dictionary<int, List<Multiply_Helper_1>>();
                {
                    foreach(int iProc in A.ReceiveLists.Keys)
                        MpiSendCollection.Add(iProc, new List<Multiply_Helper_1>());

                    for(long iBlkRow = A_iBlk0; iBlkRow < A_IBlkE; iBlkRow++) {
                        if(A.m_ExternalBlock[(int)(iBlkRow - A_iBlk0)]) {
                            var row = A.m_BlockRows[iBlkRow - A_iBlk0];

                            int BT = A._RowPartitioning.GetBlockType(iBlkRow);
                            int[] Sbkl_i0 = A._RowPartitioning.GetSubblk_i0(BT);
                            int[] SbklLen = A._RowPartitioning.GetSubblkLen(BT);
                            Debug.Assert(Sbkl_i0.Length == SbklLen.Length);
                            int NoOfSblk = Sbkl_i0.Length;
                            long RowOffset = A._RowPartitioning.GetBlockI0(iBlkRow);

                            foreach(var kv in row) {
                                long jBlkCol = kv.Key;
                                BlockEntry BE = kv.Value;
                                Debug.Assert(BE.jBlkCol == jBlkCol, "1");
                                Debug.Assert(BE.MembnkIdx.GetLength(0) == NoOfSblk, "2");
                                Debug.Assert(BE.InMembnk.GetLength(0) == NoOfSblk, "3");

                                if(jBlkCol < A_coliBlk0 || jBlkCol >= A_colIBlkE) {
                                    // external block
                                    Debug.Assert(BE.MembnkIdx.GetLength(1) == 1, "4");
                                    Debug.Assert(BE.InMembnk.GetLength(1) == 1, "5");

                                    long j0, jE;
                                    BlockMsrMatrix.GetExternalSubblockIndices(A._ColPartitioning, jBlkCol, out j0, out jE);
                                    int proc = A._ColPartitioning.FindProcess(j0);
                                    long[,] ReceiveList = A.ReceiveLists[proc];
#if DEBUG_EXTENDED
                                    Debug.Assert(A.ReceiveLists.ContainsKey(proc), "6");
#endif

                                    List<Multiply_Helper_1> copyColection = MpiSendCollection[proc];

                                    for(int iSblk = 0; iSblk < NoOfSblk; iSblk++) {
                                        if(BE.MembnkIdx[iSblk, 0] >= 0) {
                                            // this sub-block must be sent to processor 'proc'
                                            Multiply_Helper_1 H;
                                            H.i0 = RowOffset + Sbkl_i0[iSblk];
                                            H.j0 = j0;
                                            H.MembankIdx = BE.MembnkIdx[iSblk, 0];
                                            H.InMembank = BE.InMembnk[iSblk, 0];
                                            copyColection.Add(H);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // sent external blocks to other processors
                // ----------------------------------------

                int[] RecvRanks = A.SendLists.Keys.ToArray(); // MPI ranks from which this processor receives data
                int[] SendRanks = MpiSendCollection.Keys.ToArray(); //

                MiniTransceiver miniTx1 = new MiniTransceiver(comm, SendRanks, RecvRanks);

                {
                    // calculate size of send buffer
                    int GetSendBufferSize(int Rank) {
                        List<Multiply_Helper_1> copyColection = MpiSendCollection[Rank];
                        int Sz = sizeof(int);
                        foreach(var H in copyColection) {
                            Membank Mbnk = A.m_Membanks[H.MembankIdx];
                            int I = Mbnk.Mem.GetLength(1);
                            int J = Mbnk.Mem.GetLength(2);
                            Sz += sizeof(long) * 2 + sizeof(int) * 2 + sizeof(double) * I * J;
                        }
                        return Sz;
                    }

                    // assemble send buffer
                    void FillSendBuffer(int Rank, IntPtr Buffer) {
                        unsafe {
                            List<Multiply_Helper_1> copyColection = MpiSendCollection[Rank];
                            void* pBuffer = (void*)Buffer;
                            int* pI = (int*)pBuffer;
                            *pI = copyColection.Count;
                            pI++;
                            pBuffer = pI;
                            Debug.Assert((((byte*)pBuffer) - ((byte*)Buffer)) == sizeof(int), "8");
                            foreach(var H in copyColection) {
                                Membank Mbnk = A.m_Membanks[H.MembankIdx];
                                Debug.Assert(Mbnk.Occupied.Length == Mbnk.Mem.GetLength(0), "9");
                                int I = Mbnk.Mem.GetLength(1);
                                int J = Mbnk.Mem.GetLength(2);

                                long* pH1 = (long*)pBuffer;
                                pH1[0] = H.i0;
                                pH1[1] = H.j0;
                                pH1 += 2;
                                int* pH2 = (int*)pH1; pH1 = null;
                                pH2[0] = I;
                                pH2[1] = J;

                                double* pB = (double*)(pH2 + 2); pH2 = null;

                                double[] buf;
                                int offset, Ci, Cj;
                                bool isDense;
                                Mbnk.GetFastBlockAccessInfo(out buf, out offset, out Ci, out Cj, out isDense, H.InMembank);

                                for(int i = 0; i < I; i++) {
                                    for(int j = 0; j < J; j++) {
                                        Debug.Assert(i * Ci + j * Cj >= 0);
                                        Debug.Assert(i * Ci + j * Cj < I * J);

                                        pB[i * J + j] = buf[offset + i * Ci + j * Cj];
                                    }
                                }

                                double* _pB = pB;
                                pB += I * J;
                                pBuffer = pB;
                                Debug.Assert((((byte*)pBuffer) - ((byte*)_pB)) == I * J * sizeof(double), "10");
                            }
#if DEBUG_EXTENDED
                            Debug.Assert((((byte*)pBuffer) - ((byte*)Buffer)) == GetSendBufferSize(Rank), "10");
#endif
                        }
                    }

                    miniTx1.TransceiveStart(GetSendBufferSize, FillSendBuffer);
                }
                // perform local sparse matrix-matrix multiplication
                // -------------------------------------------------

                long C_iBlk0 = C._RowPartitioning.FirstBlock;
                int C_NoBlk = C._RowPartitioning.LocalNoOfBlocks;
                long C_IBlkE = C_iBlk0 + C_NoBlk;

                unsafe {
                    IntPtr tempBuf = IntPtr.Zero;
                    int tempBuf_Length = -1;

                    using(TempBlockRow rc1 = new TempBlockRow(), rc2 = new TempBlockRow(),
                        rA1 = new TempBlockRow(), rA2 = new TempBlockRow(), rA3 = new TempBlockRow(), bR = new TempBlockRow()) {

                        for(long iBlkRow = C_iBlk0; iBlkRow < C_IBlkE; iBlkRow++) { // loop over the rows of 'C'

                            // extract row from 'C'
                            long i0 = C._RowPartitioning.GetBlockI0(iBlkRow);
                            long iE = i0 + C._RowPartitioning.GetBlockLen(iBlkRow);
                            if(iE == i0)
                                // empty block row
                                continue;
                            TempBlockRow Caccu = rc1, Cnext = rc2;
                            Caccu.Clear(i0, iE);
                            Caccu.Init(C, iBlkRow, true);

                            // extract row from 'A'
                            TempBlockRow Arow = TempBlockRow.Init(A, i0, iE, false, rA1, rA2, rA3);
                            int NoRows = (int)(Arow.iE - Arow.i0);
                            Debug.Assert(NoRows == C._RowPartitioning.GetBlockLen(iBlkRow));

                            if(Arow.count > 0) {
                                //int jBlk = 0;
                                //while (jBlk < Arow.count) {

                                int NoCols = int.MinValue;
                                int _BlkCount = 0;
                                bool _inBlock = false;
                                long _j0 = long.MinValue, _jE = long.MinValue;
                                long jPointer = 0;
                                bool inBlock = false;
                                long jBlock = long.MinValue, j0 = long.MinValue, jE = long.MinValue;
                                double* pAij = null;
                                double AijSum = 0.0;
                                while(_BlkCount < Arow.count) {
                                    if(!_inBlock) {
                                        _j0 = Arow.m_j0S[_BlkCount];
                                        _jE = Arow.m_JES[_BlkCount];
                                        _inBlock = true;
                                        jPointer = _j0;

                                        Debug.Assert(_j0 >= A.ColPartition.i0, "12");
                                        Debug.Assert(_jE <= A.ColPartition.iE, "13");
                                        Debug.Assert(_j0 < _jE, "14");
                                    }
                                    if(!inBlock) {
                                        jBlock = B._RowPartitioning.GetBlockIndex(jPointer);
                                        j0 = B._RowPartitioning.GetBlockI0(jBlock);
                                        jE = B._RowPartitioning.GetBlockLen(jBlock) + j0;
                                        NoCols = (int)(jE - j0);
                                        inBlock = true;
                                    }

                                    long jNext = Math.Min(jE, _jE);

                                    int RelOffsetSource = (int)(jPointer - _j0);
                                    int RelOffsetDest = (int)(jPointer - j0);
                                    int Length = (int)(jNext - jPointer);
                                    Debug.Assert(Length >= 0, "15");
                                    Debug.Assert(RelOffsetSource >= 0, "16");
                                    Debug.Assert(RelOffsetDest >= 0, "17");
                                    //Console.WriteLine("{0} {1}--{2} off: {3}  ===>  {4} {5}--{6} off: {7}", _BlkCount, j, jNext, RelOffsetSource,
                                    //    jBlock, j, jNext, RelOffsetDest);
                                    bool ExMatch = false;
                                    if(Length > 0) {
                                        ExMatch = RelOffsetDest == 0 && RelOffsetSource == 0 && Length == (jE - j0);
                                        Debug.Assert(ExMatch == (pAij == null), "18");

                                        if(ExMatch) {
                                            int PointerOffset, ANoOfCols;
                                            Arow.AccessBlock(_BlkCount, out PointerOffset, out ANoOfCols);
                                            pAij = ((double*)(Arow.m_pMem)) + PointerOffset;
                                            Debug.Assert(ANoOfCols == NoCols, "19");

                                            double* _pAij = pAij;
                                            for(int l = NoRows * NoCols; l > 0; l--) {
                                                AijSum += Math.Abs(*_pAij);
                                                _pAij++;
                                                if(AijSum != 0.0)
                                                    break;
                                            }

                                        } else {
                                            if(pAij == null) {
                                                int ReqSize = NoRows * NoCols;
                                                if(ReqSize < tempBuf_Length) {
                                                    if(tempBuf != IntPtr.Zero)
                                                        Marshal.FreeHGlobal(tempBuf);
                                                    tempBuf_Length = ReqSize * 2;
                                                    tempBuf = Marshal.AllocHGlobal(tempBuf_Length * sizeof(double));
                                                }
                                                pAij = (double*)tempBuf;
                                                double* pD = pAij;
                                                for(int ii = NoCols * NoRows; ii > 0; ii--) {
                                                    *pD = 0;
                                                    pD++;
                                                }
                                            }

                                            int PointerOffset, C_Ai;
                                            Arow.AccessBlock(_BlkCount, out PointerOffset, out C_Ai);
                                            double* pAijSrc = (double*)(Arow.m_pMem);
                                            pAijSrc += PointerOffset;

                                            for(int i = 0; i < NoRows; i++) {
                                                for(int j = 0; j < Length; j++) {
                                                    double vl = pAijSrc[i * C_Ai + j + RelOffsetSource];
                                                    pAij[i * NoCols + j + RelOffsetDest] = vl;
                                                    AijSum += Math.Abs(vl);
                                                }
                                            }
                                        }
                                    }

                                    if(jNext == jE) {
                                        // block complete
                                        //Console.WriteLine("destination block complete");
                                        if(AijSum != 0) {
                                            bR.Clear(j0, jE);
                                            bR.Init(B, jBlock, true);
                                            Cnext.Clear(Caccu.i0, Caccu.iE);
                                            multiply_core.Start();
                                            TempBlockRow.Merge(Cnext, bR, Caccu, pAij, NoRows, NoCols);
                                            multiply_core.Stop();
                                            var _Caccu = Caccu;
                                            Caccu = Cnext;
                                            Cnext = _Caccu;
                                        }

                                        inBlock = false;
                                        pAij = null;
                                    }

                                    if(jNext == _jE) {
                                        //Console.WriteLine("source block complete");
                                        _inBlock = false;
                                        _BlkCount++;
                                    }


                                    jPointer = jNext;
                                }
                                //}
                            }


                            // save row
                            C.ClearBlockRow(iBlkRow);
                            Caccu.Save(C, iBlkRow);
                        }
                    }

                    Debug.Assert((tempBuf_Length >= 0) == (tempBuf != IntPtr.Zero), "20");
                    if(tempBuf_Length >= 0) {
                        Marshal.FreeHGlobal(tempBuf);
                    }
                } // end of local multiplication...

                // receive external blocks & perform multiplication
                // ------------------------------------------------

                // data structure for external rows:
                // 1st key: owner processor rank
                // 2nd key: row index
                // content: the block row
                Dictionary<int, Dictionary<long, TempBlockRow>> ExternalRows = new Dictionary<int, Dictionary<long, TempBlockRow>>();

                TempBlockRow CrowNext = new TempBlockRow();
                using(TempBlockRow rt1 = new TempBlockRow(), rt2 = new TempBlockRow(), rt3 = new TempBlockRow()) {

                    void OnReceive(int OriginRank, IntPtr Buffer, int _RecvSize) {
                        // create data structure for the respective rank
                        unsafe {
                            Dictionary<long, TempBlockRow> ExternalRows_rank;
                            if(!ExternalRows.TryGetValue(OriginRank, out ExternalRows_rank)) {
                                ExternalRows_rank = new Dictionary<long, TempBlockRow>();
                                ExternalRows.Add(OriginRank, ExternalRows_rank);
                            }

                            // start computing the external stuff
                            void* pBuffer = (void*)Buffer;
                            int* pI = (int*)pBuffer;
                            int NoOfBuffers = *pI;
                            pI++;
                            pBuffer = pI;

                            for(int iBuf = 0; iBuf < NoOfBuffers; iBuf++) {
                                long* pH1 = (long*)pBuffer;
                                long i0 = pH1[0];
                                long j0 = pH1[1];
                                int* pH2 = (int*)(pH1 + 2); pH1 = null;
                                int I = pH2[0];
                                int J = pH2[1];
                                Debug.Assert(B._RowPartitioning.IsInLocalRange(j0), "2nd wave: first index not in local range");
                                Debug.Assert(J <= 0 || B._RowPartitioning.IsInLocalRange(j0 + J - 1), "2nd wave: last index not in local range");
                                Debug.Assert(A._RowPartitioning.FindProcess(i0) == OriginRank, "2nd wave: receive index mismatch (1)");
                                Debug.Assert(I <= 0 || A._RowPartitioning.FindProcess(i0 + I - 1) == OriginRank, "2nd wave: receive index mismatch (2)");
                                double* pBlock = (double*)(pH2 + 2); pH2 = null;

                                double* _pAij = pBlock;
                                double AijSum = 0;
                                for(int l = I * J; l > 0; l--) {
                                    AijSum += Math.Abs(*_pAij);
                                    _pAij++;
                                    if(AijSum != 0.0)
                                        break;
                                }

                                if(AijSum != 0.0) {
                                    // extract block-row [j0 .. j0+J[ from B,
                                    // multiply with pBlock from left,
                                    // send back (somehow) later

                                    TempBlockRow Btmp = TempBlockRow.Init(B, j0, j0 + J, true, rt1, rt2, rt3);

                                    // access row
                                    TempBlockRow Crow;
                                    if(!ExternalRows_rank.TryGetValue(i0, out Crow)) {
                                        Crow = new TempBlockRow();
                                        Crow.Clear(i0, i0 + I);
                                        ExternalRows_rank.Add(i0, Crow);
                                    } else {
                                        Debug.Assert(Crow.i0 == i0, "22");
                                        Debug.Assert(Crow.iE == i0 + I, "23");
                                    }

                                    // multiply & merge
                                    TempBlockRow.Merge(CrowNext, Btmp, Crow, pBlock, I, J);

                                    // swap
                                    var _CrowNext = CrowNext;
                                    CrowNext = Crow;
                                    ExternalRows_rank[i0] = _CrowNext;
                                }

                                // move to next block
                                pBlock += I * J;
                                pBuffer = pBlock;

                            }
#if DEBUG_EXTENDED
                            Debug.Assert((((byte*)pBuffer) - ((byte*)Buffer)) == _RecvSize, "mismatch between pointer difference and receive size");
#endif
                        }
                    }


                    miniTx1.Finish(OnReceive);
                }
                CrowNext.Dispose();

                // send the block rows back to their respective owners
                // ---------------------------------------------------

                // now, we have to return data, i.e. the role of send and receive ranks flips
                {
                    var tmp = SendRanks;
                    SendRanks = RecvRanks; // processes which receive data from us
                    RecvRanks = tmp;       // processes from which we receive data
                }
                MiniTransceiver miniTx2 = new MiniTransceiver(comm, SendRanks, RecvRanks);

                {
                    int GetSendBufferSize(int DestRank) {
                        var Value = ExternalRows[DestRank];

                        int Sz = sizeof(int);
                        foreach(var _kv in Value) {
                            var r = _kv.Value;
                            Sz += r.GetSize();
                        }

                        return Sz;
                    }

                    // serialize block rows, send data
                    void FillSendBuffer(int DestRank, IntPtr _Buffer) {
                        unsafe {
                            var RowCollection = ExternalRows[DestRank];

                            int* pBuffer = (int*)_Buffer;
                            *pBuffer = RowCollection.Count;
                            pBuffer++;
                            IntPtr Buffer = (IntPtr)pBuffer;
                            foreach(var _kv in RowCollection) {
                                var r = _kv.Value;
                                Debug.Assert(C._RowPartitioning.FindProcess(r.i0) == DestRank, "24");
                                Debug.Assert(C._RowPartitioning.FindProcess(Math.Max(r.i0, r.iE - 1)) == DestRank, "25");

                                Buffer = r.Serialize(Buffer);
                            }
#if DEBUG_EXTENDED
                            Debug.Assert((byte*)Buffer - (byte*)_Buffer == GetSendBufferSize(DestRank));
#endif
                        }
                    }

                    miniTx2.TransceiveStart(GetSendBufferSize, FillSendBuffer);

                }


                // wait for the block rows to come in and accumulate them
                // ------------------------------------------------------

                using(TempBlockRow rt1 = new TempBlockRow()) {

                    void OnReceive(int OriginRank, IntPtr _Buffer, int _RecvSize) {
                        unsafe {
                            // de-serialize data & accumulate
                            int* pBuffer = (int*)_Buffer;
                            int NoOfRows = *pBuffer;
                            pBuffer++;
                            IntPtr Buffer = (IntPtr)(pBuffer);
                            for(int j = 0; j < NoOfRows; j++) {
                                Buffer = rt1.Deserialize(Buffer);
                                Debug.Assert(C._RowPartitioning.IsInLocalRange(rt1.i0), "31");
                                Debug.Assert(C._RowPartitioning.IsInLocalRange(Math.Max(rt1.i0, rt1.iE - 1)), "32");


                                // accumulation:
                                long jPointer = rt1.i0;
                                while(jPointer < rt1.iE) {
                                    long iBlock = C._RowPartitioning.GetBlockIndex(jPointer);
                                    rt1.Save(C, iBlock);
                                    int BlockLen = C._RowPartitioning.GetBlockLen(iBlock);
                                    jPointer += BlockLen;
                                }
                            }
#if DEBUG_EXTENDED
                            Debug.Assert(((byte*)Buffer - (byte*)_Buffer) == _RecvSize, "33");
#endif
                        }
                    }

                    miniTx2.Finish(OnReceive);
                }

                // dispose external rows
                // ---------------------
                {
                    foreach(var k in ExternalRows.Values)
                        foreach(var t in k.Values)
                            t.Dispose();
                }
            }
        }

        class MiniTransceiver {

            public MiniTransceiver(MPI_Comm __comm, int[] __RanksToSendTo, int[] __RanksToReceiveFrom) {
                comm = __comm;
                MPICollectiveWatchDog.Watch(comm);
                RanksToSendTo = __RanksToSendTo;
                RanksToReceiveFrom = __RanksToReceiveFrom;

#if DEBUG_EXTENDED
                csMPI.Raw.Comm_Size(comm, out int MPIsize);
                csMPI.Raw.Comm_Rank(comm, out int MPIrank);

                int[] isTarget = new int[MPIsize];
                foreach(int iDest in RanksToSendTo) {
                    isTarget[iDest] = 1;
                }

                int[] Combined = isTarget.MPIAllGather(comm);
                var Sources = new List<int>();
                for(int r = 0; r < MPIsize; r++) {
                    if(Combined[MPIrank + r*MPIsize] != 0) {
                        Sources.Add(r);
                    } 
                }

                Debug.Assert(__RanksToReceiveFrom.SetEquals(Sources), "MPI-global mismatch in sources and destinations.");
#endif
            }

            MPI_Comm comm;
            int[] RanksToSendTo;
            int[] RanksToReceiveFrom;

#if DEBUG_EXTENDED
            /// <summary>
            /// 1st index: origin process rank 'o'
            /// 2nd index: destination process rank 'r'
            /// content: number of bytes send from rank 'o' to 'r'
            /// </summary>
            int[,] SizesCheck;
#endif

            internal void Finish(Action<int, IntPtr, int> ReceiveAction) {
                MPICollectiveWatchDog.Watch(comm);
                unsafe {
                    int NoOfSend = secondWave_SendBuffers.Length;
                    int NoOfRecv = secondWave_ReceiveBuffer.Length;
                    Debug.Assert(secondWave.Length == NoOfSend + NoOfRecv, "21");

#if DEBUG_EXTENDED
                    bool[] checker = new bool[NoOfRecv + NoOfSend];
#endif
                    for(int i = 0; i < NoOfRecv + NoOfSend; i++) {
                        MPI_Status Stat;
                        int index;
                        csMPI.Raw.Waitany(NoOfSend + NoOfRecv, secondWave, out index, out Stat);

#if DEBUG_EXTENDED
                        Debug.Assert(checker[index] == false, $"2nd wave: request {index} already handled.");
                        checker[index] = true;
#endif
                        if(index >= NoOfSend) {
                            // received buffer
                            Debug.Assert(Stat.MPI_TAG == 32567, "2nd wave, receiving: tag mismatch in MPI status.");
                            Debug.Assert(Stat.MPI_SOURCE == RanksToReceiveFrom[index - NoOfSend], "2nd wave, receiving: source rank mismatch in MPI status.");
#if DEBUG_EXTENDED
                            //the count filed seeems to be wrongly mapped for OpenMPI
                            //Debug.Assert(Stat.count == _RecvSize[index - NoOfSend], $"2nd wave, receiving: receive count mismatch -- Status count {Stat.count}, expected {_RecvSize[index - NoOfSend]}.");
#endif
                            IntPtr Buffer = secondWave_ReceiveBuffer[index - NoOfSend];

                            // 
                            ReceiveAction(RanksToReceiveFrom[index - NoOfSend], Buffer, _RecvSize[index - NoOfSend]);

                            Marshal.FreeHGlobal(secondWave_ReceiveBuffer[index - NoOfSend]);
                            secondWave_ReceiveBuffer[index - NoOfSend] = IntPtr.Zero;
                        } else {
                            // send operation complete: free send buffer
                            Marshal.FreeHGlobal(secondWave_SendBuffers[index]);
                            secondWave_SendBuffers[index] = IntPtr.Zero;
                        }
                    }
                }
            }

            internal void TransceiveStart(Func<int, int> GetSentBufferSize, Action<int, IntPtr> FillSendBuffer) {
                MPICollectiveWatchDog.Watch(comm);
#if DEBUG_EXTENDED
                csMPI.Raw.Comm_Size(comm, out int MPIsize);
                csMPI.Raw.Comm_Rank(comm, out int MPIrank);
                int[] SendSizes = new int[MPIsize];
#endif
                unsafe {
                    int NoOfSend = RanksToSendTo.Length;
                    int NoOfRecv = RanksToReceiveFrom.Length;


                    MPI_Request[] firstWave = new MPI_Request[NoOfSend + NoOfRecv];
                    secondWave = new MPI_Request[NoOfSend + NoOfRecv];
                    secondWave_SendBuffers = new IntPtr[NoOfSend];


                    // initiate receiving of buffer size
                    int* RecvSize = stackalloc int[NoOfRecv];
                    for(int i = 0; i < RanksToReceiveFrom.Length; i++) {
                        RecvSize[i] = int.MinValue;
                        csMPI.Raw.Irecv((IntPtr)(RecvSize + i), 1, csMPI.Raw._DATATYPE.INT, RanksToReceiveFrom[i], 567, comm, out firstWave[NoOfSend + i]);
                    }

                    // serialize blocks & send them
                    int* SendSize = stackalloc int[NoOfSend];
                    for(int i = 0; i < RanksToSendTo.Length; i++) {
                        int DestinationProc = RanksToSendTo[i];

                        // calculate size of send buffer
                        int Sz = GetSentBufferSize(DestinationProc);
                        SendSize[i] = Sz;
#if DEBUG_EXTENDED
                        SendSizes[DestinationProc] = Sz;
#endif
                        // send size of send buffer
                        csMPI.Raw.Issend((IntPtr)(SendSize + i), 1, csMPI.Raw._DATATYPE.INT, DestinationProc, 567, comm, out firstWave[i]);

                        // allocate send buffer
                        IntPtr Buffer = Marshal.AllocHGlobal(Sz);
                        secondWave_SendBuffers[i] = Buffer;

                        // assemble send buffer
                        FillSendBuffer(DestinationProc, Buffer);

                        // send send buffer
                        csMPI.Raw.Issend(Buffer, Sz, csMPI.Raw._DATATYPE.BYTE, DestinationProc, 32567, comm, out secondWave[i]);
                    }

                    // initiate receiving of data
                    secondWave_ReceiveBuffer = new IntPtr[NoOfRecv];
#if DEBUG_EXTENDED
                    bool[] checker = new bool[NoOfRecv + NoOfSend];
                    SizesCheck = SendSizes.MPIAllGather(comm).Reshape(false, MPIsize, MPIsize);
#endif
                    for(int i = 0; i < NoOfRecv + NoOfSend; i++) {
                        MPI_Status Stat;
                        int index;
                        csMPI.Raw.Waitany(NoOfSend + NoOfRecv, firstWave, out index, out Stat);

#if DEBUG_EXTENDED
                        Debug.Assert(checker[index] == false, $"1st wave: request {index} already handled.");
                        checker[index] = true;
#endif
                        if(index >= NoOfSend) {
                            // received buffer size
#if DEBUG_EXTENDED                            
                            //Debug.Assert(Stat.count == 1 * sizeof(int), $"1st wave, receiving: receive count mismatch -- Status count {Stat.count}, expected {sizeof(int)}.");
                            Debug.Assert(RecvSize[index - NoOfSend] >= 0, "11");
                            Debug.Assert(Stat.MPI_TAG == 567, $"1st wave, receiving: tag mismatch.");
                            Debug.Assert(Stat.MPI_SOURCE == RanksToReceiveFrom[index - NoOfSend], $"1st wave, receiving: source mismatch.");

                            Debug.Assert(RecvSize[index - NoOfSend] == SizesCheck[RanksToReceiveFrom[index - NoOfSend], MPIrank], $"Buffer size check failed: received size {RecvSize[index - NoOfSend]}, check says {SizesCheck[RanksToReceiveFrom[index - NoOfSend], MPIrank]}");
#endif
                            IntPtr Buffer = Marshal.AllocHGlobal(RecvSize[index - NoOfSend]);
                            secondWave_ReceiveBuffer[index - NoOfSend] = Buffer;
                            csMPI.Raw.Irecv(Buffer, RecvSize[index - NoOfSend], csMPI.Raw._DATATYPE.BYTE, RanksToReceiveFrom[index - NoOfSend], 32567, comm, out secondWave[index]);
                        } else {
                            // nop
                        }
                    }
                    _RecvSize = new int[NoOfRecv];
                    for(int i44 = 0; i44 < NoOfRecv; i44++) {
                        _RecvSize[i44] = RecvSize[i44];
                    }
                }
#if DEBUG_EXTENDED
                
                for(int r = 0; r < MPIsize; r++) {
                    bool receiveFrom_r_check = SizesCheck[r, MPIrank] != 0;
                    bool receiveFrom_r = RanksToReceiveFrom.Contains(r);
                    Debug.Assert(receiveFrom_r == receiveFrom_r_check, "Receive/Send pattern failed 2nd check.");

                    int idx = Array.IndexOf(RanksToReceiveFrom, r);
                    int recvSize = idx >= 0 ? _RecvSize[idx] : 0;
                    Debug.Assert(recvSize == SizesCheck[r, MPIrank], $"Buffer size check failed (2): received size {recvSize}, check says {SizesCheck[r, MPIrank]}");
                }
#endif
            }
            // for checking: number of bytes we receive from the respective rank, index correlates with 'RecvRanks'
            int[] _RecvSize;

            // requests for the send *and* receive ops: [0 .. MpiSendCollection.Count[      : MPI_Issend -- operations
            //                                          [MpiSendCollection.Count .. Length[ : MPI_Irecv  -- operations
            MPI_Request[] secondWave;

            IntPtr[] secondWave_SendBuffers;
            IntPtr[] secondWave_ReceiveBuffer;
        }

        /// <summary>
        /// Sparse Matrix-Matrix multiplication.
        /// </summary>
        /// <param name="A">Left operand.</param>
        /// <param name="B">Right operand.</param>
        /// <returns>Hopefully equal to <paramref name="A"/>*<paramref name="B"/>;</returns>
        public static BlockMsrMatrix Multiply(BlockMsrMatrix A, BlockMsrMatrix B) {
            BlockMsrMatrix C = new BlockMsrMatrix(A._RowPartitioning, B._ColPartitioning);
            Multiply(C, A, B);
            return C;
        }

        #endregion

        /// <summary>
        /// Returns a matrix with inverted blocks.
        /// </summary>
        /// <param name="OnlyDiagonal">
        /// If true, only diagonal blocks resp. sub-blocks are inverted.
        /// </param>
        /// <param name="Subblocks">
        /// Whether the blocks or sub-blocks should be inverted.
        /// - false: inversion of (full) blocks, as defined by <see cref="IBlockPartitioning.GetBlockI0"/> and <see cref="IBlockPartitioning.GetBlockLen"/>.
        /// - true: sub-block inversion, sub-blocks are defined by the <see cref="IBlockPartitioning"/>, 
        ///    see especially <see cref="IBlockPartitioning.GetBlockType"/>, <see cref="IBlockPartitioning.GetSubblk_i0"/> and <see cref="IBlockPartitioning.GetSubblkLen"/>
        /// </param>
        /// <param name="ignoreEmptyBlocks"></param>
        /// <param name="SymmetricalInversion">
        /// If true, symmetrical blocks are assumed, which are inverted by <see cref="IMatrixExtensions.InvertSymmetrical{T}(T)"/>.
        /// </param>
        /// <returns></returns>
        public BlockMsrMatrix InvertBlocks(bool OnlyDiagonal = true, bool Subblocks = false, bool ignoreEmptyBlocks = false, bool SymmetricalInversion = false) {
            using(new FuncTrace()) {
                // Matrix for the result
                BlockMsrMatrix res = new BlockMsrMatrix(this._RowPartitioning, this._ColPartitioning);

                // FullMatrix to calculate inverse
                MultidimensionalArray tmp = null;
                //for (int BlockNo = 0; BlockNo < NoOfBlocks; BlockNo++) {

                Debug.Assert(m_BlockRows.Length == _RowPartitioning.LocalNoOfBlocks);
                for(int iBlockLoc = 0; iBlockLoc < m_BlockRows.Length; iBlockLoc++) {
                    var BlockRow = this.m_BlockRows[iBlockLoc];
                    if(BlockRow != null) {
                        long iBlockGlb = iBlockLoc + this._RowPartitioning.FirstBlock;
                        long i0 = this._RowPartitioning.GetBlockI0(iBlockGlb);

                        foreach(var KV in BlockRow) {
                            long jBlockGlb = KV.Key;
                            BlockEntry Block = KV.Value;
                            Debug.Assert(jBlockGlb == Block.jBlkCol);

                            if(OnlyDiagonal && (iBlockGlb != jBlockGlb)) {
                                // skip of-diagonal entries
                                continue;
                            }

                            long j0 = this._ColPartitioning.GetBlockI0(jBlockGlb);

                            if(Subblocks) {
                                // ++++++++++++++++++
                                // work on sub-blocks
                                // ++++++++++++++++++

                                int NoOfSblk_rows = Block.InMembnk.GetLength(0);
                                int NoOfSblk_cols = Block.InMembnk.GetLength(1);
                                int RowBlockType = _RowPartitioning.GetBlockType(iBlockGlb);
                                int ColBlockType = _ColPartitioning.GetBlockType(jBlockGlb);
                                int[] RowSblk_i0 = _RowPartitioning.GetSubblk_i0(RowBlockType);
                                int[] RowSblkLen = _RowPartitioning.GetSubblkLen(RowBlockType);
                                Debug.Assert(RowSblk_i0.Length == NoOfSblk_rows);
                                Debug.Assert(RowSblkLen.Length == NoOfSblk_rows);

                                int[] ColSblk_i0 = _ColPartitioning.GetSubblk_i0(ColBlockType);
                                int[] ColSblkLen = _ColPartitioning.GetSubblkLen(ColBlockType);
                                Debug.Assert(ColSblk_i0.Length == NoOfSblk_cols);
                                Debug.Assert(ColSblkLen.Length == NoOfSblk_cols);

                                for(int sblkRow = 0; sblkRow < NoOfSblk_rows; sblkRow++) {
                                    for(int sblkCol = 0; sblkCol < NoOfSblk_cols; sblkCol++) {
                                        if(OnlyDiagonal && (sblkCol != sblkRow)) {
                                            // skip of-diagonal entries
                                            continue;
                                        }
                                        if(Block.InMembnk[sblkCol, sblkRow] < 0 || Block.MembnkIdx[sblkCol, sblkRow] < 0)
                                            continue;

                                        int M = RowSblkLen[sblkRow];
                                        int N = ColSblkLen[sblkCol];

                                        if(tmp == null || tmp.GetLength(0) != M || tmp.GetLength(1) != N)
                                            tmp = MultidimensionalArray.Create(M, N);

                                        // Copy values from current block to FullMatrix
                                        this.ReadBlock(i0 + RowSblk_i0[sblkRow], j0 + ColSblk_i0[sblkCol], tmp);

                                        if(ignoreEmptyBlocks && tmp.AbsSum() < double.Epsilon) {
                                            // No inversion of empty blocks requested
                                            continue;
                                        }

                                        // Calculate inverse of FullMatrix
                                        if(SymmetricalInversion)
                                            tmp.InvertSymmetrical();
                                        else
                                            tmp.InvertInPlace();

                                        // Write inverse of FullMatrix to current block
                                        res.AccBlock(i0 + RowSblk_i0[sblkRow], j0 + ColSblk_i0[sblkCol], 1.0, tmp);
                                    }
                                }


                            } else {
                                // ++++++++++++++
                                // work on blocks
                                // ++++++++++++++

                                int M = this._RowPartitioning.GetBlockLen(iBlockGlb);
                                int N = this._ColPartitioning.GetBlockLen(jBlockGlb);


                                if(tmp == null || tmp.GetLength(0) != M || tmp.GetLength(1) != N)
                                    tmp = MultidimensionalArray.Create(M, N);

                                // Copy values from current block to FullMatrix
                                this.ReadBlock(i0, j0, tmp);

                                if(ignoreEmptyBlocks && tmp.AbsSum() < 1e-15) {
                                    // No inversion of empty blocks requested
                                    continue;
                                }

                                // Calculate inverse of FullMatrix
                                if(SymmetricalInversion)
                                    tmp.InvertSymmetrical();
                                else
                                    tmp.InvertInPlace();

                                // Write inverse of FullMatrix to current block
                                res.AccBlock(i0, j0, 1.0, tmp);
                            }
                        }
                    }
                }

                return res;
            }
        }
        /*
        public IDictionary<int,MultidimensionalArray> GetBlockRow(int iBlockRow) {
            var R = new SortedDictionary<int, MultidimensionalArray>();

            if (!this._RowPartitioning.IsLocalBlock(iBlockRow))
                throw new ArgumentException("Row Block index out of range.");

            var IntRow = this.m_BlockRows[iBlockRow - this._RowPartitioning.FirstBlock];
            if(IntRow != null) {

                int I = _RowPartitioning.GetBlockLen(iBlockRow);
                int i0 = _RowPartitioning.GetBlockI0(iBlockRow);

                foreach(int jBlockCol in IntRow.Keys) {
                    int J = _ColPartitioning.GetBlockLen(jBlockCol);
                    int j0 = _RowPartitioning.GetBlockI0(iBlockRow);

                    MultidimensionalArray Blk = MultidimensionalArray.Create(I, J);
                    this.ReadBlock(i0, j0, Blk);
                    R.Add(i0, Blk);
                }
            }
            
            return R;
        }
        */




        
    }
}
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers
{
    public class SubBlockSelector : SubBlockSelectorBase {
        public SubBlockSelector(MultigridMapping map) : base(map) { }
    }


    /// <summary>
    /// Enables creation of submatrices of some or single cells and masks vectors accordingly. Useful for linear solvers in particular (e.g. Schwarz, levelpmg, etc.)
    /// </summary>
    public class BlockMask {

        /// <summary>
        /// Auxiliary class, provides Local Block masks
        /// </summary>
        private class BlockMaskLoc : BlockMaskBase {
            public BlockMaskLoc(SubBlockSelector sbs) : base(sbs, csMPI.Raw._COMM.SELF) {
                foreach (var idx in this.m_GlobalMask) {
                    Debug.Assert(idx >= m_map.i0);
                    Debug.Assert(idx < m_map.iE);
                }
            }

            protected override int m_NoOfCells {
                get {
                    return m_map.LocalNoOfBlocks;
                }
            }

            protected override int m_CellOffset {
                get {
                    return 0;
                }
            }

            protected override int m_LocalLength {
                get {
                    return m_map.LocalLength;
                }
            }

        }

        /// <summary>
        /// Auxiliary class, provides masking of external cells: overlap that is located on other MPI processes
        /// </summary>
        private class BlockMaskExt : BlockMaskBase {

            public BlockMaskExt(SubBlockSelector SBS) : base(SBS, SBS.GetMapping.MPI_Comm) {
                foreach (int idx in this.m_GlobalMask) {
                    Debug.Assert(idx < m_map.i0 || idx >= m_map.iE);
                }

                int LL = m_map.LocalLength;
                int jMax = m_map.AggGrid.iLogicalCells.Count - 1;
                int LE = m_map.LocalUniqueIndex(0, jMax, 0) + m_map.GetLength(jMax);

                foreach (int idx in this.m_LocalMask) {
                    Debug.Assert(idx >= LL);
                    Debug.Assert(idx < LE);
                }
                m_ExtCellLen = m_map.GetLocalLength_Ext();
            }

            private static int m_ExtCellLen;

            protected override int m_NoOfCells {
                get {
                    return m_map.AggGrid.iLogicalCells.NoOfExternalCells;
                }
            }

            protected override int m_CellOffset {
                get {
                    return m_map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                }
            }

            /// <summary>
            /// These are the DOF of external ghost cells available on this proc
            /// </summary>
            protected override int m_LocalLength {
                get {
                    return m_ExtCellLen;
                }
            }

            //public BlockMsrMatrix GetExternalRows(BlockMsrMatrix target) {

            //    var ExternalRows_BlockI0 = base.GetAllSubMatrixCellOffsets();
            //    var ExternalRows_BlockN = base.GetAllSubMatrixCellLength();
            //    var ExternalRowsIndices = base.m_GlobalMask;
            //    var Mop = target;
            //    BlockPartitioning PermRow = new BlockPartitioning(ExternalRowsIndices.Count, ExternalRows_BlockI0, ExternalRows_BlockN, Mop.MPI_Comm, i0isLocal: true);

            //    // Remark: we use a permutation matrix for MPI-exchange of rows
            //    //  set   Perm[l,m] = I;
            //    //  then  ExternalRowsTemp[l,-]  =  Mop[m,-]

            //    BlockMsrMatrix Perm = new BlockMsrMatrix(PermRow, Mop._RowPartitioning);
            //    for (int iRow = 0; iRow < ExternalRowsIndices.Count; iRow++) {
            //        Debug.Assert(Mop._RowPartitioning.IsInLocalRange(ExternalRowsIndices[iRow]) == false);
            //        Perm[iRow + PermRow.i0, ExternalRowsIndices[iRow]] = 1;
            //    }

            //    return BlockMsrMatrix.Multiply(Perm, Mop);
            //}

        }

        /// <summary>
        /// generates a masking of a <see cref="BlockMsrMatrix"/> according to <see cref="SubBlockSelector"/> <paramref name="sbs"/>
        /// and provides operations exclusively executed on mask
        /// (e.g. sub matrix generation, sub vector extraction, etc.).
        /// The smallest unit are DGdegree blocks.
        /// The masking operates on MPI local blocks per default,
        /// which can be exceeded to external blocks by providing external rows: <paramref name="ExtRows"/>.
        /// Ghost cells (covert by <see cref="BoSSS.Foundation.Grid.ILogicalCellData.NoOfExternalCells"/>) can be acquired by <see cref="GetAllExternalRows"/>.
        /// </summary>
        /// <param name="sbs">selection instruction</param>
        /// <param name="ExtRows">external rows collected from other MPI-processes on this proc</param>
        public BlockMask(SubBlockSelector sbs, BlockMsrMatrix ExtRows=null) {
            m_map = sbs.GetMapping;
            m_ExtRows = ExtRows;
            m_includeExternalCells = (ExtRows != null) && m_map.MpiSize > 1;
            BMLoc = new BlockMaskLoc(sbs);
            
            if (m_includeExternalCells) {
                BMExt = new BlockMaskExt(sbs);
                SetThisShitUp(new BlockMaskBase[] { BMLoc, BMExt });
            } else {
                SetThisShitUp(new BlockMaskBase[] { BMLoc });
            }
#if Debug
            CheckIndices();
#endif
        }

        private void CheckIndices() {
            for(int i=0;i < BMLoc.m_GlobalMask.Count; i++) {
                int GlobIdx = BMLoc.m_GlobalMask.ToArray()[i];
                Debug.Assert(m_map.IsInLocalRange(GlobIdx));
            }
            for(int i = 0; i < BMExt.m_GlobalMask.Count; i++) {
                int GlobIdx = BMExt.m_GlobalMask.ToArray()[i];
                Debug.Assert(!m_map.IsInLocalRange(GlobIdx));
            }
        }

        private void SetThisShitUp(BlockMaskBase[] masks) {
            List<int> tmpOffsetList = new List<int>();
            List<int> tmpLengthList = new List<int>();
            List<extNi0[][][]> tmpNi0 = new List<extNi0[][][]>();
            foreach (var mask in masks) {
                var tmp = mask.GetAllSubMatrixCellOffsets().ToArray();
                //If there are any external cells do something different:
                if (mask.GetType() == typeof(BlockMaskExt)) {
                    int offset = BMLoc.LocalDOF;
                    for (int i = 0; i < tmp.Length; i++)
                        tmp[i] += offset;
                }
                tmpOffsetList.AddRange(tmp);
                tmpLengthList.AddRange(mask.GetAllSubMatrixCellLength());
                tmpNi0.AddRange(mask.m_StructuredNi0.ToList());

            }
            Debug.Assert(tmpOffsetList != null);
            Debug.Assert(tmpLengthList != null);
            Debug.Assert(tmpNi0 != null);
            Debug.Assert(tmpOffsetList.GroupBy(x => x).Any(g => g.Count() == 1));
            Debug.Assert(tmpNi0.GroupBy(x => x).Any(g => g.Count() == 1));

            SubMatrixOffsets = tmpOffsetList;
            SubMatrixLen = tmpLengthList;
            StructuredNi0 = tmpNi0.ToArray();
        }

        BlockMaskLoc BMLoc;
        BlockMaskExt BMExt;
        List<int> SubMatrixOffsets;
        List<int> SubMatrixLen;
        bool m_includeExternalCells;
        extNi0[][][][] StructuredNi0;
        MultigridMapping m_map;
        BlockMsrMatrix m_ExtRows;

        //public MultigridMapping GenerateMappingfromMask() {

        //}


        /// <summary>
        /// If you want nothing special. Take this one. If you want only diagonal block matrix choose one of the other methods
        /// </summary>
        /// <returns></returns>
        public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix source) {
            if (source == null)
                throw new ArgumentNullException();
            if (source.NoOfRows < BMLoc.LocalDOF)
                throw new ArgumentException();

            BlockMsrMatrix target;
            if (m_includeExternalCells) {

                BlockPartitioning targetBlocking = new BlockPartitioning(BMLoc.LocalDOF + BMExt.LocalDOF, SubMatrixOffsets, SubMatrixLen, csMPI.Raw._COMM.SELF);

                //make an extended block dummy to fit local and external blocks
                target = new BlockMsrMatrix(targetBlocking, targetBlocking);

                //get the external rows via MPI exchange
                var ExtRowsTmp = m_ExtRows;

                int offset = BMLoc.m_GlobalMask.Count;
                var extBlockRows = BMExt.m_GlobalMask.Count.ForLoop(i => i + offset);
                var extBlockCols = BMExt.m_GlobalMask.Count.ForLoop(i => i + offset);

                //ExtRowsTmp lives at the MPI-Communicator of the target, thus the global index is related to a new partitioning and has nothing to do with the partitioning of the multigrid operator ...
                
                var GlobalIdxExtRows = BMExt.m_GlobalMask.Count.ForLoop(i => BMExt.m_LocalMask[i] - m_map.LocalLength);
                for (int iGlob=0; iGlob< GlobalIdxExtRows.Count(); iGlob++) {
                    Debug.Assert(GlobalIdxExtRows[iGlob] < ExtRowsTmp._RowPartitioning.LocalLength);
                    GlobalIdxExtRows[iGlob] += ExtRowsTmp._RowPartitioning.i0;
                    Debug.Assert(ExtRowsTmp._RowPartitioning.IsInLocalRange(GlobalIdxExtRows[iGlob]));
                }
                
                //add local Block ...
                source.WriteSubMatrixTo(target, BMLoc.m_GlobalMask, default(int[]), BMLoc.m_GlobalMask, default(int[]));

                //add columns related to external rows ...
                source.AccSubMatrixTo(1.0, target, BMLoc.m_GlobalMask, default(int[]), new int[0], default(int[]), BMExt.m_GlobalMask, extBlockCols);

                //add external rows ...
                ExtRowsTmp.AccSubMatrixTo(1.0, target, GlobalIdxExtRows, extBlockRows, BMLoc.m_GlobalMask, default(int[]), BMExt.m_GlobalMask, extBlockCols);
            } else {
                BlockPartitioning localBlocking = new BlockPartitioning(BMLoc.LocalDOF, SubMatrixOffsets, SubMatrixLen, csMPI.Raw._COMM.SELF, i0isLocal: true);

                target = new BlockMsrMatrix(localBlocking);

                source.AccSubMatrixTo(1.0, target, BMLoc.m_GlobalMask, default(int[]), BMLoc.m_GlobalMask, default(int[]));
            }
            Debug.Assert(target != null);
            return target;
        }


        /// <summary>
        /// Gets number of blocks covert by mask
        /// </summary>
        public int GetNoOfMaskedCells {
            get {
                if (m_includeExternalCells) {
                    return BMLoc.m_StructuredNi0.Length + BMExt.m_StructuredNi0.Length;
                } else {
                    return BMLoc.m_StructuredNi0.Length;
                }
            }
        }

        /// <summary>
        /// Gets number of rows covert by mask
        /// </summary>
        public int GetNoOfMaskedRows {
            get {
                if (m_includeExternalCells) {
                    return BMLoc.LocalDOF + BMExt.LocalDOF;
                } else {
                    return BMLoc.LocalDOF;
                }
            }
        }

#region Testing
        // ==========
        // Stuff dedicated to testing ...
        // ==========

        protected List<int> GlobalIList_Internal {
            get {
                return BMLoc.m_GlobalMask;
            }
        }

        protected List<int> GlobalIList_External {
            get {
                return BMExt.m_GlobalMask;
            }
        }

        protected int[] GetGlobalIdxOfCell(int iCell) {
            if (iCell < m_map.LocalNoOfBlocks)
                return BMLoc.GetCellwiseGlobalidx(iCell);
            else if (iCell < m_map.AggGrid.iLogicalCells.NoOfExternalCells)
                return BMExt.GetCellwiseGlobalidx(iCell);
            else
                throw new NotSupportedException("Selected Cells are outside valid bounds. This is not supported right now");
        }

        protected int[] GetLocalIdxOfCell(int iCell) {
            if (iCell < m_map.LocalNoOfBlocks)
                return BMLoc.GetCellwiseLocalidx(iCell);
            else if (iCell < m_map.AggGrid.iLogicalCells.NoOfExternalCells)
                return BMExt.GetCellwiseLocalidx(iCell);
            else
                throw new NotSupportedException("Selected Cells are outside of valid bounds. This is not supported right now");
        }

#endregion
        /// <summary>
        /// Get Array of Cellblocks
        /// </summary>
        /// <param name="target"></param>
        /// <param name="ignoreCellCoupling"></param>
        /// <param name="ignoreVarCoupling"></param>
        /// <param name="ignoreSpecCoupling"></param>
        /// <returns></returns>
        public MultidimensionalArray[] GetSubBlocks(BlockMsrMatrix target, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {

            List<MultidimensionalArray> Sblocks = new List<MultidimensionalArray>();
            Sblocks.AddRange(GetMaskBlocks(target, BMLoc, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling));

            if (m_includeExternalCells&& BMExt!=null) {
                //var ExtMatrix = this.GetExternalSubBlockMatrixOnly(target);
                //int RowOffset = m_map.iE;
                Sblocks.AddRange(GetMaskBlocks(m_ExtRows, BMExt, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling));
            }

            Debug.Assert(Sblocks!=null);
            return Sblocks.ToArray();
        }

        //private BlockMsrMatrix GetExternalSubBlockMatrixOnly(BlockMsrMatrix target) {
        //    BlockPartitioning Blocking = new BlockPartitioning(BMExt.LocalDOF, BMExt.GetAllSubMatrixCellOffsets(), BMExt.GetAllSubMatrixCellLength(), m_map.MPI_Comm,i0isLocal:true);

        //    //make an extended block dummy to fit local and external blocks
        //    var Block = new BlockMsrMatrix(Blocking, Blocking);

        //    //get the external rows via MPI exchange
        //    BlockMsrMatrix ExtRowsTmp = m_ExtRows;

        //    int offset = BMLoc.m_GlobalMask.Count;
        //    var extBlockRows = BMExt.m_GlobalMask.Count.ForLoop(i => i );

        //    //ExtRowsTmp lives at the MPI-Communicator of the target, thus the global index is related to a new partitioning and has nothing to do with the partitioning of the multigrid operator ...

        //    var GlobalIdxExtRows = BMExt.m_SubBlockMask;
        //    for (int iGlob = 0; iGlob < GlobalIdxExtRows.Count(); iGlob++) {
        //        Debug.Assert(GlobalIdxExtRows[iGlob] < ExtRowsTmp._RowPartitioning.LocalLength);
        //        GlobalIdxExtRows[iGlob] += ExtRowsTmp._RowPartitioning.i0;
        //        Debug.Assert(ExtRowsTmp._RowPartitioning.IsInLocalRange(GlobalIdxExtRows[iGlob]));
        //    }

        //    ExtRowsTmp.WriteSubMatrixTo(Block, GlobalIdxExtRows, default(int[]), GlobalIdxExtRows, default(int[]));
        //    //ExtRowsTmp.AccSubMatrixTo(1.0, extBlock, GlobalIdxExtRows, extBlockRows, BMLoc.m_GlobalMask, default(int[]), BMExt.m_GlobalMask, extBlockCols);
        //    //public void AccSubMatrixTo<V1, V2, V3, V4, V5, V6>(
        //    //double alpha, IMutableMatrixEx Target,
        //    //V1 RowIndicesSource, V2 RowIndicesTarget,
        //    //V3 ColumnIndicesSource, V4 ColIndicesTarget,
        //    //V5 ExternalColumnIndicesSource, V6 ExternalColIndicesTarget)
        //    return Block;
        //}


        private MultidimensionalArray[] GetMaskBlocks(BlockMsrMatrix source, BlockMaskBase mask, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {
            //BlockMsrMatrix ExtRowsTmp = null;
            //if (BMExt != null) {
            //    ExtRowsTmp = BMExt.GetExternalRows(target);
            //}
            bool ListenToGlobalId = true;
            // pay attention to offset of local stuff ...
            // target is concidered to 
            if (mask.GetType() == typeof(BlockMaskExt)) {
                ListenToGlobalId = false;
            }

            int NoOfCells = mask.m_StructuredNi0.Length;
            int size = ignoreCellCoupling ? NoOfCells : NoOfCells * NoOfCells;

            MultidimensionalArray[] Sblocks = new MultidimensionalArray[size];

            int auxIdx = 0;
            for (int iLoc = 0; iLoc < mask.m_StructuredNi0.Length; iLoc++) {
                for (int jLoc = 0; jLoc < mask.m_StructuredNi0.Length; jLoc++) {
                    if (ignoreCellCoupling && jLoc != iLoc) {
                        continue;
                    }
                    int CellBlockLen = mask.GetCellwiseLength(jLoc);
                    Sblocks[auxIdx] = MultidimensionalArray.Create(CellBlockLen, CellBlockLen);
                    for (int iVar = 0; iVar < mask.m_StructuredNi0[iLoc].Length; iVar++) {
                        for (int jVar = 0; jVar < mask.m_StructuredNi0[jLoc].Length; jVar++) {
                            if (ignoreVarCoupling && jVar != iVar) {
                                continue;
                            }
                            for (int iSpc = 0; iSpc < mask.m_StructuredNi0[iLoc][iVar].Length; iSpc++) {
                                for (int jSpc = 0; jSpc < mask.m_StructuredNi0[jLoc][jVar].Length; jSpc++) {
                                    if (ignoreSpecCoupling && jSpc != iSpc) {
                                        continue;
                                    }
                                    for (int iMode = 0; iMode < mask.m_StructuredNi0[iLoc][iVar][iSpc].Length; iMode++) {
                                        for (int jMode = 0; jMode < mask.m_StructuredNi0[jLoc][jVar][jSpc].Length; jMode++) {
                                            extNi0 RowNi0 = mask.m_StructuredNi0[iLoc][iVar][iSpc][iMode];
                                            extNi0 ColNi0 = mask.m_StructuredNi0[jLoc][jVar][jSpc][jMode];
                                            int Targeti0 = ListenToGlobalId ? RowNi0.Gi0 : RowNi0.Li0 + source._RowPartitioning.i0 - m_map.LocalLength;
                                            int Targetj0 = ColNi0.Gi0;
                                            //int Targeti0 = RowNi0.Gi0;
                                            //int Targetj0 = ColNi0.Gi0;
                                            int Subi0 = mask.GetRelativeSubBlockOffset(iLoc, iVar, iSpc, iMode);
                                            int Subj0 = mask.GetRelativeSubBlockOffset(jLoc, jVar, jSpc, jMode);
                                            int Subie = Subi0 + RowNi0.N - 1;
                                            int Subje = Subj0 + ColNi0.N - 1;

                                            var tmp = Sblocks[auxIdx].ExtractSubArrayShallow(new int[] { Subi0, Subj0 }, new int[] { Subie, Subje });

                                            Debug.Assert((m_map.IsInLocalRange(Targeti0) && m_map.IsInLocalRange(Targetj0) && mask.GetType() == typeof(BlockMaskLoc)) || mask.GetType() == typeof(BlockMaskExt));


                                                try {
                                                    source.ReadBlock(Targeti0, Targetj0, tmp
                                                );
                                                } catch (Exception e) {
                                                    Console.WriteLine("row: " + Targeti0);
                                                    Console.WriteLine("col: " + Targetj0);
                                                    throw new Exception(e.Message);
                                                }
                                            //ilPSP.Environment.StdoutOnlyOnRank0 = false;
                                            //if (!ListenToGlobalId) Console.WriteLine("proc_{0}: row {1}, col {2}", m_map.MpiRank, Targeti0, Targetj0);

                                        }
                                    }
                                }
                            }
                        }
                    }
                    auxIdx++;
                }
            }
            return Sblocks;
        }

        /// <summary>
        /// Coupling can be ignored. If you want full output take the basic GetSubBlockMatrix() method, which will be faster.
        /// </summary>
        /// <param name="target"></param>
        /// <param name="ignoreCellCoupling"></param>
        /// <param name="ignoreVarCoupling"></param>
        /// <param name="ignoreSpecCoupling"></param>
        /// <returns></returns>
        public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix target, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {


            BlockMsrMatrix submatrix = null;
            if (m_includeExternalCells) {
                BlockPartitioning extBlocking = new BlockPartitioning(BMLoc.LocalDOF + BMExt.LocalDOF, SubMatrixOffsets, SubMatrixLen, csMPI.Raw._COMM.SELF);
                submatrix = new BlockMsrMatrix(extBlocking);
                var ExtRowsTmp = m_ExtRows;
                GetMaskMatrix(submatrix, target, BMLoc, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling);
                GetMaskMatrix(submatrix, ExtRowsTmp, BMExt, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling);
            } else {
                int Loclength = BMLoc.LocalDOF;
                var tmpN = BMLoc.GetAllSubMatrixCellLength();
                var tmpi0 = BMLoc.GetAllSubMatrixCellOffsets();
                BlockPartitioning localBlocking = new BlockPartitioning(Loclength, tmpi0.ToArray(), tmpN.ToArray(), csMPI.Raw._COMM.SELF, i0isLocal: true);
                submatrix = new BlockMsrMatrix(localBlocking);
                GetMaskMatrix(submatrix, target, BMLoc, ignoreCellCoupling,  ignoreVarCoupling,  ignoreSpecCoupling);
            }
            Debug.Assert(submatrix!=null);
            return submatrix;
        }


      
        private void GetMaskMatrix(BlockMsrMatrix target, BlockMsrMatrix source, BlockMaskBase mask, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {

            var SubMSR = target;

            //int SubRowIdx = 0;
            //int SubColIdx = 0;

            bool ListenToGlobalId = true;
            // pay attention to offset of local stuff ...
            // target is concidered to 
            if (mask.GetType() == typeof(BlockMaskExt)) {
                ListenToGlobalId = false;
            }

            int auxIdx = 0;
            for (int iLoc = 0; iLoc < mask.m_StructuredNi0.Length; iLoc++) {
                for (int jLoc = 0; jLoc < mask.m_StructuredNi0.Length; jLoc++) {
                    if (ignoreCellCoupling && jLoc != iLoc) {
                        continue;
                    }
                    for (int iVar = 0; iVar < mask.m_StructuredNi0[iLoc].Length; iVar++) {
                        for (int jVar = 0; jVar < mask.m_StructuredNi0[jLoc].Length; jVar++) {
                            if (ignoreVarCoupling && jVar != iVar) {
                                continue;
                            }
                            for (int iSpc = 0; iSpc < mask.m_StructuredNi0[iLoc][iVar].Length; iSpc++) {
                                for (int jSpc = 0; jSpc < mask.m_StructuredNi0[jLoc][jVar].Length; jSpc++) {
                                    if (ignoreSpecCoupling && jSpc != iSpc) {
                                        continue;
                                    }
                                    for (int iMode = 0; iMode < mask.m_StructuredNi0[iLoc][iVar][iSpc].Length; iMode++) {
                                        int SubRowIdx = mask.m_StructuredNi0[iLoc][iVar][iSpc][iMode].Si0;
                                        for (int jMode = 0; jMode < mask.m_StructuredNi0[jLoc][jVar][jSpc].Length; jMode++) {

                                            extNi0 RowNi0 = mask.m_StructuredNi0[iLoc][iVar][iSpc][iMode];
                                            extNi0 ColNi0 = mask.m_StructuredNi0[jLoc][jVar][jSpc][jMode];
                                            int Targeti0 = ListenToGlobalId? RowNi0.Gi0: RowNi0.Li0 + source._RowPartitioning.i0 - m_map.LocalLength;
                                            int Targetj0 = ColNi0.Gi0;

                                            var tmpBlock = MultidimensionalArray.Create(RowNi0.N, ColNi0.N);

                                            Debug.Assert(source.RowPartitioning.IsInLocalRange(Targeti0) && source.RowPartitioning.IsInLocalRange(Targetj0) && mask.GetType() == typeof(BlockMaskLoc) || mask.GetType() == typeof(BlockMaskExt));
                                            int SubColIdx = mask.m_StructuredNi0[jLoc][jVar][jSpc][jMode].Si0;
#if Debug
                                            
                                            SubMSR.ReadBlock(SubRowIdx, SubColIdx, tmpBlock);
                                            Debug.Assert(tmpBlock.Sum() == 0);
                                            Debug.Assert(tmpBlock.InfNorm() == 0);
#endif

                                            try {
                                                source.ReadBlock(Targeti0, Targetj0,
                                                tmpBlock);
                                            } catch (Exception e) {
                                                Console.WriteLine("row: " + Targeti0);
                                                Console.WriteLine("col: " + Targetj0);
                                                throw new Exception(e.Message);
                                            }
                                            Debug.Assert(SubRowIdx < SubMSR.RowPartitioning.LocalLength);
                                            Debug.Assert(SubColIdx < SubMSR.ColPartition.LocalLength);

                                            
                                            SubMSR.AccBlock(SubRowIdx, SubColIdx, 1.0, tmpBlock);
                                            //SubColIdx += mask.m_StructuredNi0[jLoc][jVar][jSpc][jMode].N;
                                        }
                                       
                                        //SubColIdx = 0;
                                        //SubRowIdx += mask.m_StructuredNi0[iLoc][iVar][iSpc][iMode].N;
                                    }
                                    //SubRowIdx = 0;
                                }
                            }
                        }
                    }
                    auxIdx++;
                }
                auxIdx++;
            }
        }


        /// <summary>
        /// Translates back your masked vector to the full one
        /// </summary>
        /// <typeparam name="V"></typeparam>
        /// <typeparam name="W"></typeparam>
        /// <param name="targetVector">output</param>
        /// <param name="accVector"></param>
        public void AccVecToFull<V, W>(W accVector, V targetVector)
            where V : IList<double>
            where W : IList<double> {

            if (m_includeExternalCells) {
                if (targetVector.Count() != GetLocalandExternalDOF())
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                if (accVector.Count() != BMLoc.LocalDOF + BMExt.LocalDOF)
                    throw new ArgumentException("accVector length is not equal to length of mask");
                targetVector.AccV(1.0, accVector, BMLoc.m_LocalMask, default(int[]));
                targetVector.AccV(1.0, accVector, BMExt.m_LocalMask, default(int[]),b_index_shift: BMLoc.LocalDOF);
            } else {
                if (targetVector.Count() != m_map.LocalLength)
                    throw new ArgumentException("Length of targetVector != Length of original");
                if (accVector.Count() != BMLoc.LocalDOF)
                    throw new ArgumentException("accVector length is not equal to length of mask");
                targetVector.AccV(1.0, accVector, BMLoc.m_LocalMask, default(int[]));
            }
        }

        /// <summary>
        /// Translates back your masked vector to the full one cellwise. i.e. Useful, if you solved blockwise
        /// </summary>
        /// <typeparam name="V"></typeparam>
        /// <typeparam name="W"></typeparam>
        /// <param name="targetVector">output</param>
        /// <param name="accVector"></param>
        /// <param name="iCell">starts with 0</param>
        public void AccVecCellwiseToFull<V,W>(W accVector, int iCell, V targetVector)
            where V : IList<double>
            where W : IList<double> {


            if (m_includeExternalCells) {
                if (targetVector.Count() != GetLocalandExternalDOF())
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                if (iCell >= BMLoc.m_StructuredNi0.Length+ BMExt.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                BlockMaskBase mask;
                if(iCell < m_map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells) {
                    mask = BMLoc;
                } else {
                    mask = BMExt;
                }
                AuxAcc(mask, accVector, iCell, targetVector);
            } else {
                if (targetVector.Count() != m_map.LocalLength)
                    throw new ArgumentException("Length of targetVector != Length of original");
                if (iCell >= BMLoc.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                AuxAcc(BMLoc, accVector, iCell, targetVector);
            }

        }

        private int GetLocalandExternalDOF() {
            int eCell = m_map.LocalNoOfBlocks + m_map.AggBasis[0].AggGrid.iLogicalCells.NoOfExternalCells - 1;
            int eVar = m_map.AggBasis.Length-1;
            int eN = m_map.GetLength(eCell)-1;
            return m_map.LocalUniqueIndex(eVar, eCell, eN)+1;
        }

        public double[] GetSubBlockVec( IList<double> extvector, IList<double> locvector)
        {
            if (BMExt == null)
                return GetSubBlockVec(locvector);
            int SubL = BMExt.LocalDOF + BMLoc.LocalDOF;
            var subvector = new double[SubL];
            int acc_offset = BMLoc.LocalDOF;
            int target_offset = m_map.LocalLength;
            if(BMExt.m_LocalMask!=null && BMExt.LocalDOF > 0)
                subvector.AccV(1.0, extvector, default(int[]), BMExt.m_LocalMask, acc_index_shift: acc_offset, b_index_shift: (-target_offset));
            subvector.AccV(1.0, locvector, default(int[]), BMLoc.m_LocalMask);
            return subvector;
        }

        public void AccVecToFull<U ,V, W>(W accVector, V extvector, U locvector)
            where V : IList<double>
            where W : IList<double>
            where U : IList<double> {
            locvector.AccV(1.0, accVector, BMLoc.m_LocalMask, default(int[]));
            if (BMExt != null) {
                int target_offset = m_map.LocalLength;
                int acc_offset = BMLoc.LocalDOF;
                if (BMExt.m_LocalMask != null && BMExt.LocalDOF > 0)
                    extvector.AccV(1.0, accVector, BMExt.m_LocalMask, default(int[]), acc_index_shift: (-target_offset), b_index_shift: acc_offset);
            }
        }

        private void AuxAcc<V, W>(BlockMaskBase mask, W accVector, int iCell, V targetVector)
            where V : IList<double>
            where W : IList<double> {
            int nCell = mask.GetType() == typeof(BlockMaskExt) ? iCell - m_map.LocalNoOfBlocks : iCell;
            if (nCell > mask.m_StructuredNi0.Length - 1)
                throw new ArgumentOutOfRangeException("iCell is greater than Cells in mask");
            if (accVector.Count() != mask.GetCellwiseLength(nCell))
                throw new ArgumentException("accVector length is not equal to length of mask");

            var Cidx = mask.GetCellwiseLocalidx(nCell);
            Debug.Assert(accVector.Count() == Cidx.Count());
            targetVector.AccV(1.0, accVector, Cidx, default(int[]));
        }

        /// <summary>
        /// Translates a full vector to a vector corresponding to matrix mask cellwise
        /// </summary>
        /// <param name="fullVector"></param>
        /// <param name="iCell"></param>
        /// <returns></returns>
        public double[] GetVectorCellwise(IList<double> fullVector, int iCell) {

            double[] tmp;
            
            if (m_includeExternalCells) {
                if (fullVector.Count() != GetLocalandExternalDOF())
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                if (iCell >= BMLoc.m_StructuredNi0.Length + BMExt.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                BlockMaskBase mask;
                if (iCell < m_map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells) {
                    mask = BMLoc;
                } else {
                    mask = BMExt;
                }
                tmp=GetAuxAccVec(mask, fullVector, iCell);
            } else {
                if (fullVector.Count() != m_map.LocalLength)
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                if (iCell >= BMLoc.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                tmp=GetAuxAccVec(BMLoc, fullVector, iCell);
            }

            return tmp;
        }

        private double[] GetAuxAccVec(BlockMaskBase mask, IList<double> fullVector, int iCell) {
            int nCell=mask.GetType() == typeof(BlockMaskExt)? iCell-m_map.LocalNoOfBlocks:iCell;
            double[] subVector = new double[mask.GetCellwiseLength(nCell)];
            var Cidx = mask.GetCellwiseLocalidx(nCell);
            Debug.Assert(subVector.Length == Cidx.Length);
            ArrayTools.GetSubVector<int[], int[], double>(fullVector, subVector, Cidx);
            return subVector;
        }


            /// <summary>
            /// Translates a full vector to a vector corresponding to matrix mask
            /// </summary>
            /// <param name="fullVector"></param>
        public double[] GetSubBlockVec(IList<double> fullVector) {
            double[] subVector;
            if (m_includeExternalCells) {
                if (fullVector.Count() != GetLocalandExternalDOF())
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                subVector =new double[BMLoc.LocalDOF+BMExt.LocalDOF];
                subVector.AccV(1.0, fullVector, default(int[]), BMLoc.m_LocalMask);
                subVector.AccV(1.0, fullVector,default(int[]), BMExt.m_LocalMask, acc_index_shift: BMLoc.LocalDOF);
            } else {
                if (fullVector.Count() != m_map.LocalLength)
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                subVector = new double[BMLoc.LocalDOF];
                subVector.AccV(1.0, fullVector, default(int[]), BMLoc.m_LocalMask);
            }
            return subVector;
        }




        public static void MatlabDebugging(BlockMask[] BMs, MultigridOperator Mop) {

            int NoOfSchwzBlocks = BMs.Length;
            int myMpisize = Mop.Mapping.MpiSize;
            int myMpiRank = Mop.Mapping.MpiRank;
            var MopMap = Mop.Mapping;

            ilPSP.Connectors.Matlab.BatchmodeConnector matlab;
            matlab = new ilPSP.Connectors.Matlab.BatchmodeConnector();


            List<BlockMsrMatrix> Blocks = new List<BlockMsrMatrix>();
            var BlkIdx_gI_lR = NoOfSchwzBlocks.ForLoop(b => new List<int>());
            var BlkIdx_gI_eR = NoOfSchwzBlocks.ForLoop(b => new List<int>());
            int[][] BlockIndices_Local = new int[NoOfSchwzBlocks][];
            int[][] BlockIndices_External = new int[NoOfSchwzBlocks][];

            int LocalI0 = MopMap.i0;
            for (int iPart = 0; iPart < NoOfSchwzBlocks; iPart++) {
                BlkIdx_gI_lR[iPart] = BMs[iPart].BMLoc.m_GlobalMask;
                BlkIdx_gI_eR[iPart] = BMs[iPart].BMExt.m_GlobalMask;
                var locallist = new List<int>();
                var extlist = new List<int>();
                foreach (int lIdx in BlkIdx_gI_lR[iPart]) {
                    locallist.Add(lIdx - LocalI0);
                }
                foreach (int eIdx in BlkIdx_gI_eR[iPart]) {
                    extlist.Add(eIdx - LocalI0);
                }
                BlockIndices_Local[iPart] = locallist.ToArray();
                BlockIndices_External[iPart] = extlist.ToArray();
                Blocks.Add(BMs[iPart].GetSubBlockMatrix(Mop.OperatorMatrix));
            }


            int globalBlockCounter = 0;
            for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                if (rankCounter == myMpiRank)
                    Debug.Assert(rank_NoBlks == NoOfSchwzBlocks);

                for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {
                    double[] vec;
                    if (rankCounter == myMpiRank) {
                        vec = ArrayTools.Cat(BlkIdx_gI_lR[iBlock], BlkIdx_gI_eR[iBlock]).Select(ii => ((double)(ii + 1))).ToArray();
                    } else {
                        vec = new double[0];
                    }

                    matlab.PutVector(vec, string.Format("BlockIdx{0}", globalBlockCounter));

                    globalBlockCounter++;
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                }

            }

            //matlab.PutSparseMatrix(Perm, "Perm");
            //matlab.PutSparseMatrix(ExternalRowsTemp, "ExternalRowsTemp");

            globalBlockCounter = 0;
            for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {
                    BlockMsrMatrix Block;
                    if (rankCounter == myMpiRank) {
                        Block = Blocks[iBlock];
                    } else {
                        Block = null;
                    }

                    matlab.PutSparseMatrix(Block, string.Format("Block{0}", globalBlockCounter));

                    globalBlockCounter++;
                    csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                }
            }


            Console.WriteLine("Matlab dir: " + matlab.WorkingDirectory);

            matlab.PutSparseMatrix(Mop.OperatorMatrix, "Full");
            int GlobalNoOfBlocks = NoOfSchwzBlocks.MPISum();



            for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                matlab.Cmd("BlockErr({0} + 1, 1) = norm( Block{0} - Full( BlockIdx{0}, BlockIdx{0} ), inf );", iGlbBlock);
            }

            Random rnd = new Random(myMpiRank);
            double[] testRHS = new double[MopMap.LocalLength];
            for (int i = 0; i < testRHS.Length; i++) {
                testRHS[i] = rnd.NextDouble();
            }
            matlab.PutVector(testRHS, "testRHS");

            MPIexchange<double[]> ResExchange = new MPIexchange<double[]>(Mop.Mapping, testRHS);
            ResExchange.TransceiveStartImReturn();
            ResExchange.TransceiveFinish(0.0);

            int offset = MopMap.LocalLength;

            int g = 0;
            for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {
                    double[] SubVec;
                    if (rankCounter == myMpiRank) {
                        int LL = BlockIndices_Local[iBlock].Length;
                        int LE;
                        if (BlockIndices_External[iBlock] != null) {
                            LE = BlockIndices_External[iBlock].Length;
                        } else {
                            LE = 0;
                        }
                        int L = LL + LE;

                        SubVec = new double[L];
                        for (int i = 0; i < LL; i++) {
                            SubVec[i] = testRHS[BlockIndices_Local[iBlock][i]];
                        }
                        if (LE > 0) {
                            for (int i = 0; i < LE; i++) {
                                SubVec[i + LL] = ResExchange.Vector_Ext[BlockIndices_External[iBlock][i] - offset];
                            }
                        }
                    } else {
                        SubVec = new double[0];
                    }

                    matlab.PutVector(SubVec, "SubVec" + g);

                    g++;
                }
            }

            for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                matlab.Cmd("RhsErr({0} + 1, 1) = norm( SubVec{0} - testRHS( BlockIdx{0} ), inf );", iGlbBlock);
            }

            double[] testX = new double[testRHS.Length];
            MPIexchangeInverse<double[]> XXExchange = new MPIexchangeInverse<double[]>(MopMap, testX);

            g = 0;
            for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {

                    if (rankCounter == myMpiRank) {
                        int LL = BlockIndices_Local[iBlock].Length;
                        int LE;
                        if (BlockIndices_External[iBlock] != null) {
                            LE = BlockIndices_External[iBlock].Length;
                        } else {
                            LE = 0;
                        }
                        int L = LL + LE;


                        for (int i = 0; i < LL; i++) {
                            testX[BlockIndices_Local[iBlock][i]] += (g + 1);
                        }
                        if (LE > 0) {
                            for (int i = 0; i < LE; i++) {
                                XXExchange.Vector_Ext[BlockIndices_External[iBlock][i] - offset] += (g + 1);
                            }
                        }
                    } else {
                        //nop
                    }

                    g++;
                }
            }
            XXExchange.TransceiveStartImReturn();
            XXExchange.TransceiveFinish(1.0);

            matlab.Cmd("testXref = zeros({0},1);", MopMap.TotalLength);
            for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                matlab.Cmd("testXref(BlockIdx{0},1) = testXref(BlockIdx{0},1) + ({0} + 1);", iGlbBlock);
            }

            matlab.PutVector(testX, "testX");
            matlab.Cmd("testXErr = norm(testX - testXref, inf);");

            MultidimensionalArray BlockErr = MultidimensionalArray.Create(GlobalNoOfBlocks, 1);
            MultidimensionalArray RhsErr = MultidimensionalArray.Create(GlobalNoOfBlocks, 1);
            MultidimensionalArray testXErr = MultidimensionalArray.Create(1, 1);

            matlab.GetMatrix(BlockErr, "BlockErr");
            matlab.GetMatrix(RhsErr, "RhsErr");
            matlab.GetMatrix(testXErr, "testXErr");

            matlab.Execute();

            for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                Console.WriteLine("Block #{0} Error (external? ) " + BlockErr[iGlbBlock, 0], iGlbBlock);
                Console.WriteLine("RHS #{0} Error " + RhsErr[iGlbBlock, 0], iGlbBlock);
                Debug.Assert(BlockErr[iGlbBlock, 0] == 0);
                Debug.Assert(RhsErr[iGlbBlock, 0] == 0);
            }

            Console.WriteLine("X Error " + testXErr[0, 0]);
            Debug.Assert(testXErr[0, 0] == 0.0);

            matlab.Dispose();

        }

        public static BlockMsrMatrix GetAllExternalRows(MultigridMapping map, BlockMsrMatrix M) {

            //var map = Mop.Mapping;
            //var op = Mop.OperatorMatrix;

            var extcells = map.AggGrid.iLogicalCells.NoOfExternalCells.ForLoop(i => i + map.LocalNoOfBlocks);

            var SBS = new SubBlockSelector(map);
            SBS.CellSelector(extcells, false);
            var AllExtMask = new BlockMaskExt(SBS);

            var ExternalRows_BlockI0 = AllExtMask.GetAllSubMatrixCellOffsets();
            var ExternalRows_BlockN = AllExtMask.GetAllSubMatrixCellLength();
            var ExternalRowsIndices = AllExtMask.m_GlobalMask;

            BlockPartitioning PermRow = new BlockPartitioning(ExternalRowsIndices.Count, ExternalRows_BlockI0, ExternalRows_BlockN, M.MPI_Comm, i0isLocal: true);

            BlockMsrMatrix Perm = new BlockMsrMatrix(PermRow, M._RowPartitioning);
            for (int iRow = 0; iRow < ExternalRowsIndices.Count; iRow++) {
                Debug.Assert(M._RowPartitioning.IsInLocalRange(ExternalRowsIndices[iRow]) == false);
                Perm[iRow + PermRow.i0, ExternalRowsIndices[iRow]] = 1;
            }

            Perm.SaveToTextFileSparseDebug("Perm");

            return BlockMsrMatrix.Multiply(Perm, M);
        }
    }

}

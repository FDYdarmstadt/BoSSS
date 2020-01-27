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
    public class SubBlockSelector : SubBlockSelectorBase
    {
        public SubBlockSelector(MultigridMapping map) : base(map) { }
    }

    public class BlockMask {
        private class BlockMaskLoc : BlockMaskBase {
            public BlockMaskLoc(SubBlockSelector sbs) : base(sbs, csMPI.Raw._COMM.SELF) { }

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
        }
        private class BlockMaskExt : BlockMaskBase {

            public BlockMaskExt(SubBlockSelector SBS) : base(SBS, SBS.GetMapping.MPI_Comm) { }

            protected override int m_NoOfCells {
                get {
                    return m_map.AggBasis[0].AggGrid.iLogicalCells.NoOfExternalCells;
                }
            }

            protected override int m_CellOffset {
                get {
                    return m_map.AggBasis[0].AggGrid.iLogicalCells.NoOfLocalUpdatedCells - 1;
                }
            }

            public BlockMsrMatrix GetExternalRows(BlockMsrMatrix target) {

                var ExternalRows_BlockI0 = base.GetAllSubMatrixCellOffsets();
                var ExternalRows_BlockN = base.GetAllSubMatrixCellLength();
                var ExternalRowsIndices = base.m_GlobalMask;
                var Mop = target;
                BlockPartitioning PermRow = new BlockPartitioning(ExternalRowsIndices.Count, ExternalRows_BlockI0, ExternalRows_BlockN, Mop.MPI_Comm, i0isLocal: true);

                // Remark: we use a permutation matrix for MPI-exchange of rows
                //  set   Perm[l,m] = I;
                //  then  ExternalRowsTemp[l,-]  =  Mop[m,-]

                BlockMsrMatrix Perm = new BlockMsrMatrix(PermRow, Mop._RowPartitioning);
                for (int iRow = 0; iRow < ExternalRowsIndices.Count; iRow++) {
                    Debug.Assert(Mop._RowPartitioning.IsInLocalRange(ExternalRowsIndices[iRow]) == false);
                    Perm[iRow + PermRow.i0, ExternalRowsIndices[iRow]] = 1;
                }
                return BlockMsrMatrix.Multiply(Perm, Mop);
            }
        }

        public BlockMask (SubBlockSelector sbs,bool includeExternalCells) {
            m_includeExternalCells = includeExternalCells; 
            BMLoc =new BlockMaskLoc(sbs);
            BMExt =new BlockMaskExt(sbs);
            m_map = sbs.GetMapping;
            if (includeExternalCells) {
                SetThisShitUp(new BlockMaskBase[] { BMLoc, BMExt });
            } else {
                SetThisShitUp(new BlockMaskBase[] { BMLoc });
            }
        }

        private void SetThisShitUp(BlockMaskBase[] masks) {
            List<int> tmpOffsetList = new List<int>();
            List<int> tmpLengthList = new List<int>();
            List<extNi0[][][]> tmpNi0 = null;
            foreach (var mask in masks) {
                tmpOffsetList.AddRange(mask.GetAllSubMatrixCellOffsets());
                tmpLengthList.AddRange(mask.GetAllSubMatrixCellLength());
                tmpNi0.AddRange(mask.m_StructuredNi0.ToList());
            }
            Debug.Assert(SubMatrixOffsets != null);
            Debug.Assert(SubMatrixLen != null);
            Debug.Assert(tmpNi0 != null);
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

        /// <summary>
        /// If you want nothing special. Take this one. If you want only diagonal block matrix choose one of the other methods
        /// </summary>
        /// <returns></returns>
        public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix target) {
            if (target == null)
                throw new ArgumentNullException();
            if (target.NoOfRows < BMLoc.LocalDOF)
                throw new ArgumentException();

            BlockMsrMatrix extBlock;
            if (m_includeExternalCells) {
                BlockPartitioning extBlocking = new BlockPartitioning(BMLoc.LocalDOF + BMExt.LocalDOF, SubMatrixOffsets, SubMatrixLen, csMPI.Raw._COMM.SELF);

                //make an extended block dummy to fit local and external blocks
                extBlock = new BlockMsrMatrix(extBlocking, extBlocking);

                //get the external rows via MPI exchange
                BlockMsrMatrix ExtRowsTmp = BMExt.GetExternalRows(target);

                int offset = BMLoc.m_GlobalMask.Count;
                var extBlockRows = BMExt.m_GlobalMask.Count.ForLoop(i => i + offset);
                var extBlockCols = BMExt.m_LocalMask.Count.ForLoop(i => i + offset);

                //add local Block ...
                target.WriteSubMatrixTo(extBlock, BMLoc.m_GlobalMask, default(int[]), BMLoc.m_GlobalMask, default(int[]));

                //add columns related to external rows ...
                target.AccSubMatrixTo(1.0, extBlock, BMLoc.m_GlobalMask, default(int[]), new int[0], default(int[]), BMExt.m_GlobalMask, extBlockCols);

                //add external rows ...
                ExtRowsTmp.AccSubMatrixTo(1.0, extBlock, BMExt.m_LocalMask, extBlockRows, BMLoc.m_GlobalMask, default(int[]), BMExt.m_GlobalMask, extBlockCols);
            } else {
                BlockPartitioning localBlocking = new BlockPartitioning(BMLoc.LocalDOF, SubMatrixOffsets, SubMatrixLen, csMPI.Raw._COMM.SELF, i0isLocal: true);

                extBlock = new BlockMsrMatrix(localBlocking);

                target.AccSubMatrixTo(1.0, extBlock, BMLoc.m_GlobalMask, default(int[]), BMLoc.m_GlobalMask, default(int[]));
            }
            Debug.Assert(extBlock!=null);
            return extBlock;
        }

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

            if (m_includeExternalCells) {
                var ExtRowsTmp = BMExt.GetExternalRows(target);
                Sblocks.AddRange(GetMaskBlocks(ExtRowsTmp, BMExt, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling));
            }

            Debug.Assert(Sblocks!=null);
            return Sblocks.ToArray();
        }


        private MultidimensionalArray[] GetMaskBlocks(BlockMsrMatrix target, BlockMaskBase mask, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {
            BlockMsrMatrix ExtRowsTmp = null;
            if (BMExt != null) {
                ExtRowsTmp = BMExt.GetExternalRows(target);
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
                                            int Targeti0 = RowNi0.Gi0;
                                            int Targetj0 = ColNi0.Gi0;
                                            int Subi0 = mask.GetRelativeSubBlockOffset(iLoc, iVar, iSpc, iMode);
                                            int Subj0 = mask.GetRelativeSubBlockOffset(jLoc, jVar, jSpc, jMode);
                                            int Subie = Subi0 + RowNi0.N - 1;
                                            int Subje = Subj0 + ColNi0.N - 1;

                                            var tmp = Sblocks[auxIdx].ExtractSubArrayShallow(new int[] { Subi0, Subj0 }, new int[] { Subie, Subje });

                                            Debug.Assert(m_map.IsInLocalRange(Targeti0) && m_map.IsInLocalRange(Targetj0) && mask.GetType() == typeof(BlockMaskLoc) || mask.GetType() == typeof(BlockMaskExt));


                                                try {
                                                    target.ReadBlock(Targeti0, Targetj0, tmp
                                                );
                                                } catch (Exception e) {
                                                    Console.WriteLine("row: " + Targeti0);
                                                    Console.WriteLine("col: " + Targetj0);
                                                    throw new Exception(e.Message);
                                                }

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
                var ExtRowsTmp = BMExt.GetExternalRows(target);
                GetMaskMatrix(submatrix, target, BMLoc, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling);
                GetMaskMatrix(ExtRowsTmp, target, BMExt, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling);
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


      
        private void GetMaskMatrix(BlockMsrMatrix submatrixempty, BlockMsrMatrix target, BlockMaskBase mask, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {

            var SubMSR = submatrixempty;

            int SubRowIdx = 0;
            int SubColIdx = 0;

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
                                        for (int jMode = 0; jMode < mask.m_StructuredNi0[jLoc][jVar][jSpc].Length; jMode++) {

                                            int[] RowIdx = new int[] { iLoc, iVar, iSpc, iMode };
                                            int[] ColIdx = new int[] { jLoc, jVar, jSpc, jMode };
                                            extNi0 RowNi0 = mask.m_StructuredNi0[iLoc][iVar][iSpc][iMode];
                                            extNi0 ColNi0 = mask.m_StructuredNi0[jLoc][jVar][jSpc][jMode];
                                            int Targeti0 = RowNi0.Gi0;
                                            int Targetj0 = ColNi0.Gi0;

                                            var tmpBlock = MultidimensionalArray.Create(RowNi0.N, ColNi0.N);

                                            Debug.Assert(m_map.IsInLocalRange(Targeti0) && m_map.IsInLocalRange(Targetj0) && mask.GetType() == typeof(BlockMaskLoc) || mask.GetType() == typeof(BlockMaskExt));

                                            try {
                                                target.ReadBlock(Targeti0, Targetj0,
                                                tmpBlock);
                                            } catch (Exception e) {
                                                Console.WriteLine("row: " + Targeti0);
                                                Console.WriteLine("col: " + Targetj0);
                                                throw new Exception(e.Message);
                                            }
                                            SubMSR.AccBlock(SubRowIdx, SubColIdx, 1, tmpBlock);
                                            SubColIdx += ColNi0.N;
                                        }
                                        SubRowIdx += mask.m_StructuredNi0[iLoc][iVar][iSpc][iMode].N;
                                    }
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
            if (targetVector.Count() != m_map.LocalLength)
                throw new ArgumentException("Length of targetVector != Length of original");

            if (m_includeExternalCells) {
                if (accVector.Count() != BMLoc.LocalDOF+BMExt.LocalDOF)
                    throw new ArgumentException("accVector length is not equal to length of mask");
                targetVector.AccV(1.0, accVector, BMLoc.m_LocalMask, default(int[]));
                targetVector.AccV(1.0, accVector, BMExt.m_LocalMask, default(int[]));
            } else {
                
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
        public void AccVecCellwiseToFull<V, W>(W accVector, int iCell, V targetVector)
            where V : IList<double>
            where W : IList<double> {
            if (targetVector.Count() != m_map.LocalLength)
                throw new ArgumentException("Length of targetVector not equal Length of original");

            if (m_includeExternalCells) {
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
                if (iCell >= BMLoc.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                AuxAcc(BMLoc, accVector, iCell, targetVector);
            }

        }

        private void AuxAcc<V, W>(BlockMaskBase mask, W accVector, int iCell, V targetVector)
            where V : IList<double>
            where W : IList<double> {
            if (iCell > mask.m_StructuredNi0.Length - 1)
                throw new ArgumentOutOfRangeException("iCell is greater than Cells in mask");
            if (accVector.Count() != mask.GetCellwiseLength(iCell))
                throw new ArgumentException("accVector length is not equal to length of mask");

            var Cidx = mask.GetCellwiseLocalidx(iCell);
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
            if (fullVector.Count() != m_map.LocalLength)
                throw new ArgumentException("Length of targetVector not equal Length of original");
            if (m_includeExternalCells) {
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
                if (iCell >= BMLoc.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                tmp=GetAuxAccVec(BMLoc, fullVector, iCell);
            }

            return tmp;
        }

        public double[] GetAuxAccVec(BlockMaskBase mask, IList<double> fullVector, int iCell) {
            double[] subVector = new double[mask.GetCellwiseLength(iCell)];
            var Cidx = mask.GetCellwiseLocalidx(iCell);
            Debug.Assert(subVector.Length == Cidx.Length);
            ArrayTools.GetSubVector<int[], int[], double>(fullVector, subVector, Cidx);
            return subVector;
        }


            /// <summary>
            /// Translates a full vector to a vector corresponding to matrix mask
            /// </summary>
            /// <param name="fullVector"></param>
        public double[] GetSubBlockVec(IList<double> fullVector) {
            if (fullVector.Count() != m_map.LocalLength)
                throw new ArgumentException("Length of targetVector not equal Length of original");

            double[] subVector;
            if (m_includeExternalCells) {
                subVector=new double[BMLoc.LocalDOF+BMExt.LocalDOF];
                subVector.AccV(1.0, fullVector, default(int[]), BMLoc.m_LocalMask);
                subVector.AccV(1.0, fullVector, default(int[]), BMExt.m_LocalMask);
            } else {
                subVector = new double[BMLoc.LocalDOF];
                subVector.AccV(1.0, fullVector, default(int[]), BMLoc.m_LocalMask);
            }
            return subVector;
        }










    }


}

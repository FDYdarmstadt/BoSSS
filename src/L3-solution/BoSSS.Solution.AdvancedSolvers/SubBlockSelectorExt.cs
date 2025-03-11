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
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers
{
    /// <summary>
    /// Enables the selection of sub blocks within a <see cref="MultigridMapping"/>,
    /// which can be specified by the developer.
    /// There are four selection types: cell, variable, species and dg block selection.
    /// Default: Selects all blocks.
    /// </summary>
    public class SubBlockSelector : SubBlockSelectorBase {
        
        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="map"></param>
        public SubBlockSelector(ICoordinateMapping map) : base(map) { }
    }


    /// <summary>
    /// Interface to apply masking onto matrix and vectors according to sub block selection.
    /// Used in linear solvers (e.g. <see cref="Schwarz"/>, <see cref="LevelPmg"/>, etc.).
    /// </summary>
    public class BlockMask {

        /// <summary>
        /// Auxiliary class, provides Local Block masks
        /// </summary>
        private class BlockMaskLoc : BlockMaskBase {
            public BlockMaskLoc(SubBlockSelector sbs) : base(sbs, csMPI.Raw._COMM.SELF) {
                base.GenerateAllMasks();
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

            protected override int m_SubBlockOffset {
                get {
                    return 0;
                }
            }
        }

        /// <summary>
        /// Auxiliary class, provides masking of external cells: overlap that is located on other MPI processes
        /// </summary>
        private class BlockMaskExt : BlockMaskBase {

            public BlockMaskExt(SubBlockSelector SBS, int LocMaskOffset) : base(SBS, SBS.Mapping.MPI_Comm) {
                
                // must be set before mask generation
                m_LocOffset = LocMaskOffset;
                m_extLocLen = GetLocalLength_Ext(SBS.Mapping);

                base.GenerateAllMasks();
                foreach (int idx in this.m_GlobalMask) {
                    Debug.Assert(idx < m_map.i0 || idx >= m_map.iE);
                }

                int LL = m_map.LocalLength;
                int jMax = m_map.LocalCellCount - 1;
                int LE = m_map.LocalUniqueIndex(0, jMax, 0) + m_map.GetLength(jMax);

                foreach (int idx in this.m_LocalMask) {
                    Debug.Assert(idx >= LL);
                    Debug.Assert(idx < LE);
                }
            }


            /// <summary>
            /// Gets DOF of ghost cells available on this proc
            /// </summary>
            /// <returns>DOF of ghost cells available on this proc
            /// </returns>
            static public int GetLocalLength_Ext(ICoordinateMapping cm) {
                int Jup = cm.NoOfLocalUpdatedCells;
                int Jtt = cm.LocalCellCount;
                int Len = 0;
                for(int j = Jup; j < Jtt; j++)
                    Len += cm.GetLength(j);

                //int[] LocCellIdxExt = this.AggGrid.iLogicalCells.NoOfExternalCells.ForLoop(i => i + Locoffset);
                //foreach(int jCell in LocCellIdxExt) {
                //    for(int fld = 0; fld < NoOfVariables; fld++) {
                //        Len += this.AggBasis[fld].GetLength(jCell, this.DgDegree[fld]);
                //    }
                //}
                return Len;
            }
            

            private int m_LocOffset;
            private int m_extLocLen;

            protected override int m_NoOfCells {
                get {
                    return m_map.NoOfExternalCells;
                }
            }

            protected override int m_CellOffset {
                get {
                    return m_map.NoOfLocalUpdatedCells;
                }
            }

            /// <summary>
            /// These are the DOF of external ghost cells available on this proc
            /// Note: this is method requires MPI communication
            /// and thus must be called on ALL procs
            /// </summary>
            protected override int m_LocalLength {
                get {
                    return m_extLocLen;
                }
            }

            protected override int m_SubBlockOffset {
                get {
                    return m_LocOffset;
                }
            }
        }

        /// <summary>
        /// generates a masking of subblocks within a <see cref="MultigridMapping"/> according to <see cref="SubBlockSelector"/> <paramref name="sbs"/>
        /// enables applying mask onto matrices and vectors, which comply with this <see cref="MultigridMapping"/>.
        /// (e.g. sub matrix generation, sub vector extraction, etc.).
        /// The smallest unit are dg blocks.
        /// The masking operates on local blocks available on this proc per default,
        /// which can be exceeded to external blocks by providing external rows: <paramref name="ExtRows"/>.
        /// Ghost cells (covert by <see cref="BoSSS.Foundation.Grid.ILogicalCellData.NoOfExternalCells"/>) can be acquired by <see cref="GetAllExternalRows"/>.
        /// </summary>
        /// <param name="sbs">sub block selection defined by dev</param>
        /// <param name="ExtRows">external rows collected from other MPI-processes on this proc</param>
        public BlockMask(SubBlockSelector sbs, BlockMsrMatrix ExtRows = null) {
            m_map = sbs.Mapping;
            m_ExtRows = ExtRows;
            m_ExchangeExternalRows = (ExtRows != null) && m_map.MpiSize > 1;
            BMLoc = new BlockMaskLoc(sbs);
            
            if (m_ExchangeExternalRows) {
                BMExt = new BlockMaskExt(sbs, BMLoc.LocalDOF);
                SetThisShitUp(new BlockMaskBase[] { BMLoc, BMExt });
            } else {
                SetThisShitUp(new BlockMaskBase[] { BMLoc });
            }
#if Debug
            CheckIndices();
#endif
        }

        private void CheckIndices() {
            for(int i = 0; i < BMLoc.m_GlobalMask.Count; i++) {
                long GlobIdx = BMLoc.m_GlobalMask.ToArray()[i];
                Debug.Assert(m_map.IsInLocalRange(GlobIdx));
            }
            for(int i = 0; i < BMExt.m_GlobalMask.Count; i++) {
                long GlobIdx = BMExt.m_GlobalMask.ToArray()[i];
                Debug.Assert(!m_map.IsInLocalRange(GlobIdx));
            }
        }

        /// <summary>
        /// true if no elements are selected;
        /// Typically some phatological use case, e.g. very coarse meshes
        /// </summary>
        public bool IsEmpty {
            get;
            private set;
        }

        private void SetThisShitUp(BlockMaskBase[] masks) {
            List<long> tmpOffsetList = new List<long>();
            List<int> tmpLengthList = new List<int>();
            List<extNi0[][][]> tmpNi0 = new List<extNi0[][][]>();
            foreach(var mask in masks) {
                var tmp = mask.GetAllSubMatrixCellOffsets().ToArray();
                tmpOffsetList.AddRange(tmp);
                tmpLengthList.AddRange(mask.GetAllSubMatrixCellLength());
                tmpNi0.AddRange(mask.m_StructuredNi0.ToList());
            }
            if (tmpOffsetList.Count == 0)
                IsEmpty = true; // typically some phatological use case, e.g. very coarse meshes
                //throw new ArgumentException("Nothing Selected. Mask is empty");
            Debug.Assert(tmpOffsetList != null);
            Debug.Assert(tmpLengthList != null);
            Debug.Assert(tmpNi0 != null);
            Debug.Assert(IsEmpty || tmpOffsetList.GroupBy(x => x).Any(g => g.Count() == 1));
            Debug.Assert(IsEmpty || tmpNi0.GroupBy(x => x).Any(g => g.Count() == 1));

            SubMatrixOffsets = tmpOffsetList;
            SubMatrixLen = tmpLengthList;
            StructuredNi0 = tmpNi0.ToArray();
        }

        // internal set

        /// <summary>
        /// Block mask of local entires yield by <see cref="BlockMaskLoc"/>
        /// </summary>
        BlockMaskLoc BMLoc;

        /// <summary>
        /// Block mask of external entires yield by <see cref="BlockMaskExt"/>
        /// </summary>
        BlockMaskExt BMExt;

        /// <summary>
        /// collected sub matrix offsets of <see cref="BMLoc"/> and <see cref="BMExt"/>
        /// </summary>
        List<long> SubMatrixOffsets;

        /// <summary>
        /// collected sub matrix lengths of <see cref="BMLoc"/> and <see cref="BMExt"/>
        /// </summary>
        List<int> SubMatrixLen;


        /// <summary>
        /// inter-process communication of matrix rows
        /// </summary>
        bool m_ExchangeExternalRows;

        /// <summary>
        /// structured Ni0 of <see cref="BMLoc"/> and <see cref="BMExt"/>
        /// </summary>
        extNi0[][][][] StructuredNi0;

        /// <summary>
        /// <see cref="ICoordinateMapping"/>, which this mask is based upon
        /// </summary>
        ICoordinateMapping m_map;

        /// <summary>
        /// external rows, which correspond to ghost cells, overgiven in cctor
        /// </summary>
        BlockMsrMatrix m_ExtRows;

        /// <summary>
        /// Gets number of blocks/cells covert by mask
        /// </summary>
        public int NoOfMaskedCells {
            get {
                if (m_ExchangeExternalRows) {
                    return BMLoc.m_StructuredNi0.Length + BMExt.m_StructuredNi0.Length;
                } else {
                    return BMLoc.m_StructuredNi0.Length;
                }
            }
        }

        /// <summary>
        /// Gets number of rows covert by mask
        /// </summary>
        public int NoOfMaskedRows {
            get {
                if (m_ExchangeExternalRows) {
                    return BMLoc.LocalDOF + BMExt.LocalDOF;
                } else {
                    return BMLoc.LocalDOF;
                }
            }
        }


        /// <summary>
        /// If you just want to get the <see cref="BlockMsrMatrix"/>, which corresponds to this <see cref="BlockMask"/>.
        /// This is the method to choose!
        /// </summary>
        /// <returns>sub-matrix on <see cref="IMPI_CommConstants.SELF"/></returns>
        public BlockMsrMatrix GetSubBlockMatrix_MpiSelf(BlockMsrMatrix source) {
            return GetSubBlockMatrix(source, csMPI.Raw._COMM.SELF);
        }

        /// <summary>
        /// If you just want to get the <see cref="BlockMsrMatrix"/>, which corresponds to this <see cref="BlockMask"/>.
        /// This is the method to choose! In addition, MPI communicator can be defined via <paramref name="comm"/>.
        /// </summary>
        /// <returns>sub-matrix on <paramref name="comm"/></returns>
        public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix source, MPI_Comm comm) {
            return GetSubBlockMatrix(source, this, comm);
        }

        /// <summary>
        /// If you just want to get the <see cref="BlockMsrMatrix"/>, which corresponds to this <see cref="BlockMask"/>.
        /// This is the method to choose! In addition, MPI communicator can be defined via <paramref name="comm"/>.
        /// </summary>
        /// <returns>sub-matrix on <paramref name="comm"/></returns>
        public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix source, BlockMask colMask, MPI_Comm comm) {
            if (source == null)
                throw new ArgumentNullException();
            if (source.NoOfRows < BMLoc.LocalDOF)
                throw new ArgumentException();

            BlockMsrMatrix target;
            if (m_ExchangeExternalRows) {

               BlockPartitioning targetBlocking = new BlockPartitioning(BMLoc.LocalDOF + BMExt.LocalDOF, SubMatrixOffsets, SubMatrixLen, comm);

                //make an extended block dummy to fit local and external blocks
                target = new BlockMsrMatrix(targetBlocking, targetBlocking);

                //get the external rows via MPI exchange
                var ExtRowsTmp = m_ExtRows;

                int offset = BMLoc.m_GlobalMask.Count;
                long[] extBlockRows = BMExt.m_GlobalMask.Count.ForLoop(i => (long)(i + offset));
                long[] extBlockCols = BMExt.m_GlobalMask.Count.ForLoop(i => (long)(i + offset));

                //ExtRowsTmp lives at the MPI-Communicator of the target, thus the global index is related to a new partitioning and has nothing to do with the partitioning of the multigrid operator ...
                
                long[] GlobalIdxExtRows = BMExt.m_GlobalMask.Count.ForLoop(i => (long)(BMExt.m_LocalMask[i] - m_map.LocalLength));
                for(int iGlob = 0; iGlob < GlobalIdxExtRows.Length; iGlob++) {
                    Debug.Assert(GlobalIdxExtRows[iGlob] < ExtRowsTmp._RowPartitioning.LocalLength);
                    GlobalIdxExtRows[iGlob] += ExtRowsTmp._RowPartitioning.i0;
                    Debug.Assert(ExtRowsTmp._RowPartitioning.IsInLocalRange(GlobalIdxExtRows[iGlob]));
                }

                //add local Block ...
                source.WriteSubMatrixTo(target, BMLoc.m_GlobalMask, default(long[]), BMLoc.m_GlobalMask, default(long[]));

                //add columns related to external rows ...
                source.AccSubMatrixTo(1.0, target, BMLoc.m_GlobalMask, default(long[]), new long[0], default(long[]), BMExt.m_GlobalMask, extBlockCols);

                //add external rows ...
                ExtRowsTmp.AccSubMatrixTo(1.0, target, GlobalIdxExtRows, extBlockRows, BMLoc.m_GlobalMask, default(long[]), BMExt.m_GlobalMask, extBlockCols);
            } else {
                BlockPartitioning rowPart = new BlockPartitioning(BMLoc.LocalDOF, SubMatrixOffsets, SubMatrixLen, comm, i0isLocal: true);
                BlockPartitioning colPart;
                if (colMask == this)
                    colPart = rowPart;
                else
                    colPart = new BlockPartitioning(colMask.BMLoc.LocalDOF, colMask.SubMatrixOffsets, colMask.SubMatrixLen, comm, i0isLocal: true);

                target = new BlockMsrMatrix(rowPart, colPart);

                //BMLoc.m_GlobalMask.SaveToTextFile("mask-" + BMLoc.m_GlobalMask.Count + ".txt");
                source.AccSubMatrixTo(1.0, target, BMLoc.m_GlobalMask, default(long[]), colMask.BMLoc.m_GlobalMask, default(long[]));
            }
            Debug.Assert(target != null);
            return target;
        }

        /// <summary>
        /// Get array of diagonal cell-blocks covert by <see cref="BlockMask"/>. With the ignore flags,
        /// coupling blocks can be left out (e.g. blocks containing level-set).
        /// - index i: block of i-th cell within mask (note: if some cell selection specified, i corresponds not to local cell index)
        /// - content: matrix corresponding to masking
        /// </summary>
        /// <param name="source">matrix to apply masking to</param>
        /// <param name="ignoreVarCoupling">flag to ignore variable coupling</param>
        /// <param name="ignoreSpecCoupling">flag to ignore species coupling</param>
        /// <returns></returns>
        public MultidimensionalArray[] GetDiagonalBlocks(BlockMsrMatrix source, bool ignoreVarCoupling, bool ignoreSpecCoupling) {

            bool ignoreCellCoupling = true; // consider only diagonal blocks

            List<MultidimensionalArray> Sblocks = new List<MultidimensionalArray>();
            Sblocks.AddRange(AuxGetSubBlocks(source, BMLoc, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling));

            if (m_ExchangeExternalRows && BMExt!=null) {
                //var ExtMatrix = this.GetExternalSubBlockMatrixOnly(target);
                //int RowOffset = m_map.iE;
                if (!ignoreCellCoupling)
                    throw new NotImplementedException("if overlapping is desired, only extraction of diagonal blocks is supported");
                Sblocks.AddRange(AuxGetSubBlocks(m_ExtRows, BMExt, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling));
            }

            Debug.Assert(Sblocks!=null);
            return Sblocks.ToArray();
        }

        private MultidimensionalArray[] AuxGetSubBlocks(BlockMsrMatrix source, BlockMaskBase mask, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {

            bool IsLocMask = mask.GetType() == typeof(BlockMaskLoc); // if external cells are masked, we have to consider other offsets ...

            int NoOfCells = mask.m_StructuredNi0.Length;
            int size = ignoreCellCoupling ? NoOfCells : NoOfCells * NoOfCells;

            MultidimensionalArray[] Sblocks = new MultidimensionalArray[size];

            int auxIdx = 0;
            for (int iLoc = 0; iLoc < mask.m_StructuredNi0.Length; iLoc++) {
                for (int jLoc = 0; jLoc < mask.m_StructuredNi0.Length; jLoc++) {
                    if (ignoreCellCoupling && jLoc != iLoc) {
                        continue;
                    }
                    int CellBlockLen = mask.GetLengthOfCell(jLoc);
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
                                            long Targeti0 = IsLocMask ? RowNi0.Gi0 : RowNi0.Li0 + source._RowPartitioning.i0 - m_map.LocalLength;
                                            long Targetj0 = ColNi0.Gi0;
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
        /// Get SubMatrix corresponding to this <see cref="BlockMask"/>.
        /// With the ignore flags, coupling blocks can be left out (e.g. blocks containing level-set).
        /// If <paramref name="ignoreCellCoupling"/> is set true, only diagonal blocks are considered.
        /// Probably slower than <see cref="GetSubBlockMatrix_MpiSelf(BlockMsrMatrix)"/>.
        /// </summary>
        /// <remarks>
        /// If you are using <paramref name="ignoreCellCoupling"/>, you may dismiss coupling with other cells.
        /// By comparison to variables, whose ordering is fixed in the mapping, this does not hold for species.
        /// E.g. iSpc=0 is not corresponding to the same species within every cell
        /// </remarks>
        /// <param name="source">matrix to apply masking to</param>
        /// <param name="ignoreCellCoupling">flag to ignore cell coupling</param>
        /// <param name="ignoreVarCoupling">flag to ignore variable coupling</param>
        /// <param name="ignoreSpecCoupling">flag to ignore species coupling</param>
        /// <returns>sub block matrix</returns>
        public BlockMsrMatrix GetSubBlockMatrix(BlockMsrMatrix source, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {

            BlockMsrMatrix submatrix = null;
            if (m_ExchangeExternalRows) {
                BlockPartitioning extBlocking = new BlockPartitioning(BMLoc.LocalDOF + BMExt.LocalDOF, SubMatrixOffsets, SubMatrixLen, csMPI.Raw._COMM.SELF);
                submatrix = new BlockMsrMatrix(extBlocking);
                var ExtRowsTmp = m_ExtRows;
                AuxGetSubBlockMatrix(submatrix, source, BMLoc, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling);
                AuxGetSubBlockMatrix(submatrix, ExtRowsTmp, BMExt, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling);
            } else {
                int Loclength = BMLoc.LocalDOF;
                var tmpN = BMLoc.GetAllSubMatrixCellLength();
                var tmpi0 = BMLoc.GetAllSubMatrixCellOffsets();
                BlockPartitioning localBlocking = new BlockPartitioning(Loclength, tmpi0.ToArray(), tmpN.ToArray(), csMPI.Raw._COMM.SELF, i0isLocal: true);
                submatrix = new BlockMsrMatrix(localBlocking);
                AuxGetSubBlockMatrix(submatrix, source, BMLoc, ignoreCellCoupling, ignoreVarCoupling, ignoreSpecCoupling);
            }
            Debug.Assert(submatrix != null);
            return submatrix;
        }


        private void AuxGetSubBlockMatrix(BlockMsrMatrix target, BlockMsrMatrix source, BlockMaskBase Rowmask, bool ignoreCellCoupling, bool ignoreVarCoupling, bool ignoreSpecCoupling) {

            bool IsLocalMask = Rowmask.GetType() == typeof(BlockMaskLoc);

            extNi0[][][][] RowNi0s = Rowmask.m_StructuredNi0;
            extNi0[][][][] ColNi0s = this.StructuredNi0;

            int auxIdx = 0;
            for (int iLoc = 0; iLoc < RowNi0s.Length; iLoc++) {
                for (int jLoc = 0; jLoc < ColNi0s.Length; jLoc++) {
                    if (ignoreCellCoupling && jLoc != iLoc) {
                        continue;
                    }
                    for (int iVar = 0; iVar < RowNi0s[iLoc].Length; iVar++) {
                        for (int jVar = 0; jVar < ColNi0s[jLoc].Length; jVar++) {
                            if (ignoreVarCoupling && jVar != iVar) {
                                continue;
                            }
                            for (int iSpc = 0; iSpc < RowNi0s[iLoc][iVar].Length; iSpc++) {
                                for (int jSpc = 0; jSpc < ColNi0s[jLoc][jVar].Length; jSpc++) {
                                    if (ignoreSpecCoupling && jSpc != iSpc) {
                                        continue;
                                    }
                                    for (int iMode = 0; iMode < RowNi0s[iLoc][iVar][iSpc].Length; iMode++) {
                                        int Trgi0 = RowNi0s[iLoc][iVar][iSpc][iMode].Si0;
                                        for (int jMode = 0; jMode < ColNi0s[jLoc][jVar][jSpc].Length; jMode++) {

                                            extNi0 RowNi0 = RowNi0s[iLoc][iVar][iSpc][iMode];
                                            extNi0 ColNi0 = ColNi0s[jLoc][jVar][jSpc][jMode];
                                            long Srci0 = IsLocalMask? RowNi0.Gi0: RowNi0.Li0 + source._RowPartitioning.i0 - m_map.LocalLength;
                                            long Srcj0 = ColNi0.Gi0;

                                            var tmpBlock = MultidimensionalArray.Create(RowNi0.N, ColNi0.N);
                                            int Trgj0 = ColNi0s[jLoc][jVar][jSpc][jMode].Si0;

                                            try {
                                                source.ReadBlock(Srci0, Srcj0,
                                                tmpBlock);
                                            } catch (Exception e) {
                                                Console.WriteLine("error at rank: " + m_map.MpiRank);
                                                Console.WriteLine("row: " + Srci0);
                                                Console.WriteLine("col: " + Srcj0);
                                                throw new Exception(e.Message);
                                            }
                                            Debug.Assert(Trgi0 < target.RowPartitioning.LocalLength);
                                            Debug.Assert(Trgj0 < target.ColPartition.LocalLength);

                                            target.AccBlock(Trgi0, Trgj0, 1.0, tmpBlock);
                                        } 
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
        /// accumulates the <paramref name="subVector"/>, which correspond to <see cref="BlockMask"/>, back to the <paramref name="fullVector"/>.
        /// </summary>
        /// <typeparam name="V"></typeparam>
        /// <typeparam name="W"></typeparam>
        /// <param name="fullVector">output, unmasked vector</param>
        /// <param name="subVector">input, masked vector</param>
        public void AccSubVec<V, W>(W subVector, V fullVector)
            where V : IList<double>
            where W : IList<double> {

            if (m_ExchangeExternalRows) {
                if (fullVector.Count != GetLocalAndExternalDOF(m_map))
                //    throw new ArgumentException("Length of targetVector not equal Length of original");
                if (subVector.Count != BMLoc.LocalDOF + BMExt.LocalDOF)
                    throw new ArgumentException("accVector length is not equal to length of mask");
                fullVector.AccV(1.0, subVector, BMLoc.m_LocalMask, default(int[]));
                fullVector.AccV(1.0, subVector, BMExt.m_LocalMask, default(int[]), b_index_shift: BMLoc.LocalDOF);
            } else {
                if (fullVector.Count != m_map.LocalLength)
                    throw new ArgumentException("Length of targetVector != Length of original");
                if (subVector.Count != BMLoc.LocalDOF)
                    throw new ArgumentException("accVector length is not equal to length of mask");
                fullVector.AccV(1.0, subVector, BMLoc.m_LocalMask, default(int[]));
            }
        }

        /// <summary>
        /// accumulates the <paramref name="subVector"/>, which correspond to <see cref="BlockMask"/>,
        /// back to the unmasked full vector.
        /// Entries of <paramref name="subVector"/> are distributed to the local (<paramref name="locFullVector"/>)
        /// and external (<paramref name="extFullVector"/>) part of the unmasked vector.
        /// </summary>
        /// <typeparam name="U"></typeparam>
        /// <typeparam name="V"></typeparam>
        /// <typeparam name="W"></typeparam>
        /// <param name="subVector">input, sub vector to accumulate</param>
        /// <param name="extFullVector">output, external part of unmasked vector</param>
        /// <param name="locFullVector">output, local part of unmasked vector</param>
        public void AccSubVec<U, V, W>(W subVector, V extFullVector, U locFullVector)
            where V : IList<double>
            where W : IList<double>
            where U : IList<double> {
            if (BMLoc != null && BMLoc.LocalDOF > 0)
                locFullVector.AccV(1.0, subVector, BMLoc.m_LocalMask, default(int[]));
            if (BMExt != null && BMExt.LocalDOF > 0) {
                int target_offset = m_map.LocalLength;
                int acc_offset = BMLoc.LocalDOF;
                if (BMExt.m_LocalMask != null && BMExt.LocalDOF > 0)
                    extFullVector.AccV(1.0, subVector, BMExt.m_LocalMask, default(int[]), acc_index_shift: (-target_offset), b_index_shift: acc_offset);
            }
        }

        /// <summary>
        /// accumulates masked vector <paramref name="subVector"/>
        /// of cell block <paramref name="iBlock"/>
        /// back to the full vector <paramref name="fullVector"/>. 
        /// <paramref name="iBlock"/> is the consecutive index of cellblocks,
        /// which are yield by <see cref="GetDiagonalBlocks(BlockMsrMatrix, bool, bool)"/>.
        /// note: if some cell selection specified, i corresponds not to local cell index
        /// </summary>
        /// <typeparam name="V"></typeparam>
        /// <typeparam name="W"></typeparam>
        /// <param name="fullVector">output, unmaksed vector</param>
        /// <param name="subVector">input, masked vector</param>
        /// <param name="iBlock">cell-block index</param>
        public void AccSubVecOfCell<V,W>(W subVector, int iBlock, V fullVector)
            where V : IList<double>
            where W : IList<double> {


            if (m_ExchangeExternalRows) {
                if (fullVector.Count != GetLocalAndExternalDOF(m_map))
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                if (iBlock >= BMLoc.m_StructuredNi0.Length + BMExt.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                BlockMaskBase mask;
                if(iBlock < BMLoc.m_StructuredNi0.Length) {
                    mask = BMLoc;
                } else {
                    mask = BMExt;
                }
                AuxAcc(mask, subVector, iBlock, fullVector);
            } else {
                if (fullVector.Count != m_map.LocalLength)
                    throw new ArgumentException("Length of targetVector != Length of original");
                if (iBlock >= BMLoc.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iBlock is greater than number of cellblocks in this mask");
                AuxAcc(BMLoc, subVector, iBlock, fullVector);
            }

        }



        private void AuxAcc<V, W>(BlockMaskBase mask, W subVector, int iCell, V fullVector)
            where V : IList<double>
            where W : IList<double> {
            int nCell = mask.GetType() == typeof(BlockMaskExt) ? iCell - BMLoc.m_StructuredNi0.Length : iCell;
            if (nCell >= mask.m_StructuredNi0.Length)
                throw new ArgumentOutOfRangeException("iCell is greater than Cells in mask");
            if (subVector.Count != mask.GetLengthOfCell(nCell))
                throw new ArgumentException("accVector length is not equal to length of mask");

            var Cidx = mask.GetLocalidcOfCell(nCell);
            Debug.Assert(subVector.Count == Cidx.Length);
            fullVector.AccV(1.0, subVector, Cidx, default(int[]));
        }

        /// <summary>
        /// returns the sub-vector of <paramref name="fullVector"/>
        /// corresponding to <see cref="BlockMask"/> of cell block <paramref name="iBlock"/>.
        /// <paramref name="iBlock"/> is the consecutive index of cell-blocks,
        /// which are yield by <see cref="GetDiagonalBlocks(BlockMsrMatrix, bool, bool)"/>.
        /// note: if some cell selection specified, i corresponds not to local cell index
        /// </summary>
        /// <param name="fullVector">input, unmasked vector</param>
        /// <param name="iBlock">cell-block index</param>
        /// <returns>sub-vector of cell block <paramref name="iBlock"/></returns>
        public double[] GetSubVecOfCell<T>(T fullVector, int iBlock) where T : IList<double> {

            double[] tmp;

            if(m_ExchangeExternalRows) {
                if(fullVector.Count != GetLocalAndExternalDOF(m_map))
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                if(iBlock >= BMLoc.m_StructuredNi0.Length + BMExt.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                BlockMaskBase mask;
                if(iBlock < BMLoc.m_StructuredNi0.Length) {
                    mask = BMLoc;
                } else {
                    mask = BMExt;
                }
                tmp = GetAuxAccVec(mask, fullVector, iBlock);
            } else {
                if(fullVector.Count != m_map.LocalLength)
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                if(iBlock >= BMLoc.m_StructuredNi0.Length)
                    throw new ArgumentOutOfRangeException("iCell is greater than Cellindex of mask");
                tmp = GetAuxAccVec(BMLoc, fullVector, iBlock);
            }
            return tmp;
        }

        private double[] GetAuxAccVec<T>(BlockMaskBase mask, T fullVector, int iCell) where T : IList<double> {
            int nCell = mask.GetType() == typeof(BlockMaskExt) ? iCell - BMLoc.m_StructuredNi0.Length : iCell;
            int L = mask.GetLengthOfCell(nCell);
            double[] subVector = new double[L];
            var Cidx = mask.GetLocalidcOfCell(nCell);
            Debug.Assert(subVector.Length == Cidx.Length);
            for(int l = 0; l < L; l++)
                subVector[l] = fullVector[Cidx[l]];
            return subVector;
        }


        /// <summary>
        /// returns the subvector of <paramref name="fullVector"/>
        /// corresponding to <see cref="BlockMask"/>.
        /// </summary>
        /// <param name="fullVector">input, unmasked vector</param>
        public double[] GetSubVec<T>(T fullVector) where T : IList<double> {
            double[] subVector;
            if(m_ExchangeExternalRows) {
                if(fullVector.Count != GetLocalAndExternalDOF(m_map))
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                subVector = new double[BMLoc.LocalDOF + BMExt.LocalDOF];
                subVector.AccV(1.0, fullVector, default(int[]), BMLoc.m_LocalMask);
                subVector.AccV(1.0, fullVector, default(int[]), BMExt.m_LocalMask, acc_index_shift: BMLoc.LocalDOF);
            } else {
                if(fullVector.Count != m_map.LocalLength)
                    throw new ArgumentException("Length of targetVector not equal Length of original");
                subVector = new double[BMLoc.LocalDOF];
                subVector.AccV(1.0, fullVector, default(int[]), BMLoc.m_LocalMask);
            }
            return subVector;
        }

        /// <summary>
        /// returns the subvector
        /// corresponding to <see cref="BlockMask"/>.
        /// Entries of subvector are taken from the local (<paramref name="locfullVector"/>)
        /// and external (<paramref name="extFullVector"/>) part of the unmasked vector.
        /// </summary>
        /// <param name="extFullVector">input, external part of unmasked vector</param>
        /// <param name="locfullVector">input, local part of unmasked vector</param>
        /// <returns></returns>
        public double[] GetSubVec(IList<double> extFullVector, IList<double> locfullVector) {
            if (BMExt == null)
                return GetSubVec(locfullVector);
            int SubL = BMExt.LocalDOF + BMLoc.LocalDOF;
            var subvector = new double[SubL];
            int acc_offset = BMLoc.LocalDOF;
            int target_offset = m_map.LocalLength;
            if (BMExt != null && BMExt.LocalDOF > 0)
                subvector.AccV(1.0, extFullVector, default(int[]), BMExt.m_LocalMask, acc_index_shift: acc_offset, b_index_shift: (-target_offset));
            if (BMLoc != null && BMLoc.LocalDOF > 0)
                subvector.AccV(1.0, locfullVector, default(int[]), BMLoc.m_LocalMask);
            return subvector;
        }


        /// <summary>
        /// returns all external rows of <paramref name="M"/>
        /// corresponding to ghost cells of <paramref name="map"/>,
        /// which are located on other MPI-ranks.
        /// </summary>
        /// <param name="map">Multigrid mapping</param>
        /// <param name="M">matrix distributed according to <paramref name="map"/></param>
        /// <returns></returns>
        /// <remarks>
        /// Exchange of matrix rows between MPI processors is implemented using multiplication with a permutation matrix.
        /// In this way, the MPI-communication routines of <see cref="BlockMsrMatrix.Multiply(BlockMsrMatrix, BlockMsrMatrix)"/> can be re-used.
        /// </remarks>
        public static BlockMsrMatrix GetAllExternalRows(MultigridMapping map, BlockMsrMatrix M) {
            using(new FuncTrace()) {
                var extcells = map.AggGrid.iLogicalCells.NoOfExternalCells.ForLoop(i => i + map.LocalNoOfBlocks);

                var SBS = new SubBlockSelector(map);
                SBS.CellSelector(extcells, false);
                var AllExtMask = new BlockMaskExt(SBS, 0);

                var ExternalRows_BlockI0 = AllExtMask.GetAllSubMatrixCellOffsets();
                var ExternalRows_BlockN = AllExtMask.GetAllSubMatrixCellLength();
                var ExternalRowsIndices = AllExtMask.m_GlobalMask;

                BlockPartitioning PermRow = new BlockPartitioning(ExternalRowsIndices.Count, ExternalRows_BlockI0, ExternalRows_BlockN, M.MPI_Comm, i0isLocal: true);

                BlockMsrMatrix Perm = new BlockMsrMatrix(PermRow, M._RowPartitioning);
                for(int iRow = 0; iRow < ExternalRowsIndices.Count; iRow++) {
                    Debug.Assert(M._RowPartitioning.IsInLocalRange(ExternalRowsIndices[iRow]) == false);
                    Perm[iRow + PermRow.i0, ExternalRowsIndices[iRow]] = 1;
                }

#if TEST
                Perm.SaveToTextFileSparseDebug("Perm");
#endif
                return BlockMsrMatrix.Multiply(Perm, M);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public static int GetLocalAndExternalDOF(ICoordinateMapping map) {
            int eCell = map.LocalNoOfBlocks + map.NoOfExternalCells - 1;
            return map.LocalUniqueIndex(0, eCell, 0) + map.GetLength(eCell);
        }

        #region stuff for Operator Testing
        /// <summary>
        /// Extracts blocks from matrix <paramref name="source"/>, i.e. one block per cell and the coupling between these cells.
        /// This might only work for local available cells
        /// </summary>
        public MultidimensionalArray[,] GetFullSubBlocks(BlockMsrMatrix source, bool ignoreVarCoupling = false, bool ignoreSpecCoupling = false) {
            int NoOfCells = StructuredNi0.Length;
            MultidimensionalArray[,] Sblocks = new MultidimensionalArray[NoOfCells, NoOfCells];

            for (int iLoc = 0; iLoc < StructuredNi0.Length; iLoc++) {
                for (int jLoc = 0; jLoc < StructuredNi0.Length; jLoc++) {

                    Sblocks[iLoc, jLoc] = GetBlock(source, ignoreVarCoupling, ignoreSpecCoupling, iLoc, jLoc);
                }
            }
            return Sblocks;
        }

        private MultidimensionalArray GetBlock(BlockMsrMatrix target, bool ignoreVarCoupling, bool ignoreSpecCoupling, int iLoc, int jLoc) {
            var _Sblocks = MultidimensionalArray.Create(BMLoc.GetLengthOfCell(iLoc), BMLoc.GetLengthOfCell(jLoc));


            for (int iVar = 0; iVar < StructuredNi0[iLoc].Length; iVar++) { // loop over (row/codomain/test) variables
                for (int jVar = 0; jVar < StructuredNi0[jLoc].Length; jVar++) { // loop over (column/domain/trial) variables
                    if (ignoreVarCoupling && jVar != iVar) {
                        continue;
                    }
                    for (int iSpc = 0; iSpc < StructuredNi0[iLoc][iVar].Length; iSpc++) { // loop over species 
                        for (int jSpc = 0; jSpc < StructuredNi0[jLoc][jVar].Length; jSpc++) {
                            if (ignoreSpecCoupling && jSpc != iSpc) {
                                continue;
                            }
                            for (int iMode = 0; iMode < StructuredNi0[iLoc][iVar][iSpc].Length; iMode++) {
                                for (int jMode = 0; jMode < StructuredNi0[jLoc][jVar][jSpc].Length; jMode++) {
                                    extNi0 RowNi0 = StructuredNi0[iLoc][iVar][iSpc][iMode];
                                    extNi0 ColNi0 = StructuredNi0[jLoc][jVar][jSpc][jMode];
                                    long Targeti0 = RowNi0.Gi0;
                                    long Targetj0 = ColNi0.Gi0;
                                    int Subi0 = BMLoc.GetRelativeSubBlockOffset(iLoc, iVar, iSpc, iMode);
                                    int Subj0 = BMLoc.GetRelativeSubBlockOffset(jLoc, jVar, jSpc, jMode);
                                    int Subie = Subi0 + RowNi0.N - 1;
                                    int Subje = Subj0 + ColNi0.N - 1;

                                    target.ReadBlock(Targeti0, Targetj0,
                                        _Sblocks.ExtractSubArrayShallow(new int[] { Subi0, Subj0 }, new int[] { Subie, Subje }));
                                }
                            }
                        }
                    }
                }
            }

            return _Sblocks;
        }

        /// <summary>
        /// global (i.e. across all MPI processors) indices of vector entries, resp. matrix rows and columns to select.
        /// </summary>
        public long[] GlobalIndices {
            get {
                List<long> tmp = new List<long>();
                if (m_ExchangeExternalRows) {
                    if (BMLoc != null)
                        tmp.AddRange(BMLoc.m_GlobalMask);
                    if (BMExt != null)
                        tmp.AddRange(BMExt.m_GlobalMask);
                    if (BMLoc == null && BMExt == null)
                        throw new Exception("empty mask. what a waste of time!");
                } else {
                    if (BMLoc != null)
                        tmp.AddRange(BMLoc.m_GlobalMask);
                    else
                        throw new Exception("empty mask. what a waste of time!");
                }
                return tmp.ToArray();

            }
        }

        /// <summary>
        /// local (i.e. valid only on current MPI processor) indices of vector entries, resp. matrix rows and columns to select.
        /// </summary>
        public int[] LocalIndices {
            get {
                List<int> tmp = new List<int>();
                if(m_ExchangeExternalRows) {
                    if (BMLoc != null)
                        tmp.AddRange(BMLoc.m_LocalMask);
                    if (BMExt != null)
                        tmp.AddRange(BMExt.m_LocalMask);
                    if (BMLoc == null && BMExt == null)
                        throw new Exception("empty mask. what a waste of time!");
                } else {
                    if (BMLoc != null)
                        tmp.AddRange(BMLoc.m_LocalMask);
                    else
                        throw new Exception("empty mask. what a waste of time!");
                }

                return tmp.ToArray();
            }
        }

        /// <summary>
        /// number of indices on the current processor
        /// </summary>
        public int LocalLength {
            get {
                return BMLoc.m_LocalMask.Count;
            }
        }


        // ==========
        // Stuff dedicated to testing ...
        // ==========

        /// <summary>
        /// This is provided for testing or if you know what you are doing!
        /// </summary>
        public List<long> GlobalIndices_Internal {
            get {
                return BMLoc.m_GlobalMask;
            }
        }

        /// <summary>
        /// This is provided for testing or if you know what you are doing!
        /// </summary>
        public List<long> GlobalIndices_External {
            get {
                return BMExt.m_GlobalMask;
            }
        }

#endregion
    }

}

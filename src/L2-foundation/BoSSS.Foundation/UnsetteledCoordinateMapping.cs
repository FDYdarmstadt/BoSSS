/* =======================================================================
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
using System.Linq;
using BoSSS.Foundation.Grid;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using MPI.Wrappers;
using BoSSS.Platform;
using ilPSP.Utils;

namespace BoSSS.Foundation {

   


    /// <summary>
    /// This class provides bijective mappings:
    /// <list type="bullet">
    ///   <item>
    ///    First, between <em>local field coordinate indices</em> and <em>local unique indices</em>,
    ///    by the methods <see cref="LocalUniqueCoordinateIndex"/> and <see cref="LocalFieldCoordinateIndex"/>, 
    ///    and ...
    ///   </item>
    ///   <item>
    ///    second, between <em>local unique indices</em> and <em>global unique indices</em>,
    ///    by the methods <see cref="Global2LocalIndex"/> and <see cref="Local2GlobalIndex"/>.
    ///   </item>
    /// </list>
    /// In easy words, that means it maps the DG coordinate from a list of fields into one long, one-dimensional vector
    /// and it converts between indices that are valid only on the actual MPI process (local indices)
    /// and indices that are valid among all MPI processes in the current MPI communicator (global indices).
    /// </summary>
    /// <remarks>
    /// A "<em>local index</em>" is an index that is valid only on the actual MPI process,
    /// while a "<em>global index</em>" is an index that is valid among all MPI processes 
    /// (in the corresponding MPI communicator).<br/>
    /// A "<em>local field coordinate index</em>" is a tuple (<i>f</i>,<i>j</i>,<i>n</i>), where...<br/>
    /// <list type="bullet">
    ///   <item><i>f</i> is a <see cref="DGField"/>,</item>
    ///   <item><i>j</i> is a local cell index and </item>
    ///   <item><i>n</i> is a polynomial index;</item>
    /// </list>
    /// A "<em>local unique index</em>" is a single integer that is unique for every field in this mapping, 
    /// for every local cell and for every polynomial;
    /// <br/>
    /// In this class, the DG fields are not specified. Hence, only a list of
    /// <see cref="Basis"/>-objects is provided (<see cref="BasisS"/>) to do the index computation.
    /// Any list of fields with equal DG basis can be accessed;
    /// If the mapping should be limited to a specific list of DG fields, the <see cref="CoordinateMapping"/>
    /// class can be used.
    /// </remarks>
    public class UnsetteledCoordinateMapping : Partitioning, IBlockPartitioning {

        /// <summary>
        /// Constructs an empty mapping.
        /// </summary>
        public UnsetteledCoordinateMapping(IGridData grd)
            : this(grd, new Basis[0], grd.CellPartitioning.MPI_Comm) {
        }

        /// <summary>
        /// Constructs a new mapping from an ordered list of basis functions;
        /// </summary>
        /// <param name="_basis">the list of DG basis'es that define this mapping</param>
        public UnsetteledCoordinateMapping(params Basis[] _basis)
            : this((IEnumerable<Basis>)_basis) {
        }

        static int ComputeLength(IEnumerable<Basis> _basis) {
            if (_basis == null || _basis.Count() <= 0)
                return 0;
            else 
                return _basis.Sum(b => b.MaximalLength) * _basis.ElementAt(0).GridDat.iLogicalCells.NoOfLocalUpdatedCells;
        }

        /// <summary>
        /// Constructs a new mapping from an ordered list of basis functions;
        /// </summary>
        /// <param name="_basis">the list of DG basis'es that define this mapping</param>
        public UnsetteledCoordinateMapping(IEnumerable<Basis> _basis) :
            this(_basis.ElementAt(0).GridDat, _basis, _basis.ElementAt(0).GridDat.CellPartitioning.MPI_Comm) { }

        /// <summary>
        /// Constructs a new mapping from an ordered list of basis functions;
        /// </summary>
        /// <param name="_basis">the list of DG basis'es that define this mapping</param>
        private UnsetteledCoordinateMapping(IGridData g, IEnumerable<Basis> _basis, MPI_Comm _Comm) :
            base(ComputeLength(_basis), _Comm) //
        {

            // =========================
            // initial checks
            // =========================

            Basis[] basis = _basis.ToArray();
            m_Context = g;
            for (int i = 0; i < basis.Length; i++) {
                if (!object.ReferenceEquals(m_Context, basis[i].GridDat))
                    throw new ArgumentException("all basis-objects must be associated with the same grid.");
            }
            m_BasisS = basis;

             // cell - independent index calculation
            //int[] j0CoordinateIndexMin = new int[m_BasisS.Length + 1];
            //int[] j0CoordinateIndexMax = new int[m_BasisS.Length + 1];
            m_j0CoordinateIndex = new int[m_BasisS.Length + 1];
            m_MaxTotalNoOfCoordinatesPerCell = 0;
            m_MinTotalNoOfCoordinatesPerCell = 0;
            
            for (int i = 0; i < m_BasisS.Length; i++) {
                m_j0CoordinateIndex[i] = m_MaxTotalNoOfCoordinatesPerCell;
                m_MaxTotalNoOfCoordinatesPerCell += m_BasisS[i].MaximalLength;
                m_MinTotalNoOfCoordinatesPerCell += m_BasisS[i].MinimalLength;
            }
            m_j0CoordinateIndex[m_BasisS.Length] = m_MaxTotalNoOfCoordinatesPerCell;


            // ==================
            // sub-blocking setup
            // ==================

            m_SubblkLen = new List<int[]>();
        }

        BlockPartitioning GetBlocking(bool bConstantFrame) {
            int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            Debug.Assert(J == GridDat.CellPartitioning.LocalLength);

            Basis[] basisS = this.m_BasisS;
            int NoOfBasisS = basisS.Length;

            int[] BlockType = new int[J];
            List<int[]> _SubblkLen = new List<int[]>();
            List<int[]> _Subblk_i0 = new List<int[]>();
            int[] SubblkLen_j = null;
            int[] Subblk_i0_j = null;
            for(int j = 0; j < J; j++) {
                if (SubblkLen_j == null)
                    SubblkLen_j = new int[NoOfBasisS];
                if (Subblk_i0_j == null)
                    Subblk_i0_j = new int[NoOfBasisS];

                int i0 = 0;
                for(int iBs = 0; iBs < NoOfBasisS; iBs++) {
                    Debug.Assert(i0 == m_j0CoordinateIndex[iBs]);
                    Subblk_i0_j[iBs] = m_j0CoordinateIndex[iBs];
                    SubblkLen_j[iBs] = basisS[iBs].GetLength(j);
                    i0 += basisS[iBs].MaximalLength;
                }

                Debug.Assert(_SubblkLen.Count == _Subblk_i0.Count);
                int K = _SubblkLen.Count;
                int blockType;
                bool bFound = false;
                for(blockType = 0; blockType < K; blockType++) { // loop over all block types found so far
                    int[] len_bT = _SubblkLen[blockType];
                    int[] _i0_bT = _Subblk_i0[blockType];
                    Debug.Assert(len_bT.Length == NoOfBasisS);
                    Debug.Assert(_i0_bT.Length == NoOfBasisS);
                    Debug.Assert(ArrayTools.ListEquals<int>(_i0_bT, m_j0CoordinateIndex.GetSubVector(0, NoOfBasisS)));

                    bFound = true;                    
                    for(int i = 0; i < NoOfBasisS; i++) {
                        if(len_bT[i] != SubblkLen_j[i]) {
                            bFound = false;
                            break;
                        }
                        if(_i0_bT[i] != Subblk_i0_j[i]) {
                            bFound = false;
                            break;
                        }
                    }
                                       

                    if (bFound)
                        break;
                }

                if(bFound == false) {
                    Debug.Assert(blockType == _SubblkLen.Count);
                    Debug.Assert(blockType == _Subblk_i0.Count);

                    _Subblk_i0.Add(Subblk_i0_j.CloneAs());
                    _SubblkLen.Add(SubblkLen_j.CloneAs());
                }

                BlockType[j] = blockType;
            }

            return new BlockPartitioning(
                this.LocalLength,
                bConstantFrame ? this.MaxTotalNoOfCoordinatesPerCell : -1,
                _Subblk_i0.ToArray(), _SubblkLen.ToArray(), 
                BlockType,
                this.MPI_Comm);
        }


        List<int[]> m_SubblkLen;

        public int[] GetSubblkLen(int blockType) {
            return m_SubblkLen[blockType];
        }
        
#if DEBUG
        int[] m_j0Subblk_i0_Backup = null;
#endif
        int[] m_j0Subblk_i0 = null;

        public int[] GetSubblk_i0(int blockType) {
            if(m_j0Subblk_i0 == null) {
                m_j0Subblk_i0 = ArrayTools.GetSubVector(this.m_j0CoordinateIndex, 0, this.m_BasisS.Length);
            }

#if DEBUG
            if (m_j0Subblk_i0_Backup == null) {
                m_j0Subblk_i0_Backup = m_j0Subblk_i0.CloneAs();
            } else {
                Debug.Assert(ArrayTools.ListEquals<int>(m_j0Subblk_i0, m_j0Subblk_i0_Backup), "Some moron messed with the start indices.");
            }
#endif
            return m_j0Subblk_i0;
        }

        /// <summary>
        /// Equal to total number of cells, i.e. each cell corresponds to one block, see also
        /// <see cref="IBlockPartitioning.TotalNoOfBlocks"/>.
        /// </summary>
        public int TotalNoOfBlocks {
            get {
                return m_Context.CellPartitioning.TotalLength;
            }
        }

        /// <summary>
        /// Equal to local number of cells, i.e. each cell corresponds to one block, see also
        /// <see cref="IBlockPartitioning.TotalNoOfBlocks"/>.
        /// </summary>
        public int LocalNoOfBlocks {
            get {
                Debug.Assert(m_Context.CellPartitioning.LocalLength == m_Context.iLogicalCells.NoOfLocalUpdatedCells);
                return m_Context.CellPartitioning.LocalLength;
            }
        }

        public int FirstBlock {
            get {
                return m_Context.CellPartitioning.i0;
            }
        }

        public int GetFirstBlock(int proc) {
            return m_Context.CellPartitioning.GetI0Offest(proc);
        }

        public int GetLocalNoOfBlocks(int proc) {
            return m_Context.CellPartitioning.GetLocalLength(proc);
        }

        public int FindProcessForBlock(int iBlk) {
            return m_Context.CellPartitioning.FindProcess(iBlk);
        }

        /// <summary>
        /// Always true-frame blocks are used here.
        /// </summary>
        public bool AllBlockSizesEqual {
            get {
                return true;
            }
        }

        public int GetBlockType(int iBlock) {
            this.GridDat.CellPartitioning.TestIfInLocalRange(iBlock);
            if(!CellDepLength) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // if the no. of DOFs per cell is constant, there is only one possible blocking,
                // i.e. the block type is 0
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                if(this.m_SubblkLen.Count == 0) {
                    // ensure that we have the blocking...
                    int NoOfBs = m_BasisS.Length;
                    int[] Len = new int[NoOfBs];
                    for (int iBs = 0; iBs < m_BasisS.Length; iBs++) {
                        Len[iBs] = m_BasisS[iBs].Length;
                    }

                    m_SubblkLen.Add(Len);
                }
                return 0;
            } else {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // if the no. of DOFs per cell is *not*, we have to find the blocking
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                int j = iBlock - this.GridDat.CellPartitioning.i0; // local cell index
                Basis[] Bs = m_BasisS;
                int NoOfBs = Bs.Length;
                int blockType;
                unsafe
                {
                    int* _Len_j = stackalloc int[NoOfBs];
                    for(int iBs = 0; iBs < NoOfBs; iBs++) {
                        _Len_j[iBs] = m_BasisS[iBs].GetLength(j);
                    }

                    bool bFound = false;
                    
                    int K = m_SubblkLen.Count;
                    for (blockType = 0; blockType < K; blockType++) { // loop over all block types found so far
                        int[] len_bT = m_SubblkLen[blockType];
                        Debug.Assert(len_bT.Length == NoOfBs);

                        bFound = true;
                        for (int i = 0; i < NoOfBs; i++) {
                            if (len_bT[i] != _Len_j[i]) {
                                bFound = false;
                                break;
                            }
                        }

                        if (bFound)
                            break;
                    }
                    if(bFound == false) {
                        Debug.Assert(blockType == m_SubblkLen.Count);
                        int[] __Len_j = new int[NoOfBs];
                        for(int iBs = 0; iBs < NoOfBs; iBs++) {
                            __Len_j[iBs] = _Len_j[iBs];
                        }
                        m_SubblkLen.Add(__Len_j);
                    }
                }
                return blockType;
            }

        }

        public int GetBlockI0(int iBlock) {
            this.GridDat.CellPartitioning.TestIfInLocalRange(iBlock);
            return this.m_MaxTotalNoOfCoordinatesPerCell * iBlock;
        }

        public int GetBlockLen(int iBlock) {
            return m_MaxTotalNoOfCoordinatesPerCell;
        }

        public int GetBlockIndex(int i) {
            return i / m_MaxTotalNoOfCoordinatesPerCell;
        }

        /// <summary>
        /// <see cref="IBlockPartitioning.GetImmutableBlockPartitioning"/>
        /// </summary>
        public IBlockPartitioning GetImmutableBlockPartitioning() {
            if(this.IsMutable) {
                return GetBlocking(true);
            } else {
                return this;
            }
        }

        /// <summary>
        /// maybe, it depends...
        /// </summary>
        override public bool IsMutable {
            get {
                return CellDepLength;
            }
        }
        
        /// <summary>
        /// True, if the degrees of freedom are cell-dependent.
        /// This depends on the used DG basis (see <see cref="BasisS"/>);
        /// </summary>
        public bool CellDepLength {
            get { return (m_MaxTotalNoOfCoordinatesPerCell != m_MinTotalNoOfCoordinatesPerCell) ; }
        }

        /// <summary>
        /// DG context
        /// </summary>
        protected internal IGridData m_Context;

		/// <summary>
		/// The grid on which this mapping relies on.
		/// </summary>
		public IGridData GridDat {
			get {
                return m_Context;
            }
		}


        /// <summary>
        /// local unique index of the 0th coordinate in the 0th cell for each field;
        /// index: field or basis index;
        /// </summary>
        internal protected int[] m_j0CoordinateIndex;


        /// <summary>
        /// 
        /// </summary>
        /// <param name="j">local cell index</param>
        /// <returns></returns>
        public int GetTotalNoOfCoordinatesPerCell(int j) {
            if (m_MaxTotalNoOfCoordinatesPerCell == m_MinTotalNoOfCoordinatesPerCell)
                // constant number of DG coordinates per cell
                return m_MaxTotalNoOfCoordinatesPerCell;


            // need to calculate the actual number of coordinates per cell
            int L = 0;
            for (int f = m_BasisS.Length - 1; f >= 0; f--)
                L += m_BasisS[f].GetLength(j);
            return L;
        }

        /// <summary>
        /// Local (on this MPI process) number of degrees-of-freedom; 
        /// this can be different than <see cref="Partitioning.LocalLength"/>, e.g. in the XDG case.
        /// </summary>
        public int GetLocalNoOfDOFs() {
            int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            int dofs = 0;
            for (int j = 0; j < J; j++) {
                dofs += GetTotalNoOfCoordinatesPerCell(j);
            }
            return dofs;
        }

        /// <summary>
        /// the sum of <see cref="GetLocalNoOfDOFs"/> over all MPI processors in communicator <see cref="MPI_Comm"/>.
        /// </summary>
        public int GetTotalNoOfDOFs() {
            return GetLocalNoOfDOFs().MPISum(this.MPI_Comm);
        }


        /// <summary>
        /// the total number of coordinates of all fields in this mapping in one cell.
        /// </summary>
        private int m_MaxTotalNoOfCoordinatesPerCell = 0;

        /// <summary>
        /// the minimum number of coordinates of all fields in this mapping in one cell.
        /// </summary>
        private int m_MinTotalNoOfCoordinatesPerCell = 0;

        /// <summary>
        /// the total number of coordinates of all fields in this mapping in one cell.
        /// </summary>
        public int MaxTotalNoOfCoordinatesPerCell {
            get { return m_MaxTotalNoOfCoordinatesPerCell; }
        }


        /// <summary>
        /// the total number of coordinates of all fields in this mapping in one cell.
        /// </summary>
        public int MinTotalNoOfCoordinatesPerCell {
            get { return m_MinTotalNoOfCoordinatesPerCell; }
        }

        /// <summary>
        /// Number of DG coordinates  of all fields in this mapping in one cell.
        /// </summary>
        public int NoOfCoordinatesPerCell {
            get {
                if (m_MinTotalNoOfCoordinatesPerCell != m_MaxTotalNoOfCoordinatesPerCell)
                    throw new ApplicationException();

                return m_MinTotalNoOfCoordinatesPerCell;
            }
        }

        /// <summary>
        /// Number of elements in <see cref="BasisS"/>, i.e. number of dependent variables resp. scalar fields.
        /// </summary>
        public int NoOfVariables {
            get {
                return m_BasisS.Length;
            }
        }


        /// <summary>
        /// <see cref="BasisS"/>;
        /// </summary>
        protected Basis[] m_BasisS;
               
        /// <summary>
        /// the basis object of all field in this mapping;
        /// </summary>
        /// <remarks>
        /// By this list, an <i>basis index</i>, or <i>field index</i> for each basis in this mapping is defined;
        /// </remarks>
        public IList<Basis> BasisS {
            get {
                //if(m_RdonlyBasisS == null) {
                //    m_RdonlyBasisS = m_BasisS.ToList().AsReadOnly();
                //}
                //return m_RdonlyBasisS;
                return m_BasisS;
            }
        }

        /// <summary>
        /// converts a global unique coordinate index into a local unique coordinate index
        /// </summary>
        /// <param name="iglobal"></param>
        /// <returns></returns>
        /// <see cref="Local2GlobalIndex"/>
        public int Global2LocalIndex(long iglobal) {
            long iLocal = iglobal - i0;

            if (iLocal >= 0 && iLocal < Ntotal) {
                return (int)iLocal;
            } else {
                throw new ApplicationException("unable to get local index from global index which does not belong to the local update range");
            }


        }

        /// <summary>
        /// total number of locally stored coordinates
        /// (internal, border and external elements);
        /// </summary>
        public int Ntotal {
            get {
                return (this.m_MaxTotalNoOfCoordinatesPerCell * m_Context.iLogicalCells.Count);
                //return m_Count;
            }
        }

        /// <summary>
        /// number of external, or ghost coordinates (DG coordinates in external (aka. ghost) cells);
        /// </summary>
        public int NExternal {
            get {
                return Ntotal - LocalLength;
            }
        }

        /// <summary>
        /// total number of coordinates in all processors
        /// </summary>
        public long GlobalCount {
            get {
                return m_Context.CellPartitioning.TotalLength*m_MaxTotalNoOfCoordinatesPerCell;
            }
        }

        /*
        public MPI_Comm MPI_Comm {
            get {
                return this.GridDat.;
            }
        }

        public int MpiSize {
            get {
                throw new NotImplementedException();
            }
        }

        public int MpiRank {
            get {
                throw new NotImplementedException();
            }
        }

        public int i0 {
            get {
                throw new NotImplementedException();
            }
        }

        public int iE {
            get {
                throw new NotImplementedException();
            }
        }

        public int LocalLength {
            get {
                throw new NotImplementedException();
            }
        }

        public int TotalLength {
            get {
                throw new NotImplementedException();
            }
        }

        public int GetI0Offest(int proc) {
            throw new NotImplementedException();
        }

        public IPartitioning GetImmutablePartition() {
            throw new NotImplementedException();
        }

        public int FindProcess(int index) {
            throw new NotImplementedException();
        }

        public int FindProcess(long index) {
            throw new NotImplementedException();
        }

        public bool IsInLocalRange(int i) {
            throw new NotImplementedException();
        }

        public int GetLocalLength(int proc) {
            throw new NotImplementedException();
        }
        */

        /// <summary>
        /// converts a local unique coordinate index into a global unique coordinate index
        /// </summary>
        /// <param name="iLocal"></param>
        /// <returns></returns>
        /// <see cref="Global2LocalIndex"/>
        public int Local2GlobalIndex(int iLocal) {
            if (iLocal < 0 || iLocal >= Ntotal)
                throw new IndexOutOfRangeException();

            if (iLocal < LocalLength) {
                return (i0 + iLocal);
            } else {
                int find;
                int j;
                int n;

                LocalFieldCoordinateIndex(iLocal, out find, out j, out n);

                j -= m_Context.iLogicalCells.NoOfLocalUpdatedCells;
                int jGlobal = (int)(m_Context.iParallel.GlobalIndicesExternalCells[j]);

                return jGlobal * m_MaxTotalNoOfCoordinatesPerCell + m_j0CoordinateIndex[find] + n;
            }
        }



        /// <summary>
        /// finds the owning processor (<paramref name="TargetProc"/>) and the local index on this processor 
        /// (<paramref name="iLocalTargetProc"/>) which corresponds to 
        /// a local index in the external range of this processor.
        /// </summary>
        /// <param name="iLocalExt">
        /// a local index on this processor; If this argument is not within the external range (greater or equal to
        /// <see cref="NUpdate"/> and smaller than <see cref="Ntotal"/>), this method still works
        /// </param>
        /// <param name="TargetProc">
        /// on exit, the process rank of the owner process of local external index <paramref name="iLocalExt"/>;
        /// </param>
        /// <param name="iLocalTargetProc">
        /// on exit, the local index on process <paramref name="TargetProc"/> which corresponds to
        /// <paramref name="iLocalExt"/>
        /// </param>
        public void TransformLocalExternal(int iLocalExt, out int TargetProc, out int iLocalTargetProc) {
            Debug.Assert(!(iLocalExt < 0 || iLocalExt >= Ntotal), "Index out of range.");

            if (iLocalExt < LocalLength) {
                iLocalTargetProc = iLocalExt;
                TargetProc = m_Context.CellPartitioning.MpiRank;
            } else {
                IGridData gd = m_Context;

                //Field f;
                int find;
                int j;
                int n;
                LocalFieldCoordinateIndex(iLocalExt, out find, out j, out n);

                j -= gd.iLogicalCells.NoOfLocalUpdatedCells;
                long jGlobal = gd.iParallel.GlobalIndicesExternalCells[j];
                TargetProc = gd.CellPartitioning.FindProcess(jGlobal);
                int jLocalTargetProc = (int)(jGlobal - gd.CellPartitioning.GetI0Offest(TargetProc));

                //int find = Array.IndexOf<Field>(m_Fields, f);

                iLocalTargetProc =  jLocalTargetProc * m_MaxTotalNoOfCoordinatesPerCell + m_j0CoordinateIndex[find] + n;
            }
        }

        /// <summary>
        /// computes a local unique coordinate index ("local" means local on this processor);
        /// this index is unique over all fields (in this mapping), over all cells, over all basis functions, 
        /// but it's only locally (on this processor) valid.
        /// A local index in the update range (smaller than <see cref="NUpdate"/>) can be converted into 
        /// a global index by adding <see cref="Partitioning.i0"/>.
        /// </summary>
        /// <param name="find">
        /// the field or basis index (see <see cref="BasisS"/>);
        /// </param>
        /// <param name="j">local cell index</param>
        /// <param name="n">basis index.</param>
        /// <returns>
        /// the local coordinate index for the given parameters.
        /// </returns>
        public int LocalUniqueCoordinateIndex(int find, int j, int n) {
            Debug.Assert(!(find < 0 || find > this.m_BasisS.Length),
                "Field is not in this mapping");
            Debug.Assert(!(j < 0 || j >= m_Context.iLogicalCells.Count),
                "j must be greater or equal 0 and less than number of locally updated cells");
            //Debug.Assert(!(n < 0 || n >= m_BasisS[find].GetLength(j)),
            Debug.Assert(!(n < 0 || n >= m_BasisS[find].MaximalLength),
                "n must be greater or equal 0 and less than the basis index determined by the specified interpolation order");

            return j * m_MaxTotalNoOfCoordinatesPerCell + m_j0CoordinateIndex[find] + n;
        }
        
        /// <summary>
        /// computes a global unique coordinate index ("global" means over all MPI processors)
        /// from a local field coordinate index (i.e. the tuple of field index 
        /// <paramref name="find"/>, local cell index <paramref name="j"/> and 
        /// DG coordinate index <paramref name="n"/>).
        /// </summary>
        /// <param name="find">
        /// the field or basis index (see <see cref="BasisS"/>);
        /// </param>
        /// <param name="j">local cell index</param>
        /// <param name="n">basis index.</param>
        /// <returns>
        /// the global DG coordinate index for the given parameters.
        /// </returns>
        /// <remarks>
        /// The global unique coordinate index is unique over all fields (in this mapping), 
        /// over all cells, over all basis functions, over all MPI processes in the current
        /// communicator.
        /// </remarks>
        public int GlobalUniqueCoordinateIndex(int find, int j, int n) {
            int iloc = LocalUniqueCoordinateIndex(find, j, n);
            return Local2GlobalIndex(iloc);
        }


        /// <summary>
        /// computes a global unique coordinate index ("global" means over all MPI processors)
        /// from a global field coordinate index (i.e. the tuple of field index 
        /// <paramref name="find"/>, global cell index <paramref name="jGlobal"/> and 
        /// DG coordinate index <paramref name="n"/>).
        /// this index is unique over all fields (in this mapping), over all cells, over all basis functions;
        /// </summary>
        /// <param name="find">
        /// the field or basis index (see <see cref="BasisS"/>);
        /// </param>
        /// <param name="jGlobal">global cell index</param>
        /// <param name="n">basis index.</param>
        /// <returns>
        /// the global DG coordinate index for the given parameters.
        /// </returns>
        /// <remarks>
        /// The global unique coordinate index is unique over all fields (in this mapping), 
        /// over all cells, over all basis functions, over all MPI processes in the current
        /// communicator.
        /// </remarks>
        public long GlobalUniqueCoordinateIndex_FromGlobal(int find, long jGlobal, int n) {
            if (find < 0 || find > m_BasisS.Length)
                throw new ArgumentOutOfRangeException("find", "field index");
            if (jGlobal < 0 || jGlobal >= GridDat.CellPartitioning.TotalLength)
                throw new ArgumentOutOfRangeException("j must be greater or equal 0 and less than the total number of cells over all processors");
            if (n < 0 || n >= m_BasisS[find].MaximalLength)
                throw new ArgumentOutOfRangeException("n must be greater or equal 0 and less than the basis index determined by the specified interpolation order");

            return jGlobal * m_MaxTotalNoOfCoordinatesPerCell + m_j0CoordinateIndex[find] + n;
         }

        
        /// <summary>
        /// inverse mapping of <see cref="LocalUniqueCoordinateIndex"/>;
        /// </summary>
        /// <param name="Index">local unique coordinate index</param>
        /// <param name="ifld">on exit, the field or basis index (see
        /// <see cref="BasisS"/>) of the field
        /// to which <paramref name="Index"/> belongs.</param>
        /// <param name="j">on exit, the local cell index</param>
        /// <param name="n">on exit, the basis function index</param>
        /// <returns></returns>
        public void LocalFieldCoordinateIndex(int Index, out int ifld, out int j, out int n) {
            Debug.Assert(!(Index < 0 || Index >= Ntotal), "Index out of range");


            j = Index / m_MaxTotalNoOfCoordinatesPerCell;
            int jrst = Index % m_MaxTotalNoOfCoordinatesPerCell;

            int Lf = m_BasisS.Length;
            ifld = 0;
            while (ifld < Lf) {
                if (jrst < m_j0CoordinateIndex[ifld + 1])
                    break;
                ifld++;
            }
            n = jrst - m_j0CoordinateIndex[ifld];
        }

        /// <summary>
        /// two <see cref="UnsetteledCoordinateMapping"/>, if their <see cref="BasisS"/>-properties
        /// are equal (for every entry);
        /// </summary>
        public override bool Equals(object obj) {
            UnsetteledCoordinateMapping othr = obj as UnsetteledCoordinateMapping;
            if (obj == null)
                return false;
            else 
                return EqualsUnsetteled(othr);
        }

        /// <summary>
        /// two <see cref="UnsetteledCoordinateMapping"/>, if their <see cref="BasisS"/>-properties
        /// are equal (for every entry);
        /// </summary>
        /// <remarks>
        /// This is a non-virtual version of <see cref="Equals"/>;
        /// See <see cref="CoordinateMapping.Equals"/> for explanation.
        /// </remarks>
        public bool EqualsUnsetteled(UnsetteledCoordinateMapping obj) {
            UnsetteledCoordinateMapping othr = obj;
            
            // the following calls may be not necessary, but since they are very cheap ...
            if (!object.ReferenceEquals(this.m_Context, othr.m_Context))
                return false;


            // the real test...
            if (this.m_BasisS.Length != othr.m_BasisS.Length)
                return false;
            for (int i = 0; i < this.m_BasisS.Length; i++)
                if (!this.m_BasisS[i].Equals(othr.m_BasisS[i]))
                    return false;
            
            // success
            return true;
        }
        


        /// <summary>
        /// calls base implementation
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();            
        }

        /// <summary>
        /// Returns global unique indices which correlate to a certain sub-set of this mappings's basises (<see cref="BasisS"/>).
        /// </summary>
        /// <param name="Fields">
        /// indices into <see cref="BasisS"/>
        /// </param>
        /// <param name="lo">
        /// Inclusion of unused indices:
        ///  - true: all indices up to <see cref="Basis.MaximalLength"/>;
        ///  - false: only occupied entries, up to <see cref="Basis.GetLength"/>
        /// </param>
        /// <returns>a list of global (over all MPI processes) unique indices.</returns>
        public int[] GetSubvectorIndices(bool lo, params int[] Fields) {

            // <param name="includeExternal">
            // true, if indices which correlate to external cells should be included; false if not.
            // </param>
            bool includeExternal = false;

            int L = 0;
            for (int i = 0; i < Fields.Length; i++) {
                int iField = Fields[i];
                if (iField < 0 || iField > this.m_BasisS.Length)
                    throw new IndexOutOfRangeException();

                L += m_BasisS[iField].MaximalLength;

            }
            
            int J;
            if (includeExternal)
                J = this.m_Context.iLogicalCells.Count;
            else
                J = this.m_Context.iLogicalCells.NoOfLocalUpdatedCells;
            List<int> R = new List<int>(L*J);

            for (int j = 0; j < J; j++) {

                for (int i = 0; i < Fields.Length; i++) {
                    int iField = Fields[i];
                    int N;
                    if (lo)
                        N = m_BasisS[iField].MaximalLength;
                    else
                        N = m_BasisS[iField].GetLength(j); 

                    int i0 = (int) this.GlobalUniqueCoordinateIndex(iField, j, 0);


                    for (int n = 0; n < N; n++) {
                        R.Add(i0 + n);
                    }
                }
            }

            return R.ToArray();
        }

        /// <summary>
        ///  returns global unique indices which correlate to a certain sub-set of this mappings's basis (<see cref="BasisS"/>),
        ///  and to cells in a certain sub-grid (<paramref name="sgrd"/>).
        /// </summary>
        /// <param name="sgrd"></param>
        /// <param name="Fields">
        /// indices into <see cref="BasisS"/>
        /// </param>
        /// <param name="lo">
        /// Inclusion of unused indices:
        ///  - true: all indices up to <see cref="Basis.MaximalLength"/>
        ///  - false: only occupied entries, up to <see cref="Basis.GetLength"/>
        /// </param>
        /// <returns>a list of global (over all MPI processes) unique indices.</returns>
        public int[] GetSubvectorIndices(SubGrid sgrd, bool lo, int[] Fields) {

            // <param name="includeExternal">
            // true, if indices which correlate to external cells (i.e. ghost cells) should be included; false if not.
            // </param>
            bool includeExternal = false;

            var _BasisS = this.BasisS.ToArray();
            
            int L = 0;
            for (int i = 0; i < Fields.Length; i++) {
                int iField = Fields[i];
                if (iField < 0 || iField > _BasisS.Length)
                    throw new IndexOutOfRangeException();

                L += _BasisS[iField].MaximalLength;
            }

            var jSub2Full = sgrd.SubgridIndex2LocalCellIndex;
            var R = new List<int>();

            int JSUB;
            if(includeExternal)
                JSUB = sgrd.LocalNoOfCells_WithExternal;
            else
                JSUB = sgrd.LocalNoOfCells;
            for (int jsub = 0; jsub < JSUB; jsub++) {
                int jCell = jSub2Full[jsub];
                
                for (int i = 0; i < Fields.Length; i++) {
                    int iField = Fields[i];

                    int N;
                    if (lo)
                        N = m_BasisS[iField].MaximalLength;
                    else
                        N = m_BasisS[iField].GetLength(jCell);

                    int i0 = (int)this.GlobalUniqueCoordinateIndex(iField, jCell, 0);


                    for (int n = 0; n < N; n++) {
                        R.Add(i0 + n);
                    }
                }
            }

            return R.ToArray();
        }

        
    }
}

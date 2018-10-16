using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// Mapping of an old to a new grid under adaptive mesh refinement (see <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]})"/>.
    /// </summary>
    public class GridCorrelation {

        /// <summary>
        /// Mapping from old cell index to old GlobalId.
        /// - index: local cell index into old grid
        /// - content: GlobalId of respective cell
        /// </summary>
        public long[] OldGlobalId;

        /// <summary>
        /// Mapping form old cell index to GlobalId of new cells.
        /// - 1st index: local cell index into old grid
        /// - 2nd index: enumeration; for cells which are refined, this maps to the subdivisions; for coarsened cells, or cells which are not changed, the enumeration contains only one element.
        /// - content: GlobalId of destination cell
        /// </summary>
        public long[][] DestGlobalId;

        /// <summary>
        /// Affine transformation between old and new cells.
        /// - indexing: correlates with <see cref="DestGlobalId"/>
        /// - content: index into <see cref="GeometricMapping"/>; if negative, the mapping is the identity.
        /// </summary>
        public int[][] MappingIndex;

        /// <summary>
        /// Affine transformation from cell-local coordinate system to cell-local coordinate system.
        /// - 1st index: cell reference element, correlates with <see cref="GridCommons.RefElements"/>
        /// - 2nd index: refinement leaf index
        /// </summary>
        public AffineTrafo[][] GeometricMapping;

        /// <summary>
        /// Subdivision nodes used for refinement.
        /// - 1st index: cell reference element, correlates with <see cref="GridCommons.RefElements"/>
        /// - 2nd index: refinement leaf index
        /// </summary>
        public RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves;

        /// <summary>
        /// Cache for <see cref="GetSubdivBasisTransform(int, int, int)"/>.
        /// - 1st index: reference element
        /// - 2nd index: subdivision node
        /// - Item 1: basis transformation matrix, <see cref="RefElement.SubdivisionTreeNode.GetBasisTrafo(int)"/>
        /// - Item 2: Jacobian determinate, <see cref="RefElement.SubdivisionTreeNode.Trafo2Root"/>
        /// </summary>
        Tuple<MultidimensionalArray,double>[][] Subdiv_BasisTransform;
        
        /// <summary>
        /// Transformation matrices, for cell subdivision, as computed by <see cref="RefElement.SubdivisionTreeNode.GetBasisTrafo(int)"/>,
        /// are cached by this method.
        /// </summary>
        /// <param name="iKref">Reference element index, correlates with 1st index of <see cref="KrefS_SubdivLeaves"/>.</param>
        /// <param name="iSubdiv">Subdivision leaf index, correlates with 2nd index of <see cref="KrefS_SubdivLeaves"/>.</param>
        /// <param name="p">Requested polynomial degree.</param>
        /// <returns></returns>
        public Tuple<MultidimensionalArray,double> GetSubdivBasisTransform(int iKref, int iSubdiv, int p) {
            if(Subdiv_BasisTransform == null) {
                Subdiv_BasisTransform = new Tuple<MultidimensionalArray,double>[KrefS_SubdivLeaves.Length][];
            }

            if(Subdiv_BasisTransform[iKref] == null) {
                Subdiv_BasisTransform[iKref] = new Tuple<MultidimensionalArray,double>[KrefS_SubdivLeaves[iKref].Length];
            }

            RefElement Kref = KrefS_SubdivLeaves[iKref][0].RefElement;
            int Np = Kref.GetNoOfOrthonormalPolynomialsUptoDegree(p);

            if(Subdiv_BasisTransform[iKref][iSubdiv] == null || Subdiv_BasisTransform[iKref][iSubdiv].Item1.NoOfCols < Np) {
                MultidimensionalArray R = KrefS_SubdivLeaves[iKref][iSubdiv].GetBasisTrafo(p);
                Subdiv_BasisTransform[iKref][iSubdiv] = new Tuple<MultidimensionalArray, double>(
                    R,
                    KrefS_SubdivLeaves[iKref][iSubdiv].TrafoFromRoot.Matrix.Determinant()
                    );
            }

            if(Subdiv_BasisTransform[iKref][iSubdiv].Item1.NoOfCols >= Np) {
                return new Tuple<MultidimensionalArray, double>(
                    Subdiv_BasisTransform[iKref][iSubdiv].Item1.ExtractSubArrayShallow(new[] { 0, 0 }, new[] { Np - 1, Np - 1 }),
                    Subdiv_BasisTransform[iKref][iSubdiv].Item2
                    );
            } else {
                Debug.Assert(Subdiv_BasisTransform[iKref][iSubdiv].Item1.NoOfCols == Np);
                Debug.Assert(Subdiv_BasisTransform[iKref][iSubdiv].Item1.NoOfRows == Np);
                return Subdiv_BasisTransform[iKref][iSubdiv];
            }

        }

        /// <summary>
        /// Mapping from local cell indices of the old grid to global cell indices of the new grid.
        /// - 1st index: local cell index into old grid
        /// - 2nd index: enumeration; for cells which are refined, this maps to the subdivisions; for coarsened cells, or cells which are not changed, the enumeration contains only one element.
        /// - content: global target index (index into new grid) of respective cell.
        /// </summary>
        int[][] TargetIdx;
       
        public void ComputeDataRedist(GridData NewGrid) {
            using(new FuncTrace()) {
                Debug.Assert(DestGlobalId.Length == MappingIndex.Length);
                Debug.Assert(OldGlobalId.Length == MappingIndex.Length);
                int oldJ = DestGlobalId.Length;

                Permutation invSigma = NewGrid.CurrentGlobalIdPermutation.Invert();

                List<long> _DestGlobalId = new List<long>();
                for(int j = 0; j < oldJ; j++) {
                    _DestGlobalId.AddRange(DestGlobalId[j]);
                }

                long[] _TargetIdx = new long[_DestGlobalId.Count];
                invSigma.EvaluatePermutation(_DestGlobalId, _TargetIdx);

                TargetIdx = new int[oldJ][];
                int c2 = 0;
                for(int j = 0; j < oldJ; j++) {
                    int K = DestGlobalId[j].Length;
                    TargetIdx[j] = new int[K];
                    for(int k = 0; k < K; k++) {
                        TargetIdx[j][k] = (int)(_TargetIdx[c2]);
                        c2++;
                    }
                }
                Debug.Assert(c2 == _DestGlobalId.Count);
            }
        }


        public void ApplyToVector<I>(IList<I> input, IList<I[]> output, IPartitioning outputPartitioning) {
            using(new FuncTrace()) {
                Debug.Assert(DestGlobalId.Length == MappingIndex.Length);
                Debug.Assert(OldGlobalId.Length == MappingIndex.Length);
                int oldJ = DestGlobalId.Length;

                if(input.Count != oldJ)
                    throw new ArgumentException("Mismatch between input vector length and current data length.");
                if(output.Count != outputPartitioning.LocalLength)
                    throw new ArgumentException("Length mismatch of output list and output partition.");

                int j0Dest = outputPartitioning.i0;

                // keys: processors which should receive data from this processor
                Dictionary<int, ApplyToVector_Helper<I>> AllSendData = new Dictionary<int, ApplyToVector_Helper<I>>();

                for(int j = 0; j < oldJ; j++) {
                    I data_j = input[j];

                    foreach(int jDest in TargetIdx[j]) {
                        if(outputPartitioning.IsInLocalRange(jDest)) {
                            I[] destCollection = output[jDest - j0Dest];
                            ArrayTools.AddToArray(data_j, ref destCollection);
                            output[jDest - j0Dest] = destCollection;
                        } else {
                            int targProc = outputPartitioning.FindProcess(jDest);

                            ApplyToVector_Helper<I> dataTargPrc;
                            if(!AllSendData.TryGetValue(targProc, out dataTargPrc)) {
                                dataTargPrc = new ApplyToVector_Helper<I>();
                                AllSendData.Add(targProc, dataTargPrc);
                            }

                            dataTargPrc.TargetIndices.Add(jDest);
                            dataTargPrc.Items.Add(data_j);


                        }

                    }

                }

                var AllRcvData = SerialisationMessenger.ExchangeData(AllSendData, outputPartitioning.MPI_Comm);

                foreach(var kv in AllRcvData) {
                    int rcvProc = kv.Key;
                    j0Dest = outputPartitioning.GetI0Offest(rcvProc);

                    var TIdxs = kv.Value.TargetIndices;
                    var TVals = kv.Value.Items;
                    Debug.Assert(TIdxs.Count == TVals.Count);
                    int L = TIdxs.Count;

                    for(int l = 0; l < L; l++) {
                        int idx = TIdxs[l] - j0Dest;
                        Debug.Assert(outputPartitioning.IsInLocalRange(idx));

                        I[] destCollection = output[idx];
                        ArrayTools.AddToArray(TVals[idx], ref destCollection);
                        output[idx] = destCollection;
                    }


                }

            }
        }

        /// <summary>
        /// Data structure used in <see cref="ApplyToVector{I}(IList{I}, IList{I[]}, IPartitioning)"/>.
        /// </summary>
        [Serializable]
        class ApplyToVector_Helper<I> {
            public List<int> TargetIndices = new List<int>();
            public List<I> Items = new List<I>();
        }


        int[][] m_TargetMappingIndex;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="TargetMappingIndex"></param>
        /// <param name="outputPartitioning">
        /// Partitioning of the new grid, resp the return array.
        /// </param>
        /// <returns>
        /// - 1st index: cell index in new grid, correlates with <paramref name="outputPartitioning"/>.
        /// - 2nd index: enumeration over cells (in the old grid) which are combined in the new grid. 
        ///   For cells with refinement, always one entry, for cells which are coarsened a greater number of entries.
        ///   If null, the cell is not changed.
        /// - content: Subdivision leaf index, correlates with 2nd index of <see cref="KrefS_SubdivLeaves"/>, can be used as an input to <see cref="GetSubdivBasisTransform(int, int, int)"/>.
        /// </returns>
        public int[][] GetTargetMappingIndex(IPartitioning outputPartitioning) {
            using(new FuncTrace()) {
                Debug.Assert(DestGlobalId.Length == MappingIndex.Length);
                Debug.Assert(OldGlobalId.Length == MappingIndex.Length);
                int oldJ = DestGlobalId.Length;

                // Caching
                // =======

                if(m_TargetMappingIndex != null) {
                    // caching 
                    if(m_TargetMappingIndex.Length != outputPartitioning.LocalLength)
                        throw new ArgumentException("Length mismatch of output list and output partition.");

                    return m_TargetMappingIndex;
                }

                // local evaluation, prepare communication
                // =======================================

                m_TargetMappingIndex = new int[outputPartitioning.LocalLength][];

                int j0Dest = outputPartitioning.i0;

                // keys: processors which should receive data from this processor
                Dictionary<int, GetTargetMapping_Helper> AllSendData = new Dictionary<int, GetTargetMapping_Helper>();

                for(int j = 0; j < oldJ; j++) {
                    int[] MappingIndex_j = MappingIndex[j];
                    if(MappingIndex_j != null) {
                        Debug.Assert(TargetIdx[j].Length == MappingIndex_j.Length);
                        int L = MappingIndex_j.Length;

                        for(int l = 0; l < L; l++) {
                            
                            int jDest = TargetIdx[j][l];
                            int MapIdx = MappingIndex_j[l];

                            if(outputPartitioning.IsInLocalRange(jDest)) {
                                int[] destCollection = m_TargetMappingIndex[jDest - j0Dest];
                                ArrayTools.AddToArray(MapIdx, ref destCollection);
                                m_TargetMappingIndex[jDest - j0Dest] = destCollection;
                            } else {
                                int targProc = outputPartitioning.FindProcess(jDest);

                                GetTargetMapping_Helper dataTargPrc;
                                if(!AllSendData.TryGetValue(targProc, out dataTargPrc)) {
                                    dataTargPrc = new GetTargetMapping_Helper();
                                    AllSendData.Add(targProc, dataTargPrc);
                                }

                                dataTargPrc.TargetIndices.Add(jDest);
                                dataTargPrc.Items.Add(MapIdx);
                            }

                        }
                    } else {
                        Debug.Assert(TargetIdx[j].Length == 1);
                    }
                }

                // communication
                // =============

                var AllRcvData = SerialisationMessenger.ExchangeData(AllSendData, outputPartitioning.MPI_Comm);

                foreach(var kv in AllRcvData) {
                    int rcvProc = kv.Key;
                    j0Dest = outputPartitioning.GetI0Offest(rcvProc);

                    var TIdxs = kv.Value.TargetIndices;
                    var TVals = kv.Value.Items;
                    Debug.Assert(TIdxs.Count == TVals.Count);
                    int L = TIdxs.Count;

                    for(int l = 0; l < L; l++) {
                        int idx = TIdxs[l] - j0Dest;
                        Debug.Assert(outputPartitioning.IsInLocalRange(idx));

                        int[] destCollection = m_TargetMappingIndex[idx];
                        ArrayTools.AddToArray(TVals[idx], ref destCollection);
                        m_TargetMappingIndex[idx] = destCollection;
                    }
                }

                // return
                // ======

                return m_TargetMappingIndex;
            }
        }
               


        /// <summary>
        /// Data structure used in <see cref="GetTargetMappingIndex(int[][], IPartitioning)"/>.
        /// </summary>
        [Serializable]
        class GetTargetMapping_Helper {
            public List<int> TargetIndices = new List<int>();
            public List<int> Items = new List<int>();
        }

    }
}

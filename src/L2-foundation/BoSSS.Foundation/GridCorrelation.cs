using BoSSS.Foundation.Comm;
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
        /// - 2nd index: enumeration
        /// - content: 
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
        /// </summary>
        public AffineTrafo[] GeometricMapping;



        int[][] TargetIdx;

        public void ComputeDataRedist(GridData NewGrid) {
            using(new FuncTrace()) {
                Debug.Assert(DestGlobalId.Length == MappingIndex.Length);
                Debug.Assert(OldGlobalId.Length == MappingIndex.Length);
                int J = DestGlobalId.Length;

                Permutation invSigma = NewGrid.CurrentGlobalIdPermutation.Invert();

                List<long> _DestGlobalId = new List<long>();
                for(int j = 0; j < J; j++) {
                    _DestGlobalId.AddRange(DestGlobalId[j]);
                }

                long[] _TargetIdx = new long[_DestGlobalId.Count];
                invSigma.EvaluatePermutation(_DestGlobalId, _TargetIdx);

                TargetIdx = new int[J][];
                int c2 = 0;
                for(int j = 0; j < J; j++) {
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
                int J = DestGlobalId.Length;

                if(input.Count != J)
                    throw new ArgumentException("Mismatch between input vector length and current data length.");
                if(output.Count != outputPartitioning.LocalLength)
                    throw new ArgumentException("Length mismatch of output list and output partition.");

                int j0Dest = outputPartitioning.i0;

                // keys: processors which should receive data from this processor
                Dictionary<int, ApplyToVector_Helper<I>> AllSendData = new Dictionary<int, ApplyToVector_Helper<I>>();

                for(int j = 0; j < J; j++) {
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

        [Serializable]
        class ApplyToVector_Helper<I> {
            public List<int> TargetIndices = new List<int>();
            public List<I> Items = new List<I>();
        }


    }
}

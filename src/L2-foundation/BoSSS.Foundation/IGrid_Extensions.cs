using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid {
    
    /// <summary>
    /// Various extension methods.
    /// </summary>
    static public class IGrid_Extensions {

        

        /// <summary>
        /// sets values for <see cref="Cell.CellFaceTags"/> by using a
        /// <paramref name="EdgeTagFunc"/>-function; also adds entries with empty names
        /// to the <see cref="GridCommons.EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        static public void DefineEdgeTags(this IGrid g, Func<Vector, byte> EdgeTagFunc, params string[] EdgeTagNamesToEnsure) {
            g.DefineEdgeTags(delegate (Vector v) {
                string EdgeTagName = null;
                byte et = EdgeTagFunc(v);
                if (!g.EdgeTagNames.ContainsKey(et))
                    throw new ArgumentException("unable to find EdgeTagName for EdgeTag = " + et);

                if (et >= Classic.GridCommons.FIRST_PERIODIC_BC_TAG)
                    throw new ApplicationException("edge tags greater or equal to " + Classic.GridCommons.FIRST_PERIODIC_BC_TAG + " are reserved for periodic \"boundaries\"");
                EdgeTagName = g.EdgeTagNames[et];
                return EdgeTagName;
            });
        }

        
        /// <summary>
        /// sets values for <see cref="Cell.CellFaceTags"/> by using a
        /// <paramref name="EdgeTagFunc"/>-function; also adds entries with empty names
        /// to the <see cref="IGridData.EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        static public void DefineEdgeTags(this IGrid g, Func<double[], byte> EdgeTagFunc) {
            g.DefineEdgeTags((Vector v) => EdgeTagFunc(v));
        }

        /// <summary>
        /// sets values for <see cref="Cell.CellFaceTags"/> by using a
        /// <paramref name="EdgeTagFunc"/>-function; also adds entries with empty names
        /// to the <see cref="IGridData.EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        static public void DefineEdgeTags(this IGrid g, Func<double[], string> EdgeTagFunc) {
            g.DefineEdgeTags((Vector v) => EdgeTagFunc(v));
        }

        /// <summary>
        /// Adds the edge tag name <paramref name="EdgeTagName"/> to the grid, 
        /// regardless of whether it is used or not.
        /// </summary>
        static public byte AddEdgeTag(this IGrid g, string EdgeTagName) {

            bool[] UsedEt = new bool[byte.MaxValue];
            foreach (var kv in g.EdgeTagNames) {
                UsedEt[kv.Key] = true;
                if (kv.Value.Equals(EdgeTagName))
                    return kv.Key;

            }

            byte NewEt = byte.MaxValue;
            for(int i = 1; i < UsedEt.Length; i++) {
                if (i >= GridCommons.FIRST_PERIODIC_BC_TAG)
                    throw new ArgumentOutOfRangeException("already all edge tags used");

                if (UsedEt[i] == false) {
                    NewEt = (byte)i;
                    break;
                }
            }
            Debug.Assert(g.EdgeTagNames.Keys.Contains(NewEt) == false);
            g.EdgeTagNames.Add(NewEt, EdgeTagName);

            return NewEt;
        }


        /// <summary>
        /// sets values for <see cref="Cell.CellFaceTags"/> by using a
        /// <paramref name="EdgeTagFunc"/>-function; also adds entries with empty names
        /// to the <see cref="IGridData.EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        static public void DefineEdgeTags(this IGrid g, Func<Vector, string> EdgeTagFunc, params string[] EdgeTagNamesToEnsure) {

            int D = g.SpatialDimension;
            double[] x = new double[D];
            MultidimensionalArray GlobalVerticesOut = MultidimensionalArray.CreateWrapper(x, 1, D);

            Dictionary<string, byte> EdgeTagNames_Reverse = new Dictionary<string, byte>();
            EdgeTagNames_Reverse.Add("inner edge", 0);
            int[] etUseCount = new int[byte.MaxValue];
            foreach(var kv in g.EdgeTagNames) {
                if(kv.Key == 0) {
                    
                } else {
                    EdgeTagNames_Reverse.Add(kv.Value, kv.Key);

                    if (kv.Key >= GridCommons.FIRST_PERIODIC_BC_TAG)
                        etUseCount[kv.Key] = 1;
                }
            }



            GridCommons baseGrid = null;
            if(g is GridCommons gcm) {
                baseGrid = gcm;
            } else if(g is Aggregation.AggregationGrid ag) {
                baseGrid = ag.RootGrid as GridCommons;
            }
            var GrdDatTmp = g.iGridData;
            
            // loop over edges...
            bool LoopOverEdges(Func<int,int,int,byte> GetEt) {
                bool _GridChanged = false;
                int NoOfEdges = GrdDatTmp.iGeomEdges.Count;
                for (int iEdge = 0; iEdge < NoOfEdges; ++iEdge) {
                    if (GrdDatTmp.iGeomEdges.IsEdgeBoundaryEdge(iEdge)) {

                        int jCell = GrdDatTmp.iGeomEdges.CellIndices[iEdge, 0];
                        int iFace = GrdDatTmp.iGeomEdges.FaceIndices[iEdge, 0];

                        // get edge tag
                        byte et = GetEt(iEdge, jCell, iFace);

                        // record edge tag
                        if (baseGrid != null) {
                            var _Cell = baseGrid.Cells[jCell];
                            var allCFTs = _Cell.CellFaceTags;
                            int found = 0;
                            if (allCFTs != null) {
                                for (int i = 0; i < allCFTs.Length; i++) {
                                    if (allCFTs[i].NeighCell_GlobalID < 0 && allCFTs[i].FaceIndex == iFace) {
                                        found++;
                                        if (allCFTs[i].EdgeTag == et) {
                                            // nop
                                        } else {
                                            allCFTs[i].EdgeTag = et;
                                            _GridChanged = true;
                                        }
                                    }
                                }
                            }

                            if (found > 1)
                                throw new ApplicationException(string.Format("Cell face tags inconsistent in cell (GlId={0},LocIdx={1}): found {2} boundary tags for face {3}.", _Cell.GlobalID, jCell, found, iFace));
                            if (found <= 0) {
                                CellFaceTag CFT = new CellFaceTag() {
                                    EdgeTag = et,
                                    FaceIndex = iFace,
                                    NeighCell_GlobalID = long.MinValue
                                };

                                CFT.AddToArray(ref baseGrid.Cells[jCell].CellFaceTags);
                                _GridChanged = true;
                            }
                        }
                        GrdDatTmp.iGeomEdges.EdgeTags[iEdge] = et;

                        etUseCount[et]++;
                    }
                }

                return _GridChanged;
            }

            byte RecordTag(int iedge, int jCell, int iFace) {
                var KRef = GrdDatTmp.iGeomCells.GetRefElement(jCell);

                // call edge-tag-name function
                GrdDatTmp.TransformLocal2Global(KRef.GetFaceCenter(iFace), GlobalVerticesOut, jCell);
                string EdgeTagName = EdgeTagFunc(x);

                // obtain edge tag
                if (!EdgeTagNames_Reverse.ContainsKey(EdgeTagName)) {
                    int NewTag = EdgeTagNames_Reverse.Count;
                    Debug.Assert(NewTag > 0);
                    if (NewTag >= GridCommons.FIRST_PERIODIC_BC_TAG)
                        throw new ApplicationException("To many different edge tag names; at maximum " + (Classic.GridCommons.FIRST_PERIODIC_BC_TAG - 1) + " different tags for non-periodic boundaries.");
                    EdgeTagNames_Reverse.Add(EdgeTagName, (byte)NewTag);
                }
                var et = EdgeTagNames_Reverse[EdgeTagName];
                return et;
            }


            // pass 1: assign edge tags locally
            bool GridChanged = LoopOverEdges(RecordTag);
            GridChanged = GridChanged.MPIOr();

            // pass 2: MPI syncronization
            if(GridChanged && g.Size > 1) {
                byte[] ETTranslation = SyncEdgeTagsOverMPI(EdgeTagNames_Reverse);
                if(ETTranslation != null) {
                    LoopOverEdges(delegate (int iedge, int jCell, int iFace) {
                        byte oldEt = GrdDatTmp.iGeomEdges.EdgeTags[iedge];
                        byte newEt = ETTranslation[oldEt];
                        return newEt;
                    });
                }
            }


            // store & return
            etUseCount = etUseCount.MPISum();
            g.EdgeTagNames.Clear();
            foreach(var kv in EdgeTagNames_Reverse) {
                if(kv.Value == 0 || etUseCount[kv.Value] > 0)
                    g.EdgeTagNames.Add(kv.Value, kv.Key);
            }

            if (GridChanged) {
                g.InvalidateGridData();
                Console.WriteLine("Grid Edge Tags changed.");
            }

            if (EdgeTagNamesToEnsure != null) {
                foreach (string EdgeTagName in EdgeTagNamesToEnsure) {
                    g.AddEdgeTag(EdgeTagName);
                }

                {
                    string[] allBndys = g.EdgeTagNames.Values.ToArray();
                    string[][] allBndys_allProcs = allBndys.MPIAllGatherO();

                    for (int r = 0; r < allBndys_allProcs.Length; r++) {
                        if (!allBndys_allProcs[r].SetEquals(allBndys)) {
                            throw new ApplicationException("Internal Error: mismatch in edge tag names among MPI processors.");
                        }
                    }
                }
            }
        }

        private static byte[] SyncEdgeTagsOverMPI(Dictionary<string, byte> EdgeTagNames_Reverse) {

            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MyRank);

            var To0 = new Dictionary<int, KeyValuePair<string, byte>[]>();
            if (MyRank > 0) {
                To0.Add(0, EdgeTagNames_Reverse.ToArray());
            }

            var allData = SerialisationMessenger.ExchangeData(To0);

            bool[] usedEdgeTags = new bool[byte.MaxValue + 1];
            foreach (var et in EdgeTagNames_Reverse.Values)
                usedEdgeTags[et] = true;
            byte GetNewEt() {
                for (int i = 1; i < usedEdgeTags.Length; i++) {
                    if (i >= GridCommons.FIRST_PERIODIC_BC_TAG)
                        throw new ApplicationException("Running out of edge tags.");
                    if (usedEdgeTags[i] == false) {
                        usedEdgeTags[i] = true;
                        return (byte)i;
                    }
                }
                throw new ApplicationException("Running out of edge tags.");
            }


            // collect everything on rank 0, sync and broadcast:

            if (MyRank == 0) {

                foreach (var kv in allData) {
                    var backData = kv.Value;

                    for (int i = 0; i < backData.Length; i++) {
                        if (EdgeTagNames_Reverse.ContainsKey(backData[i].Key)) {

                        } else {
                            byte sugKey = backData[i].Value;
                            if (usedEdgeTags[sugKey]) {
                                sugKey = GetNewEt();
                            }

                            EdgeTagNames_Reverse.Add(backData[i].Key, sugKey);
                        }
                    }
                }
            }

            var AllEts = EdgeTagNames_Reverse.ToArray().MPIBroadcast(0);

            byte[] EtTanslations = (byte.MaxValue + 1).ForLoop(i => (byte)i);

            bool AnyTranslation = false;
            if (MyRank > 0) {
                foreach (var kv in AllEts) {
                    if (EdgeTagNames_Reverse.ContainsKey(kv.Key)) {
                        byte oldVal = EdgeTagNames_Reverse[kv.Key];
                        byte newVal = kv.Value;
                        AnyTranslation = AnyTranslation | (oldVal != newVal);
                        EdgeTagNames_Reverse[kv.Key] = newVal;
                        EtTanslations[oldVal] = newVal;
                    } else {
                        EdgeTagNames_Reverse.Add(kv.Key, kv.Value);
                    }
                }
            }
            AnyTranslation = AnyTranslation.MPIOr();


            if (AnyTranslation)
                return EtTanslations;
            else
                return null;
        }

  



        /// <summary>
        /// 
        /// </summary>
        static public void EnsureMinimalBalance(this IGrid g) {
            using(var tr = new FuncTrace()) {
                int J = g.CellPartitioning.LocalLength;
                int MpiSz = g.CellPartitioning.MpiSize;
                int MpiRk = g.CellPartitioning.MpiRank;

                long JG = g.CellPartitioning.TotalLength;
                if(JG < MpiSz) {
                    throw new ApplicationException($"Grid contains only {JG} cells, but {MpiSz} processors are used; unable to assign at least one cell to each processor; un-able to continue.");
                }

                // determine if there is any processor which has no cell at all
                bool locallyFucked = J <= 0;
                bool globallyFucked = locallyFucked.MPIOr();

                tr.Info($"Processor {MpiRk} contains {J} cells");
                if(!globallyFucked) {
                    tr.Info("All processors contain at least one cell - no re-balancing required.");
                    return; // nothing to do
                }

                g.InvalidateGridData();

                int _i0(int rnk) {
                    return (rnk * J) / MpiSz;
                }
                int _iE(int rnk) {
                    return ((rnk + 1) * J) / MpiSz;
                }
                int FindRank(int j) {
                    for(int r = 0; r < MpiSz; r++) {
                        if(_i0(r) <= j && j < _iE(r))
                            return r;
                    }
                    throw new ApplicationException();
                }



                int[] newPart = new int[J];
                long j0G = g.CellPartitioning.i0;
                for(int j = 0; j < J; j++) {
                    newPart[j] = FindRank(j);
                }

                g.RedistributeGrid(newPart);

            }
        }
    }
}

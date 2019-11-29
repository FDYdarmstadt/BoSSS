using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Utils;
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
        /// to the <see cref="EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        /// <param name="EdgeTagFunc"></param>
        static public void DefineEdgeTags(this IGrid g, Func<Vector, byte> EdgeTagFunc) {
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
        /// to the <see cref="EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        /// <param name="EdgeTagFunc"></param>
        static public void DefineEdgeTags(this IGrid g, Func<double[], byte> EdgeTagFunc) {
            g.DefineEdgeTags((Vector v) => EdgeTagFunc(v));
        }
        

        /// <summary>
        /// sets values for <see cref="Cell.CellFaceTags"/> by using a
        /// <paramref name="EdgeTagFunc"/>-function; also adds entries with empty names
        /// to the <see cref="EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        /// <param name="EdgeTagFunc"></param>
        static public void DefineEdgeTags(this IGrid g, Func<double[], string> EdgeTagFunc) {
            g.DefineEdgeTags((Vector v) => EdgeTagFunc(v));
        }
        

        /// <summary>
        /// sets values for <see cref="Cell.CellFaceTags"/> by using a
        /// <paramref name="EdgeTagFunc"/>-function; also adds entries with empty names
        /// to the <see cref="EdgeTagNames"/>-dictionary, if the edge tag
        /// returned by the <paramref name="EdgeTagFunc"/>-function is not in
        /// the dictionary
        /// </summary>
        /// <param name="EdgeTagFunc"></param>
        static public void DefineEdgeTags(this IGrid g, Func<Vector, string> EdgeTagFunc) {

            int D = g.SpatialDimension;
            double[] x = new double[D];
            MultidimensionalArray GlobalVerticesOut = MultidimensionalArray.CreateWrapper(x, 1, D);

            Dictionary<string, byte> EdgeTagNames_Reverse = new Dictionary<string, byte>();
            EdgeTagNames_Reverse.Add("inner edge", 0);
            foreach(var kv in g.EdgeTagNames) {
                if(kv.Key == 0) {
                    
                } else {
                    EdgeTagNames_Reverse.Add(kv.Value, kv.Key);
                }
            }
            int[] etUseCount = new int[byte.MaxValue];



            GridCommons baseGrid = null;
            if(g is GridCommons gcm) {
                baseGrid = gcm;
            } else if(g is Aggregation.AggregationGrid ag) {
                baseGrid = ag.RootGrid as GridCommons;
            }

            bool GridChanged = false;

            var GrdDatTmp = g.iGridData;
            int NoOfEdges = GrdDatTmp.iGeomEdges.Count;
            for (int iEdge = 0; iEdge < NoOfEdges; ++iEdge) {
                if (GrdDatTmp.iGeomEdges.IsEdgeBoundaryEdge(iEdge)) {

                    int jCell = GrdDatTmp.iGeomEdges.CellIndices[iEdge, 0];
                    int iFace = GrdDatTmp.iGeomEdges.FaceIndices[iEdge, 0];
                    var KRef = GrdDatTmp.iGeomCells.GetRefElement(jCell);

                    GrdDatTmp.TransformLocal2Global(KRef.GetFaceCenter(iFace), GlobalVerticesOut, jCell);
                    string EdgeTagName = EdgeTagFunc(x);

                    if (!EdgeTagNames_Reverse.ContainsKey(EdgeTagName)) {
                        int NewTag = EdgeTagNames_Reverse.Count;
                        Debug.Assert(NewTag > 0);
                        if (NewTag >= GridCommons.FIRST_PERIODIC_BC_TAG)
                            throw new ApplicationException("To many different edge tag names; at maximum " + (Classic.GridCommons.FIRST_PERIODIC_BC_TAG - 1) + " different tags for non-periodic boundaries.");
                        EdgeTagNames_Reverse.Add(EdgeTagName, (byte)NewTag);
                    }
                    var et = EdgeTagNames_Reverse[EdgeTagName];

                    if (baseGrid != null) {
                        var _Cell = baseGrid.Cells[jCell];
                        var allCFTs = _Cell.CellFaceTags;
                        int found = 0;
                        for(int i = 0; i < allCFTs.Length; i++) {
                            if(allCFTs[i].NeighCell_GlobalID < 0 && allCFTs[i].FaceIndex == iFace) {
                                found++;
                                if(allCFTs[i].EdgeTag == et) {
                                    // nop
                                } else {
                                    allCFTs[i].EdgeTag = et;
                                    GridChanged = true;
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
                            GridChanged = true;
                        }
                    }

                    etUseCount[et]++;
                    GrdDatTmp.iGeomEdges.EdgeTags[iEdge] = et;
                }
            }

            g.EdgeTagNames.Clear();
            foreach(var kv in EdgeTagNames_Reverse) {
                if(kv.Value == 0 || etUseCount[kv.Value] > 0)
                    g.EdgeTagNames.Add(kv.Value, kv.Key);
            }

            if (GridChanged) {
                Console.WriteLine("Grid Edge Tags changed.");
            }
        }
    }
}

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
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// masks some edges in a <see cref="GridData"/>-object
    /// </summary>
    public class EdgeMask : ExecutionMask {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="grddat"></param>
        /// <param name="mask">
        /// a "true" entry for all edges in grid <paramref name="grddat"/> that should be in the mask;
        /// The length of this array must not exceed <see cref="GridData.EdgeData.Count"/>
        /// </param>
        public EdgeMask(IGridData grddat, BitArray mask) :
            base(grddat, mask) {
            if (mask.Length != grddat.iLogicalEdges.Count)
                throw new ArgumentException();
        }

        static BitArray MaskFromSelector(IGridData grddat, Func<double[], bool> GeomSelector) {
            int NoOfEdges = grddat.iLogicalEdges.Count;
            BitArray mask = new BitArray(NoOfEdges);
            int D = grddat.SpatialDimension;

            for(int i = 0; i < NoOfEdges; i++) {
                int iEref = grddat.iGeomEdges.GetRefElementIndex(i);
                var Eref = grddat.iGeomEdges.EdgeRefElements[iEref];
                NodeSet Center = Eref.Center;

                MultidimensionalArray GlobalCoord = grddat.GlobalNodes.GetValue_EdgeSV(Center, i, 1);
                double[] X = new double[D];
                Debug.Assert(GlobalCoord.Dimension == 3);
                Debug.Assert(GlobalCoord.GetLength(0) == 1);
                Debug.Assert(GlobalCoord.GetLength(0) == 1);
                Debug.Assert(GlobalCoord.GetLength(0) == D);

                for(int d = 0; d < D; d++) {
                    X[d] = GlobalCoord[0, 0, d];
                }

                mask[i] = GeomSelector(X);
            }

            return mask;
        }


        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="grddat"></param>
        /// <param name="GeomSelector">
        /// Retuns true, if the edge with given center coordinate should be in the mask, otherwise false.
        /// </param>
        public EdgeMask(IGridData grddat, Func<double[], bool> GeomSelector) :
            base(grddat, MaskFromSelector(grddat, GeomSelector)) {
        }


        /// <summary>
        /// ctor
        /// </summary>
        public EdgeMask(IGridData grddat, int[] Sequence) :
            base(grddat, Sequence) {
        }

        


        /// <summary>
        /// compiles an edge mask from a set of chunks
        /// </summary>
        /// <param name="parts">
        /// a list of chunks, which may overlap
        /// </param>
        /// <param name="grddat">
        /// the grid that this mask will be associated with;
        /// </param>
        public EdgeMask(IGridData grddat, params Chunk[] parts)
            : this(grddat, (IEnumerable<Chunk>)parts) {
        }

        /// <summary>
        /// Creates an edge mask for all edges with an edge tag that
        /// corresponds to the given <paramref name="edgeTagName"/>
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="edgeTagName"></param>
        public EdgeMask(IGridData gridData, string edgeTagName)
            : this(gridData, GetBitMaskForEdgeTag(gridData, edgeTagName)) {
        }

        /// <summary>
        /// compiles an quadrature execution mask from a set of chunks
        /// </summary>
        /// <param name="Parts">
        /// a list of chunks, which may overlap
        /// </param>
        /// <param name="grddat">
        /// the grid that this mask will be associated with;
        /// </param>
        protected EdgeMask(IGridData grddat, IEnumerable<Chunk> Parts)
            : this(grddat, FromChunkEnum(Parts)) {
        }

        /// <summary>
        /// Retrieves an empty edge mask;
        /// </summary>
        /// <param name="grdDat">
        /// grid that the returned mask will be assigned to
        /// </param>
        static public EdgeMask GetEmptyMask(IGridData grdDat) {
            BitArray ba = new BitArray(grdDat.iLogicalEdges.Count, false);
            return new EdgeMask(grdDat, ba);
        }

        /// <summary>
        /// complimentary edge mask
        /// </summary>
        public EdgeMask Complement() {
            return base.Complement<EdgeMask>();
        }

        /// <summary>
        /// Retrieves a mask containing all edges (i.e. returns the
        /// complement of <see cref="GetEmptyMask"/>)
        /// </summary>
        /// <param name="gridDat">
        /// Grid data that the returned mask will be assigned with
        /// </param>
        /// <returns>A full mask</returns>
        public static EdgeMask GetFullMask(IGridData gridDat) {
            return new EdgeMask(gridDat, new Chunk {
                i0 = 0,
                Len = gridDat.iLogicalEdges.Count
            });
        }

        /// <summary>
        /// like ctor;
        /// </summary>
        protected override ExecutionMask CreateInstance(IGridData grdDat, BitArray mask) {
            return new EdgeMask(grdDat, mask);
        }

        /// <summary>
        /// see <see cref="ExecutionMask.GetTotalNumberOfElements"/>
        /// </summary>
        protected override int GetTotalNumberOfElements(IGridData gridData) {
            return gridData.iLogicalEdges.Count;
        }


        /// <summary>
        /// computes a cell mask that contains all cells with an edge in this edge mask
        /// </summary>
        public CellMask GetAdjacentCells(Grid.Classic.GridData gridData) {
            int J = gridData.Cells.NoOfLocalUpdatedCells;
            int[,] AllEdges = gridData.Edges.CellIndices;
            BitArray mask = new BitArray(J);

            foreach (var ch in this) {
                int L = ch.i0 + ch.Len;
                for (int l = ch.i0; l < L; l++) {
                    int cellIn = AllEdges[l, 0];
                    int cellOt = AllEdges[l, 1];

                    if (cellIn < J)
                        mask[cellIn] = true;
                    if (cellOt >= 0 && cellOt < J)
                        mask[cellOt] = true;

                }
            }

            return new CellMask(gridData, mask);
        }

        /// <summary>
        /// a cell is in the returned cell mask if
        /// it is a neighbor cell of this edge mask and if it is also neighbor of <paramref name="X"/>
        /// </summary>
        public CellMask GetAdjacentCellsCond(Grid.Classic.GridData gridData, CellMask X) {
            int J = gridData.Cells.NoOfLocalUpdatedCells;
            int[,] AllEdges = gridData.Edges.CellIndices;
            BitArray mask = new BitArray(J);
            BitArray Xmask = X.GetBitMaskWithExternal();

            foreach (var ch in this) {
                int L = ch.i0 + ch.Len;
                for (int l = ch.i0; l < L; l++) {
                    int cellIn = AllEdges[l, 0];
                    int cellOt = AllEdges[l, 1];

                    bool Out = (cellOt > 0 && Xmask[cellOt]);
                    bool In_ = Xmask[cellIn];
                    if (!(Out || In_))
                        continue;


                    if (cellIn < J)
                        mask[cellIn] = true;
                    if (cellOt >= 0 && cellOt < J)
                        mask[cellOt] = true;

                }
            }

            return new CellMask(gridData, mask);

        }


        /// <summary>
        /// Writes the mid-points of each edge in this mask to a text file with
        /// the given name. Useful for debugging purposes.
        /// </summary>
        public override void SaveToTextFile(string fileName, bool WriteHeader = true, params ItemInfo[] infoFunc) {
            int D = GridData.SpatialDimension;
            int LI = infoFunc.Length;
            using (var file = new StreamWriter(fileName)) {
                if (WriteHeader) {

                    file.Write("Edge");
                    switch (D) {
                        case 1:
                        file.Write("\tx");
                        break;

                        case 2:
                        file.Write("\tx\ty");
                        break;

                        case 3:
                        file.Write("\tx\ty\tz");
                        break;

                        default:
                        throw new Exception();
                    }


                    for (int i = 0; i < LI; i++) {
                        file.Write("\ti(" + i + ")");
                    }
                    file.WriteLine();
                }
                double[] x = new double[D];

                MultidimensionalArray localCenterEdge = MultidimensionalArray.Create(1, D - 1);
                MultidimensionalArray localCenterVolume = MultidimensionalArray.Create(1, D);
                foreach (Chunk chunk in this) {
                    for (int i = 0; i < chunk.Len; i++) {
                        int edgeLog = chunk.i0 + i;
                        foreach (int edge in this.GridData.GetGeometricEdgeIndices(edgeLog)) {
                            int iTrafo = GridData.iGeomEdges.Edge2CellTrafoIndex[edge, 0];
                            int localEdge = GridData.iGeomEdges.FaceIndices[edge, 0];
                            int cell = GridData.iGeomEdges.CellIndices[edge, 0];
                            RefElement KrefCell = this.GridData.iGeomCells.GetRefElement(cell);

                            GridData.iGeomEdges.Edge2CellTrafos[iTrafo].Transform(localCenterEdge, localCenterVolume);

                            MultidimensionalArray globalCenter = MultidimensionalArray.Create(1, 1, D);
                            GridData.TransformLocal2Global(new NodeSet(KrefCell, localCenterVolume), cell, 1, globalCenter, 0);

                            file.Write(chunk.i0 + i);
                            for (int d = 0; d < D; d++) {
                                file.Write("\t" + globalCenter[0, 0, d].ToString("e", NumberFormatInfo.InvariantInfo));
                            }

                            if (LI > 0) {
                                for (int d = 0; d < D; d++)
                                    x[d] = globalCenter[0, 0, d];
                                for (int li = 0; li < LI; li++) {
                                    double info_i = infoFunc[li](x, edgeLog, edge);
                                    file.Write("\t" + info_i.ToString("e", NumberFormatInfo.InvariantInfo));
                                }
                            }

                            file.WriteLine();
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Constructs the bit mask for all edges with an edge tag that
        /// corresponds to the given <paramref name="edgeTagName"/>
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="edgeTagName"></param>
        /// <returns></returns>
        private static BitArray GetBitMaskForEdgeTag(IGridData gridData, string edgeTagName) {
            byte edgeTag;
            try {
                edgeTag = gridData.EdgeTagNames.First(item => item.Value.Equals(edgeTagName)).Key;
            } catch (InvalidOperationException e) {
                throw new ArgumentException("An edge tag with name \"" + edgeTagName + "\" does not exist", "edgeTagName", e);
            }

            byte[] edgeTags = gridData.iGeomEdges.EdgeTags;
            BitArray maskArray = new BitArray(edgeTags.Length);
            for (int i = 0; i < edgeTags.Length; i++) {
                if (edgeTag == edgeTags[i]) {
                    maskArray[i] = true;
                }
            }

            return maskArray;
        }
    }
}

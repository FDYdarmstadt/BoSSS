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
        /// <param name="mt">
        /// <see cref="ExecutionMask.MaskType"/>
        /// </param>
        public EdgeMask(IGridData grddat, BitArray mask, MaskType mt = MaskType.Logical) :
            base(grddat, mask, mt) //
        {
            if (mask.Length != this.GetUpperIndexBound(grddat))
                throw new ArgumentException("Mismatch in number of edges/length of input bitmask.");
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
                Debug.Assert(GlobalCoord.GetLength(1) == 1);
                Debug.Assert(GlobalCoord.GetLength(2) == D);

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
        /// Returns true, if the edge with given center coordinate should be in the mask, otherwise false.
        /// </param>
        /// <param name="mt">
        /// <see cref="ExecutionMask.MaskType"/>
        /// </param>
        public EdgeMask(IGridData grddat, Func<double[], bool> GeomSelector, MaskType mt = MaskType.Logical) :
            base(grddat, MaskFromSelector(grddat, GeomSelector), mt) {
        }


        /// <summary>
        /// ctor
        /// </summary>
        public EdgeMask(IGridData grddat, int[] Sequence, MaskType mt = MaskType.Logical) :
            base(grddat, Sequence, mt) {
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
            : this(grddat, (IEnumerable<Chunk>)parts, MaskType.Logical) {
        }

        /// <summary>
        /// compiles an edge mask from a single chunk
        /// </summary>
        /// <param name="part">
        /// a chunk
        /// </param>
        /// <param name="grddat">
        /// the grid that this mask will be associated with;
        /// </param>
        /// <param name="mt"></param>
        public EdgeMask(IGridData grddat, Chunk part, MaskType mt = MaskType.Logical)
            : this(grddat, new Chunk[] { part }, mt) {
        }

        /// <summary>
        /// Creates an edge mask for all edges with an edge tag that
        /// corresponds to the given <paramref name="edgeTagName"/>
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="edgeTagName"></param>
        public EdgeMask(IGridData gridData, string edgeTagName)
            : this(gridData, GetBitMaskForEdgeTag(gridData, edgeTagName), MaskType.Logical) {
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
        /// <param name="mt">
        /// <see cref="ExecutionMask.MaskType"/>
        /// </param>
        protected EdgeMask(IGridData grddat, IEnumerable<Chunk> Parts, MaskType mt = MaskType.Logical)
            : this(grddat, FromChunkEnum(Parts), mt) {
        }

        /// <summary>
        /// Retrieves an empty edge mask;
        /// </summary>
        /// <param name="grdDat">
        /// grid that the returned mask will be assigned to
        /// </param>
        /// <param name="mt">
        /// <see cref="ExecutionMask.MaskType"/>
        /// </param>
        static public EdgeMask GetEmptyMask(IGridData grdDat, MaskType mt = MaskType.Logical) {
            //BitArray ba = new BitArray(grdDat.iLogicalEdges.Count, false);
            return new EdgeMask(grdDat, new int[0], mt);
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
        /// <param name="mt">
        /// <see cref="ExecutionMask.MaskType"/>
        /// </param>
        public static EdgeMask GetFullMask(IGridData gridDat, MaskType mt) {

            switch (mt) {
                case MaskType.Logical: {
                    int L = gridDat.iLogicalEdges.Count;
                    return new EdgeMask(gridDat, new[]{ new Chunk {
                            i0 = 0,
                            Len = L
                        } }, mt);
                }
                case MaskType.Geometrical: {
                    return GetFullMask(gridDat, MaskType.Logical).ToGeometicalMask();
                }
                default: throw new NotImplementedException();
            }
        }

        /// <summary>
        /// like ctor;
        /// </summary>
        protected override ExecutionMask CreateInstance(BitArray mask, MaskType mt) {
            return new EdgeMask(base.GridData, mask, mt);
        }

        /// <summary>
        /// see <see cref="ExecutionMask.GetUpperIndexBound"/>
        /// </summary>
        protected override int GetUpperIndexBound(IGridData gridData) {
            switch (base.MaskType) {
                case MaskType.Logical: return gridData.iLogicalEdges.Count;
                case MaskType.Geometrical: return gridData.iGeomEdges.Count;
                default: throw new NotImplementedException();
            }
        }


        /// <summary>
        /// computes a cell mask that contains all cells with an edge in this edge mask
        /// </summary>
        public CellMask GetAdjacentCells() {
            var gridData = (Grid.Classic.GridData)(base.GridData);

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

            return new CellMask(gridData, mask, MaskType.Logical);
        }

        /// <summary>
        /// a cell is in the returned cell mask if
        /// it is a neighbor cell of this edge mask and if it is also neighbor of <paramref name="X"/>
        /// </summary>
        public CellMask GetAdjacentCellsCond(CellMask X) {
            var gridData = (Grid.Classic.GridData)(base.GridData);

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

            return new CellMask(gridData, mask, MaskType.Logical);

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
                            GridData.TransformLocal2Global(new NodeSet(KrefCell, localCenterVolume, false), cell, 1, globalCenter, 0);

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


        /// <summary>
        /// Converts this  
        /// from a logical (<see cref="IGridData.iLogicalEdges"/>) mask
        /// to a geometrical (<see cref="IGridData.iGeomEdges"/>) mask.
        /// </summary>
        /// <returns></returns>
        public EdgeMask ToGeometicalMask() {
            if(base.MaskType != MaskType.Logical)
                throw new NotSupportedException();

            if(base.GridData is Grid.Classic.GridData) 
                // logical and geometrical cells are identical - return a clone of this mask
                return new EdgeMask(base.GridData, base.Sequence, MaskType.Geometrical);

            int Jg = GridData.iGeomEdges.Count;
            int[][] jl2jg = GridData.iLogicalEdges.EdgeToParts;
            BitArray ba = new BitArray(Jg);
            foreach(Chunk c in this) { // loop over chunks of logical edges...
                for(int jl = 0; jl < c.JE; jl++) { // loop over edges in chunk...
                    foreach(int jg in jl2jg[jl]) { // loop over geometrical edges in logical cell 'jl'...
                        ba[jg] = true;
                    }
                }
            }

            return new EdgeMask(base.GridData, ba, MaskType.Geometrical);
        }

    }
}

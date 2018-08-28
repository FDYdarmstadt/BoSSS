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

using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi {
    partial class VoronoiGrid {

        public IGeometricalCellsData iGeomCells {
            get {
                return m_CellData;
            }
        }
        public ILogicalCellData iLogicalCells {
            get {
                return m_CellData;
            }
        }

        

        CellData m_CellData;

        class CellData : IGeometricalCellsData, ILogicalCellData {

            /// <summary>
            /// The reference element for the Voronoi cell parts
            /// </summary>
            public RefElement[] RefElements {
                get;
                internal set;
            }


            internal VoronoiGrid m_Owner;

            public int[][] CellNeighbours {
                get;
                internal set;
            }

            public int Count {
                get {
                    return NoOfExternalCells + NoOfLocalUpdatedCells;
                }
            }

            public int NoOfExternalCells {
                get {
                    return 0;
                }
            }

            public int NoOfLocalUpdatedCells {
                get {
                    return m_Owner.DelaunayVertices.NoOfRows;
                }
            }

            public int[][] CellVertices {
                get;
                internal set;
            }

            public int[][] AggregateCellToParts {
                get;
                internal set;
            }


            /// <summary>
            /// Linear part of transformation.
            /// - 1st index: Voronoi cell 
            /// - 2nd index, 3rd index: row and column index of transformation matrix.
            /// </summary>
            public MultidimensionalArray PartTransformation {
                get;
                internal set;
            }

            /// <summary>
            /// Affine part of transformation.
            /// </summary>
            public MultidimensionalArray PartCenter {
                get;
                internal set;
            }

            /// <summary>
            /// Linear part of transformation.
            /// - 1st index: Voronoi cell 
            /// - 2nd index, 3rd index: row and column index of transformation matrix.
            /// </summary>
            public MultidimensionalArray BoundingBoxTransformation {
                get;
                internal set;
            }

            /// <summary>
            /// Affine part of transformation.
            /// </summary>
            public MultidimensionalArray BoundingBoxCenter {
                get;
                internal set;
            }

            /// <summary>
            /// see <see cref="CellInfo"/>
            /// </summary>
            public CellInfo[] InfoFlags {
                get;
                internal set;
            }

            public int[][] Cells2Edges {
                get;
                internal set;
            }

            public MultidimensionalArray InverseTransformation {
                get {
                    throw new NotImplementedException();
                }
            }

            public MultidimensionalArray JacobiDet {
                get {
                    throw new NotImplementedException();
                }
            }

            public MultidimensionalArray Transformation {
                get {
                    throw new NotImplementedException();
                }
            }

            public MultidimensionalArray h_min {
                get {
                    throw new NotImplementedException();
                }
            }

            public MultidimensionalArray h_max {
                get {
                    throw new NotImplementedException();
                }
            }

            public int[] GeomCell2LogicalCell => throw new NotImplementedException();


            /// <summary>
            /// Computes the bounding box of cell <paramref name="j"/>
            /// </summary>
            /// <param name="j">local cell index</param>
            /// <param name="bb">
            /// on exit, the bounding box of cell j.
            /// </param>
            public void GetCellBoundingBox(int j, BoundingBox bb) {
                int D = bb.D;
                if (bb.D != m_Owner.SpatialDimension)
                    throw new ArgumentException("wrong dimension of bounding box.");
                bb.Clear();

                MultidimensionalArray Points = m_Owner.m_VertexData.Coordinates;
                double[] pt = new double[D];

                foreach (int iVtx in this.CellVertices[j]) {
                    Points.GetRow(iVtx, pt, 0, 1);
                    bb.AddPoint(pt);
                }
            }

            public int GetRefElementIndex(int jCell) {
                throw new NotImplementedException();
            }

            public bool IsCellAffineLinear(int jCell) {
                throw new NotImplementedException();
            }

            public RefElement GetRefElement(int j) {
                throw new NotImplementedException();
            }

            public double GetCellVolume(int j) {
                throw new NotImplementedException();
            }

            public CellType GetCellType(int jCell) {
                throw new NotImplementedException();
            }

            public long GetGlobalID(int j) {
                throw new NotImplementedException();
            }

            public int GetNoOfSimilarConsecutiveCells(CellInfo mask, int j0, int Lmax) {
                throw new NotImplementedException();
            }

            public CellMask GetCells4Refelement(RefElement Kref) {
                throw new NotImplementedException();
            }

            public int GetInterpolationDegree(int jCell) {
                throw new NotImplementedException();
            }
        }
    }
}

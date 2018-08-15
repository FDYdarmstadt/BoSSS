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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Caching;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using ilPSP;

namespace BoSSS.Foundation.Grid.Aggregation {
    partial class AggregationGrid {

        public IGeometricalEdgeData iGeomEdges {
            get {
                return m_GeomEdgeData;
            }
        }

        public ILogicalEdgeData iLogicalEdges {
            get {
                return m_LogEdgeData;
            }
        }

        GeomEdgeData m_GeomEdgeData;

        LogEdgeData m_LogEdgeData;

        /// <summary>
        /// Just a wrapper/proxy around the geometrical edge data (<see cref="IGridData.iGeomEdges"/>) 
        /// of the parent grid (<see cref="AggregationGrid.ParentGrid"/>).
        /// </summary>
        class GeomEdgeData : IGeometricalEdgeData {
            internal AggregationGrid m_Owner;
            
            public int[,] CellIndices {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.CellIndices;
                }
            }

            public int Count {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.Count;
                }
            }

            public int e2C_offet {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.e2C_offet;
                }
            }

            public int[,] Edge2CellTrafoIndex {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.Edge2CellTrafoIndex;
                }
            }

            public IList<AffineTrafo> Edge2CellTrafos {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.Edge2CellTrafos;
                }
            }

            public IList<int> Edge2CellTrafosRefElementIndices {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.Edge2CellTrafosRefElementIndices;
                }
            }

            public MultidimensionalArray Edge2CellTrafos_SqrtGramian {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.Edge2CellTrafos_SqrtGramian;
                }
            }

            public RefElement[] EdgeRefElements {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.EdgeRefElements;
                }
            }

            public byte[] EdgeTags {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.EdgeTags;
                }
            }

            public byte[,] FaceIndices {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.FaceIndices;
                }
            }

            public EdgeInfo[] Info {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.Info;
                }
            }

            public EdgeNormalsCacheLogic_CNsFace NormalsCache {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.NormalsCache;
                }
            }

            public MultidimensionalArray SqrtGramian {
                get {
                    return m_Owner.ParentGrid.iGeomEdges.SqrtGramian;
                }
            }


            EdgeMask[] m_Edges4RefElement;

            public EdgeMask GetEdges4RefElement(RefElement Kref) {
                if (m_Edges4RefElement == null)
                    m_Edges4RefElement = new EdgeMask[this.EdgeRefElements.Length];

                int iKref = this.EdgeRefElements.IndexOf(Kref, (a, b) => object.ReferenceEquals(a, b));

                if (m_Edges4RefElement[iKref] == null) {
                    EdgeMask parrEdg = m_Owner.ParentGrid.iGeomEdges.GetEdges4RefElement(Kref);
                    EdgeMask thisAll = EdgeMask.GetFullMask(this.m_Owner, MaskType.Geometrical);
                    if (parrEdg.MaskType != MaskType.Geometrical)
                        throw new ApplicationException("expecting a geometrical mask");
                    if (thisAll.MaskType != MaskType.Geometrical)
                        throw new ApplicationException("expecting a geometrical mask");

                    BitArray parrBitmask = parrEdg.GetBitMask().CloneAs();
                    BitArray thisBitMask = thisAll.GetBitMask();
                    Debug.Assert(parrBitmask.Length == thisBitMask.Length);

                    BitArray intersect = parrBitmask.And(thisBitMask);
                    Debug.Assert(object.ReferenceEquals(intersect, parrBitmask));

                    m_Edges4RefElement[iKref] = new EdgeMask(m_Owner, intersect, MaskType.Geometrical);
                }

                return m_Edges4RefElement[iKref];
            }

            public void GetNormalsForCell(NodeSet Nodes, int jCell, int iFace, MultidimensionalArray NormalsOut) {
                m_Owner.ParentGrid.iGeomEdges.GetNormalsForCell(Nodes, jCell, iFace, NormalsOut);
            }

            public void GetNormalsForCell(NodeSet Nodes, int jCell, int iFace, MultidimensionalArray NormalsOut, MultidimensionalArray QuadMetric, int Offset) {
                m_Owner.ParentGrid.iGeomEdges.GetNormalsForCell(Nodes, jCell, iFace, NormalsOut, QuadMetric, Offset);
            }

            public int GetRefElementIndex(int e) {
                return m_Owner.ParentGrid.iGeomEdges.GetRefElementIndex(e);
            }

            public bool IsEdgeAffineLinear(int e) {
                return m_Owner.ParentGrid.iGeomEdges.IsEdgeAffineLinear(e);
            }


        }

        class LogEdgeData : ILogicalEdgeData {
            public int[,] CellIndices {
                get;
                internal set;
            }

            public int Count {
                get {
                    return CellIndices.GetLength(0);
                }
            }

            public int[][] EdgeToParts {
                get;
                internal set;
            }

            public EdgeMask GetEdges4RefElement(RefElement Kref) {
                throw new NotImplementedException();
            }
        }

    }
}

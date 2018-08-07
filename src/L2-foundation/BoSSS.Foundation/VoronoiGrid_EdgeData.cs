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
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Foundation.Caching;

namespace BoSSS.Foundation.Grid.Voronoi {
    partial class VoronoiGrid {

        public IGeometricalEdgeData iGeomEdges {
            get {
                return m_GeomEdges;
            }
        }

        GeomEdgeData m_GeomEdges;

        class GeomEdgeData : IGeometricalEdgeData {
            public int Count {
                get {
                    throw new NotImplementedException();
                }
            }

            public int e2C_offet {
                get {
                    throw new NotImplementedException();
                }
            }

            public int[,] Edge2CellTrafoIndex {
                get;
                internal set;
            }

            public IList<AffineTrafo> Edge2CellTrafos {
                get {
                    throw new NotImplementedException();
                }
            }

            public IList<int> Edge2CellTrafosRefElementIndices {
                get {
                    throw new NotImplementedException();
                }
            }

            public RefElement[] EdgeRefElements {
                get;
                internal set;
            }

            public byte[] EdgeTags {
                get;
                internal set;
            }

            public byte[,] FaceIndices {
                get {
                    throw new NotImplementedException();
                }
            }

            public EdgeInfo[] Info {
                get;
                internal set;
            }

            public int GetRefElementIndex(int e) {
                throw new NotImplementedException();
            }

            public bool IsEdgeAffineLinear(int e) {
                throw new NotImplementedException();
            }

            public void GetNormalsForCell(NodeSet Nodes, int jCell, int iFace, MultidimensionalArray NormalsOut, MultidimensionalArray QuadMetric, int Offset) {
                throw new NotImplementedException();
            }

            public void GetNormalsForCell(NodeSet Nodes, int jCell, int iFace, MultidimensionalArray NormalsOut) {
                throw new NotImplementedException();
            }

            public EdgeMask GetEdges4RefElement(RefElement Kref) {
                throw new NotImplementedException();
            }

            public int[][] VertexIndices {
                get;
                internal set;
            }

            public int[,] CellIndices {
                get {
                    throw new NotImplementedException();
                }
            }

            public MultidimensionalArray SqrtGramian {
                get {
                    throw new NotImplementedException();
                }
            }

            public EdgeNormalsCacheLogic_CNsFace NormalsCache {
                get {
                    throw new NotImplementedException();
                }
            }

            public MultidimensionalArray Edge2CellTrafos_SqrtGramian {
                get {
                    throw new NotImplementedException();
                }
            }
        }

        public ILogicalEdgeData iLogicalEdges {
            get {
                return m_LogEdges;
            }
        }

        LogEdgeData m_LogEdges;

        class LogEdgeData : ILogicalEdgeData {
            public int[,] CellIndices {
                get;
                internal set;
            }

            public int Count {
                get {
                    throw new NotImplementedException();
                }
            }

            public int[][] EdgeToParts {
                get {
                    throw new NotImplementedException();
                }
            }

            public EdgeMask GetEdges4RefElement(RefElement Kref) {
                throw new NotImplementedException();
            }
        }

    }
}

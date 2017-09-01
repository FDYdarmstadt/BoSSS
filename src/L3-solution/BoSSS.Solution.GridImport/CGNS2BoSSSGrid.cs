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
using System.Linq;
using System.Text;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP.Utils;
using System.Diagnostics;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Solution.GridImport {
    
    
    public static class CGNS2BoSSSGrid {

        static public GridCommons ToBoSSSGrid(this Cgns cgns, int BaseNodeIndex = 0) {
            var b = cgns.base_t.ElementAt(BaseNodeIndex);
            int D = b.CellDimension;

            var grid = CreateBoSSSGrid(b);

            long GlobalID_cnt = 0;

            int[] NodesPerZone = b.zones.Select(delegate(Cgns.Zone_t z) {
                if(!(z is Cgns.UnstructuredZone_t))
                    throw new NotSupportedException("Only unstructured Grids are supported for CGNS ot BoSSS conversion.");

                Cgns.UnstructuredZone_t uz = (Cgns.UnstructuredZone_t)z;

                int Ncoords = uz.coord[0].Length;
                for(int d = 1; d < uz.coord.Length; d++) {
                    if(uz.coord[d].Length != Ncoords) {
                        throw new ApplicationException("invalid cgns data.");
                    }
                }

                return Ncoords;
            }).ToArray();


            // import cells
            // ============

            {

                List<Cell> Cells = new List<Cell>();
                int NodeOffset = 0;
                for ( int iZone = 0; iZone < b.zones.Length; iZone++) {
                    var z = b.zones[iZone];

                    if (!(z is Cgns.UnstructuredZone_t))
                        throw new NotSupportedException("Only unstructured Grids are supported for CGNS ot BoSSS conversion.");

                    Cgns.UnstructuredZone_t uz = (Cgns.UnstructuredZone_t)z;

                    int NumberOfNodes = NodesPerZone[iZone];


                    foreach (Cgns.Elements_t elm in uz.elements) {
                        if (elm.element_type.Dimension() != D)
                            continue;

                        CellType BoSSSCelType = elm.element_type.ToBoSSSCellType();

                        for (int j = 0; j < elm.ielem.GetLength(0); j++) {
                            Cell cell = new Cell();
                            cell.GlobalID = GlobalID_cnt;
                            cell.Type = BoSSSCelType;

                            int[] CgnsNodes = elm.ielem.GetRow(j);
                            foreach(int iNode in CgnsNodes)
                                if(iNode > NumberOfNodes)
                                    throw new ArgumentException("CGNS node index out of range.");
                            cell.NodeIndices = GetBoSSSConnectivityNodes(elm.element_type, CgnsNodes);
                            for(int iNode = 0; iNode < cell.NodeIndices.Length; iNode++) {
                                cell.NodeIndices[iNode] += NodeOffset;
                            }

                            cell.TransformationParams = GetTransformationParams(elm.element_type, CgnsNodes, uz.coord, D);

                            GlobalID_cnt++;

                            Cells.Add(cell);
                        }
                    }

                    NodeOffset += NumberOfNodes;
                }
                grid.Cells = Cells.ToArray();
            }

            // import Bc
            // =========
            {
                List<BCElement> BcCells = new List<BCElement>();

                byte EdgeTagCnt = 0;

                int NodeOffset = 0;
                for(int iZone = 0; iZone < b.zones.Length; iZone++) {
                    var z = b.zones[iZone];

                    if(!(z is Cgns.UnstructuredZone_t))
                        throw new NotSupportedException("Only unstructured Grids are supported.");

                    Cgns.UnstructuredZone_t uz = (Cgns.UnstructuredZone_t)z;

                    int NumberOfNodes = NodesPerZone[iZone];

                    foreach (Cgns.BC_t bc in uz.bcs) {

                        Cgns.Elements_t bc_elements = uz.elements.Single(
                            elem => (elem.Name.Equals(bc.Name)));

                        byte EdgeTag;
                        if (!(grid.EdgeTagNames.Values.Contains(bc.Name, (x, y) => x.Equals(y)))) {
                            EdgeTagCnt++;
                            grid.EdgeTagNames.Add(EdgeTagCnt, bc.Name);
                            EdgeTag = EdgeTagCnt;
                        } else {
                            var f = grid.EdgeTagNames.Single(kv => kv.Value.Equals(bc.Name));
                            EdgeTag = f.Key;
                        }

                        CellType BoSSSCelType = bc_elements.element_type.ToBoSSSCellType();
                        
                        for (int j = 0; j < bc_elements.ielem.GetLength(0); j++) {
                            BCElement cell = new BCElement();
                            cell.GlobalID = GlobalID_cnt;
                            cell.Type = BoSSSCelType;
                            cell.EdgeTag = EdgeTag;

                            int[] CgnsNodes = bc_elements.ielem.GetRow(j);
                            foreach(int iNode in CgnsNodes)
                                if(iNode > NumberOfNodes)
                                    throw new ArgumentException("CGNS node index out of range.");
                            cell.NodeIndices = GetBoSSSConnectivityNodes(bc_elements.element_type, CgnsNodes);
                            for(int iNode = 0; iNode < cell.NodeIndices.Length; iNode++) {
                                cell.NodeIndices[iNode] += NodeOffset;
                            }

                            cell.TransformationParams = GetTransformationParams(bc_elements.element_type, CgnsNodes, uz.coord, D);

                            GlobalID_cnt++;

                            BcCells.Add(cell);
                        }
                    }

                    NodeOffset += NumberOfNodes;
                }

                if (BcCells.Count > 0)
                    grid.BcCells = BcCells.ToArray();
            }


            // finalize & return
            // =================

            if(b.zones.Length > 1)
                grid.MergeAndCheckNodes();
            grid.CompressGlobalID();
            grid.CompressNodeIndices();
            grid.CheckAndFixJacobianDeterminat();

            return grid;
        }

        static GridCommons CreateBoSSSGrid(this Cgns.CGNSBase_t b) {
            int D = b.CellDimension;

            List<RefElement> RefElements = new List<RefElement>();
            List<RefElement> EdgeRefElements = new List<RefElement>();

            foreach (var z in b.zones) {
                if (!(z is Cgns.UnstructuredZone_t))
                    throw new NotSupportedException("Only unstructured CGNS Grids are supported.");
                
                Cgns.UnstructuredZone_t uz = (Cgns.UnstructuredZone_t)z;

                foreach (var elm in uz.elements) {
                    if (elm.element_type.Dimension() != D)
                        continue;

                    var Kref = elm.element_type.GetBoSSSElement();

                    if (!RefElements.Contains(Kref, (x, y) => x.GetType().Equals(y.GetType()))) {
                        RefElements.Add(Kref);
                        if (!EdgeRefElements.Contains(Kref.FaceRefElement, (x, y) => x.GetType().Equals(y.GetType()))) {
                            EdgeRefElements.Add(Kref.FaceRefElement);
                        }
                    }
                }
            }

            return new GridCommons(RefElements.ToArray(), EdgeRefElements.ToArray());
        }


        static RefElement GetBoSSSElement(this ElementType_t t) {
            switch (t) {
                case ElementType_t.ElementTypeNull: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.ElementTypeUserDefined: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.NODE: return Point.Instance;
                case ElementType_t.BAR_2: return Line.Instance;
                case ElementType_t.BAR_3: return Line.Instance;
                case ElementType_t.TRI_3: return Triangle.Instance;
                case ElementType_t.TRI_6: return Triangle.Instance;
                case ElementType_t.QUAD_4: return Square.Instance;
                case ElementType_t.QUAD_8: return Square.Instance;
                case ElementType_t.QUAD_9: return Square.Instance;
                case ElementType_t.TETRA_4: return Tetra.Instance;
                case ElementType_t.TETRA_10: return Tetra.Instance;
                case ElementType_t.PYRA_5: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.PYRA_14: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.PENTA_6: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.PENTA_15: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.PENTA_18: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.",t.ToString()));
                case ElementType_t.HEXA_8: return Cube.Instance;
                case ElementType_t.HEXA_20: return Cube.Instance;
                case ElementType_t.HEXA_27: return Cube.Instance;
                case ElementType_t.MIXED: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.NGON_n: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                default: throw new NotImplementedException();
            }
        }


        static CellType ToBoSSSCellType(this ElementType_t t) {
            switch (t) {
                case ElementType_t.ElementTypeNull: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.ElementTypeUserDefined: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.NODE: return CellType.Point;
                case ElementType_t.BAR_2: return CellType.Line_2;
                case ElementType_t.BAR_3: return CellType.Line_3;
                case ElementType_t.TRI_3: return CellType.Triangle_3;
                case ElementType_t.TRI_6: return CellType.Triangle_6;
                case ElementType_t.QUAD_4: return CellType.Square_4;
                case ElementType_t.QUAD_8: return CellType.Square_8;
                case ElementType_t.QUAD_9: return CellType.Square_9;
                case ElementType_t.TETRA_4: return CellType.Tetra_Linear;
                case ElementType_t.TETRA_10: return CellType.Tetra_10;
                case ElementType_t.PYRA_5: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.PYRA_14: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.PENTA_6: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.PENTA_15: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.PENTA_18: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.HEXA_8: return CellType.Cube_8;
                case ElementType_t.HEXA_20: return CellType.Cube_20;
                case ElementType_t.HEXA_27: return CellType.Cube_27;
                case ElementType_t.MIXED: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                case ElementType_t.NGON_n: throw new NotSupportedException(string.Format("CGNS Element {0} not supported in BoSSS.", t.ToString()));
                default: throw new NotImplementedException();
            }
        }


        static int[] GetBoSSSConnectivityNodes(ElementType_t t, int[] cgNodes) {
            if (t == ElementType_t.NODE) {
                return cgNodes.CloneAs();
            }
            if (t == ElementType_t.BAR_2 || t == ElementType_t.BAR_3) {
                return new int[] { cgNodes[0], cgNodes[1] };
            }
            if (t == ElementType_t.TRI_3 || t == ElementType_t.TRI_6) {
                //return new int[] { cgNodes[0], cgNodes[1], cgNodes[2] };
                return new int[] { cgNodes[2], cgNodes[1], cgNodes[0] };
            }
            if (t == ElementType_t.QUAD_4 || t == ElementType_t.QUAD_8 ||t == ElementType_t.QUAD_9) {
                return new int[] { cgNodes[3], cgNodes[0], cgNodes[2], cgNodes[1] };
            }
            if (t == ElementType_t.TETRA_4 || t == ElementType_t.TETRA_10) {
                return new int[] { cgNodes[3], cgNodes[1], cgNodes[0], cgNodes[2] };
            }
            if (t == ElementType_t.HEXA_8 || t == ElementType_t.HEXA_20 || t == ElementType_t.HEXA_27) {
                return new int[] { cgNodes[3], cgNodes[0], cgNodes[7], cgNodes[2], cgNodes[4], cgNodes[5], cgNodes[1], cgNodes[6] };
            }

            throw new NotSupportedException();
        }

        /// <summary>
        /// permutation between CGNS node numbering and BoSSS node numbering
        /// </summary>
        /// <remarks>
        /// Falls irgendwo ein Fehler ist, nicht bei mir beschweren, sondern reparieren:
        /// http://www.grc.nasa.gov/WWW/cgns/CGNS_docs_current/sids/sids.pdf, Seite 20 ff;
        /// </remarks>
        static int[] GetBoSSSNodes(ElementType_t t, int[] cgNodes) {
            if (t == ElementType_t.NODE) {
                return cgNodes.CloneAs();
            }
            // ---------------------
            if (t == ElementType_t.BAR_2) {
                return new int[] { cgNodes[0], cgNodes[1] };
            }
            if (t == ElementType_t.BAR_3) {
                return new int[] { cgNodes[2], cgNodes[0], cgNodes[1] };
            }
            // ---------------------
            if (t == ElementType_t.TRI_3) {
                //return new int[] { cgNodes[0], cgNodes[1], cgNodes[2] };
                return new int[] { cgNodes[2], cgNodes[1], cgNodes[0] };
            }
            if (t == ElementType_t.TRI_6) {
                return new int[] { cgNodes[3], cgNodes[4], cgNodes[5], cgNodes[0], cgNodes[1], cgNodes[2] };
            }
            // ---------------------
            if (t == ElementType_t.QUAD_4) {
                return new int[] { cgNodes[3], cgNodes[0], cgNodes[2], cgNodes[1] };
            }
            if (t == ElementType_t.QUAD_8) {
                return new int[] { cgNodes[6], cgNodes[4], cgNodes[5], cgNodes[7], cgNodes[3], cgNodes[0], cgNodes[2], cgNodes[1] };
            }
            if (t == ElementType_t.QUAD_9) {
                return new int[] { cgNodes[8], cgNodes[6], cgNodes[4], cgNodes[5], cgNodes[7], cgNodes[3], cgNodes[0], cgNodes[2], cgNodes[1] };
            }
            // ---------------------
            if (t == ElementType_t.TETRA_4) {
                return new int[] { cgNodes[3], cgNodes[1], cgNodes[0], cgNodes[2] };
            }
            if (t == ElementType_t.TETRA_10) {
                return new int[] { cgNodes[4], cgNodes[6], cgNodes[5], cgNodes[7], cgNodes[9], cgNodes[8], cgNodes[3], cgNodes[1], cgNodes[0], cgNodes[2] };
            }
            // ---------------------
            if (t == ElementType_t.HEXA_8) {
                return new int[] { cgNodes[3], cgNodes[0], cgNodes[7], cgNodes[2], cgNodes[4], cgNodes[5], cgNodes[1], cgNodes[6] };
            }
            if (t == ElementType_t.HEXA_20 ) {
                return new int[] { cgNodes[15], cgNodes[14], cgNodes[18], cgNodes[10], cgNodes[12], cgNodes[13], cgNodes[16], cgNodes[8], cgNodes[17], cgNodes[19], cgNodes[9], cgNodes[11], cgNodes[3], cgNodes[0], cgNodes[7], cgNodes[2], cgNodes[4], cgNodes[5], cgNodes[1], cgNodes[6] };
            }
            if (t == ElementType_t.HEXA_27) {
                return new int[] { cgNodes[26], cgNodes[23], cgNodes[21], cgNodes[25], cgNodes[20], cgNodes[22], cgNodes[24], cgNodes[15], cgNodes[14], cgNodes[18], cgNodes[10], cgNodes[12], cgNodes[13], cgNodes[16], cgNodes[8], cgNodes[17], cgNodes[19], cgNodes[9], cgNodes[11], cgNodes[3], cgNodes[0], cgNodes[7], cgNodes[2], cgNodes[4], cgNodes[5], cgNodes[1], cgNodes[6] };
            }

            throw new NotSupportedException();
        }

        static MultidimensionalArray GetTransformationParams(ElementType_t t, int[] cgNodes, double[][] coords, int D) {
            var R = MultidimensionalArray.Create(cgNodes.Length, D);
            int[] NodePerm = GetBoSSSNodes(t, cgNodes);
            for (int l = 0; l < NodePerm.Length; l++) {
                for (int d = 0; d < D; d++) {
                    R[l, d] = coords[d][NodePerm[l]-1];
                }
            }
            return R;
        }
    }
}

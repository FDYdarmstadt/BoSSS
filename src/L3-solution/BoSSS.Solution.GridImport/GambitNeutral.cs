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
using System.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Solution.GridImport.Gambit {

    /// <summary>
    /// a gambit neutral file (at least the sections that are important for BoSSS)
    /// </summary>
    public class GambitNeutral : IGridImporter {

        /// <summary>
        /// 1st section which is read from gambit file;
        /// </summary>
        public string name_1;

        /// <summary>
        /// 
        /// </summary>
        public class PROGRAM {

            /// <summary>
            /// guess what ?
            /// </summary>
            /// <returns></returns>
            public PROGRAM Clone() {
                PROGRAM r = new PROGRAM();
                r.NBSETS = this.NBSETS;
                r.NDFCD_NoOfDims = this.NDFCD_NoOfDims;
                r.NDFVL = this.NDFVL;
                r.NELEM_NoOfElements = this.NELEM_NoOfElements;
                r.NGRPS = this.NGRPS;
                r.NUMNP_NoOfNodes = this.NUMNP_NoOfNodes;
                return r;
            }


            /// <summary>
            /// Number of nodes
            /// </summary>
            public int NUMNP_NoOfNodes;

            /// <summary>
            /// Number of elements
            /// </summary>
            public int NELEM_NoOfElements;

            /// <summary>
            /// 
            /// </summary>
            public int NGRPS;

            /// <summary>
            /// 
            /// </summary>
            public int NBSETS;

            /// <summary>
            /// number of dimensions
            /// </summary>
            public int NDFCD_NoOfDims;

            /// <summary>
            /// 
            /// </summary>
            public int NDFVL;
        }

        /// <summary>
        /// 2nd section which is read from gambit file;
        /// </summary>
        PROGRAM program_2;

        /// <summary>
        /// 
        /// </summary>
        public class NODAL_COORDINATES {


            /// <summary>
            /// calculates the range of node numbers (<see cref="nodeNumber"/>
            /// </summary>
            /// <param name="nodeNoMin"></param>
            /// <param name="nodeNoMax"></param>
            public void GetNodeNoRange(out int nodeNoMin, out int nodeNoMax) {
                nodeNoMin = int.MaxValue;
                nodeNoMax = int.MinValue;
                for (int i = 0; i < nodeNumber.Length; i++) {
                    int id = nodeNumber[i];
                    if (id > nodeNoMax)
                        nodeNoMax = id;
                    if (id < nodeNoMin)
                        nodeNoMin = id;
                }
            }


            /// <summary>
            /// computes a mapping from node numbers (<see cref="nodeNumber"/>) to
            /// node indices;
            /// </summary>
            /// <returns>
            /// 1st index: node number;
            /// 2nd index: always 0;
            /// note: this is a non-zero based array;
            /// </returns>
            public int[,] ComuteNodeNumber2IndexMap() {
                int idmin, idmax;
                GetNodeNoRange(out idmin, out idmax);

                int[,] id2index = (int[,])Array.CreateInstance(typeof(int), new int[] { idmax - idmin + 1, 1 }, new int[] { idmin, 0 });
                ;
                for (int i = id2index.GetLowerBound(0); i <= id2index.GetUpperBound(0); i++)
                    id2index[i, 0] = int.MinValue;
                for (int i = 0; i < nodeNumber.Length; i++) {
                    int id = nodeNumber[i];
                    if (id2index[id, 0] != int.MinValue)
                        throw new ApplicationException("node number " + id + " is used more than once.");
                    id2index[id, 0] = i;
                }

                return id2index;
            }

            /// <summary>
            /// node nuberd that identify the nodes (mapping: "node index" to "node number");
            /// index: node index;
            /// </summary>
            public int[] nodeNumber;

            /// <summary>
            /// points in space;
            /// 1st index: node index;
            /// 2nd index: spatial dimension
            /// </summary>
            public double[,] coord;



            /// <summary>
            /// euclidean distance between two node indices
            /// </summary>
            /// <param name="NodeIndexA">
            /// node index (!! not to be confused with node number !!) of 1st node
            /// </param>
            /// <param name="NodeIndexB">
            /// node index (!! not to be confused with node number !!) of 2nd node
            /// </param>
            /// <returns></returns>
            public double CalcDist(int NodeIndexA, int NodeIndexB) {

                double s = 0;
                for (int d = 0; d < coord.GetLength(1); d++) {
                    double _d = coord[NodeIndexA, d] - coord[NodeIndexB, d];
                    s += _d * _d;
                }
                return Math.Sqrt(s);
            }

        }

        /// <summary>
        /// 3rd section which is read from gambit file;
        /// </summary>
        NODAL_COORDINATES nodal_coordinates_3;


        /// <summary>
        /// gambit elemnet types;
        /// </summary>
        public enum ElementType {
            /// <summary> </summary>
            Edge = 1,
            /// <summary> </summary>
            Quadrilateral = 2,
            /// <summary> </summary>
            Triangle = 3,
            /// <summary> </summary>
            Brick = 4,
            /// <summary> </summary>
            Wedge = 5,
            /// <summary> </summary>
            Tetrahedron = 6,
            /// <summary> </summary>
            Pyramid = 7,
            /// <summary> </summary>
            Undefined = 0

        }


        /// <summary>
        /// 
        /// </summary>
        public class ELEMENTS_CELLS {

            /// <summary>
            /// calculates the range of element numbers (<see cref="ne"/>
            /// </summary>
            /// <param name="cellNoMin">on exit, the lowest used cell number</param>
            /// <param name="cellNoMax">on exit, the highest used cell number</param>
            public void GetElementsNoRange(out int cellNoMin, out int cellNoMax) {
                cellNoMin = int.MaxValue;
                cellNoMax = int.MinValue;
                for (int i = 0; i < ne.Length; i++) {
                    int id = ne[i];
                    if (id > cellNoMax)
                        cellNoMax = id;
                    if (id < cellNoMin)
                        cellNoMin = id;
                }
            }

            /// <summary>
            /// computes a mapping from element numbers (<see cref="ne"/>) to
            /// element indices;
            /// </summary>
            /// <returns>
            /// 1st index: element number;
            /// 2nd index: always 0;
            /// note: this is a non-zero based array;
            /// </returns>
            public int[,] ComuteElementNumber2IndexMap() {
                int idmin, idmax;
                GetElementsNoRange(out idmin, out idmax);

                int[,] id2index = (int[,])Array.CreateInstance(typeof(int), new int[] { idmax - idmin + 1, 1 }, new int[] { idmin, 0 });
                for (int i = id2index.GetLowerBound(0); i <= id2index.GetUpperBound(0); i++)
                    id2index[i, 0] = int.MinValue;
                for (int i = 0; i < ne.Length; i++) {
                    int id = ne[i];
                    if (id2index[id, 0] != int.MinValue)
                        throw new ApplicationException("element number number " + id + " is used more than once.");
                    id2index[id, 0] = i;
                }

                return id2index;
            }


            /// <summary>
            /// element identification
            /// index: element index;
            /// </summary>
            public int[] ne;

            /// <summary>
            /// element type
            /// index: element index;
            /// </summary>
            public ElementType[] NTYPE_ElementType;

            /// <summary>
            /// nodes that define an element
            /// 1st index: element index;
            /// 2nd index: node index
            /// </summary>
            public int[][] NODE;


            /// <summary>
            /// checks wether the content is suitable for BoSSS;
            /// throws an exception if it is not the case;
            /// </summary>
            public void ValidateForBoSSS() {
                int N = ne.Length;

                ElementType et = ElementType.Undefined;
                for (int n = 0; n < N; n++) {
                    if (NTYPE_ElementType[n] == ElementType.Edge)
                        continue; // ignoring edges

                    if (et == ElementType.Undefined)
                        et = NTYPE_ElementType[n];

                    if (et != NTYPE_ElementType[n])
                        throw new ApplicationException("BoSSS supports only one type of elements/cells");

                    if (NTYPE_ElementType[n] == ElementType.Pyramid)
                        throw new ApplicationException("BoSSS doesn't support Pyramid elements.");

                    if (NTYPE_ElementType[n] == ElementType.Wedge)
                        throw new ApplicationException("BoSSS doesn't support Wedge elements");

                    if (NTYPE_ElementType[n] == ElementType.Brick
                        && NODE[n].Length != 8)
                        throw new ApplicationException("BoSSS doesn't support brickes with other than 8 nodes.");

                    if (NTYPE_ElementType[n] == ElementType.Quadrilateral
                        && NODE[n].Length != 4)
                        throw new ApplicationException("BoSSS doesn't support quadliterals with other than 4 nodes.");

                    if (NTYPE_ElementType[n] == ElementType.Tetrahedron
                        && NODE[n].Length != 4)
                        throw new ApplicationException("BoSSS doesn't support tetrahedrons with other than 4 nodes.");

                    if (NTYPE_ElementType[n] == ElementType.Triangle
                        && NODE[n].Length != 3)
                        throw new ApplicationException("BoSSS doesn't support triangles with other than 3 nodes.");
                }
            }

        }

        /// <summary>
        /// 4th section that is read from from gambit file.
        /// </summary>
        ELEMENTS_CELLS elements_cells_4;

        /// <summary>
        /// 
        /// </summary>
        public class BOUNDARY_CONDITIONS {

            /// <summary>
            /// 
            /// </summary>
            public List<BCSet_Cell> bcsets;


            /// <summary>
            /// boundary condition set for cell edges (in contrast to bc.'s for nodes);
            /// </summary>
            public class BCSet_Cell {

                /// <summary>
                /// non-shallow copy
                /// </summary>
                /// <returns></returns>
                public BCSet_Cell Clone() {
                    BCSet_Cell ret_bcset = new BCSet_Cell();
                    ret_bcset.name = this.name;
                    ret_bcset.IBCODE = (int[])this.IBCODE.Clone();
                    ret_bcset.ELEM = (int[])this.ELEM.Clone();
                    ret_bcset.Values = (double[,])this.Values.Clone();
                    ret_bcset.FaceNumber = (byte[])this.FaceNumber.Clone();
                    ret_bcset.elementType = (ElementType[])this.elementType.Clone();

                    return ret_bcset;
                }

                /// <summary>
                /// name of boundary condition set
                /// </summary>
                public string name;

                /// <summary>
                /// optioanal boundary condition codes
                /// index: 0 to 4, at maximum;
                /// </summary>
                public int[] IBCODE;

                /// <summary>
                /// element/cell number
                /// index: bc. set element index
                /// </summary>
                public int[] ELEM;

                /// <summary>
                /// 
                /// </summary>
                public ElementType[] elementType;

                /// <summary>
                /// By calling <see cref="GambitNeutral.Rename"/>,
                /// these are transformed from the Gambit convention to the
                /// BoSSS convention.
                /// </summary>
                public byte[] FaceNumber;

                /// <summary>
                /// optional element/cell values;
                /// </summary>
                public double[,] Values;
            }

        }

        BOUNDARY_CONDITIONS boundary_conditions_5;




        string ReadUntil(string s, bool opt) {
            bool found = false;
            while (!found) {
                string l = ReadLine();
                if (l == null)
                    break;
                l = l.TrimStart(' ', '\t'); // delete leading whitespaces
                if (l.StartsWith(s, StringComparison.InvariantCulture)) {
                    return l;
                }
            }

            if (opt)
                return null;
            else
                throw new ApplicationException("expecting \"" + s + "\".");

        }


        char[] whitespaces = new char[] { ' ', '\t', '\n', '\r' };

        TextReader m_rd;
        int lineNo;


        string ReadLine() {
            lineNo++;
            return m_rd.ReadLine();
        }

        /// <summary>
        /// true if numberings, etc. are BoSSS-compatible;
        /// <see cref="Rename"/>;
        /// </summary>
        public bool BoSSSCompatible {
            get {
                return m_BoSSSCompatible;
            }
        }

        bool m_BoSSSCompatible = false;

        /// <summary>
        /// brings al numberings, etc. into a boss compatible - order
        /// </summary>
        private void Rename() {

            if (m_BoSSSCompatible)
                throw new ApplicationException("Renaming has allready been done.");



            // 1st: rename node id's 
            // ---------------------
            {
                int[,] id2index = nodal_coordinates_3.ComuteNodeNumber2IndexMap();
                for (int i = 0; i < nodal_coordinates_3.nodeNumber.Length; i++) {
                    nodal_coordinates_3.nodeNumber[i] = i;
                }


                for (int j = 0; j < elements_cells_4.ne.Length; j++) {
                    int[] nodes = elements_cells_4.NODE[j];

                    for (int k = 0; k < nodes.Length; k++)
                        nodes[k] = id2index[nodes[k], 0];

                }
            }

            // 2nd: rename cell id's
            // ---------------------
            {
                int[,] id2index = (int[,])elements_cells_4.ComuteElementNumber2IndexMap();
                for (int i = 0; i < elements_cells_4.ne.Length; i++) {
                    elements_cells_4.ne[i] = i;
                }

                foreach (BOUNDARY_CONDITIONS.BCSet_Cell bcset in boundary_conditions_5.bcsets) {
                    for (int j = 0; j < bcset.ELEM.Length; j++) {
                        bcset.ELEM[j] = id2index[bcset.ELEM[j], 0];
                    }

                }
            }

            // 3rd: reorder vertex numbering
            // -----------------------------
            {
                for (int i = 0; i < elements_cells_4.ne.Length; i++) {
                    int[] NODES = elements_cells_4.NODE[i];

                    if (elements_cells_4.NTYPE_ElementType[i] == ElementType.Quadrilateral) {
                        int b = NODES[2];
                        NODES[2] = NODES[3];
                        NODES[3] = b;
                    }

                    if (elements_cells_4.NTYPE_ElementType[i] == ElementType.Brick) {
                        int _0 = NODES[0];
                        int _1 = NODES[1];
                        int _2 = NODES[2];
                        int _3 = NODES[4];
                        int _4 = NODES[3];
                        int _5 = NODES[7];
                        int _6 = NODES[5];
                        int _7 = NODES[6];

                        NODES[0] = _0;
                        NODES[1] = _1;
                        NODES[2] = _2;
                        NODES[3] = _3;
                        NODES[4] = _4;
                        NODES[5] = _5;
                        NODES[6] = _6;
                        NODES[7] = _7;
                    }
                }
            }

            // 4th: edge numbering
            // -------------------
            {
                foreach (BOUNDARY_CONDITIONS.BCSet_Cell bcset in boundary_conditions_5.bcsets) {
                    for (int i = 0; i < bcset.ELEM.Length; i++) {

                        int iEdge_Gambit = bcset.FaceNumber[i];
                        int iEdge_BoSSS = -1;
                        if (bcset.elementType[i] == ElementType.Brick) {
                            switch (iEdge_Gambit) {
                                case 1:
                                    iEdge_BoSSS = (int)Cube.Edge.Bottom;
                                    break;
                                case 2:
                                    iEdge_BoSSS = (int)Cube.Edge.Right;
                                    break;
                                case 3:
                                    iEdge_BoSSS = (int)Cube.Edge.Top;
                                    break;
                                case 4:
                                    iEdge_BoSSS = (int)Cube.Edge.Left;
                                    break;
                                case 5:
                                    iEdge_BoSSS = (int)Cube.Edge.Back;
                                    break;
                                case 6:
                                    iEdge_BoSSS = (int)Cube.Edge.Front;
                                    break;
                                default:
                                    throw new ApplicationException("unknow face number for Gambit Brick: " + iEdge_Gambit + ";");
                            }
                        } else if (bcset.elementType[i] == ElementType.Quadrilateral) {
                            switch (iEdge_Gambit) {
                                case 1:
                                    iEdge_BoSSS = (int)Square.Edge.Bottom;
                                    break;
                                case 2:
                                    iEdge_BoSSS = (int)Square.Edge.Right;
                                    break;
                                case 3:
                                    iEdge_BoSSS = (int)Square.Edge.Top;
                                    break;
                                case 4:
                                    iEdge_BoSSS = (int)Square.Edge.Left;
                                    break;
                                default:
                                    throw new ApplicationException("unknow face number for Gambit Quadliteral: " + iEdge_Gambit + ";");
                            }
                        } else if (bcset.elementType[i] == ElementType.Triangle) {
                            switch (iEdge_Gambit) {
                                case 1:
                                    iEdge_BoSSS = 0;
                                    break;
                                case 2:
                                    iEdge_BoSSS = 1;
                                    break;
                                case 3:
                                    iEdge_BoSSS = 2;
                                    break;
                                default:
                                    throw new ApplicationException("unknow face number for Gambit Triangle: " + iEdge_Gambit + ";");
                            }
                        } else {
                            throw new NotImplementedException(bcset.elementType.ToString() + " will come soon.");
                        }


                        bcset.FaceNumber[i] = (byte)iEdge_BoSSS;
                    }

                }
            }
            m_BoSSSCompatible = true;
        }



        /// <summary>
        /// Cell Neighbours;
        /// 1st index: cell index;
        /// 2nd index; edge index, BoSSS convention;
        /// The content of this member is computed by <see cref="BuildNeighbourship"/>,
        /// because this info is not present in the ".neu"-file;
        /// </summary>
        int[][] m_NeighbourCells;


        /// <summary>
        /// cannot be called before <see cref="Rename"/>;
        /// finds cell neighbours;
        /// </summary>
        private void BuildNeighbourship() {
            if (!m_BoSSSCompatible)
                throw new ApplicationException("cannot be called before \"Rename\".");


            int J = elements_cells_4.ne.Length;
            m_NeighbourCells = new int[J][];

            // build a mapping from vertex index -> cell index
            List<int>[] points2cells = new List<int>[nodal_coordinates_3.nodeNumber.Length];
            for (int i = 0; i < points2cells.Length; i++)
                points2cells[i] = new List<int>();

            for (int j = 0; j < J; j++) {
                int[] nodes = elements_cells_4.NODE[j];
                for (int k = 0; k < nodes.Length; k++) {
                    points2cells[nodes[k]].Add(j);
                }
            }


            for (int j = 0; j < J; j++) {

                int[] neigh = null;
                int[] nodes = elements_cells_4.NODE[j];

                if (elements_cells_4.NTYPE_ElementType[j] == ElementType.Edge) {
                    // nothing to do
                } else if (elements_cells_4.NTYPE_ElementType[j] == ElementType.Brick) {
                    neigh = new int[6];

                    // Neighbour at Top (= between Vertex 7,5,2,4)
                    neigh[(int)Cube.Edge.Top] = NeighInd(j, points2cells[nodes[7]], points2cells[nodes[5]], points2cells[nodes[2]], points2cells[nodes[4]]);

                    // Neighbour at Bottom (= between Vertex 0,1,3,6)
                    neigh[(int)Cube.Edge.Bottom] = NeighInd(j, points2cells[nodes[0]], points2cells[nodes[1]], points2cells[nodes[3]], points2cells[nodes[6]]);

                    // Neighbour at Left (= between Vertex 0,2,3,7)
                    neigh[(int)Cube.Edge.Left] = NeighInd(j, points2cells[nodes[0]], points2cells[nodes[2]], points2cells[nodes[3]], points2cells[nodes[7]]);

                    // Neighbour at Right (= between Vertex 1,4,5,6)
                    neigh[(int)Cube.Edge.Right] = NeighInd(j, points2cells[nodes[1]], points2cells[nodes[4]], points2cells[nodes[5]], points2cells[nodes[6]]);

                    // Neighbour at Front (= between Vertex 3,6,5,7)
                    neigh[(int)Cube.Edge.Front] = NeighInd(j, points2cells[nodes[3]], points2cells[nodes[6]], points2cells[nodes[5]], points2cells[nodes[7]]);

                    // Neighbour at Back (= between Vertex 0,1,2,4)
                    neigh[(int)Cube.Edge.Back] = NeighInd(j, points2cells[nodes[0]], points2cells[nodes[1]], points2cells[nodes[2]], points2cells[nodes[4]]);



                } else if (elements_cells_4.NTYPE_ElementType[j] == ElementType.Quadrilateral) {
                    neigh = new int[4];

                    // Neighbour at Top (= between Vertex 2 and 3)
                    neigh[(int)Square.Edge.Top] = NeighInd(j, points2cells[nodes[2]], points2cells[nodes[3]]);

                    // Neighbour at Bottom (= between Vertex 0 and 1)
                    neigh[(int)Square.Edge.Bottom] = NeighInd(j, points2cells[nodes[0]], points2cells[nodes[1]]);

                    // Neighbour at Left (= between Vertex 0 and 2)
                    neigh[(int)Square.Edge.Left] = NeighInd(j, points2cells[nodes[0]], points2cells[nodes[2]]);

                    // Neighbour at Right (= between Vertex 1 and 3)
                    neigh[(int)Square.Edge.Right] = NeighInd(j, points2cells[nodes[1]], points2cells[nodes[3]]);

                } else if (elements_cells_4.NTYPE_ElementType[j] == ElementType.Tetrahedron) {
                    neigh = new int[4];
                    throw new NotImplementedException("Tetras will come soon.");
                } else if (elements_cells_4.NTYPE_ElementType[j] == ElementType.Triangle) {
                    neigh = new int[3];

                    // Neighbour at Edge 0 (=between Vertex 0 and 1):
                    neigh[0] = NeighInd(j, points2cells[nodes[0]], points2cells[nodes[1]]);

                    // Neighbour at Edge 1 (=between Vertex 1 and 2):
                    neigh[1] = NeighInd(j, points2cells[nodes[1]], points2cells[nodes[2]]);

                    // Neighbour at Edge 2 (=between Vertex 2 and 0):
                    neigh[2] = NeighInd(j, points2cells[nodes[2]], points2cells[nodes[0]]);
                } else {
                    throw new NotSupportedException(elements_cells_4.NTYPE_ElementType[j].ToString() + " is not supported by BoSSS.");
                }


                m_NeighbourCells[j] = neigh;

            }


            // track internal edges (like infinite thin walls)
            // -----------------------------------------------

            foreach (var bcset in this.boundary_conditions_5.bcsets) {

                int JJ = bcset.FaceNumber.Length;
                for (int e = 0; e < JJ; e++) {
                    int jCell = bcset.ELEM[e];
                    int edg = bcset.FaceNumber[e];

                    m_NeighbourCells[jCell][edg] = -1;
                }
            }
        }

        /// <summary>
        /// returns: 
        /// (intersection over all <paramref name="Grids"/>) without {<paramref name="TriaIndex"/>};
        /// this should at maximum one element, if the data structures are correct.
        /// </summary>
        /// <param name="TriaIndex">cell (index) who's neighbour we are searching</param>
        /// <param name="Grids">cells (their indices) which bound on each vertex of edge to investigate</param>
        /// <returns>
        /// index of neighbour cell, or negative number if specifyed triangle has no neighbour 
        /// at specifyed edge.
        /// </returns>
        int NeighInd(int TriaIndex, params List<int>[] Grids) {
            int ret = int.MinValue;
            int fndcnt = 0;
            foreach (int iTria in Grids[0]) {
                if (iTria == TriaIndex)
                    continue;

                bool found = true;

                for (int k = 1; found && k < Grids.Length; k++) {
                    if (!Grids[k].Contains(iTria)) {
                        found = false;
                    }
                }

                if (found) {
                    fndcnt++;
                    ret = iTria;
                }

            }
            if (fndcnt > 1)
                throw new ApplicationException("more than two neighbours found for one edge - data structure corrupted.");
            return ret;
        }

        /// <summary>
        /// loads a Gambit Neutral file (.neu - file) from a specified path
        /// </summary>
        /// <param name="FilePath"></param>
        public GambitNeutral(string FilePath) {

            m_rd = new StreamReader(FilePath);
            lineNo = 0;
            try {
                {
                    ReadLine();
                    string l = ReadLine();
                    if (!l.StartsWith("** GAMBIT NEUTRAL FILE", StringComparison.InvariantCulture))
                        throw new IOException("not a gambit neutral file: missing headnut;");
                    name_1 = ReadLine();
                }
                {
                    ReadUntil("PROGRAM:", false);
                    string l = ReadUntil("NUMNP", false);
                    string[] parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                    if (parts.Length != 6) {
                        throw new ApplicationException("expecting \"NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\".");
                    }
                    if (!parts[0].Equals("NUMNP", StringComparison.InvariantCulture)
                        || !parts[1].Equals("NELEM", StringComparison.InvariantCulture)
                        || !parts[2].Equals("NGRPS", StringComparison.InvariantCulture)
                        || !parts[3].Equals("NBSETS", StringComparison.InvariantCulture)
                        || !parts[4].Equals("NDFCD", StringComparison.InvariantCulture)
                        || !parts[5].Equals("NDFVL", StringComparison.InvariantCulture)) {
                        throw new ApplicationException("expecting \"NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\".");
                    }


                    l = ReadLine();
                    parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                    if (parts.Length != 6) {
                        throw new ApplicationException("expecting 6 numbers.");
                    }
                    program_2 = new PROGRAM();
                    program_2.NUMNP_NoOfNodes = int.Parse(parts[0]);
                    program_2.NELEM_NoOfElements = int.Parse(parts[1]);
                    program_2.NGRPS = int.Parse(parts[2]);
                    program_2.NBSETS = int.Parse(parts[3]);
                    program_2.NDFCD_NoOfDims = int.Parse(parts[4]);
                    program_2.NDFVL = int.Parse(parts[5]);

                    ReadUntil("ENDOFSECTION", false);
                }
                {
                    ReadUntil("NODAL COORDINATES", false);
                    nodal_coordinates_3 = new NODAL_COORDINATES();
                    nodal_coordinates_3.nodeNumber = new int[program_2.NUMNP_NoOfNodes];
                    nodal_coordinates_3.coord = new double[program_2.NUMNP_NoOfNodes, program_2.NDFCD_NoOfDims];

                    int N = program_2.NUMNP_NoOfNodes;
                    int D = program_2.NDFCD_NoOfDims;
                    for (int i = 0; i < N; i++) {
                        string l = ReadLine();
                        string[] parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                        if (parts.Length != D + 1)
                            throw new ApplicationException("expecting node number and " + D + " spatial coordinates.");
                        nodal_coordinates_3.nodeNumber[i] = int.Parse(parts[0]);
                        for (int d = 0; d < D; d++)
                            nodal_coordinates_3.coord[i, d] = double.Parse(parts[d + 1], System.Globalization.NumberFormatInfo.InvariantInfo);
                    }

                    ReadUntil("ENDOFSECTION", false);
                }
                {
                    ReadUntil("ELEMENTS/CELLS", false);
                    int N = program_2.NELEM_NoOfElements;
                    int D = program_2.NDFCD_NoOfDims;

                    elements_cells_4 = new ELEMENTS_CELLS();
                    elements_cells_4.ne = new int[N];
                    elements_cells_4.NODE = new int[N][];
                    elements_cells_4.NTYPE_ElementType = new ElementType[N];

                    for (int n = 0; n < N; n++) {
                        string l = ReadLine();
                        string[] parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);

                        elements_cells_4.ne[n] = int.Parse(parts[0]);
                        elements_cells_4.NTYPE_ElementType[n] = (ElementType)int.Parse(parts[1]);
                        int NN = int.Parse(parts[2]);

                        elements_cells_4.NODE[n] = new int[NN];
                        for (int nn = 0; nn < NN; nn++) {
                            elements_cells_4.NODE[n][nn] = int.Parse(parts[nn + 3]);
                        }
                    }

                    ReadUntil("ENDOFSECTION", false);
                }
                {
                    boundary_conditions_5 = new BOUNDARY_CONDITIONS();
                    boundary_conditions_5.bcsets = new List<BOUNDARY_CONDITIONS.BCSet_Cell>();
                    while (ReadUntil("BOUNDARY CONDITIONS", true) != null) {

                        //while (true) {
                        string l = ReadLine();
                        string[] parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                        if (parts[0].StartsWith("ENDOFSECTION", StringComparison.InvariantCulture))
                            continue;

                        if (parts.Length < 4)
                            throw new ApplicationException("expecting at least 4 entries.");

                        BOUNDARY_CONDITIONS.BCSet_Cell bcset;

                        int itype = int.Parse(parts[1]);
                        if (itype == 0)
                            throw new NotSupportedException("node boundary conditions are not supported in BoSSS.");
                        else if (itype == 1)
                            bcset = new BOUNDARY_CONDITIONS.BCSet_Cell();
                        else
                            throw new ApplicationException("unknown data type: " + itype + " (boundary condition type).");

                        boundary_conditions_5.bcsets.Add(bcset);
                        bcset.name = parts[0];

                        int Nentry = int.Parse(parts[2]);
                        int Nvalues = int.Parse(parts[3]);

                        bcset.IBCODE = new int[parts.Length - 4];
                        for (int i = 0; i < bcset.IBCODE.Length; i++) {
                            bcset.IBCODE[i] = int.Parse(parts[i + 4]);
                        }

                        bcset.ELEM = new int[Nentry];
                        bcset.FaceNumber = new byte[Nentry];
                        bcset.elementType = new ElementType[Nentry];
                        bcset.Values = new double[Nentry, Nvalues];

                        for (int i = 0; i < Nentry; i++) {
                            l = ReadLine();
                            parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);

                            bcset.ELEM[i] = int.Parse(parts[0]);
                            bcset.elementType[i] = (ElementType)int.Parse(parts[1]);
                            bcset.FaceNumber[i] = byte.Parse(parts[2]);

                            for (int j = 0; j < Nvalues; j++)
                                bcset.Values[i, j] = double.Parse(parts[3 + j]);
                        }

                        ReadUntil("ENDOFSECTION", false);
                    }
                }


            } catch (Exception e) {
                throw new ApplicationException("Error, line: " + lineNo + ": " + e.Message + ";", e);
            } finally {
                m_rd.Close();
            }
        }

        /// <summary>
        /// empty constructor; creates uninitialized element
        /// </summary>
        private GambitNeutral() {
        }


        /// <summary>
        /// converts cell' that are not supported in BoSSS into
        /// triangles/tetras;
        /// Cannot be called after <see cref="GenerateBoSSSGrid"/>
        /// has been executed;
        /// </summary>
        /// <returns></returns>
        public GambitNeutral ToLinearElements() {
            if (m_BoSSSCompatible)
                throw new ApplicationException("allready converted to BoSSS - format;");
            GambitNeutral ret = new GambitNeutral();

            int[,] NodeNo2NodeInd = this.nodal_coordinates_3.ComuteNodeNumber2IndexMap();


            ret.m_BoSSSCompatible = false;
            ret.name_1 = this.name_1;
            ret.program_2 = this.program_2.Clone();


            ret.nodal_coordinates_3 = this.nodal_coordinates_3;

            int Jold = this.elements_cells_4.ne.Length;

            // (non-bijective) mapping form "old" cell id's (index) to "new" cell id's (content);
            // "old" cell id: this.elements_cells_4.ne
            // "new" cell id: ret.elements_cells_4.ne
            int[,] oldId2newId;       // 1st index: cell number (id); 2nd index: face number
            byte[,] oldFaceIndex2new;  // 1st index: cell number (id); 2nd index: face number
            {
                int CellIndMin = int.MaxValue;
                int CellIndMax = int.MinValue;

                for (int i = 0; i < this.elements_cells_4.ne.Length; i++) {
                    CellIndMin = Math.Min(CellIndMin, this.elements_cells_4.ne[i]);
                    CellIndMax = Math.Max(CellIndMin, this.elements_cells_4.ne[i]);
                }

                oldId2newId = (int[,])Array.CreateInstance(typeof(int), new int[] { 1 + CellIndMax - CellIndMin, 6 }, new int[] { CellIndMin, 1 });
                oldFaceIndex2new = (byte[,])Array.CreateInstance(typeof(byte), new int[] { 1 + CellIndMax - CellIndMin, 6 }, new int[] { CellIndMin, 1 });

                for (int i = CellIndMin; i <= CellIndMax; i++)
                    for (int j = 1; j <= 6; j++) {
                        oldFaceIndex2new[i, j] = byte.MaxValue;
                        oldId2newId[i, j] = int.MinValue;
                    }
            }

            // cell conversion
            {
                ret.elements_cells_4 = new ELEMENTS_CELLS();
                var ret_ne = new List<int>(Jold);
                var ret_NODE = new List<int[]>(Jold);
                var ret_NType = new List<ElementType>(Jold);

                int j_new = 0;
                for (int j = 0; j < Jold; j++) {
                    ElementType et = this.elements_cells_4.NTYPE_ElementType[j];
                    int[] Nodes = this.elements_cells_4.NODE[j];

                    if (et == ElementType.Edge)
                        continue; // we do not care about this type
                    else if (et == ElementType.Triangle) {
                        if (Nodes.Length == 3) {
                            // element can be taken as it is

                            ret_ne.Add(j_new + 1);
                            oldId2newId[this.elements_cells_4.ne[j], 1] = ret_ne[j_new];
                            oldId2newId[this.elements_cells_4.ne[j], 2] = ret_ne[j_new];
                            oldId2newId[this.elements_cells_4.ne[j], 3] = ret_ne[j_new];
                            oldFaceIndex2new[this.elements_cells_4.ne[j], 1] = 1;
                            oldFaceIndex2new[this.elements_cells_4.ne[j], 2] = 2;
                            oldFaceIndex2new[this.elements_cells_4.ne[j], 3] = 3;
                            ret_NType.Add(ElementType.Triangle);
                            ret_NODE.Add(Nodes);

                            j_new++;
                            continue;
                        } else if (Nodes.Length == 6) {
                            // conversion needed

                            throw new NotImplementedException("will come somewhen.");

                        } else if (Nodes.Length == 7) {
                            // conversion needed

                            throw new NotImplementedException("will come somewhen.");

                        } else
                            throw new ApplicationException("unknown Gambit triangle type.");

                    } else if (et == ElementType.Quadrilateral) {
                        if (Nodes.Length == 4) {

                            // distance between vertices 0 and 2 (first diagonal)
                            double dist02 = this.nodal_coordinates_3.CalcDist(NodeNo2NodeInd[Nodes[0], 0], NodeNo2NodeInd[Nodes[2], 0]);

                            // distance between vertices 3 and 1 (second diagonal)
                            double dist13 = this.nodal_coordinates_3.CalcDist(NodeNo2NodeInd[Nodes[3], 0], NodeNo2NodeInd[Nodes[1], 0]);


                            if (dist13 < dist02) {
                                ret_ne.Add(j_new + 1);
                                ret_NType.Add(ElementType.Triangle);
                                int[] NewNodes = new int[] { Nodes[0], Nodes[1], Nodes[3] };
                                ret_NODE.Add(NewNodes);
                                j_new++;

                                ret_ne.Add(j_new + 1);
                                ret_NType.Add(ElementType.Triangle);
                                NewNodes = new int[] { Nodes[1], Nodes[2], Nodes[3] };
                                ret_NODE.Add(NewNodes);
                                j_new++;

                                oldFaceIndex2new[this.elements_cells_4.ne[j], 1] = 1;
                                oldFaceIndex2new[this.elements_cells_4.ne[j], 2] = 1;
                                oldFaceIndex2new[this.elements_cells_4.ne[j], 3] = 2;
                                oldFaceIndex2new[this.elements_cells_4.ne[j], 4] = 3;

                                oldId2newId[this.elements_cells_4.ne[j], 1] = ret_ne[j_new - 2];
                                oldId2newId[this.elements_cells_4.ne[j], 2] = ret_ne[j_new - 1];
                                oldId2newId[this.elements_cells_4.ne[j], 3] = ret_ne[j_new - 1];
                                oldId2newId[this.elements_cells_4.ne[j], 4] = ret_ne[j_new - 2];
                            } else {
                                ret_ne.Add(j_new + 1);
                                ret_NType.Add(ElementType.Triangle);
                                int[] NewNodes = new int[] { Nodes[0], Nodes[1], Nodes[2] };
                                ret_NODE.Add(NewNodes);
                                j_new++;

                                ret_ne.Add(j_new + 1);
                                ret_NType.Add(ElementType.Triangle);
                                NewNodes = new int[] { Nodes[0], Nodes[2], Nodes[3] };
                                ret_NODE.Add(NewNodes);
                                j_new++;

                                oldFaceIndex2new[this.elements_cells_4.ne[j], 1] = 1;
                                oldFaceIndex2new[this.elements_cells_4.ne[j], 2] = 2;
                                oldFaceIndex2new[this.elements_cells_4.ne[j], 3] = 2;
                                oldFaceIndex2new[this.elements_cells_4.ne[j], 4] = 3;

                                oldId2newId[this.elements_cells_4.ne[j], 1] = ret_ne[j_new - 2];
                                oldId2newId[this.elements_cells_4.ne[j], 2] = ret_ne[j_new - 2];
                                oldId2newId[this.elements_cells_4.ne[j], 3] = ret_ne[j_new - 1];
                                oldId2newId[this.elements_cells_4.ne[j], 4] = ret_ne[j_new - 1];
                            }

                            continue;

                        } else if (Nodes.Length == 8) {
                            throw new NotImplementedException("will come somewhen.");
                        } else if (Nodes.Length == 9) {
                            throw new NotImplementedException("will come somewhen.");
                        } else
                            throw new ApplicationException("unknown Gambit quad type.");
                    } else {
                        throw new NotImplementedException("will come somewhen.");
                    }
                }
                ret.elements_cells_4.ne = ret_ne.ToArray();
                ret.elements_cells_4.NODE = ret_NODE.ToArray();
                ret.elements_cells_4.NTYPE_ElementType = ret_NType.ToArray();
            }

            // conversion of boundary conditions
            {
                ret.boundary_conditions_5 = new BOUNDARY_CONDITIONS();
                ret.boundary_conditions_5.bcsets = new List<BOUNDARY_CONDITIONS.BCSet_Cell>();

                foreach (var this_bcset in this.boundary_conditions_5.bcsets) {
                    BOUNDARY_CONDITIONS.BCSet_Cell ret_bcset = this_bcset.Clone();
                    ret.boundary_conditions_5.bcsets.Add(ret_bcset);

                    int E = ret_bcset.ELEM.Length;
                    for (int e = 0; e < E; e++) {

                        ret_bcset.elementType[e] = ElementType.Triangle;

                        int cellNo = this_bcset.ELEM[e];
                        int FaceNo = this_bcset.FaceNumber[e];
                        ret_bcset.FaceNumber[e] = oldFaceIndex2new[cellNo, FaceNo];

                        ret_bcset.ELEM[e] = oldId2newId[cellNo, FaceNo];
                    }
                }
            }

            return ret;
        }


        /// <summary>
        /// Detects wether the grid contains elements which require a nonlinear
        /// mapping (to be exact anyting else but affine-linear) to their 
        /// (e.g. trapezoidal quads)
        /// reference domain.
        /// </summary>
        /// <returns>
        /// true, if all elements (cells) are linear (or affine-linear)
        /// </returns>
        public bool TestForLinearElements() {
            int J = elements_cells_4.NODE.Length;


            int[,] NodeNumbers2Index = null;

            int D = nodal_coordinates_3.coord.GetLength(1);
            Vector v1 = new Vector(3), v2 = new Vector(3);

            for (int j = 0; j < J; j++) {
                ElementType et = elements_cells_4.NTYPE_ElementType[j];
                int[] nodes = elements_cells_4.NODE[j];


                if (et == ElementType.Edge)
                    // the 3-vertices-edge could be nonlinear
                    // (we do not consider it)
                    continue;

                if (et == ElementType.Triangle)
                    // the 6- or 7-vertice triangle could be nonlinear
                    // (we do not consider it)
                    continue;

                if (et == ElementType.Quadrilateral) {
                    if (NodeNumbers2Index == null)
                        NodeNumbers2Index = nodal_coordinates_3.ComuteNodeNumber2IndexMap();

                    // we're checking wether the sides are parallel

                    v1.z = 0;
                    v2.z = 0;
                    for (int d = 0; d < D; d++) {
                        v1[d] = nodal_coordinates_3.coord[NodeNumbers2Index[nodes[3], 0], d] - nodal_coordinates_3.coord[NodeNumbers2Index[nodes[0], 0], d];
                        v2[d] = nodal_coordinates_3.coord[NodeNumbers2Index[nodes[2], 0], d] - nodal_coordinates_3.coord[NodeNumbers2Index[nodes[1], 0], d];
                    }

                    v2.Sub(v1);

                    if (v2.Abs() > 1.0e-6 * v1.Abs())
                        return false;

                }

            }

            return true;
        }



        /// <summary>
        /// tests wether this grid is suitable for BoSSS or not;
        /// </summary>
        /// <returns>
        /// true, if this grid needs to be converted 
        /// (by <see cref="ToLinearElements"/>)
        /// to be compatible with BoSSS;
        /// </returns>
        public bool BoSSSConversionNeccessary() {
            if (!TestForLinearElements())
                return true;

            try {
                elements_cells_4.ValidateForBoSSS();
            } catch (Exception) {
                return true;
            }

            return false;
        }

        #region IImporter Members



        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public GridCommons GenerateBoSSSGrid() {

            elements_cells_4.ValidateForBoSSS();
            Rename();
            BuildNeighbourship();

            int J = elements_cells_4.ne.Length;
            int D = program_2.NDFCD_NoOfDims;
            int N = elements_cells_4.NODE[0].Length;


            double[, ,] vertices = new double[J, N, D];
            long[] globalID = new long[J];
            byte[,] edgeTags = new byte[J, N];
            long[,] CellNeighbours = new long[J, N];


            // fill vertices
            // -------------
            for (int j = 0; j < J; j++) {
                for (int n = 0; n < N; n++) {
                    int nodeNo = elements_cells_4.NODE[j][n];
                    for (int d = 0; d < D; d++) {
                        vertices[j, n, d] = nodal_coordinates_3.coord[nodeNo, d];
                    }
                }
            }

            // fill globalID
            // -------------
            for (int j = 0; j < J; j++) {
                globalID[j] = j;
            }


            // fill neighbourship
            // ------------------
            for (int j = 0; j < J; j++) {
                for (int n = 0; n < N; n++) {
                    CellNeighbours[j, n] = m_NeighbourCells[j][n];

                }
            }

            // create grid object
            // ------------------

            GridCommons grd = null;
            if (elements_cells_4.NTYPE_ElementType[0] == ElementType.Quadrilateral)
                //grd = new Cartesian2DGrid();
                grd = new Grid2D(Square.Instance);
            else if (elements_cells_4.NTYPE_ElementType[0] == ElementType.Triangle)
                //grd = new UnstructuredTriangleGrid();
                grd = new Grid2D(Triangle.Instance);
            else if (elements_cells_4.NTYPE_ElementType[0] == ElementType.Brick)
                throw new NotImplementedException();
            else
                throw new NotSupportedException("currently, only Quad's, Trinagles and Bricks are supported.");


            // fill edgeTags
            // -------------

            int currEdgeTag = 1;
            foreach (BOUNDARY_CONDITIONS.BCSet_Cell bcset in boundary_conditions_5.bcsets) {
                if (currEdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                    throw new ApplicationException("not more than " + (GridCommons.FIRST_PERIODIC_BC_TAG - 1) + "edge sets are allowed.");
                }

                grd.EdgeTagNames.Add((byte)currEdgeTag, bcset.name);

                int JJ = bcset.ELEM.Length;

                for (int j = 0; j < JJ; j++)
                    edgeTags[bcset.ELEM[j], bcset.FaceNumber[j]] = (byte)currEdgeTag;


                currEdgeTag++;
            }

            // set variables
            // -------------
            grd.Description = this.name_1;
            //grd.GlobalID = globalID;
            //grd.CellNeighbours = CellNeighbours;
            //grd.Vertices = vertices;
            //grd.EdgeTags = edgeTags;

            return grd;
        }

        #endregion
    }
}

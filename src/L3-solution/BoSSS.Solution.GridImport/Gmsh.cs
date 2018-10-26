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
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using BoSSS.Foundation.Grid;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Solution.GridImport {

    /// <summary>
    /// Curved element types; To do: StandardDriver Elementypes_t should only contain non curved elements!!
    /// First the base nodes are used to create the curved element key nodes. 
    /// Second the base nodes are transformed to BoSSS node numbering. 
    /// The curved element nodes are *nor* transformed!
    /// </summary>
    public enum CurvedElementType_t {
        /// <summary>2-node line</summary>
        BAR_2 = 1,
        /// <summary>3-node triangle</summary>
        TRI_3 = 2,
        /// <summary>4-node quadrangle</summary>
        QUAD_4 = 3,
        /// <summary>4-node tetrahedron</summary>
        TETRA_4 = 4,
        /// <summary>8-node hexahedron</summary>
        HEXA_8 = 5,
        /// <summary>6-node prism</summary>
        PENTA_6 = 6,
        /// <summary>5-node pyramid</summary>
        PYRA_5 = 7,
        /// <summary>3-node second order line (2 nodes associated with the vertices and 1 with the edge)</summary>
        BAR_3 = 8,
        /// <summary>6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)</summary>
        TRI_6 = 9,
        /// <summary>9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)</summary>
        QUAD_9 = 10,
        /// <summary>10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)</summary>
        TETRA_10 = 11,
        /// <summary>27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)</summary>
        HEXA_27 = 12,
        /// <summary>18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)</summary>
        PENTA_18 = 13,
        /// <summary>14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)</summary>
        PYRA_14 = 14,
        /// <summary>1-node point</summary>
        NODE = 15,
        /// <summary>8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)</summary>
        QUAD_8 = 16,
        /// <summary>20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)</summary>
        HEXA_20 = 17,
        /// <summary>15-node second order prism (6 nodes associated with the vertices and 9 with the edges)</summary>
        PENTA_15 = 18,
        /// <summary>13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)</summary>
        PYRA_13 = 19,
        /// <summary>9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)</summary>
        TRI_9 = 20,
        /// <summary>10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)</summary>
        TRI_10 = 21,
        /// <summary>12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges) </summary>
        TRI_12 = 22,
        /// <summary>15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face) </summary>
        TRI_15 = 23,
        /// <summary>15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)</summary>
        TRI_152 = 24,
        /// <summary>21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)</summary>
        TRI_21 = 25,
        /// <summary>4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)</summary>
        BAR_4 = 26,
        /// <summary>5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)</summary>
        BAR_5 = 27,
        /// <summary>6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)</summary>
        BAR_6 = 28,
        /// <summary>20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)</summary>
        TETRA_20 = 29,
        /// <summary>35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)</summary>
        TETRA_35 = 30,
        /// <summary>56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)</summary>
        TETRA_56 = 31,
        /// <summary>16-node third order quadrangle (8 nodes associated with the vertices, 4 with the edges, 4 with the face)</summary>
        QUAD_16 = 36,
        /// <summary>25-node fourth order quadrangle (8 nodes associated with the vertices, 4 with the edges, 5 with the face)</summary>
        QUAD_25 = 37,
        /// <summary>36-node fourth order quadrangle (12 nodes associated with the vertices, 4 with the edges, 16 with the face)</summary>
        QUAD_36 = 38,
        /// <summary>17-node fourth order quadrangle (12 nodes associated with the vertices, 4 with the edges, 1 with the face)</summary>
        QUAD_172 = 39, // Achtung!!!! Unbekannte Nummerierung!!!!!!!!!!!
        /// <summary>64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)</summary>
        HEXA_64 = 92,
        /// <summary>125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)</summary>
        HEXA_125 = 93
    };

    /// <summary>
    /// Importer of msh format.
    ///  - First index is circularity, 0 clockwise, 1 counterclockwise
    ///  - second index is node reordering
    ///  - third index is element type.
    /// </summary>
    public class Gmsh : IGridImporter {

        static int[][,] REORDER = {
            // <summary>2-node line</summary>
            new int[,] {{ 0, 1 }},
            // <summary>3-node triangle</summary>
            new int[,] {{ 0, 1, 2 },{ 2, 1, 0 }},
            // <summary>4-node quadrangle</summary>
            new int[,] {{ 0, 1, 3, 2 },{ 0, 3, 1, 2 }},
            // <summary>4-node tetrahedron</summary>
            new int[,] {{ 3, 2, 0, 1 },{ 3, 2, 0, 1 }},
            // <summary>8-node hexahedron</summary>
            new int[,] {{ 2, 3, 1, 0, 6, 7, 5, 4 },{ 0, 1, 3, 2, 4, 5, 7, 6 }},
            // <summary>6-node prism</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>5-node pyramid</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>3-node second order line (2 nodes associated with the vertices and 1 with the edge)</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)</summary>
            new int[,] {{ 3, 4, 5, 0, 1, 2 },{ 5, 4, 3, 0, 2, 1 }},
            // <summary>9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)</summary>
            new int[,] {{ 8, 7, 5, 6, 4, 0, 1, 3, 2 },{ 8, 5, 7, 6, 4, 1, 0, 2, 3 }},
            // <summary>10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)</summary>
            new int[,] {{ 6, 4, 5, 7, 9, 8, 3, 2, 0, 1 },{ 6, 4, 5, 7, 9, 8, 3, 2, 0, 1 }},
            // <summary>27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)</summary> // Not implemented 
            new int[,] {{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 }}, 
            // <summary>18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>1-node point</summary> // Not BoSSS element
            new int[,] {{ 0 }}, 
            // <summary>8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)</summary>
            new int[,] {{ 7, 5, 6, 4, 0, 1, 3, 2 },{ 5, 7, 6, 4, 1, 0, 2, 3 }},
            // <summary>20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)</summary>  // Not implemented 
            new int[,] {{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 },{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }},
            // <summary>15-node second order prism (6 nodes associated with the vertices and 9 with the edges)</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)</summary> 
            new int[,] {{ 4, 3, 5, 6, 7, 8, 0, 1, 2 },{ 7, 8, 6, 5, 4, 3, 0, 2, 1 }},
            // <summary>10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)</summary>
            new int[,] {{ 9, 4, 3, 5, 6, 7, 8, 0, 1, 2 },{ 9, 7, 8, 6, 5, 4, 3, 0, 2, 1 }},
            // <summary>12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges) </summary>
            new int[,] {{ 5, 4, 3, 6, 7, 8, 9, 10, 11, 0, 1, 2 },{ 11, 10, 9, 8, 7, 6, 5, 3, 4, 0, 2, 1 }},  
            // <summary>15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face) </summary>
            new int[,] {{ 13, 14, 12, 5, 4, 3, 6, 7, 8, 9, 10, 11, 0, 1, 2 },{ 14, 13, 12, 9, 10, 11, 8, 7, 6, 5, 4, 3, 0, 2, 1 }}, 
            // <summary>15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)</summary>
            new int[,] {{ 6, 5, 4, 3, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2 },{ 11, 12, 13, 14, 10, 9, 8, 7, 6, 5, 4, 3, 0, 2, 1 }}, 
            // <summary>21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)</summary>
            new int[,] {{ 16, 19, 18, 17, 20, 15, 6, 5, 4, 3, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2 },{ 16, 19, 18, 17, 20, 15, 6, 5, 4, 3, 7, 8, 9, 10, 11, 12, 13, 14, 0, 1, 2 }},
            // <summary>4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)</summary> // Not BoSSS element
            new int[,] {{}},
            // <summary>20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)</summary>
            new int[,] {{ 16, 17, 18, 19, 8, 9, 5, 4, 6, 7, 11, 10, 15, 14, 13, 12, 3, 2, 0, 1 },{ 16, 17, 18, 19, 8, 9, 5, 4, 6, 7, 11, 10, 15, 14, 13, 12, 3, 2, 0, 1 }},
            // <summary>35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)</summary> // Not implemented 
            new int[,] {{ 34, 23, 24, 22, 26, 25, 27, 30, 28, 29, 32, 33, 31, 10, 11, 12, 6, 5, 4, 7, 8, 9, 15, 14, 13, 21, 20, 19, 18, 17, 16, 3, 2, 0, 1 },{ 34, 23, 24, 22, 26, 25, 27, 30, 28, 29, 32, 33, 31, 10, 11, 12, 6, 5, 4, 7, 8, 9, 15, 14, 13, 21, 20, 19, 18, 17, 16, 3, 2, 0, 1 }},
            // <summary>56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)</summary> // Not implemented 
            new int[,] {{}},
            // <summary>Unknown</summary> // Not implemented 
            new int[,] {{}},
            // <summary>Unknown</summary> // Not implemented 
            new int[,] {{}},
            // <summary>Unknown</summary> // Not implemented 
            new int[,] {{}},
            // <summary>Unknown</summary> // Not implemented 
            new int[,] {{}},
            // <summary>16-node third order quadrangle (8 nodes associated with the vertices, 4 with the edges, 4 with the face)</summary>
            new int[,] {{ 14, 13, 15, 12, 10, 11, 7, 6, 8, 9, 5, 4, 0, 1, 3, 2 },{ 15, 12, 14, 13, 7, 6, 10, 11, 9, 8, 4, 5, 1, 0, 2, 3 }},  
            // <summary>25-node third order quadrangle (8 nodes associated with the vertices, 4 with the edges, 9 with the face)</summary>
            new int[,] {{ 18, 21, 17, 22, 24, 20, 19, 23, 16, 13, 14, 15, 9, 8, 7, 10, 11, 12, 6, 5, 4, 0, 1, 3, 2 },{ 19, 23, 16, 22, 24, 20, 18, 21, 17, 9, 8, 7, 13, 14, 15, 12, 11, 10, 4, 5, 6, 1, 0, 2, 3 }},  
            // <summary>36-node third order quadrangle (12 nodes associated with the vertices, 4 with the edges, 16 with the face)</summary>
            new int[,] {{ 22, 27, 26, 21, 28, 34, 33, 25, 29, 35, 32, 24, 23, 30, 31, 20, 16, 17, 18, 19, 11, 10, 9, 8, 12, 13, 14, 15, 7, 6, 5, 4, 0, 1, 3, 2 },{ 23, 30, 31, 20, 29, 35, 32, 24, 28, 34, 33, 25, 22, 27, 26, 21, 11, 10, 9, 8, 16, 17, 18, 19, 15, 14, 13, 12, 4, 5, 6, 7, 1, 0, 2, 3 }},                                
        };

        private List<string> physicalNameLines;
        private CurvedElementType_t[] curved_element_types;
        private ElementType_t[] element_types;
        private bool[] element_is_curved;
        private static int number_of_coord_dimensions = 2;
        private int verticesPerCell;
        private ElementType_t bosss_element_type;
        private CurvedElementType_t curved_element_type;
        private int[] msh_element_index_to_bosss_element_index;
        private int[] bosss_element_index_to_msh_element_index;
        private PointSetType_t[] pointset_types;
        private GridLocation_t[] gridlocation_types;
        private BCType_t[] bc_types;
        private int[][] bcelements; // Either elements or points! bcelements[bc_index_with_index_mapping][element_index/point_index] = element_index
        private int[] bcindex_of_element; // Either elements or points! bcindex_of_element[element_index] = bc_index_with_index_mapping
        private int bosss_elements;
        private bool binary = true;
        private bool hasbc;
        private TextReader m_rd;
        private BinaryReader m_rdb;
        private FileStream m_rdf;
        private int lineNo;
        private char[] whitespaces = new char[] { ' ', '\t', '\n', '\r' };
        private double[] x;
        private double[] y;
        private double[] z;
        private double[] ref_x;
        private double[] ref_y;
        private double[] ref_z;
        // Important: Do not directly call elements[index][node_index] = number, do always call elements[index][node_index] = node_index_mapping[number]
        // Because: Note that the element-numbers do not necessarily have to form a dense nor an ordered sequence. 
        private int[][] elements;
        private int[][] base_elements;
        private int[] ref_element;
        private int[,] permutation_array;
        // Important: Do not directly call boconames[index], do always call boconames[physical_index_mapping[[index]]
        // Because: Note that the physical-numbers do not necessarily have to form a dense nor an ordered sequence. 
        private String[] boconames;
        private int[] node_index_mapping;
        private int[] physical_index_mapping;

        private int[] chirality_index;

        private static char[] ca = System.Environment.NewLine.ToCharArray();
        private static char[] trimChars = new char[3];

        /// <summary>
        /// static constructor
        /// </summary>
        static Gmsh() {
            int i = 0;
            for (; i < ca.Length; i++) {
                trimChars[i] = ca[i];
            }
            trimChars[i] = ' ';
        }

        /// <summary>
        /// loads an msh file either binary or ASCII from or to a specified path
        /// </summary>
        /// <param name="filePath"></param>
        public Gmsh(string filePath) {
            ParseBinary(filePath);
        }

        /// <summary>
        /// loads an msh file ASCII from or to a specified path
        /// </summary>
        /// <param name="filePath"></param>
        private void ParseAscii(string filePath) {
            if (m_rdb != null) {
                m_rdb.Close();
            }
            lineNo = 0;
            int i, j;
            using (m_rd = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))) {
                string l = ReadLine();
                if (!l.StartsWith("$MeshFormat", StringComparison.InvariantCulture)) {
                    throw new IOException("not a msh file: missing header");
                }
                l = ReadLine();
                string[] parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                if (parts[0] != "2.2") {
                    throw new ApplicationException("expecting .msh version 2.2");
                }
                binary = int.Parse(parts[1]) == 1;
                int datasize = int.Parse(parts[2]);
                if (datasize != 8) {
                    throw new ApplicationException("currently only data-size = sizeof(double) is supported");
                }

                l = ReadLine();
                if (!l.StartsWith("$EndMeshFormat", StringComparison.InvariantCulture)) {
                    throw new IOException("Missing end of header");
                }

                int array_length;
                int[] temp;
                int max;

                // Read optional physical names that define boundary condition
                // tags, but process them later
                l = ReadLine();
                physicalNameLines = new List<string>();
                if (l.TrimStart(' ', '\t').StartsWith("$PhysicalNames", StringComparison.InvariantCulture)) {
                    hasbc = true;

                    l = ReadLine();
                    parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                    array_length = int.Parse(parts[0]);

                    for (i = 0; i < array_length; i++) {
                        l = ReadLine();
                        physicalNameLines.Add(l);
                    }
                } else {
                    hasbc = false;
                }

                // The mesh nodes
                if (!l.TrimStart(' ', '\t').StartsWith("$Nodes")) {
                    ReadUntil("$Nodes", false);
                }
                l = ReadLine();
                parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length != 1) {
                    throw new ApplicationException("expecting only one node number");
                }
                array_length = int.Parse(parts[0]);

                x = new double[array_length];
                y = new double[array_length];
                z = new double[array_length];

                temp = new int[array_length];
                max = 0;

                for (i = 0; i < array_length; i++) {
                    l = ReadLine();
                    parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                    if (parts.Length != 4) {
                        throw new ApplicationException("expecting 3 coordinates per node");
                    }
                    int tempi = temp[i] = int.Parse(parts[0]);
                    if (max < tempi) {
                        max = temp[i];
                    }
                    // Important: Do not directly call x,y,z[index], do always call e.g. x[node_index_mapping[[index]]
                    // Because: Note that the node-numbers do not necessarily have to form a dense nor an ordered sequence. 
                    x[i] = double.Parse(parts[1], System.Globalization.NumberFormatInfo.InvariantInfo);
                    y[i] = double.Parse(parts[2], System.Globalization.NumberFormatInfo.InvariantInfo);
                    z[i] = double.Parse(parts[3], System.Globalization.NumberFormatInfo.InvariantInfo);
                }
                max++;
                node_index_mapping = new int[max];
                for (i = 0; i < max; i++) {
                    for (j = 0; j < array_length; j++) {
                        if (temp[j] == i) {
                            node_index_mapping[i] = j;
                            continue;
                        }
                    }
                }

                // Finally, the mesh elements
                ReadUntil("$Elements", false);
                l = ReadLine();
                parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length != 1) {
                    throw new ApplicationException("expecting only one element number");
                }
                array_length = int.Parse(parts[0]);

                elements = new int[array_length][];
                chirality_index = new int[array_length];
                base_elements = new int[array_length][];
                element_is_curved = new bool[array_length];
                curved_element_types = new CurvedElementType_t[array_length];
                element_types = new ElementType_t[array_length];
                bcindex_of_element = new int[array_length];

                for (i = 0; i < array_length; i++) {
                    l = ReadLine();
                    parts = l.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);

                    int element_type = int.Parse(parts[1]);

                    curved_element_types[i] = (CurvedElementType_t)Enum.ToObject(typeof(CurvedElementType_t), element_type);

                    element_types[i] = convertCurvedElementTypeToStandardElementType(curved_element_types[i]);

                    element_is_curved[i] = !element_types[i].Equals(curved_element_types[i]);

                    if (!BoSSSCompliantElementType((ElementType_t)Enum.ToObject(typeof(ElementType_t), element_types[i]))) {
                        continue;
                    }

                    int number_of_tags = int.Parse(parts[2]);

                    if (number_of_tags >= 1) {
                        bcindex_of_element[i] = int.Parse(parts[3]);
                    }

                    int offset = number_of_tags + 3;

                    int array_length_inner = parts.Length - offset;
                    int c = NumberOfNodesPerElement(element_types[i]);

                    elements[i] = new int[array_length_inner];
                    base_elements[i] = new int[c];

                    for (j = 0; j < array_length_inner; j++) {
                        elements[i][j] = node_index_mapping[int.Parse(parts[j + offset])];
                    }
                    if (Is3D(element_types[i])) {
                        number_of_coord_dimensions = 3;
                    }
                    if (verticesPerCell < c) {
                        verticesPerCell = c;
                        bosss_element_type = element_types[i];
                        curved_element_type = curved_element_types[i];
                    }
                }
                parseBoundaryConditions();
            }
        }

        /// <summary>
        /// loads an msh file binary from or to a specified path
        /// </summary>
        /// <param name="filePath"></param>
        private void ParseBinary(string filePath) {
            lineNo = 0;
            int i, j, k;
            using (m_rdf = File.Open(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
            using (m_rdb = new BinaryReader(m_rdf, Encoding.ASCII)) {
                string l = ReadLine();
                if (!l.StartsWith("$MeshFormat", StringComparison.InvariantCulture)) {
                    throw new IOException("not a msh file: missing headnut");
                }
                if (ReadString() != "2.2") {
                    throw new ApplicationException("expecting .msh version 2.2");
                }
                binary = ReadString() == "1";
                if (!binary) {
                    ParseAscii(filePath);
                    return;
                }
                if (ReadAsciiInt() != 8) {
                    throw new ApplicationException("currently only data-size = sizeof(double) is supported");
                }
                int array_length;
                // Read optional physical names that define boundary condition
                // tags, but process them later
                l = ReadLine();
                physicalNameLines = new List<string>();
                if (l.TrimStart(' ', '\t').StartsWith("$PhysicalNames", StringComparison.InvariantCulture)) {
                    hasbc = true;

                    array_length = ReadAsciiInt();

                    for (i = 0; i < array_length; i++) {
                        l = ReadLine();
                        physicalNameLines.Add(l);
                    }
                } else {
                    hasbc = false;
                }

                ReadUntil("$Nodes", false);
                array_length = ReadAsciiInt();

                x = new double[array_length];
                y = new double[array_length];
                z = new double[array_length];

                int[] temp = new int[array_length];
                int max = 0;

                for (i = 0; i < array_length; i++) {
                    int tempi = temp[i] = ReadInt();
                    if (max < tempi) {
                        max = temp[i];
                    }
                    // Important: Do not direktly call x,y,z[index], do always call e.g. x[node_index_mapping[[index]]
                    // Because: Note that the node-numbers do not necessarily have to form a dense nor an ordered sequence. 
                    x[i] = ReadDouble();
                    y[i] = ReadDouble();
                    z[i] = ReadDouble();
                }
                max++;
                node_index_mapping = new int[max];
                for (i = 0; i < max; i++) {
                    for (j = 0; j < array_length; j++) {
                        if (temp[j] == i) {
                            node_index_mapping[i] = j;
                            continue;
                        }
                    }
                }
                ReadUntil("$Elements", false);
                array_length = ReadAsciiInt();

                elements = new int[array_length][];
                chirality_index = new int[array_length];
                base_elements = new int[array_length][];
                element_is_curved = new bool[array_length];
                curved_element_types = new CurvedElementType_t[array_length];
                element_types = new ElementType_t[array_length];
                bcindex_of_element = new int[array_length];

                for (i = 0; i < array_length; ) {
                    int element_type = ReadInt();
                    int num_elm_follow = ReadInt();
                    int number_of_tags = ReadInt();

                    for (k = 0; k < num_elm_follow; i++, k++) {
                        int tempi = ReadInt();

                        curved_element_types[i] = (CurvedElementType_t)Enum.ToObject(typeof(CurvedElementType_t), element_type);

                        element_types[i] = convertCurvedElementTypeToStandardElementType(curved_element_types[i]);

                        element_is_curved[i] = !element_types[i].Equals(curved_element_types[i]);

                        if (!BoSSSCompliantElementType((ElementType_t)Enum.ToObject(typeof(ElementType_t), element_types[i]))) {
                            continue;
                        }

                        if (number_of_tags >= 1) {
                            bcindex_of_element[i] = ReadInt();
                            number_of_tags--;
                            while (number_of_tags > 0) {
                                ReadInt();
                                number_of_tags--;
                            }
                        }

                        int array_length_inner = NumberOfNodesPerElement(curved_element_types[i]);
                        int c = NumberOfNodesPerElement(element_types[i]);

                        elements[i] = new int[array_length_inner];
                        base_elements[i] = new int[c];

                        for (j = 0; j < array_length_inner; j++) {
                            elements[i][j] = node_index_mapping[ReadInt()];
                        }
                        if (Is3D(element_types[i])) {
                            number_of_coord_dimensions = 3;
                        }
                        if (verticesPerCell < c) {
                            verticesPerCell = c;
                            bosss_element_type = element_types[i];
                            curved_element_type = curved_element_types[i];
                        }
                    }
                }
                parseBoundaryConditions();
            }
        }

        private void parseBoundaryConditions() {
            int i, j;
            msh_element_index_to_bosss_element_index = new int[elements.Length];
            bosss_element_index_to_msh_element_index = new int[elements.Length];
            for (i = 0; i < elements.Length; i++) {
                if (ChoosenBoSSSCompliantElementType(element_types[i])) {
                    bosss_element_index_to_msh_element_index[msh_element_index_to_bosss_element_index[i] = bosss_elements++] = i;
                }
            }
            if (!hasbc) {
                //Console.WriteLine("No BoundaryConditions");
                return;
            }
            for (i = 0; i < physicalNameLines.Count; i++) {
                string[] parts = physicalNameLines[i].Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                int dimension = int.Parse(parts[0]);
                if (dimension != number_of_coord_dimensions - 1) {
                    // Boundary conditions are always defined on
                    // D-1-dimensional entities; ignore others
                    physicalNameLines.RemoveAt(i);
                    continue;
                }
            }

            int array_length = physicalNameLines.Count;

            gridlocation_types = new GridLocation_t[array_length];
            pointset_types = new PointSetType_t[array_length];
            bc_types = new BCType_t[array_length];
            boconames = new String[array_length];
            bcelements = new int[array_length][];
            int[] temp = new int[array_length];
            int max = 0;

            for (i = 0; i < array_length; i++) {
                string[] parts = physicalNameLines[i].Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);

                int dimension = int.Parse(parts[0]); // What the f* is physical_dimension?

                int tempi = temp[i] = int.Parse(parts[1]);
                if (max < tempi) {
                    max = temp[i];
                }

                boconames[i] = parts[2];
                gridlocation_types[i] = GridLocation_t.CellCenter;
                pointset_types[i] = PointSetType_t.ElementList;
            }
            physical_index_mapping = new int[max + 1];
            for (i = 0; i < max + 1; i++) {
                for (j = 0; j < array_length; j++) {
                    if (temp[j] == i) {
                        physical_index_mapping[i] = j + 1;
                        continue;
                    }
                }
            }

            // Create actual boundary elements
            for (j = 0; j < bcelements.Length; j++) {
                bcelements[j] = new int[bcindex_of_element.Length];
            }

            i = 0;
            foreach (var z in boconames) {
                for (j = 0; j < bcindex_of_element.Length; j++) {
                    if (bcindex_of_element[j] <= max
                        && physical_index_mapping[bcindex_of_element[j]] == i) {
                        bcelements[i][j] = j;
                    }
                }
                i++;
            }
        }

        private string ReadLine() {
            if (!binary) {
                lineNo++;
                return m_rd.ReadLine();
            }
            long position = m_rdf.Position;
            int c;
            try {
                do {
                    c = m_rdb.ReadByte();
                    if (ca.Length == 2 && c == ca[0]) {
                        c = m_rdb.ReadByte();
                        if (c == ca[1]) {
                            break;
                        }
                    }
                    if (IsSpaceChar((char)c)) {
                        break;
                    }
                } while (c != -1);
                if (c == -1) {
                    return null;
                }
            } catch (EndOfStreamException) {
                return null;
            }
            long newPosition = m_rdf.Position;
            m_rdf.Position = position;
            string ret = System.Text.ASCIIEncoding.ASCII.GetString(m_rdb.ReadBytes((int)(newPosition - position))).Trim(trimChars);
            m_rdf.Position = newPosition;
            return ret;
        }

        private int Read(char[] buffer, int i, int j) {
            return m_rd.Read(buffer, i, j);
        }

        private int ReadAsciiInt() {
            string s = ReadString();
            return Convert.ToInt32(s);
        }

        private int ReadInt() {
            return m_rdb.ReadInt32();
        }

        private double ReadDouble() {
            return m_rdb.ReadDouble();
        }

        private string ReadString() {
            long position = m_rdf.Position;
            int c;
            try {
                do {
                    c = m_rdb.ReadByte();
                    if (ca.Length == 2 && c == ca[0]) {
                        c = m_rdb.ReadByte();
                        if (c == ca[1]) {
                            break;
                        }
                    }
                    if (IsSpaceChar((char)c)) {
                        break;
                    }
                } while (c != -1);
                if (c == -1) {
                    return null;
                }
            } catch (EndOfStreamException) {
                return null;
            }
            long newPosition = m_rdf.Position;
            m_rdf.Position = position;
            string ret = System.Text.ASCIIEncoding.ASCII.GetString(m_rdb.ReadBytes((int)(newPosition - position))).Trim(trimChars);
            m_rdf.Position = newPosition;
            return ret;
        }

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

        /// <summary>
        /// Is it a space?
        /// </summary>
        /// <param name="c"></param>
        private static bool IsSpaceChar(char c) {
            for (int i = 0; i < ca.Length; i++) {
                if (c == ca[i]) {
                    return true;
                }
            }
            if (c == ' ') {
                return true;
            }
            return false;
        }

        /// <summary>
        /// Converts an curved element type to an standard element type
        /// </summary>
        /// <param name="elm"></param>
        private static ElementType_t convertCurvedElementTypeToStandardElementType(CurvedElementType_t elm) {
            String s = elm.ToString();

            ElementType_t relm;

            int index = s.IndexOf("_");

            if (index > -1) {
                s = s.Substring(0, index);
            }

            if (s.Equals("NODE")) {
                relm = (ElementType_t)Enum.Parse(typeof(ElementType_t), s);
            } else if (s.Equals("BAR")) {
                relm = (ElementType_t)Enum.Parse(typeof(ElementType_t), s + "_2");
            } else if (s.Equals("TRI")) {
                relm = (ElementType_t)Enum.Parse(typeof(ElementType_t), s + "_3");
            } else if (s.Equals("QUAD")) {
                relm = (ElementType_t)Enum.Parse(typeof(ElementType_t), s + "_4");
            } else if (s.Equals("TETRA")) {
                relm = (ElementType_t)Enum.Parse(typeof(ElementType_t), s + "_4");
            } else if (s.Equals("PYRA")) {
                relm = (ElementType_t)Enum.Parse(typeof(ElementType_t), s + "_5");
            } else if (s.Equals("PENTA")) {
                relm = (ElementType_t)Enum.Parse(typeof(ElementType_t), s + "_5");
            } else {
                relm = (ElementType_t)Enum.Parse(typeof(ElementType_t), s + "_8");
            }
            return relm;
        }

        /// <summary>
        /// Generates parameters for the transformation from reference to physical domain.
        /// All curved elements key nodes are in Gmesh numbering order.
        /// </summary>
        private MultidimensionalArray GenerateModalTransformationParams(int numberOfNodesPerCurvedElement, int l) {
            MultidimensionalArray transformationParams = MultidimensionalArray.Create(number_of_coord_dimensions, numberOfNodesPerCurvedElement);
            int numberOfCoefficientsAndNodes = NumberOfNodesPerElement(curved_element_type);

            int j = bosss_element_index_to_msh_element_index[l];

            double[,] coords = new double[number_of_coord_dimensions, numberOfCoefficientsAndNodes];
            double[,] coordsRef = new double[number_of_coord_dimensions, numberOfCoefficientsAndNodes];

            for (int d = 0; d < number_of_coord_dimensions; d++) {
                for (int i = 0; i < numberOfCoefficientsAndNodes; i++) {
                    int node_no = elements[j][i];

                    if (d == 0) {
                        coords[d, i] = x[node_no];
                        coordsRef[d, i] = ref_x[i];
                    } else if (d == 1) {
                        coords[d, i] = y[node_no];
                        coordsRef[d, i] = ref_y[i];
                    } else {
                        coords[d, i] = z[node_no];
                        coordsRef[d, i] = ref_z[i];
                    }
                }
            }

            for (int d = 0; d < number_of_coord_dimensions; d++) {
                double[,] source = new double[numberOfCoefficientsAndNodes, numberOfCoefficientsAndNodes];

                for (int i = 0; i < numberOfCoefficientsAndNodes; i++) {
                    generatesAndLoadsModalBaseFuncs(source, coordsRef, i);
                }

                MultidimensionalArray m = MultidimensionalArray.Create(numberOfCoefficientsAndNodes, numberOfCoefficientsAndNodes);

                m.Set2DArray(source);
                m = m.GetInverse();
                loadsModalTransformationParams(transformationParams, m, coords, d);
            }
            return transformationParams;
        }

        /// <summary>
        /// Loads the base functions.
        /// All curved elements key nodes are in Gmesh numbering order.
        /// </summary>
        /// <param name="TransformationParams"></param>
        /// <param name="m"></param>
        /// <param name="coords"></param>
        /// <param name="d"></param>
        private void loadsModalTransformationParams(MultidimensionalArray TransformationParams, MultidimensionalArray m, double[,] coords, int d) {
            int numberOfNodes = coords.Length / number_of_coord_dimensions;

            for (int i = 0; i < numberOfNodes; i++) {
                double sum = 0;
                for (int j = 0; j < numberOfNodes; j++) {
                    sum += m[i, j] * coords[d, j];
                }
                TransformationParams[d, i] = sum;
            }
        }

        /// <summary>
        /// Generates and loads the base functions.
        /// All curved elements key nodes are in Gmesh numbering order.
        /// </summary>
        /// <param name="target"></param>
        /// <param name="coordsRef"></param>
        /// <param name="i"></param>
        private void generatesAndLoadsModalBaseFuncs(double[,] target, double[,] coordsRef, int i) {
            int numberOfNodes = coordsRef.Length / number_of_coord_dimensions;
            Debug.Assert(numberOfNodes == 4 || numberOfNodes == 8);
            if (numberOfNodes == 4) {
                for (int j = 0; j < numberOfNodes; j++) {
                    if (j == 0) {
                        target[i, j] = 1;
                    } else if (j == 1) {
                        target[i, j] = coordsRef[0, i];
                    } else if (j == 2) {
                        target[i, j] = coordsRef[1, i];
                    } else {
                        target[i, j] = coordsRef[0, i] * coordsRef[1, i];
                    }
                }
            } else {
                for (int j = 0; j < numberOfNodes; j++) {
                    if (j == 0) {
                        target[i, j] = 1;
                    } else if (j == 1) {
                        target[i, j] = coordsRef[0, i];
                    } else if (j == 2) {
                        target[i, j] = coordsRef[1, i];
                    } else if (j == 3) {
                        target[i, j] = coordsRef[2, i];
                    } else if (j == 4) {
                        target[i, j] = coordsRef[0, i] * coordsRef[0, i];
                    } else if (j == 5) {
                        target[i, j] = coordsRef[1, i] * coordsRef[1, i];
                    } else if (j == 6) {
                        target[i, j] = coordsRef[2, i] * coordsRef[2, i];
                    } else {
                        target[i, j] = coordsRef[0, i] * coordsRef[1, i] * coordsRef[2, i];
                    }
                }
            }
        }

        /// <summary>
        /// Check if shuld be a down strip of transformations params
        /// </summary>
        /// <param name="numberOfNodesPerCurvedElement"></param>
        /// <param name="l"></param>
        private bool CheckIfElementIsNonLinear(int numberOfNodesPerCurvedElement, int l) {
            var trafops = GenerateModalTransformationParams(numberOfNodesPerCurvedElement, l);

            int j, length;
            bool ret = false;

            double max = 0;
            double temp;

            for (j = numberOfNodesPerCurvedElement - 1, length = numberOfNodesPerCurvedElement; j >= 0; j--) {
                if (number_of_coord_dimensions == 3) {
                    if ((temp = Math.Abs(trafops[0, j]) + Math.Abs(trafops[1, j]) + Math.Abs(trafops[2, j])) > max) {
                        max = temp;
                    }
                } else {
                    if ((temp = Math.Abs(trafops[0, j]) + Math.Abs(trafops[1, j])) > max) {
                        max = temp;
                        ;
                    }
                }
            }

            double dist = 0.00000001 / max;

            for (j = numberOfNodesPerCurvedElement - 1, length = numberOfNodesPerCurvedElement; j >= 0; j--) {
                if (number_of_coord_dimensions == 3) {
                    if (Math.Abs(trafops[0, j]) > dist || Math.Abs(trafops[1, j]) > dist || Math.Abs(trafops[2, j]) > dist) {
                        break;
                    }
                } else {
                    if (Math.Abs(trafops[0, j]) > dist || Math.Abs(trafops[1, j]) > dist) {
                        break;
                    }
                }
                length--;
            }
            if (length == numberOfNodesPerCurvedElement) {
                ret = true;
            }
            return ret;
        }

        private CellType ConvertToCellType(CurvedElementType_t c, bool nonLinear) {
            switch (c) {
                case CurvedElementType_t.NODE:
                    return CellType.Point;

                case CurvedElementType_t.BAR_2:
                    return CellType.Line_2;

                case CurvedElementType_t.BAR_3:
                    return CellType.Line_3;

                case CurvedElementType_t.BAR_4:
                    return CellType.Line_4;

                case CurvedElementType_t.BAR_5:
                    return CellType.Line_5;

                case CurvedElementType_t.BAR_6:
                    return CellType.Line_6;

                case CurvedElementType_t.TRI_3:
                    return CellType.Triangle_3;

                case CurvedElementType_t.TRI_6:
                    return CellType.Triangle_6;

                case CurvedElementType_t.TRI_9:
                    return CellType.Triangle_9;

                case CurvedElementType_t.TRI_10:
                    return CellType.Triangle_10;

                case CurvedElementType_t.TRI_12:
                    return CellType.Triangle_12;

                case CurvedElementType_t.TRI_15:
                    return CellType.Triangle_15;

                case CurvedElementType_t.TRI_152:
                    return CellType.Triangle_152;

                case CurvedElementType_t.TRI_21:
                    return CellType.Triangle_21;

                case CurvedElementType_t.TETRA_4:
                    return CellType.Tetra_Linear;

                case CurvedElementType_t.TETRA_10:
                    return CellType.Tetra_10;

                case CurvedElementType_t.TETRA_20:
                    return CellType.Tetra_35;

                case CurvedElementType_t.TETRA_35:
                    return CellType.Tetra_35;

                case CurvedElementType_t.TETRA_56:
                    return CellType.Tetra_56;

                case CurvedElementType_t.QUAD_4:
                    if (nonLinear) {
                        return CellType.Square_4;
                    }
                    return CellType.Square_Linear;

                case CurvedElementType_t.QUAD_8:
                    return CellType.Square_8;

                case CurvedElementType_t.QUAD_9:
                    return CellType.Square_9;

                case CurvedElementType_t.QUAD_16:
                    return CellType.Square_16;

                case CurvedElementType_t.HEXA_8:
                    if (nonLinear) {
                        return CellType.Cube_8;
                    }
                    return CellType.Cube_Linear;

                case CurvedElementType_t.HEXA_20:
                    return CellType.Cube_20;

                case CurvedElementType_t.HEXA_27:
                    return CellType.Cube_27;

                case CurvedElementType_t.HEXA_64:
                    return CellType.Cube_64;

                case CurvedElementType_t.HEXA_125:
                    return CellType.Cube_125;

                default:
                    throw new NotSupportedException();
            }
        }

        /// <summary>
        /// Generates parameters for the transformation from reference to physical domain.
        /// All curved elements key nodes are in Gmesh numbering order.
        /// </summary>
        /// <param name="transformationParams"></param>
        /// <param name="l"></param>
        private void GenerateTransformationParams(MultidimensionalArray[] transformationParams, int l) {

            int k = bosss_element_index_to_msh_element_index[l];

            MultidimensionalArray trafo = transformationParams[l];

            Debug.Assert(trafo.GetLength(1) == number_of_coord_dimensions);

            double[,] coords = new double[number_of_coord_dimensions, trafo.GetLength(0)];

            for (int d = 0; d < number_of_coord_dimensions; d++) {
                for (int i = 0; i < trafo.GetLength(0); i++) {
                    int j = permutation_array[chirality_index[k], i];
                    if (d == 0) {
                        trafo[i, d] = x[elements[k][j]];
                    } else if (d == 1) {
                        trafo[i, d] = y[elements[k][j]];
                    } else {
                        trafo[i, d] = z[elements[k][j]];
                    }
                }
            }
        }

        private void GeneratePermutationArray(CurvedElementType_t c) {
            permutation_array = REORDER[(int)c - 1];
        }

        /// <summary>
        /// Generates the reference element.
        /// </summary>
        private RefElement GenerateReferenceElement() {
            int array_length = NumberOfNodesPerElement(curved_element_type);
            int i;

            ref_x = new double[array_length];
            ref_y = new double[array_length];
            ref_z = new double[array_length];
            ref_element = new int[array_length];

            RefElement refel = GetReferenceElement(bosss_element_type);
            MultidimensionalArray nodes = refel.GetInterpolationNodes(ConvertToCellType(curved_element_type, true)); // Always non linear
            int length = nodes.Length / number_of_coord_dimensions;
            for (i = 0; i < length; i++) {
                ref_x[i] = nodes[permutation_array[0, i], 0];
                ref_y[i] = nodes[permutation_array[0, i], 1];
                if (number_of_coord_dimensions > 2) {
                    ref_z[i] = nodes[permutation_array[0, i], 2];
                }
            }

            return refel;
        }

        /// <summary>
        /// Generates the reference element.
        /// </summary>
        private void DetermineChirality(int i, int numberOfNodes) {
            int[] indices = elements[i];

            int index1 = indices[0];
            int index2 = indices[1];

            double[] averagePoint = CalculateAveragePoint(indices, numberOfNodes);

            Vector u = new Vector(x[index1], y[index1]);
            Vector v = new Vector(x[index2], y[index2]);
            Vector w = new Vector(averagePoint[0], averagePoint[1]);

            u -= w;
            v -= w;

            Double crossProduct = CrossProduct(u, v);

            if (number_of_coord_dimensions == 3) {
                double crossProductAbs = Math.Abs(crossProduct);

                u = new Vector(x[index1], z[index1]);
                v = new Vector(x[index2], z[index2]);

                u -= w;
                v -= w;

                double crossProduct1 = CrossProduct(u, v);

                double crossProductAbs1 = Math.Abs(crossProduct);

                if (crossProductAbs1 > crossProductAbs) {
                    crossProduct = crossProduct1;
                    crossProductAbs = crossProductAbs1;
                }

                u = new Vector(y[index1], z[index1]);
                v = new Vector(y[index2], z[index2]);

                u -= w;
                v -= w;

                crossProduct1 = CrossProduct(u, v);

                crossProductAbs1 = Math.Abs(crossProduct);

                if (crossProductAbs1 > crossProductAbs) {
                    crossProduct = crossProduct1;
                }
            }

            if (crossProduct < 0) {
                chirality_index[i] = 1;
            } else {
                chirality_index[i] = 0;
            }
        }

        /// <summary>
        /// Generates the cross product.
        /// </summary>
        private static double CrossProduct(Vector u, Vector v) {
            return u.CrossProduct(v).z;
        }

        /// <summary>
        /// Helper function.
        /// </summary>
        private double[] CalculateAveragePoint(int[] indices, int length) {
            int i, j;

            double[] ret = new double[3];

            double normfactor = length;

            for (i = 0; i < normfactor; i++) {
                j = indices[i];

                ret[0] += x[j];
                ret[1] += y[j];
                if (number_of_coord_dimensions == 3) {
                    ret[2] += z[j];
                }
            }

            ret[0] /= normfactor;
            ret[1] /= normfactor;
            if (number_of_coord_dimensions == 3) {
                ret[2] /= normfactor;
            }

            return ret;
        }

        private Tuple<int, int> FindCellAndFaceIndex(int element, GridCommons grid, RefElement refElement) {
            int[] nodeIndices = elements[element];

            // Loop through all *real* elements and the one using the same
            // nodes for a face
            int cellIndex = -1;
            int faceIndex = -1;
            for (int i = 0; i < bosss_elements; i++) {
                int mshIndex = bosss_element_index_to_msh_element_index[i];
                int[] otherNodeIndices = elements[mshIndex];

                bool containsAllIndices = true;
                for (int j = 0; j < nodeIndices.Length; j++) {
                    containsAllIndices &= otherNodeIndices.Contains(nodeIndices[j]);
                }

                // Element does not contain all indices used by the BC element
                // => Not the element we're looking for;
                if (!containsAllIndices) {
                    continue;
                }

                cellIndex = i;

                // Determine face by looking for node indices (don't check
                // order, should have been checked elsewhere)
                for (int j = 0; j < refElement.FaceToVertexIndices.GetLength(0); j++) {
                    otherNodeIndices = new int[refElement.FaceToVertexIndices.GetLength(1)];
                    for (int k = 0; k < refElement.FaceToVertexIndices.GetLength(1); k++) {
                        otherNodeIndices[k] = grid.Cells[i].NodeIndices[
                            refElement.FaceToVertexIndices[j, k]];
                    }

                    containsAllIndices = true;
                    // Careful: $nodeIndices potentially contains more nodes
                    // (for higher order elements), so we have to compare this
                    // way round
                    for (int k = 0; k < otherNodeIndices.Length; k++) {
                        containsAllIndices &= nodeIndices.Contains(otherNodeIndices[k]);
                    }

                    if (containsAllIndices) {
                        faceIndex = j;
                        break;
                    }
                }

                if (faceIndex < 0) {
                    throw new Exception("This should not have happened");
                }

                break;
            }

            return Tuple.Create(cellIndex, faceIndex);
        }

        private static int NumberOfNodesPerElement(Enum type) {
            if (type.Equals(CurvedElementType_t.NODE) || type.Equals(ElementType_t.NODE)) {
                return 1;
            }
            if (type.Equals(CurvedElementType_t.TRI_152)) {
                return 15;
            }
            string s = type.ToString();
            string[] words = s.Split('_');
            return Convert.ToInt32(words[words.Length - 1]);
        }

        private static RefElement GetReferenceElement(Enum type) {
            if (type.Equals(ElementType_t.TRI_3) || type.Equals(ElementType_t.TRI_6)) {
                return Triangle.Instance;
            } else if (type.Equals(ElementType_t.TETRA_4) || type.Equals(ElementType_t.TETRA_10)) {
                return Tetra.Instance;
            } else if (type.Equals(ElementType_t.QUAD_4) || type.Equals(ElementType_t.QUAD_8) || type.Equals(ElementType_t.QUAD_9)) {
                return Square.Instance;
            } else if (type.Equals(ElementType_t.HEXA_8) || type.Equals(ElementType_t.HEXA_20) || type.Equals(ElementType_t.HEXA_27)) {
                return Cube.Instance;
            }
            throw new NotImplementedException("Wrong Elementtype");
        }

        private static bool Is3D(ElementType_t type) {
            return type.Equals(ElementType_t.TETRA_4) || type.Equals(ElementType_t.HEXA_8);
        }

        private static bool BoSSSCompliantElementType(ElementType_t type) {
            return type.Equals(ElementType_t.BAR_2) ||
                type.Equals(ElementType_t.TRI_3) ||
                type.Equals(ElementType_t.TETRA_4) ||
                type.Equals(ElementType_t.QUAD_4) ||
                type.Equals(ElementType_t.HEXA_8);
        }

        private bool ChoosenBoSSSCompliantElementType(ElementType_t type) {
            return type.Equals(bosss_element_type);
        }

        private static String CleanString(String s) {
            StringBuilder result = new StringBuilder();
            for (int i = 0; i < s.Length; i++) {
                char c = s[i];
                byte b = (byte)c;
                if (b >= 32) {
                    result.Append(c);
                }
            }
            return result.ToString();
        }

        #region IImporter Members

        /// <summary>
        /// Generates a grid commons object - this is the main function of this class. 
        /// All curved elements key nodes are in Gmesh numbering order.
        /// </summary>
        /// <returns>the grid commons object</returns>
        public GridCommons GenerateBoSSSGrid() {
            GridCommons grd = null;

            if (bosss_element_type.Equals(ElementType_t.QUAD_4)) {
                grd = new Grid2D(Square.Instance);
            } else if (bosss_element_type.Equals(ElementType_t.TRI_3)) {
                grd = new Grid2D(Triangle.Instance);
            } else if (bosss_element_type.Equals(ElementType_t.HEXA_8)) {
                grd = new Grid3D(Cube.Instance);
            } else if (bosss_element_type.Equals(ElementType_t.TETRA_4)) {
                grd = new Grid3D(Tetra.Instance);
            } else {
                throw new NotSupportedException("Unsupported BoSSS element type: " + bosss_element_type);
            }

            int numberOfNodesPerCurvedElement = NumberOfNodesPerElement(curved_element_type);
            int numberOfNodesPerElement = NumberOfNodesPerElement(bosss_element_type);

            for (int j = 0; j < bosss_elements; j++) {
                int i = bosss_element_index_to_msh_element_index[j];
                DetermineChirality(i, numberOfNodesPerCurvedElement);
            }
            
            GeneratePermutationArray(curved_element_type);

            for (int j = 0; j < bosss_elements; j++) {
                int i = bosss_element_index_to_msh_element_index[j];
                for (int k = 0; k < base_elements[bosss_element_index_to_msh_element_index[j]].Length; k++) {
                    // over coordinate indices per elements
                    int l = elements[i].Length - numberOfNodesPerElement + k;
                    base_elements[i][k] = elements[i][permutation_array[chirality_index[i], l]];
                }
            }

            RefElement refElement = GenerateReferenceElement();

            grd.Cells = new Cell[bosss_elements];

            MultidimensionalArray[] transformationParams = new MultidimensionalArray[bosss_elements];

            for (int i = 0; i < bosss_elements; i++) {
                grd.Cells[i] = new Cell();
                bool nonLinear = false;
                if (curved_element_type.Equals(CurvedElementType_t.QUAD_4) || curved_element_type.Equals(CurvedElementType_t.HEXA_8)) {
                    if (CheckIfElementIsNonLinear(numberOfNodesPerCurvedElement, i)) {
                        nonLinear = true;
                    }
                }
                if (curved_element_type.Equals(CurvedElementType_t.QUAD_4) || curved_element_type.Equals(CurvedElementType_t.HEXA_8)) {
                    if (nonLinear) {
                        transformationParams[i] = MultidimensionalArray.Create(numberOfNodesPerCurvedElement, number_of_coord_dimensions);
                    } else {
                        transformationParams[i] = MultidimensionalArray.Create(numberOfNodesPerCurvedElement - 1, number_of_coord_dimensions);
                    }
                } else {
                    transformationParams[i] = MultidimensionalArray.Create(numberOfNodesPerCurvedElement, number_of_coord_dimensions);
                }
                GenerateTransformationParams(transformationParams, i);
                grd.Cells[i].Type = ConvertToCellType(curved_element_type, nonLinear);
            }

            grd.Description = "Created from BoSSS";

            for (int j = 0; j < bosss_elements; j++) {
                grd.Cells[j].GlobalID = j;
                grd.Cells[j].NodeIndices = base_elements[bosss_element_index_to_msh_element_index[j]];
                grd.Cells[j].TransformationParams = transformationParams[j];
            }

            if (hasbc) {
                for (int i = 0; i < boconames.Length; i++) {
                    grd.EdgeTagNames.Add((byte)(i + 1), boconames[i]);
                }

                for (int element = 0; element < bcindex_of_element.Length; element++) {
                    int bcindex = bcindex_of_element[element];
                    if (bcindex >= physical_index_mapping.Length) {
                        // Not a BC element
                        continue;
                    }

                    byte edgeTag = (byte)physical_index_mapping[bcindex];
                    Tuple<int, int> cellAndFaceIndex = FindCellAndFaceIndex(element, grd, refElement);

                    CellFaceTag CFT = new CellFaceTag() {
                        EdgeTag = edgeTag,
                        FaceIndex = cellAndFaceIndex.Item2,
                        NeighCell_GlobalID = long.MinValue
                    };
                    CFT.AddToArray(ref grd.Cells[cellAndFaceIndex.Item1].CellFaceTags);
                }
            }
            grd.CompressNodeIndices();
            grd.CheckAndFixJacobianDeterminat();
            return grd;
        }

        #endregion
    }
}


using System;
using System.Collections.Generic;
using System.Text;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.GridImport {



    /// <summary>
    /// a cgns file (at least the sections that are important for BoSSS)
    /// </summary>
    public class Cgns {

        private double[][, ,] x_of_all_zones;
        private double[][, ,] y_of_all_zones;
        private double[][, ,] z_of_all_zones;
        private int[][][][] bcelements_of_all_zones;
        private int[][][][] normal_indices_of_all_zones;
        private int[][] starts_of_all_zones;
        private int[][] ends_of_all_zones;
        private int[][] nbndrys_of_all_zones;
        private PointSetType_t[][][] pointset_types_of_all_zones;
        private GridLocation_t[][][] gridlocation_types_of_all_zones;
        private BCType_t[][][] bc_types_of_all_zones;
        private ElementType_t[][] element_types_of_all_zones;
        private int[][][,] elements_of_all_zones;
        private char[][][][] boconames_of_all_zones;
        private ZoneType_t[] zone_types_of_all_zones;

        private double[] x;
        private double[] y;
        private double[] z;
        private int[][] bcelements; // Either elements or points!
        private int[][] normal_indices;
        private PointSetType_t[] pointset_types;
        private GridLocation_t[] gridlocation_types;
        private BCType_t[] bc_types;
        private ElementType_t[] element_types;
        private int[][] elements;
        private int[][][] faces;
        private int[,] bcfaces;
        private char[][] boconames;
        private ZoneType_t[] zone_types;
        private int[][] to_node_index;
        private int[][][] to_element_index;

        private Boolean need_struct_conversion;
        private int number_of_coord_dimensions = 2;
        private List<int>[] node_to_element_indices;
        private List<int>[] bc_to_element_indices;
        private int bosss_elements;
        private int number_of_face_edges;
        private int number_of_faces_per_element;
        private int number_of_faces_per_element_mask;
        private int[] cgns_element_index_to_bosss_element_index;
        private int[] bosss_element_index_to_cgns_element_index;
        private ElementType_t bosss_element_type;
        private long[,] cell_neighbours;

        ///// <summary>
        ///// Permutation tables for the conversion between Tecplot ordering and
        ///// BoSSS ordering of the nodes in an element.
        ///// </summary>
        //static public Dictionary<ElementType_t, int[]> elementTypeToPermutationTableMap = new Dictionary<ElementType_t, int[]> {
        //    {ElementType_t.BAR_2, new int[] {0, 1}},
        //    {ElementType_t.TRI_3, new int[] {0, 1, 2}},
        //    {ElementType_t.QUAD_4, new int[] {0, 1, 3, 2}},
        //    {ElementType_t.TETRA_4, new int[] {0, 1, 2, 3}},
        //    {ElementType_t.HEXA_8, new int[] {0, 1, 4, 2, 3, 6, 5, 7}}
        //};

        private void SwitchNode(ref int source, ref int target) {
            int buffer = source;

            source = target;
            target = buffer;
        }

        private void SwitchXZ() {
            double buffer;

            for (int i = 0; i < x.Length; i++) {
                buffer = x[i];
                x[i] = z[i];
                z[i] = buffer;
            }
        }

        private void SwitchYZ() {
            double buffer;

            for (int i = 0; i < x.Length; i++) {
                buffer = y[i];
                y[i] = z[i];
                z[i] = buffer;
            }
        }

        private void SwitchXY() {
            double buffer;

            for (int i = 0; i < x.Length; i++) {
                buffer = x[i];
                x[i] = y[i];
                y[i] = buffer;
            }
        }

        private void ThrowError(CgnsDriver cg) {
            throw new ApplicationException(cg.get_error());
        }

        private static void ThrowError(String s) {
            throw new ApplicationException(s);
        }

        private void ThrowError(int index_file, CgnsDriver cg) {
            String s = cg.get_error();

            if (cg.close(index_file) != (int)Error.CG_OK) {
                s = s + " " + cg.get_error();
            }

            throw new ApplicationException(s);
        }

        private void ThrowError(String s, int index_file, CgnsDriver cg) {
            if (cg.close(index_file) != (int)Error.CG_OK) {
                s = s + " " + cg.get_error();
            }

            throw new ApplicationException(s);
        }

        private int NumberOfFaceEdges(ElementType_t type) {
            if (Is3D(type)) {
                if (type.Equals(ElementType_t.HEXA_8)) {
                    return 4;
                }
                return 3;
            } else {
                if (number_of_coord_dimensions == 3) {
                    if (type.Equals(ElementType_t.TRI_3)) {
                        return 3;
                    }
                    if (type.Equals(ElementType_t.QUAD_4)) {
                        return 4;
                    }
                }
                return 2;
            }
        }

        private int NumberOfFacesPerElement(ElementType_t type) {
            if (number_of_coord_dimensions == 3) {
                if (type.Equals(ElementType_t.TRI_3) || type.Equals(ElementType_t.QUAD_4)) {
                    return 1;
                }
            }
            if (type.Equals(ElementType_t.HEXA_8)) {
                return 6;
            }
            char[] number = type.ToString().ToNullTermCharAry();
            return number[number.Length - 1] - 48;
        }

        private static int NumberOfNodesPerElement(ElementType_t type) {
            char[] number = type.ToString().ToNullTermCharAry();
            return number[number.Length - 1] - 48;
        }

        private static bool Is3D(ElementType_t type) {
            return type.Equals(ElementType_t.TETRA_4) || type.Equals(ElementType_t.HEXA_8);
        }

        private static bool BoSSSCompliantElementType(ElementType_t type) {
            return type.Equals(ElementType_t.TRI_3) || type.Equals(ElementType_t.TETRA_4) || type.Equals(ElementType_t.QUAD_4) || type.Equals(ElementType_t.HEXA_8);
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

        /// <summary>
        /// Convert a structured grid to an unstructured grid
        /// </summary>
        private void ConvertStructure2Unstructure() {
            // insert here
            need_struct_conversion = false;
        }

        int verticesPerCell;

        /// <summary>
        /// Generates zoneless arrays
        /// </summary>
        private void GenerateZonelessAndSectionlessArrays() {
            if (need_struct_conversion) {
                ThrowError("Struct conversion needed");
            }

            int i, j, k, l, m, n, a, b, c;

            zone_types = zone_types_of_all_zones;

            int array_length = 0;

            for (i = 0; i < zone_types.Length; i++) {
                array_length += x_of_all_zones[i].GetLength(0);
            }

            x = new double[array_length];
            y = new double[array_length];
            z = new double[array_length];

            node_to_element_indices = new List<int>[array_length];
            to_node_index = new int[zone_types.Length][];
            to_element_index = new int[zone_types.Length][][];
            for (i = 0, j = 0; i < zone_types.Length; i++) { // over zone indices
                a = x_of_all_zones[i].GetLength(0);
                to_node_index[i] = new int[a];
                for (k = 0; k < a; j++, k++) { // over coordinate indices per zone
                    x[j] = x_of_all_zones[i][k, 0, 0];
                    y[j] = y_of_all_zones[i][k, 0, 0];
                    z[j] = z_of_all_zones[i][k, 0, 0];
                    to_node_index[i][k] = j; // i|k index -> j index
                    node_to_element_indices[j] = new List<int>();
                }
            }

            array_length = 0;
            for (i = 0; i < elements_of_all_zones.Length; i++) {
                a = elements_of_all_zones[i].Length;
                for (k = 0; k < a; k++) {
                    array_length += elements_of_all_zones[i][k].GetLength(0);
                }
            }

            element_types = new ElementType_t[array_length];
            elements = new int[array_length][];
            for (i = 0, j = 0, m = 0; i < elements_of_all_zones.Length; i++) { // over zone indices
                a = elements_of_all_zones[i].Length;
                to_element_index[i] = new int[a][];
                for (k = 0; k < a; j++, k++) { // over section indices per zone
                    b = elements_of_all_zones[i][k].GetLength(0);
                    c = NumberOfNodesPerElement(element_types_of_all_zones[i][k]);
                    to_element_index[i][k] = new int[b];
                    for (l = 0; l < b; l++, m++) { // over element indices (per section) per zone
                        to_element_index[i][k][l] = m; // i|k|l index -> m index
                        elements[m] = new int[c];
                        element_types[m] = element_types_of_all_zones[i][k];
                        if (Is3D(element_types[m])) {
                            number_of_coord_dimensions = 3;
                        }
                        if (BoSSSCompliantElementType(element_types[m])) {
                            if (verticesPerCell < c) {
                                verticesPerCell = c;
                                bosss_element_type = element_types[m];
                            }
                            if (number_of_face_edges < NumberOfFaceEdges(element_types[m])) {
                                number_of_face_edges = NumberOfFaceEdges(element_types[m]);
                                bosss_element_type = element_types[m];
                            }
                            if (number_of_faces_per_element < NumberOfFacesPerElement(element_types[m])) {
                                number_of_faces_per_element = NumberOfFacesPerElement(element_types[m]);
                                number_of_faces_per_element_mask = (1 << number_of_faces_per_element) - 1;
                                bosss_element_type = element_types[m];
                            }
                        }
                        for (n = 0; n < c; n++) { // over coordinate indices per elements (per section) per zone
                            node_to_element_indices[elements[m][n] = to_node_index[i][elements_of_all_zones[i][k][l, n] - 1]].Add(m);
                        }
                    }
                }
            }
            for (m = 0; m < elements.Length; m++) {
                if (BoSSSCompliantElementType(element_types[m]) && (number_of_coord_dimensions == 2 || Is3D(element_types[m]))) {
                    c = NumberOfNodesPerElement(element_types[m]);
                    if (verticesPerCell > c) {
                        verticesPerCell = c;
                        bosss_element_type = element_types[m];
                    }
                    if (number_of_face_edges > NumberOfFaceEdges(element_types[m])) {
                        number_of_face_edges = NumberOfFaceEdges(element_types[m]);
                        bosss_element_type = element_types[m];
                    }
                    if (number_of_faces_per_element > NumberOfFacesPerElement(element_types[m])) {
                        number_of_faces_per_element = NumberOfFacesPerElement(element_types[m]);
                        number_of_faces_per_element_mask = (1 << number_of_faces_per_element) - 1;
                        bosss_element_type = element_types[m];
                    }
                }
            }
            cgns_element_index_to_bosss_element_index = new int[elements.Length];
            bosss_element_index_to_cgns_element_index = new int[elements.Length];
            for (m = 0; m < elements.Length; m++) {
                if (ChoosenBoSSSCompliantElementType(element_types[m])) {
                    bosss_element_index_to_cgns_element_index[cgns_element_index_to_bosss_element_index[m] = bosss_elements++] = m;
                }
            }
            array_length = 0;
            for (i = 0; i < bcelements_of_all_zones.Length; i++) { // over zone indices
                a = bcelements_of_all_zones[i].Length;
                for (k = 0; k < a; k++) { // over section indices per zone
                    array_length += bcelements_of_all_zones[i][k].Length;
                }
            }

            gridlocation_types = new GridLocation_t[array_length];
            pointset_types = new PointSetType_t[array_length];
            bc_types = new BCType_t[array_length];
            boconames = new char[array_length][];
            bcelements = new int[array_length][];
            normal_indices = new int[array_length][];
            for (i = 0, j = 0, m = 0; i < bcelements_of_all_zones.Length; i++) { // over zone indices
                a = bcelements_of_all_zones[i].Length;
                for (k = 0; k < a; j++, k++) { // over section indices per zone
                    b = bcelements_of_all_zones[i][k].Length;
                    for (l = 0; l < b; l++, m++) { // over boundary condition indices (per section) per zone
                        gridlocation_types[m] = gridlocation_types_of_all_zones[i][k][l];

                        if (!gridlocation_types[m].Equals(GridLocation_t.FaceCenter) && !gridlocation_types[m].Equals(GridLocation_t.Vertex)) {
                            ThrowError("Unsupported gridlocation");
                        }

                        pointset_types[m] = pointset_types_of_all_zones[i][k][l];
                        bc_types[m] = bc_types_of_all_zones[i][k][l];
                        boconames[m] = boconames_of_all_zones[i][k][l];
                        normal_indices[m] = normal_indices_of_all_zones[i][k][l];
                        c = bcelements_of_all_zones[i][k][l].Length;
                        bcelements[m] = new int[c];
                        for (n = 0; n < c; n++) { // over coordinate indices per boundary condition (per section) per zone
                            if (gridlocation_types[m].Equals(GridLocation_t.Vertex)) {
                                bcelements[m][n] = to_node_index[i][bcelements_of_all_zones[i][k][l][n] - 1];
                            } else {
                                bcelements[m][n] = bcelements_of_all_zones[i][k][l][n] - 1;
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Reorder the coordinates inside an element
        /// </summary>
        private void VertexIndexNumberingTransformation() {
            for (int i = 0; i < elements.Length; i++) {
                if (element_types[i].Equals(ElementType_t.QUAD_4)) {
                    SwitchNode(ref elements[i][2], ref elements[i][3]);
                }
                if (!element_types[i].Equals(ElementType_t.HEXA_8)) {
                    continue;
                }
                SwitchNode(ref elements[i][2], ref elements[i][3]);
                SwitchNode(ref elements[i][3], ref elements[i][4]);
                SwitchNode(ref elements[i][5], ref elements[i][6]);
            }
        }

        /// <summary>
        /// Reorder the coordinates inside an element
        /// </summary>
        private void CoordinateSystemTransformation() {
            // insert here
        }

        /// <summary>
        /// Builds faces from elements
        /// </summary>
        private void BuildFaces() {
            int i, j;

            faces = new int[elements.Length][][];

            bc_to_element_indices = new List<int>[bc_types.Length];

            for (i = 0; i < bc_types.Length; i++) { // over elements
                bc_to_element_indices[i] = new List<int>();
            }
            for (i = 0; i < faces.Length; i++) { // over elements
                if (element_types[i].Equals(ElementType_t.BAR_2)) {
                    faces[i] = new int[1][];
                    faces[i][0] = new int[2];
                    faces[i][0][0] = elements[i][0];
                    faces[i][0][1] = elements[i][1];
                } else {
                    faces[i] = new int[NumberOfFacesPerElement(element_types[i])][];

                    for (j = 0; j < faces[i].Length; j++) {
                        faces[i][j] = new int[NumberOfFaceEdges(element_types[i])];
                    }
                    if (element_types[i].Equals(ElementType_t.TRI_3)) {
                        if (number_of_coord_dimensions == 2) {
                            faces[i][0][0] = elements[i][0];
                            faces[i][0][1] = elements[i][1];
                            faces[i][1][0] = elements[i][1];
                            faces[i][1][1] = elements[i][2];
                            faces[i][2][0] = elements[i][2];
                            faces[i][2][1] = elements[i][0];
                        } else {
                            faces[i][0][0] = elements[i][0];
                            faces[i][0][1] = elements[i][1];
                            faces[i][0][2] = elements[i][2];
                        }
                    } else if (element_types[i].Equals(ElementType_t.TETRA_4)) {
                        faces[i][0][0] = elements[i][1];
                        faces[i][0][1] = elements[i][2];
                        faces[i][0][2] = elements[i][3];
                        faces[i][1][0] = elements[i][0];
                        faces[i][1][1] = elements[i][2];
                        faces[i][1][2] = elements[i][3];
                        faces[i][2][0] = elements[i][0];
                        faces[i][2][1] = elements[i][1];
                        faces[i][2][2] = elements[i][2];
                        faces[i][3][0] = elements[i][0];
                        faces[i][3][1] = elements[i][1];
                        faces[i][3][2] = elements[i][3];
                    } else if (element_types[i].Equals(ElementType_t.QUAD_4)) {
                        if (number_of_coord_dimensions == 2) {
                            faces[i][0][0] = elements[i][0];
                            faces[i][0][1] = elements[i][2];
                            faces[i][1][0] = elements[i][1];
                            faces[i][1][1] = elements[i][3];
                            faces[i][2][0] = elements[i][2];
                            faces[i][2][1] = elements[i][3];
                            faces[i][3][0] = elements[i][0];
                            faces[i][3][1] = elements[i][1];
                        } else {
                            faces[i][0][0] = elements[i][0];
                            faces[i][0][1] = elements[i][1];
                            faces[i][0][2] = elements[i][2];
                            faces[i][0][3] = elements[i][3];
                        }
                    } else if (element_types[i].Equals(ElementType_t.HEXA_8)) {
                        faces[i][0][0] = elements[i][0];
                        faces[i][0][1] = elements[i][3];
                        faces[i][0][2] = elements[i][7];
                        faces[i][0][3] = elements[i][2];
                        faces[i][1][0] = elements[i][1];
                        faces[i][1][1] = elements[i][4];
                        faces[i][1][2] = elements[i][5];
                        faces[i][1][3] = elements[i][6];
                        faces[i][2][0] = elements[i][2];
                        faces[i][2][1] = elements[i][4];
                        faces[i][2][2] = elements[i][5];
                        faces[i][2][3] = elements[i][7];
                        faces[i][3][0] = elements[i][0];
                        faces[i][3][1] = elements[i][1];
                        faces[i][3][2] = elements[i][6];
                        faces[i][3][3] = elements[i][3];
                        faces[i][4][0] = elements[i][3];
                        faces[i][4][1] = elements[i][7];
                        faces[i][4][2] = elements[i][5];
                        faces[i][4][3] = elements[i][6];
                        faces[i][5][0] = elements[i][0];
                        faces[i][5][1] = elements[i][2];
                        faces[i][5][2] = elements[i][4];
                        faces[i][5][3] = elements[i][1];
                    }
                }
            }
        }

        /// <summary>
        /// Builds and shares the boundary conditions of one point with all other points positioned at the same location
        /// </summary>
        private void BuildAndShareBoundaryConditions() {
            int i, j, k, l, m, n;
            int range_from;
            int range_to;
            int index;
            int index2;
            int edge_counter;
            int verify;

            bcfaces = new int[elements.Length, number_of_faces_per_element];
            for (i = 0; i < bcfaces.GetLength(0); i++) {
                for (j = 0; j < bcfaces.GetLength(1); j++) {
                    bcfaces[i, j] = -1;
                }
            }
            for (i = 0; i < bcelements.Length; i++) { // over boundary condition indices
                range_from = 0;
                range_to = bcelements[i].Length;

                if (pointset_types[i].Equals(PointSetType_t.ElementRange) || pointset_types[i].Equals(PointSetType_t.PointRange)) {
                    range_from = bcelements[i][0];
                    range_to = bcelements[i][1] + 1;
                }
                for (j = range_from; j < range_to; j++) { // over coordinate indices per boundary condition (if branch) or element indices per boundary condition (else branch)
                    if (pointset_types[i].Equals(PointSetType_t.ElementRange) || pointset_types[i].Equals(PointSetType_t.PointRange)) {
                        index = j;
                    } else {
                        index = bcelements[i][j];
                    }
                    if (gridlocation_types[i].Equals(GridLocation_t.Vertex)) {
                        foreach (int z in node_to_element_indices[index]) { // over elements
                            for (k = 0; k < faces[z].Length; k++) { // over face indices per element
                                if (k >= bcfaces.GetLength(1)) {
                                    continue;
                                }
                                for (m = edge_counter = 0; m < faces[z][k].Length; m++) { // over coordinate indices per faces
                                    for (n = range_from; n < range_to; n++) { // over coordinate indices per boundary condition
                                        if (pointset_types[i].Equals(PointSetType_t.ElementRange) || pointset_types[i].Equals(PointSetType_t.PointRange)) {
                                            index2 = n;
                                        } else {
                                            index2 = bcelements[i][n];
                                        }
                                        if (index2 == faces[z][k][m]) {
                                            edge_counter++;
                                        }
                                    }
                                }
                                if (edge_counter == number_of_face_edges) {
                                    bcfaces[z, k] = i;
                                }
                            }
                        }
                    } else {
                        for (k = 0; k < faces[j].Length; k++) { // over face indices per element
                            if (k >= bcfaces.GetLength(1)) {
                                continue;
                            }
                            bcfaces[j, k] = i;
                        }
                    }
                }
                for (j = 0; j < faces.Length; j++) { // over element indices
                    for (k = 0; k < faces[j].Length; k++) { // over face indices per element
                        foreach (int z in node_to_element_indices[faces[j][k][0]]) { // over elements
                            if (z == j || !ChoosenBoSSSCompliantElementType(element_types[z])) {
                                continue;
                            }
                            for (l = 0; l < faces[z].Length; l++) { // over face indices per element
                                if (l >= bcfaces.GetLength(1)) {
                                    continue;
                                }
                                for (m = edge_counter = 0; m < faces[z][l].Length; m++) { // over coordinate indices per faces
                                    for (n = 0; n < faces[j][k].Length; n++) { // over coordinate indices per faces
                                        if (faces[j][k][n] == faces[z][l][m]) {
                                            edge_counter++;
                                        }
                                    }
                                }
                                if (edge_counter == number_of_face_edges) {
                                    if (bcfaces[z, l] == -1) {
                                        bcfaces[z, l] = bcfaces[j, k];
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (bc_to_element_indices.Length == 0) {
                return;
            }
            for (i = 0; i < faces.Length; i++) { // over elements
                if (!ChoosenBoSSSCompliantElementType(element_types[i])) {
                    continue;
                }
                for (j = 0; j < faces[i].Length; j++) { // over faces
                    if ((verify = bcfaces[i, j]) >= 0) {
                        bc_to_element_indices[verify].Add(i << number_of_faces_per_element | j);
                    }
                }
            }
        }

        /// <summary>
        /// Build neighbor relations
        /// </summary>
        private void BuildNeighbourship() {
            int i, j, k, l, m, n;
            int edge_counter;

            cell_neighbours = new long[bosss_elements, number_of_faces_per_element];

            for (i = 0; i < bosss_elements; i++) { // over elements
                for (j = 0; j < number_of_faces_per_element; j++) { // over face indices per element
                    cell_neighbours[i, j] = -1;
                }
            }
            for (n = 0; n < bosss_elements; n++) { // over elements
                i = bosss_element_index_to_cgns_element_index[n];
                if (!ChoosenBoSSSCompliantElementType(element_types[i])) {
                    continue;
                }
                for (j = 0; j < faces[i].Length; j++) { // over face indices per element
                    if (cell_neighbours[cgns_element_index_to_bosss_element_index[i], j] > -1) {
                        continue;
                    }
                    foreach (int z in node_to_element_indices[faces[i][j][0]]) { // over elements
                        if (z == i || !ChoosenBoSSSCompliantElementType(element_types[z])) {
                            continue;
                        }
                        for (l = 0; l < faces[z].Length; l++) { // over face indices per element
                            for (m = edge_counter = 0; m < faces[z][l].Length; m++) { // over coordinate indices per faces
                                for (k = 0; k < faces[i][j].Length; k++) { // over coordinate indices per faces
                                    if (faces[i][j][k] == faces[z][l][m]) {
                                        edge_counter++;
                                    }
                                }
                            }
                            if (edge_counter == number_of_face_edges) {
                                cell_neighbours[cgns_element_index_to_bosss_element_index[i], j] = cgns_element_index_to_bosss_element_index[z];
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Replaces old edge tag number to new one in case of duplicate boundary conditions
        /// </summary>
        /// <param name="EdgeTags"></param>
        /// <param name="oldValue"></param>
        /// <param name="newValue"></param>
        private void replaceEdgeTagNumerInCaseOfDublicateBoundaryConditions(byte[,] EdgeTags, byte oldValue, byte newValue) {
            for (int i = 0; i < EdgeTags.GetLength(0); i++) {
                for (int j = 0; j < EdgeTags.GetLength(1); j++) {
                    if (EdgeTags[i, j] == oldValue) {
                        EdgeTags[i, j] = newValue;
                    }
                }
            }
        }

        /// <summary>
        /// loads or saves an adf or hdf5 file from or to a specified path
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="hdf5"></param>
        public Cgns(string filePath, bool hdf5) {
            Build(filePath, hdf5);
        }

        /// <summary>
        /// loads an adf or hdf5 file from or to a specified path
        /// </summary>
        void Build(string filePath, bool hdf5) {
            using (CgnsDriver cg = new CgnsDriver(hdf5)) {
                int index_file = 0;
                int index_base = 0;
                int index_zone = 0;
                int index_section = 0;
                int index_bc = 0;
                int nbases;
                int nzones;
                int ncoords;
                int nsections;
                int nbocos;
                char[] basename = new char[100];
                char[] zonename = new char[100];
                char[] coordnameX = new char[100];
                char[] coordnameY = new char[100];
                char[] coordnameZ = new char[100];
                char[] flow = new char[100];
                char[] secname = new char[100];
                char[] boconame;
                char[] solname = new char[100];
                int[,] isize = new int[3, 3];
                int[,] ielem;
                int[] irmin = new int[3];
                int[] irmax = new int[3];
                int celldim, physdim; // zone definition
                int start, end, nbndry, parent_flag; // element definition
                int normal_list_flag, nbcelement, ndataset; // bc definition
                DataType_t data_type;
                ElementType_t element_type;
                ZoneType_t zone_type;
                BCType_t bc_type;
                GridLocation_t gridlocation_type;
                PointSetType_t pointset_type;

                index_file = index_base = 1;

                if (cg.open(filePath.ToNullTermCharAry(), (int)Mode.MODE_READ, out index_file) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                if (cg.nbases(index_file, out nbases) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                if (nbases != 1) {
                    ThrowError("Only one base allowed", index_file, cg);
                }
                if (cg.base_read(index_file, index_base, basename, out celldim, out physdim) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                if (cg.nzones(index_file, index_base, out nzones) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }

                x_of_all_zones = new double[nzones][, ,];
                y_of_all_zones = new double[nzones][, ,];
                z_of_all_zones = new double[nzones][, ,];
                bcelements_of_all_zones = new int[nzones][][][];
                normal_indices_of_all_zones = new int[nzones][][][];
                starts_of_all_zones = new int[nzones][];
                ends_of_all_zones = new int[nzones][];
                nbndrys_of_all_zones = new int[nzones][];
                pointset_types_of_all_zones = new PointSetType_t[nzones][][];
                gridlocation_types_of_all_zones = new GridLocation_t[nzones][][];
                bc_types_of_all_zones = new BCType_t[nzones][][];
                element_types_of_all_zones = new ElementType_t[nzones][];
                elements_of_all_zones = new int[nzones][][,];
                boconames_of_all_zones = new char[nzones][][][];
                zone_types_of_all_zones = new ZoneType_t[nzones];

                for (index_zone = 1; index_zone <= nzones; index_zone++) {
                    int zone_offset = 0;

                    if (cg.zone_type(index_file, index_base, index_zone, out zone_type) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    zone_types_of_all_zones[index_zone - 1 - zone_offset] = zone_type;
                    if (cg.zone_read(index_file, index_base, index_zone, zonename, isize) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (cg.ncoords(index_file, index_base, index_zone, out ncoords) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (cg.coord_info(index_file, index_base, index_zone, 1, out data_type, coordnameX) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (cg.coord_info(index_file, index_base, index_zone, 2, out data_type, coordnameY) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (celldim > 2) {
                        if (cg.coord_info(index_file, index_base, index_zone, 3, out data_type, coordnameZ) != (int)Error.CG_OK) {
                            ThrowError(index_file, cg);
                        }
                    }

                    irmin[0] = 1;
                    irmin[1] = 1;
                    irmin[2] = 1;
                    irmax[0] = isize[0, 0];
                    irmax[1] = isize[0, 1];
                    irmax[2] = isize[0, 2];

                    double[, ,] x = null;
                    double[, ,] y = null;
                    double[, ,] z = null;

                    if ((int)zone_type == (int)ZoneType_t.Unstructured) {
                        x = new double[isize[0, 0], 1, 1];
                        y = new double[isize[0, 0], 1, 1];
                        z = new double[isize[0, 0], 1, 1];
                    } else {
                        x = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                        y = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                        z = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                    }

                    char[][][] boconames = null;
                    char[][] boconametemp = null;
                    int[][][] bcelements = null;
                    int[][] bcelementtemp = null;
                    int[] bcelement = null;
                    int[][][] normal_indices = null;
                    int[][] normal_indextemp = null;
                    int[] normal_index = null;
                    int[] starts = null;
                    int[] ends = null;
                    int[] nbndrys = null;
                    PointSetType_t[][] pointset_types = null;
                    PointSetType_t[] pointset_typtemp = null;
                    GridLocation_t[][] gridlocation_types = null;
                    GridLocation_t[] gridlocation_typtemp = null;
                    BCType_t[][] bc_types = null;
                    BCType_t[] bc_typtemp = null;
                    ElementType_t[] element_types = null;
                    int[][,] elements = null;

                    if (cg.coord_read3(index_file, index_base, index_zone, coordnameX, data_type, irmin, irmax, x) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (cg.coord_read3(index_file, index_base, index_zone, coordnameY, data_type, irmin, irmax, y) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (celldim > 2) {
                        if (cg.coord_read3(index_file, index_base, index_zone, coordnameZ, data_type, irmin, irmax, z) != (int)Error.CG_OK) {
                            ThrowError(index_file, cg);
                        }
                    }
                    if (cg.nsections(index_file, index_base, index_zone, out nsections) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }

                    if ((int)zone_type == (int)ZoneType_t.Unstructured) {


                        boconames = new char[nsections][][];
                        bcelements = new int[nsections][][];
                        normal_indices = new int[nsections][][];
                        pointset_types = new PointSetType_t[nsections][];
                        gridlocation_types = new GridLocation_t[nsections][];
                        bc_types = new BCType_t[nsections][];
                        starts = new int[nsections];
                        ends = new int[nsections];
                        nbndrys = new int[nsections];
                        element_types = new ElementType_t[nsections];
                        elements = new int[nsections][,];


                        //



                        for (index_section = 1; index_section <= nsections; index_section++) {
                            if (cg.section_read(index_file, index_base, index_zone, index_section, secname, out element_type, out start, out end, out nbndry, out parent_flag) != (int)Error.CG_OK) {
                                ThrowError(index_file, cg);
                            }

                            starts[index_section - 1] = start;
                            ends[index_section - 1] = end;
                            nbndrys[index_section - 1] = nbndry;
                            element_types[index_section - 1] = element_type;

                            char[] c = element_type.ToString().ToNullTermCharAry();

                            ielem = new int[end - start + 1, c[c.Length - 1] - 48];

                            if (cg.elements_read2(index_file, index_base, index_zone, index_section, ielem, null) != (int)Error.CG_OK) {
                                ThrowError(index_file, cg);
                            }

                            elements[index_section - 1] = ielem;

                            if (cg.nbocos(index_file, index_base, index_zone, out nbocos) != (int)Error.CG_OK) {
                                ThrowError(index_file, cg);
                            }
                            boconametemp = new char[nbocos][];
                            bcelementtemp = new int[nbocos][];
                            normal_indextemp = new int[nbocos][];
                            pointset_typtemp = new PointSetType_t[nbocos];
                            gridlocation_typtemp = new GridLocation_t[nbocos];
                            bc_typtemp = new BCType_t[nbocos];

                            int depht = 3;
                            char[] label0 = "Zone_t".ToNullTermCharAry(), label1 = "ZoneBC_t".ToNullTermCharAry(), label2 = "BC_t".ToNullTermCharAry();
                            string[] labels = new string[depht];
                            labels[0] = new String(label0);
                            labels[1] = new String(label1);
                            labels[2] = new String(label2);

                            int[] num = new int[depht];

                            num[0] = index_zone;
                            num[1] = 1;
                            for (index_bc = 1; index_bc <= nbocos; index_bc++) {
                                boconame = new char[100];
                                normal_index = new int[3];
                                if (cg.boco_info(index_file, index_base, index_zone, index_bc, boconame, out bc_type, out pointset_type, out nbcelement, normal_index, out normal_list_flag, out data_type, out ndataset) != (int)Error.CG_OK) {
                                    ThrowError(index_file, cg);
                                }
                                bcelement = new int[nbcelement];
                                if (cg.boco_read1(index_file, index_base, index_zone, index_bc, bcelement, null) != (int)Error.CG_OK) {
                                    ThrowError(index_file, cg);
                                }
                                num[2] = index_bc;
                                if (pointset_type.Equals(PointSetType_t.ElementList) || pointset_type.Equals(PointSetType_t.ElementRange)) {
                                    gridlocation_type = GridLocation_t.FaceCenter;
                                } else {
                                    if (cg.golist(index_file, index_base, depht, labels, num) != (int)Error.CG_OK) {
                                        ThrowError(index_file, cg);
                                    }
                                    if (cg.gridlocation_read(out gridlocation_type) != (int)Error.CG_OK) {
                                        ThrowError(index_file, cg);
                                    }
                                }
                                boconametemp[index_bc - 1] = boconame;
                                bcelementtemp[index_bc - 1] = bcelement;
                                normal_indextemp[index_bc - 1] = normal_index;
                                pointset_typtemp[index_bc - 1] = pointset_type;
                                bc_typtemp[index_bc - 1] = bc_type;
                                gridlocation_typtemp[index_bc - 1] = gridlocation_type;
                            }
                            boconames[index_section - 1] = boconametemp;
                            bcelements[index_section - 1] = bcelementtemp;
                            normal_indices[index_section - 1] = normal_indextemp;
                            pointset_types[index_section - 1] = pointset_typtemp;
                            bc_types[index_section - 1] = bc_typtemp;
                            gridlocation_types[index_section - 1] = gridlocation_typtemp;
                        }
                    } else {
                        need_struct_conversion = true;
                        // insert here ...
                    }
                    x_of_all_zones[index_zone - 1 - zone_offset] = x;
                    y_of_all_zones[index_zone - 1 - zone_offset] = y;
                    z_of_all_zones[index_zone - 1 - zone_offset] = z;
                    bcelements_of_all_zones[index_zone - 1 - zone_offset] = bcelements;
                    normal_indices_of_all_zones[index_zone - 1 - zone_offset] = normal_indices;
                    starts_of_all_zones[index_zone - 1 - zone_offset] = starts;
                    ends_of_all_zones[index_zone - 1 - zone_offset] = ends;
                    nbndrys_of_all_zones[index_zone - 1 - zone_offset] = nbndrys;
                    pointset_types_of_all_zones[index_zone - 1 - zone_offset] = pointset_types;
                    bc_types_of_all_zones[index_zone - 1 - zone_offset] = bc_types;
                    element_types_of_all_zones[index_zone - 1 - zone_offset] = element_types;
                    elements_of_all_zones[index_zone - 1 - zone_offset] = elements;
                    boconames_of_all_zones[index_zone - 1 - zone_offset] = boconames;
                    gridlocation_types_of_all_zones[index_zone - 1 - zone_offset] = gridlocation_types;

                }
                if (cg.close(index_file) != (int)Error.CG_OK) {
                    ThrowError(cg);
                }
            }
        }
    }
}

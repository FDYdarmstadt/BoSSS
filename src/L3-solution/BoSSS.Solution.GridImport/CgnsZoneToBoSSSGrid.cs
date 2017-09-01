using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.GridImport {
#pragma warning disable  1591
    public class CgnsZoneToBoSSSGrid {

        CgnsFile.CgnsZone m_zone;

        public CgnsZoneToBoSSSGrid(CgnsFile.CgnsZone _zone) {
            m_zone = _zone;
        }

        private double[] x;
        private double[] y;
        private double[] z;
        private ZoneType_t[] zone_types;
        private int[][] to_node_index;
        private int[][][] to_element_index;
        private List<int>[] node_to_element_indices;
        private ElementType_t[] element_types;
        private int[][] elements;
        private int number_of_coord_dimensions = 2;
        private int verticesPerCell;
        private ElementType_t bosss_element_type;
        private int number_of_faces_per_element;
        private int number_of_face_edges;
        private int number_of_faces_per_element_mask;
        private int[] cgns_element_index_to_bosss_element_index;
        private int[] bosss_element_index_to_cgns_element_index;
        private PointSetType_t[] pointset_types;
        private GridLocation_t[] gridlocation_types;
        private char[][] boconames;
        private BCType_t[] bc_types;
        private int[][] bcelements; // Either elements or points!
        private int[][] normal_indices;
        //private int[][][] faces;
        private List<int>[] bc_to_element_indices;
        private int bosss_elements;
        //private long[,] cell_neighbours;
        private int[,] bcfaces;

        /// <summary>
        /// Generates zoneless arrays
        /// </summary>
        private void GenerateZonelessAndSectionlessArrays() {

            int i, j, k, l, m, n, a, b, c;

            zone_types = new ZoneType_t[] { m_zone.zone_type }; // zone_types_of_all_zones;

            int array_length = 0;

            for (i = 0; i < zone_types.Length; i++) {
                array_length += m_zone.x.GetLength(0);    // x_of_all_zones[i].GetLength(0);
            }

            x = new double[array_length];
            y = new double[array_length];
            z = new double[array_length];

            node_to_element_indices = new List<int>[array_length];
            to_node_index = new int[zone_types.Length][];
            to_element_index = new int[zone_types.Length][][];
            for (i = 0, j = 0; i < zone_types.Length; i++) { // over zone indices
                a = m_zone.x.GetLength(0);    //x_of_all_zones[i].GetLength(0);
                to_node_index[i] = new int[a];
                for (k = 0; k < a; j++, k++) { // over coordinate indices per zone
                    x[j] = m_zone.x[k, 0, 0];  //x_of_all_zones[i][k, 0, 0];
                    y[j] = m_zone.y[k, 0, 0]; //y_of_all_zones[i][k, 0, 0];
                    z[j] = m_zone.z[k, 0, 0]; //z_of_all_zones[i][k, 0, 0];
                    to_node_index[i][k] = j; // i|k index -> j index
                    node_to_element_indices[j] = new List<int>();
                }
            }

            array_length = 0;
            //for (i = 0; i < elements_of_all_zones.Length; i++) {
            {
                a = m_zone.elements.Length; // elements_of_all_zones[i].Length;
                for (k = 0; k < a; k++) {
                    array_length += m_zone.elements[k].GetLength(0);  // elements_of_all_zones[i][k].GetLength(0);
                }
            }

            element_types = new ElementType_t[array_length];
            elements = new int[array_length][];
            //for (i = 0, j = 0, m = 0; i < elements_of_all_zones.Length; i++) // over zone indices
            j = 0;
            m = 0;
            i = 0;
            {
                a = m_zone.elements.Length; // elements_of_all_zones[i].Length;
                to_element_index[i] = new int[a][];
                for (k = 0; k < a; j++, k++) { // over section indices per zone
                    b = m_zone.elements[k].GetLength(0); //elements_of_all_zones[i][k].GetLength(0);
                    c = NumberOfNodesPerElement(m_zone.element_types[k]); // element_types_of_all_zones[i][k]);
                    to_element_index[i][k] = new int[b];
                    for (l = 0; l < b; l++, m++) { // over element indices (per section) per zone
                        to_element_index[i][k][l] = m; // i|k|l index -> m index
                        elements[m] = new int[c];
                        element_types[m] = m_zone.element_types[k]; // element_types_of_all_zones[i][k];
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
                            node_to_element_indices[elements[m][n] = to_node_index[i][m_zone.elements/*elements_of_all_zones[i]*/[k][l, n] - 1]].Add(m);
                        }
                    }
                }
            }
            // Michael: Das gehört hier raus! War wohl vergessen worden zu löschen (siehe direkt oben)
            /*
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
            }*/
            cgns_element_index_to_bosss_element_index = new int[elements.Length];
            bosss_element_index_to_cgns_element_index = new int[elements.Length];
            for (m = 0; m < elements.Length; m++) {
                if (ChoosenBoSSSCompliantElementType(element_types[m])) {
                    bosss_element_index_to_cgns_element_index[cgns_element_index_to_bosss_element_index[m] = bosss_elements++] = m;
                }
            }
            array_length = 0;
            //for (i = 0; i < bcelements_of_all_zones.Length; i++) // over zone indices
            {

                a = /*bcelements_of_all_zones[i]*/m_zone.bcelements.Length;
                for (k = 0; k < a; k++) { // over section indices per zone
                    array_length += /*bcelements_of_all_zones[i]*/m_zone.bcelements[k].Length;
                }
            }

            gridlocation_types = new GridLocation_t[array_length];
            pointset_types = new PointSetType_t[array_length];
            bc_types = new BCType_t[array_length];
            boconames = new char[array_length][];
            bcelements = new int[array_length][];
            normal_indices = new int[array_length][];
            //for (i = 0, j = 0, m = 0; i < bcelements_of_all_zones.Length; i++) // over zone indices
            i = 0;
            j = 0;
            m = 0;
            {
                a = /*bcelements_of_all_zones[i]*/m_zone.bcelements.Length;
                for (k = 0; k < a; j++, k++) { // over section indices per zone
                    b = /*bcelements_of_all_zones[i]*/m_zone.bcelements[k].Length;
                    for (l = 0; l < b; l++, m++) { // over boundary condition indices (per section) per zone
                        gridlocation_types[m] = /*gridlocation_types_of_all_zones[i]*/m_zone.gridlocation_types[k][l];

                        // Michael: Wieso ist diese Einschränkung hier? Gab es mal Probleme ohne diese Exception? Mache ich mal raus.
                        if (!gridlocation_types[m].Equals(GridLocation_t.FaceCenter) && !gridlocation_types[m].Equals(GridLocation_t.Vertex)) {
                            // ThrowError("Unsupported gridlocation");
                        }

                        pointset_types[m] = /*pointset_types_of_all_zones[i]*/m_zone.pointset_types[k][l];
                        bc_types[m] = /*bc_types_of_all_zones[i]*/m_zone.bc_types[k][l];
                        boconames[m] = /*boconames_of_all_zones[i]*/m_zone.boconames[k][l];
                        normal_indices[m] = /*normal_indices_of_all_zones[i]*/m_zone.normal_indices[k][l];
                        c = /*bcelements_of_all_zones[i]*/m_zone.bcelements[k][l].Length;
                        bcelements[m] = new int[c];
                        for (n = 0; n < c; n++) { // over coordinate indices per boundary condition (per section) per zone
                            if (gridlocation_types[m].Equals(GridLocation_t.Vertex)) {
                                bcelements[m][n] = to_node_index[i][/*bcelements_of_all_zones[i]*/m_zone.bcelements[k][l][n] - 1];
                            } else {
                                bcelements[m][n] = /*bcelements_of_all_zones[i]*/ m_zone.bcelements[k][l][n] - 1;
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        private void CoordinateSystemTransformation() {
            // insert here
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

        private void SwitchNode(ref int source, ref int target) {
            int buffer = source;

            source = target;
            target = buffer;
        }

        /*

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

        */
        /*

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

        */


        /// <summary>
        /// Builds and shares the boundary conditions of one point with all other points positioned at the same location
        /// </summary>
        private void BuildAndShareBoundaryConditions() {
            
            //int i, j, k, l, m, n;
            //int range_from;
            //int range_to;
            //int index;
            //int index2;
            //int edge_counter;
            //int verify;

            //bcfaces = new int[elements.Length, number_of_faces_per_element];
            //for (i = 0; i < bcfaces.GetLength(0); i++) {
            //    for (j = 0; j < bcfaces.GetLength(1); j++) {
            //        bcfaces[i, j] = -1;
            //    }
            //}
            //for (i = 0; i < bcelements.Length; i++) { // over boundary condition indices
            //    range_from = 0;
            //    range_to = bcelements[i].Length;

            //    if (pointset_types[i].Equals(PointSetType_t.ElementRange) || pointset_types[i].Equals(PointSetType_t.PointRange)) {
            //        range_from = bcelements[i][0];
            //        range_to = bcelements[i][1] + 1;
            //    }
            //    for (j = range_from; j < range_to; j++) { // over coordinate indices per boundary condition (if branch) or element indices per boundary condition (else branch)
            //        if (pointset_types[i].Equals(PointSetType_t.ElementRange) || pointset_types[i].Equals(PointSetType_t.PointRange)) {
            //            index = j;
            //        } else {
            //            index = bcelements[i][j];
            //        }
            //        if (gridlocation_types[i].Equals(GridLocation_t.Vertex)) {
            //            foreach (int z in node_to_element_indices[index]) { // over elements
            //                for (k = 0; k < faces[z].Length; k++) { // over face indices per element
            //                    if (k >= bcfaces.GetLength(1)) {
            //                        continue;
            //                    }
            //                    for (m = edge_counter = 0; m < faces[z][k].Length; m++) { // over coordinate indices per faces
            //                        for (n = range_from; n < range_to; n++) { // over coordinate indices per boundary condition
            //                            // Michael: Kann nicht ElementRange sein
            //                            if (/*pointset_types[i].Equals(PointSetType_t.ElementRange) ||*/ pointset_types[i].Equals(PointSetType_t.PointRange)) {
            //                                index2 = n;
            //                            } else {
            //                                index2 = bcelements[i][n];
            //                            }
            //                            if (index2 == faces[z][k][m]) {
            //                                edge_counter++;
            //                            }
            //                        }
            //                    }
            //                    if (edge_counter == number_of_face_edges) {
            //                        bcfaces[z, k] = i;
            //                    }
            //                }
            //            }
            //        } else {
            //            for (k = 0; k < faces[j].Length; k++) { // over face indices per element
            //                if (k >= bcfaces.GetLength(1)) {
            //                    continue;
            //                }
            //                bcfaces[j, k] = i;
            //            }
            //        }
            //    }
            //    for (j = 0; j < faces.Length; j++) { // over element indices
            //        for (k = 0; k < faces[j].Length; k++) { // over face indices per element
            //            foreach (int z in node_to_element_indices[faces[j][k][0]]) { // over elements
            //                if (z == j || !ChoosenBoSSSCompliantElementType(element_types[z])) {
            //                    continue;
            //                }
            //                for (l = 0; l < faces[z].Length; l++) { // over face indices per element
            //                    if (l >= bcfaces.GetLength(1)) {
            //                        continue;
            //                    }
            //                    for (m = edge_counter = 0; m < faces[z][l].Length; m++) { // over coordinate indices per faces
            //                        for (n = 0; n < faces[j][k].Length; n++) { // over coordinate indices per faces
            //                            if (faces[j][k][n] == faces[z][l][m]) {
            //                                edge_counter++;
            //                            }
            //                        }
            //                    }
            //                    if (edge_counter == number_of_face_edges) {
            //                        if (bcfaces[z, l] == -1) {
            //                            bcfaces[z, l] = bcfaces[j, k];
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //    }
            //}
            //if (bc_to_element_indices.Length == 0) {
            //    return;
            //}
            //for (i = 0; i < faces.Length; i++) { // over elements
            //    if (!ChoosenBoSSSCompliantElementType(element_types[i])) {
            //        continue;
            //    }
            //    for (j = 0; j < faces[i].Length; j++) { // over faces
            //        if ((verify = bcfaces[i, j]) >= 0) {
            //            bc_to_element_indices[verify].Add(i << number_of_faces_per_element | j);
            //        }
            //    }
            //}
        }


        /// <summary>
        /// Generates a grid commons object
        /// </summary>
        /// <returns>the grid commons object</returns>
        public GridCommons GenerateBoSSSGrid() {
            int j, k, l, i;
            int max_bc_index = -1;

            GenerateZonelessAndSectionlessArrays();
            CoordinateSystemTransformation();
            VertexIndexNumberingTransformation();
            //BuildFaces();
            //BuildAndShareBoundaryConditions();
            //BuildNeighbourship();

            double[, ,] vertices = new double[bosss_elements, verticesPerCell, number_of_coord_dimensions];
            long[] globalID = new long[bosss_elements];
            byte[,] edgeTags = new byte[bosss_elements, number_of_faces_per_element];

            for (l = 0, k = 0; l < bosss_elements; l++) {
                j = bosss_element_index_to_cgns_element_index[l];
                for (int n = 0; n < verticesPerCell; n++) {
                    k = cgns_element_index_to_bosss_element_index[j];

                    int node_no = elements[j][n];

                    vertices[k, n, 0] = x[node_no];
                    vertices[k, n, 1] = y[node_no];
                    if (number_of_coord_dimensions == 3) {
                        vertices[k, n, 2] = z[node_no];
                    }
                }
            }

            GridCommons grd = null;
            if (bosss_element_type.Equals(ElementType_t.QUAD_4)) {
                grd = new Grid2D(new Square());
            } else if (bosss_element_type.Equals(ElementType_t.TRI_3)) {
                grd = new Grid2D(new Triangle());
            } else if (bosss_element_type.Equals(ElementType_t.HEXA_8)) {
                throw new NotImplementedException();
            } else if (bosss_element_type.Equals(ElementType_t.TETRA_4)) {
                throw new NotImplementedException();
            } else {
                throw new NotSupportedException("Unsupported BoSSS element type");
            }

            int tag_number = 1;
            for (i = 0; i < bcelements.Length; i++) { // over boundary condition indices
                foreach (int z in bc_to_element_indices[i]) { // over elements
                    j = z >> number_of_faces_per_element;
                    j = cgns_element_index_to_bosss_element_index[j];
                    k = z & number_of_faces_per_element_mask;
                    edgeTags[j, k] = (byte)tag_number;
                    if (tag_number > max_bc_index) {
                        max_bc_index = tag_number;
                    }
                }
                tag_number++;
            }


            bool duplicatedBoundaryConditionNames = false;
            List<string> edgeTagNameList = new List<string>();

            for (i = 0; i < max_bc_index; i++) { // over boundary condition interfaces
                string edgeTagName = CleanString(new string(boconames[i]));

                if (edgeTagNameList.Contains(edgeTagName)) {
                    duplicatedBoundaryConditionNames = true;
                }
                edgeTagNameList.Add(edgeTagName);
            }

            tag_number = 1;
            List<string> edgeTagNameForbiddenList = new List<string>();
            List<int> numberList = new List<int>();
            if (duplicatedBoundaryConditionNames) {
                for (i = 0; i < edgeTagNameList.Count; i++) {
                    string edgeTagName = edgeTagNameList[i];

                    if (edgeTagNameForbiddenList.Contains(edgeTagName)) {
                        continue;
                    }
                    for (j = i + 1; j < edgeTagNameList.Count; j++) {
                        if (edgeTagNameList[j].Equals(edgeTagName)) {
                            replaceEdgeTagNumerInCaseOfDublicateBoundaryConditions(edgeTags, (byte)(j + 1), (byte)(i + 1));
                            max_bc_index--;
                            if (numberList.Contains(i + 1)) {
                                continue;
                            }
                            numberList.Add(i + 1);
                        }
                    }
                    edgeTagNameForbiddenList.Add(edgeTagName);
                }
                for (i = 0; i < max_bc_index; i++) { // over boundary condition interfaces
                    replaceEdgeTagNumerInCaseOfDublicateBoundaryConditions(edgeTags, (byte)numberList[i], (byte)(i + 1));
                }
                for (i = 0; i < max_bc_index; i++) { // over boundary condition interfaces
                    grd.EdgeTagsNames.Add((byte)tag_number, edgeTagNameForbiddenList[i]);
                    tag_number++;
                }
            } else {
                foreach (string edgeTagName in edgeTagNameList) {
                    grd.EdgeTagsNames.Add((byte)tag_number, edgeTagName);
                    tag_number++;
                }
            }

            /*

            // remove 'internal' boundaries ...
            {
                // pointvise may also assign boundaries for internal edges...
                // this is not allowed by bosss, so we have to clear them

                int E = cell_neighbours.GetLength(1);

                for (j = 0; j < bosss_elements; j++) {
                    for (int e = 0; e < E; e++) {
                        if (cell_neighbours[j, e] >= 0)
                            // this is an internal edge
                            if (edgeTags[j, e] != 0)
                                edgeTags[j, e] = 0;
                    }
                }

                // remove unused edgetag-names
                int[] EdgeTagUseCnt = new int[tag_number];
                for (j = 0; j < bosss_elements; j++) {
                    for (int e = 0; e < E; e++) {
                        int tag = edgeTags[j, e];
                        EdgeTagUseCnt[tag]++;
                    }
                }


                for (int tag = 0; tag < tag_number; tag++) {
                    if (EdgeTagUseCnt[tag] <= 0)
                        grd.EdgeTagsNames.Remove((byte)tag);
                }


                //// test whether all boundary edges have an assigned edge tag...
                //// (it's legal for a BoSSS grid to have unspecified edge tags)
                //for (j = 0; j < bosss_elements; j++) {
                //    for (int e = 0; e < E; e++) {
                //        int tag = edgeTags[j, e];
                //        long neigh = cell_neighbours[j,e];

                //        if (neigh < 0 && tag == 0)
                //            throw new ApplicationException("found boundary edge without assigned edge tag.");
                //    }
                //}

            } */

            grd.Description = "Created from BoSSS";

            grd.Cells = new Cell[bosss_elements];
            for (j = 0; j < bosss_elements; j++) {
                grd.Cells[j] = new Cell();
                grd.Cells[j].GlobalID = j;

                //grd.Cells[j].Neighbours = new long[number_of_faces_per_element];
                //grd.Cells[j].EdgeTags = new byte[number_of_faces_per_element];
                throw new NotImplementedException("todo");
                for (i = 0; i < number_of_faces_per_element; i++) {
                    //grd.Cells[j].Neighbours[i] = cell_neighbours[j, i];
                    //grd.Cells[j].EdgeTags[i] = edgeTags[j, i];
                }
            }

            //grd.Vertices = vertices;

            throw new NotImplementedException("todo");

            //return grd;
        }


        private static int NumberOfNodesPerElement(ElementType_t type) {
            char[] number = type.ToString().ToCharArray();
            return number[number.Length - 1] - 48;
        }


        private static bool Is3D(ElementType_t type) {
            return type.Equals(ElementType_t.TETRA_4) || type.Equals(ElementType_t.HEXA_8);
        }

        private static bool BoSSSCompliantElementType(ElementType_t type) {
            return type.Equals(ElementType_t.TRI_3) || type.Equals(ElementType_t.TETRA_4) || type.Equals(ElementType_t.QUAD_4) || type.Equals(ElementType_t.HEXA_8);
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
            char[] number = type.ToString().ToCharArray();
            return number[number.Length - 1] - 48;
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

    }

#pragma warning restore  1591
}

using System;
using System.Collections.Generic;
using System.Collections;
using System.Runtime.InteropServices;
using System.Text;

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.GridImport;
using BoSSS.Solution;

namespace BoSSS.Solution.CgnsExport {

    /// <summary>
    /// write DG fields to CGNS files
    /// </summary>
    public class CgnsExport : PlotDriver {

        /// <summary>
        /// Permutation tables for the conversion between Tecplot ordering and
        /// BoSSS ordering of the nodes in an element.
        /// </summary>
        static public Dictionary<ElementType_t, int[]> elementTypeToPermutationTableMap = new Dictionary<ElementType_t, int[]> {
            {ElementType_t.BAR_2, new int[] {0, 1}},
            {ElementType_t.TRI_3, new int[] {0, 1, 2}},
            {ElementType_t.QUAD_4, new int[] {0, 1, 3, 2}},
            {ElementType_t.TETRA_4, new int[] {0, 1, 2, 3}},
            {ElementType_t.HEXA_8, new int[] {0, 1, 4, 2, 3, 6, 5, 7}}
        };

        /// <summary>
        /// Implement this method by returning the file ending that should be
        /// appended at the end of a file name.
        /// </summary>
        protected override string FileEnding {
            get {
                if (hdf5) {
                    return "hdf5";
                } else {
                    return "cgns";
                }
            }
        }

        private IEnumerable<BoSSS.Foundation.Field> fields;
        private double time;
        private string title;

        /// <summary>
        /// Implement this methods by exporting/plotting all fields contained
        /// in <paramref name="fieldsToPlot"/> while regarding the settings
        /// <see cref="PlotDriver.superSampling"/>, <see cref="PlotDriver.ghostZone"/> and <see cref="PlotDriver.showJumps"/>.
        /// </summary>
        /// <param name="fileNameBase">
        /// The base name of the file (if any) that should be created. If
        /// called in parallel, rank information will be appended. If a file
        /// with the given name already exists, the resulting file name will
        /// contain an additional counter.
        /// </param>
        /// <param name="title">The title of the plot</param>
        /// <param name="time">
        /// The solution time represented by <paramref name="fieldsToPlot"/>
        /// </param>
        /// <param name="fieldsToPlot">
        /// A list of fields which should be plotted
        /// </param>
        public override void PlotFields(string fileNameBase, string title, double time, IEnumerable<BoSSS.Foundation.Field> fieldsToPlot) {
            string filename = GenerateFileName(fileNameBase);
            this.time = time;
            this.fields = fieldsToPlot;
            this.title = title;
            Build(filename, hdf5, base.GridDat);
        }

        /// <summary>
        /// Constructs a new PlotDriver
        /// </summary>
        /// <param name="context">The omnipresent context</param>
        /// <param name="hdf5"></param>
        /// <param name="showJumps">
        /// <param name="ghostZone"><see cref="PlotDriver"/></param>
        /// If true, the real DG data (including discontinuities) should be
        /// exported. Otherwise, the mean value of all values in one node
        /// should be calculated
        /// </param>
        /// <param name="superSampling"><see cref="PlotDriver.superSampling"/></param>
        public CgnsExport(GridData context, bool hdf5, bool showJumps, bool ghostZone, uint superSampling)
            : base(context, showJumps, ghostZone, superSampling,null) {
            this.hdf5 = hdf5;
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

        /// <summary>
        /// ture: CGNS-HDF5 is used; false: CGNS-ADF is used;
        /// </summary>
        public bool hdf5 {
            get;
            private set;
        }

        /// <summary>
        /// loads or saves an adf or hdf5 file from or to a specified path
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="hdf5"></param>
        /// <param name="context"></param>
        public void Build(string filePath, bool hdf5, GridData context) {
            using (CgnsDriver cg = new CgnsDriver(hdf5)) {
                int index_file = 0;
                int index_base = 0;
                int index_zone = 0;
                int index_section = 0;
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
                int[] irmin = new int[3];
                int[] irmax = new int[3];
                int celldim, physdim; // zone definition
                int start, end, nbndry; // element definition
                DataType_t data_type;


                // some header information
                // =======================

                basename = "Base".ToNullTermCharAry();
                zonename = "Zone".ToNullTermCharAry();
                coordnameX = "CoordinateX".ToNullTermCharAry();
                coordnameY = "CoordinateY".ToNullTermCharAry();
                coordnameZ = "CoordinateZ".ToNullTermCharAry();
                flow = "VectorMagnitude".ToNullTermCharAry();
                secname = "Section".ToNullTermCharAry();
                boconame = "BoundaryCondition".ToNullTermCharAry();
                solname = "FlowSolution".ToNullTermCharAry();
                data_type = DataType_t.RealDouble;

                int index_coord;
                int index_flow;
                int index_field;

                if (cg.open(filePath.ToNullTermCharAry(), (int)Mode.MODE_WRITE, out index_file) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }

                celldim = physdim = context.Grid.SpatialDimension;

                if (cg.base_write(index_file, basename, celldim, physdim, out index_base) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }

                cg.gopath(index_file, "/Base".ToNullTermCharAry());

                cg.descriptor_write("Information".ToNullTermCharAry(), title.ToNullTermCharAry());

                GridCommons grd = context.Grid;

                ElementType_t bosss_element_type = ElementType_t.ElementTypeNull;

                if (grd.RefElements.Length > 1)
                    throw new NotSupportedException("implementation todo");
                var RefElm = grd.RefElements[0];
                if (RefElm.GetType() == typeof(Line)) {
                    bosss_element_type = ElementType_t.BAR_2;
                } else if (RefElm.GetType() == typeof(Square)) {
                    bosss_element_type = ElementType_t.QUAD_4;
                } else if (RefElm.GetType() == typeof(Triangle)) {
                    bosss_element_type = ElementType_t.TRI_3;
                } else if (RefElm.GetType() == typeof(Cube)) {
                    bosss_element_type = ElementType_t.HEXA_8;
                } else if (RefElm.GetType() == typeof(Tetra)) {
                    bosss_element_type = ElementType_t.TETRA_4;
                } else {
                    ThrowError("Unsupported BoSSS element type", index_file, cg);
                }
                
                if (!showJumps) {
                    isize[0, 0] = totalVertices;
                } else {
                    isize[0, 0] = NoOfCells * verticesPerCell;
                }
                isize[0, 1] = totalCells;
                isize[0, 2] = 0;

                if (cg.zone_write(index_file, index_base, zonename, isize, ZoneType_t.Unstructured, out index_zone) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }

                // write grid vertices
                // ===================

                double[, ,] x = null;
                double[, ,] y = null;
                double[, ,] z = null;

                if (!showJumps) {
                    x = new double[isize[0, 0], 1, 1];
                    for (int vertex = 0; vertex < isize[0, 0]; vertex++) {
                        x[vertex, 0, 0] = vertices[vertex, 0];
                    }
                    if (cg.coord_write3(index_file, index_base, index_zone, data_type, coordnameX, x, out index_coord) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    y = new double[isize[0, 0], 1, 1];
                    for (int vertex = 0; vertex < isize[0, 0]; vertex++) {
                        y[vertex, 0, 0] = vertices[vertex, 1];
                    }
                    if (cg.coord_write3(index_file, index_base, index_zone, data_type, coordnameY, y, out index_coord) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (celldim > 2) {
                        z = new double[isize[0, 0], 1, 1];
                        for (int vertex = 0; vertex < isize[0, 0]; vertex++) {
                            z[vertex, 0, 0] = vertices[vertex, 2];
                        }
                        if (cg.coord_write3(index_file, index_base, index_zone, data_type, coordnameZ, z, out index_coord) != (int)Error.CG_OK) {
                            ThrowError(index_file, cg);
                        }
                    }
                } else {
                    x = new double[isize[0, 0], 1, 1];
                    int cnt = 0;
                    for (int j = 0; j < NoOfCells; j++) {
                        for (int i = 0; i < verticesPerCell; i++) {
                            x[cnt, 0, 0] = verticeCoordinates[j, i, 0];
                            cnt++;
                        }
                    }
                    if (cg.coord_write3(index_file, index_base, index_zone, data_type, coordnameX, x, out index_coord) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    y = new double[isize[0, 0], 1, 1];
                    cnt = 0;
                    for (int j = 0; j < NoOfCells; j++) {
                        for (int i = 0; i < verticesPerCell; i++) {
                            y[cnt, 0, 0] = verticeCoordinates[j, i, 1];
                            cnt++;
                        }
                    }
                    if (cg.coord_write3(index_file, index_base, index_zone, data_type, coordnameY, y, out index_coord) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (celldim > 2) {
                        z = new double[isize[0, 0], 1, 1];
                        cnt = 0;
                        for (int j = 0; j < NoOfCells; j++) {
                            for (int i = 0; i < verticesPerCell; i++) {
                                z[cnt, 0, 0] = verticeCoordinates[j, i, 2];
                                cnt++;
                            }
                        }
                        if (cg.coord_write3(index_file, index_base, index_zone, data_type, coordnameZ, z, out index_coord) != (int)Error.CG_OK) {
                            ThrowError(index_file, cg);
                        }
                    }
                }

                if (cg.sol_write(index_file, index_base, index_zone, solname, GridLocation_t.Vertex, out index_flow) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }

                // write DG fields
                // ===============

                foreach (Field field in fields) {
                    calculateField(field);

                    if (!showJumps) {
                        cg.field_write1(index_file, index_base, index_zone, index_flow, data_type, field.ToString().ToNullTermCharAry(), smoothedResult, out index_field);
                    } else {
                        cg.field_write1(index_file, index_base, index_zone, index_flow, data_type, field.ToString().ToNullTermCharAry(), notSmoothedResult, out index_field);
                    }
                }

                // write cell connectivity
                // =======================

                int[] permutationTable = elementTypeToPermutationTableMap[bosss_element_type];
                if (!showJumps) {
                    int[,] permutedConnectivity = new int[isize[0, 1], connectivity.GetLength(1)];

                    for (int i = 0; i < isize[0, 1]; i++) {
                        for (int j = 0; j < permutationTable.Length; j++) {
                            permutedConnectivity[i, j] = 1 + connectivity[i, permutationTable[j]];
                        }
                    }

                    start = 1;
                    end = isize[0, 1];
                    nbndry = 0;

                    if (cg.section_write2(index_file, index_base, index_zone, secname, bosss_element_type, start, end, nbndry, permutedConnectivity, out index_section) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }

                } else {
                    int[, ,] permutedConnectivity = new int[NoOfCells, subdivisionsPerCell, permutationTable.Length];
                    for (int j = 0; j < NoOfCells; j++) {
                        for (int jj = 0; jj < subdivisionsPerCell; jj++) {
                            for (int i = 0; i < permutationTable.Length; i++) {
                                permutedConnectivity[j, jj, i] = 1 + j * verticesPerCell + subdivisionTreeLeaves[jj].GlobalVerticeInd[permutationTable[i]];
                            }
                        }
                    }

                    start = 1;
                    end = isize[0, 1];
                    nbndry = 0;

                    if (cg.section_write3(index_file, index_base, index_zone, secname, bosss_element_type, start, end, nbndry, permutedConnectivity, out index_section) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                }
                if (cg.close(index_file) != (int)Error.CG_OK) {
                    ThrowError(cg);
                }
            }
        }
    }
}

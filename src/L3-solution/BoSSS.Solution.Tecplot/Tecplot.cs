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
using System.Runtime.InteropServices;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using System.Linq;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Solution.Tecplot {

    /// <summary>
    /// Tecplot plotting
    /// </summary>
    public class Tecplot : PlotDriver {

        private readonly string path = null;

        /// <summary>
        /// see <see cref="PlotDriver.PlotDriver"/>.
        /// </summary>
        public Tecplot(IGridData context, uint superSampling)
            : base(context, true, false, superSampling, null) {
        }

        /// <summary>
        /// see <see cref="PlotDriver.PlotDriver"/>.
        /// </summary>
        public Tecplot(IGridData context, bool showJumps, bool ghostZone, uint superSampling, CellMask mask = null)
            : base(context, showJumps, ghostZone, superSampling, mask) {
        }

        /// <summary>
        /// see <see cref="PlotDriver.PlotDriver"/>.
        /// </summary>
        /// <param name="path">path to output folder</param>
        public Tecplot(IGridData context, uint superSampling, string path)
            : this(context, superSampling) {
            this.path = path;
        }

        /// <summary>
        /// creates a zone driver for Tecplot;
        /// </summary>
        protected override PlotDriver.ZoneDriver CreateZoneDriver(IGridData context, int iKref, bool showJumps, bool showGhostCells, uint superSampling, CellMask __mask) {
            return new TecplotZone(context, iKref, showJumps, showGhostCells, superSampling, __mask);
        }

        /// <summary>
        /// Use ".plt" for all Tecplot files
        /// </summary>
        protected override string FileEnding {
            get {
                return "plt";
            }
        }

        /// <summary>
        /// Opens an output file by writing header information for Tecplot
        /// </summary>
        override protected void OpenFile(string filename, IEnumerable<string> fieldNames) {

            int dimension = base.gridData.SpatialDimension;

            StringWriter stw = new StringWriter();
            if (dimension == 1) {
                stw.Write("x ");
            } else if (dimension == 2) {
                stw.Write("x y ");
            } else if (dimension == 3) {
                stw.Write("x y z ");
            } else {
                throw new ApplicationException("No Tecplot for more than 3 dimensions.");
            }

            int unknowncnt = 1;
            List<string> writtenIDs = new List<string>();
            foreach (string _id in fieldNames) {
                var id = _id;
                if (id == null || id.Length <= 0) {
                    // write a name for field without identification
                    stw.Write("unknown_" + unknowncnt);
                    unknowncnt++;
                } else if (id.Contains(" ") || id.Contains("(") || id.Contains(")")) {
                    throw new Exception("A field identification contains white space or parentheses, please change the name of this field");
                } else if (id.Contains("\t") || id.Contains("\n") || id.Contains("\r")) {
                    // remove whitespaces
                    //string[] subStrg = id.Split(new char[] { ' ', '\t', '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries);

                    string[] subStrg = id.Split(new char[] { '\t', '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries);

                    string idNew = "";
                    for (int __i = 0; __i < subStrg.Length; __i++) {
                        idNew = idNew + subStrg[__i];
                        if (__i < subStrg.Length - 1)
                            idNew += "_";
                    }
                    id = idNew;
                }

                if (writtenIDs.Contains(id))
                    // ensure that name is unique
                    id = "_" + id;

                writtenIDs.Add(id);

                stw.Write(id);
                stw.Write(" ");
            }

            int Debug = 0;
            int VIsDouble = 1;
            string ScratchDir = ".";
            string Variables = stw.ToString();
            string filenameWithPath = path != null ? Path.Combine(path, filename) : filename;

            IntPtr ptrTitle, ptrVariables, ptrFName, ptrScratchDir;
            ptrTitle = Marshal.StringToHGlobalAnsi(filename);
            ptrScratchDir = Marshal.StringToHGlobalAnsi(ScratchDir);
            ptrFName = Marshal.StringToHGlobalAnsi(filenameWithPath);
            ptrVariables = Marshal.StringToHGlobalAnsi(Variables);

            int errorWhileOpening = tecini110(ptrTitle, ptrVariables, ptrFName, ptrScratchDir, ref Debug, ref VIsDouble);
            if (errorWhileOpening == -1)
            {
                throw new Exception("Tecplot could not create file. Do you have writing permission?");
            }
            Marshal.FreeHGlobal(ptrTitle);
            Marshal.FreeHGlobal(ptrScratchDir);
            Marshal.FreeHGlobal(ptrFName);
            Marshal.FreeHGlobal(ptrVariables);
        }

        /// <summary>
        /// Closes the currently open output file 
        /// </summary>
        override protected void CloseFile() {
            tecend110();
        }

        /// <summary>
        /// Tecplot Zone
        /// </summary>
        public class TecplotZone : ZoneDriver {

            /// <summary>
            /// The different finite-element zone types supported by Tecplot. The
            /// assigned byte codes correspond to the Tecplot notation
            /// </summary>
            private enum ZoneType {

                /// <summary>
                /// One-Dimensional line segment (two nodes)
                /// </summary>
                FELineSeg = 1,

                /// <summary>
                /// Two-dimensional triangular element (three nodes)
                /// </summary>
                FETriangle = 2,

                /// <summary>
                /// Two-dimensional quadrilateral element (four nodes)
                /// </summary>
                FEQuad = 3,

                /// <summary>
                /// Three-dimensional tetrahedral element (four nodes)
                /// </summary>
                FETetrahedron = 4,

                /// <summary>
                /// Three-dimensional brick element (eight nodes)
                /// </summary>
                FEBrick = 5
            }

            /// <summary>
            /// Permutation tables for the conversion between Tecplot ordering and
            /// BoSSS ordering of the nodes in an element.
            /// </summary>
            private Dictionary<ZoneType, int[]> zoneTypeToPermutationTableMap = new Dictionary<ZoneType, int[]> {
                {ZoneType.FELineSeg, new int[] {0, 1}},
                {ZoneType.FETriangle, new int[] {0, 1, 2}},
                {ZoneType.FEQuad, new int[] {0, 1, 3, 2}},
                {ZoneType.FETetrahedron, new int[] {0, 1, 2, 3}},
                {ZoneType.FEBrick, new int[] {3, 6, 1, 0, 7, 5, 4, 2}},
            };

            /// <summary>
            /// ctor.
            /// </summary>
            public TecplotZone(IGridData context, int iKref, bool showJumps, bool ghostZone, uint superSampling, CellMask mask = null)
                : base(context, iKref, showJumps, ghostZone, superSampling, mask) {
            }

            /// <summary>
            /// Plots all fields in <paramref name="fieldsToPlot"/> into a Tecplot file, see
            /// <see cref="PlotDriver.ZoneDriver.PlotZone"/>.
            /// </summary>
            public override void PlotZone(string ZoneName, double time, IEnumerable<Tuple<string, ScalarFunctionEx>> fieldsToPlot) {
                List<ScalarFunctionEx> fields = new List<ScalarFunctionEx>(
                    fieldsToPlot.Select(x => x.Item2));


                //for (int iKref = 0; iKref < GridDat.Grid.RefElements.Length; iKref++) {
                {

                    int totalVertices = vertices.GetLength(0);
                    if (showJumps) {
                        totalVertices = NoOfCells * verticesPerCell;
                    }
                    ZoneType zoneType = GetZoneType(base.Zone_Element);
                    WriteZoneInformation(time, zoneType, totalVertices, totalCells, ZoneName);

                    // write grid vertices
                    // ===================

                    if (!showJumps) {
                        // reconstructed mode 
                        // vertices that are at the same position are written only once

                        double[] globalCoordinates = new double[totalVertices];

                        for (int d = 0; d < dimension; d++) {
                            for (int vertex = 0; vertex < totalVertices; vertex++) {
                                globalCoordinates[vertex] = vertices[vertex, d];
                            }

                            int VIsDouble = 1;
                            int n = globalCoordinates.Length;
                            tecdat110(ref n, globalCoordinates, ref VIsDouble);
                        }
                    } else {
                        // each-cell-on-its-own -- mode

                        double[] globalCoordinates = new double[verticesPerCell];

                        for (int d = 0; d < dimension; d++) {
                            for (int j = 0; j < NoOfCells; j++) {
                                for (int i = 0; i < verticesPerCell; i++) {
                                    globalCoordinates[i] = verticeCoordinates[j, i, d];
                                }

                                int VIsDouble = 1;
                                int n = globalCoordinates.Length;
                                tecdat110(ref n, globalCoordinates, ref VIsDouble);
                            }
                        }
                    }

                    // write DG field data
                    // ===================
                    foreach (ScalarFunctionEx field in fields) {
                        int VIsDouble = 1;

                        SampleField(field, showJumps);

                        if (!showJumps) {
                            tecdat110(ref totalVertices, smoothedResult, ref VIsDouble);
                        } else {
                            tecdat110(ref totalVertices, notSmoothedResult, ref VIsDouble);
                        }
                    }

                    // write cells
                    // ===========

                    int[] permutationTable = zoneTypeToPermutationTableMap[zoneType];

                    if (!showJumps) {
                        int[,] permutedConnectivity = new int[totalCells, connectivity.GetLength(1)];

                        for (int i = 0; i < totalCells; i++) {
                            for (int j = 0; j < permutationTable.Length; j++) {
                                permutedConnectivity[i, j] = 1 + connectivity[i, permutationTable[j]];
                            }
                        }
                        tecnod110(permutedConnectivity);
                    } else {
                        int[,,] permutedConnectivity = new int[NoOfCells, subdivisionsPerCell, permutationTable.Length];
                        for (int j = 0; j < NoOfCells; j++) {
                            for (int jj = 0; jj < subdivisionsPerCell; jj++) {
                                for (int i = 0; i < permutationTable.Length; i++) {
                                    permutedConnectivity[j, jj, i] = 1 + j * verticesPerCell + subdivisionTreeLeaves[jj].GlobalVerticeInd[permutationTable[i]];
                                }
                            }
                        }
                        tecnod110(permutedConnectivity);
                    }

                }
            }

            /// <summary>
            /// Write the zone information. Since we always write exactly one zone,
            /// only the <paramref name="numberOfPoints"/> and the
            /// <paramref name="numberOfElements"/> differ between the
            /// reconstructed and the continuous case.
            /// </summary>
            /// <param name="time">
            /// <see cref="PlotDriver.PlotFields(string, string, double, IEnumerable{DGField})"/>
            ///</param>
            /// <param name="zoneType"><see cref="Tecplot.TecplotZone.ZoneType"/></param>
            /// <param name="numberOfPoints">
            /// The total number of vertices in the plot
            /// </param>
            /// <param name="numberOfElements">
            /// The total number of cells in the plot
            /// </param>
            /// <param name="zone_name">
            /// arbitrary naming for the zone
            /// </param>
            private void WriteZoneInformation(double time, ZoneType zoneType, int numberOfPoints, int numberOfElements, string zone_name) {
                IntPtr ptrZoneTitle = Marshal.StringToHGlobalAnsi(zone_name);
                int zoneTypeIndex = (int)zoneType;
                int KMax = 0, ICellMax = 0, JCellMax = 0, KCellMax = 0, NFConn = 0, FNMode = 0;
                int IsBlock = 1;
                int StrandID = 0, ParentZone = 0, ShrConn = 0;

                teczne110(ptrZoneTitle,
                    ref zoneTypeIndex,
                    ref numberOfPoints,
                    ref numberOfElements,
                    ref KMax,
                    ref ICellMax,
                    ref JCellMax,
                    ref KCellMax,
                    ref time,
                    ref StrandID,
                    ref ParentZone,
                    ref IsBlock,
                    ref NFConn,
                    ref FNMode,
                    null,           // No passive variables
                    null,           // All variables node-centered
                    null,           // No shared variables
                    ref ShrConn);

                Marshal.FreeHGlobal(ptrZoneTitle);
            }

            /// <summary>
            /// Defines a mapping between BoSSS elements (<see cref="RefElement"/>)
            /// and the Tecplot representations.
            /// </summary>
            /// <returns>The Tecplot zone type for the current grid</returns>
            private ZoneType GetZoneType(RefElement Kref) {

                ZoneType zoneType;
                if (Kref.GetType() == typeof(Line)) {
                    zoneType = ZoneType.FELineSeg;
                } else if (Kref.GetType() == typeof(Square)) {
                    zoneType = ZoneType.FEQuad;
                } else if (Kref.GetType() == typeof(Triangle)) {
                    zoneType = ZoneType.FETriangle;
                } else if (Kref.GetType() == typeof(Cube)) {
                    zoneType = ZoneType.FEBrick;
                } else if (Kref.GetType() == typeof(Tetra)) {
                    zoneType = ZoneType.FETetrahedron;
                } else {
                    throw new NotSupportedException(Kref.GetType().ToString() + " - simplex currently not supported in tecplot.");
                }
                return zoneType;
            }
        }

        /// <summary>
        /// Legacy interface
        /// </summary>
        public static void PlotFields(IEnumerable<DGField> _FieldsToPlot, string filename, double time, int supersampling) {
            Tecplot tecplot = new Tecplot((_FieldsToPlot.First().GridDat), true, false, (uint)supersampling, null);
            tecplot.PlotFields(filename, time, _FieldsToPlot);
        }

        /*
        /// <summary>
        /// Legacy interface
        /// </summary>
        public static void PlotFields(IEnumerable<DGField> _FieldsToPlot, GridData context, string filename, double time, int supersampling) {
            Tecplot tecplot = new Tecplot(context, true, false, (uint)supersampling, null);
            tecplot.PlotFields(filename, title, time, _FieldsToPlot);
        }
        */

        /// <summary>
        /// Initializes the process of writing a binary data file. This must be
        /// called first before any other TecIO calls are made (except
        /// TECFOREIGN110). You may write to multiple files by calling
        /// TECINI110 more than once. Each time TECINI110 is called, a new file
        /// is opened. Use TECFIL110 to switch between files. For each call to
        /// TECINI, there must be a corresponding call to TECEND110
        /// </summary>
        /// <param name="Title">
        /// Title of the data set
        /// </param>
        /// <param name="Variables">
        /// List of variable names. If a comma appears in the string it will be
        /// used as the separator between variable names, otherwise a space is
        /// used.
        /// </param>
        /// <param name="FName">Name of the file to create</param>
        /// <param name="ScratchDir">
        /// Name of the directory to put the scratch file.
        /// </param>
        /// <param name="Debug">
        /// Pointer to the integer flag for debugging. Set to 0 for no
        /// debugging or 1 to debug. When set to 1, the debug messages will be
        /// sent to the standard output (stdout).
        /// </param>
        /// <param name="VIsDouble">
        /// Pointer to the integer flag for specifying whether field data
        /// generated in future calls to TECDAT110 are to be written in single
        /// or double precision.
        /// </param>
        /// <returns>0 if successful, -1 if unsuccessful.</returns>
        [DllImport("tecio")]
        private static unsafe extern int tecini110(
            IntPtr Title,
            IntPtr Variables,
            IntPtr FName,
            IntPtr ScratchDir,
            ref int Debug,
            ref int VIsDouble);

        /// <summary>
        /// Writes header information about the next zone to be added to the
        /// data file. After TECZNE110 is called, you must call TECDAT110 one
        /// or more times. If the zone is a finite-element zone, call TECNOD110
        /// (cell-based zones) or TECPOLY111 (face-based zones) after calling
        /// TECDAT110
        /// </summary>
        /// <param name="ZoneTitle">
        /// The name of the zone
        /// </param>
        /// <param name="ZoneType">
        /// 0=ORDERED
        /// 1=FELINESEG
        /// 2=FETRIANGLE
        /// 3=FEQUADRILATERAL
        /// 4=FETETRAHEDRON
        /// 5=FEBRICK
        /// 6=FEPOLYGON
        /// 7=FEPOLYHEDRON
        /// </param>
        /// <param name="IMxOrNumPts">
        /// For ordered zones, the number of nodes in the Iindex direction. For
        /// finite-element zones (cellbased and face-based), the number of
        /// nodes.
        /// </param>
        /// <param name="JMxOrNumElements">
        /// For ordered zones, the number of nodes in the J index direction.
        /// For finite-element zones (cellbased and face-based), the number of
        /// elements.
        /// </param>
        /// <param name="KMx">
        /// For ordered zones, the number of nodes in the K index direction.
        /// For polyhedral and polygonal finite-element zones, it is the number
        /// of faces. Not used all other finite-element zone types
        /// </param>
        /// <param name="ICellMax">Reserved for future use. Set to zero</param>
        /// <param name="JCellMax">Reserved for future use. Set to zero.</param>
        /// <param name="KCellMax">Reserved for future use. Set to zero.</param>
        /// <param name="SolutionTime">
        /// Scalar double precision value specifying the time associated with
        /// the zone.
        /// </param>
        /// <param name="StrandID">
        /// Scalar integer value specifying the strand to which the zone is
        /// associated.
        /// 0 = static zone, not associated with a strand.
        /// Values greater than 0 indicate a zone is assigned to a given strand
        /// </param>
        /// <param name="ParentZone">
        /// Scalar integer value representing the relationship between this
        /// zone and its parent. With a parent zone association, Tecplot can
        /// generate a surface streamtrace on a no-slip boundary zone. A zone
        /// may not specify itself as its own parent.
        /// 0 = indicates that this zone is not associated with a parent zone.
        /// >0 = A value greater than zero is considered this zone's parent.
        /// </param>
        /// <param name="IsBlock">
        /// Indicates whether the data will be passed into TECDAT110 in BLOCK
        /// or POINT format.
        /// 0=POINT
        /// 1=BLOCK
        /// </param>
        /// <param name="NumFaceConnections">
        /// Used for cell-based finite-element and ordered zones only. The
        /// number of face connections that will be passed in routine
        /// TECFACE110
        /// </param>
        /// <param name="FaceNeighborMode">
        /// Used for cell-baseda finite-element and ordered zones only. The
        /// type of face connections that will be passed in routine TECFACE110.
        /// 0=LocalOneToOne
        /// 2=GlobalOneToOne
        /// 1=LocalOneToMany
        /// 3=GlobalOneToMany
        /// </param>
        /// <param name="PassiveVarList">
        /// Array, dimensioned by the number of variables, of 4 byte integer
        /// values specifying the active/passive nature of each variable. A
        /// value of 0 indicates the associated variable is active while a
        /// value of 1 indicates that it is passive.
        /// </param>
        /// <param name="ValueLocation">
        /// The location of each variable in the data set. ValueLocation(I)
        /// indicates the location of variable I for this zone.
        /// 0=cell-centered
        /// 1=node-centered.
        /// Pass null to indicate that all variables are nodecentered.
        /// </param>
        /// <param name="ShareVarFromZone">
        /// Indicates variable sharing. Array, dimensioned by the number of
        /// variables. ShareVarFromZone(I) indicates the zone number with which
        /// variable I will be shared. This reduces the amount of data to be
        /// passed via TECDAT110. A value of 0 indicates that the variable is
        /// not shared. Pass null to indicate no variable sharing for this
        /// zone. You must pass null for the first zone in a data set (there is
        /// no data available to share).
        /// </param>
        /// <param name="ShareConnectivityFromZone">
        /// Indicates the zone number with which connectivity is shared. Pass 0
        /// to indicate no connectivity sharing. You must pass 0 for the first
        /// zone in a data set.
        /// NOTE: Connectivity and/or face neighbors cannot be shared when the
        /// face neighbor mode is set to Global. Connectivity cannot be shared
        /// between cell-based and face-based finite-element zones.
        /// </param>
        /// <returns>0 if successful, -1 if unsuccessful.</returns>
        [DllImport("tecio")]
        private static unsafe extern int teczne110(
            IntPtr ZoneTitle,
            ref int ZoneType,
            ref int IMxOrNumPts,
            ref int JMxOrNumElements,
            ref int KMx,
            ref int ICellMax,
            ref int JCellMax,
            ref int KCellMax,
            ref double SolutionTime,
            ref int StrandID,
            ref int ParentZone,
            ref int IsBlock,
            ref int NumFaceConnections,
            ref int FaceNeighborMode,
            int[] PassiveVarList,
            int[] ValueLocation,
            int[] ShareVarFromZone,
            ref int ShareConnectivityFromZone);

        /// <summary>
        /// Writes an array of data to the data file. Data should not be passed
        /// for variables that have been indicated as passive or shared (via
        /// TECZNE110). TECDAT110 allows you to write your data in a piecemeal
        /// fashion in case it is not contained in one contiguous block in your
        /// program. TECDAT110 must be called enough times to ensure that the
        /// correct number of values are written for each zone and that the
        /// aggregate order for the data is correct.
        /// </summary>
        /// <param name="N">
        /// Pointer to an integer value specifying number of values to write.
        /// </param>
        /// <param name="Data">
        /// Array of single or double precision data values
        /// </param>
        /// <param name="IsDouble">
        /// Pointer to the integer flag stating whether the array Data is
        /// single (0) or double (1) precision.
        /// </param>
        /// <returns>0 if successful, -1 if unsuccessful.</returns>
        [DllImport("tecio")]
        private static unsafe extern int tecdat110(ref int N, double[] Data, ref int IsDouble);

        /// <summary>
        /// Must be called to close out the current data file. There must be
        /// one call to TECEND110 for each TECINI111.
        /// </summary>
        /// <returns>0 if successful, -1 if unsuccessful.</returns>
        [DllImport("tecio")]
        private static unsafe extern int tecend110();

        /// <summary>
        /// Writes an array of node data to the binary data file. This is the
        /// connectivity list for cell-based finite-element zones (line
        /// segment, triangle, quadrilateral, brick, and tetrahedral zones)
        /// </summary>
        /// <param name="nodelist">
        /// <see cref="TecplotZone.PlotZone"/>
        /// </param>
        /// <returns>0 if successful, -1 if unsuccessful.</returns>
        [DllImport("tecio")]
        private static unsafe extern int tecnod110(int[,,] nodelist);

        /// <summary>
        /// Writes an array of node data to the binary data file. This is the
        /// connectivity list for cell-based finite-element zones (line
        /// segment, triangle, quadrilateral, brick, and tetrahedral zones)
        /// </summary>
        /// <param name="nodelist">
        /// Array of integers listing the nodes for each element. This is the
        /// connectivity list, dimensioned (T, M) (T moving fastest), where M
        /// is the number of elements in the zone and T is set according to the
        /// following list:
        /// 2=Line Segment
        /// 3=Triangle
        /// 4=Tetrahedral
        /// 4=Quadrilateral
        /// 8=Brick
        /// </param>
        /// <returns>0 if successful, -1 if unsuccessful.</returns>
        [DllImport("tecio")]
        private static unsafe extern int tecnod110(int[,] nodelist);


    }
}

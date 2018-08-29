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
using System.Globalization;
using System.IO;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.ASCIIExport {

    /// <summary>
    /// Export driver for VisIts curve file format for one-dimensional data.
    /// </summary>
    /// <remarks>
    /// An exemplary file looks as follows:
    /// <example>
    /// |file format|
    /// |:-:|
    /// |# variableName1|
    /// |x1  y1|
    /// |x2  y2|
    /// | ...  |
    /// |xN  yN|
    /// |# variableName2|
    /// |...|
    /// </example>
    /// </remarks>
    public class CurveExportDriver : PlotDriver {

        /// <summary>
        /// see <see cref="PlotDriver.PlotDriver"/>.
        /// </summary>
        public CurveExportDriver(IGridData context, bool showJumps, bool ghostZone, uint superSampling, CellMask sgrd = null)
            : base(context, showJumps, ghostZone, superSampling, sgrd) //
        {
            if (context.SpatialDimension != 1) {
                throw new ArgumentException("Only supported in one dimension", "context");
            }
        }

        /// <summary>
        /// see <see cref="PlotDriver.PlotDriver"/>.
        /// </summary>
        public CurveExportDriver(IGridData context, bool showJumps, uint superSampling)
            : base(context, showJumps, false, superSampling, null) {
            if (context.SpatialDimension != 1) {
                throw new ArgumentException("Only supported in one dimension", "context");
            }
        }

        StreamWriter file;

        /// <summary>
        /// closes file, e.g disposes StreamWriter
        /// </summary>
        protected override void CloseFile() {
            file.Dispose();
        }

        /// <summary>
        /// generates StreamWriter
        /// </summary>
        protected override void OpenFile(string filename, IEnumerable<string> fieldNames) {
            file = new StreamWriter(filename);
            m_fieldNames = fieldNames.ToArray();
        }

        string[] m_fieldNames;

        /// <summary>
        /// creates new CurveExportZone
        /// </summary>
        protected override PlotDriver.ZoneDriver CreateZoneDriver(IGridData context, int iKref, bool showJumps, bool showGhostCells, uint superSampling, CellMask __sgrd) {
            return new CurveExportZone(context, this, iKref, showJumps, showGhostCells, superSampling, __sgrd);
        }

        /// <summary>
        /// File ending is "curve"
        /// </summary>
        protected override string FileEnding {
            get {
                return "curve";
            }
        }
        /// <summary>
        /// Legacy interface
        /// </summary>
        public static void PlotFields(IEnumerable<DGField> _FieldsToPlot, GridData context, string filename, string title, double time, int supersampling) {
            CurveExportDriver CurveExport = new CurveExportDriver(context, true, false, (uint)supersampling, null);
            CurveExport.PlotFields(filename, time, _FieldsToPlot);
        }

        /// <summary>
        /// CSVExport Zone
        /// </summary>
        public class CurveExportZone : ZoneDriver {

            private CurveExportDriver owner;

            /// <summary>
            /// ctor.
            /// </summary>
            public CurveExportZone(IGridData context, CurveExportDriver owner, int iKref, bool showJumps, bool ghostZone, uint superSampling, CellMask sgrd = null)
                : base(context, iKref, showJumps, ghostZone, superSampling, sgrd) {
                this.owner = owner;
            }

            /// <summary>
            /// <see cref="PlotDriver.ZoneDriver.PlotZone"/>
            /// </summary>
            public override void PlotZone(string ZoneName, double time, IEnumerable<Tuple<string, ScalarFunctionEx>> _fieldsToPlot) {
                var fieldsToPlot = _fieldsToPlot.ToArray();
                Debug.Assert(owner.m_fieldNames.Length == fieldsToPlot.Length);
                
                for(int ii = 0; ii < fieldsToPlot.Length; ii++) {
                    var field = fieldsToPlot[ii];

                    SampleField(field.Item2, showJumps);

                    owner.file.WriteLine("# " + owner.m_fieldNames[ii]);

                    Tuple<double, double>[] tuples;

                    if (showJumps) {
                        int noOfNodesPerCell = verticeCoordinates.GetLength(1);
                        tuples = new Tuple<double, double>[NoOfCells * noOfNodesPerCell];
                        for (int i = 0; i < NoOfCells; i++) {
                            for (int j = 0; j < noOfNodesPerCell; j++) {
                                int index = i * noOfNodesPerCell + j;
                                tuples[index] = Tuple.Create(verticeCoordinates[i, j, 0], notSmoothedResult[index]);
                            }
                        }
                    } else {
                        tuples = new Tuple<double, double>[totalVertices];
                        for (int i = 0; i < totalVertices; i++) {
                            tuples[i] = Tuple.Create(vertices[i, 0], smoothedResult[i]);
                        }
                    }

                    // Write tuples with increasing x values
                    foreach (var tuple in tuples.OrderBy(t => t.Item1)) {
                        owner.file.WriteLine(
                            tuple.Item1.ToString("e16", NumberFormatInfo.InvariantInfo) +
                            "\t" +
                            tuple.Item2.ToString("e16", NumberFormatInfo.InvariantInfo));
                    }
                }
            }
        }
    }
}

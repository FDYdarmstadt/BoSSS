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
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;

namespace BoSSS.Solution.ASCIIExport {

    /// <summary>
    /// Export driver for generic ASCII files. In this particular
    /// implementation, tabs will be used as separators and numbers will be
    /// formatted in the language-independent (see
    /// <see cref="NumberFormatInfo.InvariantInfo"/>) scientific format. May be
    /// useful for debugging and one-dimensional plots (e.g., using gnuplot)
    /// </summary>
    /// <remarks>
    /// An exemplary output file for the two-dimensional case looks as follows:
    /// <example>
    /// | x | y | variable1 | variable2 |
    /// |:-:|:-:|:---------:|:---------:|
    /// | 0 | 0 |     1     |     2     |
    /// | 1 | 0 |     1     |     3     |  
    /// | 0 | 1 |     2     |     2     |
    /// | 1 | 1 |     2     |     3     | 
    /// | 1 | 1 |     3     |     3     |
    /// </example>
    /// As shown in the example, identical coordinates may appear multiple times
    /// if jumps are activated (see
    /// <see cref="CSVExportDriver.CSVExportDriver(GridData, bool, bool, uint, SubGrid)"/>).
    /// </remarks>
    public class CSVExportDriver : PlotDriver {

        /// <summary>
        /// <see cref="NumberFormatInfo.InvariantInfo"/>
        /// </summary>
        private static readonly IFormatProvider formatInfo
            = NumberFormatInfo.InvariantInfo;

        /// <summary>
        /// see <see cref="PlotDriver.PlotDriver"/>.
        /// </summary>
        public CSVExportDriver(GridData context, bool showJumps, uint superSampling)
            : base(context, showJumps, false, superSampling, null) {
        }

        /// <summary>
        /// see <see cref="PlotDriver.PlotDriver"/>.
        /// </summary>
        public CSVExportDriver(GridData context, bool showJumps, bool ghostZone, uint superSampling, SubGrid sgrd = null)
            : base(context, showJumps, ghostZone, superSampling, sgrd) {
        }

        StreamWriter file;

        /// <summary>
        /// closes file, e.g disposes StreamWriter
        /// </summary>
        protected override void CloseFile() {
            file.Dispose();
        }

        /// <summary>
        /// Creates StreamWriter and writes header of the data file
        /// </summary>
        protected override void OpenFile(string filename, IEnumerable<string> fieldNames) {
            file = new StreamWriter(filename);
            switch (base.gridData.SpatialDimension) {
                case 1:
                    file.Write("x");
                    break;

                case 2:
                    file.Write("x\ty");
                    break;

                case 3:
                    file.Write("x\ty\tz");
                    break;

                default:
                    throw new Exception();
            }

            foreach (var field in fieldNames) {
                file.Write("\t" + field);
            }
            file.WriteLine();
        }

        /// <summary>
        /// creates new CSVExportZone
        /// </summary>
        protected override PlotDriver.ZoneDriver CreateZoneDriver(GridData context, int iKref, bool showJumps, bool showGhostCells, uint superSampling, SubGrid __sgrd) {
            return new CSVExportZone(context, this, iKref, showJumps, showGhostCells, superSampling, __sgrd);
        }

        /// <summary>
        /// File ending is "csv"
        /// </summary>
        protected override string FileEnding {
            get {
                return "csv";
            }
        }

        /// <summary>
        /// Legacy interface
        /// </summary>
        public static void PlotFields(IEnumerable<DGField> _FieldsToPlot, GridData context, string filename, string title, double time, int supersampling) {
            CSVExportDriver CSVExport = new CSVExportDriver(context, true, false, (uint)supersampling, null);
            CSVExport.PlotFields(filename, time, _FieldsToPlot);
        }

        /// <summary>
        /// CSVExport Zone
        /// </summary>
        public class CSVExportZone : ZoneDriver {

            private CSVExportDriver owner;

            /// <summary>
            /// ctor.
            /// </summary>
            public CSVExportZone(GridData context, CSVExportDriver owner, int iKref, bool showJumps, bool ghostZone, uint superSampling, SubGrid sgrd = null)
                : base(context, iKref, showJumps, ghostZone, superSampling, sgrd) {
                this.owner = owner;
            }

            /// <summary>
            /// <see cref="PlotDriver.ZoneDriver.PlotZone"/>
            /// </summary>
            public override void PlotZone(string ZoneName, double time, IEnumerable<Tuple<string, ScalarFunctionEx>> fieldsToPlot) {
                List<double[]> fieldValues = new List<double[]>(fieldsToPlot.Count());
                foreach (var field in fieldsToPlot) {
                    SampleField(field.Item2, showJumps);

                    if (showJumps) {
                        fieldValues.Add(notSmoothedResult.CloneAs());
                    } else {
                        fieldValues.Add(smoothedResult.CloneAs());
                    }
                }

                // Item1: Point coordinates
                // Item2: Field values
                Tuple<double[], double[]>[] tuples;

                if (showJumps) {
                    int noOfNodesPerCell = verticeCoordinates.GetLength(1);
                    tuples = new Tuple<double[], double[]>[NoOfCells * noOfNodesPerCell];
                    // Generates tuples with all field values in tuple.Item2
                    for (int i = 0; i < NoOfCells; i++) {
                        for (int j = 0; j < noOfNodesPerCell; j++) {
                            int index = i * noOfNodesPerCell + j;
                            double[] values = new double[fieldValues.Count];

                            int k = 0;
                            foreach (double[] valueList in fieldValues) {
                                values[k] = valueList[index];
                                k++;
                            }

                            double[] point = new double[dimension];
                            for (int d = 0; d < dimension; d++) {
                                point[d] = verticeCoordinates[i, j, d];
                            }
                            tuples[index] = Tuple.Create(point, values);
                        }
                    }
                } else {
                    tuples = new Tuple<double[], double[]>[totalVertices];
                    for (int i = 0; i < totalVertices; i++) {
                        double[] values = new double[fieldValues.Count];

                        int k = 0;
                        foreach (double[] valueList in fieldValues) {
                            values[k] = valueList[i];
                            k++;
                        }

                        double[] point = new double[dimension];
                        for (int d = 0; d < dimension; d++) {
                            point[d] = vertices[i, d];
                        }

                        tuples[i] = Tuple.Create(point, values);
                    }
                }

                // Write tuples with increasing x values
                foreach (var tuple in tuples.OrderBy(t => t.Item1[0])) {
                    owner.file.Write(tuple.Item1[0].ToString("e16", NumberFormatInfo.InvariantInfo));
                    for (int d = 1; d < dimension; d++) {
                        owner.file.Write("\t" + tuple.Item1[d].ToString("e16", NumberFormatInfo.InvariantInfo));
                    }
                    foreach (double value in tuple.Item2) {
                        owner.file.Write("\t" + value.ToString("e16", NumberFormatInfo.InvariantInfo));
                    }
                    owner.file.WriteLine();
                }
            }
        }
    }
}

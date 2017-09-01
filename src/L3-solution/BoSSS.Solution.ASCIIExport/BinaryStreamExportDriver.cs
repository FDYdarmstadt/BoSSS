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

using ilPSP;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.ASCIIExport {

    /// <summary>
    /// Export driver for binary streams of double values such as they are used
    /// by Matlab
    /// </summary>
    public class BinaryStreamExportDriver : PlotDriver {

        /// <summary>
        /// File writers for the different variables to be exported (cf.
        /// <see cref="GetOutputFieldNames(IEnumerable{string})"/>)
        /// </summary>
        private Dictionary<string, BinaryWriter> fileWriters;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="context"></param>
        /// <param name="showJumps"></param>
        /// <param name="showGhostCells"></param>
        /// <param name="superSampling"></param>
        /// <param name="sgrd"></param>
        public BinaryStreamExportDriver(GridData context, bool showJumps, bool showGhostCells, uint superSampling, SubGrid sgrd = null)
            : base(context, showJumps, showGhostCells, superSampling, sgrd) {
        }

        /// <summary>
        /// Constructs a zone driver with the given options
        /// </summary>
        /// <param name="gridData"></param>
        /// <param name="iKref"></param>
        /// <param name="showJumps"></param>
        /// <param name="showGhostCells"></param>
        /// <param name="superSampling"></param>
        /// <param name="sgrd"></param>
        /// <returns></returns>
        protected override PlotDriver.ZoneDriver CreateZoneDriver(GridData gridData, int iKref, bool showJumps, bool showGhostCells, uint superSampling, SubGrid sgrd) {
            return new BinaryZoneDriver(gridData, this, iKref, showJumps, showGhostCells, superSampling, sgrd);
        }

        /// <summary>
        /// Opens one file per output variable. Existing files are overwritten
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="fieldNames"></param>
        protected override void OpenFile(string filename, IEnumerable<string> fieldNames) {
            string[] outputFields = GetOutputFieldNames(fieldNames);

            fileWriters = new Dictionary<string, BinaryWriter>(outputFields.Length);
            for (int i = 0; i < outputFields.Length; i++) {
                string specificFileName = String.Format(filename, outputFields[i]);

                fileWriters[outputFields[i]] = new BinaryWriter(
                    new FileStream(specificFileName, FileMode.Create),
                    Encoding.ASCII);
            }
        }

        /// <summary>
        /// File name format is $baseName.$rank.$variableName.$fileEnding
        /// </summary>
        /// <param name="fileNameBase"></param>
        /// <returns></returns>
        protected override string GenerateFileName(string fileNameBase) {
            int size = gridData.MpiSize;

            string infix = "";
            if (size > 1) {
                infix = "." + (gridData.MpiRank + 1) + "of" + size;
            }

            return fileNameBase + infix + ".{0}." + FileEnding;
        }

        /// <summary>
        /// Returns all fields within <paramref name="fieldNames"/> plus the
        /// <see cref="Coordinates"/>
        /// </summary>
        /// <param name="fieldNames"></param>
        /// <returns></returns>
        private string[] GetOutputFieldNames(IEnumerable<string> fieldNames) {
            if (fieldNames.Any(name => Coordinates.Contains(name))) {
                throw new Exception(
                    "Field names 'x', 'y' and 'z' are forbidden");
            }

            return Coordinates.Concat(fieldNames).ToArray();
        }

        /// <summary>
        /// Closes all files within <see cref="fileWriters"/>
        /// </summary>
        protected override void CloseFile() {
            foreach (BinaryWriter writer in fileWriters.Values) {
                writer.Close();
                writer.Dispose();
            }
        }

        /// <summary>
        /// ".bin"
        /// </summary>
        protected override string FileEnding {
            get {
                return "bin";
            }
        }

        private string[] coordinates;

        /// <summary>
        /// { "x" }, { "x", "y" } or { "x", "y", "z" }, depending on the
        /// spatial dimension
        /// </summary>
        protected string[] Coordinates {
            get {
                if (coordinates == null) {
                    coordinates = new[] { "x", "y", "z" }.Take(gridData.SpatialDimension).ToArray();
                }
                return coordinates;
            }
        }

        /// <summary>
        /// Export for each individual zone
        /// </summary>
        public class BinaryZoneDriver : ZoneDriver {

            private BinaryStreamExportDriver owner;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="context"></param>
            /// <param name="owner"></param>
            /// <param name="iKref"></param>
            /// <param name="showJumps"></param>
            /// <param name="ghostZone"></param>
            /// <param name="superSampling"></param>
            /// <param name="sgrd"></param>
            public BinaryZoneDriver(GridData context, BinaryStreamExportDriver owner, int iKref, bool showJumps, bool ghostZone, uint superSampling, SubGrid sgrd = null)
                : base(context, iKref, showJumps, ghostZone, superSampling, sgrd) {
                this.owner = owner;
            }

            /// <summary>
            /// Export the values of <paramref name="fieldsToPlot"/>, while the
            /// one file is used per variable.
            /// </summary>
            /// <param name="ZoneName"></param>
            /// <param name="time"></param>
            /// <param name="fieldsToPlot"></param>
            public override void PlotZone(string ZoneName, double time, IEnumerable<Tuple<string, ScalarFunctionEx>> fieldsToPlot) {
                Dictionary<string, double[]> fieldValues = new Dictionary<string, double[]>(fieldsToPlot.Count());
                foreach (var field in fieldsToPlot) {
                    SampleField(field.Item2, showJumps);

                    if (showJumps) {
                        fieldValues.Add(field.Item1, notSmoothedResult.CloneAs());
                    } else {
                        fieldValues.Add(field.Item1, smoothedResult.CloneAs());
                    }
                }

                // Write everything into tuples first in order so support
                // sorting of the points later on
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
                            foreach (double[] valueList in fieldValues.Values) {
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
                        foreach (double[] valueList in fieldValues.Values) {
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

                // Sort tuples by x, then by y, then by z
                var sortedTuples = tuples.OrderBy(t => t.Item1[0]);
                for (int d = 1; d < dimension; d++) {
                    // Copy d to dummy because loop variable will be captured
                    // and incremented otherwise
                    int dummy = d;
                    sortedTuples = sortedTuples.ThenBy(t => t.Item1[dummy]);
                }

                // Write coordinates and variable values
                foreach (var tuple in sortedTuples) {
                    for (int d = 0; d < dimension; d++) {
                        owner.fileWriters[owner.Coordinates[d]].Write(tuple.Item1[d]);
                    }

                    for (int f = 0; f < tuple.Item2.Length; f++) {
                        owner.fileWriters[fieldValues.ElementAt(f).Key].Write(tuple.Item2[f]);
                    }
                }
            }
        }
    }
}

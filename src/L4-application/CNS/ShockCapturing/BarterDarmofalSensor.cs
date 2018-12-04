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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CNS.ShockCapturing {

    public class BarterDarmofalSensor : IShockSensor {

        private Variable sensorVariable;

        private DGField sensorValues;

        public BarterDarmofalSensor(Variable sensorVariable) {
            this.sensorVariable = sensorVariable;
        }

        public double GetSensorValue(int cellIndex) {
            return sensorValues.GetMeanValue(cellIndex);
        }

        public void UpdateSensorValues(IEnumerable<DGField> fieldSet, ISpeciesMap speciesMap, CellMask cellMask) {
            using (new FuncTrace()) {
                DGField fieldToTest = fieldSet.
                    Where(f => f.Identification == sensorVariable.Name).
                    Single();

                sensorValues = new SinglePhaseField(
                    new Basis(fieldToTest.GridDat, 0));

                var quadrature = EdgeQuadrature.GetQuadrature(
                    new int[] { 1 },
                    (GridData)fieldToTest.GridDat,
                    new EdgeQuadratureScheme().Compile(fieldToTest.GridDat, 2 * fieldToTest.Basis.Degree),
                    delegate (int e0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        MultidimensionalArray gIn = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                        MultidimensionalArray gOut = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                        fieldToTest.EvaluateEdge(
                            e0,
                            Length,
                            QR.Nodes,
                            gIn,
                            gOut,
                            MeanValueIN: null,
                            MeanValueOT: null,
                            GradientIN: null,
                            GradientOT: null,
                            ResultIndexOffset: 0,
                            ResultPreScale: 0.0);

                        for (int e = 0; e < Length; e++) {
                            int edge = e + e0;

                        // Check if "out" neighbor exist; ignore boundary edges for now
                        int outCell = fieldToTest.GridDat.iLogicalEdges.CellIndices[edge, 1];
                            if (outCell < 0) {
                                continue;
                            }

                            for (int node = 0; node < QR.NoOfNodes; node++) {
                                double jump = gOut[e, node] - gIn[e, node];
                                double mean = 0.5 * (gOut[e, node] + gIn[e, node]);
                                double s = jump / mean;
                                if (!double.IsNaN(s)) {
                                    EvalResult[e, node, 0] = s * s / (Math.Sign(s) * s + 0.001);
                                }
                            }
                        }
                    },
                    delegate (int e0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int e = 0; e < Length; e++) {
                            int edge = e + e0;
                            for (int neighbor = 0; neighbor < 2; neighbor++) { // loop over IN/OUT cells
                            int jCell = fieldToTest.GridDat.iLogicalEdges.CellIndices[edge, neighbor];
                                if (jCell < 0) {
                                    break;
                                }

                                sensorValues.Coordinates[jCell, 0] += ResultsOfIntegration[e, 0];
                            }
                        }
                    });
                quadrature.Execute();
            }
        }
    }
}

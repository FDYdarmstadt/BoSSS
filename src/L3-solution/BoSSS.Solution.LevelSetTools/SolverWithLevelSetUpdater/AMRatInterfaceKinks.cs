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
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// refinement on cells which show strong kinks (alternating) in the interface shape
    /// </summary>
    [Serializable]
    public class AMRatInterfaceKinks : AMRLevelIndicatorWithLevelset {
        public override int[] DesiredCellChanges() {

            var grddDat = LsTrk.GridDat;
            int J = grddDat.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            double[] gradDiff = new double[J];      // stores per cell the absolute diff values of each neighboring edge with interface 

            int order = LevSet.Basis.Degree;
            var metrics = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS, order);
            var factory = metrics.XQuadFactoryHelper.GetSurfaceElement_BoundaryRuleFactory(0, grddDat.Grid.RefElements[0]);
            CellMask CC = LsTrk.Regions.GetCutCellMask();
            EdgeQuadratureScheme InnerBoundaryEdgeSchemme = new EdgeQuadratureScheme(factory, CC.GetAllInnerEdgesMask());

            EdgeQuadrature.GetQuadrature(new int[] { 1 }, grddDat,
                InnerBoundaryEdgeSchemme.Compile(grddDat, order),
                delegate (int i0, int length, QuadRule QR, MultidimensionalArray EvalResult) {

                    MultidimensionalArray GradientIN = MultidimensionalArray.Create(length, QR.NoOfNodes, grddDat.SpatialDimension);
                    MultidimensionalArray GradientOUT = MultidimensionalArray.Create(length, QR.NoOfNodes, grddDat.SpatialDimension);
                    LevSet.EvaluateEdge(i0, length, QR.Nodes, null, null, null, null, GradientIN, GradientOUT);

                    for(int d = 0; d < grddDat.SpatialDimension; d++) {
                        MultidimensionalArray absDiff = MultidimensionalArray.Create(length, QR.NoOfNodes);
                        absDiff.Acc(1.0, GradientOUT.ExtractSubArrayShallow(-1, -1, d));
                        absDiff.Acc(-1.0, GradientIN.ExtractSubArrayShallow(-1, -1, d));
                        absDiff.ApplyAll(delegate (int[] index, ref double entry) { entry = entry.Abs(); });
                        EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(1.0, absDiff);
                    }

                },
                delegate (int i0, int length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < length; i++) {
                        gradDiff[grddDat.Edges.CellIndices[i0 + i, 0]] += ResultsOfIntegration[i, 0];
                        gradDiff[grddDat.Edges.CellIndices[i0 + i, 1]] += ResultsOfIntegration[i, 0];
                    }
                }
            ).Execute();

            SinglePhaseField AbsDiffField = new SinglePhaseField(new Basis(grddDat, 0), "absEdgeGradientDiff");
            for(int j = 0; j < J; j++) {
                AbsDiffField.SetMeanValue(j, gradDiff[j]);
            }
            Tecplot.Tecplot.PlotFields(new DGField[] { AbsDiffField, LevSet }, "InterfaceKinks", 0.0, 2);

            for(int j = 0; j < J; j++) {
                if(gradDiff[j] > 1e-4 && levels[j] < maxRefinementLevel) { 
                    levels[j] += 1;
                }
                if(gradDiff[j] == 0.0 && levels[j] > 0) {
                    levels[j] -= 1;
                }
            }

            return levels;
        }
    }

}

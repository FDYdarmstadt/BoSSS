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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.EquationSystem;
using CNS.IBM;
using ilPSP;
using System;

namespace CNS.Diffusion {

    /// <summary>
    /// Builder for the <see cref="SIPGFlux"/>.
    /// </summary>
    public class SIPGFluxBuilder : FluxBuilder {

        /// <summary>
        /// Information about the grid.
        /// </summary>
        private IGridData gridData;

        private Func<MultidimensionalArray> cellMetricFunc;

        /// <summary>
        /// Constructs a new flux builder.
        /// </summary>
        /// <param name="control"></param>
        /// <param name="boundaryMap"></param>
        /// <param name="speciesMap"></param>
        /// <param name="gridData"></param>
        public SIPGFluxBuilder(CNSControl control, BoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, IGridData gridData)
            : base(control, boundaryMap, speciesMap) {
            this.gridData = gridData;

            //Create Functions for calculation the cell metric, needed as Func<> because
            //LevelSet field and HMF options are not known at this point
            if (speciesMap is IBM.ImmersedSpeciesMap) {
                // IBM case
                ImmersedSpeciesMap IBMspeciesMap = speciesMap as ImmersedSpeciesMap;
                SpeciesId species = IBMspeciesMap.Tracker.GetSpeciesId(IBMspeciesMap.Control.FluidSpeciesName);
                cellMetricFunc = delegate () {
                    MultidimensionalArray cellMetric = IBMspeciesMap.CellAgglomeration.CellLengthScales[species].CloneAs();
                    cellMetric.ApplyAll(x => 1 / x);
                    // Needed, because 1/x produces NaN in void cells and can happen that penalty factor leads then to NaN
                    cellMetric.ApplyAll(delegate (double x) {
                        if (double.IsNaN(x) || double.IsInfinity(x)) {
                            return 0;
                        } else {
                            return x;
                        }
                    });
                    return cellMetric;
                };
            } else {
                // Non-IBM
                MultidimensionalArray cj = ((GridData)gridData).Cells.cj;
                cellMetricFunc = () => cj;
            }

        }

        /// <summary>
        /// Adds instances of appropriate sub-classes of <see cref="SIPGFlux"/>
        /// to the given <paramref name="op"/>
        /// </summary>
        /// <param name="op"></param>
        public override void BuildFluxes(Operator op) {

            op.DensityComponents.Add(
                new SIPGDensityFlux(control, boundaryMap, speciesMap, gridData, 0, cellMetricFunc));

            for (int i = 0; i < CompressibleEnvironment.NumberOfDimensions; i++) {
                op.MomentumComponents[i].Add(
                    new SIPGMomentumFlux(control, boundaryMap, speciesMap, gridData, i + 1, cellMetricFunc));
            }

            op.EnergyComponents.Add(
                new SIPGEnergyFlux(control, boundaryMap, speciesMap, gridData, CompressibleEnvironment.NumberOfDimensions + 1, cellMetricFunc));
        }
    }
}

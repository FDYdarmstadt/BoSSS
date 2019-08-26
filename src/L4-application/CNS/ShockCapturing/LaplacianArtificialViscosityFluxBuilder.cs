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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using CNS.EquationSystem;

namespace CNS.ShockCapturing {

    class LaplacianArtificialViscosityFluxBuilder : FluxBuilder {

        public LaplacianArtificialViscosityFluxBuilder(CNSControl control, CompressibleBoundaryCondMap boundaryMap, ISpeciesMap speciesMap)
            : base(control, boundaryMap, speciesMap) {
        }

        public override void BuildFluxes(Operator mapping) {
            GridData gridData = (GridData)this.speciesMap.GridData;

            mapping.DensityComponents.Add(new OptimizedLaplacianArtificialViscosityFlux(
                gridData,
                CompressibleVariables.Density.Name,
                penaltySafetyFactor: 1.0,
                penaltyFactor: (control.DensityDegree + 1) * (control.DensityDegree + gridData.SpatialDimension) / gridData.SpatialDimension,
                inverseLengthScales: gridData.Cells.cj
                ));

            for (int d = 0; d < gridData.SpatialDimension; d++) {
                mapping.MomentumComponents[d].Add(new OptimizedLaplacianArtificialViscosityFlux(
                    gridData,
                    CompressibleVariables.Momentum[d].Name,
                    penaltySafetyFactor: 1.0,
                    penaltyFactor: (control.MomentumDegree + 1) * (control.MomentumDegree + gridData.SpatialDimension) / gridData.SpatialDimension,
                    inverseLengthScales: gridData.Cells.cj
                    ));
            }

            mapping.EnergyComponents.Add(new OptimizedLaplacianArtificialViscosityFlux(
                 gridData,
                 CompressibleVariables.Energy.Name,
                    penaltySafetyFactor: 1.0,
                    penaltyFactor: (control.EnergyDegree + 1) * (control.EnergyDegree + gridData.SpatialDimension) / gridData.SpatialDimension,
                    inverseLengthScales: gridData.Cells.cj
                 ));


            //// Old artificial viscosity flux (hack for testing different AV boundary conditions)
            //mapping.DensityComponents.Add(new LaplacianArtificialViscosityFlux(
            //    control.DensityDegree, this.speciesMap.GridData.Cells.cj, CompressibleEnvironment.PrimalArgumentOrdering, CompressibleEnvironment.PrimalArgumentToIndexMap[Variables.Density], boundaryMap, speciesMap));

            //for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
            //    mapping.MomentumComponents[d].Add(new LaplacianArtificialViscosityFlux(
            //        control.MomentumDegree, this.speciesMap.GridData.Cells.cj, CompressibleEnvironment.PrimalArgumentOrdering, CompressibleEnvironment.PrimalArgumentToIndexMap[Variables.Momentum[d]], boundaryMap, speciesMap));
            //}
            //mapping.EnergyComponents.Add(new LaplacianArtificialViscosityFlux(
            //    control.EnergyDegree, this.speciesMap.GridData.Cells.cj, CompressibleEnvironment.PrimalArgumentOrdering, CompressibleEnvironment.PrimalArgumentToIndexMap[Variables.Energy], boundaryMap, speciesMap));


            // Old artificial viscosity flux
            //mapping.DensityComponents.Add(new LaplacianArtificialViscosityFlux(
            //    control.DensityDegree, this.speciesMap.GridData.Cells.cj, Variables.Density, boundaryMap, speciesMap));

            //for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
            //    mapping.MomentumComponents[d].Add(new LaplacianArtificialViscosityFlux(
            //        control.MomentumDegree, this.speciesMap.GridData.Cells.cj, Variables.Momentum[d], boundaryMap, speciesMap));
            //}

            //mapping.EnergyComponents.Add(new LaplacianArtificialViscosityFlux(
            //    control.EnergyDegree, this.speciesMap.GridData.Cells.cj, Variables.Energy, boundaryMap, speciesMap));
        }
    }
}

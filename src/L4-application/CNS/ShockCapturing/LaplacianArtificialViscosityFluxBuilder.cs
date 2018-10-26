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

using BoSSS.Solution.CompressibleFlowCommon;
using CNS.Boundary;
using CNS.EquationSystem;

namespace CNS.ShockCapturing {

    class LaplacianArtificialViscosityFluxBuilder : FluxBuilder {

        public LaplacianArtificialViscosityFluxBuilder(CNSControl control, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap)
            : base(control, boundaryMap, speciesMap) {
        }

        public override void BuildFluxes(Operator mapping) {
            //double hMin = this.speciesMap.GridData.Cells.h_minGlobal;

            // Optimized artificial viscosity flux
            mapping.DensityComponents.Add(new OptimizedLaplacianArtificialViscosityFlux(
                (BoSSS.Foundation.Grid.Classic.GridData) this.speciesMap.GridData, 
                control.DensityDegree, this.speciesMap.GridData.SpatialDimension, Variables.Density.Name));

            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                mapping.MomentumComponents[d].Add(new OptimizedLaplacianArtificialViscosityFlux(
                    (BoSSS.Foundation.Grid.Classic.GridData) this.speciesMap.GridData, 
                    control.MomentumDegree, this.speciesMap.GridData.SpatialDimension, Variables.Momentum[d].Name));
            }   
            
            mapping.EnergyComponents.Add(new OptimizedLaplacianArtificialViscosityFlux(
                 (BoSSS.Foundation.Grid.Classic.GridData) this.speciesMap.GridData, 
                 control.EnergyDegree, this.speciesMap.GridData.SpatialDimension, Variables.Energy.Name));


            //// Old artificial viscosity flux (hack for testing different AV boundary conditions)
            //mapping.DensityComponents.Add(new LaplacianArtificialViscosityFlux(
            //    control.DensityDegree, this.speciesMap.GridData.Cells.cj, CNSEnvironment.PrimalArgumentOrdering, CNSEnvironment.PrimalArgumentToIndexMap[Variables.Density], boundaryMap, speciesMap));

            //for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
            //    mapping.MomentumComponents[d].Add(new LaplacianArtificialViscosityFlux(
            //        control.MomentumDegree, this.speciesMap.GridData.Cells.cj, CNSEnvironment.PrimalArgumentOrdering, CNSEnvironment.PrimalArgumentToIndexMap[Variables.Momentum[d]], boundaryMap, speciesMap));
            //}
            //mapping.EnergyComponents.Add(new LaplacianArtificialViscosityFlux(
            //    control.EnergyDegree, this.speciesMap.GridData.Cells.cj, CNSEnvironment.PrimalArgumentOrdering, CNSEnvironment.PrimalArgumentToIndexMap[Variables.Energy], boundaryMap, speciesMap));


            // Old artificial viscosity flux
            //mapping.DensityComponents.Add(new LaplacianArtificialViscosityFlux(
            //    control.DensityDegree, this.speciesMap.GridData.Cells.cj, Variables.Density, boundaryMap, speciesMap));

            //for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
            //    mapping.MomentumComponents[d].Add(new LaplacianArtificialViscosityFlux(
            //        control.MomentumDegree, this.speciesMap.GridData.Cells.cj, Variables.Momentum[d], boundaryMap, speciesMap));
            //}

            //mapping.EnergyComponents.Add(new LaplacianArtificialViscosityFlux(
            //    control.EnergyDegree, this.speciesMap.GridData.Cells.cj, Variables.Energy, boundaryMap, speciesMap));
        }
    }
}

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
using CNS.Boundary;
using CNS.EquationSystem;
using System;

namespace CNS.Diffusion {

    /// <summary>
    /// Supported types of diffusive fluxes
    /// </summary>
    public enum DiffusiveFluxTypes {

        /// <summary>
        /// No diffusive flux
        /// </summary>
        None = 0,

        /// <summary>
        /// Symmetric interior penalty Discontinuous Galerkin
        /// </summary>
        SIPG,

        /// <summary>
        /// Optimized version of SIPG
        /// </summary>
        OptimizedSIPG
    }

    /// <summary>
    /// Extension methods for <see cref="DiffusiveFluxTypes"/>
    /// </summary>
    public static class DiffusiveFluxTypesExtensions {

        /// <summary>
        /// Instantiates the <see cref="FluxBuilder"/> associated with the
        /// given <paramref name="diffusiveFlux"/>.
        /// </summary>
        /// <param name="diffusiveFlux">
        /// The selected flux type.
        /// </param>
        /// <param name="control">
        /// <see cref="FluxBuilder.FluxBuilder"/>
        /// </param>
        /// <param name="boundaryMap">
        /// <see cref="FluxBuilder.FluxBuilder"/>
        /// </param>
        /// <param name="speciesMap">
        /// <see cref="FluxBuilder.FluxBuilder"/>
        /// </param>
        /// <param name="gridData">
        /// Grid information; e.g. required for the calculation of penalty
        /// parameters
        /// </param>
        /// <returns>
        /// An instance of <see cref="FluxBuilder"/> that constructs the fluxes
        /// corresponding to <paramref name="diffusiveFlux"/>.
        /// </returns>
        public static FluxBuilder GetBuilder(this DiffusiveFluxTypes diffusiveFlux, CNSControl control, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, GridData gridData) {
            int minDegree = Math.Min(
                Math.Min(control.DensityDegree, control.MomentumDegree),
                control.EnergyDegree);

            switch (diffusiveFlux) {
                case DiffusiveFluxTypes.SIPG:
                    if (minDegree < 1) {
                        throw new Exception(
                            "SIPG is only valid for DG degrees greater than 0");
                    }
                    return new SIPGFluxBuilder(control, boundaryMap, speciesMap, gridData);

                case DiffusiveFluxTypes.OptimizedSIPG:
                    if (minDegree < 1) {
                        throw new Exception(
                            "SIPG is only valid for DG degrees greater than 0");
                    }
                    return new OptimizedSIPGFluxBuilder(control, boundaryMap, speciesMap, gridData);

                case DiffusiveFluxTypes.None:
                    return NullFluxBuilder.Instance;

                default:
                    throw new Exception("Unknown flux function \"" + diffusiveFlux + "\"");
            }
        }
    }
}

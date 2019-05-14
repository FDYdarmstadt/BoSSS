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
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using CNS.EquationSystem;
using CNS.IBM;
using System;

namespace CNS.Convection {

    /// <summary>
    /// Supported types of convective fluxes
    /// </summary>
    public enum ConvectiveFluxTypes {

        /// <summary>
        /// No convective flux
        /// </summary>
        None = 0,

        /// <summary>
        /// See <see cref="RusanovFluxBuilder"/>.
        /// </summary>
        Rusanov,

        /// <summary>
        /// See <see cref="HLLFluxBuilder"/>
        /// </summary>
        HLL,

        /// <summary>
        /// See <see cref="HLLCFluxBuilder"/>
        /// </summary>
        HLLC,

        /// <summary>
        /// See <see cref="OptimizedHLLCFluxBuilder"/>
        /// </summary>
        OptimizedHLLC,

        /// <summary>
        /// See <see cref="GodunovFluxBuilder"/>
        /// </summary>
        Godunov,

        MovingFrameRusanov
    }

    /// <summary>
    /// Extension methods for <see cref="ConvectiveFluxTypes"/>
    /// </summary>
    public static class ConvectiveFluxTypesExtensions {

        /// <summary>
        /// Instantiates the correct sub-class of <see cref="FluxBuilder"/>
        /// corresponding to the selected <paramref name="flux"/>.
        /// </summary>
        /// <param name="flux">
        /// The flux for which the builder should be instantiated.
        /// </param>
        /// <param name="control">
        /// Configuration options
        /// </param>
        /// <param name="boundaryMap">
        /// Information about boundary conditions
        /// </param>
        /// <param name="speciesMap">
        /// Mapping of different species int he domain
        /// </param>
        /// <returns>
        /// An instance of a flux builder that builds fluxes
        /// corresponding to the given <paramref name="flux"/>.
        /// </returns>
        public static FluxBuilder GetBuilder(this ConvectiveFluxTypes flux, CNSControl control, BoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, IEquationOfState equationOfState, double machNumber) {
            switch (flux) {
                case ConvectiveFluxTypes.Rusanov:
                    return new RusanovFluxBuilder(control, boundaryMap, speciesMap, machNumber);

                case ConvectiveFluxTypes.HLL:
                    return new HLLFluxBuilder(control, boundaryMap, speciesMap, machNumber);

                case ConvectiveFluxTypes.HLLC:
                    return new HLLCFluxBuilder(control, boundaryMap, speciesMap, machNumber);

                case ConvectiveFluxTypes.OptimizedHLLC:
                    return new OptimizedHLLCFluxBuilder(control, boundaryMap, speciesMap, equationOfState, machNumber);

                case ConvectiveFluxTypes.Godunov:
                    return new GodunovFluxBuilder(control, boundaryMap, speciesMap, machNumber);

                case ConvectiveFluxTypes.MovingFrameRusanov:
                    return new MovingFrameRusanovFluxBuilder(control, boundaryMap, speciesMap, machNumber); 

                case ConvectiveFluxTypes.None:
                    return NullFluxBuilder.Instance;

                default:
                    throw new Exception("Unknown flux function \"" + flux + "\"");
            }
        }
    }
}

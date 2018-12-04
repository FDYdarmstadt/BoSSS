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

using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;

namespace BoSSS.Solution.CompressibleFlowCommon.Convection {

    /// <summary>
    /// Common interface for all components of the Euler equations (i.e.,
    /// density, momentum and energy)
    /// </summary>
    public interface IEulerEquationComponent {

        /// <summary>
        /// Returns the conserved variable associated to the represented
        /// component of the Euler equations (i.e., density, momentum or
        /// energy).
        /// </summary>
        /// <param name="state">
        /// The state vector from which the variable will be computed.
        /// </param>
        /// <returns>
        /// Depending on the implementation: Density, momentum or energy
        /// </returns>
        double VariableValue(StateVector state);

        /// <summary>
        /// Returns the flux (<b>not</b> the numerical flux) associated to the
        /// represented component of the Euler equations.
        /// </summary>
        /// <param name="state">
        /// The state vector from which the flux will be computed.
        /// </param>
        /// <returns>
        /// The flux part of the represented equation.
        /// </returns>
        /// <remarks>
        /// The term <i>flux</i> must not be mistaken for the <i>numerical</i>
        /// flux. This is implemented in <see cref="EulerFlux"/>.
        /// </remarks>
        Vector Flux(StateVector state);
    }
}

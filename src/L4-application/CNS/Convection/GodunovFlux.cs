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
using CNS.Boundary;
using CNS.MaterialProperty;
using System;

namespace CNS.Convection {

    /// <summary>
    /// Flux based on the exact Riemann solver by Toro.
    /// </summary>
    /// <seealso cref="ExactRiemannSolver"/>
    public class GodunovFlux : EulerFlux {

        /// <summary>
        /// <see cref="EulerFlux.EulerFlux"/>
        /// </summary>
        /// <param name="config">
        /// <see cref="EulerFlux.EulerFlux"/>
        /// </param>
        /// <param name="boundaryMap">
        /// <see cref="EulerFlux.EulerFlux"/>
        /// </param>
        /// <param name="equationComponent">
        /// <see cref="EulerFlux.EulerFlux"/>
        /// </param>
        /// <param name="speciesMap">
        /// <see cref="EulerFlux.EulerFlux"/>
        /// </param>
        public GodunovFlux(CNSControl config, IBoundaryConditionMap boundaryMap, IEulerEquationComponent equationComponent, ISpeciesMap speciesMap)
            : base(config, boundaryMap, equationComponent, speciesMap) {
            if (config.EquationOfState is IdealGas == false) {
                throw new Exception("Riemann solver currently only supports ideal gases");
            }
        }

        /// <summary>
        /// Uses <see cref="ExactRiemannSolver"/> in order to compute the exact
        /// solution of the Riemann problem.
        /// </summary>
        /// <param name="x">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector3D, int)"/>
        /// </param>
        /// <param name="time">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector3D, int)"/>
        /// </param>
        /// <param name="stateIn">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector3D, int)"/>
        /// </param>
        /// <param name="stateOut">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector3D, int)"/>
        /// </param>
        /// <param name="normal">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector3D, int)"/>
        /// </param>
        /// <param name="edgeIndex">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector3D, int)"/>
        /// </param>
        /// <returns>
        /// <see cref="ExactRiemannSolver.GetCentralState"/>
        /// </returns>
        protected internal override double InnerEdgeFlux(double[] x, double time, StateVector stateIn, StateVector stateOut, ref Vector3D normal, int edgeIndex) {
            ExactRiemannSolver riemannSolver = new ExactRiemannSolver(
                stateIn, stateOut, normal);

            StateVector stateEdge = riemannSolver.GetCentralState();
            return equationComponent.Flux(stateEdge) * normal;
        }
    }
}

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

using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using System;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;

namespace BoSSS.Solution.CompressibleFlowCommon.Convection {

    /// <summary>
    /// Flux based on the exact Riemann solver by Toro.
    /// </summary>
    /// <seealso cref="ExactRiemannSolver"/>
    public class GodunovFlux : EulerFlux {

        /// <summary>
        /// <see cref="EulerFlux.EulerFlux"/>
        /// </summary>
        public GodunovFlux(CompressibleControl config, IBoundaryConditionMap boundaryMap, IEulerEquationComponent equationComponent, Material material)
            : base(config, boundaryMap, equationComponent, material) {
            if (config.EquationOfState is IdealGas == false) {
                throw new Exception("Riemann solver currently only supports ideal gases");
            }
        }

        /// <summary>
        /// Uses <see cref="ExactRiemannSolver"/> in order to compute the exact
        /// solution of the Riemann problem.
        /// </summary>
        public override double InnerEdgeFlux(double[] x, double time, StateVector stateIn, StateVector stateOut, ref ilPSP.Vector normal, int edgeIndex) {
            ExactRiemannSolver riemannSolver = new ExactRiemannSolver(
                stateIn, stateOut, normal);

            StateVector stateEdge = riemannSolver.GetCentralState();
            return equationComponent.Flux(stateEdge) * normal;
        }

    }
}

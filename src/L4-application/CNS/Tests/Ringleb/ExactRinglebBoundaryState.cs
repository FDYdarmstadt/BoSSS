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
using BoSSS.Solution.CompressibleFlowCommon.Boundary;

namespace CNS.Tests.Ringleb {

    /// <summary>
    /// A boundary condition that represents the exact solution of the Ringleb
    /// test case.
    /// </summary>
    public class ExactRinglebBoundaryState : BoundaryCondition {

        private RinglebControl control;

        /// <summary>
        /// Constructs a new boundary state
        /// </summary>
        /// <param name="control"></param>
        public ExactRinglebBoundaryState(RinglebControl control)
            : base(control.GetMaterial()) {
            this.control = control;
        }

        /// <summary>
        /// Uses <see cref="RinglebExactSolution"/> in order to construct the
        /// boundary state corresponding to the exact solution of the Ringleb
        /// problem.
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="stateIn"></param>
        /// <returns></returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            double kappa = control.EquationOfState.HeatCapacityRatio;
            double pi = control.EquationOfState.ReferencePressure;
            double a0 = control.RinglebReferenceSpeedOfSound;
            double p0 = control.RinglebReferenceTotalPressure;
            RinglebExactSolution.FlowState ringlebState = RinglebExactSolution.GetFlowState(
                x[0], x[1], kappa, pi, a0, p0);

            double rho = ringlebState.Density;
            Vector v = new Vector(ringlebState.Velocity[0], ringlebState.Velocity[1], 0.0);
            return new StateVector(
                stateIn.Material,
                rho,
                rho * v,
                (ringlebState.Pressure + kappa * pi) / (kappa - 1.0) + 0.5 * rho * v * v);
        }
    }
}

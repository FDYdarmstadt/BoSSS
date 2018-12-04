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

using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using CNS.MaterialProperty;
using static BoSSS.Solution.CompressibleFlowCommon.Boundary.ExactRinglebBoundaryState;

namespace CNS.Tests.Ringleb {

    /// <summary>
    /// Specialized control file for the well-known Ringleb flow
    /// </summary>
    public class RinglebControl : CNSControl {

        /// <summary>
        /// The equation of state to be used; use
        /// <see cref="StiffenedGas.ReferencePressure"/> = 0 for ideal gas
        /// </summary>
        private StiffenedGas equationOfState;

        /// <summary>
        /// Type-safe access to <see cref="equationOfState"/>
        /// </summary>
        public new StiffenedGas EquationOfState {
            get {
                return equationOfState;
            }
            set {
                this.equationOfState = value;
                base.EquationOfState = value;
            }
        }

        /// <summary>
        /// The reference speed of sound \f$ a_0\f$  that
        /// has to be specified in the definition of the Ringleb problem
        /// </summary>
        public double RinglebReferenceSpeedOfSound = 1.0;

        private double? ringlebReferenceTotalPressure = null;

        /// <summary>
        /// The reference total pressure \f$ p^t_0\f$ 
        /// that has to be specified in the definition of the Ringleb problem
        /// </summary>
        public double RinglebReferenceTotalPressure {
            get {
                if (!ringlebReferenceTotalPressure.HasValue) {
                    ringlebReferenceTotalPressure =
                        RinglebReferenceSpeedOfSound * RinglebReferenceSpeedOfSound / EquationOfState.HeatCapacityRatio
                        - equationOfState.ReferencePressure;
                }
                return ringlebReferenceTotalPressure.Value;
            }
            set {
                ringlebReferenceTotalPressure = value;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        override public Material GetMaterial() {
            return new RinglebMaterial(EquationOfState, ViscosityLaw, MachNumber, ReynoldsNumber, PrandtlNumber, FroudeNumber, ViscosityRatio) {
                RinglebReferenceSpeedOfSound = this.RinglebReferenceSpeedOfSound,
                RinglebReferencePressure = this.equationOfState.ReferencePressure,
                RinglebReferenceTotalPressure = this.RinglebReferenceTotalPressure
            };
        }
    }
}
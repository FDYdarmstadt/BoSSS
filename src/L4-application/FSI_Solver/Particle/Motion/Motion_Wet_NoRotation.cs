/* =======================================================================
Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using System.Collections.Generic;

namespace BoSSS.Application.FSI_Solver {
    public class Motion_Wet_NoRotation : Motion_Wet {

        /// <summary>
        /// The dry description of motion including hydrodynamics without rotation.
        /// </summary>
        /// <param name="gravity">
        /// The gravity (volume forces) acting on the particle.
        /// </param>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        /// /// <param name="underrelaxationParam">
        /// The underrelaxation parameters (convergence limit, prefactor and a bool whether to use addaptive underrelaxation) defined in <see cref="ParticleUnderrelaxationParam"/>.
        /// </param>
        public Motion_Wet_NoRotation(double[] gravity, double density) : base(new Vector(gravity), density) {
            IncludeRotation = false;
        }

        /// <summary>
        /// Include translation?
        /// </summary>
        internal override bool IncludeRotation { get; } = false;

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected override double CalculateParticleAngle(double dt) {
            double l_Angle = GetAngle(1);
            Aux.TestArithmeticException(l_Angle, "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new particle angle after a collision
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected override double CalculateParticleAngle(double dt, double collisionTimestep) {
            double l_Angle = GetAngle(1);
            Aux.TestArithmeticException(l_Angle, "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected override double CalculateAngularVelocity(double dt) {
            double l_RotationalVelocity = 0;
            Aux.TestArithmeticException(l_RotationalVelocity, "particle rotational velocity");
            return l_RotationalVelocity;
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected override double CalculateAngularVelocity(double dt, double collisionTimestep) {
            double l_RotationalVelocity = 0;
            Aux.TestArithmeticException(l_RotationalVelocity, "particle rotational velocity");
            return l_RotationalVelocity;
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        public override Vector CalculateHydrodynamicForces(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, double fluidDensity, CellMask cutCells, double dt) {
            Vector tempForces = new Vector(hydrodynamicsIntegration.Forces(out List<double[]>[] stressToPrintOut, cutCells));
            currentStress = TransformStressToPrint(stressToPrintOut);
            Aux.TestArithmeticException(tempForces, "temporal forces during calculation of hydrodynamics");
            tempForces = Force_MPISum(tempForces);
            tempForces = CalculateGravity(fluidDensity, tempForces);
            return tempForces;
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="levelSetTracker"></param>
        /// <param name="fluidViscosity"></param>
        /// <param name="cutCells"></param>
        /// <param name="dt"></param>
        public override double CalculateHydrodynamicTorque(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, double dt) {
            return 0;
        }
    }
}

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

using BoSSS.Foundation.Grid;
using ilPSP;
using System;

namespace BoSSS.Application.XNSERO_Solver {
    [Serializable]
    public class MotionFixedWithVelocity : Motion {

        /// <summary>
        /// No motion
        /// </summary>
        /// <param name="gravity">
        /// The gravity (volume forces) acting on the particle.
        /// </param>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        public MotionFixedWithVelocity(double density = 0) : base(density) {  }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        public override Vector CalculateParticlePosition(double dt) {
            Vector position = GetPosition(1);
            Aux.TestArithmeticException(position, "particle translational velocity");
            return position;
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        public override double CalculateParticleAngle(double dt) {
            double angle = GetAngle(1);
            Aux.TestArithmeticException(angle, "particle rotational velocity");
            return angle;
        }

        public override object Clone() {
            MotionFixed clonedMotion = new MotionFixed(Density);
            clonedMotion.Volume = this.Volume;
            clonedMotion.MomentOfInertia = this.MomentOfInertia;
            return clonedMotion;
        }
    }
}

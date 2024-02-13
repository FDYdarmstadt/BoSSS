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

using ilPSP;
using System;

namespace BoSSS.Application.XNSERO_Solver {
    [Serializable]
    public class MotionOnlyVertical : Motion {

        /// <summary>
        /// The dry description of motion without hydrodynamics and rotation.
        /// </summary>
        /// <param name="gravity">
        /// The gravity (volume forces) acting on the particle.
        /// </param>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        public MotionOnlyVertical(double density) : base(density) { }

        public override bool IncludeRotation() {
            return false;
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public override double CalculateParticleAngle(double dt) {
            double l_Angle = GetAngle(1);
            Aux.TestArithmeticException(l_Angle, "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override Vector CalculateParticlePosition(double dt) {
            Vector position = Position[1] + (5 * TranslationalVelocity[0] + 8 * TranslationalVelocity[1] - TranslationalVelocity[2]) * dt / 12;
            position[0] = Position[1][0];
            position[1] = Position[1][1];
            Aux.TestArithmeticException(position, "particle position");
            return position;
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override Vector CalculateParticlePositionDuringCollision(double dt) {
            Vector position = Position[0] + (5 * TranslationalVelocity[0] + 8 * TranslationalVelocity[1] - TranslationalVelocity[2]) * dt / 12;
            position[0] = Position[0][0];
            Aux.TestArithmeticException(position, "particle position");
            return position;
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public override double CalculateRotationalAcceleration(double dt) {
            return 0;
        }

        public override object Clone() {
            MotionDryNoRotation clonedMotion = new MotionDryNoRotation(Density);
            clonedMotion.Volume = this.Volume;
            clonedMotion.MomentOfInertia = this.MomentOfInertia;
            return clonedMotion;
        }
    }
}

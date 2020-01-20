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
using ilPSP;

namespace BoSSS.Application.FSI_Solver {
    public class MotionFixedWithForces : Motion {

        /// <summary>
        /// No motion
        /// </summary>
        /// <param name="gravity">
        /// The gravity (volume forces) acting on the particle.
        /// </param>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        public MotionFixedWithForces(Vector gravity, double density = 0) : base(new Vector(gravity), density) {
            IncludeRotation = false;
            IncludeTranslation = false;
        }

        /// <summary>
        /// Include rotation?
        /// </summary>
        internal override bool IncludeRotation { get; } = false;

        /// <summary>
        /// Include translation?
        /// </summary>
        internal override bool IncludeTranslation { get; } = false;

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected override Vector CalculateParticlePosition(double dt, double collisionTimestep) {
            Vector l_Position = GetPosition(1);
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected override double CalculateParticleAngle(double dt, double collisionTimestep = 0) {
            double l_Angle = GetAngle(1);
            Aux.TestArithmeticException(l_Angle, "particle angle");
            return l_Angle;
        }

        public override object Clone() {
            Motion clonedMotion = new MotionFixedWithForces(Gravity, Density);
            clonedMotion.GetParticleArea(ParticleArea);
            clonedMotion.GetParticleMomentOfInertia(MomentOfInertia);
            return clonedMotion;
        }
    }
}

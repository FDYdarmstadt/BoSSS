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

namespace BoSSS.Application.FSI_Solver {
    [Serializable]
    public class MotionWetNoTranslation : Motion {

        /// <summary>
        /// The dry description of motion including hydrodynamics without translation.
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
        public MotionWetNoTranslation(Vector gravity, double density) : base(new Vector(gravity), density) {
            IncludeTranslation = false;
        }

        /// <summary>
        /// Include translation?
        /// </summary>
        internal override bool IncludeTranslation { get; } = false;

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected override Vector CalculateParticlePosition(double dt) {
            Vector l_Position = GetPosition(1);
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override Vector CalculateTranslationalVelocity(double dt) {
            Vector l_TranslationalVelocity = new Vector(SpatialDim);
            Aux.TestArithmeticException(l_TranslationalVelocity, "particle translational velocity");
            return l_TranslationalVelocity;
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        public override Vector CalculateHydrodynamicForces(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, double fluidDensity, CellMask cutCells, double dt) {
            return new Vector(SpatialDim); ;
        }

        public override object Clone() {
            Motion clonedMotion = new MotionWetNoTranslation(Gravity, Density);
            clonedMotion.SetParticleArea(ParticleArea);
            clonedMotion.SetParticleMomentOfInertia(MomentOfInertia);
            return clonedMotion;
        }
    }
}

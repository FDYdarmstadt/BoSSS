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
        public MotionWetNoTranslation(double density) : base(density) {  }

        public override bool IncludeTranslation() {
            return false;
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override Vector CalculateParticlePosition(double dt) {
            Vector l_Position = GetPosition(1);
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        public override Vector CalculateTranslationalVelocity(double dt) {
            Vector l_TranslationalVelocity = new Vector(SpatialDim);
            Aux.TestArithmeticException(l_TranslationalVelocity, "particle translational velocity");
            return l_TranslationalVelocity;
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        //public override Vector CalculateHydrodynamics(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, string[] FluidSpecies, double dt) {
        //    int NoOfFluidSpecies = FluidSpecies.Length;
        //    Vector forcesAndTorque = new Vector(SpatialDim + 1);
        //    for(int i = 0; i < NoOfFluidSpecies; i++) {
        //        forcesAndTorque += hydrodynamicsIntegration.Main(GetPosition(0), CutCells, FluidSpecies[i]);
        //    }
        //    forcesAndTorque[0] = 0;
        //    forcesAndTorque[1] = 0;
        //    Aux.TestArithmeticException(forcesAndTorque, "during calculation of hydrodynamics");
        //    return forcesAndTorque;
        //}

        public override object Clone() {
            MotionWetNoTranslation clonedMotion = new MotionWetNoTranslation(Density);
            clonedMotion.Volume = this.Volume;
            clonedMotion.MomentOfInertia = this.MomentOfInertia;
            return clonedMotion;
        }

        /// <summary>
        /// Calculates the gravitational forces.
        /// </summary>
        /// <param name="fluidDensity"></param>
        /// <param name="tempForces"></param>
        public override Vector GravityForce(Vector Gravity) {
            return new Vector(0, 0);
        }
    }
}

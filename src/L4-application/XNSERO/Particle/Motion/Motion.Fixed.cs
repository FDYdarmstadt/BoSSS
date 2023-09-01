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
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSERO_Solver {
    [Serializable]
    public class MotionFixed : Motion {

        /// <summary>
        /// No motion
        /// </summary>
        /// <param name="gravity">
        /// The gravity (volume forces) acting on the particle.
        /// </param>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        public MotionFixed(double density = 0) : base(density) { }

        public override bool IncludeRotation() {
            return false;
        }
        public override bool IncludeTranslation() {
            return false;
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
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        public override double CalculateAngularVelocity(double dt) {
            double l_RotationalVelocity = 0;
            Aux.TestArithmeticException(l_RotationalVelocity, "particle rotational velocity");
            return l_RotationalVelocity;
        }

        /// <summary>
        /// Calculates the new translational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        public override Vector CalculateTranslationalAcceleration(double dt = 0) {
            return new Vector(SpatialDim);
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public override double CalculateRotationalAcceleration(double dt) {
            return 0;
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        //public override Vector CalculateHydrodynamics(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, string[] FluidSpecies, double dt) {
        //    return new Vector(SpatialDim + 1);
        //}

        public override object Clone() {
            MotionFixed clonedMotion = new MotionFixed(Density);
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

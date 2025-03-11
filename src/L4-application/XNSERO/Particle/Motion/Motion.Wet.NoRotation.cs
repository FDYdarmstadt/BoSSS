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
using System.Collections.Generic;
using System;

namespace BoSSS.Application.XNSERO_Solver {
    [Serializable]
    public class MotionWetNoRotation : Motion {

        /// <summary>
        /// The dry description of motion including hydrodynamics without rotation.
        /// </summary>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        public MotionWetNoRotation(double density) : base(density) { }

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
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        //public override Vector CalculateHydrodynamics(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, string[] FluidSpecies, double dt) {
        //    int NoOfFluidSpecies = FluidSpecies.Length;
        //    Vector tempForces = new Vector(SpatialDim + 1);
        //    for(int i = 0; i < NoOfFluidSpecies; i++) {
        //        tempForces += hydrodynamicsIntegration.Main(GetPosition(0), CutCells, FluidSpecies[i]);
        //    }
        //    Aux.TestArithmeticException(tempForces, "temporal forces during calculation of hydrodynamics");
        //    tempForces[SpatialDim] = 0;
        //    return tempForces;
        //}

        public override object Clone() {
            MotionWetNoRotation clonedMotion = new MotionWetNoRotation(Density);
            clonedMotion.Volume = this.Volume;
            clonedMotion.MomentOfInertia = this.MomentOfInertia;
            return clonedMotion;
        }
    }
}

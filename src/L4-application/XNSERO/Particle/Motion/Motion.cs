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
using ilPSP.Tracing;
using System;

namespace BoSSS.Application.XNSERO_Solver {
    [Serializable]
    public class Motion : MotionBase {

        /// <summary>
        /// The standard description of motion including hydrodynamics.
        /// </summary>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        public Motion(double density) : base(density) { }

        public override bool IncludeRotation() {
            return true;
        }

        public override bool IncludeTranslation() {
            return true;
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override Vector CalculateParticlePosition(double dt) {
            using (new FuncTrace()) {
                Vector position = Position[1] + (5 * TranslationalVelocity[0] + 8 * TranslationalVelocity[1] - TranslationalVelocity[2]) * dt / 12;
                Aux.TestArithmeticException(position, "particle position");
                return position;
            }
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public override double CalculateParticleAngle(double dt) {
            using (new FuncTrace()) {
                double angle = Angle[1] + (5 * RotationalVelocity[0] + 8 * RotationalVelocity[1] - RotationalVelocity[2]) * dt / 12;
                Aux.TestArithmeticException(angle, "particle angle");
                return angle;
            }
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override Vector CalculateParticlePositionDuringCollision(double dt) {
            using (new FuncTrace()) {
                Vector position = Position[0] + (5 * TranslationalVelocity[0] + 8 * TranslationalVelocity[1] - TranslationalVelocity[2]) * dt / 12;
                Aux.TestArithmeticException(position, "particle position");
                return position;
            }
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public override double CalculateParticleAngleDuringCollision(double dt) {
            using (new FuncTrace()) {
                double angle = Angle[0] + (5 * RotationalVelocity[0] + 8 * RotationalVelocity[1] - RotationalVelocity[2]) * dt / 12;
                Aux.TestArithmeticException(angle, "particle angle");
                return angle;
            }
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Time-step</param>
        public override Vector CalculateTranslationalVelocity(double dt) {
            using (new FuncTrace()) {
                Vector translationalVelocity = TranslationalVelocity[1] + (5 * TranslationalAcceleration[0] + 8 * TranslationalAcceleration[1] - TranslationalAcceleration[2]) * dt / 12;
                Aux.TestArithmeticException(translationalVelocity, "particle translational velocity");
                return translationalVelocity;
            }
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Time-step</param>
        public override double CalculateAngularVelocity(double dt) {
            using (new FuncTrace()) {
                double rotationalVelocity = RotationalVelocity[1] + (5 * RotationalAcceleration[0] + 8 * RotationalAcceleration[1] - RotationalAcceleration[2]) * dt / 12;
                Aux.TestArithmeticException(rotationalVelocity, "particle rotational velocity");
                return rotationalVelocity;
            }
        }

        /// <summary>
        /// Calculate the new translational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        public override Vector CalculateTranslationalAcceleration(double dt) {
            using (new FuncTrace()) {
                Vector l_Acceleration = HydrodynamicForces[0] / (Density * Volume);
                Aux.TestArithmeticException(l_Acceleration, "particle translational acceleration");
                return l_Acceleration;
            }
        }

        /// <summary>
        /// Calculate the new rotational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        public override double CalculateRotationalAcceleration(double dt) {
            using (new FuncTrace()) {
                double l_Acceleration = HydrodynamicTorque[0] / MomentOfInertia;
                Aux.TestArithmeticException(l_Acceleration, "particle rotational acceleration");
                return l_Acceleration;
            }
        }

        /// <summary>
        /// Calls the integration of the hydrodynamic stress at this particles level-set
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        public override Vector CalculateHydrodynamics(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, string[] FluidSpecies, double dt = 0) {
            Vector returnVec;

            int NoOfFluidSpecies = FluidSpecies.Length;
            Vector forcesAndTorque = new Vector(SpatialDim + 1);
            for (int i = 0; i < NoOfFluidSpecies; i++) {
                forcesAndTorque += hydrodynamicsIntegration.Main(Position[0], cutCells, FluidSpecies[i]);
            }
            Aux.TestArithmeticException(forcesAndTorque, "during calculation of hydrodynamics");

            if (cutCells.IsEmptyOnRank)
                returnVec = new Vector(SpatialDim + 1);
            else {
                returnVec = forcesAndTorque;
            }
            return returnVec;
        }

        public override object Clone() {
            Motion clonedMotion = new Motion(Density);
            clonedMotion.Volume = this.Volume;
            clonedMotion.MomentOfInertia = this.MomentOfInertia;
            return clonedMotion;
        }
    }
}

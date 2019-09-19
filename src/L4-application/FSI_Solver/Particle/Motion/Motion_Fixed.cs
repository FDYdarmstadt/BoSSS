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

namespace BoSSS.Application.FSI_Solver {
    public class Motion_Fixed : Motion_Wet {

        public Motion_Fixed(double[] gravity = null, double density = 0) : base(gravity, density) {
            IncludeRotation = false;
            IncludeTranslation = false;
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateParticlePosition(double dt, double collisionTimestep = 0) {
            position[0] = position[1];
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateParticleAngle(double dt, double collisionTimestep = 0) {
            angle[0] = angle[1];
            Aux.TestArithmeticException(angle[0], "particle angle");
        }


        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        protected override void CalculateTranslationalVelocity(double dt) {
            CalculateTranslationalAcceleration();
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            CalculateTranslationalAcceleration();
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        protected override void CalculateAngularVelocity(double dt = 0) {
            rotationalVelocity[0] = 0;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override void CalculateAngularVelocity(double dt = 0, double collisionTimestep = 0) {
            rotationalVelocity[0] = 0;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculates the new translational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateTranslationalAcceleration(double dt = 0) {
            for (int d = 0; d < spatialDim; d++) {
                translationalAcceleration[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle translational acceleration");
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateRotationalAcceleration(double dt) {
            rotationalAcceleration[0] = 0;
            Aux.TestArithmeticException(rotationalAcceleration[0], "particle rotational acceleration");
        }

        /// <summary>
        /// Overrides the calculation of hydrodynamics for fixed particles, so that nothing happens.
        /// </summary>
        public override void UpdateForcesAndTorque(VectorField<SinglePhaseField> U = null, SinglePhaseField P = null, LevelSetTracker LsTrk = null, CellMask CutCells_P = null, double fluidViscosity = 0, double fluidDensity = 0, bool firstIteration = false, double dt = 0) {
            for (int d = 0; d < spatialDim; d++) {
                HydrodynamicForces[0][d] = 0;
            }
            HydrodynamicTorque[0] = 0;
        }
    }
}

﻿/* =======================================================================
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
    public class Motion_Dry : Motion_Wet {
        public Motion_Dry(double[] gravity) : base(gravity) {
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        protected override void CalculateTranslationalVelocity(double dt) {
            CalculateTranslationalAcceleration(dt);
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + translationalAcceleration[0][d] * dt;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            CalculateTranslationalAcceleration(dt);
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + translationalAcceleration[0][d] * (dt - collisionTimestep);
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        protected override void CalculateAngularVelocity(double dt = 0) {
            rotationalVelocity[0] = rotationalVelocity[1];
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override void CalculateAngularVelocity(double dt = 0, double collisionTimestep = 0) {
            rotationalVelocity[0] = rotationalVelocity[1];
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
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
                HydrodynamicForces[0][d] = m_Gravity[d] * particleDensity * particleArea;
            }
            HydrodynamicTorque[0] = 0;
        }
    }
}

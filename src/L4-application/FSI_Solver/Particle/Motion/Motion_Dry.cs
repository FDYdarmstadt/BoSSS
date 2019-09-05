using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using FSI_Solver;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.FSI_Solver {
    public class Motion_Dry : ParticleMotion {
        public Motion_Dry(double[] gravity) : base(gravity) {
        }

        public override void CheckCorrectInit(bool IsDry) {
            // Do nothing
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        protected override void CalculateTranslationalVelocity(double dt) {
            CalculateTranslationalAcceleration();
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
            CalculateTranslationalAcceleration();
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
        /// Calculates the new translational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateTranslationalAcceleration(double dt = 0) {
            for (int d = 0; d < spatialDim; d++) {
                translationalAcceleration[0][d] = m_Gravity[d];
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle translational acceleration");
        }
    }
}

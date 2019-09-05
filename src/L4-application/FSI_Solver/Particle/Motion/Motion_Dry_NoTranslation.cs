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
    public class Motion_Dry_NoTranslation : Motion_Dry {
        public Motion_Dry_NoTranslation(double[] gravity) : base(gravity) {
        }

        public override void CheckCorrectInit(bool IsDry) {
            // Do nothing
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateParticlePosition(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d];
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateParticlePosition(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d];
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        protected override void CalculateTranslationalVelocity(double dt) {
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
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculates the new translational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateTranslationalAcceleration(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                translationalAcceleration[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle translational acceleration");
        }
    }
}

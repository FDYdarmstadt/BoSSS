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

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticlePosition(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d] + translationalVelocity[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticleAngle(double dt, double collisionTimestep) {
            angle[0] = angle[1] + rotationalVelocity[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public override void CalculateTranslationalVelocity(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + m_Gravity[d] * dt;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        public override void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + m_Gravity[d] * (dt - collisionTimestep);
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public override void CalculateAngularVelocity(double dt = 0) {
            rotationalVelocity[0] = rotationalVelocity[1];
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        public override void CalculateAngularVelocity(double dt = 0, double collisionTimestep = 0) {
            rotationalVelocity[0] = rotationalVelocity[1];
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }
    }
}

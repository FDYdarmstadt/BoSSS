using System.Collections.Generic;

namespace BoSSS.Application.FSI_Solver {
    public class Motion_Dry : ParticleMotion {
        public Motion_Dry(double[] gravity) : base() {
            m_Gravity = gravity;
        }

        private double[] m_Gravity;

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public override void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            int spatialDim = translationalVelocity[0].Length;
            double[] newTranslationalVelocity = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + m_Gravity[d] * dt;
            }
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public override void CalculateAngularVelocity(double dt, double collisionTimestep) {
            rotationalVelocity[0] = rotationalVelocity[1];
        }
    }
}

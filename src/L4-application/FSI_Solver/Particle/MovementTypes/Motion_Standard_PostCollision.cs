using BoSSS.Foundation.XDG;
using FSI_Solver;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.FSI_Solver {
    public class Motion_Standard_PostCollision : ParticleMotion {
        public Motion_Standard_PostCollision() : base() {
        }
        readonly private FSI_Auxillary Aux = new FSI_Auxillary();

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticlePosition(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d] + translationalVelocity[0][d] * (dt - collisionTimestep) / 6;
            }
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticleAngle(double dt, double collisionTimestep) {
            int clearPreCollisionValues = collisionTimestep != 0 ? 0 : 1;
            angle[0] = angle[1] + rotationalVelocity[0] * (dt - collisionTimestep) / 6;
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateAcceleration(double dt, double[] force, double torque) {
            for (int d = 0; d < spatialDim; d++) {
                translationalAcceleration[0][d] = force[d] / particleMass;
            }

            rotationalAcceleration[0] = torque / particleMomentOfInertia;
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public override void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            int clearPreCollisionValues = collisionTimestep != 0 ? 0 : 1;
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + translationalAcceleration[0][d] * (dt - collisionTimestep) / 6;
            }
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public override void CalculateAngularVelocity(double dt, double collisionTimestep) {
            int clearPreCollisionValues = collisionTimestep != 0 ? 0 : 1;
            rotationalVelocity[0] = rotationalVelocity[1] + rotationalAcceleration[0] * (dt - collisionTimestep) / 6;
        }
    }
}
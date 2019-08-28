using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver {
    public class ParticleMotion {

        private const int historyLength = 4;
        protected static int spatialDim = 2;

        public ParticleMotion() {
            for (int i = 0; i < historyLength; i++) {
                position.Add(new double[spatialDim]);
                angle.Add(new double());
                translationalVelocity.Add(new double[spatialDim]);
                translationalAcceleration.Add(new double[spatialDim]);
                rotationalVelocity.Add(new double());
                rotationalAcceleration.Add(new double());
            }
        }

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public double particleMass = new double();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public double particleMomentOfInertia = new double();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> position = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> angle = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> translationalVelocity = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> rotationalVelocity = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> translationalAcceleration = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> rotationalAcceleration = new List<double>();

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public virtual void CalculateParticlePosition(double dt, double collisionTimestep) {
            int clearPreCollisionValues = collisionTimestep != 0 ? 0 : 1;
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d] + (translationalVelocity[0][d] + clearPreCollisionValues * (4 * translationalVelocity[1][d] + translationalVelocity[2][d])) * (dt - collisionTimestep) / 6;
            }
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public virtual void CalculateParticleAngle(double dt, double collisionTimestep) {
            int clearPreCollisionValues = collisionTimestep != 0 ? 0 : 1;
            angle[0] = angle[1] + (rotationalVelocity[0] + clearPreCollisionValues * (4 * rotationalVelocity[1] + rotationalVelocity[2])) * (dt - collisionTimestep) / 6;
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public virtual void CalculateAcceleration(double dt, double[] force, double torque) {
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
        public virtual void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            int clearPreCollisionValues = collisionTimestep != 0 ? 0 : 1;
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + (translationalAcceleration[0][d] + clearPreCollisionValues * (4 * translationalAcceleration[1][d] + translationalAcceleration[2][d])) * dt / 6;
            }
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public virtual void CalculateAngularVelocity(double dt, double collisionTimestep) {
            int clearPreCollisionValues = collisionTimestep != 0 ? 0 : 1;
            rotationalVelocity[0] = rotationalVelocity[1] + dt * (rotationalAcceleration[0] + clearPreCollisionValues * (4 * rotationalAcceleration[1] + rotationalAcceleration[2])) / 6;
        }
    }
}

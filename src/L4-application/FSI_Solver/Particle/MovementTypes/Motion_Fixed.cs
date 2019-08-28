using System.Collections.Generic;

namespace BoSSS.Application.FSI_Solver {
    public class Motion_Fixed : ParticleMotion {
        public Motion_Fixed() : base() { }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticlePosition(double dt, double collisionTimestep = 0) {
            position[0] = position[1];
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticleAngle(double dt, double collisionTimestep = 0) {
            angle[0] = angle[1];
        }
    }
}

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
    public class Motion_Dry_NoRotation : Motion_Dry {
        public Motion_Dry_NoRotation(double[] gravity) : base(gravity) {
        }

        public override void CheckCorrectInit(bool IsDry) {
            // Do nothing
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateParticleAngle(double dt) {
            angle[0] = angle[1];
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected override void CalculateParticleAngle(double dt, double collisionTimestep) {
            angle[0] = angle[1];
            Aux.TestArithmeticException(angle[0], "particle angle");
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
    }
}

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
    public class Motion_Wet_NoRotation : ParticleMotion {
        public Motion_Wet_NoRotation(double[] gravity) : base(gravity) {
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
        /// Calculate the new particle angle after a collision
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected override void CalculateParticleAngle(double dt, double collisionTimestep) {
            angle[0] = angle[1];
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override void CalculateAngularVelocity(double dt) {
            rotationalVelocity[0] = 0;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected override void CalculateAngularVelocity(double dt, double collisionTimestep) {
            rotationalVelocity[0] = 0;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calls the calculation of the hydrodynamics
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        public override void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double relativeParticleMass, double dt = 0) {
            double[] tempForces = CalculateHydrodynamicForces(U, P, LsTrk, CutCells_P, muA, relativeParticleMass);
            HydrodynamicsPostprocessing(tempForces);
        }

        protected override void HydrodynamicsPostprocessing(double[] tempForces, double tempTorque = 0) {
            for (int d = 0; d < spatialDim; d++) {
                hydrodynamicForces[0][d] = tempForces[d];
            }
            hydrodynamicTorque[0] = 0;
            Aux.TestArithmeticException(hydrodynamicForces[0], "hydrodynamic forces");
        }
    }
}

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
    public class Motion_Wet_NoTranslation : ParticleMotion {
        public Motion_Wet_NoTranslation(double[] gravity) : base(gravity) {
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
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected override void CalculateParticlePosition(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d];
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override void CalculateTranslationalVelocity(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected override void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calls the calculation of the hydrodynamics
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        public override void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double relativeParticleMass, double dt = 0) {
            double[] tempForces = new double[spatialDim];
            double tempTorque = CalculateHydrodynamicTorque(U, P, LsTrk, CutCells_P, muA);
            HydrodynamicsPostprocessing(tempForces, tempTorque);
        }

        protected override void HydrodynamicsPostprocessing(double[] tempForces, double tempTorque) {
            for (int d = 0; d < spatialDim; d++) {
                hydrodynamicForces[0][d] = 0;
            }
            hydrodynamicTorque[0] = tempTorque;
            Aux.TestArithmeticException(hydrodynamicTorque[0], "hydrodynamic torque");
        }
    }
}

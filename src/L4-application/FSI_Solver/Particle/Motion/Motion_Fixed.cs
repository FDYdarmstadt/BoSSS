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
    public class Motion_Fixed : ParticleMotion {
        public Motion_Fixed(double[] gravity) : base(gravity) {
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticlePosition(double dt, double collisionTimestep = 0) {
            position[0] = position[1];
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticleAngle(double dt, double collisionTimestep = 0) {
            angle[0] = angle[1];
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Overrides the calculation of hydrodynamics for fixed particles, so that nothing happens.
        /// </summary>
        public override void UpdateForcesAndTorque(VectorField<SinglePhaseField> U = null, SinglePhaseField P = null, LevelSetTracker LsTrk = null, CellMask CutCells_P = null, double muA = 0, double relativeParticleMass = 0, double dt = 0) {
            for (int d = 0; d < spatialDim; d++) {
                hydrodynamicForces[0][d] = 0;
            }
            hydrodynamicTorque[0] = 0;
        }
    }
}

using BoSSS.Foundation.XDG;
using FSI_Solver;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.FSI_Solver {
    public class Motion_AddedDamping : ParticleMotion {
        public Motion_AddedDamping() : base() {
            m_StartingAngle = angle[0];
        }
        readonly private FSI_Auxillary Aux = new FSI_Auxillary();
        readonly private ParticleAddedDamping AddedDamping = new ParticleAddedDamping();
        readonly private ParticleAcceleration Acceleration = new ParticleAcceleration();

        /// <summary>
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        public double[,] addedDampingTensor = new double[6, 6];

        /// <summary>
        /// Added damping coefficient, should be between 0.5 and 1.5, for reference: Banks et.al. 2017
        /// </summary>
        readonly int addedDampingCoefficient = 1;

        /// <summary>
        /// Saving the initial angle of the particle for <see cref="UpdateDampingTensors()"/>
        /// </summary>
        double m_StartingAngle;

        /// <summary>
        /// Calculate tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public void CalculateDampingTensor(Particle particle, LevelSetTracker LsTrk, double muA, double rhoA, double dt) {
            addedDampingTensor = AddedDamping.IntegrationOverLevelSet(particle, LsTrk, muA, rhoA, dt, position[0]);
            Aux.TestArithmeticException(addedDampingTensor, "particle added damping tensor");
        }

        /// <summary>
        /// Update in every timestep tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public void UpdateDampingTensors() {
            addedDampingTensor = AddedDamping.RotateTensor(angle[0], m_StartingAngle, addedDampingTensor);
            Aux.TestArithmeticException(addedDampingTensor, "particle added damping tensor");
        }

        /// <summary>
        /// Predict the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public void PredictAcceleration() {

            for (int d = 0; d < spatialDim; d++) {
                translationalAcceleration[0][d] = (translationalAcceleration[1][d] + 4 * translationalAcceleration[2][d] + translationalAcceleration[3][d]) / 8;
                if (Math.Abs(translationalAcceleration[0][d]) < 1e-20)
                    translationalAcceleration[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle acceleration");

            rotationalAcceleration[0] = (rotationalAcceleration[1] + 4 * rotationalAcceleration[2] + rotationalAcceleration[3]) / 8;
            if (Math.Abs(rotationalAcceleration[0]) < 1e-20)
                rotationalAcceleration[0] = 0;
            Aux.TestArithmeticException(rotationalAcceleration[0], "particle angular acceleration");
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateAcceleration(double dt, double[] force, double torque) {
            double[,] coefficientMatrix = Acceleration.CalculateCoefficientMatrix(addedDampingTensor, particleMass, particleMomentOfInertia, dt, addedDampingCoefficient);
            Aux.TestArithmeticException(coefficientMatrix, "particle acceleration coefficients");
            double Denominator = Acceleration.CalculateDenominator(coefficientMatrix);
            Aux.TestArithmeticException(Denominator, "particle acceleration denominator");

            translationalAcceleration[0] = Acceleration.Translational(coefficientMatrix, Denominator, force, torque);

            for (int d = 0; d < spatialDim; d++) {
                if (Math.Abs(translationalAcceleration[0][d]) < 1e-20)
                    translationalAcceleration[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle acceleration");

            rotationalAcceleration[0] = Acceleration.Rotational(coefficientMatrix, Denominator, force, torque);
            if (Math.Abs(rotationalAcceleration[0]) < 1e-20)
                rotationalAcceleration[0] = 0;
        }
    }
}

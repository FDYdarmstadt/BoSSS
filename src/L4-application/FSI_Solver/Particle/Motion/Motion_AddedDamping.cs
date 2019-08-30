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
    public class Motion_AddedDamping : ParticleMotion {
        public Motion_AddedDamping(
            double[] gravity,
            double maxParticleLengthScale,
            double forceAndTorqueConvergence,
            double underrelaxationFactor,
            bool useAddaptiveUnderrelaxation = true,
            double addedDampingCoefficient = 1) : base(gravity) {
            m_StartingAngle = angle[0];
            m_AddedDampingCoefficient = addedDampingCoefficient;
            m_ForceAndTorqueConvergence = forceAndTorqueConvergence;
            m_UnderrelaxationFactor = underrelaxationFactor;
            m_UseAddaptiveUnderrelaxation = useAddaptiveUnderrelaxation;
            m_MaxParticleLengthScale = maxParticleLengthScale;
        }

        readonly ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation();
        readonly ParticleAddedDamping AddedDamping = new ParticleAddedDamping();
        readonly ParticleAcceleration Acceleration = new ParticleAcceleration();

        /// <summary>
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        public double[,] addedDampingTensor = new double[6, 6];

        /// <summary>
        /// Added damping coefficient, should be between 0.5 and 1.5, for reference: Banks et.al. 2017
        /// </summary>
        protected double m_AddedDampingCoefficient;

        /// <summary>
        /// Saving the initial angle of the particle for <see cref="UpdateDampingTensors()"/>
        /// </summary>
       protected double m_StartingAngle;

        /// <summary>
        /// Convergence criterion.
        /// </summary>
        protected double m_ForceAndTorqueConvergence;

        /// <summary>
        /// Force and torque underrelaxation.
        /// </summary>
        protected double m_UnderrelaxationFactor;

        /// <summary>
        /// 
        /// </summary>
        protected bool m_UseAddaptiveUnderrelaxation;

        /// <summary>
        /// 
        /// </summary>
        protected double m_MaxParticleLengthScale;

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
        /// Calculate the new particle position after a collision
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticlePosition(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d] + translationalVelocity[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle angle after a collision
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateParticleAngle(double dt, double collisionTimestep) {
            angle[0] = angle[1] + rotationalVelocity[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public override void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + translationalAcceleration[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <returns></returns>
        public override void CalculateAngularVelocity(double dt, double collisionTimestep) {
            rotationalVelocity[0] = rotationalVelocity[1] + rotationalAcceleration[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public override void CalculateAcceleration(double dt, double[] force, double torque) {
            // Added damping
            double[,] coefficientMatrix = Acceleration.CalculateCoefficientMatrix(addedDampingTensor, particleMass, particleMomentOfInertia, dt, m_AddedDampingCoefficient);
            Aux.TestArithmeticException(coefficientMatrix, "particle acceleration coefficients");
            double Denominator = Acceleration.CalculateDenominator(coefficientMatrix);
            Aux.TestArithmeticException(Denominator, "particle acceleration denominator");

            // Translation
            translationalAcceleration[0] = Acceleration.Translational(coefficientMatrix, Denominator, force, torque);
            for (int d = 0; d < spatialDim; d++) {
                if (Math.Abs(translationalAcceleration[0][d]) < 1e-20)
                    translationalAcceleration[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle translational acceleration");

            // Rotation
            rotationalAcceleration[0] = Acceleration.Rotational(coefficientMatrix, Denominator, force, torque);
            Aux.TestArithmeticException(translationalAcceleration[0], "particle rotational acceleration");
            if (Math.Abs(rotationalAcceleration[0]) < 1e-20)
                rotationalAcceleration[0] = 0;
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        public void CalculateAcceleration(double dt, double collisionTimestep, double[] force, double torque) {
            // Added damping
            double[,] coefficientMatrix = Acceleration.CalculateCoefficientMatrix(addedDampingTensor, particleMass, particleMomentOfInertia, (dt - collisionTimestep), m_AddedDampingCoefficient);
            Aux.TestArithmeticException(coefficientMatrix, "particle acceleration coefficients");
            double Denominator = Acceleration.CalculateDenominator(coefficientMatrix);
            Aux.TestArithmeticException(Denominator, "particle acceleration denominator");

            // Translation
            translationalAcceleration[0] = Acceleration.Translational(coefficientMatrix, Denominator, force, torque);
            for (int d = 0; d < spatialDim; d++) {
                if (Math.Abs(translationalAcceleration[0][d]) < 1e-20)
                    translationalAcceleration[0][d] = 0;
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle translational acceleration");

            // Rotation
            rotationalAcceleration[0] = Acceleration.Rotational(coefficientMatrix, Denominator, force, torque);
            Aux.TestArithmeticException(translationalAcceleration[0], "particle rotational acceleration");
            if (Math.Abs(rotationalAcceleration[0]) < 1e-20)
                rotationalAcceleration[0] = 0;
        }

        /// <summary>
        /// Predicts the new acceleration (translational and rotational)
        /// </summary>
        public void PredictForceAndTorque() {
            for (int d = 0; d < spatialDim; d++) {
                hydrodynamicForces[0][d] = (hydrodynamicForces[1][d] + 4 * hydrodynamicForces[2][d] + hydrodynamicForces[3][d]) / 6;
                if (Math.Abs(hydrodynamicForces[0][d]) < 1e-20)
                    hydrodynamicForces[0][d] = 0;
            }
            Aux.TestArithmeticException(hydrodynamicForces[0], "hydrodynamic forces");
            hydrodynamicTorque[0] = (hydrodynamicTorque[1] + 4 * hydrodynamicTorque[2] + hydrodynamicTorque[3]) / 6;
            if (Math.Abs(hydrodynamicTorque[0]) < 1e-20)
                hydrodynamicTorque[0] = 0;
            Aux.TestArithmeticException(hydrodynamicTorque[0], "hydrodynamic torque");
        }

        /// <summary>
        /// Calls the calculation of the hydrodynamics
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        public override void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double relativeParticleMass, double dt) {
            double[] tempForces = CalculateHydrodynamicForces(U, P, LsTrk, CutCells_P, muA, relativeParticleMass, dt);
            double tempTorque = CalculateHydrodynamicTorque(U, P, LsTrk, CutCells_P, muA, dt);
            HydrodynamicsPostprocessing(tempForces, tempTorque);
        }

        /// <summary>
        /// Calls the calculation of the hydrodynamics
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        public void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double relativeParticleMass, double dt, double collisionTimestep) {
            double[] tempForces = CalculateHydrodynamicForces(U, P, LsTrk, CutCells_P, muA, relativeParticleMass, (dt - collisionTimestep));
            double tempTorque = CalculateHydrodynamicTorque(U, P, LsTrk, CutCells_P, muA, (dt - collisionTimestep));
            HydrodynamicsPostprocessing(tempForces, tempTorque);
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        protected override double[] CalculateHydrodynamicForces(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double relativeParticleMass, double dt) {
            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;
            double[] tempForces = ForcesIntegration(UA, pA, LsTrk, CutCells_P, RequiredOrder, muA);
            Force_MPISum(ref tempForces);
            for (int d = 0; d < spatialDim; d++) {
                tempForces[d] += relativeParticleMass * m_Gravity[d];
            }
            ForceAddedDamping(ref tempForces, dt);
            return tempForces;
        }

        private void ForceAddedDamping(ref double[] forces, double dt) {
            for (int d = 0; d < spatialDim; d++) {
                forces[d] += m_AddedDampingCoefficient * dt * (addedDampingTensor[0, d] * translationalAcceleration[0][0] + addedDampingTensor[1, d] * translationalAcceleration[0][1] + addedDampingTensor[d, 2] * rotationalAcceleration[0]);
            }
        }

        protected override double CalculateHydrodynamicTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double dt) {
            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;
            double tempTorque = TorqueIntegration(UA, pA, LsTrk, CutCells_P, RequiredOrder, muA);
            Torque_MPISum(ref tempTorque);
            TorqueAddedDamping(ref tempTorque, dt);
            return tempTorque;
        }

        private void TorqueAddedDamping(ref double torque, double dt) {
            torque += m_AddedDampingCoefficient * dt * (addedDampingTensor[2, 0] * translationalAcceleration[0][0] + addedDampingTensor[2, 1] * translationalAcceleration[0][1] + addedDampingTensor[2, 2] * rotationalAcceleration[0]);
        }

        protected override void HydrodynamicsPostprocessing(double[] tempForces, double tempTorque) {
            Underrelaxation.CalculateAverageForces(tempForces, tempTorque, m_MaxParticleLengthScale, out double averagedForces);
            Underrelaxation.Forces(ref tempForces, forcesPrevIteration, m_ForceAndTorqueConvergence, m_UnderrelaxationFactor, m_UseAddaptiveUnderrelaxation, averagedForces);
            Underrelaxation.Torque(ref tempTorque, torquePrevIteration, m_ForceAndTorqueConvergence, m_UnderrelaxationFactor, m_UseAddaptiveUnderrelaxation, averagedForces);
            ForceClearSmallValues(tempForces);
            TorqueClearSmallValues(tempTorque);
            Aux.TestArithmeticException(hydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(hydrodynamicTorque[0], "hydrodynamic torque");
        }

        private void ForceClearSmallValues(double[] tempForces) {
            for (int d = 0; d < spatialDim; d++) {
                hydrodynamicForces[0][d] = 0;
                if (Math.Abs(tempForces[d]) > m_ForceAndTorqueConvergence * 1e-2)
                    hydrodynamicForces[0][d] = tempForces[d];
            }
        }

        private void TorqueClearSmallValues(double tempTorque) {
            hydrodynamicTorque[0] = 0;
            if (Math.Abs(tempTorque) > m_ForceAndTorqueConvergence * 1e-2)
                hydrodynamicTorque[0] = tempTorque;
        }
    }
}

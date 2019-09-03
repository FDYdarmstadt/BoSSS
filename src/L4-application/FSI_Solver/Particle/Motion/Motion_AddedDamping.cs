﻿using BoSSS.Foundation;
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
            double underrelaxationFactor = 1,
            double forceAndTorqueConvergence = 1e-10,
            bool useAddaptiveUnderrelaxation = true,
            double addedDampingCoefficient = 1) : base(gravity, underrelaxationFactor, forceAndTorqueConvergence, useAddaptiveUnderrelaxation) {
            m_StartingAngle = angle[0];
            m_AddedDampingCoefficient = addedDampingCoefficient;
        }

        readonly ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation();
        readonly ParticleAddedDamping AddedDamping = new ParticleAddedDamping();

        /// <summary>
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        public double[,] addedDampingTensor = new double[6, 6];

        /// <summary>
        /// Added damping coefficient, should be between 0.5 and 1.5, for reference: Banks et.al. 2017
        /// </summary>
        private double m_AddedDampingCoefficient;

        /// <summary>
        /// Saving the initial angle of the particle for <see cref="UpdateDampingTensors()"/>
        /// </summary>
        private double m_StartingAngle;

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

        protected override void CalculateTranslationalAcceleration(double dt) {
            double[,] coefficientMatrix = CalculateCoefficientMatrix(dt);
            double denominator = CalculateDenominator(coefficientMatrix);

            double[] tempAcceleration = new double[2];
            tempAcceleration[0] = hydrodynamicForces[0][0] * (coefficientMatrix[1, 1] * coefficientMatrix[2, 2] - coefficientMatrix[1, 2] * coefficientMatrix[2, 1]);
            tempAcceleration[0] += hydrodynamicForces[0][1] * (-coefficientMatrix[0, 1] * coefficientMatrix[2, 2] + coefficientMatrix[0, 2] * coefficientMatrix[2, 1]);
            tempAcceleration[0] += hydrodynamicTorque[0] * (coefficientMatrix[0, 1] * coefficientMatrix[1, 2] - coefficientMatrix[0, 2] * coefficientMatrix[1, 1]);
            tempAcceleration[0] = tempAcceleration[0] / denominator;

            tempAcceleration[1] = hydrodynamicForces[0][0] * (-coefficientMatrix[1, 0] * coefficientMatrix[2, 2] + coefficientMatrix[1, 2] * coefficientMatrix[2, 0]);
            tempAcceleration[1] += hydrodynamicForces[0][1] * (coefficientMatrix[0, 0] * coefficientMatrix[2, 2] - coefficientMatrix[0, 2] * coefficientMatrix[2, 0]);
            tempAcceleration[1] += hydrodynamicTorque[0] * (-coefficientMatrix[0, 0] * coefficientMatrix[1, 2] + coefficientMatrix[0, 2] * coefficientMatrix[1, 0]);
            tempAcceleration[1] = tempAcceleration[1] / denominator;
            translationalAcceleration[0] = tempAcceleration.CloneAs();
        }

        protected override void CalculateRotationalAcceleration(double dt) {
            double[,] coefficientMatrix = CalculateCoefficientMatrix(dt);
            double denominator = CalculateDenominator(coefficientMatrix);

            double tempAcceleration = hydrodynamicForces[0][0] * (coefficientMatrix[1, 0] * coefficientMatrix[2, 1] - coefficientMatrix[1, 1] * coefficientMatrix[2, 0]);
            tempAcceleration += hydrodynamicForces[0][1] * (coefficientMatrix[0, 1] * coefficientMatrix[2, 0] - coefficientMatrix[0, 0] * coefficientMatrix[2, 1]);
            tempAcceleration += hydrodynamicTorque[0] * (coefficientMatrix[0, 0] * coefficientMatrix[1, 1] - coefficientMatrix[0, 1] * coefficientMatrix[1, 0]);
            rotationalAcceleration[0] = tempAcceleration / denominator;
        }

        private double[,] CalculateCoefficientMatrix(double Timestep) {
            double[,] massMatrix = GetMassMatrix();
            double[,] coefficientMatrix = massMatrix.CloneAs();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    coefficientMatrix[i, j] = massMatrix[i, j] + Timestep * m_AddedDampingCoefficient * addedDampingTensor[i, j];
                }
            }
            return coefficientMatrix;
        }
        private double[,] GetMassMatrix() {
            double[,] MassMatrix = new double[3, 3];
            MassMatrix[0, 0] = MassMatrix[1, 1] = particleMass;
            MassMatrix[2, 2] = momentOfInertia;
            return MassMatrix;
        }

        private double CalculateDenominator(double[,] coefficientMatrix) {
            double denominator = coefficientMatrix[0, 0] * coefficientMatrix[1, 1] * coefficientMatrix[2, 2];
            denominator -= coefficientMatrix[0, 0] * coefficientMatrix[1, 2] * coefficientMatrix[2, 1];
            denominator -= coefficientMatrix[0, 1] * coefficientMatrix[1, 0] * coefficientMatrix[2, 2];
            denominator += coefficientMatrix[0, 1] * coefficientMatrix[1, 2] * coefficientMatrix[2, 0];
            denominator += coefficientMatrix[0, 2] * coefficientMatrix[1, 0] * coefficientMatrix[2, 1];
            denominator -= coefficientMatrix[0, 2] * coefficientMatrix[1, 1] * coefficientMatrix[2, 0];
            return denominator;
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

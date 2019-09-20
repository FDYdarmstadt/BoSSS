/* =======================================================================
Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.FSI_Solver {
    public class Motion_AddedDamping : Motion_Wet {
        public Motion_AddedDamping(
            double[] gravity,
            double density,
            ParticleUnderrelaxationParam underrelaxationParam,
            double addedDampingCoefficient = 1) : base(gravity, density, underrelaxationParam) {
            m_StartingAngle = Angle[0];
            m_AddedDampingCoefficient = addedDampingCoefficient;
            UseAddedDamping = true;
        }

        /// <summary>
        /// Use added damping?, for reference: Banks et.al. 2017
        /// </summary>
        public override bool UseAddedDamping { get; } = true;

        private readonly ParticleAddedDamping AddedDamping = new ParticleAddedDamping();

        /// <summary>
        /// Added damping coefficient, should be between 0.5 and 1.5, for reference: Banks et.al. 2017
        /// </summary>
        private readonly double m_AddedDampingCoefficient;

        /// <summary>
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        private double[,] m_AddedDampingTensor;

        /// <summary>
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        public override double[,] AddedDampingTensor { get => m_AddedDampingTensor; } 

        /// <summary>
        /// Saving the initial angle of the particle for <see cref="UpdateDampingTensors()"/>
        /// </summary>
        private readonly double m_StartingAngle;

        /// <summary>
        /// Calculate tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public override void CalculateDampingTensor(Particle particle, LevelSetTracker LsTrk, double muA, double rhoA, double dt) {
            m_AddedDampingTensor = AddedDamping.IntegrationOverLevelSet(particle, LsTrk, muA, rhoA, dt, Position[0]);
            Aux.TestArithmeticException(m_AddedDampingTensor, "particle added damping tensor");
        }

        /// <summary>
        /// Update in every timestep tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public override void UpdateDampingTensors() {
            m_AddedDampingTensor = AddedDamping.RotateTensor(Angle[0], m_StartingAngle, AddedDampingTensor);
            Aux.TestArithmeticException(m_AddedDampingTensor, "particle added damping tensor");
        }

        /// <summary>
        /// MPI exchange of the damping tensors
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles
        /// </param>
        public override void ExchangeAddedDampingTensors() {
            int NoOfVars = 3;
            double[] StateBuffer = new double[NoOfVars * NoOfVars];
            for (int i = 0; i < NoOfVars; i++) {
                for (int j = 0; j < NoOfVars; j++) {
                    StateBuffer[i + NoOfVars * j] = m_AddedDampingTensor[i, j];
                }
            }
            double[] GlobalStateBuffer = StateBuffer.MPISum();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    m_AddedDampingTensor[i, j] = GlobalStateBuffer[i + NoOfVars * j];
                }
            }
        }

        protected override double[] CalculateTranslationalAcceleration(double dt) {
            double[,] coefficientMatrix = CalculateCoefficientMatrix(dt);
            double denominator = CalculateDenominator(coefficientMatrix);

            double[] l_Acceleration = new double[2];
            l_Acceleration[0] = HydrodynamicForces[0][0] * (coefficientMatrix[1, 1] * coefficientMatrix[2, 2] - coefficientMatrix[1, 2] * coefficientMatrix[2, 1]);
            l_Acceleration[0] += HydrodynamicForces[0][1] * (-coefficientMatrix[0, 1] * coefficientMatrix[2, 2] + coefficientMatrix[0, 2] * coefficientMatrix[2, 1]);
            l_Acceleration[0] += HydrodynamicTorque[0] * (coefficientMatrix[0, 1] * coefficientMatrix[1, 2] - coefficientMatrix[0, 2] * coefficientMatrix[1, 1]);
            l_Acceleration[0] = l_Acceleration[0] / denominator;

            l_Acceleration[1] = HydrodynamicForces[0][0] * (-coefficientMatrix[1, 0] * coefficientMatrix[2, 2] + coefficientMatrix[1, 2] * coefficientMatrix[2, 0]);
            l_Acceleration[1] += HydrodynamicForces[0][1] * (coefficientMatrix[0, 0] * coefficientMatrix[2, 2] - coefficientMatrix[0, 2] * coefficientMatrix[2, 0]);
            l_Acceleration[1] += HydrodynamicTorque[0] * (-coefficientMatrix[0, 0] * coefficientMatrix[1, 2] + coefficientMatrix[0, 2] * coefficientMatrix[1, 0]);
            l_Acceleration[1] = l_Acceleration[1] / denominator;
            Aux.TestArithmeticException(l_Acceleration, "particle translational acceleration");
            return l_Acceleration;
        }

        protected override double CalculateRotationalAcceleration(double dt) {
            double[,] coefficientMatrix = CalculateCoefficientMatrix(dt);
            double denominator = CalculateDenominator(coefficientMatrix);

            double l_Acceleration = HydrodynamicForces[0][0] * (coefficientMatrix[1, 0] * coefficientMatrix[2, 1] - coefficientMatrix[1, 1] * coefficientMatrix[2, 0]);
            l_Acceleration += HydrodynamicForces[0][1] * (coefficientMatrix[0, 1] * coefficientMatrix[2, 0] - coefficientMatrix[0, 0] * coefficientMatrix[2, 1]);
            l_Acceleration += HydrodynamicTorque[0] * (coefficientMatrix[0, 0] * coefficientMatrix[1, 1] - coefficientMatrix[0, 1] * coefficientMatrix[1, 0]);
            l_Acceleration /= denominator;
            Aux.TestArithmeticException(l_Acceleration, "particle rotational acceleration");
            return l_Acceleration;
        }

        private double[,] CalculateCoefficientMatrix(double Timestep) {
            double[,] massMatrix = GetMassMatrix();
            double[,] coefficientMatrix = massMatrix.CloneAs();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    coefficientMatrix[i, j] = massMatrix[i, j] + Timestep * m_AddedDampingCoefficient * AddedDampingTensor[i, j];
                }
            }
            return coefficientMatrix;
        }
        private double[,] GetMassMatrix() {
            double[,] MassMatrix = new double[3, 3];
            MassMatrix[0, 0] = MassMatrix[1, 1] = ParticleArea * Density;
            MassMatrix[2, 2] = MomentOfInertia;
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
        /// Calls the calculation of the hydrodynamics
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="fluidViscosity"></param>
        public override void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double fluidViscosity, double fluidDensity, bool firstIteration, double dt) {
            double[] tempForces = CalculateHydrodynamicForces(U, P, LsTrk, CutCells_P, fluidViscosity, fluidDensity, dt);
            double tempTorque = CalculateHydrodynamicTorque(U, P, LsTrk, CutCells_P, fluidViscosity, dt);
            HydrodynamicsPostprocessing(tempForces, tempTorque, firstIteration);
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        protected override double[] CalculateHydrodynamicForces(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double fluidDensity, double dt) {
            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;
            double[] tempForces = ForcesIntegration(UA, pA, LsTrk, CutCells_P, RequiredOrder, muA);
            Force_MPISum(ref tempForces);
            for (int d = 0; d < spatialDim; d++) {
                tempForces[d] += (Density - fluidDensity) * ParticleArea * Gravity[d];
            }
            ForceAddedDamping(ref tempForces, dt);
            return tempForces;
        }

        private void ForceAddedDamping(ref double[] forces, double dt) {
            for (int d = 0; d < spatialDim; d++) {
                forces[d] += m_AddedDampingCoefficient * dt * (AddedDampingTensor[0, d] * TranslationalAcceleration[0][0] + AddedDampingTensor[1, d] * TranslationalAcceleration[0][1] + AddedDampingTensor[d, 2] * RotationalAcceleration[0]);
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
            torque += m_AddedDampingCoefficient * dt * (AddedDampingTensor[2, 0] * TranslationalAcceleration[0][0] + AddedDampingTensor[2, 1] * TranslationalAcceleration[0][1] + AddedDampingTensor[2, 2] * RotationalAcceleration[0]);
        }
    }
}

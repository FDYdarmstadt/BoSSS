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
using System;

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

        private readonly ParticleAddedDamping AddedDamping = new ParticleAddedDamping();

        /// <summary>
        /// Saving the initial angle of the particle for <see cref="UpdateDampingTensors()"/>
        /// </summary>
        private readonly double m_StartingAngle;

        /// <summary>
        /// Calculate tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public override void CalculateDampingTensor(Particle particle, LevelSetTracker LsTrk, double muA, double rhoA, double dt) {
            addedDampingTensor = AddedDamping.IntegrationOverLevelSet(particle, LsTrk, muA, rhoA, dt, Position[0]);
            Aux.TestArithmeticException(addedDampingTensor, "particle added damping tensor");
        }

        /// <summary>
        /// Update in every timestep tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public override void UpdateDampingTensors() {
            addedDampingTensor = AddedDamping.RotateTensor(Angle[0], m_StartingAngle, addedDampingTensor);
            Aux.TestArithmeticException(addedDampingTensor, "particle added damping tensor");
        }

        protected override void CalculateTranslationalAcceleration(double dt) {
            double[,] coefficientMatrix = CalculateCoefficientMatrix(dt);
            double denominator = CalculateDenominator(coefficientMatrix);

            double[] tempAcceleration = new double[2];
            tempAcceleration[0] = HydrodynamicForces[0][0] * (coefficientMatrix[1, 1] * coefficientMatrix[2, 2] - coefficientMatrix[1, 2] * coefficientMatrix[2, 1]);
            tempAcceleration[0] += HydrodynamicForces[0][1] * (-coefficientMatrix[0, 1] * coefficientMatrix[2, 2] + coefficientMatrix[0, 2] * coefficientMatrix[2, 1]);
            tempAcceleration[0] += HydrodynamicTorque[0] * (coefficientMatrix[0, 1] * coefficientMatrix[1, 2] - coefficientMatrix[0, 2] * coefficientMatrix[1, 1]);
            tempAcceleration[0] = tempAcceleration[0] / denominator;

            tempAcceleration[1] = HydrodynamicForces[0][0] * (-coefficientMatrix[1, 0] * coefficientMatrix[2, 2] + coefficientMatrix[1, 2] * coefficientMatrix[2, 0]);
            tempAcceleration[1] += HydrodynamicForces[0][1] * (coefficientMatrix[0, 0] * coefficientMatrix[2, 2] - coefficientMatrix[0, 2] * coefficientMatrix[2, 0]);
            tempAcceleration[1] += HydrodynamicTorque[0] * (-coefficientMatrix[0, 0] * coefficientMatrix[1, 2] + coefficientMatrix[0, 2] * coefficientMatrix[1, 0]);
            tempAcceleration[1] = tempAcceleration[1] / denominator;
            translationalAcceleration[0] = tempAcceleration.CloneAs();
            Aux.TestArithmeticException(translationalAcceleration[0], "particle translational acceleration");
        }

        protected override void CalculateRotationalAcceleration(double dt) {
            double[,] coefficientMatrix = CalculateCoefficientMatrix(dt);
            double denominator = CalculateDenominator(coefficientMatrix);

            double tempAcceleration = HydrodynamicForces[0][0] * (coefficientMatrix[1, 0] * coefficientMatrix[2, 1] - coefficientMatrix[1, 1] * coefficientMatrix[2, 0]);
            tempAcceleration += HydrodynamicForces[0][1] * (coefficientMatrix[0, 1] * coefficientMatrix[2, 0] - coefficientMatrix[0, 0] * coefficientMatrix[2, 1]);
            tempAcceleration += HydrodynamicTorque[0] * (coefficientMatrix[0, 0] * coefficientMatrix[1, 1] - coefficientMatrix[0, 1] * coefficientMatrix[1, 0]);
            rotationalAcceleration[0] = tempAcceleration / denominator;
            Aux.TestArithmeticException(rotationalAcceleration[0], "particle rotational acceleration");
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
            MassMatrix[0, 0] = MassMatrix[1, 1] = particleArea * Density;
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
                tempForces[d] += (Density - fluidDensity) * particleArea * m_Gravity[d];
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
    }
}

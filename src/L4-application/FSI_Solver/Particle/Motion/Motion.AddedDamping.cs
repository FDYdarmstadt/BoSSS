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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using FSI_Solver;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    [Serializable]
    public class MotionAddedDamping : Motion {

        /// <summary>
        /// The added damping description of motion including hydrodynamics, for reference: Banks et.al. 2017.
        /// </summary>
        /// <param name="gravity">
        /// The gravity (volume forces) acting on the particle.
        /// </param>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        /// <param name="underrelaxationParam">
        /// The underrelaxation parameters (convergence limit, prefactor and a bool whether to use addaptive underrelaxation) defined in <see cref="ParticleUnderrelaxationParam"/>.
        /// </param>
        /// <param name="addedDampingCoefficient">
        /// The added damping coefficient is a scaling factor for the model. Should be between 0.5 and 1.5, for reference: Banks et.al. 2017.
        /// </param>
        public MotionAddedDamping(Vector gravity, double density, double addedDampingCoefficient = 1) : base(gravity, density) {
            m_StartingAngle = GetAngle(0);
            m_AddedDampingCoefficient = addedDampingCoefficient;    
            UseAddedDamping = true;
        }

        [NonSerialized]
        private ParticleAddedDamping AddedDamping = new ParticleAddedDamping();
        [DataMember]
        private double m_AddedDampingCoefficient = 1;
        [DataMember]
        private readonly double m_StartingAngle;

        /// <summary>
        /// We are using the added damping model, for reference: Banks et.al. 2017.
        /// </summary>
        [DataMember]
        internal override bool UseAddedDamping { get; } = true;
        
        /// <summary>
        /// Calculate the tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        /// <param name="particle"></param>
        /// <param name="levelSetTracker"></param>
        /// <param name="fluidViscosity"></param>
        /// <param name="fluidDensity"></param>
        /// <param name="dt"></param>
        internal override void CalculateDampingTensor(Particle particle, LevelSetTracker levelSetTracker, double fluidViscosity, double fluidDensity, double dt) {
            using (new FuncTrace()) {
                AddedDamping = new ParticleAddedDamping();
                Aux = new FSIAuxillary();
                AddedDampingTensor = AddedDamping.IntegrationOverLevelSet(particle, levelSetTracker, fluidViscosity, fluidDensity, dt, GetPosition(0));
                Aux.TestArithmeticException(AddedDampingTensor, "particle added damping tensor");
            }
        }

        /// <summary>
        /// Update in every timestep tensors to implement the added damping model (Banks et.al. 2017).
        /// </summary>
        internal override void UpdateDampingTensors() {
            using (new FuncTrace()) {
                AddedDampingTensor = AddedDamping.RotateTensor(GetAngle(0), m_StartingAngle, AddedDampingTensor);
                Aux.TestArithmeticException(AddedDampingTensor, "particle added damping tensor");
            }
        }

        /// <summary>
        /// MPI exchange of the damping tensors
        /// </summary>
        internal override void ExchangeAddedDampingTensors() {
            using (new FuncTrace()) {
                int NoOfVars = 3;
                double[] StateBuffer = new double[NoOfVars * NoOfVars];
                for (int i = 0; i < NoOfVars; i++) {
                    for (int j = 0; j < NoOfVars; j++) {
                        StateBuffer[i + NoOfVars * j] = AddedDampingTensor[i, j];
                    }
                }
                double[] GlobalStateBuffer = StateBuffer.MPISum();
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        AddedDampingTensor[i, j] = GlobalStateBuffer[i + NoOfVars * j];
                    }
                }
            }
        }

        /// <summary>
        /// Calculates the translational acceleration of the particle using the added damping model.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override Vector CalculateTranslationalAcceleration(double dt) {
            using (new FuncTrace()) {
                double[,] coefficientMatrix = CalculateCoefficientMatrix(dt);
                double denominator = CalculateDenominator(coefficientMatrix);

                Vector l_Acceleration = new Vector(SpatialDim);
                l_Acceleration[0] = GetHydrodynamicForces(0)[0] * (coefficientMatrix[1, 1] * coefficientMatrix[2, 2] - coefficientMatrix[1, 2] * coefficientMatrix[2, 1]);
                l_Acceleration[0] += GetHydrodynamicForces(0)[1] * (-coefficientMatrix[0, 1] * coefficientMatrix[2, 2] + coefficientMatrix[0, 2] * coefficientMatrix[2, 1]);
                l_Acceleration[0] += GetHydrodynamicTorque(0) * (coefficientMatrix[0, 1] * coefficientMatrix[1, 2] - coefficientMatrix[0, 2] * coefficientMatrix[1, 1]);
                l_Acceleration[0] = l_Acceleration[0] / denominator;

                l_Acceleration[1] = GetHydrodynamicForces(0)[0] * (-coefficientMatrix[1, 0] * coefficientMatrix[2, 2] + coefficientMatrix[1, 2] * coefficientMatrix[2, 0]);
                l_Acceleration[1] += GetHydrodynamicForces(0)[1] * (coefficientMatrix[0, 0] * coefficientMatrix[2, 2] - coefficientMatrix[0, 2] * coefficientMatrix[2, 0]);
                l_Acceleration[1] += GetHydrodynamicTorque(0) * (-coefficientMatrix[0, 0] * coefficientMatrix[1, 2] + coefficientMatrix[0, 2] * coefficientMatrix[1, 0]);
                l_Acceleration[1] = l_Acceleration[1] / denominator;
                Aux.TestArithmeticException(l_Acceleration, "particle translational acceleration");
                return l_Acceleration;
            }
        }

        /// <summary>
        /// Calculates the rotational acceleration of the particle using the added damping model.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected override double CalculateRotationalAcceleration(double dt) {
            using (new FuncTrace()) {
                double[,] coefficientMatrix = CalculateCoefficientMatrix(dt);
                double denominator = CalculateDenominator(coefficientMatrix);

                double l_Acceleration = GetHydrodynamicForces(0)[0] * (coefficientMatrix[1, 0] * coefficientMatrix[2, 1] - coefficientMatrix[1, 1] * coefficientMatrix[2, 0]);
                l_Acceleration += GetHydrodynamicForces(0)[1] * (coefficientMatrix[0, 1] * coefficientMatrix[2, 0] - coefficientMatrix[0, 0] * coefficientMatrix[2, 1]);
                l_Acceleration += GetHydrodynamicTorque(0) * (coefficientMatrix[0, 0] * coefficientMatrix[1, 1] - coefficientMatrix[0, 1] * coefficientMatrix[1, 0]);
                l_Acceleration /= denominator;
                Aux.TestArithmeticException(l_Acceleration, "particle rotational acceleration");
                return l_Acceleration;
            }
        }

        /// <summary>
        /// Calculates the coefficient matrix for the acceleration constituted of the mass matrix and the added damping tensor.
        /// </summary>
        /// <param name="dt">Timestep</param>
        private double[,] CalculateCoefficientMatrix(double dt) {
            using (new FuncTrace()) {
                double[,] massMatrix = GetMassMatrix();
                double[,] coefficientMatrix = massMatrix.CloneAs();
                if (m_AddedDampingCoefficient == 0)// somehow the coefficient isnt safed during restart
                    m_AddedDampingCoefficient = 1;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        coefficientMatrix[i, j] = massMatrix[i, j] + dt * m_AddedDampingCoefficient * AddedDampingTensor[i, j];
                    }
                }
                return coefficientMatrix;
            }
        }

        /// <summary>
        /// Calculates the mass matrix of the particle.
        /// </summary>
        private double[,] GetMassMatrix() {
            using (new FuncTrace()) {
                double[,] MassMatrix = new double[3, 3];
                MassMatrix[0, 0] = MassMatrix[1, 1] = ParticleArea * Density;
                MassMatrix[2, 2] = MomentOfInertia;
                return MassMatrix;
            }
        }

        /// <summary>
        /// Calculates the denominator necessary for the calculation of the acceleration of the particle.
        /// </summary>
        /// <param name="coefficientMatrix">The matrix calculated in <see cref="CalculateCoefficientMatrix"></see>/></param>
        private double CalculateDenominator(double[,] coefficientMatrix) {
            using (new FuncTrace()) {
                double denominator = coefficientMatrix[0, 0] * coefficientMatrix[1, 1] * coefficientMatrix[2, 2];
                denominator -= coefficientMatrix[0, 0] * coefficientMatrix[1, 2] * coefficientMatrix[2, 1];
                denominator -= coefficientMatrix[0, 1] * coefficientMatrix[1, 0] * coefficientMatrix[2, 2];
                denominator += coefficientMatrix[0, 1] * coefficientMatrix[1, 2] * coefficientMatrix[2, 0];
                denominator += coefficientMatrix[0, 2] * coefficientMatrix[1, 0] * coefficientMatrix[2, 1];
                denominator -= coefficientMatrix[0, 2] * coefficientMatrix[1, 1] * coefficientMatrix[2, 0];
                return denominator;
            }
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        public override Vector CalculateHydrodynamicForces(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, double dt) {
            using (new FuncTrace()) {
                Vector tempForces = new Vector(hydrodynamicsIntegration.Forces(cutCells));
                Aux.TestArithmeticException(tempForces, "temporal forces during calculation of hydrodynamics");
                tempForces = ForceAddedDamping(tempForces, dt);
                return tempForces;
            }
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="levelSetTracker"></param>
        /// <param name="fluidViscosity"></param>
        /// <param name="cutCells"></param>
        /// <param name="dt"></param>
        public override double CalculateHydrodynamicTorque(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, double dt) {
            using (new FuncTrace()) {
                double tempTorque = hydrodynamicsIntegration.Torque(GetPosition(0), cutCells);
                Aux.TestArithmeticException(tempTorque, "temporal torque during calculation of hydrodynamics");
                TorqueAddedDamping(ref tempTorque, dt);
                return tempTorque;
            }
        }

        /// <summary>
        /// Calculates the added damping effects on the hydrodynamic forces
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="forces"></param>
        private Vector ForceAddedDamping(Vector forces, double dt) {
            using (new FuncTrace()) {
                for (int d = 0; d < SpatialDim; d++) {
                    forces[d] += m_AddedDampingCoefficient * dt * (AddedDampingTensor[0, d] * GetTranslationalAcceleration(0)[0] + AddedDampingTensor[1, d] * GetTranslationalAcceleration(0)[1] + AddedDampingTensor[d, 2] * GetRotationalAcceleration(0));
                }
                return forces;
            }
        }

        /// <summary>
        /// Calculates the added damping effects on the hydrodynamic torque.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="torque"></param>
        private void TorqueAddedDamping(ref double torque, double dt) {
            using (new FuncTrace()) {
                torque += m_AddedDampingCoefficient * dt * (AddedDampingTensor[2, 0] * GetTranslationalAcceleration(0)[0] + AddedDampingTensor[2, 1] * GetTranslationalAcceleration(0)[1] + AddedDampingTensor[2, 2] * GetRotationalAcceleration(0));
            }
        }

        public override object Clone() {
            using (new FuncTrace()) {
                Motion clonedMotion = new MotionAddedDamping(Gravity, Density, m_AddedDampingCoefficient);
                clonedMotion.SetParticleArea(ParticleArea);
                clonedMotion.SetParticleMomentOfInertia(MomentOfInertia);
                return clonedMotion;
            }
        }
    }
}

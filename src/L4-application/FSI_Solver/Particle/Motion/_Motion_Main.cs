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
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using FSI_Solver;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {

    public class Motion_Wet {

        private const int historyLength = 4;
        protected static int spatialDim = 2;
        [NonSerialized]
        readonly internal FSI_Auxillary Aux = new FSI_Auxillary();
        [NonSerialized]
        readonly private ParticleForceIntegration ForceIntegration = new ParticleForceIntegration();
        readonly ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation();

        public Motion_Wet(double[] gravity,
            ParticleUnderrelaxationParam underrelaxationParam = null) {
            m_Gravity = gravity;
            m_UnderrelaxationParam = underrelaxationParam;

            for (int i = 0; i < historyLength; i++) {
                position.Add(new double[spatialDim]);
                angle.Add(new double());
                translationalVelocity.Add(new double[spatialDim]);
                translationalAcceleration.Add(new double[spatialDim]);
                rotationalVelocity.Add(new double());
                rotationalAcceleration.Add(new double());
                hydrodynamicForces.Add(new double[spatialDim]);
                hydrodynamicTorque.Add(new double());
            }
        }

        ParticleUnderrelaxationParam m_UnderrelaxationParam = null;



        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        protected double[] m_Gravity;

        /// <summary>
        /// Added damping coefficient, should be between 0.5 and 1.5, for reference: Banks et.al. 2017
        /// </summary>
        public bool includeRotation = true;

        /// <summary>
        /// Added damping coefficient, should be between 0.5 and 1.5, for reference: Banks et.al. 2017
        /// </summary>
        public bool includeTranslation = true;

        /// <summary>
        /// Added damping coefficient, should be between 0.5 and 1.5, for reference: Banks et.al. 2017
        /// </summary>
        public bool useAddedDamping = false;

        /// <summary>
        /// Added damping coefficient, should be between 0.5 and 1.5, for reference: Banks et.al. 2017
        /// </summary>
        public double m_AddedDampingCoefficient = -1;

        /// <summary>
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        public double[,] addedDampingTensor = new double[6, 6];

        /// <summary>
        /// Density of the particle.
        /// </summary>
        [DataMember]
        public double particleDensity = 1;

        /// <summary>
        /// Active velocity (alternative to active stress) on the current particle.
        /// </summary>
        [DataMember]
        public double ActiveVelocity = 0;

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        protected double particleArea = new double();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        protected double momentOfInertia = new double();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> position = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> angle = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> translationalVelocity = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> rotationalVelocity = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        public double[] PreCollisionVelocity;

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        protected List<double[]> translationalAcceleration = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        protected List<double> rotationalAcceleration = new List<double>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public List<double[]> hydrodynamicForces = new List<double[]>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public double[] forcesPrevIteration = new double[spatialDim];

        /// <summary>
        /// The Torque acting on the particle in the current time step.
        /// </summary>
        public List<double> hydrodynamicTorque = new List<double>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public double torquePrevIteration = new double();


        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public double collisionTimestep = new double();

        /// <summary>
        /// The maximum lenghtscale of the particle.
        /// </summary>
        protected double m_MaxParticleLengthScale;

        /// <summary>
        /// Saves position and angle of the last timestep
        /// </summary>
        public void SavePositionAndAngleOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(position);
            Aux.SaveValueOfLastTimestep(angle);
        }

        /// <summary>
        /// Saves translational and rotational velocities of the last timestep
        /// </summary>
        public void SaveVelocityOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(translationalVelocity);
            Aux.SaveValueOfLastTimestep(rotationalVelocity);
            Aux.SaveMultidimValueOfLastTimestep(translationalAcceleration);
            Aux.SaveValueOfLastTimestep(rotationalAcceleration);
        }

        /// <summary>
        /// Saves force and torque of the previous timestep
        /// </summary>
        public void SaveHydrodynamicsOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(hydrodynamicForces);
            Aux.SaveValueOfLastTimestep(hydrodynamicTorque);
        }

        public void InitializeParticlePositionAndAngle(double[] positionP, double angleP) {
            for (int i = 0; i < historyLength; i++) {
                position[i] = positionP.CloneAs();
                angle[i] = angleP * 2 * Math.PI / 360;
            }
        }

        public void InitializeParticleVelocity(double[] translationalVelocityP, double rotationalVelocityP) {
            if (translationalVelocityP == null)
                translationalVelocity[0] = new double[spatialDim];
            else
                translationalVelocity[0] = translationalVelocityP.CloneAs();

            rotationalVelocity[0] = rotationalVelocityP;
        }

        /// <summary>
        /// Init of the particle mass.
        /// </summary>
        /// <param name="mass"></param>
        public void GetParticleArea(double area) {
            particleArea = area;
        }

        /// <summary>
        /// Init of the particle mass.
        /// </summary>
        /// <param name="mass"></param>
        public void GetParticleDensity(double density) {
            particleDensity = density;
            Aux.TestArithmeticException(particleArea, "particle area");
            Aux.TestArithmeticException(particleDensity, "particle density");
        }

        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void GetParticleLengthscale(double lengthscale) {
            m_MaxParticleLengthScale = lengthscale;
        }

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        public double ParticleMass {
            get {
                Aux.TestArithmeticException(particleArea, "particle area");
                Aux.TestArithmeticException(particleDensity, "particle density");
                return particleArea * particleDensity;
            }
        }

        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void GetParticleMomentOfInertia(double moment) {
            momentOfInertia = moment;
        }

        /// <summary>
        /// Calls the calculation of the position and angle.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticlePositionAndAngle(double dt) {
            if (collisionTimestep == 0) {
                SavePositionAndAngleOfPreviousTimestep();
                CalculateParticlePosition(dt);
                CalculateParticleAngle(dt);
            }
            else {
                CalculateParticlePosition(dt, collisionTimestep);
                CalculateParticleAngle(dt, collisionTimestep);
            }
        }

        /// <summary>
        /// Calls the calculation of the position and angle.
        /// </summary>
        /// <param name="dt"></param>
        public void CollisionParticlePositionAndAngle(double collisionDynamicTimestep) {
            CalculateParticlePosition(collisionDynamicTimestep, collisionProcedure: true);
            CalculateParticleAngle(collisionDynamicTimestep, collisionProcedure: true);
        }

        /// <summary>
        /// Calls the calculation of the velocity.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticleVelocity(double dt) {
            if (collisionTimestep == 0) {
                CalculateTranslationalVelocity(dt);
                CalculateAngularVelocity(dt);
            }
            else {
                CalculateTranslationalVelocity(dt, collisionTimestep);// ich glaube das kann weg
                CalculateAngularVelocity(dt, collisionTimestep);
            }
        }

        /// <summary>
        /// Calls the calculation of the hydrodynamics
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="fluidViscosity"></param>
        public virtual void UpdateForcesAndTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double fluidViscosity, double fluidDensity, bool firstIteration, double dt = 0) {
            double[] tempForces = CalculateHydrodynamicForces(U, P, LsTrk, CutCells_P, fluidViscosity, fluidDensity);
            double tempTorque = CalculateHydrodynamicTorque(U, P, LsTrk, CutCells_P, fluidViscosity);
            HydrodynamicsPostprocessing(tempForces, tempTorque, firstIteration);
        }

        public virtual void PredictForceAndTorque(double activeStress, int TimestepInt) {
            if (TimestepInt == 1) {
                hydrodynamicForces[0][0] = 10 * Math.Cos(angle[0]) * activeStress + m_Gravity[1] * particleDensity * particleArea / 10;
                hydrodynamicForces[0][1] = 10 * Math.Sin(angle[0]) * activeStress + m_Gravity[1] * particleDensity * particleArea / 10;
                hydrodynamicTorque[0] = 0;
            }
            else {
                for (int d = 0; d < spatialDim; d++) {
                    hydrodynamicForces[0][d] = (hydrodynamicForces[1][d] + 4 * hydrodynamicForces[2][d] + hydrodynamicForces[3][d]) / 6;
                    if (Math.Abs(hydrodynamicForces[0][d]) < 1e-20)
                        hydrodynamicForces[0][d] = 0;
                }
                hydrodynamicTorque[0] = (hydrodynamicTorque[1] + 4 * hydrodynamicTorque[2] + hydrodynamicTorque[3]) / 6;
                if (Math.Abs(hydrodynamicTorque[0]) < 1e-20)
                    hydrodynamicTorque[0] = 0;
            }
            Aux.TestArithmeticException(hydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(hydrodynamicTorque[0], "hydrodynamic torque");
        }

        public virtual void CalculateDampingTensor(Particle particle, LevelSetTracker LsTrk, double muA, double rhoA, double dt) {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        public virtual void UpdateDampingTensors() {
            //throw new Exception("Predict force and torque should only be called if added damping is active.");
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected virtual void CalculateParticlePosition(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d] + (translationalVelocity[0][d] + 4 * translationalVelocity[1][d] + translationalVelocity[2][d]) * dt / 6;
            }

            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected virtual void CalculateParticlePosition(double dt, bool collisionProcedure) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[0][d] + (translationalVelocity[0][d] + 4 * translationalVelocity[1][d] + translationalVelocity[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual void CalculateParticlePosition(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                position[0][d] = position[1][d] + translationalVelocity[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(position[0], "particle position");
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected virtual void CalculateParticleAngle(double dt) {
            angle[0] = angle[1] + (rotationalVelocity[0] + 4 * rotationalVelocity[1] + rotationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected virtual void CalculateParticleAngle(double dt, bool collisionProcedure) {
            angle[0] = angle[0] + (rotationalVelocity[0] + 4 * rotationalVelocity[1] + rotationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new particle angle after a collision
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual void CalculateParticleAngle(double dt, double collisionTimestep) {
            angle[0] = angle[1] + rotationalVelocity[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(angle[0], "particle angle");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected virtual void CalculateTranslationalVelocity(double dt) {
            CalculateTranslationalAcceleration(dt);
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + (translationalAcceleration[0][d] + 4 * translationalAcceleration[1][d] + translationalAcceleration[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            CalculateTranslationalAcceleration(dt - collisionTimestep);
            for (int d = 0; d < spatialDim; d++) {
                translationalVelocity[0][d] = translationalVelocity[1][d] + translationalAcceleration[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(translationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected virtual void CalculateAngularVelocity(double dt) {
            CalculateRotationalAcceleration(dt);
            rotationalVelocity[0] = rotationalVelocity[1] + (rotationalAcceleration[0] + 4 * rotationalAcceleration[1] + rotationalAcceleration[2]) * dt / 6;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual void CalculateAngularVelocity(double dt, double collisionTimestep) {
            CalculateRotationalAcceleration(dt - collisionTimestep);
            rotationalVelocity[0] = rotationalVelocity[1] + rotationalAcceleration[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(rotationalVelocity[0], "particle rotational velocity");
        }

        internal void CalculateNormalAndTangentialVelocity(double[] NormalVector) {
            double[] Velocity = translationalVelocity[0];
            double[] TangentialVector = new double[] { -NormalVector[1], NormalVector[0] };
            PreCollisionVelocity = new double[] { Velocity[0] * NormalVector[0] + Velocity[1] * NormalVector[1], Velocity[0] * TangentialVector[0] + Velocity[1] * TangentialVector[1] };
            Aux.TestArithmeticException(PreCollisionVelocity, "particle velocity before collision");
        }

        /// <summary>
        /// Calculate the new acceleration
        /// </summary>
        /// <param name="dt"></param>
        protected virtual void CalculateTranslationalAcceleration(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                translationalAcceleration[0][d] = hydrodynamicForces[0][d] / (particleDensity * particleArea);
            }
            Aux.TestArithmeticException(translationalAcceleration[0], "particle translational acceleration");
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        protected virtual void CalculateRotationalAcceleration(double dt) {
            rotationalAcceleration[0] = hydrodynamicTorque[0] / momentOfInertia;
            Aux.TestArithmeticException(rotationalAcceleration[0], "particle rotational acceleration");
        }
                
        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="LsTrk"></param>
        /// <param name="muA"></param>
        protected virtual double[] CalculateHydrodynamicForces(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double fluidDensity, double dt = 0) {
            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;
            double[] tempForces = ForcesIntegration(UA, pA, LsTrk, CutCells_P, RequiredOrder, muA);
            Force_MPISum(ref tempForces);
            for (int d = 0; d < spatialDim; d++) {
                tempForces[d] += (particleDensity - fluidDensity) * particleArea * m_Gravity[d];
            }
            return tempForces;
        }

        protected double[] ForcesIntegration(SinglePhaseField[] UA, ConventionalDGField pA, LevelSetTracker LsTrk, CellMask CutCells_P, int RequiredOrder, double FluidViscosity) {
            double[] tempForces = new double[spatialDim];
            double[] IntegrationForces = tempForces.CloneAs();
            for (int d = 0; d < spatialDim; d++) {
                void ErrFunc(int CurrentCellID, int Length, NodeSet Ns, MultidimensionalArray result) {

                    int K = result.GetLength(1);
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Length, K, spatialDim, spatialDim);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Length, K);
                    MultidimensionalArray Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, CurrentCellID, Length);
                    for (int i = 0; i < spatialDim; i++) {
                        UA[i].EvaluateGradient(CurrentCellID, Length, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }
                    pA.Evaluate(CurrentCellID, Length, Ns, pARes);
                    for (int j = 0; j < Length; j++) {
                        for (int k = 0; k < K; k++) {
                            result[j, k] = ForceIntegration.CalculateStressTensor(Grad_UARes, pARes, Normals, FluidViscosity, k, j, spatialDim, d);
                        }
                    }
                }
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, CutCells_P);
                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat, cqs.Compile(LsTrk.GridDat, RequiredOrder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        IntegrationForces[d] = ForceTorqueSummationWithNeumaierArray(IntegrationForces[d], ResultsOfIntegration, Length);
                    }
                ).Execute();
            }
            return tempForces = IntegrationForces.CloneAs();
        }

        protected void Force_MPISum(ref double[] forces) {
            double[] stateBuffer = forces.CloneAs();
            double[] globalStateBuffer = stateBuffer.MPISum();
            for (int d = 0; d < spatialDim; d++) {
                forces[d] = globalStateBuffer[d];
            }
        }

        protected virtual double CalculateHydrodynamicTorque(VectorField<SinglePhaseField> U, SinglePhaseField P, LevelSetTracker LsTrk, CellMask CutCells_P, double muA, double dt = 0) {
            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;
            double tempTorque = TorqueIntegration(UA, pA, LsTrk, CutCells_P, RequiredOrder, muA);
            Torque_MPISum(ref tempTorque);
            return tempTorque;
        }

        protected double TorqueIntegration(SinglePhaseField[] UA, ConventionalDGField pA, LevelSetTracker LsTrk, CellMask CutCells_P, int RequiredOrder, double FluidViscosity) {
            double tempTorque = new double();
            void ErrFunc2(int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                int K = result.GetLength(1);
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, spatialDim, spatialDim); ;
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);
                for (int i = 0; i < spatialDim; i++) {
                    UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                }
                pA.Evaluate(j0, Len, Ns, pARes);
                for (int j = 0; j < Len; j++) {
                    MultidimensionalArray Ns_Global = Ns.CloneAs();
                    LsTrk.GridDat.TransformLocal2Global(Ns, Ns_Global, j0 + j);
                    for (int k = 0; k < K; k++) {
                        result[j, k] = ForceIntegration.CalculateTorqueFromStressTensor2D(Grad_UARes, pARes, Normals, Ns_Global, FluidViscosity, k, j, position[0]);
                    }
                }
            }
            var SchemeHelper2 = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs2 = SchemeHelper2.GetLevelSetquadScheme(0, CutCells_P);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat, cqs2.Compile(LsTrk.GridDat, RequiredOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc2(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    tempTorque = ForceTorqueSummationWithNeumaierArray(tempTorque, ResultsOfIntegration, Length);
                }
            ).Execute();
            return tempTorque;
        }

        protected void Torque_MPISum(ref double torque) {
            double stateBuffer = torque;
            double globalStateBuffer = stateBuffer.MPISum();
            torque = globalStateBuffer;
        }

        protected virtual void HydrodynamicsPostprocessing(double[] tempForces, double tempTorque, bool firstIteration) {
            double forceAndTorqueConvergence = 0;
            if (m_UnderrelaxationParam != null && !firstIteration) {
                forceAndTorqueConvergence = m_UnderrelaxationParam.hydroDynConvergenceLimit;
                double underrelaxationFactor = m_UnderrelaxationParam.underrelaxationFactor;
                bool useAddaptiveUnderrelaxation = m_UnderrelaxationParam.useAddaptiveUnderrelaxation;
                Underrelaxation.CalculateAverageForces(tempForces, tempTorque, m_MaxParticleLengthScale, out double averagedForces);
                Underrelaxation.Forces(ref tempForces, forcesPrevIteration, forceAndTorqueConvergence, underrelaxationFactor, useAddaptiveUnderrelaxation, averagedForces);
                Underrelaxation.Torque(ref tempTorque, torquePrevIteration, forceAndTorqueConvergence, underrelaxationFactor, useAddaptiveUnderrelaxation, averagedForces);
            }
            ForceClearSmallValues(tempForces, forceAndTorqueConvergence);
            TorqueClearSmallValues(tempTorque, forceAndTorqueConvergence);
            Aux.TestArithmeticException(hydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(hydrodynamicTorque[0], "hydrodynamic torque");
        }

        private void ForceClearSmallValues(double[] tempForces, double forceAndTorqueConvergence) {
            for (int d = 0; d < spatialDim; d++) {
                hydrodynamicForces[0][d] = 0;
                if (Math.Abs(tempForces[d]) > forceAndTorqueConvergence * 1e-2)
                    hydrodynamicForces[0][d] = tempForces[d];
            }
        }

        private void TorqueClearSmallValues(double tempTorque, double forceAndTorqueConvergence) {
            hydrodynamicTorque[0] = 0;
            if (Math.Abs(tempTorque) > forceAndTorqueConvergence * 1e-2)
                hydrodynamicTorque[0] = tempTorque;
        }

        /// <summary>
        /// This method performs the Neumaier algorithm form the sum of the entries of an array.
        /// </summary>
        /// <param name="ResultVariable">
        /// The variable where  the sum will be saved.
        /// </param>
        /// <param name="Summands">
        /// The array of summands
        /// </param>
        /// <param name="Length">
        /// The number of summands.
        /// </param>
        private double ForceTorqueSummationWithNeumaierArray(double ResultVariable, MultidimensionalArray Summands, double Length) {
            double sum = ResultVariable;
            double naiveSum;
            double c = 0.0;
            for (int i = 0; i < Length; i++) {
                naiveSum = sum + Summands[i, 0];
                if (Math.Abs(sum) >= Math.Abs(Summands[i, 0])) {
                    c += (sum - naiveSum) + Summands[i, 0];
                }
                else {
                    c += (Summands[i, 0] - naiveSum) + sum;
                }
                sum = naiveSum;
            }
            return sum + c;
        }
    }
}

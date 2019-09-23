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

        public Motion_Wet(double[] gravity,
            double density,
            ParticleUnderrelaxationParam underrelaxationParam = null) {
            Gravity = gravity;
            m_UnderrelaxationParam = underrelaxationParam;
            Density = density;

            for (int i = 0; i < historyLength; i++) {
                position.Add(new double[spatialDim]);
                Angle.Add(new double());
                TranslationalVelocity.Add(new double[spatialDim]);
                TranslationalAcceleration.Add(new double[spatialDim]);
                RotationalVelocity.Add(new double());
                RotationalAcceleration.Add(new double());
                HydrodynamicForces.Add(new double[spatialDim]);
                HydrodynamicTorque.Add(new double());
            }
        }

        private const int historyLength = 4;
        protected static int spatialDim = 2;
        [NonSerialized]
        readonly internal FSI_Auxillary Aux = new FSI_Auxillary();
        [NonSerialized]
        readonly private ParticleForceIntegration ForceIntegration = new ParticleForceIntegration();
        readonly ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation();
        private readonly ParticleUnderrelaxationParam m_UnderrelaxationParam = null;

        /// <summary>
        /// Gravity (volume force) acting on the particle.
        /// </summary>
        protected double[] Gravity { get; private set; }

        /// <summary>
        /// Density of the particle.
        /// </summary>
        public double Density { get; private set; } = 1;

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        protected double ParticleArea { get; private set; }

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        protected double MomentOfInertia { get; private set; } 

        /// <summary>
        /// The maximum lenghtscale of the particle.
        /// </summary>
        protected double MaxParticleLengthScale { get; private set; }
        
        /// <summary>
        /// Include rotation?
        /// </summary>
        public virtual bool IncludeRotation { get; } = true;

        /// <summary>
        /// Include translation?
        /// </summary>
        public virtual bool IncludeTranslation { get; } = true;

        /// <summary>
        /// Use added damping?, for reference: Banks et.al. 2017
        /// </summary>
        public virtual bool UseAddedDamping { get; } = false;

        /// <summary>
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        public virtual double[,] AddedDampingTensor { get; } = new double[6, 6];

        /// <summary>
        /// The position of the particle in the current time step.
        /// </summary>
        public List<double[]> position = new List<double[]>();
        
        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> Angle { get; private set; } = new List<double>();

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        public List<double[]> TranslationalVelocity { get; private set; } = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        public List<double> RotationalVelocity { get; private set; } = new List<double>(); 

        /// <summary>
        /// The translational velocity of the particle in the current time step.
        /// </summary>
        protected List<double[]> TranslationalAcceleration { get; private set; } = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step.
        /// </summary>
        protected List<double> RotationalAcceleration { get; private set; } = new List<double>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public List<double[]> HydrodynamicForces { get; private set; } = new List<double[]>();

        /// <summary>
        /// The Torque acting on the particle in the current time step.
        /// </summary>
        public List<double> HydrodynamicTorque { get; private set; } = new List<double>();

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public double[] ForcesPrevIteration { get; private set; }

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public double TorquePrevIteration { get; private set; }

        /// <summary>
        /// The force acting on the particle in the current time step.
        /// </summary>
        public double CollisionTimestep { get; set; }

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        public double[] PreCollisionVelocity { get; private set; }

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double[]> CollisionTranslationalVelocity = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double[]> collisionNormalVector = new List<double[]>();

        /// <summary>
        /// The translational velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double[]> collisionTangentialVector = new List<double[]>();

        /// <summary>
        /// The angular velocity of the particle in the current time step. This list is used by the momentum conservation model.
        /// </summary>
        [DataMember]
        public List<double> CollisionRotationalVelocity = new List<double>();

        /// <summary>
        /// Saves position and angle of the last timestep
        /// </summary>
        public void SavePositionAndAngleOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(position);
            Aux.SaveValueOfLastTimestep(Angle);
        }

        /// <summary>
        /// Saves translational and rotational velocities of the last timestep
        /// </summary>
        public void SaveVelocityOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(TranslationalVelocity);
            Aux.SaveValueOfLastTimestep(RotationalVelocity);
            Aux.SaveMultidimValueOfLastTimestep(TranslationalAcceleration);
            Aux.SaveValueOfLastTimestep(RotationalAcceleration);
        }

        /// <summary>
        /// Saves force and torque of the previous timestep
        /// </summary>
        public void SaveHydrodynamicsOfPreviousIteration() {
            ForcesPrevIteration = HydrodynamicForces[0].CloneAs();
            TorquePrevIteration = HydrodynamicTorque[0];
        }

        /// <summary>
        /// Saves force and torque of the previous timestep
        /// </summary>
        public void SaveHydrodynamicsOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(HydrodynamicForces);
            Aux.SaveValueOfLastTimestep(HydrodynamicTorque);
        }

        public void InitializeParticlePositionAndAngle(double[] positionP, double angleP) {
            for (int i = 0; i < historyLength; i++) {
                position[i] = positionP.CloneAs();
                Angle[i] = angleP * 2 * Math.PI / 360;
            }
        }

        public void InitializeParticleVelocity(double[] translationalVelocityP, double rotationalVelocityP) {
            if (translationalVelocityP == null)
                TranslationalVelocity[0] = new double[spatialDim];
            else
                TranslationalVelocity[0] = translationalVelocityP.CloneAs();

            RotationalVelocity[0] = rotationalVelocityP;
        }

        /// <summary>
        /// Init of the particle mass.
        /// </summary>
        /// <param name="mass"></param>
        public void GetParticleArea(double area) => ParticleArea = area;

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        public double Mass_P {
            get {
                Aux.TestArithmeticException(ParticleArea, "particle area");
                Aux.TestArithmeticException(Density, "particle density");
                return ParticleArea * Density;
            }
        }

        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void GetParticleLengthscale(double lengthscale) => MaxParticleLengthScale = lengthscale;

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        public double ParticleMass {
            get {
                Aux.TestArithmeticException(ParticleArea, "particle area");
                Aux.TestArithmeticException(Density, "particle density");
                return ParticleArea * Density;
            }
        }
        
        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void GetParticleMomentOfInertia(double moment) {
            MomentOfInertia = moment;
        }

        /// <summary>
        /// Calls the calculation of the position and angle.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticlePositionAndAngle(double dt) {
            if (CollisionTimestep == 0) {
                SavePositionAndAngleOfPreviousTimestep();
                position[0] = CalculateParticlePosition(dt);
                Angle[0] = CalculateParticleAngle(dt);
            }
            else {
                if (CollisionTimestep < 0)
                    CollisionTimestep = 0;
                double[] tempPos = position[0].CloneAs();
                double tempAngle = Angle[0];
                SavePositionAndAngleOfPreviousTimestep();
                position[0] = tempPos.CloneAs();
                Angle[0] = tempAngle;
                position[0] = CalculateParticlePosition(dt, CollisionTimestep);
                Angle[0] = CalculateParticleAngle(dt, CollisionTimestep);
                if (CollisionTimestep > dt) { CollisionTimestep -= dt; }
                else CollisionTimestep = 0;
            }
        }

        /// <summary>
        /// Calls the calculation of the position and angle.
        /// </summary>
        /// <param name="dt"></param>
        public void CollisionParticlePositionAndAngle(double collisionDynamicTimestep) {
            position[0] = CalculateParticlePosition(collisionDynamicTimestep, collisionProcedure: true);
            Angle[0] = CalculateParticleAngle(collisionDynamicTimestep, collisionProcedure: true);
        }

        /// <summary>
        /// Calls the calculation of the velocity.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticleVelocity(double dt) {
            TranslationalAcceleration[0] = CalculateTranslationalAcceleration(dt - CollisionTimestep);
            RotationalAcceleration[0] = CalculateRotationalAcceleration(dt - CollisionTimestep);
            if (CollisionTimestep == 0) {
                CalculateTranslationalVelocity(dt);
                CalculateAngularVelocity(dt);
            }
            else {
                CalculateTranslationalVelocity(dt, CollisionTimestep);
                CalculateAngularVelocity(dt, CollisionTimestep);
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
                HydrodynamicForces[0][0] = MaxParticleLengthScale * Math.Cos(Angle[0]) * activeStress / 2 + Gravity[0] * Density * ParticleArea / 10;
                HydrodynamicForces[0][1] = MaxParticleLengthScale * Math.Sin(Angle[0]) * activeStress / 2 + Gravity[1] * Density * ParticleArea / 10;
                HydrodynamicTorque[0] = 0;
            }
            else {
                for (int d = 0; d < spatialDim; d++) {
                    HydrodynamicForces[0][d] = (HydrodynamicForces[1][d] + 4 * HydrodynamicForces[2][d] + HydrodynamicForces[3][d]) / 6;
                    if (Math.Abs(HydrodynamicForces[0][d]) < 1e-20)
                        HydrodynamicForces[0][d] = 0;
                }
                HydrodynamicTorque[0] = (HydrodynamicTorque[1] + 4 * HydrodynamicTorque[2] + HydrodynamicTorque[3]) / 6;
                if (Math.Abs(HydrodynamicTorque[0]) < 1e-20)
                    HydrodynamicTorque[0] = 0;
            }
            Aux.TestArithmeticException(HydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(HydrodynamicTorque[0], "hydrodynamic torque");
        }

        public virtual void CalculateDampingTensor(Particle particle, LevelSetTracker LsTrk, double muA, double rhoA, double dt) {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        public virtual void UpdateDampingTensors() {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        public virtual void ExchangeAddedDampingTensors() {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double[] CalculateParticlePosition(double dt) {
            double[] l_Position = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_Position[d] = position[1][d] + (TranslationalVelocity[0][d] + 4 * TranslationalVelocity[1][d] + TranslationalVelocity[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double[] CalculateParticlePosition(double dt, bool collisionProcedure) {
            double[] l_Position = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_Position[d] = position[0][d] + (TranslationalVelocity[0][d] + 4 * TranslationalVelocity[1][d] + TranslationalVelocity[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual double[] CalculateParticlePosition(double dt, double collisionTimestep) {
            double[] l_Position = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_Position[d] = position[0][d] + TranslationalVelocity[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateParticleAngle(double dt) {
            double l_Angle = Angle[1] + (RotationalVelocity[0] + 4 * RotationalVelocity[1] + RotationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(Angle[0], "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateParticleAngle(double dt, bool collisionProcedure) {
            double l_Angle = Angle[0] + (RotationalVelocity[0] + 4 * RotationalVelocity[1] + RotationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(Angle[0], "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new particle angle after a collision
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual double CalculateParticleAngle(double dt, double collisionTimestep) {
            double l_Angle = Angle[1] + RotationalVelocity[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(Angle[0], "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected virtual void CalculateTranslationalVelocity(double dt) {
            for (int d = 0; d < spatialDim; d++) {
                TranslationalVelocity[0][d] = TranslationalVelocity[1][d] + (TranslationalAcceleration[0][d] + 4 * TranslationalAcceleration[1][d] + TranslationalAcceleration[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(TranslationalVelocity[0], "particle translational velocity");
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle using a Crank Nicolson scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual void CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            for (int d = 0; d < spatialDim; d++) {
                TranslationalVelocity[0][d] = TranslationalVelocity[1][d] + TranslationalAcceleration[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(TranslationalVelocity[0], "particle translational velocity");
        }

        public virtual void WritePostCollisionVelocity(double[] postCollisionVelocity) {
            TranslationalVelocity[0] = postCollisionVelocity;
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected virtual void CalculateAngularVelocity(double dt) {
            RotationalVelocity[0] = RotationalVelocity[1] + (RotationalAcceleration[0] + 4 * RotationalAcceleration[1] + RotationalAcceleration[2]) * dt / 6;
            Aux.TestArithmeticException(RotationalVelocity[0], "particle rotational velocity");
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle using explicit Euler scheme.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual void CalculateAngularVelocity(double dt, double collisionTimestep) {
            RotationalVelocity[0] = RotationalVelocity[1] + RotationalAcceleration[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(RotationalVelocity[0], "particle rotational velocity");
        }

        public virtual void WritePostCollisionAngularVelocity(double postCollisionAngularVelocity) {
            RotationalVelocity[0] = postCollisionAngularVelocity;
        }

        internal void CalculateNormalAndTangentialVelocity(double[] NormalVector) {
            double[] Velocity = TranslationalVelocity[0];
            double[] TangentialVector = new double[] { -NormalVector[1], NormalVector[0] };
            PreCollisionVelocity = new double[] { Velocity[0] * NormalVector[0] + Velocity[1] * NormalVector[1], Velocity[0] * TangentialVector[0] + Velocity[1] * TangentialVector[1] };
            Aux.TestArithmeticException(PreCollisionVelocity, "particle velocity before collision");
        }

        /// <summary>
        /// Calculate the new acceleration
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double[] CalculateTranslationalAcceleration(double dt) {
            double[] l_Acceleration = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_Acceleration[d] = HydrodynamicForces[0][d] / (Density * ParticleArea);
            }
            Aux.TestArithmeticException(l_Acceleration, "particle translational acceleration");
            return l_Acceleration;
        }

        /// <summary>
        /// Calculate the new acceleration (translational and rotational)
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateRotationalAcceleration(double dt) {
            double l_Acceleration = HydrodynamicTorque[0] / MomentOfInertia;
            Aux.TestArithmeticException(l_Acceleration, "particle rotational acceleration");
            return l_Acceleration;
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
            SinglePhaseField[] UA = U.ToArray();
            ConventionalDGField pA = P;
            double[] tempForces = ForcesIntegration(UA, pA, LsTrk, CutCells_P, RequiredOrder, muA);
            Force_MPISum(ref tempForces);
            for (int d = 0; d < spatialDim; d++) {
                tempForces[d] += (Density - fluidDensity) * ParticleArea * Gravity[d];
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
            if (m_UnderrelaxationParam != null && !firstIteration) {
                double forceAndTorqueConvergence = m_UnderrelaxationParam.hydroDynConvergenceLimit;
                double underrelaxationFactor = m_UnderrelaxationParam.underrelaxationFactor;
                bool useAddaptiveUnderrelaxation = m_UnderrelaxationParam.useAddaptiveUnderrelaxation;
                Underrelaxation.CalculateAverageForces(tempForces, tempTorque, MaxParticleLengthScale, out double averagedForces);
                Underrelaxation.Forces(ref tempForces, ForcesPrevIteration, forceAndTorqueConvergence, underrelaxationFactor, useAddaptiveUnderrelaxation, averagedForces);
                Underrelaxation.Torque(ref tempTorque, TorquePrevIteration, forceAndTorqueConvergence, underrelaxationFactor, useAddaptiveUnderrelaxation, averagedForces);
            }
            for (int d = 0; d < spatialDim; d++) {
                HydrodynamicForces[0][d] = 0;
                HydrodynamicForces[0][d] = tempForces[d];
            }
            HydrodynamicTorque[0] = tempTorque;
            Aux.TestArithmeticException(HydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(HydrodynamicTorque[0], "hydrodynamic torque");
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

        /// <summary>
        /// Calculating the particle Reynolds number
        /// </summary>
        public double ComputeParticleRe(double ViscosityFluid) {
            return TranslationalVelocity[0].L2Norm() * MaxParticleLengthScale / ViscosityFluid;
        }

        /// <summary>
        /// Calculating the particle Stokes number
        /// </summary>
        public double ComputeParticleSt(double ViscosityFluid, double DensityFluid) {
            return ComputeParticleRe(ViscosityFluid) * Density / (9 * DensityFluid);
        }

        public double[] CalculateParticleMomentum() {
            double[] temp = new double[spatialDim + 1];
            for (int d = 0; d < spatialDim; d++) {
                temp[d] = Mass_P * TranslationalVelocity[0][d];
            }
            temp[spatialDim] = MomentOfInertia * RotationalVelocity[0];
            return temp;
        }

        public double[] CalculateParticleKineticEnergy() {
            double[] temp = new double[spatialDim + 1];
            for (int d = 0; d < spatialDim; d++) {
                temp[d] = 0.5 * Mass_P * TranslationalVelocity[0][d].Pow2();
            }
            temp[spatialDim] = 0.5 * MomentOfInertia * RotationalVelocity[0].Pow2();
            return temp;
        }

        public void ClearParticleHistoryTranslation() {
            for (int i = 0; i < TranslationalVelocity.Count; i++) {
                for (int d = 0; d < spatialDim; d++) {
                    TranslationalVelocity[i][d] = 0;
                    TranslationalAcceleration[i][d] = 0;
                }
            }
        }

        public void ClearParticleHistoryRotational() {
            for (int i = 0; i < RotationalVelocity.Count; i++) {
                RotationalAcceleration[i] = 0;
                RotationalVelocity[i] = 0;
            }
        }

        public double[] BuildSendArray() {
            double[] dataSend = new double[13];
            dataSend[0] = RotationalVelocity[0];
            dataSend[1] = TranslationalVelocity[0][0];
            dataSend[2] = TranslationalVelocity[0][1];
            dataSend[3] = Angle[0];
            dataSend[4] = position[0][0];
            dataSend[5] = position[0][1];
            dataSend[6] = CollisionTimestep;
            dataSend[7] = RotationalVelocity[1];
            dataSend[8] = TranslationalVelocity[1][0];
            dataSend[9] = TranslationalVelocity[1][1];
            dataSend[10] = Angle[1];
            dataSend[11] = position[1][0];
            dataSend[12] = position[1][1];
            return dataSend;
        }

        public void WriteReceiveArray(double[] dataReceive, int Offset) {
            RotationalVelocity[0] = dataReceive[0 + Offset];
            TranslationalVelocity[0][0] = dataReceive[1 + Offset];
            TranslationalVelocity[0][1] = dataReceive[2 + Offset];
            Angle[0] = dataReceive[3 + Offset];
            position[0][0] = dataReceive[4 + Offset];
            position[0][1] = dataReceive[5 + Offset];
            CollisionTimestep = dataReceive[6 + Offset];
            RotationalVelocity[1] = dataReceive[7 + Offset];
            TranslationalVelocity[1][0] = dataReceive[8 + Offset];
            TranslationalVelocity[1][1] = dataReceive[9 + Offset];
            Angle[1] = dataReceive[10 + Offset];
            position[1][0] = dataReceive[11 + Offset];
            position[1][1] = dataReceive[12 + Offset];
        }
    }
}

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
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using FSI_Solver;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    [Serializable]
    public class Motion : ICloneable {

        /// <summary>
        /// The standard description of motion including hydrodynamics.
        /// </summary>
        /// <param name="gravity">
        /// The gravity (volume forces) acting on the particle.
        /// </param>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        public Motion(Vector gravity, double density) {
            if (gravity.IsNullOrEmpty())
                gravity = new Vector(0, 0);
            Gravity = new Vector(gravity);
            Density = density;
            for (int i = 0; i < NumberOfHistoryEntries; i++) {
                Position.Add(new Vector(SpatialDim));
                TranslationalVelocity.Add(new Vector(SpatialDim));
                TranslationalAcceleration.Add(new Vector(SpatialDim));
                HydrodynamicForces.Add(new Vector(SpatialDim));
                Angle.Add(new double());
                RotationalVelocity.Add(new double());
                RotationalAcceleration.Add(new double());
                HydrodynamicTorque.Add(new double());
            }
        }
        internal double omega = 1;
        [NonSerialized]
        internal FSI_Auxillary Aux = new FSI_Auxillary();
        [DataMember]
        private const int NumberOfHistoryEntries = 3;
        [DataMember]
        protected static int SpatialDim = 2;
        [DataMember]
        private readonly List<Vector> Position = new List<Vector>();
        [DataMember]
        private readonly List<double> Angle = new List<double>();
        [DataMember]
        private readonly List<Vector> TranslationalVelocity = new List<Vector>();
        [DataMember]
        private readonly List<double> RotationalVelocity = new List<double>();
        [DataMember]
        private readonly List<Vector> TranslationalAcceleration = new List<Vector>();
        [DataMember]
        private readonly List<double> RotationalAcceleration = new List<double>();
        [DataMember]
        private readonly List<Vector> HydrodynamicForces = new List<Vector>();
        [DataMember]
        private readonly List<double> HydrodynamicTorque = new List<double>();
        [DataMember]
        private double CollisionTimestep = 0;
        [DataMember]
        public readonly double Density;
        [DataMember]
        public double ParticleArea;
        [DataMember]
        public Vector Gravity = new Vector(0 , 0);
        [DataMember]
        public double MomentOfInertia;
        [DataMember]
        public double MaxParticleLengthScale;
        [DataMember]
        internal double[,] AddedDampingTensor = new double[6, 6];

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        [DataMember]
        public double ParticleMass => ParticleArea * Density;
        
        /// <summary>
        /// Include rotation?
        /// </summary>
        [DataMember]
        internal virtual bool IncludeRotation { get; } = true;

        /// <summary>
        /// Include translation?
        /// </summary>
        [DataMember]
        internal virtual bool IncludeTranslation { get; } = true;

        /// <summary>
        /// Use added damping?, for reference: Banks et.al. 2017
        /// </summary>
        [DataMember]
        internal virtual bool UseAddedDamping { get; } = false;

        /// <summary>
        /// Returns the position of the particle.
        /// </summary>
        /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetPosition(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return Position[historyPosition];
        }

        /// <summary>
        /// Returns the angle of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetAngle(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return Angle[historyPosition];
        }

        /// <summary>
        /// Returns the translational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetTranslationalVelocity(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return TranslationalVelocity[historyPosition];
        }

        /// <summary>
        /// Returns the rotational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetRotationalVelocity(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return RotationalVelocity[historyPosition];
        }

        /// <summary>
        /// Returns the translational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetTranslationalAcceleration(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return TranslationalAcceleration[historyPosition];
        }

        /// <summary>
        /// Returns the rotational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetRotationalAcceleration(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return RotationalAcceleration[historyPosition];
        }

        /// <summary>
        /// Returns the force acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for 4 timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetHydrodynamicForces(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return HydrodynamicForces[historyPosition];
        }

        /// <summary>
        /// Returns the torque acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetHydrodynamicTorque(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return HydrodynamicTorque[historyPosition];
        }

        /// <summary>
        /// Saves position and angle of the last timestep.
        /// </summary>
        public void SavePositionAndAngleOfPreviousTimestep() {
            using (new FuncTrace()) {
                Aux.SaveVectorOfLastTimestep(Position);
                Aux.SaveValueOfLastTimestep(Angle);
            }
        }

        /// <summary>
        /// Saves translational and rotational velocities of the last timestep.
        /// </summary>
        public void SaveVelocityOfPreviousTimestep() {
            using (new FuncTrace()) {
                Aux.SaveVectorOfLastTimestep(TranslationalVelocity);
                Aux.SaveValueOfLastTimestep(RotationalVelocity);
                Aux.SaveVectorOfLastTimestep(TranslationalAcceleration);
                Aux.SaveValueOfLastTimestep(RotationalAcceleration);
            }
        }

        /// <summary>
        /// Saves force and torque of the previous timestep.
        /// </summary>
        public void SaveHydrodynamicsOfPreviousTimestep() {
            using (new FuncTrace()) {
                Aux.SaveVectorOfLastTimestep(HydrodynamicForces);
                for (int i = 0; i < NumberOfHistoryEntries; i++)
                    Aux.SaveValueOfLastTimestep(HydrodynamicTorque);
            }
        }

        /// <summary>
        /// Used during init of the particle. Sets the position and the angle.
        /// </summary>
        /// <param name="initialPosition">
        /// The initial position.
        /// </param>
        /// <param name="initialAngle">
        /// The initial angle.
        /// </param>
        internal void InitializeParticlePositionAndAngle(double[] initialPosition, double initialAngle, int historyLength = 0) {
            using (new FuncTrace()) {
                if (historyLength == 0)
                    historyLength = NumberOfHistoryEntries;
                Aux = new FSI_Auxillary();
                for (int i = 0; i < historyLength; i++) {
                    Position[i] = new Vector(initialPosition);
                    Angle[i] = initialAngle * 2 * Math.PI / 360;
                    Aux.TestArithmeticException(Position[i], "initial particle position");
                    Aux.TestArithmeticException(Angle[i], "initial particle angle");
                }
                orientationVector = new Vector(Math.Cos(Angle[0]), Math.Sin(Angle[0]));
            }
        }

        /// <summary>
        /// Used during init of the particle. Sets the translational and rotational velocity.
        /// </summary>
        /// <param name="initalTranslation">
        /// The initial translational velocity.
        /// </param>
        /// <param name="initalRotation">
        /// The initial rotational velocity.
        /// </param>
        internal void InitializeParticleVelocity(double[] initalTranslation, double initalRotation, int historyLength = 0) {
            using (new FuncTrace()) {
                if (historyLength == 0)
                    historyLength = NumberOfHistoryEntries;
                for (int i = 0; i < historyLength; i++) {
                    TranslationalVelocity[i] = initalTranslation == null ? new Vector(SpatialDim) : new Vector(initalTranslation);
                    RotationalVelocity[i] = initalRotation;
                    Aux.TestArithmeticException(TranslationalVelocity[i], "initial particle translational velocity");
                    Aux.TestArithmeticException(RotationalVelocity[i], "initial particle rotational velocity");
                }
            }
        }

        /// <summary>
        /// Used during init of the particle. Sets the translational and rotational velocity.
        /// </summary>
        /// <param name="initalTranslation">
        /// The initial translational velocity.
        /// </param>
        /// <param name="initalRotation">
        /// The initial rotational velocity.
        /// </param>
        internal void InitializeParticleAcceleration(double[] initalTranslationAcceleration, double initalRotationAcceleration, int historyLength = 0) {
            using (new FuncTrace()) {
                if (historyLength == 0)
                    historyLength = NumberOfHistoryEntries;
                for (int i = 0; i < historyLength; i++) {
                    TranslationalAcceleration[i] = new Vector(initalTranslationAcceleration);
                    RotationalAcceleration[i] = initalRotationAcceleration;
                    Aux.TestArithmeticException(TranslationalVelocity[i], "initial particle translational velocity");
                    Aux.TestArithmeticException(RotationalVelocity[i], "initial particle rotational velocity");
                }
            }
        }

        /// <summary>
        /// Used during init of the particle. Sets the translational and rotational velocity.
        /// </summary>
        /// <param name="initalTranslation">
        /// The initial translational velocity.
        /// </param>
        /// <param name="initalRotation">
        /// The initial rotational velocity.
        /// </param>
        internal void InitializeParticleForceAndTorque(double[] initalForce, double initalTorque, int historyLength = 0, double dt = 0) {
            using (new FuncTrace()) {
                if (historyLength == 0)
                    historyLength = NumberOfHistoryEntries;
                for (int i = 0; i < historyLength; i++) {
                    HydrodynamicForces[i] = new Vector(initalForce);
                    HydrodynamicTorque[i] = initalTorque;
                    if (dt != 0) {
                        TranslationalAcceleration[i] = CalculateTranslationalAcceleration(dt);
                        RotationalAcceleration[i] = CalculateRotationalAcceleration(dt);
                    }
                    Aux.TestArithmeticException(TranslationalVelocity[i], "initial particle translational velocity");
                    Aux.TestArithmeticException(RotationalVelocity[i], "initial particle rotational velocity");
                }
            }
        }

        /// <summary>
        /// Transfer data form motionDataToTransfer to the current particle.Motion
        /// </summary>
        /// <param name="motionDataToTransfer">
        /// The initial translational velocity.
        /// </param>
        internal void TransferDataFromOtherParticle(Motion motionDataToTransfer) {
            using (new FuncTrace()) {
                for (int h = 0; h < NumberOfHistoryEntries; h++) {
                    Position[h] = new Vector(motionDataToTransfer.GetPosition(h));
                    TranslationalVelocity[h] = new Vector(motionDataToTransfer.GetTranslationalVelocity(h));
                    TranslationalAcceleration[h] = new Vector(motionDataToTransfer.GetTranslationalAcceleration(h));
                    HydrodynamicForces[h] = new Vector(motionDataToTransfer.GetHydrodynamicForces(h));
                    Angle[h] = motionDataToTransfer.GetAngle(h);
                    RotationalVelocity[h] = motionDataToTransfer.GetRotationalVelocity(h);
                    RotationalAcceleration[h] = motionDataToTransfer.GetRotationalAcceleration(h);
                    HydrodynamicTorque[h] = motionDataToTransfer.GetHydrodynamicTorque(h);
                }
            }
        }

        /// <summary>
        /// Adapt particle history to new time-step size. Necessary for adaptive time-stepping
        /// </summary>
        /// <param name="newTimestep">
        /// New time-step size
        /// </param>
        /// <param name="oldTimestep">
        /// Old time-step size
        /// </param>
        internal void AdaptParticleHistoryToNewTimeStep(double newTimestep, double oldTimestep) {
            using (new FuncTrace()) {
                AdaptParticleHistoryToNewTimeStep(Position, newTimestep, oldTimestep);
                AdaptParticleHistoryToNewTimeStep(TranslationalVelocity, newTimestep, oldTimestep);
                AdaptParticleHistoryToNewTimeStep(TranslationalAcceleration, newTimestep, oldTimestep);
                AdaptParticleHistoryToNewTimeStep(Angle, newTimestep, oldTimestep);
                AdaptParticleHistoryToNewTimeStep(RotationalVelocity, newTimestep, oldTimestep);
                AdaptParticleHistoryToNewTimeStep(RotationalAcceleration, newTimestep, oldTimestep);
            }
        }

        private void AdaptParticleHistoryToNewTimeStep(List<Vector> variable, double newTimestep, double oldTimestep) {
            using (new FuncTrace()) {
                if (oldTimestep == 0)// for restart
                    oldTimestep = 1e-6;
                List<Vector> stuetzstelle = new List<Vector>();
                for (int h = 0; h < NumberOfHistoryEntries; h++) {
                    stuetzstelle.Add(new Vector(variable[h]));
                }
                List<Vector> newVariable = new List<Vector> { new Vector(variable[0]) };
                for (int i = 1; i < NumberOfHistoryEntries; i++) {
                    double currentNewTimestep = -i * newTimestep;
                    double[] langrangePoly = CalculateLangrangePolynom(currentNewTimestep, oldTimestep);
                    newVariable.Add(new Vector(langrangePoly[0] * variable[0]));
                    for (int j = 1; j < NumberOfHistoryEntries; j++) {
                        newVariable[i] = newVariable[i] + langrangePoly[j] * stuetzstelle[j];
                    }
                }
                variable.Clear();
                for (int h = 0; h < NumberOfHistoryEntries; h++) {
                    variable.Add(new Vector(newVariable[h]));
                }
            }
        }

        private void AdaptParticleHistoryToNewTimeStep(List<double> variable, double newTimestep, double oldTimestep) {
            using (new FuncTrace()) {
                if (oldTimestep == 0)// for restart
                    oldTimestep = 1e-6;
                List<double> stuetzstelle = new List<double>();
                for (int h = 0; h < NumberOfHistoryEntries; h++) {
                    stuetzstelle.Add(variable[h]);
                }
                List<double> newVariable = new List<double>();
                newVariable.Add(variable[0]);
                for (int i = 1; i < NumberOfHistoryEntries; i++) {
                    double currentNewTimestep = -i * newTimestep;
                    double[] langrangePoly = CalculateLangrangePolynom(currentNewTimestep, oldTimestep);
                    newVariable.Add(langrangePoly[0] * variable[i]);
                    for (int j = 1; j < NumberOfHistoryEntries; j++) {
                        newVariable[i] = newVariable[i] + langrangePoly[j] * stuetzstelle[j];
                    }
                }
                variable.Clear();
                for (int h = 0; h < NumberOfHistoryEntries; h++) {
                    variable.Add(newVariable[h]);
                }
            }
        }

        private double[] CalculateLangrangePolynom(double time, double oldTimestep) {
            using (new FuncTrace()) {
                double[] lPoly = new double[NumberOfHistoryEntries];
                for (int i = 0; i < NumberOfHistoryEntries; i++) {
                    if (i != 0)
                        lPoly[i] = time / (-i * oldTimestep);
                    for (int j = 1; j < NumberOfHistoryEntries; j++) {
                        if (j == i)
                            continue;
                        if (i == 0 && j == 1)
                            lPoly[i] = (time + j * oldTimestep) / ((j - i) * oldTimestep);
                        lPoly[i] *= (time + j * oldTimestep) / ((j - i) * oldTimestep);
                    }
                }
                return lPoly;
            }
        }

        /// <summary>
        /// Init of the particle area.
        /// </summary>
        /// <param name="area"></param>
        public void SetParticleArea(double area) => ParticleArea = area;

        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void SetParticleMaxLengthscale(double lengthscale) => MaxParticleLengthScale = lengthscale;

        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void SetParticleMomentOfInertia(double moment) => MomentOfInertia = moment;

        /// <summary>
        /// Sets the collision timestep.
        /// </summary>
        /// <param name="collisionTimestep">
        /// The physical time consumend by the collision procedure.
        /// </param>
        public void SetCollisionTimestep(double collisionTimestep) => CollisionTimestep = collisionTimestep;

        public Vector orientationVector = new Vector(2);

        /// <summary>
        /// Calls the calculation of the position and angle.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticlePositionAndAngle(double dt) {
            using (new FuncTrace()) {
                if (CollisionTimestep < 0)
                    CollisionTimestep = 0;
                SavePositionAndAngleOfPreviousTimestep();
                if (CollisionTimestep > dt) {
                    throw new Exception("Collision timestep: " + CollisionTimestep);
                }
                Position[0] = CalculateParticlePosition(dt - CollisionTimestep);
                Angle[0] = CalculateParticleAngle(dt - CollisionTimestep);
                orientationVector = new Vector(Math.Cos(Angle[0]), Math.Sin(Angle[0]));
                CollisionTimestep = 0;
            }
        }

        /// <summary>
        /// Calls the calculation of the position and angle during the calculation of the collisions.
        /// </summary>
        /// <param name="dt"></param>
        internal void CollisionParticlePositionAndAngle(double collisionDynamicTimestep) {
            using (new FuncTrace()) {
                Position[0] = CalculateParticlePositionDuringCollision(collisionDynamicTimestep);
                Angle[0] = CalculateParticleAngleDuringCollision(collisionDynamicTimestep);
            }
        }

        /// <summary>
        /// Calls the calculation of the velocity.
        /// </summary>
        /// <param name="dt"></param>
        internal void UpdateParticleVelocity(double dt) {
            using (new FuncTrace()) {
                TranslationalAcceleration[0] = CalculateTranslationalAcceleration(dt - CollisionTimestep);
                RotationalAcceleration[0] = CalculateRotationalAcceleration(dt - CollisionTimestep);
                TranslationalVelocity[0] = CalculateTranslationalVelocity(dt - CollisionTimestep);
                RotationalVelocity[0] = CalculateAngularVelocity(dt - CollisionTimestep);
            }
        }

        internal virtual void SetGhostPosition(Vector position) {
            using (new FuncTrace()) {
                for (int h = 0; h < NumberOfHistoryEntries; h++) {
                    Position[h] = new Vector(position);
                }
            }
        }

        internal virtual void CopyNewVelocity(Vector translational, double rotational) {
            using (new FuncTrace()) {
                TranslationalVelocity[0] = new Vector(translational);
                RotationalVelocity[0] = rotational;
            }
        }

        internal void UpdateForcesAndTorque(int particleID, double[] fullListHydrodynamics) {
            using (new FuncTrace()) {
                Vector forces = new Vector(SpatialDim);
                for (int d = 0; d < SpatialDim; d++) {
                    if (Math.Abs(fullListHydrodynamics[particleID * 3 + d]) > 1e-12)
                        forces[d] = fullListHydrodynamics[particleID * 3 + d];

                }
                HydrodynamicForces[0] = forces;
                if (Math.Abs(fullListHydrodynamics[particleID * 3 + SpatialDim]) > 1e-12)
                    HydrodynamicTorque[0] = fullListHydrodynamics[particleID * 3 + SpatialDim];
                Aux.TestArithmeticException(HydrodynamicForces[0], "hydrodynamic forces");
                Aux.TestArithmeticException(HydrodynamicTorque[0], "hydrodynamic torque");
            }
        }

        /// <summary>
        /// Predicts the hydrodynamics at the beginning of the iteration loop in each timestep.
        /// </summary>
        /// <param name="activeStress"></param>
        /// <param name="timestepID">
        /// The timestep ID. Used to distinguish between the first timestep and all other steps.
        /// </param>
        internal virtual void PredictForceAndTorque(double activeStress, double circumference, int timestepID, double fluidViscosity, double fluidDensity, double dt) {
            using (new FuncTrace()) {
                HydrodynamicForces[0] = new Vector(HydrodynamicForces[1]);
                HydrodynamicTorque[0] = HydrodynamicTorque[1];
                Aux.TestArithmeticException(HydrodynamicForces[0], "hydrodynamic forces");
                Aux.TestArithmeticException(HydrodynamicTorque[0], "hydrodynamic torque");
            }
        }

        /// <summary>
        /// Calculate the tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        /// <param name="particle"></param>
        /// <param name="levelSetTracker"></param>
        /// <param name="fluidViscosity"></param>
        /// <param name="fluidDensity"></param>
        /// <param name="dt"></param>
        internal virtual void CalculateDampingTensor(Particle particle, LevelSetTracker levelSetTracker, double fluidViscosity, double fluidDensity, double dt) {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        /// <summary>
        /// Updates the tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        internal virtual void UpdateDampingTensors() {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        /// <summary>
        /// Exchange of the added damping tensors between the MPI-processes.
        /// </summary>
        internal virtual void ExchangeAddedDampingTensors() {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected virtual Vector CalculateParticlePosition(double dt) {
            using (new FuncTrace()) {
                Vector position = Position[1] + (TranslationalVelocity[0] + 4 * TranslationalVelocity[1] + TranslationalVelocity[2]) * dt / 3;
                Aux.TestArithmeticException(position, "particle position");
                return position;
            }
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateParticleAngle(double dt) {
            using (new FuncTrace()) {
                double angle = Angle[1] + (RotationalVelocity[0] + 4 * RotationalVelocity[1] + RotationalVelocity[2]) * dt / 3;
                Aux.TestArithmeticException(angle, "particle angle");
                return angle;
            }
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected virtual Vector CalculateParticlePositionDuringCollision(double dt) {
            using (new FuncTrace()) {
                Vector position = Position[0] + (TranslationalVelocity[0] + 4 * TranslationalVelocity[1] + TranslationalVelocity[2]) * dt / 3;
                Aux.TestArithmeticException(position, "particle position");
                return position;
            }
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateParticleAngleDuringCollision(double dt) {
            using (new FuncTrace()) {
                double angle = Angle[0] + (RotationalVelocity[0] + 4 * RotationalVelocity[1] + RotationalVelocity[2]) * dt / 3;
                Aux.TestArithmeticException(angle, "particle angle");
                return angle;
            }
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected virtual Vector CalculateTranslationalVelocity(double dt) {
            using (new FuncTrace()) {
                Vector translationalVelocity = TranslationalVelocity[1] + (TranslationalAcceleration[0] + 4 * TranslationalAcceleration[1] + TranslationalAcceleration[2]) * dt / 3;
                Aux.TestArithmeticException(translationalVelocity, "particle translational velocity");
                return translationalVelocity;
            }
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected virtual double CalculateAngularVelocity(double dt) {
            using (new FuncTrace()) {
                double rotationalVelocity = RotationalVelocity[1] + (RotationalAcceleration[0] + 4 * RotationalAcceleration[1] + RotationalAcceleration[2]) * dt / 3;
                Aux.TestArithmeticException(rotationalVelocity, "particle rotational velocity");
                return rotationalVelocity;
            }
        }

        /// <summary>
        /// Calculate the new tranlational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        protected virtual Vector CalculateTranslationalAcceleration(double dt) {
            using (new FuncTrace()) {
                Vector l_Acceleration = HydrodynamicForces[0] / (Density * ParticleArea);
                Aux.TestArithmeticException(l_Acceleration, "particle translational acceleration");
                return l_Acceleration;
            }
        }

        /// <summary>
        /// Calculate the new rotational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateRotationalAcceleration(double dt) {
            using (new FuncTrace()) {
                double l_Acceleration = HydrodynamicTorque[0] / MomentOfInertia;
                Aux.TestArithmeticException(l_Acceleration, "particle rotational acceleration");
                return l_Acceleration;
            }
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        public virtual Vector CalculateHydrodynamicForces(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, double dt = 0) {
            using (new FuncTrace()) {
                Vector tempForces = new Vector(hydrodynamicsIntegration.Forces(cutCells));
                Aux.TestArithmeticException(tempForces, "temporal forces during calculation of hydrodynamics");
                return tempForces;
            }
        }

        /// <summary>
        /// Calculates the gravitational forces.
        /// </summary>
        /// <param name="fluidDensity"></param>
        /// <param name="tempForces"></param>
        public Vector GetGravityForces(double fluidDensity) {
            return (Density - fluidDensity) * ParticleArea * Gravity;
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="cutCells"></param>
        /// <param name="dt"></param>
        public virtual double CalculateHydrodynamicTorque(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, double dt = 0) {
            using (new FuncTrace()) {
                double tempTorque = hydrodynamicsIntegration.Torque(Position[0], cutCells);
                Aux.TestArithmeticException(tempTorque, "temporal torque during calculation of hydrodynamics");
                return tempTorque;
            }
        }

        /// <summary>
        /// Calculating the particle Reynolds number
        /// </summary>
        /// <param name="fluidViscosity"></param>
        public double ComputeParticleReynoldsNumber(double fluidViscosity) {
            using (new FuncTrace()) {
                return 2 * TranslationalVelocity[0].L2Norm() * MaxParticleLengthScale / fluidViscosity;
            }
        }

        /// <summary>
        /// Calculating the particle Stokes number
        /// </summary>
        /// <param name="fluidViscosity"></param>
        /// <param name="fluidDensity"></param>
        public double ComputeParticleStokesNumber(double fluidViscosity, double fluidDensity) {
            using (new FuncTrace()) {
                return ComputeParticleReynoldsNumber(fluidViscosity) * Density / (9 * fluidDensity);
            }
        }

        /// <summary>
        /// Calculating the particle momentum
        /// </summary>
        public double[] CalculateParticleMomentum() {
            using (new FuncTrace()) {
                double[] temp = new double[SpatialDim + 1];
                for (int d = 0; d < SpatialDim; d++) {
                    temp[d] = ParticleMass * TranslationalVelocity[0][d];
                }
                temp[SpatialDim] = MomentOfInertia * RotationalVelocity[0];
                return temp;
            }
        }

        /// <summary>
        /// Calculating the particle kinetic energy
        /// </summary>
        public double[] CalculateParticleKineticEnergy() {
            using (new FuncTrace()) {
                double[] temp = new double[SpatialDim + 1];
                for (int d = 0; d < SpatialDim; d++) {
                    temp[d] = 0.5 * ParticleMass * TranslationalVelocity[0][d].Pow2();
                }
                temp[SpatialDim] = 0.5 * MomentOfInertia * RotationalVelocity[0].Pow2();
                return temp;
            }
        }

        /// <summary>
        /// Builds the array for the post-collision communication between MPI-processes.
        /// </summary>
        public double[] BuildSendArray() {
            using (new FuncTrace()) {
                double[] dataSend = new double[19];
                dataSend[0] = RotationalVelocity[0];
                dataSend[1] = RotationalAcceleration[0];
                dataSend[2] = TranslationalVelocity[0][0];
                dataSend[3] = TranslationalVelocity[0][1];
                dataSend[4] = TranslationalAcceleration[0][0];
                dataSend[5] = TranslationalAcceleration[0][1];
                dataSend[6] = Angle[0];
                dataSend[7] = Position[0][0];
                dataSend[8] = Position[0][1];
                dataSend[9] = CollisionTimestep;
                dataSend[10] = RotationalVelocity[1];
                dataSend[11] = TranslationalVelocity[1][0];
                dataSend[12] = TranslationalVelocity[1][1];
                dataSend[13] = Angle[1];
                dataSend[14] = Position[1][0];
                dataSend[15] = Position[1][1];
                dataSend[16] = RotationalAcceleration[1];
                dataSend[17] = TranslationalAcceleration[1][0];
                dataSend[18] = TranslationalAcceleration[1][1];
                return dataSend;
            }
        }

        /// <summary>
        /// Overwrites the particles parameters with the values received during the post-collision MPI-communication.
        /// </summary>
        public void WriteReceiveArray(double[] dataReceive, int offset) {
            using (new FuncTrace()) {
                RotationalVelocity[0] = dataReceive[0 + offset];
                RotationalAcceleration[0] = dataReceive[1 + offset];
                TranslationalVelocity[0] = new Vector(dataReceive[2 + offset], dataReceive[3 + offset]);
                TranslationalAcceleration[0] = new Vector(dataReceive[4 + offset], dataReceive[5 + offset]);
                Angle[0] = dataReceive[6 + offset];
                Position[0] = new Vector(dataReceive[7 + offset], dataReceive[8 + offset]);
                CollisionTimestep = dataReceive[9 + offset];
                RotationalVelocity[1] = dataReceive[10 + offset];
                TranslationalVelocity[1] = new Vector(dataReceive[11 + offset], dataReceive[12 + offset]);
                Angle[1] = dataReceive[13 + offset];
                Position[1] = new Vector(dataReceive[14 + offset], dataReceive[15 + offset]);
                RotationalAcceleration[1] = dataReceive[16 + offset];
                TranslationalAcceleration[1] = new Vector(dataReceive[17 + offset], dataReceive[18 + offset]);
            }
        }

        public virtual object Clone() {
            Motion clonedMotion = new Motion(Gravity, Density);
            clonedMotion.SetParticleArea(ParticleArea);
            clonedMotion.SetParticleMomentOfInertia(MomentOfInertia);
            return clonedMotion;
        }
    }
}

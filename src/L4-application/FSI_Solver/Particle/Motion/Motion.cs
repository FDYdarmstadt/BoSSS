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
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
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
            // creating list for motion parameters (to save the history)
            for (int i = 0; i < m_HistoryLength; i++) {
                m_Position.Add(new Vector(m_Dim));
                m_Angle.Add(new double());
                m_TranslationalVelocity.Add(new Vector(m_Dim));
                m_RotationalVelocity.Add(new double());
                m_TranslationalAcceleration.Add(new Vector(m_Dim));
                m_RotationalAcceleration.Add(new double());
                m_HydrodynamicForces.Add(new Vector(m_Dim));
                m_HydrodynamicTorque.Add(new double());
            }
        }
        internal double omega = 1;
        [NonSerialized]
        internal readonly FSI_Auxillary Aux = new FSI_Auxillary();
        [DataMember]
        private const int m_HistoryLength = 3;
        [DataMember]
        protected static int m_Dim = 2;
        [DataMember]
        private readonly List<Vector> m_Position = new List<Vector>();
        [DataMember]
        private readonly List<double> m_Angle = new List<double>();
        [DataMember]
        private readonly List<Vector> m_TranslationalVelocity = new List<Vector>();
        [DataMember]
        private readonly List<double> m_RotationalVelocity = new List<double>();
        [DataMember]
        private readonly List<Vector> m_TranslationalAcceleration = new List<Vector>();
        [DataMember]
        private readonly List<double> m_RotationalAcceleration = new List<double>();
        [DataMember]
        private readonly List<Vector> m_HydrodynamicForces = new List<Vector>();
        [DataMember]
        private readonly List<double> m_HydrodynamicTorque = new List<double>();
        [DataMember]
        private double m_CollisionTimestep = 0;
        [DataMember]
        private Vector NormalAndTangetialVelocityPreCollision = new Vector(m_Dim);
        [DataMember]
        private readonly List<Vector> m_CollisionTranslationalVelocity = new List<Vector>();
        [DataMember]
        private readonly List<double> m_CollisionRotationalVelocity = new List<double>();
        [DataMember]
        private readonly List<Vector> m_CollisionNormalVector = new List<Vector>();
        [DataMember]
        private readonly List<Vector> m_CollisionTangentialVector = new List<Vector>();
        [DataMember]
        private double m_Density;
        [DataMember]
        private double m_Area;
        [DataMember]
        private Vector m_Gravity = new Vector(0 , 0);
        [DataMember]
        private double m_MomentOfInertia;
        [DataMember]
        private double m_LengthScaleMax;
        [DataMember]
        private bool HasGhost = false;
        [DataMember]
        private int GhostID;

        /// <summary>
        /// Ghost particles are used to mirror a particle at a periodic boundary.
        /// </summary>
        internal void SetGhost(int ghostID) {
            HasGhost = true;
            GhostID = ghostID;
        }

        internal virtual int GetMasterID() => 0;

        internal bool GetHasGhost() => HasGhost;

        internal int GetGhostID() => GhostID;

        internal int GetHistoryLength() => m_HistoryLength;

        /// <summary>
        /// Gravity (volume force) acting on the particle.
        /// </summary>
        [DataMember]
        protected Vector Gravity {
            get {
                Aux.TestArithmeticException(m_Gravity, "particle density");
                return m_Gravity;
            }
            private set => m_Gravity = value;
        }

        /// <summary>
        /// Density of the particle.
        /// </summary>
        [DataMember]
        public double Density {
            get {
                Aux.TestArithmeticException(m_Density, "particle density");
                return m_Density;
            }
            private set => m_Density = value;
        }

        /// <summary>
        /// The area occupied by the particle. Calculated by <see cref="Particle"/>.
        /// </summary>
        [DataMember]
        protected double ParticleArea {
            get {
                Aux.TestArithmeticException(m_Area, "particle area");
                return m_Area;
            }
            private set => m_Area = value;
        }

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        [DataMember]
        public double Mass_P => ParticleArea * Density;

        /// <summary>
        /// The moment of inertia. Calculated by <see cref="Particle"/>.
        /// </summary>
        [DataMember]
        protected double MomentOfInertia {
            get {
                Aux.TestArithmeticException(m_MomentOfInertia, "particle moment of inertia");
                return m_MomentOfInertia;
            }
            private set => m_MomentOfInertia = value;
        }

        /// <summary>
        /// The maximum lenghtscale of the particle. Calculated by <see cref="Particle"/>.
        /// </summary>
        [DataMember]
        public double MaxParticleLengthScale {
            get {
                Aux.TestArithmeticException(m_LengthScaleMax, "particle lengthscale");
                return m_LengthScaleMax;
            }
            private set => m_LengthScaleMax = value;
        }

        /// <summary>
        /// Include rotation?
        /// </summary>
        [DataMember]
        internal virtual bool IsGhost { get; } = false;

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
        /// Complete added damping tensor, for reference: Banks et.al. 2017
        /// </summary>
        [DataMember]
        internal virtual double[,] AddedDampingTensor { get; } = new double[6, 6];

        /// <summary>
        /// Returns the position of the particle.
        /// </summary>
        /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetPosition(int historyPosition = 0) {
            if (historyPosition >= m_HistoryLength)
                throw new Exception("Error in Particle.Motion: Only " + m_HistoryLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_Position[historyPosition];
        }

        /// <summary>
        /// Returns the angle of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetAngle(int historyPosition = 0) {
            if (historyPosition >= m_HistoryLength)
                throw new Exception("Error in Particle.Motion: Only " + m_HistoryLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_Angle[historyPosition];
        }

        /// <summary>
        /// Returns the translational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetTranslationalVelocity(int historyPosition = 0) {
            if (historyPosition >= m_HistoryLength)
                throw new Exception("Error in Particle.Motion: Only " + m_HistoryLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_TranslationalVelocity[historyPosition];
        }

        /// <summary>
        /// Returns the rotational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetRotationalVelocity(int historyPosition = 0) {
            if (historyPosition >= m_HistoryLength)
                throw new Exception("Error in Particle.Motion: Only " + m_HistoryLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_RotationalVelocity[historyPosition];
        }

        /// <summary>
        /// Returns the translational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetTranslationalAcceleration(int historyPosition = 0) {
            if (historyPosition >= m_HistoryLength)
                throw new Exception("Error in Particle.Motion: Only " + m_HistoryLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_TranslationalAcceleration[historyPosition];
        }

        /// <summary>
        /// Returns the rotational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetRotationalAcceleration(int historyPosition = 0) {
            if (historyPosition >= m_HistoryLength)
                throw new Exception("Error in Particle.Motion: Only " + m_HistoryLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_RotationalAcceleration[historyPosition];
        }

        /// <summary>
        /// Returns the force acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for 4 timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetHydrodynamicForces(int historyPosition = 0) {
            if (historyPosition >= m_HistoryLength)
                throw new Exception("Error in Particle.Motion: Only " + m_HistoryLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_HydrodynamicForces[historyPosition];
        }

        /// <summary>
        /// Returns the torque acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetHydrodynamicTorque(int historyPosition = 0) {
            if (historyPosition >= m_HistoryLength)
                throw new Exception("Error in Particle.Motion: Only " + m_HistoryLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_HydrodynamicTorque[historyPosition];
        }

        /// <summary>
        /// The translational velocity of the particle before a colllision is triggered. Decomposed into an normal and tangential part. This value is used by the momentum conservation model.
        /// </summary>
        internal Vector GetNormalAndTangetialVelocityPreCollision() {
            return NormalAndTangetialVelocityPreCollision;
        }

        /// <summary>
        /// The tangential vector of the last calculated collsion.
        /// </summary>
        internal Vector GetLastCollisionTangentialVector() {
            return m_CollisionTangentialVector.Last();
        }


        /// <summary>
        /// Saves position and angle of the last timestep.
        /// </summary>
        public void SavePositionAndAngleOfPreviousTimestep() {
            Aux.SaveVectorOfLastTimestep(m_Position);
            Aux.SaveValueOfLastTimestep(m_Angle);
        }

        /// <summary>
        /// Saves translational and rotational velocities of the last timestep.
        /// </summary>
        public void SaveVelocityOfPreviousTimestep() {
            Aux.SaveVectorOfLastTimestep(m_TranslationalVelocity);
            Aux.SaveValueOfLastTimestep(m_RotationalVelocity);
            Aux.SaveVectorOfLastTimestep(m_TranslationalAcceleration);
            Aux.SaveValueOfLastTimestep(m_RotationalAcceleration);
        }

        /// <summary>
        /// Saves force and torque of the previous timestep.
        /// </summary>
        public void SaveHydrodynamicsOfPreviousTimestep() {
            Aux.SaveVectorOfLastTimestep(m_HydrodynamicForces);
            for (int i = 0; i < m_HistoryLength; i++)
            Aux.SaveValueOfLastTimestep(m_HydrodynamicTorque);
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
        public void InitializeParticlePositionAndAngle(double[] initialPosition, double initialAngle) {
            for (int i = 0; i < m_HistoryLength; i++) {
                m_Position[i] = new Vector(initialPosition);
                m_Angle[i] = initialAngle * 2 * Math.PI / 360;
                Aux.TestArithmeticException(m_Position[i], "initial particle position");
                Aux.TestArithmeticException(m_Angle[i], "initial particle angle");
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
        public void InitializeParticleVelocity(double[] initalTranslation, double initalRotation) {
            for (int i = 0; i < m_HistoryLength; i++) {
                m_TranslationalVelocity[i] = initalTranslation == null ? new Vector(m_Dim) : new Vector(initalTranslation);
                m_RotationalVelocity[i] = initalRotation;
                Aux.TestArithmeticException(m_TranslationalVelocity[i], "initial particle translational velocity");
                Aux.TestArithmeticException(m_RotationalVelocity[i], "initial particle rotational velocity");
            }
        }

        public void TransferPhysicalData(Motion motionDataToTransfer) {
            for(int h = 0; h < m_HistoryLength; h++) {
                m_Position[h] = new Vector(motionDataToTransfer.GetPosition(h));
                m_TranslationalVelocity[h] = new Vector(motionDataToTransfer.GetTranslationalVelocity(h));
                m_TranslationalAcceleration[h] = new Vector(motionDataToTransfer.GetTranslationalAcceleration(h));
                m_HydrodynamicForces[h] = new Vector(motionDataToTransfer.GetHydrodynamicForces(h));
                m_Angle[h] = motionDataToTransfer.GetAngle(h);
                m_RotationalVelocity[h] = motionDataToTransfer.GetRotationalVelocity(h);
                m_RotationalAcceleration[h] = motionDataToTransfer.GetRotationalAcceleration(h);
                m_HydrodynamicTorque[h] = motionDataToTransfer.GetHydrodynamicTorque(h);
            }
        }

        public void AdaptToNewTimestep(double newTimestep, double oldTimestep) {
            AdaptNewTimestep(m_Position, newTimestep, oldTimestep);
            AdaptNewTimestep(m_TranslationalVelocity, newTimestep, oldTimestep);
            AdaptNewTimestep(m_TranslationalAcceleration, newTimestep, oldTimestep);
            AdaptNewTimestep(m_Angle, newTimestep, oldTimestep);
            AdaptNewTimestep(m_RotationalVelocity, newTimestep, oldTimestep);
            AdaptNewTimestep(m_RotationalAcceleration, newTimestep, oldTimestep);
        }

        private void AdaptNewTimestep(List<Vector> variable, double newTimestep, double oldTimestep) {
            List<Vector> stuetzstelle = new List<Vector>();
            for (int h = 0; h <  m_HistoryLength; h++) {
                stuetzstelle.Add(new Vector(variable[h]));
            }
            List<Vector> newVariable = new List<Vector> { new Vector(variable[0]) };
            for (int i = 1; i < m_HistoryLength; i++) {
                double currentNewTimestep = -i * newTimestep;
                double[] langrangePoly = CalculateLangrangePolynom(currentNewTimestep, oldTimestep);
                newVariable.Add(new Vector(langrangePoly[0] * variable[0]));
                for (int j = 1; j < m_HistoryLength; j++) {
                    newVariable[i] = newVariable[i] + langrangePoly[j] * stuetzstelle[j];
                }
            }
            variable.Clear();
            for (int h = 0; h < m_HistoryLength; h++) {
                variable.Add(new Vector(newVariable[h]));
            }
        }

        private void AdaptNewTimestep(List<double> variable, double newTimestep, double oldTimestep) {
            List<double> stuetzstelle = new List<double>();
            for (int h = 0; h < m_HistoryLength; h++) {
                stuetzstelle.Add(variable[h]);
            }
            List<double> newVariable = new List<double>();
            newVariable.Add(variable[0]);
            for (int i = 1; i < m_HistoryLength; i++) {
                double currentNewTimestep = -i * newTimestep;
                double[] langrangePoly = CalculateLangrangePolynom(currentNewTimestep, oldTimestep);
                newVariable.Add(langrangePoly[0] * variable[i]);
                for (int j = 1; j < m_HistoryLength; j++) {
                    newVariable[i] = newVariable[i] + langrangePoly[j] * stuetzstelle[j];
                }
            }
            variable.Clear();
            for (int h = 0; h < m_HistoryLength; h++) {
                variable.Add(newVariable[h]);
            }
        }

        private double[] CalculateLangrangePolynom(double time, double oldTimestep) {
            double[] lPoly = new double[m_HistoryLength];
            for (int i = 0; i < m_HistoryLength; i++) {
                if (i != 0)
                    lPoly[i] = time / (-i * oldTimestep);
                for (int j = 1; j < m_HistoryLength; j++) {
                    if (j == i)
                        continue;
                    if (i == 0 && j == 1)
                        lPoly[i] = (time + j * oldTimestep) / ((j - i) * oldTimestep);
                    lPoly[i] *= (time + j * oldTimestep) / ((j - i) * oldTimestep);
                }
            }
            return lPoly;
        }

        /// <summary>
        /// Init of the particle area.
        /// </summary>
        /// <param name="area"></param>
        public void GetParticleArea(double area) => ParticleArea = area;

        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void GetParticleLengthscale(double lengthscale) => MaxParticleLengthScale = lengthscale;

        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void GetParticleMomentOfInertia(double moment) => MomentOfInertia = moment;

        /// <summary>
        /// Sets the collision timestep.
        /// </summary>
        /// <param name="collisionTimestep">
        /// The physical time consumend by the collision procedure.
        /// </param>
        public void SetCollisionTimestep(double collisionTimestep) => m_CollisionTimestep = collisionTimestep;

        /// <summary>
        /// Calls the calculation of the position and angle.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticlePositionAndAngle(double dt) {
            if (m_CollisionTimestep == 0) {
                SavePositionAndAngleOfPreviousTimestep();
                m_Position[0] = CalculateParticlePosition(dt);
                m_Angle[0] = CalculateParticleAngle(dt);
            }
            else {
                if (m_CollisionTimestep < 0)
                    m_CollisionTimestep = 0;
                Vector tempPos = new Vector(m_Position[0]);
                double tempAngle = m_Angle[0];
                SavePositionAndAngleOfPreviousTimestep();
                m_Position[0] = new Vector(tempPos);
                m_Angle[0] = tempAngle;
                m_Position[0] = CalculateParticlePosition(dt, m_CollisionTimestep);
                m_Angle[0] = CalculateParticleAngle(dt, m_CollisionTimestep);
                if (m_CollisionTimestep > dt) { m_CollisionTimestep -= dt; }
                else m_CollisionTimestep = 0;
            }
        }

        /// <summary>
        /// Calls the calculation of the position and angle.
        /// </summary>
        /// <param name="dt"></param>
        public void CollisionParticlePositionAndAngle(double collisionDynamicTimestep) {
            m_Position[0] = CalculateParticlePosition(collisionDynamicTimestep, collisionProcedure: true);
            m_Angle[0] = CalculateParticleAngle(collisionDynamicTimestep, collisionProcedure: true);
        }

        /// <summary>
        /// Calls the calculation of the velocity.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticleVelocity(double dt) {
            m_TranslationalAcceleration[0] = CalculateTranslationalAcceleration(dt - m_CollisionTimestep);
            m_RotationalAcceleration[0] = CalculateRotationalAcceleration(dt - m_CollisionTimestep);
            if (m_CollisionTimestep == 0) {
                m_TranslationalVelocity[0] = CalculateTranslationalVelocity(dt);
                m_RotationalVelocity[0] = CalculateAngularVelocity(dt);
            }
            else {
                m_TranslationalVelocity[0] = CalculateTranslationalVelocity(dt, m_CollisionTimestep);
                m_RotationalVelocity[0] = CalculateAngularVelocity(dt, m_CollisionTimestep);
            }
        }

        internal virtual void CopyNewPosition(Vector position, double angle) { }

        internal virtual void SetGhostPosition(Vector position) {
            for (int h = 0; h<m_HistoryLength; h++) {
                m_Position[h] = new Vector(position);
            }
        }

        internal virtual void CopyNewVelocity(Vector translational, double rotational) { }

        public void UpdateForcesAndTorque(int particleID, double[] fullListHydrodynamics) {
            double[] tempForces = new double[m_Dim];
            for (int d = 0; d < m_Dim; d++) {
                if (Math.Abs(fullListHydrodynamics[particleID * 3 + d]) > 1e-8)
                    tempForces[d] = fullListHydrodynamics[particleID * 3 + d];
            }
            m_HydrodynamicForces[0] = new Vector(tempForces);
            if (Math.Abs(fullListHydrodynamics[particleID * 3 + m_Dim]) > 1e-8)
                m_HydrodynamicTorque[0] = fullListHydrodynamics[particleID * 3 + m_Dim];
            Aux.TestArithmeticException(m_HydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(m_HydrodynamicTorque[0], "hydrodynamic torque");
        }

        /// <summary>
        /// Predicts the hydrodynamics at the beginning of the iteration loop in each timestep.
        /// </summary>
        /// <param name="activeStress"></param>
        /// <param name="timestepID">
        /// The timestep ID. Used to distinguish between the first timestep and all other steps.
        /// </param>
        public virtual void PredictForceAndTorque(double activeStress, double circumference, int timestepID, double fluidViscosity, double fluidDensity, double dt) {
            m_HydrodynamicForces[0] = new Vector(m_HydrodynamicForces[1]);
            m_HydrodynamicTorque[0] = m_HydrodynamicTorque[1];
            Aux.TestArithmeticException(m_HydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(m_HydrodynamicTorque[0], "hydrodynamic torque");
        }

        /// <summary>
        /// Calculate the tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        /// <param name="particle"></param>
        /// <param name="levelSetTracker"></param>
        /// <param name="fluidViscosity"></param>
        /// <param name="fluidDensity"></param>
        /// <param name="dt"></param>
        public virtual void CalculateDampingTensor(Particle particle, LevelSetTracker levelSetTracker, double fluidViscosity, double fluidDensity, double dt) {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        /// <summary>
        /// Updates the tensors to implement the added damping model (Banks et.al. 2017)
        /// </summary>
        public virtual void UpdateDampingTensors() {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        /// <summary>
        /// Exchange of the added damping tensors between the MPI-processes.
        /// </summary>
        public virtual void ExchangeAddedDampingTensors() {
            throw new Exception("Added damping tensors should only be used if added damping is active.");
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        protected virtual Vector CalculateParticlePosition(double dt) {
            Vector l_Position = m_Position[1] + (m_TranslationalVelocity[0] + 4 * m_TranslationalVelocity[1] + m_TranslationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle position during the collision procedure.
        /// </summary>
        /// <param name="dt"></param>
        protected virtual Vector CalculateParticlePosition(double dt, bool collisionProcedure) {
            Vector l_Position = m_Position[0] + (m_TranslationalVelocity[0] + 4 * m_TranslationalVelocity[1] + m_TranslationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle position after a collision
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual Vector CalculateParticlePosition(double dt, double collisionTimestep) {
            Vector l_Position = m_Position[0] + m_TranslationalVelocity[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateParticleAngle(double dt) {
            double l_Angle = m_Angle[1] + (m_RotationalVelocity[0] + 4 * m_RotationalVelocity[1] + m_RotationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(l_Angle, "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new particle angle during the collision procedure.
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateParticleAngle(double dt, bool collisionProcedure) {
            double l_Angle = m_Angle[0] + (m_RotationalVelocity[0] + 4 * m_RotationalVelocity[1] + m_RotationalVelocity[2]) * dt / 6;
            Aux.TestArithmeticException(l_Angle, "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new particle angle after a collision
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual double CalculateParticleAngle(double dt, double collisionTimestep) {
            double l_Angle = m_Angle[1] + m_RotationalVelocity[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(l_Angle, "particle angle");
            return l_Angle;
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected virtual Vector CalculateTranslationalVelocity(double dt) {
            Vector l_TranslationalVelocity = m_TranslationalVelocity[1] + (m_TranslationalAcceleration[0] + 4 * m_TranslationalAcceleration[1] + m_TranslationalAcceleration[2]) * dt / 6;
            for (int d = 0; d < l_TranslationalVelocity.Dim; d++) {
                if (Math.Abs(l_TranslationalVelocity[d]) < 1e-8)
                    l_TranslationalVelocity[d] = 0;
            }
            Aux.TestArithmeticException(l_TranslationalVelocity, "particle translational velocity");
            return l_TranslationalVelocity;
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle after a collision.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual Vector CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            Vector l_TranslationalVelocity = m_TranslationalVelocity[1] + m_TranslationalAcceleration[0] * (dt - collisionTimestep) / 6;
            for (int d = 0; d < l_TranslationalVelocity.Dim; d++) {
                if (Math.Abs(l_TranslationalVelocity[d]) < 1e-8)
                    l_TranslationalVelocity[d] = 0;
            }
            Aux.TestArithmeticException(l_TranslationalVelocity, "particle translational velocity");
            return l_TranslationalVelocity;
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        protected virtual double CalculateAngularVelocity(double dt) {
            double l_RotationalVelocity = m_RotationalVelocity[1] + (m_RotationalAcceleration[0] + 4 * m_RotationalAcceleration[1] + m_RotationalAcceleration[2]) * dt / 6;
            Aux.TestArithmeticException(l_RotationalVelocity, "particle rotational velocity");
            if (Math.Abs(l_RotationalVelocity) > 1e-8)
                return l_RotationalVelocity;
            else
                return 0;
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle after a collision.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual double CalculateAngularVelocity(double dt, double collisionTimestep) {
            double l_RotationalVelocity = m_RotationalVelocity[1] + m_RotationalAcceleration[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(l_RotationalVelocity, "particle rotational velocity");
            if (l_RotationalVelocity > 1e-8)
                return l_RotationalVelocity;
            else
                return 0;
        }

        /// <summary>
        /// Calculates the velocities along the normal and the tangential vector of a collision.
        /// </summary>
        internal void CalculateNormalAndTangentialVelocity() {
            Vector normalVector = m_CollisionNormalVector.Last();
            Vector TangentialVector = new Vector(- normalVector[1], normalVector[0]);
            Vector Velocity = m_TranslationalVelocity[0];
            NormalAndTangetialVelocityPreCollision = new Vector(Velocity * normalVector, Velocity * TangentialVector);
            Aux.TestArithmeticException(NormalAndTangetialVelocityPreCollision, "particle velocity before collision");
        }

        /// <summary>
        /// Calculate the new tranlational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        protected virtual Vector CalculateTranslationalAcceleration(double dt) {
            Vector l_Acceleration = m_HydrodynamicForces[0] / (Density * ParticleArea);
            Aux.TestArithmeticException(l_Acceleration, "particle translational acceleration");
            return l_Acceleration;
        }

        protected double[][] TransformStressToPrint(List<double[]>[] stressToPrintOut) {
            if (stressToPrintOut[0].Count() != stressToPrintOut[1].Count())
                throw new Exception("Something strange happend!");
            double[][] output = new double[stressToPrintOut[0].Count()][];
            for (int d = 0; d < m_Dim; d++) {
                for (int i = stressToPrintOut[d].Count() - 1; i > 0; i--) {
                    for (int j = 0; j < i - 1; j++) {
                        if (stressToPrintOut[d][j][0] > stressToPrintOut[d][j + 1][0]) {
                            double[] temp = stressToPrintOut[d][j].CloneAs();
                            stressToPrintOut[d][j] = stressToPrintOut[d][j + 1].CloneAs();
                            stressToPrintOut[d][j + 1] = temp;
                        }
                    }
                }
            }
            for (int i = 0; i < output.Length; i++) {
                if (Math.Abs(stressToPrintOut[0][i][0] - stressToPrintOut[1][i][0]) > 1e-15)
                    throw new Exception("Something strange happend!");
                double surfaceParam = stressToPrintOut[0][i][0];
                double normalStress = Math.Cos(surfaceParam) * stressToPrintOut[0][i][1] + Math.Sin(stressToPrintOut[0][i][0]) * stressToPrintOut[1][i][1];
                double tangentialStress = -Math.Sin(surfaceParam) * stressToPrintOut[0][i][1] + Math.Cos(stressToPrintOut[0][i][0]) * stressToPrintOut[1][i][1];
                surfaceParam = Math.PI * (1 - Math.Sign(-Math.Sin(m_Angle[0]) * Math.Cos(surfaceParam) + Math.Cos(m_Angle[0]) * Math.Sin(surfaceParam))) / 2 + Math.Acos(Math.Cos(m_Angle[0]) * Math.Cos(surfaceParam) + Math.Sin(m_Angle[0]) * Math.Sin(surfaceParam));
                double[] insert = new double[] { surfaceParam, normalStress, tangentialStress };
                output[i] = insert;
            }
            return output;
        }

        /// <summary>
        /// Calculate the new rotational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double CalculateRotationalAcceleration(double dt) {
            double l_Acceleration = m_HydrodynamicTorque[0] / MomentOfInertia;
            Aux.TestArithmeticException(l_Acceleration, "particle rotational acceleration");
            return l_Acceleration;
        }

        private TextWriter logStress;

        /// <summary>
        /// Creates a log file for the residum of the hydrodynamic forces.
        /// </summary>
        public void CreateStressLogger(SessionInfo CurrentSessionInfo, IDatabaseDriver DatabaseDriver, double phystime, int particleID) {
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MPIRank);
            if ((MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                string name = "stress_Time_" + phystime.ToString() + "_particle_" + particleID.ToString();
                logStress = DatabaseDriver.FsDriver.GetNewLog(name, CurrentSessionInfo.ID);
                logStress.WriteLine(string.Format("{0},{1},{2},{3}", "Time", "surfaceParam", "stressNormal", "stressTangential"));
            }
        }

        /// <summary>
        /// Creates a log file for the residum of the hydrodynamic forces.
        /// </summary>
        public void LogStress(double phystime) {
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MPIRank);
            if ((MPIRank == 0) && (logStress != null)) {
                for (int i = 0; i <  currentStress.Length; i++) {
                    logStress.WriteLine(string.Format("{0},{1},{2},{3}", phystime, currentStress[i][0], currentStress[i][1], currentStress[i][2]));
                    logStress.Flush();
                }
            }
        }

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        public virtual Vector CalculateHydrodynamicForces(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, double fluidDensity, CellMask cutCells, double dt = 0) {
            Vector tempForces = new Vector(hydrodynamicsIntegration.Forces(out List<double[]>[] stressToPrintOut, cutCells));
            currentStress = TransformStressToPrint(stressToPrintOut);
            Aux.TestArithmeticException(tempForces, "temporal forces during calculation of hydrodynamics");
            tempForces = Force_MPISum(tempForces);
            tempForces = CalculateGravity(fluidDensity, tempForces);
            return tempForces;
        }

        public bool UseConstantUnderrelaxation() {
            bool useConstantUnderrelaxation = false;
            for (int d = 0; d < m_Dim; d++) {
                double averageTrans = 0;
                for (int i = 1; i < m_HistoryLength; i++) {
                    averageTrans += m_TranslationalVelocity[i][d];
                }
                averageTrans /= m_HistoryLength;
                double squaresTrans = 0;
                for (int i =  1; i < m_HistoryLength; i++) {
                    squaresTrans += (m_TranslationalVelocity[i][d] - averageTrans).Pow2() / (m_TranslationalVelocity[i][d]).Pow2();
                }
                squaresTrans /= m_HistoryLength - 1;
                if(d == 0)
                Console.WriteLine("squareTrans " + squaresTrans);
                if (squaresTrans < 0.01563 && averageTrans != 0)
                    useConstantUnderrelaxation = true;
            }
            double averageRot = 0;
            for (int i = 1; i < m_HistoryLength; i++) {
                averageRot += m_RotationalVelocity[i];
            }
            averageRot /= m_HistoryLength;
            double squaresRot = 0;
            for (int i = 1; i < m_HistoryLength; i++) {
                squaresRot += (m_RotationalVelocity[i] - averageRot).Pow2() / m_RotationalVelocity[i].Pow2();
            }
            squaresRot /= m_HistoryLength - 1;
            if (squaresRot < 0.01563 && averageRot != 0)
                useConstantUnderrelaxation = true;
            return useConstantUnderrelaxation;
        }

        protected double[][] currentStress;

        /// <summary>
        /// Calculates the gravitational forces.
        /// </summary>
        /// <param name="fluidDensity"></param>
        /// <param name="tempForces"></param>
        protected Vector CalculateGravity(double fluidDensity, Vector tempForces) {
            tempForces += (Density - fluidDensity) * ParticleArea * Gravity;
            Aux.TestArithmeticException(tempForces, "temporal forces during calculation of hydrodynamics after adding gravity");
            return tempForces;
        }

        /// <summary>
        /// Summation of the hydrodynamic forces over all MPI-processes
        /// </summary>
        /// <param name="forces"></param>
        protected Vector Force_MPISum(Vector forces) {
            double[] stateBuffer = ((double[])forces).CloneAs();
            double[] globalStateBuffer = stateBuffer.MPISum();
            forces = new Vector(globalStateBuffer);
            Aux.TestArithmeticException(forces, "temporal forces during calculation of hydrodynamics after mpi-summation");
            return forces;
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
        public virtual double CalculateHydrodynamicTorque(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, double dt = 0) {
            double tempTorque = hydrodynamicsIntegration.Torque(m_Position[0], cutCells);
            Aux.TestArithmeticException(tempTorque, "temporal torque during calculation of hydrodynamics");
            Torque_MPISum(ref tempTorque);
            return tempTorque;
        }

        /// <summary>
        /// Summation of the hydrodynamic torque over all MPI-processes
        /// </summary>
        /// <param name="torque"></param>
        protected void Torque_MPISum(ref double torque) {
            double stateBuffer = torque;
            double globalStateBuffer = stateBuffer.MPISum();
            torque = globalStateBuffer;
            Aux.TestArithmeticException(torque, "temporal torque during calculation of hydrodynamics after mpi-summation");
        }

        /// <summary>
        /// Calculating the particle Reynolds number
        /// </summary>
        /// <param name="fluidViscosity"></param>
        public double ComputeParticleRe(double fluidViscosity) {
            return 2 * m_TranslationalVelocity[0].L2Norm() * MaxParticleLengthScale / fluidViscosity;
        }

        /// <summary>
        /// Calculating the particle Stokes number
        /// </summary>
        /// <param name="fluidViscosity"></param>
        /// <param name="fluidDensity"></param>
        public double ComputeParticleSt(double fluidViscosity, double fluidDensity) {
            return ComputeParticleRe(fluidViscosity) * Density / (9 * fluidDensity);
        }

        /// <summary>
        /// Calculating the particle momentum
        /// </summary>
        public double[] CalculateParticleMomentum() {
            double[] temp = new double[m_Dim + 1];
            for (int d = 0; d < m_Dim; d++) {
                temp[d] = Mass_P * m_TranslationalVelocity[0][d];
            }
            temp[m_Dim] = MomentOfInertia * m_RotationalVelocity[0];
            return temp;
        }

        /// <summary>
        /// Calculating the particle kinetic energy
        /// </summary>
        public double[] CalculateParticleKineticEnergy() {
            double[] temp = new double[m_Dim + 1];
            for (int d = 0; d < m_Dim; d++) {
                temp[d] = 0.5 * Mass_P * m_TranslationalVelocity[0][d].Pow2();
            }
            temp[m_Dim] = 0.5 * MomentOfInertia * m_RotationalVelocity[0].Pow2();
            return temp;
        }

        /// <summary>
        /// Deletes the complete history of the translational velocity and acceleration.
        /// </summary>
        public void ClearParticleHistoryTranslation() {
            for (int i = 0; i < m_TranslationalVelocity.Count; i++) {
                m_TranslationalVelocity[i] = new Vector(m_Dim);
                m_TranslationalAcceleration[i] = new Vector(m_Dim);
            }
        }

        /// <summary>
        /// Deletes the complete history of the rotational velocity and acceleration.
        /// </summary>
        public void ClearParticleHistoryRotational() {
            for (int i = 0; i < m_RotationalVelocity.Count; i++) {
                m_RotationalAcceleration[i] = 0;
                m_RotationalVelocity[i] = 0;
            }
        }

        /// <summary>
        /// Sets the normal and tangential vectors related to the current collision.
        /// </summary>
        /// <param name="normalVector"></param>
        /// <param name="tangentialVector"></param>
        public void SetCollisionVectors(Vector normalVector, Vector tangentialVector) {
            m_CollisionNormalVector.Add(normalVector);
            m_CollisionTangentialVector.Add(tangentialVector);
        }

        /// <summary>
        /// Sets the normal, tangential and rotational velocity after the current collision.
        /// </summary>
        /// <param name="normalVelocity"></param>
        /// <param name="tangentialVelocity"></param>
        /// <param name="rotationalVelocity"></param>
        public void SetCollisionVelocities(double normalVelocity, double tangentialVelocity, double rotationalVelocity) {
            m_CollisionTranslationalVelocity.Add(new Vector(normalVelocity, tangentialVelocity));
            m_CollisionRotationalVelocity.Add(rotationalVelocity);
        }

        /// <summary>
        /// Collision post-processing. Sums up the results of the multiple binary collisions of one timestep
        /// </summary>
        public void PostProcessCollisionTranslation() {
            if (m_CollisionTranslationalVelocity.Count >= 1) {
                Vector normal = new Vector(m_Dim);
                Vector tangential = new Vector(m_Dim);
                for (int t = 0; t < m_CollisionTranslationalVelocity.Count; t++) {
                    normal += m_CollisionNormalVector[t];
                    tangential += m_CollisionTangentialVector[t];
                }
                normal.Normalize();
                tangential.Normalize();
                double temp_NormalVel = 0;
                double temp_TangentialVel = 0;
                for (int t = 0; t < m_CollisionTranslationalVelocity.Count; t++) {
                    double cos = normal * m_CollisionNormalVector[t];
                    double sin = cos == 1 ? 0 : m_CollisionNormalVector[t][0] > normal[0] ? Math.Sqrt(1 + 1e-15 - cos.Pow2()) : -Math.Sqrt(1 + 1e-15 - cos.Pow2());
                    temp_NormalVel += m_CollisionTranslationalVelocity[t][0] * cos - m_CollisionTranslationalVelocity[t][1] * sin;
                    temp_TangentialVel += m_CollisionTranslationalVelocity[t][0] * sin + m_CollisionTranslationalVelocity[t][1] * cos;

                }
                temp_NormalVel /= m_CollisionTranslationalVelocity.Count;
                temp_TangentialVel /= m_CollisionTranslationalVelocity.Count;

                ClearParticleHistoryTranslation();
                m_TranslationalVelocity[0] = normal * temp_NormalVel + tangential * temp_TangentialVel;
                m_CollisionTranslationalVelocity.Clear();
                m_CollisionNormalVector.Clear();
                m_CollisionTangentialVector.Clear();
            }
        }

        /// <summary>
        /// Collision post-processing. Sums up the results for the rotatinal velocity of the multiple binary collisions of one timestep
        /// </summary>
        public void PostProcessCollisionRotation() {
            if (m_CollisionRotationalVelocity.Count >= 1) {
                ClearParticleHistoryRotational();
                m_RotationalVelocity[0] = m_CollisionRotationalVelocity.Sum() / m_CollisionRotationalVelocity.Count;
                m_CollisionRotationalVelocity.Clear();
                if (double.IsNaN(m_RotationalVelocity[0]) || double.IsInfinity(m_RotationalVelocity[0]))
                    throw new ArithmeticException("Error trying to update particle angular velocity during collision post-processing. The angular velocity is:  " + m_RotationalVelocity[0]);
            }
        }

        /// <summary>
        /// Builds the array for the post-collision communication between MPI-processes.
        /// </summary>
        public double[] BuildSendArray() {
            double[] dataSend = new double[19];
            dataSend[0] = m_RotationalVelocity[0];
            dataSend[1] = m_RotationalAcceleration[0];
            dataSend[2] = m_TranslationalVelocity[0][0];
            dataSend[3] = m_TranslationalVelocity[0][1];
            dataSend[4] = m_TranslationalAcceleration[0][0];
            dataSend[5] = m_TranslationalAcceleration[0][1];
            dataSend[6] = m_Angle[0];
            dataSend[7] = m_Position[0][0];
            dataSend[8] = m_Position[0][1];
            dataSend[9] = m_CollisionTimestep;
            dataSend[10] = m_RotationalVelocity[1];
            dataSend[11] = m_TranslationalVelocity[1][0];
            dataSend[12] = m_TranslationalVelocity[1][1];
            dataSend[13] = m_Angle[1];
            dataSend[14] = m_Position[1][0];
            dataSend[15] = m_Position[1][1];
            dataSend[16] = m_RotationalAcceleration[1];
            dataSend[17] = m_TranslationalAcceleration[1][0];
            dataSend[18] = m_TranslationalAcceleration[1][1];
            return dataSend;
        }

        /// <summary>
        /// Overwrites the particles parameters with the values received during the post-collision MPI-communication.
        /// </summary>
        public void WriteReceiveArray(double[] dataReceive, int offset) {
            m_RotationalVelocity[0] = dataReceive[0 + offset];
            m_RotationalAcceleration[0] = dataReceive[1 + offset];
            m_TranslationalVelocity[0] = new Vector(dataReceive[2 + offset], dataReceive[3 + offset]);
            m_TranslationalAcceleration[0] = new Vector(dataReceive[4 + offset], dataReceive[5 + offset]);
            m_Angle[0] = dataReceive[6 + offset];
            m_Position[0] = new Vector(dataReceive[7 + offset], dataReceive[8 + offset]);
            m_CollisionTimestep = dataReceive[9 + offset];
            m_RotationalVelocity[1] = dataReceive[10 + offset];
            m_TranslationalVelocity[1] = new Vector(dataReceive[11 + offset], dataReceive[12 + offset]);
            m_Angle[1] = dataReceive[13 + offset];
            m_Position[1] = new Vector(dataReceive[14 + offset], dataReceive[15 + offset]);
            m_RotationalAcceleration[1] = dataReceive[16 + offset];
            m_TranslationalAcceleration[1] = new Vector(dataReceive[17 + offset], dataReceive[18 + offset]);
        }

        public virtual object Clone() {
            Motion clonedMotion = new Motion(Gravity, Density);
            clonedMotion.GetParticleArea(ParticleArea);
            clonedMotion.GetParticleMomentOfInertia(MomentOfInertia);
            return clonedMotion;
        }
    }
}

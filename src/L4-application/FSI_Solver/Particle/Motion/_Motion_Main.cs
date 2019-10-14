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
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.FSI_Solver {
    public class Motion_Wet {

        /// <summary>
        /// The standard description of motion including hydrodynamics.
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
        public Motion_Wet(double[] gravity, double density, ParticleUnderrelaxationParam underrelaxationParam = null) {
            Gravity = gravity;
            m_UnderrelaxationParam = underrelaxationParam;
            Density = density;
            // creating list for motion parameters (to save the history)
            for (int i = 0; i < historyLength; i++) {
                m_Position.Add(new double[spatialDim]);
                m_Angle.Add(new double());
                m_TranslationalVelocity.Add(new double[spatialDim]);
                m_RotationalVelocity.Add(new double());
                m_TranslationalAcceleration.Add(new double[spatialDim]);
                m_RotationalAcceleration.Add(new double());
                m_HydrodynamicForces.Add(new double[spatialDim]);
                m_HydrodynamicTorque.Add(new double());
            }
        }

        [NonSerialized]
        internal readonly FSI_Auxillary Aux = new FSI_Auxillary();
        [NonSerialized]
        private readonly ParticleUnderrelaxationParam m_UnderrelaxationParam = null;

        [DataMember]
        private const int historyLength = 4;
        [DataMember]
        protected static int spatialDim = 2;
        [DataMember]
        private readonly List<double[]> m_Position = new List<double[]>();
        [DataMember]
        private readonly List<double> m_Angle = new List<double>();
        [DataMember]
        private readonly List<double[]> m_TranslationalVelocity = new List<double[]>();
        [DataMember]
        private readonly List<double> m_RotationalVelocity = new List<double>();
        [DataMember]
        private readonly List<double[]> m_TranslationalAcceleration = new List<double[]>();
        [DataMember]
        private readonly List<double> m_RotationalAcceleration = new List<double>();
        [DataMember]
        private readonly List<double[]> m_HydrodynamicForces = new List<double[]>();
        [DataMember]
        private readonly List<double> m_HydrodynamicTorque = new List<double>();
        [DataMember]
        private double[] m_ForcesPreviousIteration = new double[spatialDim];
        [DataMember]
        private double m_TorquePreviousIteration = new double();
        [DataMember]
        private double m_CollisionTimestep = 0;
        [DataMember]
        private double[] m_PreCollisionVelocity = new double[spatialDim];
        [DataMember]
        private double particleArea;
        [DataMember]
        private readonly List<double[]> m_CollisionTranslationalVelocity = new List<double[]>();
        [DataMember]
        private readonly List<double> m_CollisionRotationalVelocity = new List<double>();
        [DataMember]
        private readonly List<double[]> m_CollisionNormalVector = new List<double[]>();
        [DataMember]
        private readonly List<double[]> m_CollisionTangentialVector = new List<double[]>();
        [DataMember]
        private double density;
        [DataMember]
        private double[] gravity;
        [DataMember]
        private double momentOfInertia;
        [DataMember]
        private double maxParticleLengthScale;

        /// <summary>
        /// Gravity (volume force) acting on the particle.
        /// </summary>
        [DataMember]
        protected double[] Gravity {
            get {
                Aux.TestArithmeticException(gravity, "particle density");
                return gravity;
            }
            private set => gravity = value;
        }

        /// <summary>
        /// Density of the particle.
        /// </summary>
        [DataMember]
        public double Density {
            get {
                Aux.TestArithmeticException(density, "particle density");
                return density;
            }
            private set => density = value;
        }

        /// <summary>
        /// The area occupied by the particle. Calculated by <see cref="Particle"/>.
        /// </summary>
        [DataMember]
        protected double ParticleArea {
            get {
                Aux.TestArithmeticException(particleArea, "particle area");
                return particleArea;
            }
            private set => particleArea = value;
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
                Aux.TestArithmeticException(momentOfInertia, "particle moment of inertia");
                return momentOfInertia;
            }
            private set => momentOfInertia = value;
        }

        /// <summary>
        /// The maximum lenghtscale of the particle. Calculated by <see cref="Particle"/>.
        /// </summary>
        [DataMember]
        protected double MaxParticleLengthScale {
            get {
                Aux.TestArithmeticException(maxParticleLengthScale, "particle moment of inertia");
                return maxParticleLengthScale;
            }
            private set => maxParticleLengthScale = value;
        }

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
        internal double[] GetPosition(int historyPosition = 0) {
            if (historyPosition >= historyLength)
                throw new Exception("Error in Particle.Motion: Only " + historyLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_Position[historyPosition];
        }

        /// <summary>
        /// Returns the angle of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetAngle(int historyPosition = 0) {
            if (historyPosition >= historyLength)
                throw new Exception("Error in Particle.Motion: Only " + historyLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_Angle[historyPosition];
        }

        /// <summary>
        /// Returns the translational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double[] GetTranslationalVelocity(int historyPosition = 0) {
            if (historyPosition >= historyLength)
                throw new Exception("Error in Particle.Motion: Only " + historyLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_TranslationalVelocity[historyPosition];
        }

        /// <summary>
        /// Returns the rotational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetRotationalVelocity(int historyPosition = 0) {
            if (historyPosition >= historyLength)
                throw new Exception("Error in Particle.Motion: Only " + historyLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_RotationalVelocity[historyPosition];
        }

        /// <summary>
        /// Returns the translational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double[] GetTranslationalAcceleration(int historyPosition = 0) {
            if (historyPosition >= historyLength)
                throw new Exception("Error in Particle.Motion: Only " + historyLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_TranslationalAcceleration[historyPosition];
        }

        /// <summary>
        /// Returns the rotational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetRotationalAcceleration(int historyPosition = 0) {
            if (historyPosition >= historyLength)
                throw new Exception("Error in Particle.Motion: Only " + historyLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_RotationalAcceleration[historyPosition];
        }

        /// <summary>
        /// Returns the force acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for 4 timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double[] GetHydrodynamicForces(int historyPosition = 0) {
            if (historyPosition >= historyLength)
                throw new Exception("Error in Particle.Motion: Only " + historyLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_HydrodynamicForces[historyPosition];
        }

        /// <summary>
        /// Returns the torque acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four timesteps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetHydrodynamicTorque(int historyPosition = 0) {
            if (historyPosition >= historyLength)
                throw new Exception("Error in Particle.Motion: Only " + historyLength + " timesteps are saved. The requested value is " + historyPosition + " steps in the past!");
            return m_HydrodynamicTorque[historyPosition];
        }

        /// <summary>
        /// Returns the force of the previous iteration.
        /// </summary>
        internal double[] GetForcesPreviousIteration() {
            return m_ForcesPreviousIteration;
        }

        /// <summary>
        /// Returns the torque of the previous iteration.
        /// </summary>
        internal double GetTorquePreviousIteration() {
            return m_TorquePreviousIteration;
        }

        /// <summary>
        /// The translational velocity of the particle before a colllision is triggered. This value is used by the momentum conservation model.
        /// </summary>
        internal double[] GetPreCollisionVelocity() {
            return m_PreCollisionVelocity;
        }

        /// <summary>
        /// The tangential vector of the last calculated collsion.
        /// </summary>
        internal double[] GetLastCollisionTangentialVector() {
            return m_CollisionTangentialVector.Last();
        }


        /// <summary>
        /// Saves position and angle of the last timestep.
        /// </summary>
        public void SavePositionAndAngleOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(m_Position);
            Aux.SaveValueOfLastTimestep(m_Angle);
        }

        /// <summary>
        /// Saves translational and rotational velocities of the last timestep.
        /// </summary>
        public void SaveVelocityOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(m_TranslationalVelocity);
            Aux.SaveValueOfLastTimestep(m_RotationalVelocity);
            Aux.SaveMultidimValueOfLastTimestep(m_TranslationalAcceleration);
            Aux.SaveValueOfLastTimestep(m_RotationalAcceleration);
        }

        /// <summary>
        /// Saves force and torque of the previous iteration.
        /// </summary>
        public void SaveHydrodynamicsOfPreviousIteration() {
            m_ForcesPreviousIteration = m_HydrodynamicForces[0].CloneAs();
            m_TorquePreviousIteration = m_HydrodynamicTorque[0];
        }

        /// <summary>
        /// Saves force and torque of the previous timestep.
        /// </summary>
        public void SaveHydrodynamicsOfPreviousTimestep() {
            Aux.SaveMultidimValueOfLastTimestep(m_HydrodynamicForces);
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
            for (int i = 0; i < historyLength; i++) {
                m_Position[i] = initialPosition.CloneAs();
                Aux.TestArithmeticException(m_Position[i], "initial particle position");
                m_Angle[i] = initialAngle * 2 * Math.PI / 360;
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
            for (int i = 0; i < historyLength; i++) {
                m_TranslationalVelocity[i] = initalTranslation == null ? (new double[spatialDim]) : initalTranslation.CloneAs();
                Aux.TestArithmeticException(m_TranslationalVelocity[i], "initial particle translational velocity");
                m_RotationalVelocity[i] = initalRotation;
                Aux.TestArithmeticException(m_RotationalVelocity[i], "initial particle rotational velocity");
            }
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
                double[] tempPos = m_Position[0].CloneAs();
                double tempAngle = m_Angle[0];
                SavePositionAndAngleOfPreviousTimestep();
                m_Position[0] = tempPos.CloneAs();
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

        /// <summary>
        /// Calls the calculation of the hydrodynamics
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="levelSetTracker"></param>
        /// <param name="fluidViscosity"></param>
        public virtual void UpdateForcesAndTorque(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, double fluidDensity, bool firstIteration, double dt = 0) {
            double[] tempForces = CalculateHydrodynamicForces(hydrodynamicsIntegration, fluidDensity);
            double tempTorque = CalculateHydrodynamicTorque(hydrodynamicsIntegration);
            HydrodynamicsPostprocessing(tempForces, tempTorque, firstIteration);
        }

        /// <summary>
        /// Predicts the hydrodynamics at the beginning of the iteration loop in each timestep.
        /// </summary>
        /// <param name="activeStress"></param>
        /// <param name="timestepID">
        /// The timestep ID. Used to distinguish between the first timestep and all other steps.
        /// </param>
        public virtual void PredictForceAndTorque(double activeStress, double circumference, int timestepID) {
            if (timestepID == 1) {
                m_HydrodynamicForces[0][0] = circumference * Math.Cos(m_Angle[0]) * activeStress / 10 + Gravity[0] * Density * ParticleArea / 10;
                m_HydrodynamicForces[0][1] = circumference * Math.Sin(m_Angle[0]) * activeStress / 10 + Gravity[1] * Density * ParticleArea / 10;
                m_HydrodynamicTorque[0] = 0;
            }
            else {
                for (int d = 0; d < spatialDim; d++) {
                    m_HydrodynamicForces[0][d] = 2 * m_HydrodynamicForces[1][d] - m_HydrodynamicForces[2][d];
                    if (Math.Abs(m_HydrodynamicForces[0][d]) < 1e-20)
                        m_HydrodynamicForces[0][d] = 0;
                }
                m_HydrodynamicTorque[0] = 2 * m_HydrodynamicTorque[1] - m_HydrodynamicTorque[2];
                if (Math.Abs(m_HydrodynamicTorque[0]) < 1e-20)
                    m_HydrodynamicTorque[0] = 0;
            }
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
        protected virtual double[] CalculateParticlePosition(double dt) {
            double[] l_Position = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_Position[d] = m_Position[1][d] + (m_TranslationalVelocity[0][d] + 4 * m_TranslationalVelocity[1][d] + m_TranslationalVelocity[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle position during the collision procedure.
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double[] CalculateParticlePosition(double dt, bool collisionProcedure) {
            double[] l_Position = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_Position[d] = m_Position[0][d] + (m_TranslationalVelocity[0][d] + 4 * m_TranslationalVelocity[1][d] + m_TranslationalVelocity[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(l_Position, "particle position");
            return l_Position;
        }

        /// <summary>
        /// Calculate the new particle position after a collision
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual double[] CalculateParticlePosition(double dt, double collisionTimestep) {
            double[] l_Position = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_Position[d] = m_Position[0][d] + m_TranslationalVelocity[0][d] * (dt - collisionTimestep) / 6;
            }
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
        protected virtual double[] CalculateTranslationalVelocity(double dt) {
            double[] l_TranslationalVelocity = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_TranslationalVelocity[d] = m_TranslationalVelocity[1][d] + (4 * m_TranslationalAcceleration[0][d] + m_TranslationalAcceleration[1][d] + m_TranslationalAcceleration[2][d]) * dt / 6;
            }
            Aux.TestArithmeticException(l_TranslationalVelocity, "particle translational velocity");
            return l_TranslationalVelocity;
        }

        /// <summary>
        /// Calculate the new translational velocity of the particle after a collision.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual double[] CalculateTranslationalVelocity(double dt, double collisionTimestep) {
            double[] l_TranslationalVelocity = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_TranslationalVelocity[d] = m_TranslationalVelocity[1][d] + m_TranslationalAcceleration[0][d] * (dt - collisionTimestep) / 6;
            }
            Aux.TestArithmeticException(l_TranslationalVelocity, "particle translational velocity");
            return l_TranslationalVelocity;
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual double CalculateAngularVelocity(double dt) {
            double l_RotationalVelocity = m_RotationalVelocity[1] + (m_RotationalAcceleration[0] + 4 * m_RotationalAcceleration[1] + m_RotationalAcceleration[2]) * dt / 6;
            Aux.TestArithmeticException(l_RotationalVelocity, "particle rotational velocity");
            return l_RotationalVelocity;
        }

        /// <summary>
        /// Calculate the new angular velocity of the particle after a collision.
        /// </summary>
        /// <param name="dt">Timestep</param>
        /// <param name="collisionTimestep">The time consumed during the collision procedure</param>
        protected virtual double CalculateAngularVelocity(double dt, double collisionTimestep) {
            double l_RotationalVelocity = m_RotationalVelocity[1] + m_RotationalAcceleration[0] * (dt - collisionTimestep) / 6;
            Aux.TestArithmeticException(l_RotationalVelocity, "particle rotational velocity");
            return l_RotationalVelocity;
        }

        /// <summary>
        /// Calculates the velocities along the normal and the tangential vector of a collision.
        /// </summary>
        internal void CalculateNormalAndTangentialVelocity() {
            double[] normalVector = m_CollisionNormalVector.Last();
            double[] TangentialVector = new double[] { -normalVector[1], normalVector[0] };
            double[] Velocity = m_TranslationalVelocity[0];
            m_PreCollisionVelocity = new double[] { Velocity[0] * normalVector[0] + Velocity[1] * normalVector[1], Velocity[0] * TangentialVector[0] + Velocity[1] * TangentialVector[1] };
            Aux.TestArithmeticException(m_PreCollisionVelocity, "particle velocity before collision");
        }

        /// <summary>
        /// Calculate the new tranlational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        protected virtual double[] CalculateTranslationalAcceleration(double dt) {
            double[] l_Acceleration = new double[spatialDim];
            for (int d = 0; d < spatialDim; d++) {
                l_Acceleration[d] = m_HydrodynamicForces[0][d] / (Density * ParticleArea);
            }
            Aux.TestArithmeticException(l_Acceleration, "particle translational acceleration");
            return l_Acceleration;
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

        /// <summary>
        /// Update Forces and Torque acting from fluid onto the particle
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        /// <param name="fluidDensity"></param>
        protected virtual double[] CalculateHydrodynamicForces(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, double fluidDensity) {
            double[] tempForces = hydrodynamicsIntegration.Forces();
            Aux.TestArithmeticException(tempForces, "temporal forces during calculation of hydrodynamics");
            Force_MPISum(ref tempForces);
            CalculateGravity(fluidDensity, tempForces);
            return tempForces;
        }

        /// <summary>
        /// Calculates the gravitational forces.
        /// </summary>
        /// <param name="fluidDensity"></param>
        /// <param name="tempForces"></param>
        private void CalculateGravity(double fluidDensity, double[] tempForces) {
            for (int d = 0; d < spatialDim; d++) {
                tempForces[d] += (Density - fluidDensity) * ParticleArea * Gravity[d];
            }
            Aux.TestArithmeticException(tempForces, "temporal forces during calculation of hydrodynamics after adding gravity");
        }

        /// <summary>
        /// Summation of the hydrodynamic forces over all MPI-processes
        /// </summary>
        /// <param name="forces"></param>
        private void Force_MPISum(ref double[] forces) {
            double[] stateBuffer = forces.CloneAs();
            double[] globalStateBuffer = stateBuffer.MPISum();
            for (int d = 0; d < spatialDim; d++) {
                forces[d] = globalStateBuffer[d];
            }
            Aux.TestArithmeticException(forces, "temporal forces during calculation of hydrodynamics after mpi-summation");
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
        protected virtual double CalculateHydrodynamicTorque(ParticleHydrodynamicsIntegration hydrodynamicsIntegration) {
            double tempTorque = hydrodynamicsIntegration.Torque(m_Position[0]);
            Aux.TestArithmeticException(tempTorque, "temporal torque during calculation of hydrodynamics");
            Torque_MPISum(ref tempTorque);
            return tempTorque;
        }

        /// <summary>
        /// Summation of the hydrodynamic torque over all MPI-processes
        /// </summary>
        /// <param name="torque"></param>
        private void Torque_MPISum(ref double torque) {
            double stateBuffer = torque;
            double globalStateBuffer = stateBuffer.MPISum();
            torque = globalStateBuffer;
            Aux.TestArithmeticException(torque, "temporal torque during calculation of hydrodynamics after mpi-summation");
        }

        /// <summary>
        /// Post-processing of the hydrodynamics. If desired the underrelaxation is applied to the forces and torque.
        /// </summary>
        /// <param name="tempForces"></param>
        /// <param name="tempTorque"></param>
        /// <param name="firstIteration"></param>
        protected void HydrodynamicsPostprocessing(double[] tempForces, double tempTorque, bool firstIteration) {
            if (m_UnderrelaxationParam != null && !firstIteration) {
                ParticleUnderrelaxation Underrelaxation = new ParticleUnderrelaxation(m_UnderrelaxationParam, CalculateAverageForces(tempForces, tempTorque));
                Underrelaxation.Forces(ref tempForces, m_ForcesPreviousIteration);
                Underrelaxation.Torque(ref tempTorque, m_TorquePreviousIteration);
            }
            for (int d = 0; d < spatialDim; d++) {
                //if(tempForces[d] > 1e-15)
                    m_HydrodynamicForces[0][d] = tempForces[d];
            }
            //if (tempTorque > 1e-15)
                m_HydrodynamicTorque[0] = tempTorque;
            Aux.TestArithmeticException(m_HydrodynamicForces[0], "hydrodynamic forces");
            Aux.TestArithmeticException(m_HydrodynamicTorque[0], "hydrodynamic torque");
        }

        /// <summary>
        /// Does what it says.
        /// </summary>
        /// <param name="forces">
        /// The hydrodynamic forces.
        /// </param>
        /// <param name="torque">
        /// The hydrodynamic torque.
        /// </param>
        /// <param name="averageDistance">
        /// The average Lengthscale of the particle.
        /// </param>
        private double CalculateAverageForces(double[] forces, double torque) {
            double averageForces = Math.Abs(torque) / MaxParticleLengthScale;
            for (int d = 0; d < forces.Length; d++) {
                averageForces += forces[d];
            }
            return averageForces /= 3;
        }

        /// <summary>
        /// Calculating the particle Reynolds number
        /// </summary>
        /// <param name="fluidViscosity"></param>
        public double ComputeParticleRe(double fluidViscosity) {
            return m_TranslationalVelocity[0].L2Norm() * MaxParticleLengthScale / fluidViscosity;
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
            double[] temp = new double[spatialDim + 1];
            for (int d = 0; d < spatialDim; d++) {
                temp[d] = Mass_P * m_TranslationalVelocity[0][d];
            }
            temp[spatialDim] = MomentOfInertia * m_RotationalVelocity[0];
            return temp;
        }

        /// <summary>
        /// Calculating the particle kinetic energy
        /// </summary>
        public double[] CalculateParticleKineticEnergy() {
            double[] temp = new double[spatialDim + 1];
            for (int d = 0; d < spatialDim; d++) {
                temp[d] = 0.5 * Mass_P * m_TranslationalVelocity[0][d].Pow2();
            }
            temp[spatialDim] = 0.5 * MomentOfInertia * m_RotationalVelocity[0].Pow2();
            return temp;
        }

        /// <summary>
        /// Deletes the complete history of the translational velocity and acceleration.
        /// </summary>
        public void ClearParticleHistoryTranslation() {
            for (int i = 0; i < m_TranslationalVelocity.Count; i++) {
                for (int d = 0; d < spatialDim; d++) {
                    m_TranslationalVelocity[i][d] = 0;
                    m_TranslationalAcceleration[i][d] = 0;
                }
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
        public void SetCollisionVectors(double[] normalVector, double[] tangentialVector) {
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
            m_CollisionTranslationalVelocity.Add(new double[] { normalVelocity, tangentialVelocity });
            m_CollisionRotationalVelocity.Add(rotationalVelocity);
        }

        /// <summary>
        /// Collision post-processing. Sums up the results of the multiple binary collisions of one timestep
        /// </summary>
        public void PostProcessCollisionTranslation() {
            if (m_CollisionTranslationalVelocity.Count >= 1) {
                double[] normal = new double[spatialDim];
                double[] tangential = new double[spatialDim];
                for (int t = 0; t < m_CollisionTranslationalVelocity.Count; t++) {
                    for (int d = 0; d < spatialDim; d++) {
                        normal[d] += m_CollisionNormalVector[t][d];
                        tangential[d] += m_CollisionTangentialVector[t][d];
                    }
                }

                normal.ScaleV(1 / Math.Sqrt(normal[0].Pow2() + normal[1].Pow2()));
                tangential.ScaleV(1 / Math.Sqrt(tangential[0].Pow2() + tangential[1].Pow2()));
                double temp_NormalVel = 0;
                double temp_TangentialVel = 0;
                for (int t = 0; t < m_CollisionTranslationalVelocity.Count; t++) {
                    double cos = new double();
                    for (int d = 0; d < spatialDim; d++) {
                        cos += normal[d] * m_CollisionNormalVector[t][d];
                    }
                    double sin = cos == 1 ? 0 : m_CollisionNormalVector[t][0] > normal[0] ? Math.Sqrt(1 + 1e-15 - cos.Pow2()) : -Math.Sqrt(1 + 1e-15 - cos.Pow2());
                    temp_NormalVel += m_CollisionTranslationalVelocity[t][0] * cos - m_CollisionTranslationalVelocity[t][1] * sin;
                    temp_TangentialVel += m_CollisionTranslationalVelocity[t][0] * sin + m_CollisionTranslationalVelocity[t][1] * cos;

                }
                temp_NormalVel /= m_CollisionTranslationalVelocity.Count;
                temp_TangentialVel /= m_CollisionTranslationalVelocity.Count;

                ClearParticleHistoryTranslation();
                for (int d = 0; d < spatialDim; d++) {
                    m_TranslationalVelocity[0][d] = normal[d] * temp_NormalVel + tangential[d] * temp_TangentialVel;
                }
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
            double[] dataSend = new double[13];
            dataSend[0] = m_RotationalVelocity[0];
            dataSend[1] = m_TranslationalVelocity[0][0];
            dataSend[2] = m_TranslationalVelocity[0][1];
            dataSend[3] = m_Angle[0];
            dataSend[4] = m_Position[0][0];
            dataSend[5] = m_Position[0][1];
            dataSend[6] = m_CollisionTimestep;
            dataSend[7] = m_RotationalVelocity[1];
            dataSend[8] = m_TranslationalVelocity[1][0];
            dataSend[9] = m_TranslationalVelocity[1][1];
            dataSend[10] = m_Angle[1];
            dataSend[11] = m_Position[1][0];
            dataSend[12] = m_Position[1][1];
            return dataSend;
        }

        /// <summary>
        /// Overwrites the particles parameters with the values received during the post-collision MPI-communication.
        /// </summary>
        public void WriteReceiveArray(double[] dataReceive, int offset) {
            m_RotationalVelocity[0] = dataReceive[0 + offset];
            m_TranslationalVelocity[0][0] = dataReceive[1 + offset];
            m_TranslationalVelocity[0][1] = dataReceive[2 + offset];
            m_Angle[0] = dataReceive[3 + offset];
            m_Position[0][0] = dataReceive[4 + offset];
            m_Position[0][1] = dataReceive[5 + offset];
            m_CollisionTimestep = dataReceive[6 + offset];
            m_RotationalVelocity[1] = dataReceive[7 + offset];
            m_TranslationalVelocity[1][0] = dataReceive[8 + offset];
            m_TranslationalVelocity[1][1] = dataReceive[9 + offset];
            m_Angle[1] = dataReceive[10 + offset];
            m_Position[1][0] = dataReceive[11 + offset];
            m_Position[1][1] = dataReceive[12 + offset];
        }
    }
}

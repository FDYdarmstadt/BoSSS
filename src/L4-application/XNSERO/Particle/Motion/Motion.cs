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
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSERO_Solver {
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
        public Motion(double density) {
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
            PeriodicBoundaryPosition = new double[SpatialDim][];
        }

        /// <summary>
        /// First index: The spatial direction of the periodic boundary.
        /// Second index: The ID of the periodic boundary. Left and lower boundary has ID=0, right and upper boundary ID=1.
        /// </summary>
        [DataMember]
        private readonly double[][] PeriodicBoundaryPosition;

        /// <summary>
        /// The origin of the virtual domain at the periodic boundary.
        /// </summary>
        [DataMember]
        public List<Vector> OriginInVirtualPeriodicDomain = new List<Vector>();

        /// <summary>
        /// This method provides information about periodic boundaries for the particle. 
        /// </summary>
        /// <param name="periodicBoundaryPosition">Relative periodic boundary position towards the origin.</param>
        /// <param name="dimension"></param>
        /// <remarks>At each periodic boundary a virtual domain is created with its own origin, which coordinates are given with respect to the main origin.
        /// Considering all edges and vertices of a rectangular domain this leads to eight additional virtual domains.
        /// </remarks>
        public void SetPeriodicBoundary(double[] periodicBoundaryPosition, int dimension) {
            Aux = new Auxillary();
            Aux.TestArithmeticException(periodicBoundaryPosition, ("periodic boundary position, dimension " + dimension));
            PeriodicBoundaryPosition[dimension] = periodicBoundaryPosition.CloneAs();

            switch (dimension) {
                case 0:
                xPeriodic = true;
                for (int d = 0; d < PeriodicBoundaryPosition[dimension].Length; d++) {
                    OriginInVirtualPeriodicDomain.Add(new Vector(2 * PeriodicBoundaryPosition[0][d], 0));
                }
                break;
                case 1:
                yPeriodic = true;
                for (int d = 0; d < PeriodicBoundaryPosition[dimension].Length; d++) {
                    OriginInVirtualPeriodicDomain.Add(new Vector(0, 2 * PeriodicBoundaryPosition[1][d]));
                }
                break;
            }

            if(xPeriodic && yPeriodic) {
                for(int d1 = 0; d1 < 2; d1++) {
                    for(int d2 = 0; d2 < 2; d2++) {
                        OriginInVirtualPeriodicDomain.Add(new Vector(OriginInVirtualPeriodicDomain[d1][0], OriginInVirtualPeriodicDomain[2 + d2][1]));
                    }
                }
            }
        }

        /// <summary>
        /// Periodic boundary in x-direction?
        /// </summary>
        [DataMember]
        private bool xPeriodic = false;

        /// <summary>
        /// Periodic boundary in y-direction?
        /// </summary>
        [DataMember]
        private bool yPeriodic = false;

        /// <summary>
        /// Checks whether a point is inside of the domain or outside in the virtual domain at a periodic boundary.
        /// </summary>
        /// <param name="Point"></param>
        /// <param name="Tolerance"></param>
        /// <returns></returns>
        public bool IsInsideOfPeriodicDomain(Vector Point, double Tolerance) {
            for (int d = 0; d < PeriodicBoundaryPosition.Length; d++) {
                if (!PeriodicBoundaryPosition[d].IsNullOrEmpty()) {
                    Vector wallNormal = new Vector(1 - d, d);
                    for (int wallID = 0; wallID < PeriodicBoundaryPosition[d].Length; wallID++) {
                        double tolerance = Tolerance;
                        if (wallID == 1) {
                            tolerance = Tolerance * (-1);
                            wallNormal *= -1;
                        }
                        double wallWithTolerance = PeriodicBoundaryPosition[d][wallID] - tolerance;
                        Vector wallToPoint = d == 0 ? new Vector(Point[0] - (wallWithTolerance), Point[1]) : new Vector(Point[0], Point[1] - (wallWithTolerance));
                        if (wallNormal * wallToPoint < 0)
                            return false;
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// Parameter used for Aitken-Relaxation.
        /// </summary>
        internal double omega = 1;

        /// <summary>
        /// Provides additional methods for testing.
        /// </summary>
        [NonSerialized]
        internal Auxillary Aux = new Auxillary();

        /// <summary>
        /// The number of time-steps in the history of the particle to be safed.
        /// </summary>
        [DataMember]
        private const int NumberOfHistoryEntries = 3;

        /// <summary>
        /// Dimension.
        /// </summary>
        [DataMember]
        protected static int SpatialDim = 2;

        /// <summary>
        /// The history of the position.
        /// </summary>
        [DataMember]
        private readonly List<Vector> Position = new List<Vector>();

        /// <summary>
        /// The history of the Angle.
        /// </summary>
        [DataMember]
        private readonly List<double> Angle = new List<double>();

        /// <summary>
        /// The history of the translational velocity.
        /// </summary>
        [DataMember]
        private readonly List<Vector> TranslationalVelocity = new List<Vector>();

        /// <summary>
        /// The history of the rotational velocity.
        /// </summary>
        [DataMember]
        private readonly List<double> RotationalVelocity = new List<double>();

        /// <summary>
        /// The history of the translational acceleration.
        /// </summary>
        [DataMember]
        private readonly List<Vector> TranslationalAcceleration = new List<Vector>();

        /// <summary>
        /// The history of the rotational acceleration.
        /// </summary>
        [DataMember]
        private readonly List<double> RotationalAcceleration = new List<double>();

        /// <summary>
        /// The history of the hydrodynamic forces.
        /// </summary>
        [DataMember]
        private readonly List<Vector> HydrodynamicForces = new List<Vector>();

        /// <summary>
        /// The history of the hydrodynamic torque.
        /// </summary>
        [DataMember]
        private readonly List<double> HydrodynamicTorque = new List<double>();

        /// <summary>
        /// During the collision algorithm the particle is moved with the current velocity using "safe" time-steps without the collision.
        /// Hence, for a position update after the collision on need to subtract the already used time from the overall time-step
        /// </summary>
        [DataMember]
        private double CollisionTimestep = 0;

        /// <summary>
        /// Density
        /// </summary>
        [DataMember]
        public readonly double Density;

        /// <summary>
        /// Particle volume
        /// </summary>
        [DataMember]
        public double Volume;

        /// <summary>
        /// Moment of inertia (scalar in 2D)
        /// </summary>
        [DataMember]
        public double MomentOfInertia;

        /// <summary>
        /// The maximum length of the particle.
        /// </summary>
        [DataMember]
        public double MaxLength;

        /// <summary>
        /// If used, the added damping tensor is saved here.
        /// </summary>
        [DataMember]
        internal double[,] AddedDampingTensor = new double[6, 6];

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        [DataMember]
        public double Mass => Volume * Density;
        
        /// <summary>
        /// Include rotation?
        /// </summary>
        [DataMember]
        internal bool IncludeRotation = true;

        /// <summary>
        /// Include translation?
        /// </summary>
        [DataMember]
        internal bool IncludeTranslation = true;

        /// <summary>
        /// Use added damping?, for reference: Banks et.al. 2017
        /// </summary>
        [DataMember]
        internal virtual bool UseAddedDamping { get; } = false;

        /// <summary>
        /// Returns the position of the particle.
        /// </summary>
        /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetPosition(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return Position[historyPosition];
        }

        /// <summary>
        /// Returns the angle of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetAngle(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return Angle[historyPosition];
        }

        /// <summary>
        /// Returns the translational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetTranslationalVelocity(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return TranslationalVelocity[historyPosition];
        }

        /// <summary>
        /// Returns the rotational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetRotationalVelocity(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return RotationalVelocity[historyPosition];
        }

        /// <summary>
        /// Returns the translational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetTranslationalAcceleration(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return TranslationalAcceleration[historyPosition];
        }

        /// <summary>
        /// Returns the rotational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetRotationalAcceleration(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return RotationalAcceleration[historyPosition];
        }

        /// <summary>
        /// Returns the force acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for 4 time-steps. historyPosition=0 returns the newest value.
        /// </param>
        internal Vector GetHydrodynamicForces(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return HydrodynamicForces[historyPosition];
        }

        /// <summary>
        /// Returns the torque acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        internal double GetHydrodynamicTorque(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return HydrodynamicTorque[historyPosition];
        }

        /// <summary>
        /// Saves position and angle of the last time-step.
        /// </summary>
        public void SavePositionAndAngleOfPreviousTimestep() {
            using (new FuncTrace()) {
                Aux.SaveVectorOfLastTimestep(Position);
                Aux.SaveValueOfLastTimestep(Angle);
            }
        }

        /// <summary>
        /// Saves translational and rotational velocities of the last time-step.
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
        /// Saves force and torque of the previous time-step.
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
                Aux = new Auxillary();
                for (int i = 0; i < historyLength; i++) {
                    Position[i] = new Vector(initialPosition);
                    Angle[i] = initialAngle * 2 * Math.PI / 360;
                    Aux.TestArithmeticException(Position[i], "initial particle position");
                    Aux.TestArithmeticException(Angle[i], "initial particle angle");
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
        /// Init of the particle area.
        /// </summary>
        /// <param name="Volume"></param>
        public void SetVolume(double Volume) => this.Volume = Volume;

        /// <summary>
        /// Init of the lengthscale.
        /// </summary>
        /// <param name="moment"></param>
        public void SetMaxLength(double lengthscale) => MaxLength = lengthscale;

        /// <summary>
        /// Init of the moment of inertia.
        /// </summary>
        /// <param name="moment"></param>
        public void SetMomentOfInertia(double moment) => MomentOfInertia = moment;

        /// <summary>
        /// Sets the collision time-step.
        /// </summary>
        /// <param name="collisionTimestep">
        /// The physical time consumend by the collision procedure.
        /// </param>
        public void SetCollisionTimestep(double collisionTimestep) => CollisionTimestep = collisionTimestep;


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
                    throw new Exception("Collision time-step: " + CollisionTimestep);
                }
                Position[0] = CalculateParticlePosition(dt - CollisionTimestep);
                Angle[0] = CalculateParticleAngle(dt - CollisionTimestep);
                if (Angle[0] > 2 * Math.PI)
                    Angle[0] -= 2 * Math.PI;
                CollisionTimestep = 0;
                if(!IsInsideOfPeriodicDomain(Position[0], 0)) {
                    for(int i = 0; i < OriginInVirtualPeriodicDomain.Count(); i++) {
                        if(IsInsideOfPeriodicDomain(Position[0] + OriginInVirtualPeriodicDomain[i], 0)) {
                            for(int h = 0; h < Position.Count(); h++) {
                                Position[h] = new Vector(Position[h] + OriginInVirtualPeriodicDomain[i]);
                            }
                            break;
                        }
                    }
                }
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

        /// <summary>
        /// Sets the newly calculated hydrodynamic to this particle
        /// </summary>
        /// <param name="particleID">This particle ID</param>
        /// <param name="AllParticleHydrodynamics">Hydrodynamics of all particles.</param>
        internal void UpdateForcesAndTorque(int particleID, double[] AllParticleHydrodynamics) {
            using (new FuncTrace()) {
                Vector forces = new Vector(SpatialDim);
                for (int d = 0; d < SpatialDim; d++) {
                    if (Math.Abs(AllParticleHydrodynamics[particleID * 3 + d]) > 1e-12)
                        forces[d] = AllParticleHydrodynamics[particleID * 3 + d];

                }
                HydrodynamicForces[0] = forces;
                if (Math.Abs(AllParticleHydrodynamics[particleID * 3 + SpatialDim]) > 1e-12)
                    HydrodynamicTorque[0] = AllParticleHydrodynamics[particleID * 3 + SpatialDim];
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
        /// <param name="dt">Time-step</param>
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
        /// <param name="dt">Time-step</param>
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
                Vector l_Acceleration = HydrodynamicForces[0] / (Density * Volume);
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
        /// Calls the integration of the hydrodynamic stress at this particles level-set
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        public virtual Vector CalculateHydrodynamics(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, string[] FluidSpecies, double dt = 0) {
            if(cutCells.IsEmptyOnRank)// no cutCells on this process.
                return new Vector(SpatialDim + 1);
            int NoOfFluidSpecies = FluidSpecies.Length;
            Vector forcesAndTorque = new Vector(SpatialDim + 1);
            for(int i = 0; i < NoOfFluidSpecies; i++) {
                forcesAndTorque += hydrodynamicsIntegration.Main(Position[0], cutCells, FluidSpecies[i]);
            }
            Aux.TestArithmeticException(forcesAndTorque, "during calculation of hydrodynamics");
            return forcesAndTorque;
        }

        /// <summary>
        /// Calculates the gravitational forces.
        /// </summary>
        /// <param name="fluidDensity"></param>
        /// <param name="tempForces"></param>
        public Vector GetGravityForces(Vector Gravity) {
            return Density * Volume * Gravity;
        }

        /// <summary>
        /// Calculating the particle Reynolds number
        /// </summary>
        /// <param name="fluidViscosity"></param>
        public double ComputeParticleReynoldsNumber(double fluidViscosity) {
            using (new FuncTrace()) {
                return 2 * TranslationalVelocity[0].L2Norm() * MaxLength / fluidViscosity;
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
                    temp[d] = Mass * TranslationalVelocity[0][d];
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
                    temp[d] = 0.5 * Mass * TranslationalVelocity[0][d].Pow2();
                }
                temp[SpatialDim] = 0.5 * MomentOfInertia * RotationalVelocity[0].Pow2();
                return temp;
            }
        }

        public virtual object Clone() {
            Motion clonedMotion = new Motion(Density);
            clonedMotion.SetVolume(Volume);
            clonedMotion.SetMomentOfInertia(MomentOfInertia);
            return clonedMotion;
        }
    }
}

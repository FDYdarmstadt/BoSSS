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
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSERO_Solver {
    /// <summary>
    /// This class provides basic functionality of any motion model, e.g. saving of previous time-steps, information about periodic boundaries etc.
    /// It does not contain the calculation of the physical variables, which have to be implemented by derived classes.
    /// </summary>
    public abstract class MotionBase : IMotion{
        /// <summary>
        /// Ensures basic functionality of any motion model.
        /// </summary>
        /// <param name="density">
        /// The density of the particle.
        /// </param>
        public MotionBase(double density) {
            this.density = density;
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

        [DataMember]
        public List<Vector> originInVirtualPeriodicDomain = new();
        /// <summary>
        /// The origin of the virtual domain at the periodic boundary.
        /// </summary>
        public List<Vector> OriginInVirtualPeriodicDomain {
            get { return originInVirtualPeriodicDomain; }
        }

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
                    originInVirtualPeriodicDomain.Add(new Vector(2 * PeriodicBoundaryPosition[0][d], 0));
                }
                break;
                case 1:
                yPeriodic = true;
                for (int d = 0; d < PeriodicBoundaryPosition[dimension].Length; d++) {
                    originInVirtualPeriodicDomain.Add(new Vector(0, 2 * PeriodicBoundaryPosition[1][d]));
                }
                break;
                default:
                throw new NotImplementedException("Periodic boundaries for particles only in 2D");
            }

            if (xPeriodic && yPeriodic) {
                for (int d1 = 0; d1 < 2; d1++) {
                    for (int d2 = 0; d2 < 2; d2++) {
                        originInVirtualPeriodicDomain.Add(new Vector(originInVirtualPeriodicDomain[d1][0], originInVirtualPeriodicDomain[2 + d2][1]));
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
                    Vector wallNormal = new(1 - d, d);
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
        private double relaxationParameter = 100;

        /// <summary>
        /// Parameter used for Aitken-Relaxation.
        /// </summary>
        public double RelaxationParameter {
            get { return relaxationParameter; }
            set { relaxationParameter = value; }
        }

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
        protected readonly List<Vector> Position = new();

        /// <summary>
        /// The history of the Angle.
        /// </summary>
        [DataMember]
        protected readonly List<double> Angle = new();

        /// <summary>
        /// The history of the translational velocity.
        /// </summary>
        [DataMember]
        protected readonly List<Vector> TranslationalVelocity = new();

        /// <summary>
        /// The history of the rotational velocity.
        /// </summary>
        [DataMember]
        protected readonly List<double> RotationalVelocity = new();

        /// <summary>
        /// The history of the translational acceleration.
        /// </summary>
        [DataMember]
        protected readonly List<Vector> TranslationalAcceleration = new();

        /// <summary>
        /// The history of the rotational acceleration.
        /// </summary>
        [DataMember]
        protected readonly List<double> RotationalAcceleration = new();

        /// <summary>
        /// The history of the hydrodynamic forces.
        /// </summary>
        [DataMember]
        protected readonly List<Vector> HydrodynamicForces = new();

        /// <summary>
        /// The history of the hydrodynamic torque.
        /// </summary>
        [DataMember]
        protected readonly List<double> HydrodynamicTorque = new();

        
        [DataMember]
        private double collisionTimestep = 0;
        /// <summary>
        /// During the collision algorithm the particle is moved with the current velocity using "safe" time-steps without the collision.
        /// Hence, for a position update after the collision on need to subtract the already used time from the overall time-step
        /// </summary>
        public double CollisionTimestep {
            get { return collisionTimestep; }
            set { collisionTimestep = value; }
        }

        /// <summary>
        /// Density
        /// </summary>
        [DataMember]
        private readonly double density;

        public double Density {
            get { return density; }
        }

        /// <summary>
        /// Particle volume
        /// </summary>
        [DataMember]
        private double volume;

        /// <summary>
        /// Particle volume
        /// </summary>
        public double Volume {
            get { return volume; }
            set { volume = value; }
        }

        /// <summary>
        /// Moment of inertia (scalar in 2D)
        /// </summary>
        [DataMember]
        private double momentOfInertia;

        public double MomentOfInertia {
            get { return momentOfInertia; }
            set { momentOfInertia = value; }
        }

        [DataMember]
        private double characteristicLength;

        /// <summary>
        /// The maximum length of the particle.
        /// </summary>
        public double CharacteristicLength { 
            get { return characteristicLength; }
            set { characteristicLength = value; }
        }

        /// <summary>
        /// Mass of the current particle.
        /// </summary>
        [DataMember]
        public double Mass => volume * density;

        public abstract bool IncludeRotation();

        public abstract bool IncludeTranslation();

        /// <summary>
        /// Returns the position of the particle.
        /// </summary>
        /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public Vector GetPosition(int historyPosition = 0) {
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
        public double GetAngle(int historyPosition = 0) {
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
        public Vector GetTranslationalVelocity(int historyPosition = 0) {
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
        public double GetRotationalVelocity(int historyPosition = 0) {
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
        public Vector GetTranslationalAcceleration(int historyPosition = 0) {
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
        public double GetRotationalAcceleration(int historyPosition = 0) {
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
        public Vector GetHydrodynamicForces(int historyPosition = 0) {
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
        public double GetHydrodynamicTorque(int historyPosition = 0) {
            if (historyPosition >= NumberOfHistoryEntries)
                throw new Exception("Error in Particle.Motion: Only " + NumberOfHistoryEntries + " time-steps are saved. The requested value is " + historyPosition + " steps in the past!");
            return HydrodynamicTorque[historyPosition];
        }

        /// <summary>
        /// Saves position and angle of the last time-step.
        /// </summary>
        public void SavePositionAndAngleOfPreviousTimestep() {
            using (new FuncTrace()) {
                SaveVectorOfLastTimestep(Position);
                SaveValueOfLastTimestep(Angle);
            }
        }

        /// <summary>
        /// Saves translational and rotational velocities of the last time-step.
        /// </summary>
        public void SaveVelocityOfPreviousTimestep() {
            using (new FuncTrace()) {
                SaveVectorOfLastTimestep(TranslationalVelocity);
                SaveValueOfLastTimestep(RotationalVelocity);
                SaveVectorOfLastTimestep(TranslationalAcceleration);
                SaveValueOfLastTimestep(RotationalAcceleration);
            }
        }

        /// <summary>
        /// Saves force and torque of the previous time-step.
        /// </summary>
        public void SaveHydrodynamicsOfPreviousTimestep() {
            SaveVectorOfLastTimestep(HydrodynamicForces);
            for (int i = 0; i < NumberOfHistoryEntries; i++)
                SaveValueOfLastTimestep(HydrodynamicTorque);
        }

        private static void SaveValueOfLastTimestep(List<double> variable) {
            variable.Insert(0, new double());
            variable[0] = 0;
            variable.RemoveAt(variable.Count - 1);
        }

        private static void SaveVectorOfLastTimestep(List<Vector> variable) {
            int dim = variable[0].Dim;
            variable.Insert(0, new Vector(dim));
            variable.RemoveAt(variable.Count - 1);
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
        public void InitializeParticlePositionAndAngle(double[] initialPosition, double initialAngle, int historyLength = 0, int currentHistoryPos = 0) {
            using (new FuncTrace()) {
                if (historyLength == 0)
                    historyLength = NumberOfHistoryEntries;
                Aux = new Auxillary();
                for (int i = currentHistoryPos; i < historyLength; i++) {
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
        public void InitializeParticleVelocity(double[] initalTranslation, double initalRotation, int historyLength = 0, int currentHistoryPos = 0) {
            using (new FuncTrace()) {
                if (historyLength == 0)
                    historyLength = NumberOfHistoryEntries;
                for (int i = currentHistoryPos; i < historyLength; i++) {
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
        public void InitializeParticleAcceleration(double[] initalTranslationAcceleration, double initalRotationAcceleration, int historyLength = 0, int currentHistoryPos = 0) {
            using (new FuncTrace()) {
                if (historyLength == 0)
                    historyLength = NumberOfHistoryEntries;
                for (int i = currentHistoryPos; i < historyLength; i++) {
                    TranslationalAcceleration[i] = new Vector(initalTranslationAcceleration);
                    RotationalAcceleration[i] = initalRotationAcceleration;
                    Aux.TestArithmeticException(TranslationalVelocity[i], "initial particle translational velocity");
                    Aux.TestArithmeticException(RotationalVelocity[i], "initial particle rotational velocity");
                }
            }
        }

        /// <summary>
        /// Used for testing
        /// </summary>
        /// <param name="initalTranslation">
        /// The initial translational velocity.
        /// </param>
        /// <param name="initalRotation">
        /// The initial rotational velocity.
        /// </param>
        public void PrescribeHydrodynamicForcesAndTorque(Vector Forces, double Torque, int HistoryPosition) {
            HydrodynamicForces[HistoryPosition] = Forces;
            HydrodynamicTorque[HistoryPosition] = Torque;
        }

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
                if (!IsInsideOfPeriodicDomain(Position[0], 0)) {
                    for (int i = 0; i < originInVirtualPeriodicDomain.Count; i++) {
                        if (IsInsideOfPeriodicDomain(Position[0] + originInVirtualPeriodicDomain[i], 0)) {
                            for (int h = 0; h < Position.Count(); h++) {
                                Position[h] = new Vector(Position[h] + originInVirtualPeriodicDomain[i]);
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
        public void MoveParticleDuringCollision(double collisionDynamicTimestep) {
            using (new FuncTrace()) {
                Position[0] = CalculateParticlePositionDuringCollision(collisionDynamicTimestep);
                Angle[0] = CalculateParticleAngleDuringCollision(collisionDynamicTimestep);
            }
        }

        /// <summary>
        /// Calls the calculation of the velocity.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticleVelocity(double dt) {
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
        public void UpdateForcesAndTorque(int particleID, double[] AllParticleHydrodynamics) {
            using (new FuncTrace()) {
                Vector forces = new(SpatialDim);
                for (int d = 0; d < SpatialDim; d++) {
                    forces[d] = AllParticleHydrodynamics[particleID * 3 + d];
                }
                HydrodynamicForces[0] = forces;
                HydrodynamicTorque[0] = AllParticleHydrodynamics[particleID * 3 + SpatialDim];
                Aux.TestArithmeticException(HydrodynamicForces[0], "hydrodynamic forces");
                Aux.TestArithmeticException(HydrodynamicTorque[0], "hydrodynamic torque");
            }
        }

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public abstract Vector CalculateParticlePosition(double dt);

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public abstract double CalculateParticleAngle(double dt);

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public abstract Vector CalculateParticlePositionDuringCollision(double dt);

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public abstract double CalculateParticleAngleDuringCollision(double dt);

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Time-step</param>
        public abstract Vector CalculateTranslationalVelocity(double dt);

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Time-step</param>
        public abstract double CalculateAngularVelocity(double dt);

        /// <summary>
        /// Calculate the new tranlational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        public abstract Vector CalculateTranslationalAcceleration(double dt);

        /// <summary>
        /// Calculate the new rotational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        public abstract double CalculateRotationalAcceleration(double dt);

        /// <summary>
        /// Calls the integration of the hydrodynamic stress at this particles level-set
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        public abstract Vector CalculateHydrodynamics(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, string[] FluidSpecies, double dt = 0);

        /// <summary>
        /// Calculates the gravitational forces.
        /// </summary>
        /// <param name="fluidDensity"></param>
        /// <param name="tempForces"></param>
        public virtual Vector GravityForce(Vector Gravity) {
            return density * volume * Gravity;
        }

        /// <summary>
        /// Calculating the particle Reynolds number
        /// </summary>
        /// <param name="fluidViscosity"></param>
        public double ComputeParticleReynoldsNumber(double fluidViscosity) {
            return 2 * TranslationalVelocity[0].L2Norm() * CharacteristicLength / fluidViscosity;
        }

        /// <summary>
        /// Calculating the particle Stokes number
        /// </summary>
        /// <param name="fluidViscosity"></param>
        /// <param name="fluidDensity"></param>
        public double ComputeParticleStokesNumber(double fluidViscosity, double fluidDensity) {
            return ComputeParticleReynoldsNumber(fluidViscosity) * density / (9 * fluidDensity);
        }

        /// <summary>
        /// Calculating the particle momentum
        /// </summary>
        public double[] CalculateParticleMomentum() {
            double[] temp = new double[SpatialDim + 1];
            for (int d = 0; d < SpatialDim; d++) {
                temp[d] = Mass * TranslationalVelocity[0][d];
            }
            temp[SpatialDim] = momentOfInertia * RotationalVelocity[0];
            return temp;
        }

        /// <summary>
        /// Calculating the particle kinetic energy
        /// </summary>
        public double[] CalculateParticleKineticEnergy() {
            double[] temp = new double[SpatialDim + 1];
            for (int d = 0; d < SpatialDim; d++) {
                temp[d] = 0.5 * Mass * TranslationalVelocity[0][d].Pow2();
            }
            temp[SpatialDim] = 0.5 * momentOfInertia * RotationalVelocity[0].Pow2();
            return temp;
        }

        public virtual object Clone() => throw new NotImplementedException();
    }
}

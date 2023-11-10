/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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
using System;
using System.Collections.Generic;

namespace BoSSS.Application.XNSERO_Solver {
    public interface IMotion : ICloneable {

        /// <summary>
        /// This method provides information about periodic boundaries for the particle. 
        /// </summary>
        /// <param name="periodicBoundaryPosition">Relative periodic boundary position towards the origin.</param>
        /// <param name="dimension"></param>
        /// <remarks>At each periodic boundary a virtual domain is created with its own origin, which coordinates are given with respect to the main origin.
        /// Considering all edges and vertices of a rectangular domain in two dimensions this leads to eight additional virtual domains.
        /// </remarks>
        public void SetPeriodicBoundary(double[] periodicBoundaryPosition, int dimension);

        /// <summary>
        /// Flag for translational motion
        /// </summary>
        /// <returns></returns>
        public bool IncludeTranslation();

        /// <summary>
        /// Flag for rotational motion
        /// </summary>
        /// <returns></returns>
        public bool IncludeRotation();

        /// <summary>
        /// The origin of the virtual domain at the periodic boundary.
        /// </summary>
        public List<Vector> OriginInVirtualPeriodicDomain { get; }

        /// <summary>
        /// Checks whether a point is inside of the domain or outside in the virtual domain at a periodic boundary.
        /// </summary>
        /// <param name="Point"></param>
        /// <param name="Tolerance"></param>
        /// <returns></returns>
        public bool IsInsideOfPeriodicDomain(Vector Point, double Tolerance);

        /// <summary>
        /// Returns the position of the particle.
        /// </summary>
        /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public Vector GetPosition(int historyPosition = 0);

        /// <summary>
        /// Returns the angle of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public double GetAngle(int historyPosition = 0);

        /// <summary>
        /// Returns the translational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public Vector GetTranslationalVelocity(int historyPosition = 0);

        /// <summary>
        /// Returns the rotational velocity of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public double GetRotationalVelocity(int historyPosition = 0);

        /// <summary>
        /// Returns the translational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public Vector GetTranslationalAcceleration(int historyPosition = 0);

        /// <summary>
        /// Returns the rotational acceleration of the particle.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public double GetRotationalAcceleration(int historyPosition = 0);

        /// <summary>
        /// Returns the force acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for 4 time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public Vector GetHydrodynamicForces(int historyPosition = 0);

        /// <summary>
        /// Returns the torque acting on the particle in the current time step.
        /// </summary>
        /// /// <param name="historyPosition">
        /// The history of the particle is saved for four time-steps. historyPosition=0 returns the newest value.
        /// </param>
        public double GetHydrodynamicTorque(int historyPosition = 0);

        /// <summary>
        /// Saves position and angle of the last time-step.
        /// </summary>
        public void SavePositionAndAngleOfPreviousTimestep();

        /// <summary>
        /// Saves translational and rotational velocities of the last time-step.
        /// </summary>
        public void SaveVelocityOfPreviousTimestep();

        /// <summary>
        /// Used during initialization of the particle. Sets the position and the angle.
        /// </summary>
        /// <param name="initialPosition">
        /// The initial position.
        /// </param>
        /// <param name="initialAngle">
        /// The initial angle.
        /// </param>
        public void InitializeParticlePositionAndAngle(double[] initialPosition, double initialAngle, int historyLength = 0, int currentHistoryPos = 0);

        /// <summary>
        /// Used during initialization of the particle. Sets the translational and rotational velocity.
        /// </summary>
        /// <param name="initalTranslation">
        /// The initial translational velocity.
        /// </param>
        /// <param name="initalRotation">
        /// The initial rotational velocity.
        /// </param>
        public void InitializeParticleVelocity(double[] initalTranslation, double initalRotation, int historyLength = 0, int currentHistoryPos = 0);

        /// <summary>
        /// Used during initialization of the particle. Sets the translational and rotational acceleration.
        /// </summary>
        /// <param name="initalTranslationAcceleration"></param>
        /// <param name="initalRotationAcceleration"></param>
        /// <param name="historyLength"></param>
        /// <param name="currentHistoryPos"></param>
        public void InitializeParticleAcceleration(double[] initalTranslationAcceleration, double initalRotationAcceleration, int historyLength = 0, int currentHistoryPos = 0);

        /// <summary>
        /// Volume (3D) or area (2D) of a particle
        /// </summary>
        /// <param name="Volume"></param>
        public double Volume { get; set; }

        /// <summary>
        /// Particle density
        /// </summary>
        public double Density { get; }

        /// <summary>
        /// The characteristic length of the particle.
        /// </summary>
        public double CharacteristicLength { get; set; }

        /// <summary>
        /// Particle moment of inertia.
        /// </summary>
        public double MomentOfInertia { get; set; }

        /// <summary>
        /// Parameter used for Aitken-Relaxation.
        /// </summary>
        public double RelaxationParameter { get; set; }

        /// <summary>
        /// During the collision algorithm the particle is moved with the current velocity using "safe" time-steps without the collision.
        /// Hence, for a position update after the collision on need to subtract the already used time from the overall time-step
        /// </summary>
        public double CollisionTimestep { get; set; }

        /// <summary>
        /// Calls the calculation of the position and angle.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticlePositionAndAngle(double dt);

        /// <summary>
        /// Calls the calculation of the position and angle during the calculation of the collisions.
        /// </summary>
        /// <param name="dt"></param>
        public void MoveParticleDuringCollision(double collisionDynamicTimestep);

        /// <summary>
        /// Calls the calculation of the velocity.
        /// </summary>
        /// <param name="dt"></param>
        public void UpdateParticleVelocity(double dt);

        /// <summary>
        /// Sets the newly calculated hydrodynamic to this particle
        /// </summary>
        /// <param name="particleID">This particle ID</param>
        /// <param name="AllParticleHydrodynamics">Hydrodynamics of all particles.</param>
        public void UpdateForcesAndTorque(int particleID, double[] AllParticleHydrodynamics);

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public Vector CalculateParticlePosition(double dt);

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public double CalculateParticleAngle(double dt);

        /// <summary>
        /// Calculate the new particle position
        /// </summary>
        /// <param name="dt"></param>
        public Vector CalculateParticlePositionDuringCollision(double dt);

        /// <summary>
        /// Calculate the new particle angle
        /// </summary>
        /// <param name="dt"></param>
        public double CalculateParticleAngleDuringCollision(double dt);

        /// <summary>
        /// Calculate the new translational velocity of the particle.
        /// </summary>
        /// <param name="dt">Time-step</param>
        public Vector CalculateTranslationalVelocity(double dt);

        /// <summary>
        /// Calculate the new angular velocity of the particle.
        /// </summary>
        /// <param name="dt">Time-step</param>
        public double CalculateAngularVelocity(double dt);

        /// <summary>
        /// Calculate the new translational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        public Vector CalculateTranslationalAcceleration(double dt);

        /// <summary>
        /// Calculate the new rotational acceleration.
        /// </summary>
        /// <param name="dt"></param>
        public double CalculateRotationalAcceleration(double dt);

        /// <summary>
        /// Calls the integration of the hydrodynamic stress at this particles level-set
        /// </summary>
        /// <param name="hydrodynamicsIntegration"></param>
        public Vector CalculateHydrodynamics(ParticleHydrodynamicsIntegration hydrodynamicsIntegration, CellMask cutCells, string[] FluidSpecies, double dt = 0);

        /// <summary>
        /// Calculates the gravitational forces.
        /// </summary>
        /// <param name="fluidDensity"></param>
        /// <param name="tempForces"></param>
        public Vector GravityForce(Vector Gravity);

        /// <summary>
        /// Used for testing
        /// </summary>
        /// <param name="initalTranslation">
        /// The initial translational velocity.
        /// </param>
        /// <param name="initalRotation">
        /// The initial rotational velocity.
        /// </param>
        public void PrescribeHydrodynamicForcesAndTorque(Vector Forces, double Torque, int HistoryPosition);
    }
}

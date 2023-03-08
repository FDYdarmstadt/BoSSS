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

using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Application.XNSERO_Solver {

    /// <summary>
    /// Handles the collision between particles
    /// </summary>
    class ParticleCollision {
        private readonly double TimestepSize;
        private readonly double GridLengthScale;
        private readonly double CoefficientOfRestitution;
        private double AccumulatedCollisionTimestep = 0;
        private double[][] AccumulatedLocalSaveTimestep;
        private bool[][] Overlapping;
        private readonly bool ExceptionWhenOverlap;
        private Vector[][] DistanceVector;
        private Vector[][] ClosestPoints;
        private readonly double[][] WallCoordinates;
        private readonly bool[] IsPeriodicBoundary;
        private readonly double MinDistance;
        private double[][] TemporaryVelocity;
        private Particle[] Particles;
        private int[][] CollisionCluster;
        private int[][] ParticleCollidedWith;

        /// <summary>
        /// Constructor for the collision model.
        /// </summary>
        /// <param name="gridLenghtscale">
        /// Characteristic length of the grid.
        /// </param>
        /// <param name="coefficientOfRestitution">
        /// Static coefficient of restitution
        /// </param>
        /// <param name="dt">
        /// time step size
        /// </param>
        /// <param name="wallCoordinates">
        /// Contains the position of the wall. First index defines vertical/horizontal walls, [0]: vertical, [1]: horizontal. Second index, [0]: left/upper wall [1]: right/lower wall
        /// </param>
        /// <param name="IsPeriodicBoundary">
        /// Determines whether a boundary is periodic. [0]: vertical boundary, [1]: horizontal boundary.
        /// </param>
        /// <param name="minDistance">
        /// Min. distance threshold.
        /// </param>
        public ParticleCollision(double gridLenghtscale, double coefficientOfRestitution, double dt, double[][] wallCoordinates, bool[] IsPeriodicBoundary, double minDistance, bool DetermineOnlyOverlap) {
            CoefficientOfRestitution = coefficientOfRestitution;
            TimestepSize = dt;
            GridLengthScale = gridLenghtscale;
            WallCoordinates = wallCoordinates;
            MinDistance = minDistance;
            this.IsPeriodicBoundary = IsPeriodicBoundary;
            this.ExceptionWhenOverlap = DetermineOnlyOverlap;
        }

        /// <summary>
        /// Update collision forces between two arbitrary particles and add them to forces acting on the corresponding particle
        /// </summary>
        /// <param name="particles">
        /// List of all particles
        /// </param>
        public void Calculate(Particle[] particles, RTree Tree) {
            Particles = particles;
            int[][] potentialCollisionPartners = new int[Particles.Length][];
            double subTimeStepWithoutCollision = 0;
            int ParticleOffset = Particles.Length;
            double distanceThreshold = MinDistance == 0 ? GridLengthScale / 2 : MinDistance;

            // Step 1
            // Loop over time until the particles collide.
            // =======================================================
            while (AccumulatedCollisionTimestep < TimestepSize) {
                CreateCollisionArrarys(Particles.Length);
                double globalMinimalDistance = double.MaxValue;
                // Step 2.1
                // Loop over the distance until it is smaller  
                // than a predefined threshold
                // -------------------------------------------------------
                while (globalMinimalDistance > distanceThreshold) {
                    // Step 2.1.1
                    // Move the particle with the current save time-step.
                    // -------------------------------------------------------
                    MoveParticlesWithSubTimestep(Particles, subTimeStepWithoutCollision);
                    subTimeStepWithoutCollision = double.MaxValue;
                    bool[][] AlreadyAnalysed = new bool[Particles.Length][];
                    for (int p0 = 0; p0 < Particles.Length; p0++) {
                        // Step 2.1.2
                        // Test for particle-wall collisions 
                        // -------------------------------------------------------
                        Vector[] nearFieldWallPoints = GetNearFieldWall(Particles[p0]);
                        for (int w = 0; w < nearFieldWallPoints.Length; w++) {
                            if (!nearFieldWallPoints[w].IsNullOrEmpty()) {
                                MinimumDistance minimumDistance = new(new Particle[] { Particles[p0] }, GridLengthScale);
                                minimumDistance.CalculateParticleWallDistance(new Vector(nearFieldWallPoints[w]));
                                int wallID = ParticleOffset + w; 
                                globalMinimalDistance = CalculateGlobalMinimumDistance(globalMinimalDistance, p0, wallID, minimumDistance.DistanceVector, new Vector[] { minimumDistance.ClosestPoints[0] }, minimumDistance.Overlapping);
                                subTimeStepWithoutCollision = CalculateSubTimeStepWithoutCollision(subTimeStepWithoutCollision, p0, wallID);
                            }
                        }
                        // Step 2.1.3
                        // Test for particle-particle collisions 
                        // -------------------------------------------------------
                        potentialCollisionPartners[p0] = Tree.SearchForOverlap(Particles[p0], p0, TimestepSize - AccumulatedCollisionTimestep).ToArray();
                        AlreadyAnalysed[p0] = new bool[Particles.Length];
                        for (int p1 = 0; p1 < potentialCollisionPartners[p0].Length; p1++) {
                            int secondParticleID = potentialCollisionPartners[p0][p1];
                            if (!AlreadyAnalysed[p0][secondParticleID] && secondParticleID > p0) {
                                AlreadyAnalysed[p0][secondParticleID] = true;
                                MinimumDistance minimumDistance = new(new Particle[] { Particles[p0], Particles[secondParticleID] }, GridLengthScale);
                                minimumDistance.CalculateTwoParticleDistance();
                                globalMinimalDistance = CalculateGlobalMinimumDistance(globalMinimalDistance, p0, secondParticleID, minimumDistance.DistanceVector, minimumDistance.ClosestPoints, minimumDistance.Overlapping);
                                subTimeStepWithoutCollision = CalculateSubTimeStepWithoutCollision(subTimeStepWithoutCollision, p0, secondParticleID);
                            }
                        }
                    }
                    if (subTimeStepWithoutCollision >= 0)
                        AccumulatedCollisionTimestep += subTimeStepWithoutCollision;
                    if (AccumulatedCollisionTimestep >= TimestepSize)
                        break;
                    if (ExceptionWhenOverlap)
                        break;
                }
                if (ExceptionWhenOverlap)
                    break;
                if (AccumulatedCollisionTimestep == double.MaxValue)
                    break;

                ParticleCollidedWith = CreateArrayWithCollidedParticles(distanceThreshold, potentialCollisionPartners);
                CollisionCluster = ClusterCollisionsContainingSameParticles();
                CalculateCollision(distanceThreshold);
                SetParticleVariables(particles, subTimeStepWithoutCollision);
                Tree.UpdateTree(Particles, AccumulatedCollisionTimestep);
            }
        }

        private void CalculateCollision(double distanceThreshold) {
            for (int c = 0; c < CollisionCluster.Length; c++) {
                for (int i = 0; i < CollisionCluster[c].Length; i++) {
                    int currentParticleID = CollisionCluster[c][i];
                    if (IsParticle(currentParticleID)) {
                        for (int j = 0; j < ParticleCollidedWith[currentParticleID].Length; j++) {
                            int secondObjectID = ParticleCollidedWith[currentParticleID][j];
                            CalculateBinaryCollision(currentParticleID, secondObjectID, distanceThreshold);
                        }
                    }
                }
            }
        }

        private void CreateCollisionArrarys(int noOfParticles) {
            AccumulatedLocalSaveTimestep = new double[noOfParticles][];
            Overlapping = new bool[noOfParticles][];
            DistanceVector = new Vector[noOfParticles][];
            ClosestPoints = new Vector[noOfParticles][];
            TemporaryVelocity = new double[noOfParticles][];
            for (int p = 0; p < noOfParticles; p++) {
                AccumulatedLocalSaveTimestep[p] = new double[noOfParticles + 4];
                Overlapping[p] = new bool[noOfParticles + 4];
                DistanceVector[p] = new Vector[noOfParticles + 4];
                ClosestPoints[p] = new Vector[noOfParticles + 4];
                TemporaryVelocity[p] = new double[3];
                TemporaryVelocity[p][0] = Particles[p].Motion.GetTranslationalVelocity()[0];
                TemporaryVelocity[p][1] = Particles[p].Motion.GetTranslationalVelocity()[1];
                TemporaryVelocity[p][2] = Particles[p].Motion.GetRotationalVelocity();
            }
        }

        private void SetParticleVariables(Particle[] particles, double subTimeStepWithoutCollision) {
            for (int p = 0; p < Particles.Length; p++) {
                if (particles[p].IsCollided) {
                    particles[p].Motion.InitializeParticleVelocity(new double[] { TemporaryVelocity[p][0], TemporaryVelocity[p][1] }, TemporaryVelocity[p][2]);
                    particles[p].Motion.InitializeParticleAcceleration(new double[] { 0, 0 }, 0);
                }
                particles[p].Motion.CollisionTimestep = AccumulatedCollisionTimestep - subTimeStepWithoutCollision;
                CollisionCluster.Clear();
                ParticleCollidedWith.Clear();
            }
        }

        private double CalculateSubTimeStepWithoutCollision(double subTimeStepWithoutCollision, int firstObjectID, int secondObjectID) {
            double temp_SaveTimeStep = DynamicTimestep(firstObjectID, secondObjectID);
            AccumulatedLocalSaveTimestep[firstObjectID][secondObjectID] += temp_SaveTimeStep;
            if (IsParticle(secondObjectID))
                AccumulatedLocalSaveTimestep[secondObjectID][firstObjectID] += temp_SaveTimeStep;
            if (temp_SaveTimeStep < subTimeStepWithoutCollision && temp_SaveTimeStep > 0) {
                subTimeStepWithoutCollision = temp_SaveTimeStep;
            }
            return subTimeStepWithoutCollision;
        }

        private double CalculateGlobalMinimumDistance(double globalMinimalDistance, int firstObjectID, int secondObjectID, Vector temp_DistanceVector, Vector[] temp_ClosestPoints, bool temp_Overlapping) {
            Overlapping[firstObjectID][secondObjectID] = temp_Overlapping;
            ClosestPoints[firstObjectID][secondObjectID] = temp_ClosestPoints[0];
            DistanceVector[firstObjectID][secondObjectID] = new Vector(temp_DistanceVector);
            if (IsParticle(secondObjectID)) {
                Overlapping[secondObjectID][firstObjectID] = temp_Overlapping;
                ClosestPoints[secondObjectID][firstObjectID] = temp_ClosestPoints[1];
                temp_DistanceVector.ScaleInPlace(-1);
                DistanceVector[secondObjectID][firstObjectID] = new Vector(temp_DistanceVector);
            }
            if (DistanceVector[firstObjectID][secondObjectID].Abs() < globalMinimalDistance) {
                globalMinimalDistance = DistanceVector[firstObjectID][secondObjectID].Abs();
            }
            if (temp_Overlapping) {
                if (ExceptionWhenOverlap)
                    throw new Exception("Static particles overlap");
                DistanceVector[firstObjectID][secondObjectID] = new Vector(Particles[firstObjectID].Motion.GetPosition(0) - Particles[secondObjectID].Motion.GetPosition(0));
                DistanceVector[secondObjectID][firstObjectID] = new Vector(Particles[secondObjectID].Motion.GetPosition(0) - Particles[firstObjectID].Motion.GetPosition(0));
                globalMinimalDistance = 0;
            }

            return globalMinimalDistance;
        }

        private int[][] ClusterCollisionsContainingSameParticles() {
            List<int[]> globalParticleCluster = new List<int[]>();
            bool[] partOfCollisionCluster = new bool[Particles.Length + 4];
            for (int p0 = 0; p0 < Particles.Length; p0++) {
                if (!partOfCollisionCluster[p0]) {
                    List<int> currentParticleCluster = new List<int> { p0 };
                    partOfCollisionCluster[p0] = true;
                    for (int p1 = 1; p1 < ParticleCollidedWith[p0].Count(); p1++) {
                        currentParticleCluster.Add(ParticleCollidedWith[p0][p1]);
                        partOfCollisionCluster[ParticleCollidedWith[p0][p1]] = true;
                        FindCollisionClusterRecursive(ParticleCollidedWith, ParticleCollidedWith[p0][p1], currentParticleCluster, partOfCollisionCluster);
                    }
                    globalParticleCluster.Add(currentParticleCluster.ToArray());
                }
            }
            return globalParticleCluster.ToArray();
        }

        private void FindCollisionClusterRecursive(int[][] ParticleCollidedWith, int p0, List<int> currentParticleCluster, bool[] PartOfCollisionCluster) {
            if (IsParticle(p0)) {
                for (int p1 = 1; p1 < ParticleCollidedWith[p0].Count(); p1++) {
                    if (!PartOfCollisionCluster[ParticleCollidedWith[p0][p1]]) {
                        currentParticleCluster.Add(ParticleCollidedWith[p0][p1]);
                        PartOfCollisionCluster[p0] = true;
                        PartOfCollisionCluster[ParticleCollidedWith[p0][p1]] = true;
                        FindCollisionClusterRecursive(ParticleCollidedWith, p1, currentParticleCluster, PartOfCollisionCluster);
                    }
                }
            }
        }

        private int[][] CreateArrayWithCollidedParticles(double distanceThreshold, int[][] PotentialCollisionPartners) {
            int ParticleOffset = Particles.Length;
            int[][] particleCollidedWith = new int[Particles.Length][];
            for (int p0 = 0; p0 < Particles.Length; p0++) {
                List<int> currentParticleCollidedWith = new List<int>();
                for (int w = 0; w < 4; w++) {
                    FindCollisionPartners(p0, ParticleOffset + w, currentParticleCollidedWith, distanceThreshold);
                }
                for (int p1 = 0; p1 < PotentialCollisionPartners[p0].Length; p1++) {
                    FindCollisionPartners(p0, PotentialCollisionPartners[p0][p1], currentParticleCollidedWith, distanceThreshold);
                }
                particleCollidedWith[p0] = currentParticleCollidedWith.ToArray();
            }
            return particleCollidedWith;
        }

        private void FindCollisionPartners(int particle, int potentialCollisionPartner, List<int> currentParticleCollidedWith, double distanceThreshold) {
            if ((DistanceVector[particle][potentialCollisionPartner].Abs() <= distanceThreshold && AccumulatedCollisionTimestep < TimestepSize && AccumulatedLocalSaveTimestep[particle][potentialCollisionPartner] > 0) || Overlapping[particle][potentialCollisionPartner]) {
                int insertAtIndex = currentParticleCollidedWith.Count();
                for (int i = 1; i < currentParticleCollidedWith.Count(); i++) {
                    if (DistanceVector[particle][currentParticleCollidedWith[i]].Abs() > DistanceVector[particle][potentialCollisionPartner].Abs())
                        insertAtIndex = i;
                }
                currentParticleCollidedWith.Insert(insertAtIndex, potentialCollisionPartner);
            }
        }

        private double DynamicTimestep(int FirstParticleID, int SecondParticleID) {
            double detectCollisionVn_P0 = 0;
            double detectCollisionVn_P1 = 0;
            Vector normalVector = DistanceVector[FirstParticleID][SecondParticleID];
            double distance = normalVector.Abs();
            normalVector /= distance;
            if (Particles[FirstParticleID].Motion.IncludeTranslation() || Particles[FirstParticleID].Motion.IncludeRotation())
                detectCollisionVn_P0 = CalculateNormalSurfaceVelocity(FirstParticleID, normalVector, ClosestPoints[FirstParticleID][SecondParticleID]);
            if (IsParticle(SecondParticleID)) {
                if (Particles[FirstParticleID].Motion.IncludeTranslation() || Particles[FirstParticleID].Motion.IncludeRotation())
                    detectCollisionVn_P1 = CalculateNormalSurfaceVelocity(SecondParticleID, normalVector, ClosestPoints[SecondParticleID][FirstParticleID]);
            } else {
                detectCollisionVn_P0 = 10 * detectCollisionVn_P0;
            }
            return (detectCollisionVn_P1 - detectCollisionVn_P0 == 0) ? double.MaxValue : 0.01 * distance / (detectCollisionVn_P1 - detectCollisionVn_P0);
        }

        private double CalculateNormalSurfaceVelocity(int particleID, Vector normalVector, Vector closestPoint) {
            Vector pointVelocity = new Vector(2);
            Vector radialVector = Particles[particleID].CalculateRadialVector(closestPoint);
            pointVelocity[0] = TemporaryVelocity[particleID][0] - TemporaryVelocity[particleID][2] * radialVector[1];
            pointVelocity[1] = TemporaryVelocity[particleID][1] + TemporaryVelocity[particleID][2] * radialVector[0];
            return pointVelocity * normalVector;
        }

        private static void MoveParticlesWithSubTimestep(Particle[] particles, double dynamicTimestep) {
            for (int p = 0; p < particles.Length; p++) {
                Particle currentParticle = particles[p];
                if (dynamicTimestep != 0) {
                    currentParticle.Motion.CollisionParticlePositionAndAngle(dynamicTimestep);
                }
            }
        }


        /// <summary>
        /// Physical collision model
        /// </summary>
        /// <param name="particleID"></param>
        /// <param name="secondObjectID"></param>
        /// <param name="threshold"></param>
        private void CalculateBinaryCollision(int particleID, int secondObjectID, double threshold) {
            if (DistanceVector[particleID][secondObjectID].Count != 2)
                throw new NotImplementedException("Physical collision model only implemented for 2D");

            double distance = DistanceVector[particleID][secondObjectID].Abs();
            Vector normalVector = DistanceVector[particleID][secondObjectID];
            normalVector.NormalizeInPlace();
            Vector tangentialVector = new(-normalVector[1], normalVector[0]);

            if (distance <= threshold || Overlapping[particleID][secondObjectID]) {
                double[] normalSurfaceVelocity = new double[2];
                if (Particles[0].Motion.IncludeTranslation() || Particles[0].Motion.IncludeRotation())
                    normalSurfaceVelocity[0] = CalculateNormalSurfaceVelocity(particleID, normalVector, ClosestPoints[particleID][secondObjectID]);
                if (IsParticle(secondObjectID)) {
                    if (Particles[1].Motion.IncludeTranslation() || Particles[1].Motion.IncludeRotation())
                        normalSurfaceVelocity[1] = CalculateNormalSurfaceVelocity(secondObjectID, normalVector, ClosestPoints[secondObjectID][particleID]);
                }
                if (normalSurfaceVelocity[1] - normalSurfaceVelocity[0] <= 0)
                    return;

                Particles[particleID].CalculateEccentricity(normalVector, ClosestPoints[particleID][secondObjectID]);
                if (IsParticle(secondObjectID))
                    Particles[secondObjectID].CalculateEccentricity(normalVector, ClosestPoints[secondObjectID][particleID]);

                double[] collisionCoefficient = CalculateCollisionCoefficient(particleID, secondObjectID, normalVector);

                Vector velocityP0 = CalculateNormalAndTangentialVelocity(particleID, normalVector);
                Vector radialVectorP0 = Particles[particleID].CalculateRadialVector(ClosestPoints[particleID][secondObjectID]);
                Vector tempVel0 = Particles[particleID].Motion.IncludeTranslation()
                    ? (velocityP0[0] + collisionCoefficient[0] / Particles[particleID].Mass) * normalVector + (velocityP0[1] + collisionCoefficient[1] / Particles[particleID].Mass) * tangentialVector
                    : new Vector(0, 0);
                TemporaryVelocity[particleID][0] = tempVel0[0];
                TemporaryVelocity[particleID][1] = tempVel0[1];
                TemporaryVelocity[particleID][2] = Particles[particleID].Motion.IncludeRotation() ? TemporaryVelocity[particleID][2] + (radialVectorP0[0] * normalVector[1] - radialVectorP0[1] * normalVector[0]) * collisionCoefficient[0] / Particles[particleID].MomentOfInertia + (radialVectorP0[0] * tangentialVector[1] - radialVectorP0[1] * tangentialVector[0]) * collisionCoefficient[1] / Particles[particleID].MomentOfInertia : 0;
                Particles[particleID].IsCollided = true;
                Overlapping[particleID][secondObjectID] = false;

                if (IsParticle(secondObjectID)) {
                    Vector velocityP1 = CalculateNormalAndTangentialVelocity(secondObjectID, normalVector);
                    Vector tempVel1 = Particles[secondObjectID].Motion.IncludeTranslation()
                        ? (velocityP1[0] - collisionCoefficient[0] / Particles[secondObjectID].Mass) * normalVector + (velocityP1[1] - collisionCoefficient[1] / Particles[secondObjectID].Mass) * tangentialVector
                        : new Vector(0, 0);
                    Vector radialVectorP1 = Particles[secondObjectID].CalculateRadialVector(ClosestPoints[secondObjectID][particleID]);
                    TemporaryVelocity[secondObjectID][0] = tempVel1[0];
                    TemporaryVelocity[secondObjectID][1] = tempVel1[1];
                    TemporaryVelocity[secondObjectID][2] = Particles[secondObjectID].Motion.IncludeRotation() ? TemporaryVelocity[secondObjectID][2] - (radialVectorP1[0] * normalVector[1] - radialVectorP1[1] * normalVector[0]) * collisionCoefficient[0] / Particles[secondObjectID].MomentOfInertia - (radialVectorP1[0] * tangentialVector[1] - radialVectorP1[1] * tangentialVector[0]) * collisionCoefficient[1] / Particles[secondObjectID].MomentOfInertia : 0;
                    Overlapping[secondObjectID][particleID] = false;
                    Particles[secondObjectID].IsCollided = true;
                    Console.WriteLine("Particle " + particleID + " and particle " + secondObjectID + " collided");
                } else {
                    Console.WriteLine("Particle " + particleID + " and wall " + secondObjectID + " collided");
                }
            }
        }

        private Vector CalculateNormalAndTangentialVelocity(int particleID, Vector normalVector) {
            Vector tangentialVector = new(-normalVector[1], normalVector[0]);
            Vector velocity = new(TemporaryVelocity[particleID][0], TemporaryVelocity[particleID][1]);
            return new Vector(velocity * normalVector, velocity * tangentialVector);
        }

        private double[] CalculateCollisionCoefficient(int particleID, int secondObjectID, Vector normalVector) {
            Vector tangentialVector = new(-normalVector[1], normalVector[0]);
            Vector[] translationalVelocity = new Vector[2];
            translationalVelocity[0] = CalculateNormalAndTangentialVelocity(particleID, normalVector);
            translationalVelocity[1] = new Vector(0, 0);
            double[] massReciprocal = new double[] { 0, 0 };
            double[] normalMomentOfInertiaReciprocal = new double[] { 0, 0 };
            double[] tangentialMomentOfInertiaReciprocal = new double[] { 0, 0 };
            double[] normalEccentricity = new double[] { 0, 0 };
            double[] tangentialEccentricity = new double[] { 0, 0 };

            massReciprocal[0] = Particles[particleID].Motion.IncludeTranslation() ? 1 / Particles[particleID].Mass : 0;
            normalMomentOfInertiaReciprocal[0] = Particles[particleID].Motion.IncludeRotation() ? Particles[particleID].CalculateSecondOrderEccentricity(normalVector, ClosestPoints[particleID][secondObjectID]) / Particles[particleID].MomentOfInertia : 0;
            tangentialMomentOfInertiaReciprocal[0] = Particles[particleID].Motion.IncludeRotation() ? Particles[particleID].CalculateSecondOrderEccentricity(tangentialVector, ClosestPoints[particleID][secondObjectID]) / Particles[particleID].MomentOfInertia : 0;
            normalEccentricity[0] = Particles[particleID].CalculateEccentricity(normalVector, ClosestPoints[particleID][secondObjectID]);
            tangentialEccentricity[0] = Particles[particleID].CalculateEccentricity(tangentialVector, ClosestPoints[particleID][secondObjectID]);

            if (IsParticle(secondObjectID)) {
                translationalVelocity[1] = CalculateNormalAndTangentialVelocity(secondObjectID, normalVector);
                massReciprocal[1] = Particles[secondObjectID].Motion.IncludeTranslation() ? 1 / Particles[secondObjectID].Mass : 0;
                normalMomentOfInertiaReciprocal[1] = Particles[secondObjectID].Motion.IncludeRotation() ? Particles[secondObjectID].CalculateSecondOrderEccentricity(normalVector, ClosestPoints[secondObjectID][particleID]) / Particles[secondObjectID].MomentOfInertia : 0;
                tangentialMomentOfInertiaReciprocal[1] = Particles[secondObjectID].Motion.IncludeRotation() ? Particles[secondObjectID].CalculateSecondOrderEccentricity(tangentialVector, ClosestPoints[secondObjectID][particleID]) / Particles[secondObjectID].MomentOfInertia : 0;
                normalEccentricity[1] = Particles[secondObjectID].CalculateEccentricity(normalVector, ClosestPoints[secondObjectID][particleID]);
                tangentialEccentricity[1] = Particles[secondObjectID].CalculateEccentricity(tangentialVector, ClosestPoints[secondObjectID][particleID]);
            }
            double[] collisionCoefficient = new double[2];
            collisionCoefficient[0] = -(1 + CoefficientOfRestitution) * ((translationalVelocity[0][0] - translationalVelocity[1][0]) / (massReciprocal[0] + massReciprocal[1] + normalMomentOfInertiaReciprocal[0] + normalMomentOfInertiaReciprocal[1]));
            collisionCoefficient[1] = 0;// -(translationalVelocity[0][1] - translationalVelocity[1][1]) / (massReciprocal[0] + massReciprocal[1] + tangentialMomentOfInertiaReciprocal[0] + tangentialMomentOfInertiaReciprocal[1]);
            double tempRotVelocity2 = 0;
            if (IsParticle(secondObjectID))
                tempRotVelocity2 = TemporaryVelocity[secondObjectID][2];
            collisionCoefficient[0] -= (1 + CoefficientOfRestitution) * ((normalEccentricity[0] * TemporaryVelocity[particleID][2] - normalEccentricity[1] * tempRotVelocity2) / (massReciprocal[0] + massReciprocal[1] + normalMomentOfInertiaReciprocal[0] + normalMomentOfInertiaReciprocal[1]));
            collisionCoefficient[1] -= 0;// (tangentialEccentricity[0] * TemporaryVelocity[p0][2] - tangentialEccentricity[1] * tempRotVelocity2) / (massReciprocal[0] + massReciprocal[1] + tangentialMomentOfInertiaReciprocal[0] + tangentialMomentOfInertiaReciprocal[1]);
            return collisionCoefficient;
        }

        private Vector[] GetNearFieldWall(Particle particle) {
            Vector[] nearFieldWallPoint = new Vector[4];
            Vector particlePosition = particle.Motion.GetPosition(0);
            double particleMaxLengthscale = particle.GetLengthScales().Max();
            for (int w0 = 0; w0 < WallCoordinates.Length; w0++) {
                for (int w1 = 0; w1 < WallCoordinates[w0].Length; w1++) {
                    if (WallCoordinates[w0][w1] != 0 && !IsPeriodicBoundary[w0]) {
                        double minDistance = Math.Abs(particlePosition[w0] - WallCoordinates[w0][w1]) - particleMaxLengthscale;
                        if (minDistance < 5 * GridLengthScale) {
                            if (w0 == 0)
                                nearFieldWallPoint[w0 + w1] = new Vector(WallCoordinates[0][w1], particlePosition[1]);
                            else
                                nearFieldWallPoint[w0 * 2 + w1] = new Vector(particlePosition[0], WallCoordinates[1][w1]);
                        }
                    }
                }
            }
            return nearFieldWallPoint;
        }

        /// <summary>
        /// Returns true if <paramref name="Object"/> is a particle. Returns <see langword="false"/> if <paramref name="Object"/> is a wall.
        /// </summary>
        /// <param name="Object"></param>
        /// <returns></returns>
        private static bool IsParticle(Particle Object) {
            return Object != null;
        }

        private bool IsParticle(int objectID) {
            return objectID < Particles.Length;
        }
    }
}

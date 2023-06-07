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
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.XNSERO_Solver {

    /// <summary>
    /// Handles the collision between particles
    /// </summary>
    class ParticleCollision {
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
        public ParticleCollision(double gridLenghtscale, double coefficientOfRestitution, double dt, double[][] wallCoordinates, bool[] IsPeriodicBoundary) {
            CoefficientOfRestitution = coefficientOfRestitution;
            TimestepSize = dt;
            GridLengthScale = gridLenghtscale;
            WallCoordinates = wallCoordinates;
            this.IsPeriodicBoundary = IsPeriodicBoundary;
            DistanceThreshold = GridLengthScale / 2;
        }

        private readonly double TimestepSize;
        private readonly double GridLengthScale;
        private readonly double CoefficientOfRestitution;
        private double AccumulatedTimeWithoutCollision = 0;
        private double[][] AccumulatedLocalSaveTimestep;
        private bool[][] Overlapping;
        private Vector[][] DistanceVector;
        private Vector[][] ClosestPoints;
        private readonly double[][] WallCoordinates;
        private readonly bool[] IsPeriodicBoundary;
        private double[][] TemporaryVelocity;
        private Particle[] Particles;
        private int[][] CollisionCluster;
        private int[][] ParticleCollidedWith;
        private readonly double DistanceThreshold;

        /// <summary>
        /// Update collision forces between two arbitrary particles and add them to forces acting on the corresponding particle
        /// </summary>
        /// <param name="Particles">
        /// List of all particles
        /// </param>
        public void Calculate(Particle[] Particles, RTree Tree) {
            this.Particles = Particles;
            int[][] potentialCollisionPartners = new int[this.Particles.Length][];
            double timeWithoutCollision = 0;
            int noOfParticles = this.Particles.Length;

            // Step 1
            // Loop over time until the particles collide.
            // =======================================================
            while (AccumulatedTimeWithoutCollision < TimestepSize) {
                CreateCollisionArrarys();
                double globalMinimalDistance = double.MaxValue;
                // Step 2.1
                // Loop over the distance until it is smaller  
                // than a predefined threshold
                // -------------------------------------------------------
                while (globalMinimalDistance > DistanceThreshold) {
                    // Step 2.1.1
                    // Move the particle with the current save time-step.
                    // -------------------------------------------------------
                    if (timeWithoutCollision != 0)
                        foreach (Particle p in this.Particles)
                            p.Motion.MoveParticleDuringCollision(timeWithoutCollision);
                    timeWithoutCollision = double.MaxValue;
                    bool[][] alreadyAnalysed = new bool[this.Particles.Length][];
                    for (int p0 = 0; p0 < this.Particles.Length; p0++) {
                        // Step 2.1.2
                        // Test for particle-wall collisions 
                        // -------------------------------------------------------
                        Vector[] nearFieldWallPoints = GetNearFieldWall(this.Particles[p0]);
                        for (int w = 0; w < nearFieldWallPoints.Length; w++) {
                            if (!nearFieldWallPoints[w].IsNullOrEmpty()) {
                                DistanceAlgorithm distanceAlg = new(new Particle[] { this.Particles[p0] }, GridLengthScale);
                                distanceAlg.CalculateParticleWallDistance(new Vector(nearFieldWallPoints[w]));
                                SaveDistanceProperties(p0, noOfParticles + w, distanceAlg.DistanceVector, distanceAlg.ClosestPoints, distanceAlg.Overlapping);
                                if (distanceAlg.Overlapping)
                                    globalMinimalDistance = 0;
                                else if (distanceAlg.DistanceVector.Abs() < globalMinimalDistance) 
                                    globalMinimalDistance = distanceAlg.DistanceVector.Abs();
                                timeWithoutCollision = CalculateTimeWithoutCollision(timeWithoutCollision, p0, noOfParticles + w);
                            }
                        }
                        // Step 2.1.3
                        // Test for particle-particle collisions 
                        // -------------------------------------------------------
                        potentialCollisionPartners[p0] = Tree.SearchForOverlap(this.Particles[p0], p0, TimestepSize - AccumulatedTimeWithoutCollision).ToArray();
                        alreadyAnalysed[p0] = new bool[this.Particles.Length];
                        for (int p1 = 0; p1 < potentialCollisionPartners[p0].Length; p1++) {
                            int secondParticleID = potentialCollisionPartners[p0][p1];
                            if (!alreadyAnalysed[p0][secondParticleID] && secondParticleID > p0) {
                                alreadyAnalysed[p0][secondParticleID] = true;
                                DistanceAlgorithm distanceAlg = new(new Particle[] { this.Particles[p0], this.Particles[secondParticleID] }, GridLengthScale);
                                distanceAlg.CalculateTwoParticleDistance();
                                SaveDistanceProperties(p0, secondParticleID, distanceAlg.DistanceVector, distanceAlg.ClosestPoints, distanceAlg.Overlapping);
                                if (distanceAlg.Overlapping)
                                    globalMinimalDistance = 0;
                                else if (distanceAlg.DistanceVector.Abs() < globalMinimalDistance) 
                                    globalMinimalDistance = distanceAlg.DistanceVector.Abs();
                                timeWithoutCollision = CalculateTimeWithoutCollision(timeWithoutCollision, p0, secondParticleID);
                            }
                        }
                    }
                    if (timeWithoutCollision >= 0)
                        AccumulatedTimeWithoutCollision += timeWithoutCollision;
                    if (AccumulatedTimeWithoutCollision >= TimestepSize)
                        break;
                }
                if (AccumulatedTimeWithoutCollision == double.MaxValue)// Particles might be close but move away from each other
                    break;

                ParticleCollidedWith = CollisionPartners(potentialCollisionPartners);
                CollisionCluster = FindCollisionCluster();
                for (int c = 0; c < CollisionCluster.Length; c++) 
                    foreach (int currentParticleID in CollisionCluster[c]) 
                        if (IsParticle(currentParticleID)) 
                            foreach (int secondObjectID in ParticleCollidedWith[currentParticleID]) 
                                CalculateCollision(currentParticleID, secondObjectID);
                SetResultsToParticles(timeWithoutCollision);
                Tree.UpdateTree(this.Particles, AccumulatedTimeWithoutCollision);
            }
        }

        /// <summary>
        /// Create all arrays necessary for the calculation of the distance and the collision.
        /// </summary>
        /// <param name="noOfParticles"></param>
        private void CreateCollisionArrarys() {
            int noOfParticles = Particles.Length;
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

        /// <summary>
        /// Returns a point at each wall close to the particle. If the particle is too far away from the wall, no point is defined.
        /// </summary>
        /// <param name="particle"></param>
        /// <returns></returns>
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
        /// Calculates the time-span in which no collision happens.
        /// </summary>
        /// <param name="TimeWithoutCollision"></param>
        /// <param name="FirstParticleID"></param>
        /// <param name="SecondParticleID"></param>
        /// <returns></returns>
        private double CalculateTimeWithoutCollision(double TimeWithoutCollision, int FirstParticleID, int SecondParticleID) {
            double[] normalVelocity = new double[] { 0, 0 };
            Vector distanceVector = DistanceVector[FirstParticleID][SecondParticleID];
            double distance = distanceVector.Abs();
            Vector normalVector = distanceVector.Normalize();

            if (IsMoving(FirstParticleID))
                normalVelocity[0] = CalculateNormalSurfaceVelocity(FirstParticleID, normalVector, ClosestPoints[FirstParticleID][SecondParticleID]);
            if (IsParticle(SecondParticleID)) {
                if (IsMoving(SecondParticleID))
                    normalVelocity[1] = CalculateNormalSurfaceVelocity(SecondParticleID, normalVector, ClosestPoints[SecondParticleID][FirstParticleID]);
            } else {
                normalVelocity[0] = 10 * normalVelocity[0];
            }
            double tempTimeWithoutCollision = (normalVelocity[1] - normalVelocity[0] == 0) ? double.MaxValue : 0.01 * distance / (normalVelocity[1] - normalVelocity[0]);
            AccumulatedLocalSaveTimestep[FirstParticleID][SecondParticleID] += tempTimeWithoutCollision;
            if (IsParticle(SecondParticleID))
                AccumulatedLocalSaveTimestep[SecondParticleID][FirstParticleID] += tempTimeWithoutCollision;

            if (tempTimeWithoutCollision < TimeWithoutCollision && tempTimeWithoutCollision > 0) {
                TimeWithoutCollision = tempTimeWithoutCollision;
            }
            return TimeWithoutCollision;
        }

        /// <summary>
        /// Set results to particles.
        /// </summary>
        /// <param name="TimeWithoutCollision"></param>
        private void SetResultsToParticles(double TimeWithoutCollision) {
            for (int p = 0; p < Particles.Length; p++) {
                if (Particles[p].IsCollided) {
                    Particles[p].Motion.InitializeParticleVelocity(new double[] { TemporaryVelocity[p][0], TemporaryVelocity[p][1] }, TemporaryVelocity[p][2]);
                    Particles[p].Motion.InitializeParticleAcceleration(new double[] { 0, 0 }, 0);
                }
                Particles[p].Motion.CollisionTimestep = AccumulatedTimeWithoutCollision - TimeWithoutCollision;
                CollisionCluster.Clear();
                ParticleCollidedWith.Clear();
            }
        }

        /// <summary>
        /// Save individual results for a single particle pair to array.
        /// </summary>
        /// <param name="firstObjectID"></param>
        /// <param name="secondObjectID"></param>
        /// <param name="LocalDistanceVector"></param>
        /// <param name="LocalClosestPoints"></param>
        /// <param name="LocalOverlapping"></param>
        private void SaveDistanceProperties(int firstObjectID, int secondObjectID, Vector LocalDistanceVector, Vector[] LocalClosestPoints, bool LocalOverlapping) {
            Overlapping[firstObjectID][secondObjectID] = LocalOverlapping;
            ClosestPoints[firstObjectID][secondObjectID] = LocalClosestPoints[0];
            DistanceVector[firstObjectID][secondObjectID] = new Vector(LocalDistanceVector);
            if (IsParticle(secondObjectID)) {
                Overlapping[secondObjectID][firstObjectID] = LocalOverlapping;
                ClosestPoints[secondObjectID][firstObjectID] = LocalClosestPoints[1];
                LocalDistanceVector.ScaleInPlace(-1);
                DistanceVector[secondObjectID][firstObjectID] = new Vector(LocalDistanceVector);
            }
        }

        /// <summary>
        /// Write all collision partners (index 2) of a specific particle (index 1) to an multidimensional array.
        /// </summary>
        /// <param name="PotentialCollisionPartners"></param>
        /// <returns></returns>
        private int[][] CollisionPartners(int[][] PotentialCollisionPartners) {
            int noOfParticles = Particles.Length;
            int[][] collisionPartners = new int[Particles.Length][];
            for (int p0 = 0; p0 < Particles.Length; p0++) {
                List<int> collisionPartnersList = new();
                for (int w = 0; w < 4; w++)
                    FindCollisionPartners(p0, noOfParticles + w, collisionPartnersList);
                for (int p1 = 0; p1 < PotentialCollisionPartners[p0].Length; p1++)
                    FindCollisionPartners(p0, PotentialCollisionPartners[p0][p1], collisionPartnersList);
                collisionPartners[p0] = collisionPartnersList.ToArray();
            }
            return collisionPartners;
        }

        /// <summary>
        /// Checks whether an object <paramref name="SecondObjectID"/> is close enough towards <paramref name="ParticleID"/> for a collision.
        /// </summary>
        /// <param name="ParticleID"></param>
        /// <param name="SecondObjectID"></param>
        /// <param name="CollisionPartners"></param>
        private void FindCollisionPartners(int ParticleID, int SecondObjectID, List<int> CollisionPartners) {
            if ((DistanceVector[ParticleID][SecondObjectID].Abs() <= DistanceThreshold && AccumulatedTimeWithoutCollision < TimestepSize && AccumulatedLocalSaveTimestep[ParticleID][SecondObjectID] > 0) || Overlapping[ParticleID][SecondObjectID]) {
                int insertAtIndex = CollisionPartners.Count;
                for (int i = 1; i < CollisionPartners.Count; i++) {
                    if (DistanceVector[ParticleID][CollisionPartners[i]].Abs() > DistanceVector[ParticleID][SecondObjectID].Abs())
                        insertAtIndex = i;
                }
                CollisionPartners.Insert(insertAtIndex, SecondObjectID);
            }
        }
        /// <summary>
        /// Forms clusters, which contain particles which are close and might collide.
        /// </summary>
        /// <returns></returns>
        private int[][] FindCollisionCluster() {
            List<int[]> globalParticleCluster = new();
            bool[] partOfCollisionCluster = new bool[Particles.Length + 4];
            for (int p0 = 0; p0 < Particles.Length; p0++) {
                if (!partOfCollisionCluster[p0]) {
                    List<int> currentParticleCluster = new() { p0 };
                    partOfCollisionCluster[p0] = true;
                    for (int p1 = 1; p1 < ParticleCollidedWith[p0].Length; p1++) {
                        currentParticleCluster.Add(ParticleCollidedWith[p0][p1]);
                        partOfCollisionCluster[ParticleCollidedWith[p0][p1]] = true;
                        FindCollisionClusterRecursive(ParticleCollidedWith, ParticleCollidedWith[p0][p1], currentParticleCluster, partOfCollisionCluster);
                    }
                    globalParticleCluster.Add(currentParticleCluster.ToArray());
                }
            }
            return globalParticleCluster.ToArray();
        }

        /// <summary>
        /// Recursive call to find the collision clusters, see <see cref="FindCollisionCluster"/>.
        /// </summary>
        /// <param name="ParticleCollidedWith"></param>
        /// <param name="id"></param>
        /// <param name="currentParticleCluster"></param>
        /// <param name="PartOfCollisionCluster"></param>
        private void FindCollisionClusterRecursive(int[][] ParticleCollidedWith, int id, List<int> currentParticleCluster, bool[] PartOfCollisionCluster) {
            if (IsParticle(id)) {
                for (int p1 = 1; p1 < ParticleCollidedWith[id].Length; p1++) {
                    if (!PartOfCollisionCluster[ParticleCollidedWith[id][p1]]) {
                        currentParticleCluster.Add(ParticleCollidedWith[id][p1]);
                        PartOfCollisionCluster[id] = true;
                        PartOfCollisionCluster[ParticleCollidedWith[id][p1]] = true;
                        FindCollisionClusterRecursive(ParticleCollidedWith, p1, currentParticleCluster, PartOfCollisionCluster);
                    }
                }
            }
        }

        /// <summary>
        /// Physical collision model
        /// </summary>
        /// <param name="ParticleID"></param>
        /// <param name="SecondObjectID"></param>
        /// <param name="threshold"></param>
        private void CalculateCollision(int ParticleID, int SecondObjectID) {
            if (DistanceVector[ParticleID][SecondObjectID].Count != 2)
                throw new NotImplementedException("Physical collision model only implemented for 2D");

            double distance = DistanceVector[ParticleID][SecondObjectID].Abs();
            Vector normalVector = DistanceVector[ParticleID][SecondObjectID];
            normalVector.NormalizeInPlace();

            if ((distance <= DistanceThreshold || Overlapping[ParticleID][SecondObjectID]) && WillCollide(ParticleID, SecondObjectID, normalVector)) { 
                Particles[ParticleID].CalculateEccentricity(normalVector, ClosestPoints[ParticleID][SecondObjectID]);
                if (IsParticle(SecondObjectID))
                    Particles[SecondObjectID].CalculateEccentricity(normalVector, ClosestPoints[SecondObjectID][ParticleID]);
                double collisionCoefficient = CalculateCollisionMomentum(ParticleID, SecondObjectID, normalVector);
                CalculateVelocityAfterCollision(ParticleID, SecondObjectID, normalVector, collisionCoefficient);
                if (IsParticle(SecondObjectID)) {
                    CalculateVelocityAfterCollision(SecondObjectID, ParticleID, normalVector, collisionCoefficient);
                    Console.WriteLine("Particle " + ParticleID + " and particle " + SecondObjectID + " collided");
                } else {
                    Console.WriteLine("Particle " + ParticleID + " and wall " + SecondObjectID + " collided");
                }
            }
        }

        /// <summary>
        /// Check whether two close objects will collide based on the relative velocity in the direcetion of  <paramref name="normalVector"/>.
        /// </summary>
        /// <param name="ParticleID"></param>
        /// <param name="SecondObjectID"></param>
        /// <param name="normalVector"></param>
        /// <returns></returns>
        private bool WillCollide(int ParticleID, int SecondObjectID, Vector normalVector) {
            double[] normalSurfaceVelocity = new double[2];
            if (Particles[0].Motion.IncludeTranslation() || Particles[0].Motion.IncludeRotation())
                normalSurfaceVelocity[0] = CalculateNormalSurfaceVelocity(ParticleID, normalVector, ClosestPoints[ParticleID][SecondObjectID]);
            if (IsParticle(SecondObjectID)) {
                if (Particles[1].Motion.IncludeTranslation() || Particles[1].Motion.IncludeRotation())
                    normalSurfaceVelocity[1] = CalculateNormalSurfaceVelocity(SecondObjectID, normalVector, ClosestPoints[SecondObjectID][ParticleID]);
            }
            return normalSurfaceVelocity[1] - normalSurfaceVelocity[0] > 0;
        }

        /// <summary>
        /// Calculates the velocity of the particle <paramref name="ParticleID"/> in the direction of <paramref name="NormalVector"/> at <paramref name="SurfacePoint"/>.
        /// </summary>
        private double CalculateNormalSurfaceVelocity(int ParticleID, Vector NormalVector, Vector SurfacePoint) {
            Vector radialVector = Particles[ParticleID].CalculateRadialVector(SurfacePoint);
            Vector pointVelocity = new(TemporaryVelocity[ParticleID][0] - TemporaryVelocity[ParticleID][2] * radialVector[1],
                                       TemporaryVelocity[ParticleID][1] + TemporaryVelocity[ParticleID][2] * radialVector[0]);
            return pointVelocity * NormalVector;
        }

        /// <summary>
        /// Calculate the momentum exchanged during the collision between <paramref name="ParticleID"/> and <paramref name="secondObjectID"/>.
        /// </summary>
        /// <param name="ParticleID"></param>
        /// <param name="secondObjectID"></param>
        /// <param name="NormalVector"></param>
        /// <returns></returns>
        private double CalculateCollisionMomentum(int ParticleID, int secondObjectID, Vector NormalVector) {
            Vector tangentialVector = new(-NormalVector[1], NormalVector[0]);
            Vector translationalVelocity = new(TemporaryVelocity[ParticleID][0], TemporaryVelocity[ParticleID][1]);
            double[] normalVelocity = new double[] { translationalVelocity * NormalVector, 0 };

            //Arrays for different physical properties, first entry: first particle, second entry: second particle or wall
            double[] massReciprocal = new double[] { 0, 0 };
            double[] normalMomentOfInertiaReciprocal = new double[] { 0, 0 };
            double[] tangentialMomentOfInertiaReciprocal = new double[] { 0, 0 };
            double[] normalEccentricity = new double[] { 0, 0 };
            double[] tangentialEccentricity = new double[] { 0, 0 };

            {// first particle
                massReciprocal[0] = Particles[ParticleID].Motion.IncludeTranslation() ? 1 / Particles[ParticleID].Mass : 0;
                normalMomentOfInertiaReciprocal[0] = Particles[ParticleID].Motion.IncludeRotation() ? Particles[ParticleID].CalculateSecondOrderEccentricity(NormalVector, ClosestPoints[ParticleID][secondObjectID]) / Particles[ParticleID].MomentOfInertia : 0;
                tangentialMomentOfInertiaReciprocal[0] = Particles[ParticleID].Motion.IncludeRotation() ? Particles[ParticleID].CalculateSecondOrderEccentricity(tangentialVector, ClosestPoints[ParticleID][secondObjectID]) / Particles[ParticleID].MomentOfInertia : 0;
                normalEccentricity[0] = Particles[ParticleID].CalculateEccentricity(NormalVector, ClosestPoints[ParticleID][secondObjectID]);
                tangentialEccentricity[0] = Particles[ParticleID].CalculateEccentricity(tangentialVector, ClosestPoints[ParticleID][secondObjectID]);
            }

            if (IsParticle(secondObjectID)) {
                massReciprocal[1] = Particles[secondObjectID].Motion.IncludeTranslation() ? 1 / Particles[secondObjectID].Mass : 0;
                normalMomentOfInertiaReciprocal[1] = Particles[secondObjectID].Motion.IncludeRotation() ? Particles[secondObjectID].CalculateSecondOrderEccentricity(NormalVector, ClosestPoints[secondObjectID][ParticleID]) / Particles[secondObjectID].MomentOfInertia : 0;
                tangentialMomentOfInertiaReciprocal[1] = Particles[secondObjectID].Motion.IncludeRotation() ? Particles[secondObjectID].CalculateSecondOrderEccentricity(tangentialVector, ClosestPoints[secondObjectID][ParticleID]) / Particles[secondObjectID].MomentOfInertia : 0;
                normalEccentricity[1] = Particles[secondObjectID].CalculateEccentricity(NormalVector, ClosestPoints[secondObjectID][ParticleID]);
                tangentialEccentricity[1] = Particles[secondObjectID].CalculateEccentricity(tangentialVector, ClosestPoints[secondObjectID][ParticleID]);
            }

            double collisionCoefficient = -(1 + CoefficientOfRestitution) * ((normalVelocity[0] - normalVelocity[1]) / (massReciprocal[0] + massReciprocal[1] + normalMomentOfInertiaReciprocal[0] + normalMomentOfInertiaReciprocal[1]));
            double tempRotVelocity2 = IsParticle(secondObjectID) ? TemporaryVelocity[secondObjectID][2] : 0;
            collisionCoefficient -= (1 + CoefficientOfRestitution) * ((normalEccentricity[0] * TemporaryVelocity[ParticleID][2] - normalEccentricity[1] * tempRotVelocity2) / (massReciprocal[0] + massReciprocal[1] + normalMomentOfInertiaReciprocal[0] + normalMomentOfInertiaReciprocal[1]));
            return collisionCoefficient;
        }

        /// <summary>
        /// Calculates the velocity of <paramref name="ParticleID"/> and <paramref name="SecondObjectID"/> based on <paramref name="CollisionCoefficient"/>.
        /// </summary>
        /// <param name="ParticleID"></param>
        /// <param name="SecondObjectID"></param>
        /// <param name="NormalVector"></param>
        /// <param name="CollisionCoefficient"></param>
        private void CalculateVelocityAfterCollision(int ParticleID, int SecondObjectID, Vector NormalVector, double CollisionCoefficient) {
            if (NormalVector.Dim != 2)
                throw new NotImplementedException("Physical collision model only implemented for 2D");

            Vector tangentialVector = new(-NormalVector[1], NormalVector[0]);
            Vector velocity = new(TemporaryVelocity[ParticleID][0], TemporaryVelocity[ParticleID][1]);
            double normalVelocity = velocity * NormalVector;
            double tangentialVelocity = velocity * tangentialVector;
            Vector velocityAfterCollision = Particles[ParticleID].Motion.IncludeTranslation()
                ? (normalVelocity + CollisionCoefficient / Particles[ParticleID].Mass) * NormalVector + tangentialVelocity * tangentialVector
                : new Vector(0, 0);
            TemporaryVelocity[ParticleID][0] = velocityAfterCollision[0];
            TemporaryVelocity[ParticleID][1] = velocityAfterCollision[1];

            Vector radialVector = Particles[ParticleID].CalculateRadialVector(ClosestPoints[ParticleID][SecondObjectID]);
            TemporaryVelocity[ParticleID][2] = Particles[ParticleID].Motion.IncludeRotation()
                ? TemporaryVelocity[ParticleID][2] + (radialVector[0] * NormalVector[1] - radialVector[1] * NormalVector[0]) * CollisionCoefficient / Particles[ParticleID].MomentOfInertia : 0;

            Particles[ParticleID].IsCollided = true;
            Overlapping[ParticleID][SecondObjectID] = false;
        }

        /// <summary>
        /// True = Particle, False = Wall
        /// </summary>
        /// <param name="objectID"></param>
        /// <returns></returns>
        private bool IsParticle(int objectID) {
            return objectID < Particles.Length;
        }

        /// <summary>
        /// True: Particles move either translational or rotational or both, False: Particle does not move at all.
        /// </summary>
        /// <param name="FirstParticleID"></param>
        /// <returns></returns>
        private bool IsMoving(int FirstParticleID) {
            return Particles[FirstParticleID].Motion.IncludeTranslation() || Particles[FirstParticleID].Motion.IncludeRotation();
        }
    }
}

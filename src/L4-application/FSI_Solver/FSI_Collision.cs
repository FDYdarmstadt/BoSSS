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

using BoSSS.Application.FSI_Solver;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FSI_Solver {

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
        /// Reduced version of the collision model, used to check for periodic boundaries.
        /// </summary>
        /// <param name="gridLenghtscale">
        /// Characteristic lenght of the grid.
        /// </param>
        public ParticleCollision(double gridLenghtscale) {
            CoefficientOfRestitution = 0;
            TimestepSize = 0;
            GridLengthScale = gridLenghtscale;
        }

        /// <summary>
        /// Constructor for the collision model.
        /// </summary>
        /// <param name="gridLenghtscale">
        /// Characteristic lenght of the grid.
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
        public ParticleCollision(double gridLenghtscale, double coefficientOfRestitution, double dt, double[][] wallCoordinates, bool[] IsPeriodicBoundary, double minDistance) {
            CoefficientOfRestitution = coefficientOfRestitution;
            TimestepSize = dt;
            GridLengthScale = gridLenghtscale;
            WallCoordinates = wallCoordinates;
            MinDistance = minDistance;
            this.IsPeriodicBoundary = IsPeriodicBoundary;
        }

        private readonly FSI_Auxillary Aux = new FSI_Auxillary();

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

        /// <summary>
        /// Update collision forces between two arbitrary particles and add them to forces acting on the corresponding particle
        /// </summary>
        /// <param name="particles">
        /// List of all particles
        /// </param>
        public void CalculateCollision(Particle[] particles) {
            Particles = particles;
            // Step 1
            // Some var definintion
            // =======================================================
            double saveTimestep = 0;
            int ParticleOffset = Particles.Length;
            double distanceThreshold = GridLengthScale / 10;
            if (MinDistance != 0)
                distanceThreshold = MinDistance;
            // Step 2
            // Loop over time until the particles hit.
            // =======================================================
            while (AccumulatedCollisionTimestep < TimestepSize)// the collision needs to take place within the current timestep dt.
            {
                CreateCollisionArrarys(Particles.Length);
                double globalMinimalDistance = double.MaxValue;

                // Step 2.1
                // Loop over the distance until a predefined criterion is 
                // met.
                // -------------------------------------------------------
                while (globalMinimalDistance > distanceThreshold) {
                    // Step 2.1.1
                    // Move the particle with the current save timestep.
                    // -------------------------------------------------------
                    MoveParticlesWithSaveTimestep(Particles, saveTimestep);
                    saveTimestep = double.MaxValue;
                    for (int p0 = 0; p0 < Particles.Length; p0++) {
                        // Step 2.1.2
                        // Test for wall collisions for all particles 
                        // of the current color.
                        // -------------------------------------------------------
                        Vector[] nearFieldWallPoints = GetNearFieldWall(Particles[p0]);
                        for (int w = 0; w < nearFieldWallPoints.Length; w++) {
                            Particles[p0].ClosestPointOnOtherObjectToThis = new Vector(Particles[p0].Motion.GetPosition(0));
                            if (nearFieldWallPoints[w].IsNullOrEmpty())
                                continue;
                            else
                                particles[p0].ClosestPointOnOtherObjectToThis = new Vector(nearFieldWallPoints[w]);
                            CalculateMinimumDistance(Particles[p0], out Vector temp_DistanceVector, out Vector temp_ClosestPoint_p0, out bool temp_Overlapping);
                            int wallID = ParticleOffset + w;
                            DistanceVector[p0][wallID] = new Vector(temp_DistanceVector);
                            ClosestPoints[p0][wallID] = new Vector(temp_ClosestPoint_p0);
                            double temp_SaveTimeStep = DynamicTimestep(p0, ParticleOffset + w);
                            AccumulatedLocalSaveTimestep[p0][wallID] += temp_SaveTimeStep;
                            if (temp_SaveTimeStep < saveTimestep && temp_SaveTimeStep > 0) {
                                saveTimestep = temp_SaveTimeStep;
                            }
                            if (DistanceVector[p0][wallID].Abs() < globalMinimalDistance) {
                                globalMinimalDistance = DistanceVector[p0][wallID].Abs();
                            }
                            if (temp_Overlapping) {
                                double overlappingTimestep = -TimestepSize * 0.5; // reset time to find a particle state before they overlap.
                                saveTimestep = 0;
                                Particle[] overlappingParticles = new Particle[] { Particles[p0] };
                                Console.WriteLine("Particle " + p0 + " and wall " + w + " overlap");
                                MoveParticlesWithSaveTimestep(overlappingParticles, overlappingTimestep);
                            }
                        }

                        // Step 2.1.3
                        // Test for particle-particle collisions for all particles 
                        // of the current color.
                        // -------------------------------------------------------
                        for (int p1 = p0 + 1; p1 < Particles.Length; p1++) {
                            Particle[] currentParticles = new Particle[] { Particles[p0], Particles[p1] };
                            CalculateMinimumDistance(currentParticles, out Vector temp_DistanceVector, out Vector[] temp_ClosestPoints, out bool temp_Overlapping);
                            Overlapping[p0][p1] = temp_Overlapping;
                            Overlapping[p1][p0] = temp_Overlapping;
                            ClosestPoints[p0][p1] = temp_ClosestPoints[0];
                            ClosestPoints[p1][p0] = temp_ClosestPoints[1];
                            DistanceVector[p0][p1] = new Vector(temp_DistanceVector);
                            temp_DistanceVector.Scale(-1);
                            DistanceVector[p1][p0] = new Vector(temp_DistanceVector);
                            double temp_SaveTimeStep = DynamicTimestep(p0, p1);
                            AccumulatedLocalSaveTimestep[p0][p1] += temp_SaveTimeStep;
                            AccumulatedLocalSaveTimestep[p1][p0] += temp_SaveTimeStep;
                            if (temp_SaveTimeStep < saveTimestep && temp_SaveTimeStep > 0) {
                                saveTimestep = temp_SaveTimeStep;
                            }
                            if (DistanceVector[p0][p1].Abs() < globalMinimalDistance) {
                                globalMinimalDistance = DistanceVector[p0][p1].Abs();
                            }
                            if (temp_Overlapping) {
                                double overlappingTimestep = -TimestepSize * 0.5; // reset time to find a particle state before they overlap.
                                saveTimestep = 0;
                                Particle[] overlappingParticles = new Particle[] { Particles[p0], Particles[p1] };
                                Console.WriteLine("Particle " + p0 + " and particle " + p1 + " overlap");
                                MoveParticlesWithSaveTimestep(overlappingParticles, overlappingTimestep);
                            }
                        }
                    }
                    if (saveTimestep <= double.MaxValue && globalMinimalDistance != double.MaxValue) {
                        Console.WriteLine("Minimal distance " + globalMinimalDistance + ", threshold " + distanceThreshold + ", current save time+step " + saveTimestep + " accumulated timestep " + AccumulatedCollisionTimestep);
                    }
                    if (saveTimestep >= 0)
                        AccumulatedCollisionTimestep += saveTimestep;
                    if (AccumulatedCollisionTimestep >= TimestepSize)
                        break;
                }
                if (AccumulatedCollisionTimestep == double.MaxValue)
                    break;

                ParticleCollidedWith = CreateArrayWithCollidedParticles(distanceThreshold);
                CollisionCluster = ClusterCollisionsContainingSameParticles();

                for (int c = 0; c < CollisionCluster.Count(); c++) {
                    for (int i = 0; i < CollisionCluster[c].Count(); i++) {
                        int currentParticleID = CollisionCluster[c][i];
                        if (IsParticle(currentParticleID)) {
                            for (int j = 0; j < ParticleCollidedWith[currentParticleID].Count(); j++) {
                                int secondObjectID = ParticleCollidedWith[currentParticleID][j];
                                ComputeMomentumBalanceCollision(currentParticleID, secondObjectID, distanceThreshold);
                                TransferResultsToGhostParticles(currentParticleID);
                                if(IsParticle(secondObjectID))
                                    TransferResultsToGhostParticles(secondObjectID);
                                if (Particles[currentParticleID].IsCollided)
                                    Console.WriteLine("Particle " + currentParticleID + " and particle / wall " + secondObjectID + " collided");
                            }
                        }
                    }
                }

                for (int p = 0; p < Particles.Length; p++) {
                    if (particles[p].IsCollided) {
                        particles[p].Motion.InitializeParticleVelocity(new double[] { TemporaryVelocity[p][0], TemporaryVelocity[p][1] }, TemporaryVelocity[p][2]);
                        particles[p].Motion.InitializeParticleAcceleration(new double[] { 0, 0 }, 0);
                    }
                    particles[p].Motion.SetCollisionTimestep(AccumulatedCollisionTimestep - saveTimestep);
                    CollisionCluster.Clear();
                    ParticleCollidedWith.Clear();
                }
            }
        }

        private void TransferResultsToGhostParticles(int currentParticleID) {
            if (Particles[currentParticleID].MasterGhostIDs[0] > 0 && Particles[currentParticleID].IsCollided) {
                for (int k = 0; k < Particles[currentParticleID].MasterGhostIDs.Length; k++) {
                    if(Particles[currentParticleID].MasterGhostIDs[k] - 1 != currentParticleID && Particles[currentParticleID].MasterGhostIDs[k] > 0) {
                        Particles[Particles[currentParticleID].MasterGhostIDs[k] - 1].IsCollided = true;
                        TemporaryVelocity[Particles[currentParticleID].MasterGhostIDs[k] - 1] = TemporaryVelocity[currentParticleID].CloneAs();
                    }

                }
            }
        }

        private int[][] ClusterCollisionsContainingSameParticles() {
            List<int[]> globalParticleCluster = new List<int[]>();
            bool[] partOfCollisionCluster = new bool[Particles.Length];
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

        private int[][] CreateArrayWithCollidedParticles(double distanceThreshold) {
            int ParticleOffset = Particles.Length;
            int[][] particleCollidedWith = new int[Particles.Length][];
            for (int p0 = 0; p0 < Particles.Length; p0++) {
                List<int> currentParticleCollidedWith = new List<int>();
                for (int w = 0; w < 4; w++) {
                    FindCollisionPartners(p0, ParticleOffset + w, currentParticleCollidedWith, distanceThreshold);
                }
                for (int p1 = p0 + 1; p1 < Particles.Length; p1++) {
                    FindCollisionPartners(p0, p1, currentParticleCollidedWith, distanceThreshold);
                }
                particleCollidedWith[p0] = currentParticleCollidedWith.ToArray();
            }
            return particleCollidedWith;
        }

        private void FindCollisionPartners(int particle, int potentialCollisionPartner, List<int> currentParticleCollidedWith, double distanceThreshold) {
            if (DistanceVector[particle][potentialCollisionPartner].Abs() <= distanceThreshold && AccumulatedCollisionTimestep < TimestepSize && AccumulatedLocalSaveTimestep[particle][potentialCollisionPartner] > 0) {
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
            if (Particles[FirstParticleID].Motion.IncludeTranslation || Particles[FirstParticleID].Motion.IncludeRotation) 
                detectCollisionVn_P0 = CalculateNormalSurfaceVelocity(FirstParticleID, normalVector, ClosestPoints[FirstParticleID][SecondParticleID]);
            if (IsParticle(SecondParticleID)) {
                if (Particles[FirstParticleID].Motion.IncludeTranslation || Particles[FirstParticleID].Motion.IncludeRotation)
                    detectCollisionVn_P1 = CalculateNormalSurfaceVelocity(SecondParticleID, normalVector, ClosestPoints[SecondParticleID][FirstParticleID]);
            }
            return (detectCollisionVn_P1 - detectCollisionVn_P0 == 0) ? double.MaxValue : 0.1 * distance / (detectCollisionVn_P1 - detectCollisionVn_P0);
        }

        private double CalculateNormalSurfaceVelocity(int particleID, Vector normalVector, Vector closestPoint) {
            Vector pointVelocity = new Vector(2);
            Vector radialVector = Particles[particleID].CalculateRadialVector(closestPoint);
            pointVelocity[0] = TemporaryVelocity[particleID][0] - TemporaryVelocity[particleID][2] * radialVector[1];
            pointVelocity[1] = TemporaryVelocity[particleID][1] + TemporaryVelocity[particleID][2] * radialVector[0];
            return pointVelocity * normalVector;
        }

        private void MoveParticlesWithSaveTimestep(Particle[] particles, double dynamicTimestep) {
            for (int p = 0; p < particles.Length; p++) {
                Particle currentParticle = particles[p];
                if (dynamicTimestep != 0) {
                    currentParticle.Motion.CollisionParticlePositionAndAngle(dynamicTimestep);
                }
            }
        }

        private void CalculateMinimumDistance(Particle[] Particles, out Vector DistanceVector, out Vector[] ClosestPoints, out bool Overlapping) {
            int spatialDim = Particles[0].Motion.GetPosition(0).Dim;
            double distance = double.MaxValue;
            DistanceVector = new Vector(spatialDim);
            ClosestPoints = new Vector[2];
            ClosestPoints[0] = new Vector(spatialDim);
            ClosestPoints[1] = new Vector(spatialDim);
            Overlapping = false;
            int NoOfSubParticles1 = Particles[1] == null ? 1 : Particles[1].NoOfSubParticles;
;
            for (int i = 0; i < Particles[0].NoOfSubParticles; i++) {
                for (int j = 0; j < NoOfSubParticles1; j++) {
                    GJK_DistanceAlgorithm(Particles[0], i, Particles[1], j, out Vector temp_DistanceVector, out Vector[] temp_ClosestPoints, out Overlapping);
                    if (Overlapping)
                        break;
                    if (temp_DistanceVector.Abs() < distance) {
                        distance = temp_DistanceVector.Abs();
                        DistanceVector = new Vector(temp_DistanceVector);
                        ClosestPoints = temp_ClosestPoints.CloneAs();
                    }
                }
            }
        }

        /// <summary>
        /// Computes the minimal distance between a particle and the wall.
        /// </summary>
        /// <param name="particle">
        /// The first particle.
        /// </param>
        /// <param name="DistanceVector">
        /// The vector of the minimal distance between the two objects.
        /// </param>
        /// <param name="ClosestPoint">
        /// The point on the first object closest to the second one.
        /// </param>
        /// <param name="Overlapping">
        /// Is true if the two particles are overlapping.
        /// </param>
        internal void CalculateMinimumDistance(Particle particle, out Vector DistanceVector, out Vector ClosestPoint, out bool Overlapping) {
            int spatialDim = particle.Motion.GetPosition(0).Dim;
            double distance = double.MaxValue;
            DistanceVector = new Vector(spatialDim);
            ClosestPoint = new Vector(spatialDim);
            Overlapping = false;

            for (int i = 0; i < particle.NoOfSubParticles; i++) {
                GJK_DistanceAlgorithm(particle, i, null, 1, out Vector temp_DistanceVector, out Vector[] temp_ClosestPoints, out Overlapping);
                if (Overlapping)
                    break;
                if (temp_DistanceVector.Abs() < distance) {
                    distance = temp_DistanceVector.Abs();
                    DistanceVector = new Vector(temp_DistanceVector);
                    ClosestPoint = new Vector(temp_ClosestPoints[0]);
                }
            }
        }

        /// <summary>
        /// Computes the distance between two objects (particles or walls). Algorithm based on
        /// E.G.Gilbert, D.W.Johnson, S.S.Keerthi.
        /// </summary>
        /// <param name="Particle0">
        /// The first particle.
        /// </param>
        /// <param name="SubParticleID0">
        /// In case of concave particles the particle is devided into multiple convex subparticles. Each of them has its one ID and needs to be tested as if it was a complete particle.
        /// </param>
        ///  <param name="Particle1">
        /// The second particle, if Particle1 == null it is assumed to be a wall.
        /// </param>
        /// <param name="SubParticleID1">
        /// In case of concave particles the particle is devided into multiple convex subparticles. Each of them has its one ID and needs to be tested as if it was a complete particle.
        /// </param>
        /// <param name="DistanceVec">
        /// The vector of the minimal distance between the two objects.
        /// </param>
        /// <param name="closestPoints">
        /// The point on one object closest to the other one.
        /// </param>
        /// <param name="Overlapping">
        /// Is true if the two particles are overlapping.
        /// </param>
        private void GJK_DistanceAlgorithm(Particle Particle0, int SubParticleID0, Particle Particle1, int SubParticleID1, out Vector DistanceVec, out Vector[] closestPoints, out bool Overlapping) {

            // Step 1
            // Initialize the algorithm with the particle position
            // =======================================================
            int spatialDim = Particle0.Motion.GetPosition(0).Dim;
            Vector[] positionVectors = new Vector[2];
            positionVectors[0] = new Vector(Particle0.Motion.GetPosition(0));
            positionVectors[1] = new Vector(Particle1 == null ? Particle0.ClosestPointOnOtherObjectToThis : Particle1.Motion.GetPosition(0));
            Vector supportVector = positionVectors[0] - positionVectors[1];
            Aux.TestArithmeticException(supportVector, "support vector");

            // Define the simplex, which contains all points to be tested for their distance (max. 3 points in 2D)
            List<Vector> Simplex = new List<Vector> { new Vector(supportVector) };

            closestPoints = new Vector[2];
            closestPoints[0] = new Vector(spatialDim);
            closestPoints[1] = new Vector(spatialDim);
            Overlapping = false;
            int maxNoOfIterations = 50;

            // Step 2
            // Start the iteration
            // =======================================================
            for (int i = 0; i <= maxNoOfIterations; i++) {
                Vector negativeSupportVector = new Vector(spatialDim);
                negativeSupportVector.Sub(supportVector);

                // Calculate the support point of the two particles, 
                // which are the closest points if the algorithm is finished.
                // -------------------------------------------------------
                CalculateSupportPoint(Particle0, SubParticleID0, negativeSupportVector, out closestPoints[0]);
                // Particle-Particle collision
                if (Particle1 != null) {
                    CalculateSupportPoint(Particle1, SubParticleID1, supportVector, out closestPoints[1]);
                }
                // Particle-wall collision
                else {
                    closestPoints[1] = new Vector(spatialDim);
                    if (positionVectors[0][0] == positionVectors[1][0])
                        closestPoints[1] = new Vector(closestPoints[0][0], positionVectors[1][1]);
                    else
                        closestPoints[1] = new Vector(positionVectors[1][0], closestPoints[0][1]);
                }
                Aux.TestArithmeticException(closestPoints[0], "closest point on particle 0");
                Aux.TestArithmeticException(closestPoints[1], "closest point on particle 1");

                // The current support point can be found by forming 
                // the difference of the support points on the two particles
                // -------------------------------------------------------
                Vector supportPoint = closestPoints[0] - closestPoints[1];
                Aux.TestArithmeticException(supportPoint, "support point");

                // If the condition is true
                // we have found the closest points!
                // -------------------------------------------------------
                if (((supportVector * negativeSupportVector) - (supportPoint * negativeSupportVector)) >= -1e-12 && i > 1)
                    break;

                // Add new support point to simplex
                // -------------------------------------------------------
                Simplex.Insert(0, new Vector(supportPoint));

                // Calculation the new vector v with the distance
                // algorithm
                // -------------------------------------------------------
                supportVector = DistanceAlgorithm(Simplex, out Overlapping);

                // End algorithm if the two objects are overlapping.
                // -------------------------------------------------------
                if (Overlapping)
                    break;

                // Could not find the closest points... crash!
                // -------------------------------------------------------
                if (i == maxNoOfIterations)
                    throw new Exception("No convergence in GJK-algorithm, reached iteration #" + i);
            }

            // Step 3
            // Return min distance and distance vector.
            // =======================================================
            DistanceVec = new Vector(supportVector);
        }

        private void CalculateSupportPoint(Particle particle, int SubParticleID, Vector supportVector, out Vector supportPoint) {
            int spatialDim = particle.Motion.GetPosition(0).Dim;
            supportPoint = new Vector(spatialDim);
            // A direct formulation of the support function for a sphere exists, thus it is also possible to map it to an ellipsoid.
            if (particle is Particle_Ellipsoid || particle is Particle_Sphere || particle is Particle_Rectangle || particle is Particle_Shell) {
                supportPoint = particle.GetSupportPoint(supportVector, SubParticleID);
            }
            // Interpolated binary search in all other cases.
            else {
                double angle = particle.Motion.GetAngle(0);
                Vector particleDirection = new Vector(Math.Cos(angle), Math.Sin(angle));
                double crossProductDirectionSupportVector = particleDirection[0] * supportVector[1] - particleDirection[1] * supportVector[0];
                double searchStartAngle = (1 - Math.Sign(crossProductDirectionSupportVector)) * Math.PI / 2 + Math.Acos((supportVector * particleDirection) / supportVector.L2Norm());
                double L = searchStartAngle - Math.PI;
                double R = searchStartAngle + Math.PI;
                while (L < R && Math.Abs(L-R) > 1e-15) {
                    searchStartAngle = (L + R) / 2;
                    double dAngle = 1e-8;
                    MultidimensionalArray SurfacePoints = particle.GetSurfacePoints(dAngle, searchStartAngle, SubParticleID);
                    Vector RightNeighbour = new Vector(spatialDim);
                    Vector LeftNeighbour = new Vector(spatialDim);
                    for (int d = 0; d < spatialDim; d++) {
                        supportPoint[d] = SurfacePoints[1, d];
                        LeftNeighbour[d] = SurfacePoints[0, d];
                        RightNeighbour[d] = SurfacePoints[2, d];
                    }
                    if ((supportPoint * supportVector) > (RightNeighbour * supportVector) && (supportPoint * supportVector) > (LeftNeighbour * supportVector))
                        break; // The current temp_supportPoint is the actual support point.
                    else if ((RightNeighbour * supportVector) > (LeftNeighbour * supportVector))
                        L = searchStartAngle; // Search on the right side of the current point.
                    else
                        R = searchStartAngle; // Search on the left side.
                }
                Vector position = new Vector(particle.Motion.GetPosition(0));
                supportPoint.Acc(position);
            }
        }

        private Vector DistanceAlgorithm(List<Vector> simplex, out bool overlapping) {
            Vector supportVector = new Vector(simplex[0].Dim);
            overlapping = false;

            // Step 1
            // Test for multiple Simplex-points 
            // and remove the duplicates
            // =======================================================
            for (int s1 = 0; s1 < simplex.Count(); s1++) {
                for (int s2 = s1 + 1; s2 < simplex.Count(); s2++) {
                    if ((simplex[s1] - simplex[s2]).Abs() < 1e-8) {
                        simplex.RemoveAt(s2);
                    }
                }
            }

            // Step 2
            // Calculate dot product between all position vectors and 
            // save to an 2D-array.
            // =======================================================
            double[][] dotProductSimplex = new double[simplex.Count()][];
            for (int s1 = 0; s1 < simplex.Count(); s1++) {
                dotProductSimplex[s1] = new double[simplex.Count()];
                for (int s2 = s1; s2 < simplex.Count(); s2++) {
                    dotProductSimplex[s1][s2] = simplex[s1] * simplex[s2];
                }
            }

            // Step 3
            // Main routine to determine the relatve position of
            // the simplex towards the origin.
            // =======================================================
            // The simplex contains only one element, which must be
            // the closest point of this simplex to the origin
            // -------------------------------------------------------
            if (simplex.Count() == 1) {
                supportVector = new Vector(simplex[0]);
                Aux.TestArithmeticException(supportVector, "support vector");
            }

            // The simplex contains two elements, lets test which is
            // closest to the origin
            // -------------------------------------------------------
            else if (simplex.Count() == 2) {
                // One of the simplex point is closest to the origin, 
                // choose this and delete the other one.
                // -------------------------------------------------------
                bool continueAlgorithm = true;
                for (int s = 0; s < simplex.Count(); s++) {
                    if (dotProductSimplex[s][s] - dotProductSimplex[0][1] <= 0) {
                        supportVector = new Vector(simplex[s]);
                        simplex.RemoveAt(Math.Abs(s - 1));
                        Aux.TestArithmeticException(supportVector, "support vector");
                        continueAlgorithm = false;
                        break;
                    }
                }
                // A point at the line between the two simplex points is
                // closest to the origin, thus we need to keep both points.
                // -------------------------------------------------------
                if (continueAlgorithm) {
                    Vector simplexDistanceVector = simplex[1] - simplex[0];
                    double lambda = Math.Abs(simplex[1].CrossProduct2D(simplexDistanceVector)) / simplexDistanceVector.AbsSquare();
                    if (lambda == 0) // if the origin lies on the line between the two simplex points, the two objects are overlapping in one point
                        overlapping = true;
                    supportVector[0] = -lambda * simplexDistanceVector[1];
                    supportVector[1] = lambda * simplexDistanceVector[0];
                    Aux.TestArithmeticException(supportVector, "support vector");
                }
            }

            // The simplex contains three elements, lets test which is
            // closest to the origin
            // -------------------------------------------------------
            else if (simplex.Count() == 3) {
                bool continueAlgorithm = true;
                // Test whether one of the simplex points is closest to 
                // the origin
                // -------------------------------------------------------
                for (int s1 = 0; s1 < simplex.Count(); s1++) {
                    int s2 = s1 == 2 ? 2 : 1;
                    int s3 = s1 == 0 ? 0 : 1;
                    if (dotProductSimplex[s1][s1] - dotProductSimplex[0][s2] <= 0 && dotProductSimplex[s1][s1] - dotProductSimplex[s3][2] <= 0) {
                        supportVector = new Vector(simplex[s1]);
                        // Delete the complete simplex and add back the point closest to the origin
                        simplex.Clear();
                        simplex.Add(new Vector(supportVector));
                        continueAlgorithm = false;
                        Aux.TestArithmeticException(supportVector, "support vector");
                        break;
                    }
                }
                // None of the simplex points was the closest point, 
                // thus, it has to be any point at the edges
                // -------------------------------------------------------
                if (continueAlgorithm) {
                    for (int s1 = simplex.Count() - 1; s1 >= 0; s1--) {
                        int s2 = s1 == 0 ? 1 : 2;
                        int s3 = s1 == 2 ? 1 : 0;
                        // Calculate a crossproduct of the form (BC x BA) x BA * BX
                        double crossProduct = new double();
                        switch (s1) {
                            case 0:
                                double temp1 = dotProductSimplex[1][2] - dotProductSimplex[0][2] - dotProductSimplex[1][1] + dotProductSimplex[0][1];
                                double temp2 = dotProductSimplex[0][1] - dotProductSimplex[0][0] - dotProductSimplex[1][2] + dotProductSimplex[0][2];
                                double temp3 = dotProductSimplex[1][1] - 2 * dotProductSimplex[0][1] + dotProductSimplex[0][0];
                                crossProduct = dotProductSimplex[0][1] * temp1 + dotProductSimplex[1][1] * temp2 + dotProductSimplex[1][2] * temp3;
                                break;
                            case 1:
                                temp1 = -dotProductSimplex[2][2] + dotProductSimplex[0][2] + dotProductSimplex[1][2] - dotProductSimplex[0][1];
                                temp2 = dotProductSimplex[2][2] - 2 * dotProductSimplex[0][2] + dotProductSimplex[0][0];
                                temp3 = dotProductSimplex[0][2] - dotProductSimplex[0][0] - dotProductSimplex[1][2] + dotProductSimplex[0][1];
                                crossProduct = dotProductSimplex[0][2] * temp1 + dotProductSimplex[1][2] * temp2 + dotProductSimplex[2][2] * temp3;
                                break;
                            case 2:
                                temp1 = dotProductSimplex[2][2] - 2 * dotProductSimplex[1][2] + dotProductSimplex[1][1];
                                temp2 = -dotProductSimplex[2][2] + dotProductSimplex[1][2] + dotProductSimplex[0][2] - dotProductSimplex[0][1];
                                temp3 = dotProductSimplex[1][2] - dotProductSimplex[1][1] - dotProductSimplex[0][2] + dotProductSimplex[0][1];
                                crossProduct = dotProductSimplex[0][2] * temp1 + dotProductSimplex[1][2] * temp2 + dotProductSimplex[2][2] * temp3;
                                break;
                        }
                        // A point on one of the edges is closest to the origin.
                        if (dotProductSimplex[s3][s3] - dotProductSimplex[s3][s2] >= 0 && dotProductSimplex[s2][s2] - dotProductSimplex[s3][s2] >= 0 && crossProduct >= 0 && continueAlgorithm) {
                            Vector simplexDistanceVector = simplex[s2] - simplex[s3];
                            double Lambda = Math.Abs(simplex[s2].CrossProduct2D(simplexDistanceVector)) / simplexDistanceVector.AbsSquare();
                            supportVector[0] = -Lambda * simplexDistanceVector[1];
                            supportVector[1] = Lambda * simplexDistanceVector[0];
                            // save the two remaining simplex points and clear the simplex.
                            Vector tempSimplex1 = new Vector(simplex[s2]);
                            Vector tempSimplex2 = new Vector(simplex[s3]);
                            simplex.Clear();
                            // Re-add the remaining points
                            simplex.Add(tempSimplex1);
                            simplex.Add(tempSimplex2);
                            continueAlgorithm = false;
                            Aux.TestArithmeticException(supportVector, "support vector");
                            break;
                        }
                    }
                }
                // None of the conditions above are true, 
                // thus, the simplex must contain the origin and 
                // the two particles do overlap.
                // -------------------------------------------------------
                if (continueAlgorithm)
                    overlapping = true;
            }
            return supportVector;
        }

        private void ComputeMomentumBalanceCollision(int particleID, int secondObjectID, double threshold) {
            double distance = DistanceVector[particleID][secondObjectID].Abs();
            Vector normalVector = DistanceVector[particleID][secondObjectID];
            normalVector.Normalize();
            if (distance > threshold && !Overlapping[particleID][secondObjectID])
                return;
            if (Overlapping[particleID][secondObjectID])
                Console.WriteLine("Overlapping " + particleID + " " + secondObjectID);

            double detectCollisionVn_P0 = 0;
            double detectCollisionVn_P1 = 0;
            if (Particles[0].Motion.IncludeTranslation || Particles[0].Motion.IncludeRotation) 
                detectCollisionVn_P0 = CalculateNormalSurfaceVelocity(particleID, normalVector, ClosestPoints[particleID][secondObjectID]);
            if (IsParticle(secondObjectID)) {
                if (Particles[1].Motion.IncludeTranslation || Particles[1].Motion.IncludeRotation)
                    detectCollisionVn_P1 = CalculateNormalSurfaceVelocity(secondObjectID, normalVector, ClosestPoints[secondObjectID][particleID]);
            }
            if (detectCollisionVn_P1 - detectCollisionVn_P0 <= 0) 
                return;

            Vector tangentialVector = new Vector(-normalVector[1], normalVector[0]);
            Particles[particleID].CalculateEccentricity(normalVector, ClosestPoints[particleID][secondObjectID]);
            if (IsParticle(secondObjectID))
                Particles[secondObjectID].CalculateEccentricity(normalVector, ClosestPoints[secondObjectID][particleID]);
            double collisionCoefficient = CalculateCollisionCoefficient(particleID, secondObjectID, normalVector);
            Vector velocityP0 = CalculateNormalAndTangentialVelocity(particleID, normalVector);
            Vector radialVectorP0 = Particles[particleID].CalculateRadialVector(ClosestPoints[particleID][secondObjectID]);
            Vector tempVel0 = Particles[particleID].Motion.IncludeTranslation 
                ? (velocityP0[0] + collisionCoefficient / Particles[particleID].Motion.ParticleMass) * normalVector + CoefficientOfRestitution * velocityP0[1] * tangentialVector 
                : new Vector(0, 0);
            TemporaryVelocity[particleID][0] = tempVel0[0];
            TemporaryVelocity[particleID][1] = tempVel0[1];
            TemporaryVelocity[particleID][2] = Particles[particleID].Motion.IncludeRotation ? TemporaryVelocity[particleID][2] + (radialVectorP0[0] * normalVector[1] - radialVectorP0[1] * normalVector[0]) * collisionCoefficient / Particles[particleID].MomentOfInertia : 0;
            Particles[particleID].IsCollided = true;
            Overlapping[particleID][secondObjectID] = false;

            if (IsParticle(secondObjectID)) {
                Vector velocityP1 = CalculateNormalAndTangentialVelocity(secondObjectID, normalVector);
                Vector tempVel1 = Particles[secondObjectID].Motion.IncludeTranslation
                    ? (velocityP1[0] - collisionCoefficient / Particles[secondObjectID].Motion.ParticleMass) * normalVector + CoefficientOfRestitution * velocityP1[1] * tangentialVector
                    : new Vector(0, 0);
                Vector radialVectorP1 = Particles[secondObjectID].CalculateRadialVector(ClosestPoints[secondObjectID][particleID]);
                TemporaryVelocity[secondObjectID][0] = tempVel1[0];
                TemporaryVelocity[secondObjectID][1] = tempVel1[1];
                TemporaryVelocity[secondObjectID][2] = Particles[secondObjectID].Motion.IncludeRotation ? TemporaryVelocity[secondObjectID][2] - (radialVectorP1[0] * normalVector[1] - radialVectorP1[1] * normalVector[0]) * collisionCoefficient / Particles[secondObjectID].MomentOfInertia : 0;
                Overlapping[secondObjectID][particleID] = false;
                Particles[secondObjectID].IsCollided = true;
            }
        }

        private Vector CalculateNormalAndTangentialVelocity(int particleID, Vector normalVector) {
            Vector tangentialVector = new Vector(-normalVector[1], normalVector[0]);
            Vector velocity = new Vector(TemporaryVelocity[particleID][0], TemporaryVelocity[particleID][1]);
            return new Vector(velocity * normalVector, velocity * tangentialVector);
        }

        private double CalculateCollisionCoefficient(int p0, int p1, Vector normalVector) {
            Vector[] translationalVelocity = new Vector[2];
            translationalVelocity[0] = CalculateNormalAndTangentialVelocity(p0, normalVector);
            translationalVelocity[1] = new Vector(0, 0);
            double[] massReciprocal = new double[] { 0, 0 };
            double[] momentOfInertiaReciprocal = new double[] { 0, 0 };
            double[] eccentricity = new double[] { 0, 0 };

            massReciprocal[0] = Particles[p0].Motion.IncludeTranslation ? 1 / Particles[p0].Motion.ParticleMass : 0;
            momentOfInertiaReciprocal[0] = Particles[p0].Motion.IncludeRotation ? Particles[p0].CalculateSecondOrderEccentricity(normalVector, ClosestPoints[p0][p1]) / Particles[p0].MomentOfInertia : 0;
            eccentricity[0] = Particles[p0].CalculateEccentricity(normalVector, ClosestPoints[p0][p1]);

            if (IsParticle(p1)) {
                translationalVelocity[1] = CalculateNormalAndTangentialVelocity(p1, normalVector);
                massReciprocal[1] = Particles[p1].Motion.IncludeTranslation ? 1 / Particles[p1].Motion.ParticleMass : 0;
                momentOfInertiaReciprocal[1] = Particles[p1].Motion.IncludeRotation ? Particles[p1].CalculateSecondOrderEccentricity(normalVector, ClosestPoints[p1][p0]) / Particles[p1].MomentOfInertia : 0;
                eccentricity[1] = Particles[p1].CalculateEccentricity(normalVector, ClosestPoints[p1][p0]);
            }

            double collisionCoefficient = -(1 + CoefficientOfRestitution) * ((translationalVelocity[0][0] - translationalVelocity[1][0]) / (massReciprocal[0] + massReciprocal[1] + momentOfInertiaReciprocal[0] + momentOfInertiaReciprocal[1]));
            double tempRotVelocity2 = 0;
            if (IsParticle(p1))
                tempRotVelocity2 = TemporaryVelocity[p1][2];
            collisionCoefficient -= (1 + CoefficientOfRestitution) * ((eccentricity[0] * TemporaryVelocity[p0][2] - eccentricity[1] * tempRotVelocity2) / (massReciprocal[0] + massReciprocal[1] + momentOfInertiaReciprocal[0] + momentOfInertiaReciprocal[1]));
            return collisionCoefficient;
        }

        private Vector[] GetNearFieldWall(Particle particle) {
            Vector[] nearFieldWallPoint = new Vector[4];
            Vector particlePosition = particle.Motion.GetPosition(0);
            double particleMaxLengthscale = particle.GetLengthScales().Max();
            for (int w0 = 0; w0 < WallCoordinates.Length; w0++) {
                for(int w1 = 0; w1 < WallCoordinates[w0].Length; w1++) {
                    if (WallCoordinates[w0][w1] == 0 || IsPeriodicBoundary[w0])
                        continue;
                    double minDistance = Math.Abs(particlePosition[w0] - WallCoordinates[w0][w1]) - particleMaxLengthscale;
                    if(minDistance < 5 * GridLengthScale || minDistance < 1e-2)
                    {
                        if (w0 == 0)
                            nearFieldWallPoint[w0 + w1] = new Vector(WallCoordinates[0][w1], particlePosition[1]);
                        else
                            nearFieldWallPoint[w0 * 2 + w1] = new Vector(particlePosition[0], WallCoordinates[1][w1]);
                    }
                }
            }
            return nearFieldWallPoint;
        }

        private bool IsParticle(int objectID) {
            return objectID < Particles.Length;
        }
    }
}

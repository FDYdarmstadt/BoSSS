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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Platform.LinAlg;

namespace FSI_Solver {
    class FSI_Collision {
        private readonly double m_dt;
        private readonly double m_hMin;
        private readonly int m_CurrentColor;
        private readonly double m_CoefficientOfRestitution;
        private readonly LevelSetTracker m_LevelSetTracker;

        private double AccDynamicTimestep = 0;

        private double[][] SaveTimeStepArray;
        private double[][] Distance;
        private Vector[][] DistanceVector;
        private Vector[][] ClosestPoints;

        public FSI_Collision(LevelSetTracker levelSetTracker, int currentColor, double CoefficientOfRestitution, double dt) {
            m_CoefficientOfRestitution = CoefficientOfRestitution;
            m_dt = dt;
            m_hMin = levelSetTracker.GridDat.Cells.h_minGlobal;
            m_CurrentColor = currentColor;
            m_LevelSetTracker = levelSetTracker;
        }

        private readonly FSI_Auxillary Aux = new FSI_Auxillary();

        private void CreateCollisionArrarys(int noOfParticles) {
            SaveTimeStepArray = new double[noOfParticles][];
            Distance = new double[noOfParticles][];
            DistanceVector = new Vector[noOfParticles][];
            ClosestPoints = new Vector[noOfParticles][];
            for (int p = 0; p < noOfParticles; p++) {
                SaveTimeStepArray[p] = new double[noOfParticles + 4];
                Distance[p] = new double[noOfParticles + 4];
                DistanceVector[p] = new Vector[noOfParticles + 4];
                ClosestPoints[p] = new Vector[noOfParticles + 4];
            }
        }

        /// <summary>
        /// Update collision forces between two arbitrary particles and add them to forces acting on the corresponding particle
        /// </summary>
        /// <param name="particles">
        /// List of all particles
        /// </param>
        public void CalculateCollision(Particle[] particles, IGridData gridData, int[] cellColor, LevelSetTracker LsTrk) {
            // Step 1
            // Some var definintion
            // =======================================================
            int ParticleOffset = particles.Length;
            double distanceThreshold = m_hMin / 5;

            // Step 2
            // Loop over time until the particles hit.
            // =======================================================
            while (AccDynamicTimestep < m_dt)// the collision needs to take place within the current timestep dt.
            {
                CreateCollisionArrarys(particles.Length);
                double minimalDistance = double.MaxValue;
                double SaveTimeStep = 0;// the timestep size without any collision

                // Step 2.1
                // Loop over the distance until a predefined criterion is 
                // met.
                // -------------------------------------------------------
                while (minimalDistance > distanceThreshold) {
                    // Step 2.1.1
                    // Move the particle with the current save timestep.
                    // -------------------------------------------------------
                    UpdateParticleState(particles, SaveTimeStep);
                    SaveTimeStep = double.MaxValue;
                    for (int p0 = 0; p0 < particles.Length; p0++) {
                        // Step 2.1.2
                        // Test for wall collisions for all particles 
                        // of the current color.
                        // -------------------------------------------------------
                        CellMask ParticleCutCells = particles[p0].CutCells_P(m_LevelSetTracker);
                        CellMask ParticleBoundaryCells = gridData.GetBoundaryCells().Intersect(ParticleCutCells);
                        GetWall(gridData, ParticleBoundaryCells, out double[][] wallPoints);
                        for (int w = 0; w < wallPoints.GetLength(0); w++) {
                            particles[p0].ClosestPointOnOtherObjectToThis = new Vector(particles[p0].Motion.GetPosition(0));
                            if (wallPoints[w] == null)
                                continue;
                            else if (wallPoints[w][0] != 0) {
                                particles[p0].ClosestPointOnOtherObjectToThis = new Vector(wallPoints[w][0], particles[p0].Motion.GetPosition(0)[1]);
                            }
                            else if (wallPoints[w][1] != 0) {
                                particles[p0].ClosestPointOnOtherObjectToThis = new Vector(particles[p0].Motion.GetPosition(0)[0], wallPoints[w][1]);
                            }
                            else
                                continue;
                            CalculateMinimumDistance(particles[p0], out double temp_Distance, out Vector temp_DistanceVector, out Vector temp_ClosestPoint_p0, out bool temp_Overlapping);
                            Distance[p0][ParticleOffset + w] = temp_Distance;
                            Vector normalVector = new Vector(temp_DistanceVector) / temp_DistanceVector.Abs();
                            double temp_SaveTimeStep = DynamicTimestep(particles[p0], temp_ClosestPoint_p0, normalVector, Distance[p0][ParticleOffset + w]);
                            SaveTimeStepArray[p0][ParticleOffset + w] = temp_SaveTimeStep;
                            DistanceVector[p0][ParticleOffset + w] = new Vector(temp_DistanceVector);
                            ClosestPoints[p0][ParticleOffset + w] = new Vector(temp_ClosestPoint_p0);
                            if (temp_SaveTimeStep < SaveTimeStep && temp_SaveTimeStep > 0) {
                                SaveTimeStep = temp_SaveTimeStep;
                                minimalDistance = Distance[p0][ParticleOffset + w];
                            }
                            if (temp_Overlapping) {
                                SaveTimeStep = -m_dt * 0.25; // reset time to find a particle state before they overlap.
                                minimalDistance = double.MaxValue;
                            }
                        }

                        // Step 2.1.3
                        // Test for particle-particle collisions for all particles 
                        // of the current color.
                        // -------------------------------------------------------
                        for (int p1 = p0 + 1; p1 < particles.Length; p1++) {
                            Particle[] currentParticles = new Particle[] { particles[p0], particles[p1] };
                            Console.WriteLine("Particle " + p0 + " Particle " + p1);
                            CalculateMinimumDistance(currentParticles, out double temp_Distance,
                                                     out Vector temp_DistanceVector,
                                                     out Vector[] temp_ClosestPoints,
                                                     out bool temp_Overlapping);
                            Distance[p0][p1] = temp_Distance;
                            Vector normalVector = new Vector(temp_DistanceVector) / temp_DistanceVector.Abs();
                            double temp_SaveTimeStep = DynamicTimestep(currentParticles, temp_ClosestPoints, normalVector, Distance[p0][p1]);
                            SaveTimeStepArray[p0][p1] = temp_SaveTimeStep;
                            DistanceVector[p0][p1] = temp_DistanceVector;
                            ClosestPoints[p0][p1] = temp_ClosestPoints[0];
                            ClosestPoints[p1][p0] = temp_ClosestPoints[1];
                            if (temp_SaveTimeStep < SaveTimeStep && temp_SaveTimeStep > 0) {
                                SaveTimeStep = temp_SaveTimeStep;
                                minimalDistance = Distance[p0][p1];
                            }
                            if (temp_Overlapping) {
                                SaveTimeStep = -m_dt * 0.25; // reset time to find a particle state before they overlap.
                                minimalDistance = double.MaxValue;
                            }
                        }
                    }

                    // Step 2.1.2
                    // Accumulate the current save timestep.
                    // -------------------------------------------------------
                    if (SaveTimeStep >= 0)
                        AccDynamicTimestep += SaveTimeStep;
                    if (AccDynamicTimestep >= m_dt) {
                        break;
                    }
                }
                if (AccDynamicTimestep >= m_dt) {
                    break;
                }
                // Step 3
                // Main collision routine
                // =======================================================
                for (int p0 = 0; p0 < particles.Length; p0++) {

                    // Step 3.1
                    // Particle-wall collisions
                    // -------------------------------------------------------
                    for (int w = 0; w < 4; w++) {
                        if (Distance[p0][ParticleOffset + w] < distanceThreshold && SaveTimeStepArray[p0][ParticleOffset + w] > 0) {
                            Vector CurrentDistanceVector = new Vector(DistanceVector[p0][ParticleOffset + w]);
                            particles[p0].ClosestPointToOtherObject = new Vector(ClosestPoints[p0][ParticleOffset + w]);
                            particles[p0].Motion.SetCollisionTimestep(AccDynamicTimestep);
                            Vector normalVector = new Vector(CurrentDistanceVector) / CurrentDistanceVector.Abs();
                            Vector tangentialVector = new Vector(-normalVector[1], normalVector[0]);
                            particles[p0].Motion.SetCollisionVectors(normalVector, tangentialVector);
                            particles[p0].IsCollided = true;
                            ComputeMomentumBalanceCollision(particles[p0]);
                        }
                    }

                    // Step 3.2
                    // Particle-particle collisions
                    // -------------------------------------------------------
                    for (int p1 = p0 + 1; p1 < particles.Length; p1++) {
                        if (Distance[p0][p1] < distanceThreshold && SaveTimeStepArray[p0][p1] > 0) {
                            List<Particle> collidedParticles = new List<Particle> { particles[p0], particles[p1] };
                            Vector currentDistanceVector = new Vector(DistanceVector[p0][p1]);
                            Vector normalVector = new Vector(currentDistanceVector) / currentDistanceVector.Abs();
                            Vector tangentialVector = new Vector(-normalVector[1], normalVector[0]);
                            particles[p0].ClosestPointToOtherObject = new Vector(ClosestPoints[p0][p1]);
                            particles[p1].ClosestPointToOtherObject = new Vector(ClosestPoints[p1][p0]);
                            foreach (Particle p in collidedParticles) {
                                p.Motion.SetCollisionTimestep(AccDynamicTimestep);
                                p.Motion.SetCollisionVectors(normalVector, tangentialVector);
                                p.IsCollided = true;
                            }
                            ComputeMomentumBalanceCollision(collidedParticles);
                        }
                    }
                }

                // Step 4
                // Add multiple binary collisions together to enable 
                // multiple particle collisions
                // =======================================================
                for (int p = 0; p < particles.Length; p++) {
                    particles[p].Motion.PostProcessCollisionTranslation();
                    particles[p].Motion.PostProcessCollisionRotation();
                }
            }
        }

        /// <summary>
        /// Calculates the dynamic save timestep for a particle-particle interaction.
        /// </summary>
        /// <param name="particle0"></param>
        /// <param name="particle1"></param>
        /// <param name="closestPoint0"></param>
        /// <param name="closestPoint1"></param>
        ///  <param name="normalVector"></param>
        /// <param name="distance"></param>
        private double DynamicTimestep(Particle[] particles, Vector[] closestPoints, Vector normalVector, double distance) {
            double detectCollisionVn_P0;
            double detectCollisionVn_P1;
            if (particles[0].Motion.IncludeTranslation || particles[0].Motion.IncludeRotation) {
                CalculatePointVelocity(particles[0], closestPoints[0], out Vector pointVelocity0);
                detectCollisionVn_P0 = normalVector * pointVelocity0;
            }
            else
                detectCollisionVn_P0 = 0;
            if (particles[1].Motion.IncludeTranslation || particles[1].Motion.IncludeRotation) {
                CalculatePointVelocity(particles[1], closestPoints[1], out Vector pointVelocity1);
                detectCollisionVn_P1 = normalVector * pointVelocity1;
            }
            else
                detectCollisionVn_P1 = 0;
            return (detectCollisionVn_P1 - detectCollisionVn_P0 == 0) ? double.MaxValue : 0.9 * distance / (detectCollisionVn_P1 - detectCollisionVn_P0);
        }

        /// <summary>
        /// Calculates the dynamic save timestep for a particle-wall interaction.
        /// </summary>
        /// <param name="particle"></param>
        /// <param name="closestPoint"></param>
        ///  <param name="normalVector"></param>
        /// <param name="distance"></param>
        private double DynamicTimestep(Particle particle, Vector closestPoint, Vector normalVector, double distance) {
            CalculatePointVelocity(particle, closestPoint, out Vector pointVelocity0);
            double detectCollisionVn_P0 = normalVector * pointVelocity0;
            return detectCollisionVn_P0 == 0 ? double.MaxValue : 0.9 * distance / (-detectCollisionVn_P0);
        }

        /// <summary>
        /// Calculates the velocity of a single point on the surface of a particle.
        /// </summary>
        /// <param name="particle"></param>
        /// <param name="closestPoint"></param>
        ///  <param name="pointVelocity"></param>
        private void CalculatePointVelocity(Particle particle, Vector closestPoint, out Vector pointVelocity) {
            pointVelocity = new Vector(2);
            particle.CalculateRadialVector(closestPoint, out Vector radialVector, out double radialLength);
            pointVelocity[0] = particle.Motion.GetTranslationalVelocity(0)[0] - particle.Motion.GetRotationalVelocity(0) * radialLength * radialVector[1];
            pointVelocity[1] = particle.Motion.GetTranslationalVelocity(0)[1] + particle.Motion.GetRotationalVelocity(0) * radialLength * radialVector[0];
        }

        /// <summary>
        /// Updates the state of the current particles with the dynamic timestep
        /// </summary>
        /// <param name="particles"></param>
        ///  <param name="dynamicTimestep"></param>
        private void UpdateParticleState(Particle[] particles, double dynamicTimestep) {
            for (int p = 0; p < particles.Length; p++) {
                Particle currentParticle = particles[p];
                if (dynamicTimestep != 0) {
                    currentParticle.Motion.CollisionParticlePositionAndAngle(dynamicTimestep);
                }
            }
        }

        /// <summary>
        /// Computes the minimal distance between two particles.
        /// </summary>
        /// <param name="Particle0">
        /// The first particle.
        /// </param>
        ///  <param name="Particle1">
        /// The second particle.
        /// </param>
        /// <param name="Distance">
        /// The minimal distance between the two objects.
        /// </param>
        /// <param name="DistanceVector">
        /// The vector of the minimal distance between the two objects.
        /// </param>
        /// <param name="ClosestPoint_P0">
        /// The point on the first object closest to the second one.
        /// </param>
        /// <param name="ClosestPoint_P1">
        /// The point on the second object closest to the first one.
        /// </param>
        /// <param name="Overlapping">
        /// Is true if the two particles are overlapping.
        /// </param>
        internal void CalculateMinimumDistance(Particle[] Particles, out double Distance, out Vector DistanceVector, out Vector[] ClosestPoints,out bool Overlapping) {
            int spatialDim = Particles[0].Motion.GetPosition(0).Dim;
            Distance = double.MaxValue;
            DistanceVector = new Vector(spatialDim);
            ClosestPoints = new Vector[2];
            Overlapping = false;
            int NoOfSubParticles1 = Particles[1] == null ? 1 : Particles[1].NoOfSubParticles;
;
            for (int i = 0; i < Particles[0].NoOfSubParticles; i++) {
                for (int j = 0; j < NoOfSubParticles1; j++) {
                    GJK_DistanceAlgorithm(Particles[0], i, Particles[1], j, out Vector temp_DistanceVector, out Vector[] temp_ClosestPoints, out Overlapping);
                    if (Overlapping)
                        break;
                    if (temp_DistanceVector.Abs() < Distance) {
                        Distance = temp_DistanceVector.Abs();
                        DistanceVector = new Vector(temp_DistanceVector);
                        ClosestPoints = temp_ClosestPoints.CloneAs();
                    }
                }
            }
        }

        /// <summary>
        /// Computes the minimal distance between a particle and the wall.
        /// </summary>
        /// <param name="Particle0">
        /// The first particle.
        /// </param>
        /// <param name="Distance">
        /// The minimal distance between the two objects.
        /// </param>
        /// <param name="DistanceVector">
        /// The vector of the minimal distance between the two objects.
        /// </param>
        /// <param name="ClosestPoint_P0">
        /// The point on the first object closest to the second one.
        /// </param>
        /// <param name="Overlapping">
        /// Is true if the two particles are overlapping.
        /// </param>
        internal void CalculateMinimumDistance(Particle particle, out double Distance, out Vector DistanceVector, out Vector ClosestPoint, out bool Overlapping) {
            int spatialDim = particle.Motion.GetPosition(0).Dim;
            Distance = double.MaxValue;
            DistanceVector = new Vector(spatialDim);
            ClosestPoint = new Vector(spatialDim);
            Overlapping = false;

            for (int i = 0; i < particle.NoOfSubParticles; i++) {
                GJK_DistanceAlgorithm(particle, i, null, 1, out Vector temp_DistanceVector, out Vector[] temp_ClosestPoints, out Overlapping);
                if (Overlapping)
                    break;
                if (temp_DistanceVector.Abs() < Distance) {
                    Distance = temp_DistanceVector.Abs();
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
        internal void GJK_DistanceAlgorithm(Particle Particle0, int SubParticleID0, Particle Particle1, int SubParticleID1, out Vector DistanceVec, out Vector[] closestPoints, out bool Overlapping) {

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
            Overlapping = false;
            int maxNoOfIterations = 10000;

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

        /// <summary>
        /// Calculates the support point on a single particle.
        /// </summary>
        /// <param name="particle">
        /// Current particle.
        /// </param>
        /// <param name="Vector">
        /// The vector in which direction the support point is searched.
        /// </param>
        /// <param name="supportPoint">
        /// The support point (Cpt. Obvious)
        /// </param>
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
                double testSign = particleDirection[0] * supportVector[1] - particleDirection[1] * supportVector[0];
                double searchStartAngle = (1 - Math.Sign(testSign)) * Math.PI / 2 + Math.Acos((supportVector * particleDirection) / supportVector.L2Norm());
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

        /// <summary>
        /// The core of the GJK-algorithm. Calculates the minimum distance between the current 
        /// simplex and the origin.
        /// </summary>
        /// <param name="simplex">
        /// A list of all support points constituting the simplex.
        /// </param>
        /// <param name="v">
        /// The distance vector.
        /// </param>
        /// <param name="overlapping">
        /// Is true if the simplex contains the origin
        /// </param>
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
                            // Readd the remaining points
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

        /// <summary>
        /// Computes the post-collision velocities of two particles.
        /// </summary>
        /// <param name="collidedParticles">
        /// List of the two colliding particles
        /// </param>
        internal void ComputeMomentumBalanceCollision(List<Particle> collidedParticles) {
            for (int p = 0; p < collidedParticles.Count(); p++) {
                collidedParticles[p].Motion.CalculateNormalAndTangentialVelocity();
                collidedParticles[p].CalculateEccentricity();
            }
            double collisionCoefficient = CalculateCollisionCoefficient(collidedParticles);
            for (int p = 0; p < collidedParticles.Count(); p++) {
                double tempCollisionVn = collidedParticles[p].Motion.IncludeTranslation ? collidedParticles[p].Motion.GetPreCollisionVelocity()[0] + Math.Pow(-1, p + 1) * collisionCoefficient / collidedParticles[p].Motion.Mass_P : 0;
                double tempCollisionVt = collidedParticles[p].Motion.IncludeTranslation ? collidedParticles[p].Motion.GetPreCollisionVelocity()[1] * m_CoefficientOfRestitution : 0;
                double tempCollisionRot = collidedParticles[p].Motion.IncludeRotation ? collidedParticles[p].Motion.GetRotationalVelocity(0) + Math.Pow(-1, p) * collidedParticles[p].Eccentricity * collisionCoefficient / collidedParticles[p].MomentOfInertia : 0;
                collidedParticles[p].Motion.SetCollisionVelocities(tempCollisionVn, tempCollisionVt, tempCollisionRot);
            }
        }

        /// <summary>
        /// Computes the post-collision velocities of one particle after the collision with the wall.
        /// </summary>
        /// <param name="particle"></param>
        internal void ComputeMomentumBalanceCollision(Particle particle) {
            particle.Motion.CalculateNormalAndTangentialVelocity();
            particle.CalculateEccentricity();
            double collisionCoefficient = CalculateCollisionCoefficient(particle);
            double tempCollisionVn = particle.Motion.IncludeTranslation ? particle.Motion.GetPreCollisionVelocity()[0] - collisionCoefficient / particle.Motion.Mass_P : 0;
            double tempCollisionVt = particle.Motion.IncludeTranslation ? particle.Motion.GetPreCollisionVelocity()[1] : 0;
            double tempCollisionRot = particle.Motion.IncludeRotation ? particle.Motion.GetRotationalVelocity(0) + particle.Eccentricity * collisionCoefficient / particle.MomentOfInertia : 0;
            particle.Motion.SetCollisionVelocities(tempCollisionVn, tempCollisionVt, tempCollisionRot);
        }

        /// <summary>
        /// Computes the collision coefficient of two particle after they collided.
        /// </summary>
        /// <param name="collidedParticles"></param>
        private double CalculateCollisionCoefficient(List<Particle> collidedParticles) {
            double[] massReciprocal = new double[2];
            double[] momentOfInertiaReciprocal = new double[2];
            for (int p = 0; p < collidedParticles.Count(); p++) {
                massReciprocal[p] = collidedParticles[p].Motion.IncludeTranslation ? 1 / collidedParticles[p].Motion.Mass_P : 0;
                momentOfInertiaReciprocal[p] = collidedParticles[p].Motion.IncludeRotation ? collidedParticles[p].Eccentricity.Pow2() / collidedParticles[p].MomentOfInertia : 0;
            }
            double collisionCoefficient = (1 + m_CoefficientOfRestitution) * ((collidedParticles[0].Motion.GetPreCollisionVelocity()[0] - collidedParticles[1].Motion.GetPreCollisionVelocity()[0]) / (massReciprocal[0] + massReciprocal[1] + momentOfInertiaReciprocal[0] + momentOfInertiaReciprocal[1]));
            collisionCoefficient += (1 + m_CoefficientOfRestitution) * ((-collidedParticles[0].Eccentricity * collidedParticles[0].Motion.GetRotationalVelocity(0) + collidedParticles[1].Eccentricity * collidedParticles[1].Motion.GetRotationalVelocity(0)) / (massReciprocal[0] + massReciprocal[1] + momentOfInertiaReciprocal[0] + momentOfInertiaReciprocal[1]));
            return collisionCoefficient;
        }

        /// <summary>
        /// Computes the collision coefficient of a particle after the collision with the wall.
        /// </summary>
        /// <param name="particle"></param>
        private double CalculateCollisionCoefficient(Particle particle) {
            double collisionCoefficient = (1 + m_CoefficientOfRestitution) * (particle.Motion.GetPreCollisionVelocity()[0] / (1 / particle.Motion.Mass_P + particle.Eccentricity.Pow2() / particle.MomentOfInertia));
            collisionCoefficient += -(1 + m_CoefficientOfRestitution) * particle.Eccentricity * particle.Motion.GetRotationalVelocity(0) / (1 / particle.Motion.Mass_P + particle.Eccentricity.Pow2() / particle.MomentOfInertia);
            return collisionCoefficient;
        }

        /// <summary>
        /// Searches for all possible walls the particle might collide with
        /// </summary>
        /// <param name="GridData"></param>
        /// <param name="ParticleBoundaryCells">
        /// Cells which have the particle color and are boundary cells
        /// </param>
        /// <param name="WallPoints">
        /// The coordinate of one point on each wall found.
        /// </param>
        internal void GetWall(IGridData GridData, CellMask ParticleBoundaryCells, out double[][] WallPoints) {
            int SpatialDim = ParticleBoundaryCells.GridData.SpatialDimension;
            int NoOfMaxWallEdges = 4;
            WallPoints = new double[NoOfMaxWallEdges][];
            int[][] Cells2Edges = GridData.iLogicalCells.Cells2Edges;
            IList<AffineTrafo> trafo = GridData.iGeomEdges.Edge2CellTrafos;
            foreach (Chunk cnk in ParticleBoundaryCells) {
                for (int i = cnk.i0; i < cnk.JE; i++) {
                    foreach (int e in Cells2Edges[i]) {
                        int eId = (e < 0) ? -e - 1 : e - 1;
                        byte et = GridData.iGeomEdges.EdgeTags[eId];
                        if (GridData.EdgeTagNames[et].Contains("wall") || GridData.EdgeTagNames[et].Contains("Wall")) {
                            int jCell = GridData.iGeomEdges.CellIndices[eId, 0];
                            int iKref = GridData.iGeomEdges.GetRefElementIndex(jCell);

                            NodeSet[] refNodes = GridData.iGeomEdges.EdgeRefElements.Select(Kref2 => Kref2.GetQuadratureRule(5 * 2).Nodes).ToArray();
                            NodeSet Nodes = refNodes.ElementAt(iKref);

                            int trafoIdx = GridData.iGeomEdges.Edge2CellTrafoIndex[eId, 0];
                            MultidimensionalArray transFormed = trafo[trafoIdx].Transform(Nodes);
                            MultidimensionalArray WallVerticies = transFormed.CloneAs();
                            GridData.TransformLocal2Global(transFormed, WallVerticies, jCell);
                            double[] WallPoint1 = WallVerticies.GetRow(0);
                            double[] WallPoint2 = WallVerticies.GetRow(1);
                            if (Math.Abs(WallPoint1[0] - WallPoint2[0]) < 1e-12) {
                                if (WallPoints[0] == null || Math.Abs(WallPoint1[0] - WallPoints[0][0]) < 1e-12)
                                    WallPoints[0] = new double[] { WallPoint1[0], 0 };
                                else if (WallPoints[1] == null || Math.Abs(WallPoint1[0] - WallPoints[1][0]) < 1e-12)
                                    WallPoints[1] = new double[] { WallPoint1[0], 0 };
                                else
                                    throw new ArithmeticException("Error trying to get wall position. Please use horizontal/vertical boudaries");
                            }
                            if (Math.Abs(WallPoint1[1] - WallPoint2[1]) < 1e-12) {
                                if (WallPoints[2] == null || Math.Abs(WallPoint1[1] - WallPoints[2][1]) < 1e-12)
                                    WallPoints[2] = new double[] { 0, WallPoint1[1] };
                                else if (WallPoints[3] == null || Math.Abs(WallPoint1[1] - WallPoints[3][1]) < 1e-12)
                                    WallPoints[3] = new double[] { 0, WallPoint1[1] };
                                else
                                    throw new ArithmeticException("Error trying to get wall position. Please use horizontal/vertical boudaries");
                            }
                        }
                    }
                }
            }
        }
    }
}

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
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSERO_Solver {
    class DistanceAlgorithm {
        /// <summary>
        /// Constructor for  the distance algorithm. Provides the GJK-algorithm to calculate the minimal distance between two particles or a particle and a wall.
        /// </summary>
        /// <param name="Particles">An array containing one or two particles.</param>
        /// <param name="Tolerance">A tolerance parameter</param>
        public DistanceAlgorithm(Particle[] Particles, double Tolerance) {
            this.Particles = Particles;
            this.Tolerance = Tolerance;
            int spatialDimension = this.Particles[0].Motion.GetPosition(0).Dim;
            DistanceVector = new Vector(spatialDimension);
            ClosestPoints = new Vector[2];
            ClosestPoints[0] = new Vector(spatialDimension);
            ClosestPoints[1] = new Vector(spatialDimension);
            Overlapping = false;
        }

        /// <summary>
        /// The minimal distance vector between two particles or a particle and a wall.
        /// </summary>
        [DataMember]
        public Vector DistanceVector { get; private set; }

        /// <summary>
        /// The closest point on one particle towards the other one.
        /// </summary>
        [DataMember]
        public Vector[] ClosestPoints { get; private set; }

        /// <summary>
        /// Flag determining overlapping particles
        /// </summary>
        [DataMember]
        public bool Overlapping { get; private set; }

        [DataMember]
        private readonly double Tolerance;
        [DataMember]
        private readonly Particle[] Particles;


        /// <summary>
        /// Calculates the minimum distance between two particles
        /// </summary>
        public void CalculateTwoParticleDistance() {
            using (new FuncTrace()) {
                if (Particles.Length != 2)
                    throw new ArgumentOutOfRangeException("Number of particles in the distance calculation needs to be equal to two!");
                if (Particles[0] == null || Particles[1] == null)
                    throw new ArgumentNullException("Empty particle definition");
                for (int i = 0; i < Particles[0].NoOfSubParticles; i++) {
                    for (int j = 0; j < Particles[1].NoOfSubParticles; j++) {
                        GJK_DistanceAlgorithm(Particles[0], i, Particles[1], j);
                        if (Overlapping) {
                            DistanceVector = new Vector(Particles[0].Motion.GetPosition(0) - Particles[1].Motion.GetPosition(1));
                            return;
                        }
                    }
                }
            }
        }


        /// <summary>
        /// Computes the minimal distance between a particle and the wall.
        /// </summary>
        /// /// <param name="WallPoint">
        /// A point on the wall.
        /// </param>
        public void CalculateParticleWallDistance(Vector WallPoint) {
            using (new FuncTrace()) {
                if (Particles.Length != 1)
                    throw new ArgumentOutOfRangeException("Number of particles in the distance calculation towards a wall needs to be equal to one!");
                Particle particle = Particles[0];
                if (particle == null)
                    throw new ArgumentNullException("Empty particle definition");

                particle.ClosestPointOnOtherObjectToThis = new Vector(WallPoint);
                for (int i = 0; i < particle.NoOfSubParticles; i++) {
                    GJK_DistanceAlgorithm(particle, i, null, 1);
                    if (Overlapping) {
                        DistanceVector = new Vector(Particles[0].Motion.GetPosition(0) - WallPoint);
                        return;
                    }
                }
            }
        }

        /// <summary>
        /// Computes the distance between two objects (particles or walls). Algorithm based on
        /// E.G.Gilbert, D.W.Johnson, S.S.Keerthi.
        /// </summary>
        /// <param name="Particle">
        /// The first particle.
        /// </param>
        /// <param name="SubParticleID0">
        /// In case of concave particles the particle is divided into multiple convex subparticles. Each of them has its one ID and needs to be tested as if it was a complete particle.
        /// </param>
        ///  <param name="SecondObject">
        /// The second particle, if Particle1 == null it is assumed to be a wall.
        /// </param>
        /// <param name="SubParticleID1">
        /// In case of concave particles the particle is divided into multiple convex subparticles. Each of them has its one ID and needs to be tested as if it was a complete particle.
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
        private void GJK_DistanceAlgorithm(Particle Particle, int SubParticleID0, Particle SecondObject, int SubParticleID1) {
            int NoOfVirtualDomainsP0 = Particle.Motion.OriginInVirtualPeriodicDomain.Count;
            int NoOfVirtualDomainsP1 = IsParticle(SecondObject) ? SecondObject.Motion.OriginInVirtualPeriodicDomain.Count : 0;
            Debug.Assert(NoOfVirtualDomainsP0 == NoOfVirtualDomainsP1);
            int spatialDim = Particle.Motion.GetPosition(0).Dim;
            Overlapping = false;
            DistanceVector = spatialDim switch {
                2 => new(int.MaxValue, int.MaxValue),
                3 => new(int.MaxValue, int.MaxValue, int.MaxValue),
                _ => throw new NotImplementedException("GJK-Algorithm only for 2D or 3D"),
            };
            ClosestPoints = new Vector[2];
            Vector[] tempClosestPoints = new Vector[2];
            for (int d1 = 0; d1 < NoOfVirtualDomainsP0 + 1; d1++) {
                for (int d2 = 0; d2 < NoOfVirtualDomainsP1 + 1; d2++) {
                    // Step 1
                    // Initialize the algorithm with the particle position
                    // =======================================================
                    Vector[] positionVectors = new Vector[2];
                    positionVectors[0] = d1 == NoOfVirtualDomainsP0
                        ? new Vector(Particle.Motion.GetPosition(0))
                        : new Vector(Particle.Motion.GetPosition(0) + Particle.Motion.OriginInVirtualPeriodicDomain[d1]);
                    if (IsParticle(SecondObject))
                        positionVectors[1] = d2 == NoOfVirtualDomainsP0
                                                ? SecondObject.Motion.GetPosition(0)
                                                : SecondObject.Motion.GetPosition(0) + SecondObject.Motion.OriginInVirtualPeriodicDomain[d2];
                    else positionVectors[1] = Particle.ClosestPointOnOtherObjectToThis;

                    Vector[] orientationAngle = new Vector[2];
                    orientationAngle[0] = new Vector(Particle.Motion.GetAngle(0));
                    if (IsParticle(SecondObject))
                        orientationAngle[1] = new Vector(SecondObject.Motion.GetAngle(0));

                    Vector supportVector = positionVectors[0] - positionVectors[1];
                    if (d1 == d2 && d1 != NoOfVirtualDomainsP0)
                        continue;
                    if (SecondObject == null) {
                        if (supportVector.Abs() > 1.5 * Particle.GetLengthScales().Max() && d1 != NoOfVirtualDomainsP0)
                            continue;
                    } else if (supportVector.Abs() > 1.5 * (Particle.GetLengthScales().Max() + SecondObject.GetLengthScales().Max()) && (d1 != NoOfVirtualDomainsP0 || d2 != NoOfVirtualDomainsP1))
                        continue;

                    if (supportVector.Abs() == 0)
                        throw new ArgumentOutOfRangeException("Support vector cannot have zero length");
                    supportVector.CheckForNanOrInfV();

                    // Define the simplex, which contains all points to be tested for their distance (max. 3 points in 2D)
                    List<Vector> Simplex = new() { new Vector(supportVector) };

                    tempClosestPoints[0] = new Vector(spatialDim);
                    tempClosestPoints[1] = new Vector(spatialDim);
                    int maxNoOfIterations = 1000;

                    // Step 2
                    // Start the iteration
                    // =======================================================
                    for (int i = 0; i <= maxNoOfIterations; i++) {
                        Vector negativeSupportVector = new Vector(spatialDim) - supportVector;//-= not possible with vectors?

                        // Calculate the support point of the two particles, 
                        // which are the closest points if the algorithm is finished.
                        // -------------------------------------------------------
                        tempClosestPoints[0] = Particle.GetSupportPoint(negativeSupportVector, positionVectors[0], orientationAngle[0], SubParticleID0, Tolerance);

                        if (IsParticle(SecondObject)) // Particle-Particle collision
                            tempClosestPoints[1] = SecondObject.GetSupportPoint(supportVector, positionVectors[1], orientationAngle[1], SubParticleID1, Tolerance);
                        else {// Particle-wall collision
                            if (positionVectors[0][0] == positionVectors[1][0])
                                tempClosestPoints[1] = new Vector(tempClosestPoints[0][0], positionVectors[1][1]);
                            else
                                tempClosestPoints[1] = new Vector(positionVectors[1][0], tempClosestPoints[0][1]);
                        }

                        // The current support point can be found by forming 
                        // the difference of the support points on the two particles
                        // -------------------------------------------------------
                        Vector supportPoint = tempClosestPoints[0] - tempClosestPoints[1];
                        supportPoint.CheckForNanOrInfV();
                        if (d1 < NoOfVirtualDomainsP0)
                            tempClosestPoints[0] = new Vector(tempClosestPoints[0] - Particle.Motion.OriginInVirtualPeriodicDomain[d1]);
                        if (d2 < NoOfVirtualDomainsP1)
                            tempClosestPoints[1] = new Vector(tempClosestPoints[1] - SecondObject.Motion.OriginInVirtualPeriodicDomain[d2]);

                        // If the condition is true
                        // we have found the closest points!
                        // -------------------------------------------------------
                        if (((supportVector * negativeSupportVector) - (supportPoint * negativeSupportVector)) >= -1e-12 && i > 1)
                            break;

                        // Add new support point to simplex
                        // -------------------------------------------------------
                        Simplex.Insert(0, new Vector(supportPoint));

                        // Calculation the new support vector with the distance
                        // algorithm
                        // -------------------------------------------------------
                        supportVector = SimplexDistance(Simplex);

                        // End algorithm if the two objects are overlapping.
                        // -------------------------------------------------------
                        if (Overlapping)
                            return;

                        // Could not find the closest points... crash!
                        // -------------------------------------------------------
                        if (i > maxNoOfIterations)
                            throw new Exception("No convergence in GJK-algorithm, reached iteration #" + i + ". Note: GJK can only handle convex particles.");
                    }
                    if (supportVector.Abs() < DistanceVector.Abs()) {
                        DistanceVector = new Vector(supportVector);
                        ClosestPoints[0] = new Vector(tempClosestPoints[0]);
                        ClosestPoints[1] = new Vector(tempClosestPoints[1]);
                    }
                }
            }
        }

        /// <summary>
        /// Calculating the distance between the origin and the <paramref name="simplex"/>. 
        /// See Ericson, Christer. Real-Time Collision Detection. 2nd ed. San Francisco, CA: Morgan Kaufmann Publishers, 2005.
        /// </summary>
        /// <param name="simplex"></param>
        /// <param name="overlapping"></param>
        /// <returns></returns>
        private Vector SimplexDistance(List<Vector> simplex) {
            int spatialDimension = simplex[0].Dim;
            if (spatialDimension != 2)
                throw new NotImplementedException("Distance algorithm currently only for 2D.");
            Vector supportVector = new(spatialDimension);
            Overlapping = false;

            // Step 1
            // Test for multiple Simplex-points 
            // and remove the duplicates
            // =======================================================
            for (int s1 = 0; s1 < simplex.Count; s1++) {
                for (int s2 = s1 + 1; s2 < simplex.Count; s2++) {
                    if ((simplex[s1] - simplex[s2]).Abs() < 1e-12) {
                        simplex.RemoveAt(s2);
                    }
                }
            }

            // Step 2
            // Calculate dot product between all simplex vectors and 
            // save to an 2D-array.
            // =======================================================
            double[][] dotProductSimplex = new double[simplex.Count][];
            for (int s1 = 0; s1 < simplex.Count; s1++) {
                dotProductSimplex[s1] = new double[simplex.Count];
                for (int s2 = s1; s2 < simplex.Count; s2++) {
                    dotProductSimplex[s1][s2] = simplex[s1] * simplex[s2];
                }
            }

            // Step 3
            // Main routine to determine the relative position of
            // the simplex towards the origin.
            // =======================================================
            // The simplex contains only one element, which must be
            // the closest point of this simplex to the origin
            // -------------------------------------------------------
            if (simplex.Count == 1)
                supportVector = new Vector(simplex[0]);

            // The simplex contains two elements, lets test which is
            // closest to the origin
            // -------------------------------------------------------
            else if (simplex.Count == 2) {
                // One of the simplex point is closest to the origin, 
                // choose this and delete the other one.
                // -------------------------------------------------------
                bool continueAlgorithmFlag = true;
                for (int s = 0; s < simplex.Count; s++) {
                    if (dotProductSimplex[s][s] - dotProductSimplex[0][1] <= 0) {
                        supportVector = new Vector(simplex[s]);
                        simplex.RemoveAt(Math.Abs(s - 1));
                        continueAlgorithmFlag = false;
                        break;
                    }
                }
                // A point at the line between the two simplex points is
                // closest to the origin, thus we need to keep both points.
                // -------------------------------------------------------
                if (continueAlgorithmFlag) {
                    Vector simplexDistanceVector = simplex[1] - simplex[0];
                    double lambda = spatialDimension switch {
                        2 => Math.Abs(simplex[1].CrossProduct2D(simplexDistanceVector)) / simplexDistanceVector.AbsSquare(),
                        3 => (simplex[1].CrossProduct(simplexDistanceVector)).Abs() / simplexDistanceVector.AbsSquare(),
                        _ => throw new ArgumentOutOfRangeException("Irregular spatial dimension"),
                    };
                    if (lambda == 0) // if the origin lies on the line between the two simplex points, the two objects are overlapping in one point
                        Overlapping = true;
                    supportVector[0] = -lambda * simplexDistanceVector[1];
                    supportVector[1] = lambda * simplexDistanceVector[0];
                }
            }

            // The simplex contains three elements, lets test which is
            // closest to the origin
            // -------------------------------------------------------
            else if (simplex.Count == 3) {
                bool continueAlgorithmFlag = true;
                // Test whether one of the simplex points is closest to 
                // the origin
                // -------------------------------------------------------
                for (int s1 = 0; s1 < simplex.Count; s1++) {
                    int s2 = s1 == 2 ? 2 : 1;
                    int s3 = s1 == 0 ? 0 : 1;
                    if (dotProductSimplex[s1][s1] - dotProductSimplex[0][s2] <= 0 && dotProductSimplex[s1][s1] - dotProductSimplex[s3][2] <= 0) {
                        supportVector = new Vector(simplex[s1]);
                        // Delete the complete simplex and add back the point closest to the origin
                        simplex.Clear();
                        simplex.Add(new Vector(supportVector));
                        continueAlgorithmFlag = false;
                        break;
                    }
                }
                // None of the simplex points was the closest point, 
                // thus, it has to be any point at the edges
                // -------------------------------------------------------
                if (continueAlgorithmFlag) {
                    for (int s1 = simplex.Count() - 1; s1 >= 0; s1--) {
                        int s2 = s1 == 0 ? 1 : 2;
                        int s3 = s1 == 2 ? 1 : 0;
                        // Calculate a triple crossproduct and dotproduct ("quadrupelproduct") of the form (BC x BA) x BA * BX = (BA(BA*BC)-BC(BA*BA))*BX
                        double quadrupelProduct = new();
                        switch (s1) {
                            case 0:
                            double temp1 = dotProductSimplex[1][2] - dotProductSimplex[0][2] - dotProductSimplex[1][1] + dotProductSimplex[0][1];
                            double temp2 = dotProductSimplex[0][1] - dotProductSimplex[0][0] - dotProductSimplex[1][2] + dotProductSimplex[0][2];
                            double temp3 = dotProductSimplex[1][1] - 2 * dotProductSimplex[0][1] + dotProductSimplex[0][0];
                            quadrupelProduct = dotProductSimplex[0][1] * temp1 + dotProductSimplex[1][1] * temp2 + dotProductSimplex[1][2] * temp3;
                            break;
                            case 1:
                            temp1 = -dotProductSimplex[2][2] + dotProductSimplex[0][2] + dotProductSimplex[1][2] - dotProductSimplex[0][1];
                            temp2 = dotProductSimplex[2][2] - 2 * dotProductSimplex[0][2] + dotProductSimplex[0][0];
                            temp3 = dotProductSimplex[0][2] - dotProductSimplex[0][0] - dotProductSimplex[1][2] + dotProductSimplex[0][1];
                            quadrupelProduct = dotProductSimplex[0][2] * temp1 + dotProductSimplex[1][2] * temp2 + dotProductSimplex[2][2] * temp3;
                            break;
                            case 2:
                            temp1 = dotProductSimplex[2][2] - 2 * dotProductSimplex[1][2] + dotProductSimplex[1][1];
                            temp2 = -dotProductSimplex[2][2] + dotProductSimplex[1][2] + dotProductSimplex[0][2] - dotProductSimplex[0][1];
                            temp3 = dotProductSimplex[1][2] - dotProductSimplex[1][1] - dotProductSimplex[0][2] + dotProductSimplex[0][1];
                            quadrupelProduct = dotProductSimplex[0][2] * temp1 + dotProductSimplex[1][2] * temp2 + dotProductSimplex[2][2] * temp3;
                            break;
                        }
                        // A point on one of the edges is closest to the origin.
                        if (dotProductSimplex[s3][s3] - dotProductSimplex[s3][s2] >= 0 && dotProductSimplex[s2][s2] - dotProductSimplex[s3][s2] >= 0 && quadrupelProduct >= 0 && continueAlgorithmFlag) {
                            Vector simplexDistanceVector = simplex[s2] - simplex[s3];
                            double lambda = spatialDimension switch {
                                2 => Math.Abs(simplex[s2].CrossProduct2D(simplexDistanceVector)) / simplexDistanceVector.AbsSquare(),
                                3 => (simplex[s2].CrossProduct(simplexDistanceVector)).Abs() / simplexDistanceVector.AbsSquare(),
                                _ => throw new ArgumentOutOfRangeException("Irregular spatial dimension"),
                            };
                            supportVector[0] = -lambda * simplexDistanceVector[1];
                            supportVector[1] = lambda * simplexDistanceVector[0];
                            // save the two remaining simplex points and clear the simplex.
                            Vector tempSimplex1 = new(simplex[s2]);
                            Vector tempSimplex2 = new(simplex[s3]);
                            simplex.Clear();
                            // Re-add the remaining points
                            simplex.Add(tempSimplex1);
                            simplex.Add(tempSimplex2);
                            continueAlgorithmFlag = false;
                            break;
                        }
                    }
                }
                // None of the conditions above are true, 
                // thus, the simplex must contain the origin and 
                // the two particles do overlap.
                // -------------------------------------------------------
                if (continueAlgorithmFlag)
                    Overlapping = true;
            } else { //In case of 3D 3rd order simplices might arise, necessary to add them here!
                throw new NotImplementedException("Distance algorithm is only implemented for max 2nd order simplex (a triangle)");
            }
            return supportVector;
        }

        /// <summary>
        /// Returns true if <paramref name="Object"/> is a particle. Returns <see langword="false"/> if <paramref name="Object"/> is a wall.
        /// </summary>
        /// <param name="Object"></param>
        /// <returns></returns>
        private static bool IsParticle(Particle Object) {
            return Object != null;
        }
    }
}

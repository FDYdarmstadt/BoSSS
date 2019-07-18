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
using BoSSS.Application.IBM_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FSI_Solver
{
    class FSI_Collision
    {
        private readonly double m_FluidViscosity;
        private readonly double m_FluidDensity;
        private readonly double m_dt;
        private readonly double m_hMin;
        private readonly int m_CurrentColor;
        private double m_CoefficientOfRestitution;


        private double AccDynamicTimestep = 0;

        private MultidimensionalArray SaveTimeStepArray;
        private MultidimensionalArray Distance;
        private MultidimensionalArray DistanceVector;
        private MultidimensionalArray ClosestPoint_P0;
        private MultidimensionalArray ClosestPoint_P1;

        public FSI_Collision(int currentColor, double fluidViscosity, double fluidDensity, double CoefficientOfRestitution, double dt, double hMin)
        {
            m_FluidViscosity = fluidViscosity;
            m_FluidDensity = fluidDensity;
            m_CoefficientOfRestitution = CoefficientOfRestitution;
            m_dt = dt;
            m_hMin = hMin;
            m_CurrentColor = currentColor;
        }

        public FSI_Collision(){ }

        private readonly FSI_Auxillary Aux = new FSI_Auxillary();
        private readonly FSI_LevelSetUpdate LevelSetUpdate = new FSI_LevelSetUpdate();

        private void CreateCollisionArrarys(int noOfParticles, int spatialDim)
        {
            SaveTimeStepArray = MultidimensionalArray.Create(noOfParticles, noOfParticles + 4);
            Distance = MultidimensionalArray.Create(noOfParticles, noOfParticles + 4);
            DistanceVector = MultidimensionalArray.Create(noOfParticles, noOfParticles + 4, spatialDim);
            ClosestPoint_P0 = MultidimensionalArray.Create(noOfParticles, noOfParticles + 4, spatialDim);
            ClosestPoint_P1 = MultidimensionalArray.Create(noOfParticles, noOfParticles + 4, spatialDim);
        }

        /// <summary>
        /// Update collision forces between two arbitrary particles and add them to forces acting on the corresponding particle
        /// </summary>
        /// <param name="particles">
        /// List of all particles
        /// </param>
        public void CalculateCollision(List<Particle> particles, IGridData gridData, int[] cellColor)
        {
            // =======================================================
            // Step 1
            // Some var definintion
            // =======================================================
            int spatialDim = particles[0].position[0].Length;
            int ParticleOffset = particles.Count();
            double MaxDistance = m_dt * 1e-4;
            int J = gridData.iLogicalCells.NoOfLocalUpdatedCells;
            List<int[]> ColoredCellsSorted = LevelSetUpdate.ColoredCellsFindAndSort(cellColor);
            CellMask ParticleCutCells = LevelSetUpdate.CellsOneColor(gridData, ColoredCellsSorted, m_CurrentColor, J, false);

            // =======================================================
            // Step 2
            // Loop over time until the particles hit.
            // =======================================================
            while (AccDynamicTimestep < m_dt)// the collision needs to take place within the current timestep dt.
            {
                CreateCollisionArrarys(particles.Count(), spatialDim);
                double MinDistance = double.MaxValue;
                double SaveTimeStep = 0;// the timestep size without any collision

                // =======================================================
                // Step 2.1
                // Loop over the distance until a predefined criterion is 
                // met.
                // =======================================================
                while (MinDistance > MaxDistance)
                {
                    // =======================================================
                    // Step 2.1.1
                    // Move the particle with the current save timestep.
                    // =======================================================
                    UpdateParticleState(particles, SaveTimeStep, spatialDim);
                    SaveTimeStep = double.MaxValue;
                    for (int p0 = 0; p0 < particles.Count(); p0++)
                    {
                        // =======================================================
                        // Step 2.1.2
                        // Test for wall collisions for all particles 
                        // of the current color.
                        // =======================================================
                        CellMask ParticleBoundaryCells = gridData.GetBoundaryCells().Intersect(ParticleCutCells);
                        GetWall(gridData, ParticleBoundaryCells, out double[][] wallPoints);
                        for (int w = 0; w < wallPoints.GetLength(0); w++)
                        {
                            particles[p0].closestPointOnOtherObjectToThis = particles[p0].position[0].CloneAs();
                            if (wallPoints[w] == null)
                                continue;
                            else if (wallPoints[w][0] != 0)
                            {
                                particles[p0].closestPointOnOtherObjectToThis[0] = wallPoints[w][0];
                            }
                            else if (wallPoints[w][1] != 0)
                            {
                                particles[p0].closestPointOnOtherObjectToThis[1] = wallPoints[w][1];
                            }
                            else
                                continue;
                            CalculateMinimumDistance(particles[p0], 
                                                     out double temp_Distance,
                                                     out MultidimensionalArray temp_DistanceVector,
                                                     out MultidimensionalArray temp_ClosestPoint_p0,
                                                     out bool temp_Overlapping);
                            Distance[p0, ParticleOffset + w] = temp_Distance;
                            CalculateNormalVector(temp_DistanceVector.To1DArray(), out double[] NormalVector);
                            double temp_SaveTimeStep = DynamicTimestep(particles[p0], temp_ClosestPoint_p0.To1DArray(), NormalVector, Distance[p0, ParticleOffset + w]);
                            SaveTimeStepArray[p0, ParticleOffset + w] = temp_SaveTimeStep;
                            DistanceVector.SetSubArray(temp_DistanceVector, new int[] { p0, ParticleOffset + w, -1 });
                            ClosestPoint_P0.SetSubArray(temp_ClosestPoint_p0, new int[] { p0, ParticleOffset + w, -1 });
                            particles[p0].skipForceIntegration = temp_Distance < 0.2 * m_hMin;
                            if (temp_SaveTimeStep < SaveTimeStep && temp_SaveTimeStep > 0)
                            {
                                SaveTimeStep = temp_SaveTimeStep;
                                MinDistance = Distance[p0, ParticleOffset + w];
                            }
                            if (temp_Overlapping)
                            {
                                SaveTimeStep = -m_dt * 0.25; // reset time to find a particle state before they overlap.
                                MinDistance = double.MaxValue;
                            }
                        }

                        // =======================================================
                        // Step 2.1.3
                        // Test for particle-particle collisions for all particles 
                        // of the current color.
                        // =======================================================
                        for (int p1 = p0 + 1; p1 < particles.Count(); p1++)
                        {
                            CalculateMinimumDistance(particles[p0], particles[p1], out double temp_Distance,
                                                     out MultidimensionalArray temp_DistanceVector,
                                                     out MultidimensionalArray temp_ClosestPoint_p0,
                                                     out MultidimensionalArray temp_ClosestPoint_p1,
                                                     out bool temp_Overlapping);
                            Distance[p0, p1] = temp_Distance;
                            CalculateNormalVector(temp_DistanceVector.To1DArray(), out double[] NormalVector);
                            double temp_SaveTimeStep = DynamicTimestep(particles[p0], particles[p1], temp_ClosestPoint_p0.To1DArray(), temp_ClosestPoint_p1.To1DArray(), NormalVector, Distance[p0, p1]);
                            SaveTimeStepArray[p0, p1] = temp_SaveTimeStep;
                            DistanceVector.SetSubArray(temp_DistanceVector, new int[] { p0, p1, -1 });
                            ClosestPoint_P0.SetSubArray(temp_ClosestPoint_p0, new int[] { p0, p1, -1 });
                            ClosestPoint_P1.SetSubArray(temp_ClosestPoint_p1, new int[] { p0, p1, -1 });
                            double testMinDistance = double.MaxValue;
                            if (testMinDistance > temp_Distance)
                                testMinDistance = temp_Distance;
                            if (temp_SaveTimeStep < SaveTimeStep && temp_SaveTimeStep > 0)
                            {
                                SaveTimeStep = temp_SaveTimeStep;
                                MinDistance = Distance[p0, p1];
                            }
                            if (temp_Overlapping)
                            {
                                SaveTimeStep = -m_dt * 0.25; // reset time to find a particle state before they overlap.
                                MinDistance = double.MaxValue;
                            }
                            particles[p0].skipForceIntegration = testMinDistance < 0.2 * m_hMin;
                            particles[p1].skipForceIntegration = testMinDistance < 0.2 * m_hMin;
                        }
                    }

                    // =======================================================
                    // Step 2.1.2
                    // Accumulate the current save timestep.
                    // =======================================================
                    if (SaveTimeStep != -m_dt)
                        AccDynamicTimestep += SaveTimeStep;
                    if (AccDynamicTimestep > m_dt)
                        break;
                }
                if (AccDynamicTimestep > m_dt)
                    break;

                // =======================================================
                // Step 3
                // Main collision routine
                // =======================================================
                for (int p0 = 0; p0 < particles.Count(); p0++)
                {
                    // =======================================================
                    // Step 3.1
                    // Particle-wall collisions
                    // =======================================================
                    for (int w = 0; w < 4; w++)
                    {
                        if (Distance[p0, ParticleOffset + w] < MaxDistance && SaveTimeStepArray[p0, ParticleOffset + w] > 0)
                        {
                            particles[p0].CollisionTimestep = AccDynamicTimestep;
                            particles[p0].Collided = true;
                            double[] CurrentDistanceVector = DistanceVector.ExtractSubArrayShallow(new int[] { p0, ParticleOffset + w, -1 }).To1DArray();
                            particles[p0].closestPointToOtherObject = ClosestPoint_P0.ExtractSubArrayShallow(new int[] { p0, ParticleOffset + w, -1 }).To1DArray();
                            CalculateNormalVector(CurrentDistanceVector, out double[] NormalVector);
                            CalculateTangentialVector(NormalVector, out double[] TangentialVector);
                            particles[p0].CollisionNormalVector.Add(NormalVector);
                            particles[p0].CollisionTangentialVector.Add(TangentialVector);
                            ComputeMomentumBalanceCollision(particles[p0]);
                        }
                    }

                    // =======================================================
                    // Step 3.2
                    // Particle-particle collisions
                    // =======================================================
                    for (int p1 = p0 + 1; p1 < particles.Count(); p1++)
                    {
                        if (Distance[p0, p1] < MaxDistance && SaveTimeStepArray[p0, p1] > 0)
                        {
                            particles[p0].CollisionTimestep = AccDynamicTimestep;
                            particles[p1].CollisionTimestep = AccDynamicTimestep;
                            particles[p0].Collided = true;
                            particles[p1].Collided = true;
                            List<Particle> collidedParticles = new List<Particle> { particles[p0], particles[p1] };
                            double[] CurrentDistanceVector = DistanceVector.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                            particles[p0].closestPointToOtherObject = ClosestPoint_P0.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                            particles[p1].closestPointToOtherObject = ClosestPoint_P1.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                            CalculateNormalVector(CurrentDistanceVector, out double[] NormalVector);
                            CalculateTangentialVector(NormalVector, out double[] TangentialVector);
                            particles[p0].CollisionNormalVector.Add(NormalVector);
                            particles[p1].CollisionNormalVector.Add(NormalVector);
                            particles[p0].CollisionTangentialVector.Add(TangentialVector);
                            particles[p1].CollisionTangentialVector.Add(TangentialVector);
                            ComputeMomentumBalanceCollision(collidedParticles);
                        }
                    }
                }

                // =======================================================
                // Step 4
                // Add multiple binary collisions together to enable 
                // multiple particle collisions
                // =======================================================
                for (int p = 0; p < particles.Count(); p++)
                {
                    PostProcessCollisionTranslation(particles[p]);
                    PostProcessCollisionRotation(particles[p]);
                }

                // =======================================================
                // Step 5
                // Write already used time to particle.cs.
                // =======================================================
                for (int p = 0; p < particles.Count(); p++)
                {
                    particles[p].CollisionTimestep = AccDynamicTimestep;
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
        private double DynamicTimestep(Particle particle0, Particle particle1, double[] closestPoint0, double[] closestPoint1, double[] normalVector, double distance)
        {
            CalculatePointVelocity(particle0, closestPoint0, out double[] pointVelocity0);
            ProjectVelocityOnVector(normalVector, pointVelocity0, out double detectCollisionVn_P0);
            CalculatePointVelocity(particle1, closestPoint1, out double[] pointVelocity1);
            ProjectVelocityOnVector(normalVector, pointVelocity1, out double detectCollisionVn_P1);
            return (detectCollisionVn_P1 - detectCollisionVn_P0 == 0) ? double.MaxValue : 0.9 * distance / (detectCollisionVn_P1 - detectCollisionVn_P0);
        }

        /// <summary>
        /// Calculates the dynamic save timestep for a particle-wall interaction.
        /// </summary>
        /// <param name="particle"></param>
        /// <param name="closestPoint"></param>
        ///  <param name="normalVector"></param>
        /// <param name="distance"></param>
        private double DynamicTimestep(Particle particle, double[] closestPoint, double[] normalVector, double distance)
        {
            CalculatePointVelocity(particle, closestPoint, out double[] pointVelocity0);
            ProjectVelocityOnVector(normalVector, pointVelocity0, out double detectCollisionVn_P0);
            return detectCollisionVn_P0 == 0 ? double.MaxValue : 0.9 * distance / (-detectCollisionVn_P0);
        }

        /// <summary>
        /// Calculates the velocity of a single point on the surface of a particle.
        /// </summary>
        /// <param name="particle"></param>
        /// <param name="closestPoint"></param>
        ///  <param name="pointVelocity"></param>
        private void CalculatePointVelocity(Particle particle, double[] closestPoint, out double[] pointVelocity)
        {
            pointVelocity = new double[closestPoint.Length];
            particle.CalculateRadialVector(closestPoint, out double[] radialVector, out double radialLength);
            pointVelocity[0] = particle.translationalVelocity[0][0] - particle.rotationalVelocity[0] * radialLength * radialVector[1];
            pointVelocity[1] = particle.translationalVelocity[0][1] + particle.rotationalVelocity[0] * radialLength * radialVector[0];
        }

        /// <summary>
        /// Calculates the scalar product between the velocity and a vector.
        /// </summary>
        /// <param name="vector">
        /// </param>
        /// <param name="VelocityVector">
        /// </param>
        /// <param name="VelocityComponent">
        /// </param>
        internal void ProjectVelocityOnVector(double[] vector, double[] VelocityVector, out double VelocityComponent)
        {
            VelocityComponent = VelocityVector[0] * vector[0] + VelocityVector[1] * vector[1];
        }

        /// <summary>
        /// Updates the state of the current particles with the dynamic timestep
        /// </summary>
        /// <param name="particles"></param>
        ///  <param name="dynamicTimestep"></param>
        /// <param name="spatialDim"></param>
        private void UpdateParticleState(List<Particle> particles, double dynamicTimestep, int spatialDim)
        {
            for (int p = 0; p < particles.Count(); p++)
            {
                Particle currentParticle = particles[p];
                if (dynamicTimestep != 0)
                {
                    for (int d = 0; d < spatialDim; d++)
                    {
                        currentParticle.translationalVelocity[0][d] = currentParticle.translationalVelocity[1][d] + (currentParticle.translationalAcceleration[0][d] + 4 * currentParticle.translationalAcceleration[1][d] + currentParticle.translationalAcceleration[2][d]) * dynamicTimestep / 6;
                        currentParticle.position[0][d] = currentParticle.position[0][d] + (currentParticle.translationalVelocity[0][d] + 4 * currentParticle.translationalVelocity[1][d] + currentParticle.translationalVelocity[2][d]) * dynamicTimestep / 6;
                        Aux.TestArithmeticException(currentParticle.position[0], "particle position");
                    }
                    currentParticle.rotationalVelocity[0] = currentParticle.rotationalVelocity[1] + dynamicTimestep * (currentParticle.rotationalAcceleration[0] + 4 * currentParticle.rotationalAcceleration[1] + currentParticle.rotationalAcceleration[2]) / 6;
                    currentParticle.angle[0] = currentParticle.angle[0] + (currentParticle.rotationalVelocity[0] + 4 * currentParticle.rotationalVelocity[1] + currentParticle.rotationalVelocity[2]) * dynamicTimestep / 6;
                    Aux.TestArithmeticException(currentParticle.angle[0], "particle angle");
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
        internal void CalculateMinimumDistance(Particle Particle0, Particle Particle1, out double Distance, out MultidimensionalArray DistanceVector, out MultidimensionalArray ClosestPoint_P0, out MultidimensionalArray ClosestPoint_P1, out bool Overlapping)
        {
            int SpatialDim = Particle0.position[0].Length;
            Distance = double.MaxValue;
            DistanceVector = MultidimensionalArray.Create(SpatialDim);
            ClosestPoint_P0 = MultidimensionalArray.Create(SpatialDim);
            ClosestPoint_P1 = MultidimensionalArray.Create(SpatialDim);
            Overlapping = false;
            int NoOfSubParticles1 = Particle1 == null ? 1 : Particle1.NoOfSubParticles();

            for (int i = 0; i < Particle0.NoOfSubParticles(); i++)
            {
                for (int j = 0; j < NoOfSubParticles1; j++)
                {
                    GJK_DistanceAlgorithm(Particle0, i, Particle1, j, out double temp_Distance, out double[] temp_DistanceVector, out double[] temp_ClosestPoint_P0, out double[] temp_ClosestPoint_P1, out Overlapping);
                    if (Overlapping)
                        break;
                    if (temp_Distance < Distance)
                    {
                        Distance = temp_Distance;
                        for (int d = 0; d < SpatialDim; d++)
                        {
                            DistanceVector[d] = temp_DistanceVector[d];
                            ClosestPoint_P0[d] = temp_ClosestPoint_P0[d];
                            ClosestPoint_P1[d] = temp_ClosestPoint_P1[d];
                        }
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
        internal void CalculateMinimumDistance(Particle Particle0, out double Distance, out MultidimensionalArray DistanceVector, out MultidimensionalArray ClosestPoint_P0, out bool Overlapping)
        {
            int SpatialDim = Particle0.position[0].Length;
            Distance = double.MaxValue;
            DistanceVector = MultidimensionalArray.Create(SpatialDim);
            ClosestPoint_P0 = MultidimensionalArray.Create(SpatialDim);
            Overlapping = false;

            for (int i = 0; i < Particle0.NoOfSubParticles(); i++)
            {
                GJK_DistanceAlgorithm(Particle0, i, null, 1, out double temp_Distance, out double[] temp_DistanceVector, out double[] temp_ClosestPoint_P0, out _, out Overlapping);
                if (Overlapping)
                    break;
                if (temp_Distance < Distance)
                {
                    Distance = temp_Distance;
                    for (int d = 0; d < SpatialDim; d++)
                    {
                        DistanceVector[d] = temp_DistanceVector[d];
                        ClosestPoint_P0[d] = temp_ClosestPoint_P0[d];
                    }
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
        ///  <param name="Particle1">
        /// The second particle, if Particle1 == null it is assumed to be a wall.
        /// </param>
        /// <param name="Min_Distance">
        /// The minimal distance between the two objects.
        /// </param>
        /// <param name="DistanceVec">
        /// The vector of the minimal distance between the two objects.
        /// </param>
        /// <param name="ClosestPoint0">
        /// The point on the first object closest to the second one.
        /// </param>
        /// <param name="ClosestPoint1">
        /// The point on the second object closest to the first one.
        /// </param>
        /// <param name="Overlapping">
        /// Is true if the two particles are overlapping.
        /// </param>
        internal void GJK_DistanceAlgorithm(Particle Particle0, int SubParticleID0, Particle Particle1, int SubParticleID1, out double Min_Distance, out double[] DistanceVec, out double[] ClosestPoint0, out double[] ClosestPoint1, out bool Overlapping)
        {
            double[] Position0 = Particle0.position[0].CloneAs();
            int SpatialDim = Position0.Length;
            double[] Position1 = Particle1 == null ? Particle0.closestPointOnOtherObjectToThis.CloneAs() : Particle1.position[0].CloneAs();
            double[] v = Aux.VectorDiff(Position0, Position1);
            Aux.TestArithmeticException(v, nameof(v));
            List<double[]> Simplex = new List<double[]> { v.CloneAs() };

            ClosestPoint0 = new double[SpatialDim];
            ClosestPoint1 = new double[SpatialDim];
            Overlapping = false;

            for (int i = 0; i < 10000; i++)
            {
                double[] vt = v.CloneAs();
                for (int d = 0; d < SpatialDim; d++)
                {
                    vt[d] = -v[d];
                }

                // =======================================================
                // Step 1
                // Calculate the support point of the minkowski difference
                // and of the two particles (which are the closest points
                // if the algorithm is finished.
                // =======================================================
                CalculateSupportPoint(Particle0, SubParticleID0, vt, out ClosestPoint0);
                Aux.TestArithmeticException(ClosestPoint0, nameof(ClosestPoint0));
                if (Particle1 != null)
                    CalculateSupportPoint(Particle1, SubParticleID1, v, out ClosestPoint1);
                else
                {
                    ClosestPoint1 = ClosestPoint0.CloneAs();    
                    if (Position0[0] == Position1[0])
                        ClosestPoint1[1] = Position1[1];
                    else
                        ClosestPoint1[0] = Position1[0];
                }
                Aux.TestArithmeticException(ClosestPoint1, nameof(ClosestPoint1));

                double[] SupportPoint = Aux.VectorDiff(ClosestPoint0, ClosestPoint1);
                Aux.TestArithmeticException(SupportPoint, nameof(SupportPoint));

                // =======================================================
                // Step 2
                // Check max(x dot vt)
                // =======================================================
                if ((Aux.DotProduct(v,vt) - Aux.DotProduct(SupportPoint,vt)) >= -1e-12 && i != 0)
                    break;

                // =======================================================
                // Step 3
                // Add new support point to simplex
                // =======================================================
                Simplex.Insert(0, SupportPoint.CloneAs());

                // =======================================================
                // Step 4
                // Calculation the new vector v with the distance
                // algorithm
                // =======================================================
                DistanceAlgorithm(Simplex, out v, out Overlapping);

                // End algorithm if the two objects are overlapping.
                if (Overlapping)
                    break;
            }
            // =======================================================
            // Step 5
            // Return min distance and distance vector.
            // =======================================================
            Min_Distance = v.L2Norm();
            DistanceVec = v.CloneAs();
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
        /// <param name="SupportPoint">
        /// The support point (Cpt. Obvious)
        /// </param>
        private void CalculateSupportPoint(Particle particle, int SubParticleID, double[] Vector, out double[] SupportPoint)
        {
            int SpatialDim = particle.position[0].Length;
            SupportPoint = new double[SpatialDim];
            // A direct formulation of the support function for a sphere exists, thus it is possible to map it to an ellipsoid.
            if (particle is Particle_Ellipsoid || particle is Particle_Sphere)
            {
                particle.GetSupportPoint(SpatialDim, Vector, out SupportPoint);
            }
            // Binary search in all other cases.
            else
            {
                MultidimensionalArray SurfacePoints = particle.GetSurfacePoints(m_hMin);
                MultidimensionalArray SurfacePointsSubParticle = SurfacePoints.ExtractSubArrayShallow(new int[]{ SubParticleID, -1, -1});
                int L = 1;
                int R = SurfacePointsSubParticle.GetLength(0) - 2;
                while (L <= R && L > 0 && R < SurfacePointsSubParticle.GetLength(0) - 1)
                {
                    int Index = (L + R) / 2;
                    double[] RightNeighbour = new double[2];
                    double[] LeftNeighbour = new double[2];
                    for (int d = 0; d < 2; d++)
                    {
                        SupportPoint[d] = SurfacePointsSubParticle[Index, d];
                        LeftNeighbour[d] = SurfacePointsSubParticle[Index - 1, d];
                        RightNeighbour[d] = SurfacePointsSubParticle[Index + 1, d];
                    }
                    if (Aux.DotProduct(SupportPoint, Vector) > Aux.DotProduct(RightNeighbour, Vector) && Aux.DotProduct(SupportPoint, Vector) > Aux.DotProduct(LeftNeighbour, Vector))
                        break; // The current temp_supportPoint is the actual support point.
                    else if (Aux.DotProduct(RightNeighbour, Vector) > Aux.DotProduct(LeftNeighbour, Vector))
                        L = Index + 1; // Search on the right side of the current point.
                    else
                        R = Index - 1; // Search on the left side.
                }
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
        private void DistanceAlgorithm(List<double[]> simplex, out double[] v, out bool overlapping)
        {
            int spatialDim = simplex[0].Length;
            v = new double[spatialDim];
            overlapping = false;

            // =======================================================
            // Step 1
            // Test for multiple Simplex-points 
            // and remove the duplicates
            // =======================================================
            for (int s1 = 0; s1 < simplex.Count(); s1++)
            {
                for (int s2 = s1 + 1; s2 < simplex.Count(); s2++)
                {
                    if (Math.Abs(simplex[s1][0] - simplex[s2][0]) < 1e-8 && Math.Abs(simplex[s1][1] - simplex[s2][1]) < 1e-8)
                    {
                        simplex.RemoveAt(s2);
                    }
                }
            }

            // =======================================================
            // Step 2
            // Calculate dot product between all position vectors and 
            // save to an 2D-array.
            // =======================================================
            List<double[]> dotProductSimplex = new List<double[]>();
            for (int s1 = 0; s1 < simplex.Count(); s1++)
            {
                dotProductSimplex.Add(new double[simplex.Count()]);
                for (int s2 = s1; s2 < simplex.Count(); s2++)
                {
                    dotProductSimplex[s1][s2] = Aux.DotProduct(simplex[s1], simplex[s2]);
                }
            }

            // =======================================================
            // Step 3
            // Main routine to determine the relatve position of
            // the simplex towards the origin.
            // =======================================================
            if (simplex.Count() == 1)
            {
                v = simplex[0];// the only possibility
                Aux.TestArithmeticException(v, nameof(v));
            }
            else if (simplex.Count() == 2)
            {
                // The first simplex point is closest to the origin, choose this and delete the other one.
                if (dotProductSimplex[0][0] - dotProductSimplex[0][1] <= 0)
                {
                    v = simplex[0].CloneAs();
                    simplex.RemoveAt(1);
                    Aux.TestArithmeticException(v, nameof(v));
                }
                // The second simplex point is closest to the origin, choose this and delete the other one.
                else if (dotProductSimplex[1][1] - dotProductSimplex[0][1] <= 0)
                {
                    v = simplex[1].CloneAs();
                    simplex.RemoveAt(0);
                    Aux.TestArithmeticException(v, nameof(v));
                }
                // A point at the line between the two simplex points is closest to the origin, thus we need to keep both points.
                else
                {
                    double[] simplexDistanceVector = Aux.VectorDiff(simplex[1], simplex[0]); 
                    double lambda = Math.Abs(-simplex[1][1] * simplexDistanceVector[0] + simplex[1][0] * simplexDistanceVector[1]) / (simplexDistanceVector[0].Pow2() + simplexDistanceVector[1].Pow2());
                    if(lambda == 0) // if the origin lies on the line between the two simplex points, the two particles are overlapping in one point
                        overlapping = true;
                    v[0] = -lambda * simplexDistanceVector[1];
                    v[1] = lambda * simplexDistanceVector[0];
                    Aux.TestArithmeticException(v, nameof(v));
                }
            }
            else if (simplex.Count() == 3)
            {
                bool continueAlgorithm = true;
                // Test whether one of the simplex points is closest to the origin
                for (int s1 = 0; s1 < simplex.Count(); s1++)
                {
                    int s2 = s1 == 2 ? 2 : 1;
                    int s3 = s1 == 0 ? 0 : 1;
                    if (dotProductSimplex[s1][s1] - dotProductSimplex[0][s2] <= 0 && dotProductSimplex[s1][s1] - dotProductSimplex[s3][2] <= 0)
                    {
                        v = simplex[s1].CloneAs();
                        // Delete the complete simplex and add back the point closest to the origin
                        simplex.Clear();
                        simplex.Add(v.CloneAs());
                        continueAlgorithm = false;
                        Aux.TestArithmeticException(v, nameof(v));
                        break;
                    }
                }
                if (continueAlgorithm)
                {
                    for (int s1 = simplex.Count() - 1; s1 >= 0; s1--)
                    {
                        int s2 = s1 == 0 ? 1 : 2;
                        int s3 = s1 == 2 ? 1 : 0;
                        // Calculate a crossproduct of the form (BC x BA) x BA * BX
                        double crossProduct = new double();
                        switch (s1)
                        {
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
                        if (dotProductSimplex[s3][s3] - dotProductSimplex[s3][s2] >= 0 && dotProductSimplex[s2][s2] - dotProductSimplex[s3][s2] >= 0 && crossProduct >= 0 && continueAlgorithm)
                        {
                            double[] simplexDistanceVector = new double[2];
                            for (int d = 0; d < 2; d++)
                            {
                                simplexDistanceVector[d] = simplex[s2][d] - simplex[s3][d];
                            }
                            double Lambda = (simplex[s2][1] * simplexDistanceVector[0] - simplex[s2][0] * simplexDistanceVector[1]) / (simplexDistanceVector[0].Pow2() + simplexDistanceVector[1].Pow2());
                            v[0] = -Lambda * simplexDistanceVector[1];
                            v[1] = Lambda * simplexDistanceVector[0];
                            // save the two remaining simplex points and clear the simplex.
                            double[] tempSimplex1 = simplex[s2].CloneAs();
                            double[] tempSimplex2 = simplex[s3].CloneAs();
                            simplex.Clear();
                            // Readd the remaining points
                            simplex.Add(tempSimplex1.CloneAs());
                            simplex.Add(tempSimplex2.CloneAs());
                            continueAlgorithm = false;
                            Aux.TestArithmeticException(v, nameof(v));
                            break;
                        }
                    }
                }
                // None of the conditions above are true, thus, the simplex must contain the origin and the two particles do overlap.
                if (continueAlgorithm)
                    overlapping = true;
            }
        }

        /// <summary>
        /// Computes the coefficient of restitution. At small Stokes-No. the coefficient of restitution becomes zero.
        /// </summary>
        /// <param name="stokesNumber"></param>
        /// <param name="coefficientOfRestitution"></param>
        private void ModelCoefficientOfRestitution(double stokesNumber, out double coefficientOfRestitution)
        {
            coefficientOfRestitution = m_CoefficientOfRestitution;
            if (stokesNumber < 9 && stokesNumber != 0)// the Stokes number is zero in case of a dry simulation (no fluid)
                coefficientOfRestitution = 0;
        }

        /// <summary>
        /// Computes the post-collision velocities of two particles.
        /// </summary>
        /// <param name="collidedParticles">
        /// List of the two colliding particles
        /// </param>
        internal void ComputeMomentumBalanceCollision(List<Particle> collidedParticles)
        {
            ModelCoefficientOfRestitution(collidedParticles[0].ComputeParticleSt(m_FluidViscosity, m_FluidDensity, Aux.VectorDiff(collidedParticles[0].translationalVelocity[0], collidedParticles[1].translationalVelocity[0])), out m_CoefficientOfRestitution);

            for (int p = 0; p < collidedParticles.Count(); p++)
            {
                collidedParticles[p].CalculateNormalAndTangentialVelocity();
                collidedParticles[p].CalculateEccentricity();
            }

            CalculateCollisionCoefficient(collidedParticles, out double collisionCoefficient);

            for (int p = 0; p < collidedParticles.Count(); p++)
            {
                double tempCollisionVn = collidedParticles[p].IncludeTranslation ? collidedParticles[p].PreCollisionVelocity[0] + Math.Pow(-1, p + 1) * collisionCoefficient / collidedParticles[p].Mass_P : 0;
                double tempCollisionVt = collidedParticles[p].IncludeTranslation ? collidedParticles[p].PreCollisionVelocity[1] * m_CoefficientOfRestitution : 0;
                double tempCollisionRot = collidedParticles[p].IncludeRotation ? collidedParticles[p].rotationalVelocity[0] + Math.Pow(-1, p) * collidedParticles[p].eccentricity * collisionCoefficient / collidedParticles[p].MomentOfInertia_P : 0;
                collidedParticles[p].CollisionTranslationalVelocity.Add(new double[] { tempCollisionVn, tempCollisionVt });
                collidedParticles[p].CollisionRotationalVelocity.Add(tempCollisionRot);
            }
        }

        /// <summary>
        /// Computes the post-collision velocities of one particle after the collision with the wall.
        /// </summary>
        /// <param name="particle"></param>
        internal void ComputeMomentumBalanceCollision(Particle particle)
        {
            ModelCoefficientOfRestitution(particle.ComputeParticleSt(m_FluidViscosity, m_FluidDensity), out m_CoefficientOfRestitution);

            particle.CalculateNormalAndTangentialVelocity();
            particle.CalculateEccentricity();

            CalculateCollisionCoefficient(particle, out double collisionCoefficient);

            double tempCollisionVn = particle.IncludeTranslation ? particle.PreCollisionVelocity[0] - collisionCoefficient / particle.Mass_P : 0;
            double tempCollisionVt = particle.IncludeTranslation ? particle.PreCollisionVelocity[1] : 0;
            double tempCollisionRot = particle.IncludeRotation ? particle.rotationalVelocity[0] + particle.eccentricity * collisionCoefficient / particle.MomentOfInertia_P : 0;
            particle.CollisionTranslationalVelocity.Add(new double[] { tempCollisionVn, tempCollisionVt });
            particle.CollisionRotationalVelocity.Add(tempCollisionRot);
        }

        /// <summary>
        /// Computes the collision coefficient of two particle after they collided.
        /// </summary>
        /// <param name="collidedParticles"></param>
        /// <param name="collisionCoefficient"></param>
        private void CalculateCollisionCoefficient(List<Particle> collidedParticles, out double collisionCoefficient)
        {
            double[] massReciprocal = new double[2];
            double[] momentOfInertiaReciprocal = new double[2];
            for (int p = 0; p < collidedParticles.Count(); p++)
            {
                massReciprocal[p] = collidedParticles[p].IncludeTranslation ? 1 / collidedParticles[p].Mass_P : 0;
                momentOfInertiaReciprocal[p] = collidedParticles[p].IncludeRotation ? collidedParticles[p].eccentricity.Pow2() / collidedParticles[p].MomentOfInertia_P : 0;
            }
            collisionCoefficient = (1 + m_CoefficientOfRestitution) * ((collidedParticles[0].PreCollisionVelocity[0] - collidedParticles[1].PreCollisionVelocity[0]) / (massReciprocal[0] + massReciprocal[1] + momentOfInertiaReciprocal[0] + momentOfInertiaReciprocal[1]));
            collisionCoefficient += (1 + m_CoefficientOfRestitution) * ((-collidedParticles[0].eccentricity * collidedParticles[0].rotationalVelocity[0] + collidedParticles[1].eccentricity * collidedParticles[1].rotationalVelocity[0]) / (massReciprocal[0] + massReciprocal[1] + momentOfInertiaReciprocal[0] + momentOfInertiaReciprocal[1]));
        }

        /// <summary>
        /// Computes the collision coefficient of a particle after the collision with the wall.
        /// </summary>
        /// <param name="particle"></param>
        /// <param name="collisionCoefficient"></param>
        private void CalculateCollisionCoefficient(Particle particle, out double collisionCoefficient)
        {
            collisionCoefficient = (1 + m_CoefficientOfRestitution) * (particle.PreCollisionVelocity[0] / (1 / particle.Mass_P + particle.eccentricity.Pow2() / particle.MomentOfInertia_P));
            collisionCoefficient += -(1 + m_CoefficientOfRestitution) * particle.eccentricity * particle.rotationalVelocity[0] / (1 / particle.Mass_P + particle.eccentricity.Pow2() / particle.MomentOfInertia_P);
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
        internal void GetWall(IGridData GridData, CellMask ParticleBoundaryCells, out double[][] WallPoints)
        {
            int SpatialDim = ParticleBoundaryCells.GridData.SpatialDimension;
            int NoOfMaxWallEdges = 4;
            WallPoints = new double[NoOfMaxWallEdges][];
            int[][] Cells2Edges = GridData.iLogicalCells.Cells2Edges;
            IList<BoSSS.Platform.LinAlg.AffineTrafo> trafo = GridData.iGeomEdges.Edge2CellTrafos;
            foreach (Chunk cnk in ParticleBoundaryCells)
            {
                for (int i = cnk.i0; i < cnk.JE; i++)
                {
                    foreach (int e in Cells2Edges[i])
                    {
                        int eId = (e < 0) ? -e - 1 : e - 1;
                        byte et = GridData.iGeomEdges.EdgeTags[eId];
                        if (GridData.EdgeTagNames[et].Contains("wall") || GridData.EdgeTagNames[et].Contains("Wall"))
                        {
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
                            if (Math.Abs(WallPoint1[0] - WallPoint2[0]) < 1e-12)
                            {
                                if (WallPoints[0] == null || Math.Abs(WallPoint1[0] - WallPoints[0][0]) < 1e-12)
                                    WallPoints[0] = new double[] { WallPoint1[0], 0 };
                                else if (WallPoints[1] == null || Math.Abs(WallPoint1[0] - WallPoints[1][0]) < 1e-12)
                                    WallPoints[1] = new double[] { WallPoint1[0], 0 };
                                else
                                    throw new ArithmeticException("Error trying to get wall position. Please use horizontal/vertical boudaries");
                            }
                            if (Math.Abs(WallPoint1[1] - WallPoint2[1]) < 1e-12)
                            {
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

        /// <summary>
        /// Collision post-processing. Sums up the results of the multiple binary collisions of one timestep
        /// </summary>
        /// <param name="particle">
        /// The particle to be processed
        /// </param>
        private void PostProcessCollisionTranslation(Particle particle)
        {
            int SpatialDim = particle.position[0].Length;
            if (particle.CollisionTranslationalVelocity.Count() >= 1)
            {
                double[] Normal = new double[SpatialDim];
                double[] Tangential = new double[SpatialDim];
                for (int t = 0; t < particle.CollisionTranslationalVelocity.Count(); t++)
                {
                    for (int d = 0; d < SpatialDim; d++)
                    {
                        Normal[d] += particle.CollisionNormalVector[t][d];
                        Tangential[d] += particle.CollisionTangentialVector[t][d];
                    }
                }

                Normal.ScaleV(1 / Math.Sqrt(Normal[0].Pow2() + Normal[1].Pow2()));
                Tangential.ScaleV(1 / Math.Sqrt(Tangential[0].Pow2() + Tangential[1].Pow2()));
                double[] Cos = new double[particle.CollisionTranslationalVelocity.Count()];
                double[] Sin = new double[particle.CollisionTranslationalVelocity.Count()];
                double temp_NormalVel = 0;
                double temp_TangentialVel = 0;
                for (int t = 0; t < particle.CollisionTranslationalVelocity.Count(); t++)
                {
                    for (int d = 0; d < SpatialDim; d++)
                    {
                        Cos[t] += Normal[d] * particle.CollisionNormalVector[t][d];
                    }
                    Sin[t] = Cos[t] == 1 ? 0 : particle.CollisionNormalVector[t][0] > Normal[0] ? Math.Sqrt(1 + 1e-15 - Cos[t].Pow2()) : -Math.Sqrt(1 + 1e-15 - Cos[t].Pow2());
                    temp_NormalVel += particle.CollisionTranslationalVelocity[t][0] * Cos[t] - particle.CollisionTranslationalVelocity[t][1] * Sin[t];
                    temp_TangentialVel += particle.CollisionTranslationalVelocity[t][0] * Sin[t] + particle.CollisionTranslationalVelocity[t][1] * Cos[t];

                }
                temp_NormalVel /= particle.CollisionTranslationalVelocity.Count();
                temp_TangentialVel /= particle.CollisionTranslationalVelocity.Count();

                particle.translationalVelocity.Insert(0, new double[2]);
                for (int d = 0; d < SpatialDim; d++)
                {
                    particle.translationalVelocity[0][d] = Normal[d] * temp_NormalVel + Tangential[d] * temp_TangentialVel;
                }

                particle.CollisionTranslationalVelocity.Clear();
                particle.CollisionNormalVector.Clear();
                particle.CollisionTangentialVector.Clear();
            }
        }

        /// <summary>
        /// Collision post-processing. Sums up the results for the angular velocity of the multiple binary collisions of one timestep
        /// </summary>
        /// <param name="particle">
        /// The particle to be processed
        /// </param>
        private void PostProcessCollisionRotation(Particle particle)
        {
            if (particle.CollisionRotationalVelocity.Count() >= 1)
            {
                particle.rotationalVelocity[0] = particle.CollisionRotationalVelocity.Sum() / particle.CollisionRotationalVelocity.Count();
                particle.CollisionRotationalVelocity.Clear();

                if (double.IsNaN(particle.rotationalVelocity[0]) || double.IsInfinity(particle.rotationalVelocity[0]))
                    throw new ArithmeticException("Error trying to update particle angular velocity during collision post-processing. The angular velocity is:  " + particle.rotationalVelocity[0]);
            }
        }

        /// <summary>
        /// Calculates the normal vector from the min distance vector between the
        /// two particles.
        /// </summary>
        /// <param name="distanceVec">
        /// The min distance vector between the two particles.
        /// </param>
        /// <param name="NormalVector">
        /// The collision normal vector.
        /// </param>
        private void CalculateNormalVector(double[] distanceVec, out double[] NormalVector)
        {
            NormalVector = distanceVec.CloneAs();
            NormalVector.ScaleV(1 / Math.Sqrt(distanceVec[0].Pow2() + distanceVec[1].Pow2()));
        }

        /// <summary>
        /// Calculates the tangential vector from the normal vector between the
        /// two particles.
        /// </summary>
        /// <param name="NormalVector">
        /// The collision normal vector.
        /// </param>
        /// <param name="TangentialVector">
        /// The collision tangential vector.
        /// </param>
        private void CalculateTangentialVector(double[] NormalVector, out double[] TangentialVector)
        {
            TangentialVector = new double[] { -NormalVector[1], NormalVector[0] };
        }
    }
}

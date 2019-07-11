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
        private double m_CoefficientOfRestitution;

        public FSI_Collision(double fluidViscosity, double fluidDensity, double CoefficientOfRestitution, double dt, double hMin)
        {
            m_FluidViscosity = fluidViscosity;
            m_FluidDensity = fluidDensity;
            m_CoefficientOfRestitution = CoefficientOfRestitution;
            m_dt = dt;
            m_hMin = hMin;
        }

        public FSI_Collision(){ }

        private readonly FSI_Auxillary Aux = new FSI_Auxillary();
        private readonly FSI_LevelSetUpdate LevelSetUpdate = new FSI_LevelSetUpdate();


        /// <summary>
        /// Update collision forces between two arbitrary particles and add them to forces acting on the corresponding particle
        /// </summary>
        /// <param name="Particles">
        /// List of all particles
        /// </param>
        /// <param name="dt"></param>
        public void CalculateCollision(List<Particle> Particles, IGridData GridData, int[] CellColor)
        {
            LevelSetUpdate.DetermineGlobalParticleColor(GridData, CellColor, Particles, out int[] GlobalParticleColor);

            for (int i = 0; i < GlobalParticleColor.Length; i++)
            {
                int CurrentColor = GlobalParticleColor[i];
                int[] ParticlesOfCurrentColor = LevelSetUpdate.FindParticlesOneColor(GlobalParticleColor, CurrentColor);

                if (ParticlesOfCurrentColor.Length >= 1 && CurrentColor != 0)
                {
                    int spatialDim = 2;
                    int ParticleOffset = ParticlesOfCurrentColor.Length;
                    double MaxDistance = m_dt * 1e-4;
                    double AccDynamicTimestep = 0;
                    double PostCollisionTimestep = 0;
                    bool Overlapping = false;
                    int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
                    List<int[]> ColoredCellsSorted = LevelSetUpdate.ColoredCellsFindAndSort(CellColor);
                    CellMask ParticleCutCells = LevelSetUpdate.CellsOneColor(GridData, ColoredCellsSorted, CurrentColor, J, false);

                    while (AccDynamicTimestep < m_dt)
                    {
                        MultidimensionalArray SaveTimeStepArray = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length + 4);
                        MultidimensionalArray Distance = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length + 4);
                        MultidimensionalArray DistanceVector = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length + 4, spatialDim);
                        MultidimensionalArray ClosestPoint_P0 = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length + 4, spatialDim);
                        MultidimensionalArray ClosestPoint_P1 = MultidimensionalArray.Create(ParticlesOfCurrentColor.Length, ParticlesOfCurrentColor.Length + 4, spatialDim);
                        double MinDistance = double.MaxValue;
                        double SaveTimeStep = 0;
                        while (MinDistance > MaxDistance)
                        {
                            for (int p = 0; p < ParticlesOfCurrentColor.Length; p++)
                            {
                                UpdateParticleState(Particles[ParticlesOfCurrentColor[p]], SaveTimeStep, spatialDim);
                            }
                            Overlapping = false;
                            double LastSaveTimeStep = SaveTimeStep;
                            SaveTimeStep = double.MaxValue;
                            for (int p0 = 0; p0 < ParticlesOfCurrentColor.Length; p0++)
                            {
                                Particle Particle0 = Particles[ParticlesOfCurrentColor[p0]];
                                CellMask ParticleBoundaryCells = GridData.GetBoundaryCells().Intersect(ParticleCutCells);
                                GetWall(GridData, ParticleBoundaryCells, out double[][] wallPoints);
                                for (int w = 0; w < wallPoints.GetLength(0); w++)
                                {
                                    Particle0.closestPointOnOtherObjectToThis = Particle0.Position[0].CloneAs();
                                    if (wallPoints[w] == null)
                                        continue;
                                    else if (wallPoints[w][0] != 0)
                                    {
                                        Particle0.closestPointOnOtherObjectToThis[0] = wallPoints[w][0];
                                    }
                                    else if (wallPoints[w][1] != 0)
                                    {
                                        Particle0.closestPointOnOtherObjectToThis[1] = wallPoints[w][1];
                                    }
                                    else
                                        continue;
                                    CalculateMinimumDistance(Particle0, null,
                                                           out double temp_Distance,
                                                           out MultidimensionalArray temp_DistanceVector,
                                                           out MultidimensionalArray temp_ClosestPoint_p0,
                                                           out MultidimensionalArray temp_ClosestPoint_p1,
                                                           out bool temp_Overlapping);
                                    Distance[p0, ParticleOffset + w] = temp_Distance;
                                    CalculateNormalVector(temp_DistanceVector.To1DArray(), out double[] NormalVector);
                                    double temp_SaveTimeStep = DynamicTimestep(Particle0, null, temp_ClosestPoint_p0.To1DArray(), temp_ClosestPoint_p1.To1DArray(), NormalVector, Distance[p0, ParticleOffset + w]);
                                    SaveTimeStepArray[p0, ParticleOffset + w] = temp_SaveTimeStep;
                                    DistanceVector.SetSubArray(temp_DistanceVector, new int[] { p0, ParticleOffset + w, -1 });
                                    ClosestPoint_P0.SetSubArray(temp_ClosestPoint_p0, new int[] { p0, ParticleOffset + w, -1 });
                                    ClosestPoint_P1.SetSubArray(temp_ClosestPoint_p1, new int[] { p0, ParticleOffset + w, -1 });
                                    double testMinDistance = double.MaxValue;
                                    if (testMinDistance > temp_Distance)
                                        testMinDistance = temp_Distance;
                                    if (temp_SaveTimeStep < SaveTimeStep && temp_SaveTimeStep > 0)
                                    {
                                        SaveTimeStep = temp_SaveTimeStep;
                                        MinDistance = Distance[p0, ParticleOffset + w];
                                    }
                                    if (temp_Overlapping)
                                    {
                                        if (LastSaveTimeStep > m_dt * 1e-5)
                                            SaveTimeStep = -LastSaveTimeStep / 2;
                                        else
                                            SaveTimeStep = -m_dt * 0.25;
                                        MinDistance = double.MaxValue;
                                        Overlapping = true;
                                    }
                                    Particle0.skipForceIntegration = testMinDistance < 0.2 * m_hMin;
                                }
                                for (int p1 = p0 + 1; p1 < ParticlesOfCurrentColor.Length; p1++)
                                {
                                    Particle Particle1 = Particles[ParticlesOfCurrentColor[p1]];
                                    CalculateMinimumDistance(Particle0, Particle1, out double temp_Distance, out MultidimensionalArray temp_DistanceVector, out MultidimensionalArray temp_ClosestPoint_p0, out MultidimensionalArray temp_ClosestPoint_p1, out bool temp_Overlapping);
                                    Distance[p0, p1] = temp_Distance;
                                    CalculateNormalVector(temp_DistanceVector.To1DArray(), out double[] NormalVector);
                                    double temp_SaveTimeStep = DynamicTimestep(Particle0, Particle1, temp_ClosestPoint_p0.To1DArray(), temp_ClosestPoint_p1.To1DArray(), NormalVector, Distance[p0, p1]);
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
                                        if (LastSaveTimeStep > m_dt * 1e-5)
                                            SaveTimeStep = -LastSaveTimeStep / 2;
                                        else
                                            SaveTimeStep = -m_dt * 0.1;
                                        MinDistance = double.MaxValue;
                                        Overlapping = true;
                                    }
                                    Particle0.skipForceIntegration = testMinDistance < 0.2 * m_hMin;
                                    Particle1.skipForceIntegration = testMinDistance < 0.2 * m_hMin;
                                }
                            }
                            if (SaveTimeStep != -m_dt)
                                AccDynamicTimestep += SaveTimeStep;
                            if (AccDynamicTimestep > m_dt)
                                break;
                        }
                        for (int p = 0; p < ParticlesOfCurrentColor.Length; p++)
                        {
                            Particles[ParticlesOfCurrentColor[p]].CollisionPreviousTimestep = AccDynamicTimestep;
                        }
                        if (AccDynamicTimestep > m_dt)
                            break;
                        if (MinDistance > MaxDistance && !Overlapping)
                            break;
                        for (int p0 = 0; p0 < ParticlesOfCurrentColor.Length; p0++)
                        {
                            Particle Particle0 = Particles[ParticlesOfCurrentColor[p0]];
                            for (int w = 0; w < 4; w++)
                            {
                                if (Distance[p0, ParticleOffset + w] < MaxDistance && SaveTimeStepArray[p0, ParticleOffset + w] > 0)
                                {
                                    double[] CurrentDistanceVector = DistanceVector.ExtractSubArrayShallow(new int[] { p0, ParticleOffset + w, -1 }).To1DArray();
                                    Particle0.closestPointToOtherObject = ClosestPoint_P0.ExtractSubArrayShallow(new int[] { p0, ParticleOffset + w, -1 }).To1DArray();
                                    CalculateNormalVector(CurrentDistanceVector, out double[] NormalVector);
                                    CalculateTangentialVector(NormalVector, out double[] TangentialVector);
                                    Particle0.CollisionNormalVector.Add(NormalVector);
                                    Particle0.CollisionTangentialVector.Add(TangentialVector);
                                    Particle0.CollisionTimestep = AccDynamicTimestep;
                                    Particle0.Collided = true;
                                    ComputeMomentumBalanceCollision(Particle0);
                                    PrintCollisionResults(Particles, Particle0, Distance[p0, ParticleOffset + w], AccDynamicTimestep);
                                }
                            }
                            for (int p1 = p0 + 1; p1 < ParticlesOfCurrentColor.Length; p1++)
                            {
                                Particle Particle1 = Particles[ParticlesOfCurrentColor[p1]];
                                if (Distance[p0, p1] < MaxDistance && SaveTimeStepArray[p0, p1] > 0)
                                {
                                    List<Particle> collidedParticles = new List<Particle>();
                                    collidedParticles.Add(Particle0);
                                    collidedParticles.Add(Particle1);
                                    double[] CurrentDistanceVector = DistanceVector.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                                    Particle0.closestPointToOtherObject = ClosestPoint_P0.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                                    Particle1.closestPointToOtherObject = ClosestPoint_P1.ExtractSubArrayShallow(new int[] { p0, p1, -1 }).To1DArray();
                                    CalculateNormalVector(CurrentDistanceVector, out double[] NormalVector);
                                    CalculateTangentialVector(NormalVector, out double[] TangentialVector);
                                    Particle0.CollisionNormalVector.Add(NormalVector);
                                    Particle1.CollisionNormalVector.Add(NormalVector);
                                    Particle0.CollisionTangentialVector.Add(TangentialVector);
                                    Particle1.CollisionTangentialVector.Add(TangentialVector);
                                    Particle0.CollisionTimestep = AccDynamicTimestep;
                                    Particle1.CollisionTimestep = AccDynamicTimestep;
                                    Particle0.Collided = true;
                                    Particle1.Collided = true;
                                    ComputeMomentumBalanceCollision(collidedParticles);
                                    PrintCollisionResults(Particles, Particle0, Particle1, Distance[p0, p1], AccDynamicTimestep);
                                }
                            }
                        }
                        for (int p = 0; p < ParticlesOfCurrentColor.Length; p++)
                        {
                            PostProcessCollisionTranslation(Particles[ParticlesOfCurrentColor[p]]);
                            PostProcessCollisionRotation(Particles[ParticlesOfCurrentColor[p]]);
                        }
                        PostCollisionTimestep = (m_dt - AccDynamicTimestep) * 2;
                        bool PostCollisionOverlapping = true;
                        while (PostCollisionOverlapping)
                        {
                            PostCollisionTimestep /= 2;
                            Console.WriteLine("PostCollisionTimestep " + PostCollisionTimestep);
                            for (int p = 0; p < ParticlesOfCurrentColor.Length; p++)
                            {
                                UpdateParticleState(Particles[ParticlesOfCurrentColor[p]], PostCollisionTimestep, spatialDim);
                            }
                            if (PostCollisionTimestep == 0)
                                break;
                            for (int p0 = 0; p0 < ParticlesOfCurrentColor.Length; p0++)
                            {
                                PostCollisionOverlapping = false;
                                Particle Particle0 = Particles[ParticlesOfCurrentColor[p0]];
                                CellMask ParticleBoundaryCells = GridData.GetBoundaryCells().Intersect(ParticleCutCells);
                                GetWall(GridData, ParticleBoundaryCells, out double[][] wallPoints);
                                for (int w = 0; w < wallPoints.GetLength(0); w++)
                                {
                                    Particle0.closestPointOnOtherObjectToThis = Particle0.Position[0].CloneAs();
                                    if (wallPoints[w] == null)
                                        continue;
                                    else if (wallPoints[w][0] != 0)
                                    {
                                        Particle0.closestPointOnOtherObjectToThis[0] = wallPoints[w][0];
                                    }
                                    else if (wallPoints[w][1] != 0)
                                    {
                                        Particle0.closestPointOnOtherObjectToThis[1] = wallPoints[w][1];
                                    }
                                    else
                                        continue;
                                    CalculateMinimumDistance(Particle0, null, out _, out _, out _, out _, out PostCollisionOverlapping);
                                }
                                if (PostCollisionOverlapping)
                                    break;
                                for (int p1 = p0 + 1; p1 < ParticlesOfCurrentColor.Length; p1++)
                                {
                                    Particle Particle1 = Particles[ParticlesOfCurrentColor[p1]];
                                    CalculateMinimumDistance(Particle0, Particle1, out _, out _, out _, out _, out PostCollisionOverlapping);
                                }
                                if (PostCollisionOverlapping)
                                    break;
                            }
                        }
                        AccDynamicTimestep += PostCollisionTimestep;
                        for (int p = 0; p < ParticlesOfCurrentColor.Length; p++)
                        {
                            Particles[ParticlesOfCurrentColor[p]].CollisionTimestep = AccDynamicTimestep;
                        }
                    }
                }
                for (int j = 0; j < GlobalParticleColor.Length; j++)
                {
                    if (GlobalParticleColor[j] == CurrentColor)
                        GlobalParticleColor[j] = 0;
                }
            }
        }

        internal void PrintCollisionResults(List<Particle> Particles, Particle particle0, Particle particle1, double Distance, double AccDynamicTimestep)
        {
            Console.WriteLine("-------------------------------------------------------");
            Console.WriteLine("Collision between particle #" + Particles.IndexOf(particle0) + " and particle #" + Particles.IndexOf(particle1));
            Console.WriteLine("The minimal distance is: " + Distance);
            Console.WriteLine("The dynamic timestep is: " + AccDynamicTimestep);
            Console.WriteLine("-------------------------------------------------------");
        }

        internal void PrintCollisionResults(List<Particle> Particles, Particle particle, double Distance, double AccDynamicTimestep)
        {
            Console.WriteLine("-------------------------------------------------------");
            Console.WriteLine("Collision between particle #" + Particles.IndexOf(particle) + " and particle the wall.");
            Console.WriteLine("The minimal distance is: " + Distance);
            Console.WriteLine("The dynamic timestep is: " + AccDynamicTimestep);
            Console.WriteLine("-------------------------------------------------------");
        }

        private void ModelCoefficientOfRestitution(double StokesNumber, out double CoefficientOfRestitution)
        {
            CoefficientOfRestitution = m_CoefficientOfRestitution;
            if (StokesNumber < 5 && StokesNumber != 0)
                CoefficientOfRestitution = 0;
        }

        /// <summary>
        /// Computes the post-collision velocities of two particles.
        /// </summary>
        /// <param name="collidedParticles">
        /// List of the two colliding particles
        /// </param>
        internal void ComputeMomentumBalanceCollision(List<Particle> collidedParticles)
        {
            ModelCoefficientOfRestitution(collidedParticles[0].ComputeParticleSt(m_FluidViscosity, m_FluidDensity, Aux.VectorDiff(collidedParticles[0].TranslationalVelocity[0], collidedParticles[1].TranslationalVelocity[0])), out m_CoefficientOfRestitution);

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
                double tempCollisionRot = collidedParticles[p].IncludeRotation ? collidedParticles[p].RotationalVelocity[0] + Math.Pow(-1, p) * collidedParticles[p].eccentricity * collisionCoefficient / collidedParticles[p].MomentOfInertia_P : 0;
                collidedParticles[p].CollisionTranslationalVelocity.Add(new double[] { tempCollisionVn, tempCollisionVt });
                collidedParticles[p].CollisionRotationalVelocity.Add(tempCollisionRot);
            }
        }

        /// <summary>
        /// Computes the post-collision velocities of one particle after the collision with a wall.
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
            double tempCollisionRot = particle.IncludeRotation ? particle.RotationalVelocity[0] + particle.eccentricity * collisionCoefficient / particle.MomentOfInertia_P : 0;
            particle.CollisionTranslationalVelocity.Add(new double[] { tempCollisionVn, tempCollisionVt });
            particle.CollisionRotationalVelocity.Add(tempCollisionRot);
        }

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
            collisionCoefficient += (1 + m_CoefficientOfRestitution) * ((-collidedParticles[0].eccentricity * collidedParticles[0].RotationalVelocity[0] + collidedParticles[1].eccentricity * collidedParticles[1].RotationalVelocity[0]) / (massReciprocal[0] + massReciprocal[1] + momentOfInertiaReciprocal[0] + momentOfInertiaReciprocal[1]));
        }

        private void CalculateCollisionCoefficient(Particle particle, out double collisionCoefficient)
        {
            collisionCoefficient = (1 + m_CoefficientOfRestitution) * (particle.PreCollisionVelocity[0] / (1 / particle.Mass_P  + particle.eccentricity.Pow2() / particle.MomentOfInertia_P));
            collisionCoefficient += -(1 + m_CoefficientOfRestitution) * particle.eccentricity * particle.RotationalVelocity[0] / (1 / particle.Mass_P + particle.eccentricity.Pow2() / particle.MomentOfInertia_P);
        }

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
                                    WallPoints[2] = new double[] { 0, WallPoint1[1]};
                                else if (WallPoints[3] == null || Math.Abs(WallPoint1[1] - WallPoints[3][1]) < 1e-12)
                                    WallPoints[3] = new double[] { 0, WallPoint1[1]};
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
            int SpatialDim = particle.Position[0].Length;
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

                particle.TranslationalVelocity.Insert(0, new double[2]);
                for (int d = 0; d < SpatialDim; d++)
                {
                    particle.TranslationalVelocity[0][d] = Normal[d] * temp_NormalVel + Tangential[d] * temp_TangentialVel;
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
                particle.RotationalVelocity[0] = particle.CollisionRotationalVelocity.Sum() / particle.CollisionRotationalVelocity.Count();
                particle.CollisionRotationalVelocity.Clear();

                if (double.IsNaN(particle.RotationalVelocity[0]) || double.IsInfinity(particle.RotationalVelocity[0]))
                    throw new ArithmeticException("Error trying to update particle angular velocity during collision post-processing. The angular velocity is:  " + particle.RotationalVelocity[0]);
            }
        }

        /// <summary>
        /// Computes the minimal distance between two particles or one particle and the wall.
        /// </summary>
        /// <param name="Particle0">
        /// The first particle.
        /// </param>
        ///  <param name="Particle1">
        /// The second particle, if Particle1 == null it is assumed to be a wall.
        /// </param>
        /// <param name="LsTrk">
        /// The level set tracker.
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
            int SpatialDim = Particle0.Position[0].Length;
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
            double[] Position0 = Particle0.Position[0].CloneAs();
            int SpatialDim = Position0.Length;
            double[] Position1 = Particle1 == null ? Particle0.closestPointOnOtherObjectToThis.CloneAs() : Particle1.Position[0].CloneAs();
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
        /// <param name="Position">
        /// Position vector of the particle.
        /// </param>
        /// <param name="Angle">
        /// Angle of the particle.
        /// </param>
        /// <param name="Vector">
        /// The vector in which direction the support point is searched.
        /// </param>
        /// <param name="lsTrk">
        /// The level set tracker.
        /// </param>
        /// <param name="SupportPoint">
        /// The support point (Cpt. Obvious)
        /// </param>
        private void CalculateSupportPoint(Particle particle, int SubParticleID, double[] Vector, out double[] SupportPoint)
        {
            int SpatialDim = particle.Position[0].Length;
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
        /// <param name="Simplex">
        /// A list of all support points constituting the simplex.
        /// </param>
        /// <param name="v">
        /// The distance vector.
        /// </param>
        /// <param name="Overlapping">
        /// Is true if the simplex contains the origin
        /// </param>
        private void DistanceAlgorithm(List<double[]> Simplex, out double[] v, out bool Overlapping)
        {
            int spatialDim = Simplex[0].Length;
            v = new double[spatialDim];
            Overlapping = false;

            // =======================================================
            // Step 1
            // Test for multiple Simplex-points 
            // and remove the duplicates
            // =======================================================
            for (int s1 = 0; s1 < Simplex.Count(); s1++)
            {
                for (int s2 = s1 + 1; s2 < Simplex.Count(); s2++)
                {
                    if (Math.Abs(Simplex[s1][0] - Simplex[s2][0]) < 1e-8 && Math.Abs(Simplex[s1][1] - Simplex[s2][1]) < 1e-8)
                    {
                        Simplex.RemoveAt(s2);
                    }
                }
            }

            // =======================================================
            // Step 2
            // Calculate dot product between all position vectors
            // =======================================================
            List<double[]> DotProd_Simplex = new List<double[]>();
            for (int s1 = 0; s1 < Simplex.Count(); s1++)
            {
                DotProd_Simplex.Add(new double[Simplex.Count()]);
                for (int s2 = s1; s2 < Simplex.Count(); s2++)
                {
                    DotProd_Simplex[s1][s2] = Aux.DotProduct(Simplex[s1], Simplex[s2]);
                }
            }

            // =======================================================
            // Step 3
            // Main test to determine the relatve position of
            // the simplex towards the origin.
            // =======================================================
            if (Simplex.Count() == 1)
            {
                v = Simplex[0];// the only possibility
                Aux.TestArithmeticException(v, nameof(v));
            }
            else if (Simplex.Count() == 2)
            {
                // The first simplex point is closest to the origin, choose this and delete the other one.
                if (DotProd_Simplex[0][0] - DotProd_Simplex[0][1] <= 0)
                {
                    v = Simplex[0].CloneAs();
                    Simplex.RemoveAt(1);
                    Aux.TestArithmeticException(v, nameof(v));
                }
                // The second simplex point is closest to the origin, choose this and delete the other one.
                else if (DotProd_Simplex[1][1] - DotProd_Simplex[0][1] <= 0)
                {
                    v = Simplex[1].CloneAs();
                    Simplex.RemoveAt(0);
                    Aux.TestArithmeticException(v, nameof(v));
                }
                // A point at the line between the two simplex point is closest to the origin, thus we need to keep both points.
                else
                {
                    double[] simplexDistanceVector = Aux.VectorDiff(Simplex[1], Simplex[0]); 
                    double Lambda = Math.Abs(-Simplex[1][1] * simplexDistanceVector[0] + Simplex[1][0] * simplexDistanceVector[1]) / (simplexDistanceVector[0].Pow2() + simplexDistanceVector[1].Pow2());
                    if(Lambda == 0)
                    {
                        Overlapping = true;
                    }
                    v[0] = -Lambda * simplexDistanceVector[1];
                    v[1] = Lambda * simplexDistanceVector[0];
                    Aux.TestArithmeticException(v, nameof(v));
                }
            }
            else if (Simplex.Count() == 3)
            {
                bool Return = false;
                for (int s1 = 0; s1 < Simplex.Count(); s1++)
                {
                    int s2 = s1 == 2 ? 2 : 1;
                    int s3 = s1 == 0 ? 0 : 1;
                    double test1 = DotProd_Simplex[s1][s1] - DotProd_Simplex[0][s2];
                    double test2 = DotProd_Simplex[s1][s1] - DotProd_Simplex[s3][2];
                    if (DotProd_Simplex[s1][s1] - DotProd_Simplex[0][s2] <= 0 && DotProd_Simplex[s1][s1] - DotProd_Simplex[s3][2] <= 0)
                    {
                        v = Simplex[s1].CloneAs();
                        Simplex.Clear();
                        Simplex.Add(v.CloneAs());
                        Return = true;
                        if (double.IsNaN(v[0]) || double.IsNaN(v[1]))
                            throw new ArithmeticException("Error trying to calculate v Value:  " + v[0] + " v " + v[1] + " Simplex count == 3.0");
                        break;
                    }
                }
                if (!Return)
                {
                    int counter = 0;
                    for (int s1 = Simplex.Count() - 1; s1 >= 0; s1--)
                    {
                        int s2 = s1 == 0 ? 1 : 2;
                        int s3 = s1 == 0 ? 2 : 0;
                        int s4 = s1 == 2 ? 1 : 0;
                        double CrossProd = new double();
                        switch (s1)
                        {
                            case 0:
                                double temp1 = DotProd_Simplex[1][2] - DotProd_Simplex[0][2] - DotProd_Simplex[1][1] + DotProd_Simplex[0][1];
                                double temp2 = DotProd_Simplex[0][1] - DotProd_Simplex[0][0] - DotProd_Simplex[1][2] + DotProd_Simplex[0][2];
                                double temp3 = DotProd_Simplex[1][1] - 2 * DotProd_Simplex[0][1] + DotProd_Simplex[0][0];
                                CrossProd = DotProd_Simplex[0][1] * temp1 + DotProd_Simplex[1][1] * temp2 + DotProd_Simplex[1][2] * temp3;
                                break;
                            case 1:
                                temp1 = -DotProd_Simplex[2][2] + DotProd_Simplex[0][2] + DotProd_Simplex[1][2] - DotProd_Simplex[0][1];
                                temp2 = DotProd_Simplex[2][2] - 2 * DotProd_Simplex[0][2] + DotProd_Simplex[0][0];
                                temp3 = DotProd_Simplex[0][2] - DotProd_Simplex[0][0] - DotProd_Simplex[1][2] + DotProd_Simplex[0][1];
                                CrossProd = DotProd_Simplex[0][2] * temp1 + DotProd_Simplex[1][2] * temp2 + DotProd_Simplex[2][2] * temp3;
                                break;
                            case 2:
                                temp1 = DotProd_Simplex[2][2] - 2 * DotProd_Simplex[1][2] + DotProd_Simplex[1][1];
                                temp2 = -DotProd_Simplex[2][2] + DotProd_Simplex[1][2] + DotProd_Simplex[0][2] - DotProd_Simplex[0][1];
                                temp3 = DotProd_Simplex[1][2] - DotProd_Simplex[1][1] - DotProd_Simplex[0][2] + DotProd_Simplex[0][1];
                                CrossProd = DotProd_Simplex[0][2] * temp1 + DotProd_Simplex[1][2] * temp2 + DotProd_Simplex[2][2] * temp3;
                                break;
                        }
                        CrossProd *= 1;
                        double test1 = DotProd_Simplex[s4][s4] - DotProd_Simplex[s4][s2];
                        double test2 = DotProd_Simplex[s2][s2] - DotProd_Simplex[s4][s2];
                        counter += 1;
                        if (DotProd_Simplex[s4][s4] - DotProd_Simplex[s4][s2] >= 0 && DotProd_Simplex[s2][s2] - DotProd_Simplex[s4][s2] >= 0 && CrossProd >= 0 && !Return)
                        {
                            double[] AB = new double[2];
                            for (int d = 0; d < 2; d++)
                            {
                                AB[d] = Simplex[s2][d] - Simplex[s4][d];
                            }
                            double Lambda = (Simplex[s2][1] * AB[0] - Simplex[s2][0] * AB[1]) / (AB[0].Pow2() + AB[1].Pow2());
                            v[0] = -Lambda * AB[1];
                            v[1] = Lambda * AB[0];
                            double[] test = new double[2];
                            test[1] = Simplex[s2][1] - Simplex[s4][1];
                            test[0] = (Simplex[s2][0] - Simplex[s4][0]);
                            double[] test12 = new double[2];
                            test12[0] = -test[1];
                            test12[1] = test[0];
                            double[] tempSimplex1 = Simplex[s2].CloneAs();
                            double[] tempSimplex2 = Simplex[s4].CloneAs();
                            Simplex.Clear();
                            Simplex.Add(tempSimplex1.CloneAs());
                            Simplex.Add(tempSimplex2.CloneAs());
                            Return = true;
                            if (double.IsNaN(v[0]) || double.IsNaN(v[1]))
                                throw new ArithmeticException("Error trying to calculate v Value:  " + v[0] + " v " + v[1] + " Simplex count == 3." + s1);
                            break;
                        }
                    }
                }
                if (!Return)
                {
                    Overlapping = true;
                }
            }
        }

        /// <summary>
        /// Ensures the communication between the processes after a collision
        /// </summary>
        /// <param name="Particles">
        /// A list of all particles.
        /// </param>
        /// <param name="CurrentParticle">
        /// The current particle.
        /// </param>
        /// <param name="MPISize">
        /// Number of mpi processes
        /// </param>
        /// <param name="WallCollision">
        /// If the collision was a wall collision set this true.
        /// </param>
        internal void Collision_MPICommunication(List<Particle> Particles, Particle CurrentParticle, int MPISize)
        {
            int NoOfVars = 13;
            double[] BoolSend = new double[1];
            double[] MPI_SkipForceIntegration = new double[1];
            bool NoCurrentCollision = true;
            BoolSend[0] = CurrentParticle.Collided ? 1 : 0;
            MPI_SkipForceIntegration[0] = CurrentParticle.skipForceIntegration ? 1 : 0;

            double[] Receive_SkipForceIntegration = new double[MPISize];
            unsafe
            {
                fixed (double* pCheckSend = MPI_SkipForceIntegration, pCheckReceive = Receive_SkipForceIntegration)
                {
                    csMPI.Raw.Allgather((IntPtr)pCheckSend, BoolSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, MPI_SkipForceIntegration.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                }
            }

            for(int i = 0; i < Receive_SkipForceIntegration.Length; i++)
            {
                if (Receive_SkipForceIntegration[i] == 1)
                    CurrentParticle.skipForceIntegration = true;
            }

            double[] BoolReceive = new double[MPISize];
            unsafe
            {
                fixed (double* pCheckSend = BoolSend, pCheckReceive = BoolReceive)
                {
                    csMPI.Raw.Allgather((IntPtr)pCheckSend, BoolSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, BoolSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                }
            }

            for (int i = 0; i < BoolReceive.Length; i++)
            {
                if (BoolReceive[i] != 0)
                {
                    double[] CheckSend = new double[NoOfVars];
                    CheckSend[0] = CurrentParticle.RotationalVelocity[0];
                    CheckSend[1] = CurrentParticle.TranslationalVelocity[0][0];
                    CheckSend[2] = CurrentParticle.TranslationalVelocity[0][1];
                    CheckSend[3] = CurrentParticle.Angle[0];
                    CheckSend[4] = CurrentParticle.Position[0][0];
                    CheckSend[5] = CurrentParticle.Position[0][1];
                    CheckSend[6] = CurrentParticle.CollisionTimestep;
                    CheckSend[7] = CurrentParticle.RotationalVelocity[1];
                    CheckSend[8] = CurrentParticle.TranslationalVelocity[1][0];
                    CheckSend[9] = CurrentParticle.TranslationalVelocity[1][1];
                    CheckSend[10] = CurrentParticle.Angle[1];
                    CheckSend[11] = CurrentParticle.Position[1][0];
                    CheckSend[12] = CurrentParticle.Position[1][1];

                    double[] CheckReceive = new double[NoOfVars * MPISize];
                    unsafe
                    {
                        fixed (double* pCheckSend = CheckSend, pCheckReceive = CheckReceive)
                        {
                            csMPI.Raw.Allgather((IntPtr)pCheckSend, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                        }
                    }
                    CurrentParticle.RotationalVelocity[0] = CheckReceive[0 + i * NoOfVars];
                    CurrentParticle.TranslationalVelocity[0][0] = CheckReceive[1 + i * NoOfVars];
                    CurrentParticle.TranslationalVelocity[0][1] = CheckReceive[2 + i * NoOfVars];
                    CurrentParticle.Angle[0] = CheckReceive[3 + i * NoOfVars];
                    CurrentParticle.Position[0][0] = CheckReceive[4 + i * NoOfVars];
                    CurrentParticle.Position[0][1] = CheckReceive[5 + i * NoOfVars];
                    CurrentParticle.CollisionTimestep = CheckReceive[6 + i * NoOfVars];
                    CurrentParticle.RotationalVelocity[1] = CheckReceive[7 + i * NoOfVars];
                    CurrentParticle.TranslationalVelocity[1][0] = CheckReceive[8 + i * NoOfVars];
                    CurrentParticle.TranslationalVelocity[1][1] = CheckReceive[9 + i * NoOfVars];
                    CurrentParticle.Angle[1] = CheckReceive[10 + i * NoOfVars];
                    CurrentParticle.Position[1][0] = CheckReceive[11 + i * NoOfVars];
                    CurrentParticle.Position[1][1] = CheckReceive[12 + i * NoOfVars];
                    //if (BoolReceive[i] != 0)
                    {
                        CurrentParticle.Collided = true;
                        NoCurrentCollision = false;
                    }
                }
            }
            if (NoCurrentCollision)
            {
                CurrentParticle.Collided = false;
            }
        }

        /// <summary>
        /// Obsolete, please use GJK-algorithm
        /// </summary>
        internal void FindClosestPoint(double hmin, int SpatialDim, MultidimensionalArray interfacePoints_P0, MultidimensionalArray interfacePoints_P1, ref double[] distanceVec, ref double distance, out double[] tempPoint_P0, out double[] tempPoint_P1, out bool Overlapping)
        {
            tempPoint_P0 = new double[SpatialDim];
            tempPoint_P1 = new double[SpatialDim];
            Overlapping = false;
            if (interfacePoints_P0 != null && interfacePoints_P1 != null)
            {
                for (int i = 0; i < interfacePoints_P0.NoOfRows; i++)
                {
                    for (int g = 0; g < interfacePoints_P1.NoOfRows; g++)
                    {
                        double tempDistance = Math.Sqrt((interfacePoints_P0.GetRow(i)[0] - interfacePoints_P1.GetRow(g)[0]).Pow2() + (interfacePoints_P0.GetRow(i)[1] - interfacePoints_P1.GetRow(g)[1]).Pow2());
                        if (tempDistance < distance)
                        {
                            distanceVec = interfacePoints_P0.GetRow(i).CloneAs();
                            distanceVec.AccV(-1, interfacePoints_P1.GetRow(g));
                            tempPoint_P0 = interfacePoints_P0.GetRow(i);
                            tempPoint_P1 = interfacePoints_P1.GetRow(g);
                            distance = tempDistance;
                        }
                        if (tempDistance < 0.1 * hmin)
                        {
                            Overlapping = true;
                        }
                    }
                }
            }
            Console.WriteLine("Overlapping = " + Overlapping);
        }

        internal double DynamicTimestep(Particle particle0, Particle particle1, double[] ClosestPoint0, double[] ClosestPoint1, double[] NormalVector, double Distance)
        {
            particle0.CalculateRadialVector(ClosestPoint0, out double[] RadialVector0, out double RadialLength0);
            double[] PointVelocity0 = new double[2];
            PointVelocity0[0] = particle0.TranslationalVelocity[0][0] - particle0.RotationalVelocity[0] * RadialLength0 * RadialVector0[1];
            PointVelocity0[1] = particle0.TranslationalVelocity[0][1] + particle0.RotationalVelocity[0] * RadialLength0 * RadialVector0[0];
            ProjectVelocityOnVector(NormalVector, PointVelocity0, out double DetectCollisionVn_P0);
            double Dynamic_dt;
            if (particle1 != null)
            {
                double[] PointVelocity1 = new double[2];
                particle1.CalculateRadialVector(ClosestPoint1, out double[] RadialVector1, out double RadialLength1);
                PointVelocity1[0] = particle1.TranslationalVelocity[0][0] - particle1.RotationalVelocity[0] * RadialLength1 * RadialVector1[1];
                PointVelocity1[1] = particle1.TranslationalVelocity[0][1] + particle1.RotationalVelocity[0] * RadialLength1 * RadialVector1[0];
                ProjectVelocityOnVector(NormalVector, PointVelocity1, out double DetectCollisionVn_P1);
                if (DetectCollisionVn_P1 - DetectCollisionVn_P0 == 0)
                    return double.MaxValue;
                Dynamic_dt = 0.9 * Distance / (DetectCollisionVn_P1 - DetectCollisionVn_P0);
            }
            else if (DetectCollisionVn_P0 == 0)
                return double.MaxValue;
            else
                Dynamic_dt = 0.9 * Distance / (-DetectCollisionVn_P0);

            return Dynamic_dt;
        }

        internal void UpdateParticleState(Particle particle, double Dynamic_dt, int SpatialDim)
        {
            if (Dynamic_dt != 0)
            {
                for (int d = 0; d < SpatialDim; d++)
                {
                    particle.Position[0][d] = particle.Position[0][d] + (particle.TranslationalVelocity[1][d] * 0 + 2 * particle.TranslationalVelocity[0][d]) * Dynamic_dt / 2;
                    if (double.IsNaN(particle.Position[0][d]) || double.IsInfinity(particle.Position[0][d]))
                        throw new ArithmeticException("Error trying to update particle position. Value:  " + particle.Position[0][d]);
                }
                particle.Angle[0] = particle.Angle[0] + (particle.RotationalVelocity[1] * 0 + 2 * particle.RotationalVelocity[0]) * Dynamic_dt / 2;
                if (double.IsNaN(particle.Angle[0]) || double.IsInfinity(particle.Angle[0]))
                    throw new ArithmeticException("Error trying to update particle position. Value:  " + particle.Angle[0]);
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
        /// Calculates the point velocity due to the rotaional velocity of the particle.
        /// </summary>
        /// <param name="RotationalVelocity">
        /// </param>
        /// <param name="RadialLength">
        /// </param>
        /// <param name="RadialNormalVector">
        /// Vector normal to the radial vector.
        /// </param>
        /// <param name="PointVelocityDueToRotation">
        /// </param>
        internal void TransformRotationalVelocity(double RotationalVelocity, double RadialLength, double[] RadialNormalVector, out double[] PointVelocityDueToRotation)
        {
            PointVelocityDueToRotation = new double[RadialNormalVector.Length];
            for (int d = 0; d < RadialNormalVector.Length; d++)
            {
                PointVelocityDueToRotation[d] = RadialLength * RotationalVelocity * RadialNormalVector[d];
            }
        }

    }
}

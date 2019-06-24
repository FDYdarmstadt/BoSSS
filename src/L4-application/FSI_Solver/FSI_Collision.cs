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
        private FSI_Auxillary Aux = new FSI_Auxillary();

        /// <summary>
        /// Computes the minimal distance between two particles or one particle and the wall.
        /// </summary>
        /// <param name="Particle0">
        /// The first particle.
        /// </param>
        ///  <param name="Particle1">
        /// The second particle, if Particle1 == null it is assumed to be a wall.
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
        /// <param name="CoefficientOfRestitution">
        /// Coefficient of restitution.
        /// </param>
        internal void ComputeMomentumBalanceCollision(Particle Particle0, Particle Particle1, double[] DistanceVector, double[] ClosestPoint_P0, double[] ClosestPoint_P1, double CoefficientOfRestitution)
        {
            CalculateNormalAndTangentialVector(DistanceVector, out double[] NormalVector, out double[] TangentialVector);
            ProjectVelocity(NormalVector, TangentialVector, Particle0.TranslationalVelocity[0], out double collisionVn_P0, out double collisionVt_P0);
            ProjectVelocity(NormalVector, TangentialVector, Particle1.TranslationalVelocity[0], out double collisionVn_P1, out double collisionVt_P1);

            // Bool if collided
            //Particle0.m_collidedWithParticle[m_Particles.IndexOf(Particle1)] = true;
            //Particle1.m_collidedWithParticle[m_Particles.IndexOf(Particle0)] = true;

            // Bool if force integration should be skipped
            Particle0.skipForceIntegration = true;
            Particle1.skipForceIntegration = true;

            // coefficient of restitution (e=0 pastic; e=1 elastic)
            double e = CoefficientOfRestitution;

            // Calculate excentric parameter
            double[] RadialDistance0 = Aux.VectorDiff(ClosestPoint_P0, Particle0.Position[0]);
            double[] RadialDistance1 = Aux.VectorDiff(ClosestPoint_P1, Particle1.Position[0]);
            double a0 = Particle0 is Particle_Sphere ? 0.0 : RadialDistance0[0] * TangentialVector[0] + RadialDistance0[1] * TangentialVector[1];
            double a1 = Particle1 is Particle_Sphere ? 0.0 : RadialDistance1[0] * TangentialVector[0] + RadialDistance1[1] * TangentialVector[1];

            // Calculate post collision velocities.
            double Fx;
            double Fxrot;
            double tempCollisionVn_P0;
            double tempCollisionVn_P1;
            double tempCollisionRot_P0 = 0;
            double tempCollisionRot_P1 = 0;
            if (!Particle0.IncludeTranslation && !Particle0.IncludeRotation)
            {
                Fx = (1 + e) * ((collisionVn_P1) / (1 / Particle1.Mass_P + a1.Pow2() / Particle1.MomentOfInertia_P));
                Fxrot = (1 + e) * ((a1 * Particle1.RotationalVelocity[0]) / (1 / Particle1.Mass_P + a1.Pow2() / Particle1.MomentOfInertia_P));
                tempCollisionVn_P0 = collisionVn_P0;
                tempCollisionVn_P1 = collisionVn_P1 + (Fx + Fxrot) / Particle1.Mass_P;
                tempCollisionRot_P1 = Particle1.RotationalVelocity[0] - a1 * (Fx + Fxrot) / Particle1.MomentOfInertia_P;
            }
            else if (!Particle1.IncludeTranslation && !Particle1.IncludeRotation)
            {
                Fx = (1 + e) * ((collisionVn_P0) / (1 / Particle0.Mass_P + a0.Pow2() / Particle0.MomentOfInertia_P));
                Fxrot = (1 + e) * ((-a0 * Particle0.RotationalVelocity[0]) / (1 / Particle0.Mass_P + a0.Pow2() / Particle0.MomentOfInertia_P));
                tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / Particle0.Mass_P;
                tempCollisionRot_P0 = Particle0.RotationalVelocity[0] + a0 * (Fx + Fxrot) / Particle0.MomentOfInertia_P;
                tempCollisionVn_P1 = collisionVn_P1;
            }
            else
            {
                Fx = (1 + e) * ((collisionVn_P0 - collisionVn_P1) / (1 / Particle0.Mass_P + 1 / Particle1.Mass_P + a0.Pow2() / Particle0.MomentOfInertia_P + a1.Pow2() / Particle1.MomentOfInertia_P));
                Fxrot = (1 + e) * ((-a0 * Particle0.RotationalVelocity[0] + a1 * Particle1.RotationalVelocity[0]) / (1 / Particle0.Mass_P + 1 / Particle1.Mass_P + a0.Pow2() / Particle0.MomentOfInertia_P + a1.Pow2() / Particle1.MomentOfInertia_P));
                tempCollisionVn_P0 = collisionVn_P0 - (Fx + Fxrot) / Particle0.Mass_P;
                tempCollisionRot_P0 = Particle0.RotationalVelocity[0] + a0 * (Fx + Fxrot) / Particle0.MomentOfInertia_P;
                tempCollisionVn_P1 = collisionVn_P1 + (Fx + Fxrot) / Particle1.Mass_P;
                tempCollisionRot_P1 = Particle1.RotationalVelocity[0] - a1 * (Fx + Fxrot) / Particle1.MomentOfInertia_P;
            }

            double tempCollisionVt_P0 = collisionVt_P0 * e;
            double tempCollisionVt_P1 = collisionVt_P1 * e;
            Console.WriteLine("a0:    " + a0 + "   Fx:    " + (-Fx) + "      Fxrot:    " + (-Fxrot));
            Console.WriteLine("a1:    " + a1 + "   Fx:    " + Fx + "      Fxrot:    " + Fxrot);

            Particle0.CollisionNormal.Add(NormalVector);
            Particle1.CollisionNormal.Add(NormalVector);
            Particle0.CollisionTangential.Add(TangentialVector);
            Particle1.CollisionTangential.Add(TangentialVector);
            Particle0.CollisionRotationalVelocity.Add(tempCollisionRot_P0);
            Particle1.CollisionRotationalVelocity.Add(tempCollisionRot_P1);
            Particle0.CollisionTranslationalVelocity.Add(new double[] { tempCollisionVn_P0, tempCollisionVt_P0 });
            Particle1.CollisionTranslationalVelocity.Add(new double[] { tempCollisionVn_P1, tempCollisionVt_P1 });


            for (int d = 0; d < 2; d++)
            {
                Particle0.TranslationalAcceleration[1][d] = 0;
                Particle0.RotationalAcceleration[1] = 0;
                Particle0.RotationalVelocity[1] = 0;
                Particle1.TranslationalAcceleration[1][d] = 0;
                Particle1.RotationalAcceleration[1] = 0;
                Particle1.RotationalVelocity[1] = 0;
            }
        }

        /// <summary>
        /// Collision post-processing. Sums up the results of the multiple binary collisions of one timestep
        /// </summary>
        /// <param name="_Particle">
        /// The particle to be processed
        /// </param>
        internal void SumOverCollisionVelocities(Particle _Particle)
        {
            int SpatialDim = _Particle.Position[0].Length;
            if (_Particle.CollisionRotationalVelocity.Count() >= 1)
            {
                _Particle.RotationalVelocity[0] = 0;
                for (int r = 0; r < _Particle.CollisionRotationalVelocity.Count(); r++)
                {
                    _Particle.RotationalVelocity[0] += _Particle.CollisionRotationalVelocity[r];
                }
            }
            if (_Particle.CollisionTranslationalVelocity.Count() >= 1)
            {
                double[] Normal = new double[SpatialDim];
                double[] Tangential = new double[SpatialDim];
                for (int t = 0; t < _Particle.CollisionTranslationalVelocity.Count(); t++)
                {
                    for (int d = 0; d < SpatialDim; d++)
                    {
                        Normal[d] += _Particle.CollisionNormal[t][d];
                        Tangential[d] += _Particle.CollisionTangential[t][d];
                    }
                }

                Normal.ScaleV(1 / Math.Sqrt(Normal[0].Pow2() + Normal[1].Pow2()));
                Tangential.ScaleV(1 / Math.Sqrt(Tangential[0].Pow2() + Tangential[1].Pow2()));
                double[] Cos = new double[_Particle.CollisionTranslationalVelocity.Count()];
                double[] Sin = new double[_Particle.CollisionTranslationalVelocity.Count()];
                double temp_NormalVel = 0;
                double temp_TangentialVel = 0;
                for (int t = 0; t < _Particle.CollisionTranslationalVelocity.Count(); t++)
                {
                    for (int d = 0; d < SpatialDim; d++)
                    {
                        Cos[t] += Normal[d] * _Particle.CollisionNormal[t][d];
                    }
                    Sin[t] = Cos[t] == 1 ? 0 : _Particle.CollisionNormal[t][0] > Normal[0] ? Math.Sqrt(1 + 1e-15 - Cos[t].Pow2()) : -Math.Sqrt(1 + 1e-15 - Cos[t].Pow2());
                    temp_NormalVel += _Particle.CollisionTranslationalVelocity[t][0] * Cos[t] - _Particle.CollisionTranslationalVelocity[t][1] * Sin[t];
                    temp_TangentialVel += _Particle.CollisionTranslationalVelocity[t][0] * Sin[t] + _Particle.CollisionTranslationalVelocity[t][1] * Cos[t];

                }
                temp_NormalVel /= _Particle.CollisionTranslationalVelocity.Count();
                temp_TangentialVel /= _Particle.CollisionTranslationalVelocity.Count();
                _Particle.TranslationalVelocity.Insert(0, new double[2]);
                for (int d = 0; d < SpatialDim; d++)
                {
                    _Particle.TranslationalVelocity[0][d] = Normal[d] * temp_NormalVel + Tangential[d] * temp_TangentialVel;
                }
                _Particle.CollisionTranslationalVelocity.Clear();
                _Particle.CollisionNormal.Clear();
                _Particle.CollisionTangential.Clear();
                _Particle.CollisionPositionCorrection.Clear();
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
        internal void ComputeMinimalDistance(Particle Particle0, Particle Particle1, LevelSetTracker LsTrk, out double Distance, out MultidimensionalArray DistanceVector, out MultidimensionalArray ClosestPoint_P0, out MultidimensionalArray ClosestPoint_P1, out bool Overlapping)
        {
            int SpatialDim = Particle0.Position[0].Length;
            Distance = double.MaxValue;
            DistanceVector = MultidimensionalArray.Create(SpatialDim);
            ClosestPoint_P0 = MultidimensionalArray.Create(SpatialDim);
            ClosestPoint_P1 = MultidimensionalArray.Create(SpatialDim);
            Overlapping = false;

            for (int i = 0; i < Particle0.NoOfSubParticles(); i++)
            {
                for (int j = 0; j < Particle1.NoOfSubParticles(); j++)
                {
                    GJK_DistanceAlgorithm(Particle0, i, Particle1, j, LsTrk, Particle0.Position[0], Particle1.Position[0], Particle0.Angle[0], Particle1.Angle[0], out double temp_Distance, out double[] temp_DistanceVector, out double[] temp_ClosestPoint_P0, out double[] temp_ClosestPoint_P1, out Overlapping);
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
        /// <param name="lsTrk">
        /// The level set tracker.
        /// </param>
        /// <param name="Position0">
        /// Position vector of the first object.
        /// </param>
        /// <param name="Position1">
        /// Position vector of the second object.
        /// </param>
        /// <param name="Angle0">
        /// Angle of the first object.
        /// </param>
        /// <param name="Angle1">
        /// Angle of the second object.
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
        internal void GJK_DistanceAlgorithm(Particle Particle0, int SubParticleID0, Particle Particle1, int SubParticleID1, LevelSetTracker lsTrk, double[] Position0, double[] Position1, double Angle0, double Angle1, out double Min_Distance, out double[] DistanceVec, out double[] ClosestPoint0, out double[] ClosestPoint1, out bool Overlapping)
        {
            int SpatialDim = Position0.Length;
            ClosestPoint0 = new double[SpatialDim];
            ClosestPoint1 = new double[SpatialDim];
            double[] SupportPoint = new double[SpatialDim];
            DistanceVec = new double[SpatialDim];
            Overlapping = false;

            Initialize_GJK(Position0, Position1, out double[] v, out List<double[]> Simplex);

            for (int i = 0; i < 10000; i++)
            {
                double[] vt = v.CloneAs();
                for (int d = 0; d < SpatialDim; d++)
                {
                    vt[d] = -v[d];
                }
                //double[] vt = Aux.VectorDiff(null, v);
                if (double.IsNaN(vt[0]) || double.IsNaN(vt[1]))
                    throw new ArithmeticException("Error trying to calculate point0 Value:  " + vt[0] + " point1 " + vt[1]);

                // =======================================================
                // Step 1
                // Calculate the support point of the minkowski difference
                // and of the two particles (which are the closest points
                // if the algorithm is finished.
                // =======================================================
                CalculateSupportPoint(Particle0, SubParticleID0, Position0, Angle0, vt, lsTrk, out ClosestPoint0);
                if (double.IsNaN(ClosestPoint0[0]) || double.IsNaN(ClosestPoint0[1]))
                    throw new ArithmeticException("Error trying to calculate point0 Value:  " + ClosestPoint0[0] + " point1 " + ClosestPoint0[1]);
                if (Particle1 != null)
                    CalculateSupportPoint(Particle1, SubParticleID1, Position1, Angle1, v, lsTrk, out ClosestPoint1);
                else
                {
                    ClosestPoint1 = ClosestPoint0.CloneAs();    
                    if (Position0[0] == Position1[0])
                        ClosestPoint1[1] = Position1[1];
                    else
                        ClosestPoint1[0] = Position1[0];
                }
                SupportPoint = Aux.VectorDiff(ClosestPoint0, ClosestPoint1);

                // =======================================================
                // Step 2
                // Check max(x dot vt)
                // =======================================================
                if ((Aux.DotProduct(v,vt) - Aux.DotProduct(SupportPoint,vt)) >= -1e-12 && i != 0)
                {
                    Console.WriteLine("No of steps for distance algorithm: " + i);
                    break;
                }

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

                // End algorithm if the two object are overlapping.
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
        /// Inititalizes the GJK-algorithm
        /// </summary>
        /// <param name="Position0">
        /// Position vector of the first object.
        /// </param>
        /// <param name="Position1">
        /// Position vector of the second object.
        /// </param>
        /// <param name="v0">
        /// Initial guess for the min distance vector (diff. between the two position vectors)
        /// </param>
        /// <param name="Simplex">
        /// List of all support points defining the simplex. Initially it contains only v0.
        /// </param>
        private void Initialize_GJK(double[] Position0, double[] Position1, out double[] v0, out List<double[]> Simplex)
        {
            Simplex = new List<double[]>();
            v0 = Aux.VectorDiff(Position0, Position1);
            Simplex.Add(v0.CloneAs());
        }

        /// <summary>
        /// Calculates the support point on a single particle.
        /// </summary>
        /// <param name="_Particle">
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
        private void CalculateSupportPoint(Particle _Particle, int SubParticleID, double[] Position, double Angle, double[] Vector, LevelSetTracker lsTrk, out double[] SupportPoint)
        {
            int SpatialDim = Position.Length;
            SupportPoint = new double[SpatialDim];
            // A direct formulation of the support function for a sphere exists, thus it is possible to map it to an ellipsoid.
            if (_Particle is Particle_Ellipsoid || _Particle is Particle_Sphere)
            {
                _Particle.GetSupportPoint(SpatialDim, Vector, Position, Angle, out SupportPoint);
            }
            // Binary search in all other cases.
            else
            {
                MultidimensionalArray SurfacePoints = _Particle.GetSurfacePoints(lsTrk, Position, Angle);
                MultidimensionalArray SurfacePointsSubParticle = SurfacePoints.ExtractSubArrayShallow(new int[]{ SubParticleID, -1, -1});
                int L = 1;
                int R = SurfacePointsSubParticle.GetLength(0) - 2;
                int Counter = 0;
                while (L <= R && L > 0 && R < SurfacePointsSubParticle.GetLength(0) - 1)
                {
                    int Index = (L + R) / 2;
                    Counter = Counter + 1;
                    GetPointAndNeighbours(SurfacePointsSubParticle, Index, out SupportPoint, out double[] RightNeighbour, out double[] LeftNeighbour);
                    double DotSupportPoint = SupportPoint[0] * Vector[0] + SupportPoint[1] * Vector[1];
                    double DotRight = RightNeighbour[0] * Vector[0] + RightNeighbour[1] * Vector[1];
                    double DotLeft = LeftNeighbour[0] * Vector[0] + LeftNeighbour[1] * Vector[1];
                    if (DotSupportPoint > DotRight && DotSupportPoint > DotLeft)
                        break;
                    else if (DotRight > DotLeft)
                        L = Index + 1;
                    else
                        R = Index - 1;
                }
            }
        }

        /// <summary>
        /// Searchs for a specific point and its neighbours on a particle surface.
        /// </summary>
        /// <param name="SurfacePoints">
        /// All surface points of the current particle
        /// </param>
        /// <param name="Index">
        /// Index of the current point.
        /// </param>
        /// <param name="Point">
        /// Coordinates of the current point.
        /// </param>
        /// <param name="RightNeighbour">
        /// Its right neighbour (Index + 1)
        /// </param>
        /// <param name="LeftNeighbour">
        /// Its left neighbour (Index - 1)
        /// </param>
        private void GetPointAndNeighbours(MultidimensionalArray SurfacePoints, int Index, out double[] Point, out double[] RightNeighbour, out double[] LeftNeighbour)
        {
            Point = new double[2];
            RightNeighbour = new double[2];
            LeftNeighbour = new double[2];
            for (int d = 0; d < 2; d++)
            {
                Point[d] = SurfacePoints[Index, d];
                LeftNeighbour[d] = SurfacePoints[Index - 1, d];
                RightNeighbour[d] = SurfacePoints[Index + 1, d];
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
            v = new double[2];
            Overlapping = false;
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
            List<double[]> DotProd_Simplex = new List<double[]>();
            for (int s1 = 0; s1 < Simplex.Count(); s1++)
            {
                DotProd_Simplex.Add(new double[Simplex.Count()]);
                for (int s2 = s1; s2 < Simplex.Count(); s2++)
                {
                    DotProd_Simplex[s1][s2] = Simplex[s1][0] * Simplex[s2][0] + Simplex[s1][1] * Simplex[s2][1];
                }
            }
            if (Simplex.Count() == 1)
            {
                v = Simplex[0];
            }
            else if (Simplex.Count() == 2)
            {
                if (DotProd_Simplex[0][0] - DotProd_Simplex[0][1] <= 0)
                {
                    v = Simplex[0].CloneAs();
                    Simplex.RemoveAt(1);
                }

                else if (DotProd_Simplex[1][1] - DotProd_Simplex[0][1] <= 0)
                {
                    v = Simplex[1].CloneAs();
                    Simplex.RemoveAt(0);
                }
                else
                {
                    double[] AB = new double[2];
                    for (int d = 0; d < 2; d++)
                    {
                        AB[d] = Simplex[1][d] - Simplex[0][d];
                    }
                    double Lambda = (Simplex[1][1] * AB[0] - Simplex[1][0] * AB[1]) / (AB[0].Pow2() + AB[1].Pow2());
                    v[0] = -Lambda * AB[1];
                    v[1] = Lambda * AB[0];
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
                            break;
                        }
                        if (counter == 3)
                        {
                            Console.WriteLine("Warning");
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
        internal void Collision_MPICommunication(List<Particle> Particles, Particle CurrentParticle, int MPISize, bool WallCollision = false)
        {
            int NoOfVars = 3;
            double[] BoolSend = new double[1];
            bool NoCurrentCollision = true;
            if (CurrentParticle.m_collidedWithWall[0] && WallCollision)
                BoolSend[0] = -1;
            else
            {
                for (int p = 0; p < Particles.Count(); p++)
                {
                    if (CurrentParticle.m_collidedWithParticle[p])
                        BoolSend[0] = p + 1;
                }
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

                    double[] CheckReceive = new double[NoOfVars * MPISize];
                    unsafe
                    {
                        fixed (double* pCheckSend = CheckSend, pCheckReceive = CheckReceive)
                        {
                            csMPI.Raw.Allgather((IntPtr)pCheckSend, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, CheckSend.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                        }
                    }
                    CurrentParticle.RotationalVelocity[0] = CheckReceive[0 + i * 3];
                    CurrentParticle.TranslationalVelocity[0][0] = CheckReceive[1 + i * 3];
                    CurrentParticle.TranslationalVelocity[0][1] = CheckReceive[2 + i * 3];
                    if (!WallCollision)
                    {
                        int p = Convert.ToInt32(BoolReceive[i]);
                        CurrentParticle.m_collidedWithParticle[p - 1] = true;
                        CurrentParticle.skipForceIntegration = true;
                    }
                    NoCurrentCollision = false;
                }
            }
            if (NoCurrentCollision)
            {
                //CurrentParticle.skipForceIntegration = false;
                CurrentParticle.m_collidedWithWall[0] = false;
                for (int p = 0; p < CurrentParticle.m_collidedWithParticle.Length; p++)
                {
                    CurrentParticle.m_collidedWithParticle[p] = false;
                }
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

        /// <summary>
        /// Calculates the dynamic collision threshold based on the normal velocity of the two 
        /// particles.
        /// </summary>
        /// <param name="particle0">
        /// </param>
        /// <param name="particle1">
        /// </param>
        /// <param name="ClosestPoint0">
        /// The closest point on particle0 to particle1.
        /// </param>
        /// <param name="ClosestPoint1">
        /// The closest point on particle1 to particle0.
        /// </param>
        /// <param name="NormalVector">
        /// The collision normal vector.
        /// </param>
        /// <param name="Distance">
        /// The minimum distance between the two particles.
        /// </param>
        /// <param name="dt">
        /// The timestep.
        /// </param>
        /// <param name="Threshold">
        /// </param>
        internal void CalculateDynamicCollisionThreshold(Particle particle0, Particle particle1, double[] ClosestPoint0, double[] ClosestPoint1, double[] NormalVector, double Distance, double dt, out double Threshold)
        {
            Threshold = 0;
            CalculateRadialVector(particle0.Position[0], ClosestPoint0, out _, out double RadialLength0, out double[] RadialNormalVector0);
            TransformRotationalVelocity(particle0.RotationalVelocity[0], RadialLength0, RadialNormalVector0, out double[] PointVelocityDueToRotation0);
            double[] PointVelocity0 = new double[2];
            for (int d = 0; d < 2; d++)
                PointVelocity0[d] = particle0.TranslationalVelocity[0][d] + PointVelocityDueToRotation0[d];
            ProjectVelocityOnVector(NormalVector, PointVelocity0, out double DetectCollisionVn_P0);
            if (particle1 != null)
            {
                CalculateRadialVector(particle1.Position[0], ClosestPoint1, out _, out double RadialLength1, out double[] RadialNormalVector1);
                TransformRotationalVelocity(particle1.RotationalVelocity[0], RadialLength1, RadialNormalVector1, out double[] PointVelocityDueToRotation1);
                double[] PointVelocity1 = new double[2];
                for (int d = 0; d < 2; d++)
                {
                    PointVelocity1[d] = particle1.TranslationalVelocity[0][d] + PointVelocityDueToRotation1[d];
                }
                ProjectVelocityOnVector(NormalVector, PointVelocity1, out double DetectCollisionVn_P1);
                if (Distance <= (-DetectCollisionVn_P0 + DetectCollisionVn_P1) * dt || Math.Abs((-DetectCollisionVn_P0 + DetectCollisionVn_P1)) <= 1e-14)
                {
                    Threshold = Math.Abs(-DetectCollisionVn_P0 + DetectCollisionVn_P1) * dt;
                }
            }
            else if (Distance <= (-DetectCollisionVn_P0) * dt || Math.Abs(DetectCollisionVn_P0) <= 1e-14)
            {
                Threshold = Math.Abs(DetectCollisionVn_P0) * dt;
            }
        }

        internal double DynamicTimestep(Particle particle0, Particle particle1, double[] ClosestPoint0, double[] ClosestPoint1, double[] NormalVector, double Distance)
        {
            double Dynamic_dt = 213;
            double rMax_0 = particle0.GetLengthScales().Max();
            double[] PointVelocity0 = new double[2];
            for (int d = 0; d < 2; d++)
                PointVelocity0[d] = particle0.TranslationalVelocity[0][d] + particle0.RotationalVelocity[0] * rMax_0;
            ProjectVelocityOnVector(NormalVector, PointVelocity0, out double DetectCollisionVn_P0);
            if (particle1 != null)
            {
                double rMax_1 = particle1.GetLengthScales().Max();
                double[] PointVelocity1 = new double[2];
                for (int d = 0; d < 2; d++)
                    PointVelocity1[d] = particle1.TranslationalVelocity[0][d] + particle1.RotationalVelocity[0] * rMax_1;
                ProjectVelocityOnVector(NormalVector, PointVelocity1, out double DetectCollisionVn_P1);
                Console.WriteLine("NormalVector0: " + NormalVector[0]);
                Console.WriteLine("NormalVector1: " + NormalVector[1]);
                Console.WriteLine("DetectCollisionVn_P0: " + DetectCollisionVn_P0);
                Console.WriteLine("DetectCollisionVn_P1: " + DetectCollisionVn_P1);
                if (DetectCollisionVn_P1 - DetectCollisionVn_P0 == 0)
                    return double.MaxValue;
                Dynamic_dt = 0.99 * Distance / (DetectCollisionVn_P1 - DetectCollisionVn_P0);
                Console.WriteLine("Dynamic_dt: " + Dynamic_dt);
            }
            else if(DetectCollisionVn_P0 == 0)
                return double.MaxValue;
            else
                Dynamic_dt = 0.99 * Distance / (-DetectCollisionVn_P0);
                
            return Dynamic_dt;
        }

        internal void UpdateParticleState(Particle particle, double Dynamic_dt, int SpatialDim)
        {
            if (Dynamic_dt < 0 || Dynamic_dt == double.MaxValue)
            {
                for (int d = 0; d < SpatialDim; d++)
                {
                    particle.Position[0][d] = particle.Position[1][d];
                }
                particle.Angle[0] = particle.Angle[1];
            }
            else
            {
                for (int d = 0; d < SpatialDim; d++)
                {
                    particle.Position[0][d] = particle.Position[0][d] + (particle.TranslationalVelocity[1][d] + particle.TranslationalVelocity[0][d]) * Dynamic_dt / 2 + (particle.TranslationalAcceleration[1][d] + particle.TranslationalAcceleration[0][d]) * Dynamic_dt.Pow2() / 4;
                    if (double.IsNaN(particle.Position[0][d]) || double.IsInfinity(particle.Position[0][d]))
                        throw new ArithmeticException("Error trying to update particle position. Value:  " + particle.Position[0][d]);
                }
                particle.Angle[0] = particle.Angle[0] + (particle.RotationalVelocity[1] + particle.RotationalVelocity[0]) * Dynamic_dt / 2 + (particle.RotationalAcceleration[0] + particle.RotationalAcceleration[1]) * Dynamic_dt.Pow2() / 4;
            }
        }

        /// <summary>
        /// Calculates the normal and tangential vector from the min distance vector between the
        /// two particles.
        /// </summary>
        /// <param name="distanceVec">
        /// The min distance vector between the two particles.
        /// </param>
        /// <param name="NormalVector">
        /// The collision normal vector.
        /// </param>
        /// <param name="TangentialVector">
        /// The collision tangential vector.
        /// </param>
        internal void CalculateNormalAndTangentialVector(double[] distanceVec, out double[] NormalVector, out double[] TangentialVector)
        {
            NormalVector = distanceVec.CloneAs();
            NormalVector.ScaleV(1 / Math.Sqrt(distanceVec[0].Pow2() + distanceVec[1].Pow2()));
            TangentialVector = new double[] { -NormalVector[1], NormalVector[0] };
        }

        /// <summary>
        /// Calculates the normal and tangential components of the translational velocity.
        /// </summary>
        /// <param name="NormalVector">
        /// </param>
        /// <param name="TangentialVector">
        /// </param>
        /// <param name="TranslationalVelocity">
        /// </param>
        /// <param name="VelocityNormal">
        /// </param>
        /// <param name="VelocityTangential">
        /// </param>
        internal void ProjectVelocity(double[] NormalVector, double[] TangentialVector, double[] TranslationalVelocity, out double VelocityNormal, out double VelocityTangential)
        {
            ProjectVelocityOnVector(NormalVector, TranslationalVelocity, out VelocityNormal);
            ProjectVelocityOnVector(TangentialVector, TranslationalVelocity, out VelocityTangential);
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
        /// Calculates the radial vector (SurfacePoint-ParticlePosition)
        /// </summary>
        /// <param name="ParticlePosition">
        /// </param>
        /// <param name="SurfacePoint">
        /// </param>
        /// <param name="RadialVector">
        /// </param>
        /// <param name="RadialLength">
        /// </param>
        /// <param name="RadialNormalVector">
        /// Vector normal to the radial vector.
        /// </param>
        internal void CalculateRadialVector(double[] ParticlePosition, double[] SurfacePoint, out double[] RadialVector, out double RadialLength, out double[] RadialNormalVector)
        {
            RadialVector = new double[ParticlePosition.Length];
            for (int d = 0; d < ParticlePosition.Length; d++)
            {
                RadialVector[d] = -SurfacePoint[d] + ParticlePosition[d];
            }
            RadialVector.ScaleV(1 / Math.Sqrt(RadialVector[0].Pow2() + RadialVector[1].Pow2()));
            RadialNormalVector = new double[] { -RadialVector[1], RadialVector[0] };
            RadialLength = Math.Sqrt(RadialNormalVector[0].Pow2() + RadialNormalVector[1].Pow2());
            RadialNormalVector.ScaleV(1 / Math.Sqrt(RadialNormalVector[0].Pow2() + RadialNormalVector[1].Pow2()));
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

        /// <summary>
        /// Predicts the particle state at the next timestep.
        /// </summary>
        /// <param name="particle">
        /// </param>
        /// <param name="SpatialDim">
        /// </param>
        /// <param name="dt">
        /// The timestep.
        /// </param>
        /// <param name="Position">
        /// New position of the particle.
        /// </param>
        /// <param name="TranslationalVelocity">
        /// New velocity of the particle.
        /// </param>
        /// <param name="Angle">
        /// New angle of the particle.
        /// </param>
        /// <param name="RotationalVelocity">
        /// New rotational velocity of the particle.
        /// </param>
        internal void PredictParticleNextTimestep(Particle particle, int SpatialDim, double dt, out double[] Position, out double[] TranslationalVelocity, out double Angle, out double RotationalVelocity)
        {
            Position = new double[SpatialDim];
            TranslationalVelocity = new double[SpatialDim];
            Angle = particle.Angle[0] + particle.RotationalVelocity[0] * dt + (particle.RotationalAcceleration[1] + particle.RotationalAcceleration[0]) * dt.Pow2() / 4;
            RotationalVelocity = particle.RotationalVelocity[0] + (particle.RotationalAcceleration[1] + particle.RotationalAcceleration[0]) * dt / 2;

            for (int d = 0; d < SpatialDim; d++)
            {
                Position[d] = particle.Position[0][d] + particle.TranslationalVelocity[0][d] * dt + (particle.TranslationalAcceleration[1][d] + particle.TranslationalAcceleration[0][d]) * dt.Pow2() / 4;
                TranslationalVelocity[d] = particle.TranslationalVelocity[0][d] + (particle.TranslationalAcceleration[1][d] + particle.TranslationalAcceleration[0][d]) * dt / 2;
            }
        }
    }
}

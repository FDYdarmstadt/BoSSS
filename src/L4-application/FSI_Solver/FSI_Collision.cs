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

        /// ====================================================================================
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
        /// ====================================================================================
        internal void GJK_DistanceAlgorithm(Particle Particle0, Particle Particle1, LevelSetTracker lsTrk, double[] Position0, double[] Position1, double Angle0, double Angle1, out double Min_Distance, out double[] DistanceVec, out double[] ClosestPoint0, out double[] ClosestPoint1, out bool Overlapping)
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
                double[] vt = Aux.VectorDiff(null, v);
                if (double.IsNaN(vt[0]) || double.IsNaN(vt[1]))
                    throw new ArithmeticException("Error trying to calculate point0 Value:  " + vt[0] + " point1 " + vt[1]);

                // =======================================================
                // Step 1
                // Calculate the support point of the minkowski difference
                // and of the two particles (which are the closest points
                // if the algorithm is finished.
                // =======================================================
                CalculateSupportPoint(Particle0, Position0, Angle0, vt, lsTrk, out ClosestPoint0);
                if (double.IsNaN(ClosestPoint0[0]) || double.IsNaN(ClosestPoint0[1]))
                    throw new ArithmeticException("Error trying to calculate point0 Value:  " + ClosestPoint0[0] + " point1 " + ClosestPoint0[1]);
                if (Particle1 != null)
                    CalculateSupportPoint(Particle1, Position1, Angle1, v, lsTrk, out ClosestPoint1);
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
                if ((Aux.DotProduct(v,vt) - Aux.DotProduct(SupportPoint,vt)) >= -1e-12)
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

        /// ====================================================================================
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
        /// ====================================================================================
        private void Initialize_GJK(double[] Position0, double[] Position1, out double[] v0, out List<double[]> Simplex)
        {
            Simplex = new List<double[]>();
            v0 = Aux.VectorDiff(Position0, Position1);
            Simplex.Add(v0.CloneAs());
        }

        /// ====================================================================================
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
        /// ====================================================================================
        private void CalculateSupportPoint(Particle _Particle, double[] Position, double Angle, double[] Vector, LevelSetTracker lsTrk, out double[] SupportPoint)
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
                int L = 1;
                int R = SurfacePoints.GetLength(0) - 2;
                int Counter = 0;
                while (L <= R && L > 0 && R < SurfacePoints.GetLength(0) - 1)
                {
                    int Index = (L + R) / 2;
                    Counter = Counter + 1;
                    GetPointAndNeighbours(SurfacePoints, Index, out SupportPoint, out double[] RightNeighbour, out double[] LeftNeighbour);
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

        /// ====================================================================================
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
        /// ====================================================================================
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

        /// ====================================================================================
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
        /// ====================================================================================
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
                        if (v[0] == 0 && v[1] == 0)
                            Console.WriteLine("Stupid");
                        if (double.IsNaN(v[0]) || double.IsNaN(v[1]))
                            Console.WriteLine("Stupid");
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

        internal void ProjectVelocity(double[] NormalVector, double[] TangentialVector, double[] TranslationalVelocity, out double collisionVn, out double collisionVt)
        {
            ProjectVelocityOnVector(NormalVector, TranslationalVelocity, out collisionVn);
            ProjectVelocityOnVector(TangentialVector, TranslationalVelocity, out collisionVt);
        }

        internal void GetParticleState(Particle particle, int SpatialDim, out double[] ParticleState)
        {
            ParticleState = new double[2 * SpatialDim + 2];
            for (int d = 0; d < SpatialDim; d++)
            {
                ParticleState[d] = particle.Position[0][d];
                ParticleState[d + SpatialDim] = particle.TranslationalVelocity[0][d];
            }
            ParticleState[2 * SpatialDim] = particle.Angle[0];
            ParticleState[2*SpatialDim +1] = particle.RotationalVelocity[0];
        }

        internal void CalculateDynamicCollisionThreshold(Particle particle0, Particle particle1, double[] tempPoint_P0, double[] tempPoint_P1, double[] NormalVector, double Distance, double dt, out double Threshold)
        {
            Threshold = 0;
            FindRadialVector(particle0.Position[0], tempPoint_P0, out _, out double RadialLength0, out double[] RadialNormalVector0);
            TransformRotationalVelocity(particle0.RotationalVelocity[0], RadialLength0, RadialNormalVector0, out double[] PointVelocityDueToRotation0);
            double[] PointVelocity0 = new double[2];
            for (int d = 0; d < 2; d++)
                PointVelocity0[d] = particle0.TranslationalVelocity[0][d] + PointVelocityDueToRotation0[d];
            ProjectVelocityOnVector(NormalVector, PointVelocity0, out double DetectCollisionVn_P0);
            if (particle1 != null)
            {
                FindRadialVector(particle1.Position[0], tempPoint_P1, out _, out double RadialLength1, out double[] RadialNormalVector1);
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

        internal void FindNormalAndTangentialVector(double[] distanceVec, out double[] normal, out double[] tangential)
        {
            normal = distanceVec.CloneAs();
            normal.ScaleV(1 / Math.Sqrt(distanceVec[0].Pow2() + distanceVec[1].Pow2()));
            tangential = new double[] { -normal[1], normal[0] };
        }

        internal void ProjectVelocityOnVector(double[] vector, double[] VelocityVector, out double VelocityComponent)
        {
            VelocityComponent = VelocityVector[0] * vector[0] + VelocityVector[1] * vector[1];
        }

        internal void FindRadialVector(double[] ParticlePosition, double[] SurfacePoint, out double[] RadialVector, out double RadialLength, out double[] RadialNormalVector)
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

        internal void TransformRotationalVelocity(double RotationalVelocity, double RadialLength, double[] RadialNormalVector, out double[] PointVelocityDueToRotation)
        {
            PointVelocityDueToRotation = new double[RadialNormalVector.Length];
            for (int d = 0; d < RadialNormalVector.Length; d++)
            {
                PointVelocityDueToRotation[d] = RadialLength * RotationalVelocity * RadialNormalVector[d];
            }
        }

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

        internal void SetParticleToLastTimestep(Particle particle, int SpatialDim, out double[] Position, out double Angle)
        {
            Position = new double[SpatialDim];
            Angle = particle.Angle[1];
            for (int d = 0; d < SpatialDim; d++)
            {
                Position[d] = particle.Position[1][d];
            }
        }
    }
}

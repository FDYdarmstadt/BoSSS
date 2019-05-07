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
            FindRadialVector(particle0.Position[0], tempPoint_P0, out _, out double[] RadialNormalVector0);
            FindRadialVector(particle1.Position[0], tempPoint_P1, out _, out double[] RadialNormalVector1);
            TransformRotationalVelocity(particle0.RotationalVelocity[0], RadialNormalVector0, out double[] PointVelocityDueToRotation0);
            TransformRotationalVelocity(particle1.RotationalVelocity[0], RadialNormalVector1, out double[] PointVelocityDueToRotation1);

            //general definitions of normal and tangential components
            double[] PointVelocity0 = new double[2];
            double[] PointVelocity1 = new double[2];
            for (int d = 0; d < 2; d++)
            {
                PointVelocity0[d] = particle0.TranslationalVelocity[0][d] + PointVelocityDueToRotation0[d];
                PointVelocity1[d] = particle1.TranslationalVelocity[0][d] + PointVelocityDueToRotation1[d];
            }
            ProjectVelocityOnVector(NormalVector, PointVelocity0, out double DetectCollisionVn_P0);
            ProjectVelocityOnVector(NormalVector, PointVelocity1, out double DetectCollisionVn_P1);
            if (Distance <= ((-DetectCollisionVn_P0 + DetectCollisionVn_P1) * dt))
            {
                Threshold = (-DetectCollisionVn_P0 + DetectCollisionVn_P1) * dt;
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

        internal void FindRadialVector(double[] ParticlePosition, double[] SurfacePoint, out double[] RadialVector, out double[] RadialNormalVector)
        {
            RadialVector = new double[ParticlePosition.Length];
            for (int d = 0; d < ParticlePosition.Length; d++)
            {
                RadialVector[d] = SurfacePoint[d] - ParticlePosition[d];
            }
            RadialNormalVector = new double[] { -RadialVector[1], RadialVector[0] };
            RadialNormalVector.ScaleV(1 / Math.Sqrt(RadialNormalVector[0].Pow2() + RadialNormalVector[1].Pow2()));
        }

        internal void TransformRotationalVelocity(double RotationalVelocity, double[] RadialNormalVector, out double[] PointVelocityDueToRotation)
        {
            PointVelocityDueToRotation = new double[RadialNormalVector.Length];
            for (int d = 0; d < RadialNormalVector.Length; d++)
            {
                PointVelocityDueToRotation[d] = RotationalVelocity * RadialNormalVector[d];
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
            particle.RotationalVelocity[0] = RotationalVelocity;
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

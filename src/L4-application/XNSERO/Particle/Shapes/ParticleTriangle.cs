/* =======================================================================
Copyright 2024 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using MathNet.Numerics.Statistics.Mcmc;
using NUnit.Options;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;
using NUnit.Framework.Constraints;
using System.Runtime.Serialization;


namespace BoSSS.Application.XNSERO_Solver {

    [DataContract]
    [Serializable]
    public class ParticleTriangle : Particle {

        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private ParticleTriangle() : base() {

        }

        /// <summary>
        /// Initializes a triangle particle using a polygon of three points. the points of the triangle needs to be in clockwise order. So far only implemented for noMotion
        /// </summary>
        public ParticleTriangle(IMotion motion, MultidimensionalArray points, bool extrudeZ = false, double[] initVelocity = null, double initRotVelocity = 0.0, double activeStress = 0) 
            : base(motion, activeStress){

            SpatialDim = extrudeZ ? 3 : 2;

            if (points.Lengths[0] != 3)
                throw new ArgumentException("no three points defined");
            if (points.Lengths[1] != 2)
                throw new ArgumentException("special implementation only using 2D points");

            m_points = points;

            // check private methods
            //double dist = DistanceToLineSegment(new double[] { 0.0, 0.0 }, new double[] { 1.0, 0.0 }, new double[] { 1.25, 0.0 });
            //double sign = SignToLineSegment(new double[] { 0.0, 0.0 }, new double[] { 1.0, 0.0 }, new double[] { 1.25, 0.5 });

            //double test1 = ParticleLevelSetFunction(new double[] { -0.5, -0.5 }, new Vector(SpatialDim));
            //double test2 = ParticleLevelSetFunction(new double[] { 0.0, 0.0 }, new Vector(SpatialDim));
            //double test3 = ParticleLevelSetFunction(new double[] { 1.5, -1.0 }, new Vector(SpatialDim));
            //double test4 = ParticleLevelSetFunction(new double[] { 0.25, -0.75 }, new Vector(SpatialDim));
            //double test5 = ParticleLevelSetFunction(new double[] { 0.25, 0.5 }, new Vector(SpatialDim));

            InitializeMotionComponent(motion, GetCentroid(), 0.0, initVelocity == null ? new double[] { 0.0, 0.0} : initVelocity, initRotVelocity);

            //Aux = new Auxillary();

            Circumference = GetCircumference();
            Volume = GetArea();
            Mass = Volume * Density;
            MomentOfInertia = ComputeMomentOfInertia(); // currently just 0.0

        }


        [DataMember]
        MultidimensionalArray m_points;


        private double[] GetPoint(int pInd) {
            if (pInd < 0 || pInd > m_points.Lengths[0])
                throw new ArgumentOutOfRangeException();

            return m_points.ExtractSubArrayShallow(pInd, -1).To1DArray();
        }


        private (double[] pA, double[] pB) GetSegment(int pInd) {
            if (pInd < 0 || pInd > m_points.Lengths[0])
                throw new ArgumentOutOfRangeException();

            double[] pA = GetPoint(pInd);
            double[] pB = pInd == m_points.Lengths[0] - 1 ? GetPoint(0) : GetPoint(pInd + 1);

            return (pA, pB);
        }


        private double[] GetCentroid() {

            //int D = points[0].Length;
            //foreach (double[] p in points) {
            //    if (p.Length != D)
            //        throw new ArgumentException("different spatial dimensions for given points");
            //}
            double[] centroid2D = new double[2];

            for (int p = 0; p < m_points.Lengths[0]; p++) {
                centroid2D.AccV(1.0, GetPoint(p));
            }
            centroid2D.ScaleV(1.0 / m_points.Lengths[0]);

            return SpatialDim == 2 ? centroid2D : new double[] { centroid2D[0], centroid2D[1], 0.0 };
        }


        private double GetCircumference() { 
        
            double circ = 0.0;

            for (int p = 0; p < m_points.Lengths[0] - 1; p++) {
                    circ += GetPoint(p).L2Dist(GetPoint(p+1));
            }

            return circ;

        }


        private double ComputeMomentOfInertia() { 
        
            double mom = 0.0;

            return mom;
        }


        /// <summary>
        /// for triangles by Herons Formula
        /// </summary>
        /// <returns></returns>
        private double GetArea() { 
           
            double s = 0.5 * GetCircumference();
            double sMult = 0.0;
            for (int p = 0; p < m_points.Lengths[0]; p++) {
                (double[] pA, double[] pB) = GetSegment(p);
                double side = pA.L2Dist(pB);
                sMult *= (s - side);
            }

            double area = Math.Sqrt(s * sMult);

            return area;
        }



        //public override double Circumference => GetCircumference();

        //public override double MomentOfInertia => ComputeMomentOfInertia();

        //public override double Volume => GetArea();


        protected override bool ParticleContains(Vector point, Vector Position, double tolerance = 0.0) {

            double[] point2D = new double[] { point.x, point.y };   // 3D only for extrusion in z-direction

            for (int p = 0; p < m_points.Lengths[0]; p++) {
                (double[] pA, double[] pB) = GetSegment(p);
                (double[] proj, bool isWithin) = ProjectionOnSegment(pA, pB, point2D);
                double distSegment = point2D.L2Dist(proj);
                double signSegment = (double)SignToLine(pA, pB, point2D);
                if (signSegment < 0.0 && distSegment > tolerance)
                    return false;
            }

            return true;
        }


        protected override double ParticleLevelSetFunction(double[] X, Vector Postion) {

            double[] X2D = new double[] { X[0], X[1] };   // 3D only for extrusion in z-direction

            double dist = double.MaxValue;
            double sign = -1.0;
            for (int p = 0; p < m_points.Lengths[0] - 1; p++) {
                (double[] pA, double[] pB) = GetSegment(p);
                (double[] proj, bool isWithin) = ProjectionOnSegment(pA, pB, X2D);
                double distSegment = X2D.L2Dist(proj);
                double signSegment = (double)SignToLine(pA, pB, X2D);
                if (distSegment < dist) {
                    dist = distSegment;
                    sign = isWithin ? signSegment : -1.0;
                }
            }

            return sign * dist;
        }


        /// <summary>
        /// returns the projection of point P onto the segment between point A and B
        /// </summary>
        /// <param name="pA"></param>
        /// <param name="pB"></param>
        /// <param name="pP"></param>
        /// <returns></returns>
        private static (double[] proj, bool isWithin) ProjectionOnSegment(double[] pA, double[] pB, double[] pP) {

            double l2 = pA.L2DistPow2(pB);

            //const float t = max(0, min(1, dot(p - v, w - v) / l2));
            //const vec2 projection = v + t * (w - v);  // Projection falls on the segment
            //return distance(p, projection);

            double[] AP = pP.CloneAs();
            AP.AccV(-1.0, pA);
            double[] AB = pB.CloneAs();
            AB.AccV(-1.0, pA);

            double dotP = AP.InnerProd(AB);
            double t = dotP / l2;
            double t01 = Math.Max(0.0, Math.Min(1.0, t));

            double[] proj = pA.CloneAs();
            proj.AccV(t01, AB);

            return (proj, t >= 0.0 && t <= 1.0);
        }


        //private static double DistanceToLineSegment(double[] pA, double[] pB, double[] pP) {

        //    (double[] proj, bool isWithin) = ProjectionOnSegment(pA, pB, pP);

        //    return pP.L2Dist(proj);
        //}


        /// <summary>
        /// returns on which side on the segment AB lies the point P.
        /// -1 left of the extended segment AB and +1 right or on the extended segment AB
        /// </summary>
        /// <returns></returns>
        private static int SignToLine(double[] pA, double[] pB, double[] pP) {

            //position = sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax))

            double det = (pB[0] - pA[0]) * (pP[1] - pA[1]) - (pB[1] - pA[1]) * (pP[0] - pA[0]);

            return Math.Sign(det) == 0 ? -1 : -Math.Sign(det);

        }


        /// <summary>
        /// Returns surface points (for distance calc)
        /// </summary>
        /// <param name="hMin"></param>
        /// <param name="searchAngle"></param>
        /// <param name="subParticleID">between 0 and <see cref="NoOfSubParticles"/>, i guess</param>
        public override MultidimensionalArray GetSurfacePoints(double hMin, double searchAngle, int subParticleID) {
            // so far not needed within no motion calculations
            return MultidimensionalArray.Create(3, 2);
        }


        public override double[] GetLengthScales() {
            // returning side lengths
            double[] lsides = new double[3];
            for (int p = 0; p < m_points.Lengths[0] - 1; p++) {
                lsides[0] = GetPoint(p).L2Dist(GetPoint(p + 1));
            }
            return lsides;
        }
    }
}

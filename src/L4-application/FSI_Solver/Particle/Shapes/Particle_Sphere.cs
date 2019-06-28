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

using System;
using System.Runtime.Serialization;
using BoSSS.Foundation.XDG;
using ilPSP;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using System.Linq;
using System.Collections.Generic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Grid.Classic;
using System.Text;
using BoSSS.Platform;
using ilPSP.Tracing;
using System.Diagnostics;
using ilPSP.Utils;
using System.IO;
using System.Collections;
using BoSSS.Platform.Utils.Geom;

namespace BoSSS.Application.FSI_Solver
{
    [DataContract]
    [Serializable]
    public class Particle_Sphere : Particle
    {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Sphere() : base()
        {

        }

        public Particle_Sphere(double[] startPos = null, double startAngl = 0) : base(2, startPos, startAngl) {


        }

        /// <summary>
        /// Radius of the particle. Not necessary for particles defined by their length and thickness
        /// </summary>
        [DataMember]
        public double radius_P;

        /// <summary>
        /// %
        /// </summary>
        protected override double AverageDistance {
            get {
                return radius_P;
            }
        }

        public override double Area_P
        {
            get
            {
                return Math.PI * radius_P * radius_P;
            }
        }
        protected override double Circumference_P
        {
            get
            {
                return 2 * Math.PI * radius_P;
            }
        }
        override public double MomentOfInertia_P
        {
            get
            {
                return (1 / 2.0) * (Mass_P * radius_P * radius_P);
            }
        }
        //override public void UpdateLevelSetFunction()
        //{
        //    double alpha = -(Angle[0]);
        //    Phi_P = (X, t) => -(X[0] - Position[0][0]).Pow2() + -(X[1] - Position[0][1]).Pow2() + radius_P.Pow2();
        //}
        public override double Phi_P(double[] X) {
            double x0 = Position[0][0];
            double y0 = Position[0][1];
            return -(X[0] - x0).Pow2() + -(X[1] - y0).Pow2() + radius_P.Pow2();
        }

        public override bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false)
        {
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double radiusTolerance = !WithoutTolerance ? radius_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : radius_P;
            var distance = point.L2Distance(Position[0]);
            if (distance < (radiusTolerance))
            {
                return true;
            }
            return false;
        }

        override public MultidimensionalArray GetSurfacePoints(LevelSetTracker lsTrk, double[] Position, double Angle)
        {
            int SpatialDim = lsTrk.GridDat.SpatialDimension;
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            double hMin = lsTrk.GridDat.iGeomCells.h_min.Min();
            int NoOfSurfacePoints = Convert.ToInt32(10 * Circumference_P / hMin) + 1;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSubParticles(), NoOfSurfacePoints, SpatialDim);
            double[] InfinitisemalAngle = GenericBlas.Linspace(0, 2 * Math.PI, NoOfSurfacePoints + 1);
            if (Math.Abs(10 * Circumference_P / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");
            
            for (int j = 0; j < NoOfSurfacePoints; j++)
            {
                SurfacePoints[0, j, 0] = Math.Cos(InfinitisemalAngle[j]) * radius_P + Position[0];
                SurfacePoints[0, j, 1] = Math.Sin(InfinitisemalAngle[j]) * radius_P + Position[1];
            }
            return SurfacePoints;
        }

        override public void GetSupportPoint(int SpatialDim, double[] Vector, double[] Position, double Angle, out double[] SupportPoint)
        {
            double length = Math.Sqrt(Vector[0].Pow2() + Vector[1].Pow2());
            double CosT = Vector[0] / length;
            double SinT = Vector[1] / length;
            SupportPoint = new double[SpatialDim];
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");
            SupportPoint[0] = CosT * radius_P + Position[0];
            SupportPoint[1] = SinT * radius_P + Position[1];
            if (double.IsNaN(SupportPoint[0]) || double.IsNaN(SupportPoint[1]))
                throw new ArithmeticException("Error trying to calculate point0 Value:  " + SupportPoint[0] + " point1 " + SupportPoint[1]);
        }

        override public double[] GetLengthScales()
        {
            return new double[] { radius_P, radius_P };
        }
    }
}


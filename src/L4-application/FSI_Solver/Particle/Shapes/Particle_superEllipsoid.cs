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
using MathNet.Numerics;
using System.Diagnostics;
using ilPSP.Utils;

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_superEllipsoid : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_superEllipsoid() : base() {

        }

        /// <summary>
        /// ctor
        /// </summary>
        public Particle_superEllipsoid(double[] startPos = null, double startAngl = 0) : base(2, startPos, startAngl) {
        }

        /// <summary>
        /// Length of an elliptic particle.
        /// </summary>
        [DataMember]
        public double length_P;

        /// <summary>
        /// Thickness of an elliptic particle.
        /// </summary>
        [DataMember]
        public double thickness_P;

        /// <summary>
        /// Exponent of the super ellipsoid. Higher exponent leads to a more "squary" appearance.
        /// </summary>
        [DataMember]
        public double superEllipsoidExponent;

        protected override double Circumference_P {
            get {
                return (2 * length_P + 2 * thickness_P + 2 * Math.PI * thickness_P) / 2;
            }
        }

        public override double Area_P {
            get {
                return 4 * length_P * thickness_P * (SpecialFunctions.Gamma(1 + 1 / superEllipsoidExponent)).Pow2() / SpecialFunctions.Gamma(1 + 2 / superEllipsoidExponent);
            }

        }
        override public double MomentOfInertia_P {
            get {
                return (1 / 4.0) * Mass_P * (length_P * length_P + thickness_P * thickness_P);
            }
        }
        public override double Phi_P(double[] X) {
            double alpha = -(Angle[0]);
            double r;
            r = -Math.Pow(
                        ((X[0] - Position[0][0]) * Math.Cos(alpha) - (X[1] - Position[0][1]) * Math.Sin(alpha)) / length_P,
                        superEllipsoidExponent)
                - Math.Pow(
                    ((X[0] - Position[0][0]) * Math.Sin(alpha) + (X[1] - Position[0][1]) * Math.Cos(alpha)) / thickness_P,
                    superEllipsoidExponent)
                + 1;
            if (double.IsNaN(r) || double.IsInfinity(r))
                throw new ArithmeticException();
            return r;
        }
        public override bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false)
        {
            WithoutTolerance = false;
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double radiusTolerance = 1;
            double a = !WithoutTolerance ? length_P + h_min : length_P;
            double b = !WithoutTolerance ? thickness_P + h_min : thickness_P;
            double Superellipsoid = Math.Pow(((point[0] - Position[0][0]) * Math.Cos(Angle[0]) + (point[1] - Position[0][1]) * Math.Sin(Angle[0])) / a, superEllipsoidExponent) + (Math.Pow((-(point[0] - Position[0][0]) * Math.Sin(Angle[0]) + (point[1] - Position[0][1]) * Math.Cos(Angle[0])) / b,superEllipsoidExponent));
            if (Superellipsoid < radiusTolerance)
                return true;
            else
                return false;
        }
        override public double[] GetLengthScales()
        {
            return new double[] { length_P, thickness_P };
        }

        override public MultidimensionalArray GetSurfacePoints(LevelSetTracker lsTrk, double[] Position, double Angle)
        {
            int SpatialDim = lsTrk.GridDat.SpatialDimension;
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            double hMin = lsTrk.GridDat.iGeomCells.h_min.Min();
            int NoOfSurfacePoints = Convert.ToInt32(10 * Circumference_P / hMin);
            int QuarterSurfacePoints = NoOfSurfacePoints / 4;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSubParticles(), 4 * QuarterSurfacePoints - 2, SpatialDim);
            double[] InfinitisemalAngle = GenericBlas.Linspace(0, Math.PI / 2, QuarterSurfacePoints + 2);
            if (Math.Abs(10 * Circumference_P / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");
            for (int j = 0; j < QuarterSurfacePoints; j++)
            {
                SurfacePoints[0, j, 0] = (Math.Pow(Math.Cos(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * length_P * Math.Cos(Angle) - Math.Pow(Math.Sin(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * thickness_P * Math.Sin(Angle)) + Position[0];
                SurfacePoints[0, j, 1] = (Math.Pow(Math.Cos(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * length_P * Math.Sin(Angle) + Math.Pow(Math.Sin(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * thickness_P * Math.Cos(Angle)) + Position[1]; 
                SurfacePoints[0, 2 * QuarterSurfacePoints + j - 1, 0] = (-(Math.Pow(Math.Cos(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * length_P) * Math.Cos(Angle) + Math.Pow(Math.Sin(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * thickness_P * Math.Sin(Angle)) + Position[0];
                SurfacePoints[0, 2 * QuarterSurfacePoints + j - 1, 1] = (-(Math.Pow(Math.Cos(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * length_P) * Math.Sin(Angle) - Math.Pow(Math.Sin(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * thickness_P * Math.Cos(Angle)) + Position[1];;
            }
            for (int j = 1; j < QuarterSurfacePoints; j++)
            {
                SurfacePoints[0, 2 * QuarterSurfacePoints - j - 1, 0] = (-(Math.Pow(Math.Cos(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * length_P) * Math.Cos(Angle) - Math.Pow(Math.Sin(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * thickness_P * Math.Sin(Angle)) + Position[0];
                SurfacePoints[0, 2 * QuarterSurfacePoints - j - 1, 1] = (-(Math.Pow(Math.Cos(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * length_P) * Math.Sin(Angle) + Math.Pow(Math.Sin(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * thickness_P * Math.Cos(Angle)) + Position[1];
                SurfacePoints[0, 4 * QuarterSurfacePoints - j - 2, 0] = (Math.Pow(Math.Cos(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * length_P * Math.Cos(Angle) + Math.Pow(Math.Sin(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * thickness_P * Math.Sin(Angle)) + Position[0];
                SurfacePoints[0, 4 * QuarterSurfacePoints - j - 2, 1] = (Math.Pow(Math.Cos(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * length_P * Math.Sin(Angle) - Math.Pow(Math.Sin(InfinitisemalAngle[j]), 2 / superEllipsoidExponent) * thickness_P * Math.Cos(Angle)) + Position[1];
            }
            return SurfacePoints;
        }

        //override public void GetSupportPoint(int SpatialDim, double[] Vector, out double[] SupportPoint)
        //{
        //    SupportPoint = new double[SpatialDim];
        //    double VectorLength = Math.Sqrt(Vector[0].Pow2() + Vector[1].Pow2());
        //    for (int d = 0; d < SpatialDim; d++)
        //    {
        //        Vector[d] = Vector[d] / VectorLength;
        //    }
        //    double[] VectorBodyCOS = Vector.CloneAs();
        //    VectorBodyCOS[0] = Vector[0] * Math.Cos(Angle[0]) - Vector[1] * Math.Sin(Angle[0]);
        //    VectorBodyCOS[1] = Vector[0] * Math.Sin(Angle[0]) - Vector[1] * Math.Cos(Angle[1]);
        //    double VecBodyLength = Math.Pow(1 / (Math.Pow(VectorBodyCOS[0] / length_P, superEllipsoidExponent) + Math.Pow(VectorBodyCOS[1] / thickness_P, superEllipsoidExponent)), 1 / superEllipsoidExponent);
        //    double[] StartPoint = new double[SpatialDim];
        //    for (int d = 0; d < SpatialDim; d++)
        //    {
        //        StartPoint[d] = Position[0][d] + VecBodyLength * Vector[d];
        //    }
            
        //}
    }
}


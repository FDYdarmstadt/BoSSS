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
using ilPSP;
using System.Linq;
using ilPSP.Utils;
using MathNet.Numerics;

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
        public Particle_superEllipsoid(ParticleMotionInit motionInit, double length, double thickness, int superEllipsoidExponent, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            length_P = length;
            thickness_P = thickness;
            m_SuperEllipsoidExponent = superEllipsoidExponent;
            Motion.GetParticleLengthscale(GetLengthScales().Max());
            Motion.GetParticleArea(Area_P());
            Motion.GetParticleMomentOfInertia(MomentOfInertia_P);
        }

        /// <summary>
        /// Length of an elliptic particle.
        /// </summary>
        private readonly double length_P;

        /// <summary>
        /// Thickness of an elliptic particle.
        /// </summary>
        private readonly double thickness_P;

        /// <summary>
        /// Exponent of the super ellipsoid. Higher exponent leads to a more "squary" appearance.
        /// </summary>
        [DataMember]
        private readonly double m_SuperEllipsoidExponent;

        protected override double Circumference_P {
            get {
                return (2 * length_P + 2 * thickness_P + 2 * Math.PI * thickness_P) / 2;
            }
        }
        public override double Area_P() {
            return 4 * length_P * thickness_P * (SpecialFunctions.Gamma(1 + 1 / m_SuperEllipsoidExponent)).Pow2() / SpecialFunctions.Gamma(1 + 2 / m_SuperEllipsoidExponent);
        }

        override public double MomentOfInertia_P {
            get {
                return (1 / 4.0) * Mass_P * (length_P * length_P + thickness_P * thickness_P);
            }
        }

        public override double LevelSetFunction(double[] X) {
            double alpha = -(Motion.GetAngle(0));
            double r;
            r = -Math.Pow(((X[0] - Motion.GetPosition(0)[0]) * Math.Cos(alpha) - (X[1] - Motion.GetPosition(0)[1]) * Math.Sin(alpha)) / length_P, m_SuperEllipsoidExponent)
                - Math.Pow(((X[0] - Motion.GetPosition(0)[0]) * Math.Sin(alpha) + (X[1] - Motion.GetPosition(0)[1]) * Math.Cos(alpha)) / thickness_P, m_SuperEllipsoidExponent)
                + 1;
            if (double.IsNaN(r) || double.IsInfinity(r))
                throw new ArithmeticException();
            return r;
        }

        public override bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false) {
            WithoutTolerance = false;
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double radiusTolerance = 1;
            double a = !WithoutTolerance ? length_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : length_P;
            double b = !WithoutTolerance ? thickness_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : thickness_P;
            double Superellipsoid = Math.Pow(((point[0] - Motion.GetPosition(0)[0]) * Math.Cos(Motion.GetAngle(0)) + (point[1] - Motion.GetPosition(0)[1]) * Math.Sin(Motion.GetAngle(0))) / a, m_SuperEllipsoidExponent) + (Math.Pow((-(point[0] - Motion.GetPosition(0)[0]) * Math.Sin(Motion.GetAngle(0)) + (point[1] - Motion.GetPosition(0)[1]) * Math.Cos(Motion.GetAngle(0))) / b, m_SuperEllipsoidExponent));
            if (Superellipsoid < radiusTolerance)
                return true;
            else
                return false;
        }

        override public double[] GetLengthScales() {
            return new double[] { length_P, thickness_P };
        }

        override public MultidimensionalArray GetSurfacePoints(double hMin) {
            if (spatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            int NoOfSurfacePoints = Convert.ToInt32(10 * Circumference_P / hMin);
            int QuarterSurfacePoints = NoOfSurfacePoints / 4;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSubParticles, 4 * QuarterSurfacePoints - 2, spatialDim);
            double[] Infinitisemalangle = GenericBlas.Linspace(0, Math.PI / 2, QuarterSurfacePoints + 2);
            if (Math.Abs(10 * Circumference_P / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");
            for (int j = 0; j < QuarterSurfacePoints; j++) {
                SurfacePoints[0, j, 0] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * length_P * Math.Cos(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * thickness_P * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
                SurfacePoints[0, j, 1] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * length_P * Math.Sin(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * thickness_P * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
                SurfacePoints[0, 2 * QuarterSurfacePoints + j - 1, 0] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * length_P) * Math.Cos(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * thickness_P * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
                SurfacePoints[0, 2 * QuarterSurfacePoints + j - 1, 1] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * length_P) * Math.Sin(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * thickness_P * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1]; ;
            }
            for (int j = 1; j < QuarterSurfacePoints; j++) {
                SurfacePoints[0, 2 * QuarterSurfacePoints - j - 1, 0] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * length_P) * Math.Cos(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * thickness_P * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
                SurfacePoints[0, 2 * QuarterSurfacePoints - j - 1, 1] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * length_P) * Math.Sin(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * thickness_P * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
                SurfacePoints[0, 4 * QuarterSurfacePoints - j - 2, 0] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * length_P * Math.Cos(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * thickness_P * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
                SurfacePoints[0, 4 * QuarterSurfacePoints - j - 2, 1] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * length_P * Math.Sin(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_SuperEllipsoidExponent) * thickness_P * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
            }
            return SurfacePoints;
        }
    }
}


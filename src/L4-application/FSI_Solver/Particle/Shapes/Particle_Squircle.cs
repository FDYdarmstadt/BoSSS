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
using ilPSP.Utils;

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_Squircle : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Squircle() : base() {

        }

        /// <summary>
        /// Constructor for a squircle (symmetric superellipsoid).
        /// </summary>
        /// <param name="motionInit">
        /// Initializes the motion parameters of the particle (which model to use, whether it is a dry simulation etc.)
        /// </param>
        /// <param name="radius">
        /// The radius.
        /// </param>
        /// <param name="startPos">
        /// The initial position.
        /// </param>
        /// <param name="startAngl">
        /// The inital anlge.
        /// </param>
        /// <param name="activeStress">
        /// The active stress excerted on the fluid by the particle. Zero for passive particles.
        /// </param>
        /// <param name="startTransVelocity">
        /// The inital translational velocity.
        /// </param>
        /// <param name="startRotVelocity">
        /// The inital rotational velocity.
        /// </param>
        public Particle_Squircle(ParticleMotionInit motionInit, double radius, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            m_Radius = radius;
            Aux.TestArithmeticException(radius, "Particle radius");

            Motion.GetParticleLengthscale(radius);
            Motion.GetParticleArea(Area);
            Motion.GetParticleMomentOfInertia(MomentOfInertia);
        }

        [DataMember]
        private readonly double m_Radius;
        [DataMember]
        private int m_Exponent;

        /// <summary>
        /// Circumference. 
        /// </summary>
        protected override double Circumference => 4 * m_Radius;

        /// <summary>
        /// Area occupied by the particle.
        /// </summary>
        public override double Area => 4 * m_Radius.Pow2() * (SpecialFunctions.Gamma(1 + 1 / m_Exponent)).Pow2() / SpecialFunctions.Gamma(1 + 2 / m_Exponent);

        /// <summary>
        /// Moment of inertia. 
        /// </summary>
        override public double MomentOfInertia => (1 / 2.0) * (Mass_P * m_Radius * m_Radius);

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        public override double LevelSetFunction(double[] X) {
            double alpha = -(Motion.GetAngle(0));
            return -((((X[0] - Motion.GetPosition(0)[0]) * Math.Cos(alpha) - (X[1] - Motion.GetPosition(0)[1]) * Math.Sin(alpha)).Pow(4) + ((X[0] - Motion.GetPosition(0)[0]) * Math.Sin(alpha) + (X[1] - Motion.GetPosition(0)[1]) * Math.Cos(alpha)).Pow(4)) - m_Radius.Pow(4));
        }

        /// <summary>
        /// Returns true if a point is withing the particle.
        /// </summary>
        /// <param name="point">
        /// The point to be tested.
        /// </param>
        /// <param name="minTolerance">
        /// Minimum tolerance length.
        /// </param>
        /// <param name="maxTolerance">
        /// Maximal tolerance length. Equal to h_min if not specified.
        /// </param>
        /// <param name="WithoutTolerance">
        /// No tolerance.
        /// </param>
        public override bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false) {
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double radiusTolerance = !WithoutTolerance ? 1.0 + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : 1;
            if (-((((point[0] - Motion.GetPosition(0)[0]) * Math.Cos(Motion.GetAngle(0)) - (point[1] - Motion.GetPosition(0)[1]) * Math.Sin(Motion.GetAngle(0))).Pow(4) + ((point[0] - Motion.GetPosition(0)[0]) * Math.Sin(Motion.GetAngle(0)) + (point[1] - Motion.GetPosition(0)[1]) * Math.Cos(Motion.GetAngle(0))).Pow(4)) - radiusTolerance.Pow(4)) > 0) {
                return true;
            }
            return false;
        }

        /// <summary>
        /// Returns an array with points on the surface of the particle.
        /// </summary>
        /// <param name="hMin">
        /// Minimal cell length. Used to specify the number of surface points.
        /// </param>
        override public MultidimensionalArray GetSurfacePoints(double hMin) {
            if (spatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            int NoOfSurfacePoints = Convert.ToInt32(10 * Circumference / hMin);
            int QuarterSurfacePoints = NoOfSurfacePoints / 4;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSubParticles, 4 * QuarterSurfacePoints - 2, spatialDim);
            double[] Infinitisemalangle = GenericBlas.Linspace(0, Math.PI / 2, QuarterSurfacePoints + 2);
            if (Math.Abs(10 * Circumference / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");
            for (int j = 0; j < QuarterSurfacePoints; j++) {
                SurfacePoints[0, j, 0] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Cos(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
                SurfacePoints[0, j, 1] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Sin(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
                SurfacePoints[0, 2 * QuarterSurfacePoints + j - 1, 0] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius) * Math.Cos(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
                SurfacePoints[0, 2 * QuarterSurfacePoints + j - 1, 1] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius) * Math.Sin(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1]; ;
            }
            for (int j = 1; j < QuarterSurfacePoints; j++) {
                SurfacePoints[0, 2 * QuarterSurfacePoints - j - 1, 0] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius) * Math.Cos(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
                SurfacePoints[0, 2 * QuarterSurfacePoints - j - 1, 1] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius) * Math.Sin(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
                SurfacePoints[0, 4 * QuarterSurfacePoints - j - 2, 0] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Cos(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
                SurfacePoints[0, 4 * QuarterSurfacePoints - j - 2, 1] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Sin(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Radius * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
            }
            return SurfacePoints;
        }

        /// <summary>
        /// Returns the legnthscales of a particle.
        /// </summary>
        override public double[] GetLengthScales() {
            return new double[] { m_Radius, m_Radius };
        }
    }
}


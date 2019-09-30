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
using ilPSP.Utils;

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_Pentagone : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Pentagone() : base() {

        }

        /// <summary>
        /// Constructor for a pentagone.
        /// </summary>
        /// <param name="motionInit">
        /// Initializes the motion parameters of the particle (which model to use, whether it is a dry simulation etc.)
        /// </param>
        /// <param name="radius">
        /// The main lengthscale
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
        public Particle_Pentagone(ParticleMotionInit motionInit, double radius, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            m_Radius = radius;
            Aux.TestArithmeticException(radius, "Particle radius");

            Motion.GetParticleLengthscale(radius);
            Motion.GetParticleArea(Area);
            Motion.GetParticleMomentOfInertia(MomentOfInertia);

        }

        [DataMember]
        private readonly double m_Radius;

        /// <summary>
        /// Area occupied by the particle.
        /// </summary>
        public override double Area => (5 * m_Radius.Pow2()) / 4;

        /// <summary>
        /// Circumference. 
        /// </summary>
        protected override double Circumference => m_Radius * 4.41421356;

        /// <summary>
        /// Moment of inertia. 
        /// </summary>
        override public double MomentOfInertia => Math.Pow(m_Radius, 4) * 0.2887963;

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        public override double LevelSetFunction(double[] X) { // attention: no dependency on angle...
            double r;
            r = Math.Max(X[0] - Motion.GetPosition(0)[0] - m_Radius, Motion.GetPosition(0)[0] - m_Radius - X[0]);
            r = Math.Max(r, Motion.GetPosition(0)[1] - 0.5 * m_Radius - X[1]);
            r = Math.Max(r, Motion.GetPosition(0)[0] - m_Radius - X[1] - 1.5 * X[0]) + Math.Max(r, X[1] - Motion.GetPosition(0)[1] - 0.5 * m_Radius);
            r = -r;
            return r;
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
            double radiusTolerance = !WithoutTolerance ? m_Radius + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : 1;
            var distance = point.L2Distance(Motion.GetPosition(0));
            if (distance < (radiusTolerance)) {
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
        public override MultidimensionalArray GetSurfacePoints(double hMin) {
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            int NoOfSurfacePoints = Convert.ToInt32(20 * Circumference / hMin) + 1;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSurfacePoints, SpatialDim);
            double[] InfinitisemalAngle = GenericBlas.Linspace(0, 2 * Math.PI, NoOfSurfacePoints + 1);
            if (Math.Abs(10 * Circumference / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");

            for (int j = 0; j < NoOfSurfacePoints; j++) {
                SurfacePoints[j, 0] = Math.Cos(InfinitisemalAngle[j]) * m_Radius + Motion.GetPosition(0)[0];
                SurfacePoints[j, 1] = Math.Sin(InfinitisemalAngle[j]) * m_Radius + Motion.GetPosition(0)[1];
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

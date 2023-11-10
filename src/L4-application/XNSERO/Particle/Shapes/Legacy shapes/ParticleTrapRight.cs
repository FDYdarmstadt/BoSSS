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

namespace BoSSS.Application.XNSERO_Solver {
    [DataContract]
    [Serializable]
    public class ParticleTrapRight : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private ParticleTrapRight() : base() {

        }

        /// <summary>
        /// Constructor for the trap used in the masters thesis if E. Deriabina (2019)
        /// </summary>
        /// <param name="motionInit">
        /// Initializes the motion parameters of the particle (which model to use, whether it is a dry simulation etc.)
        /// </param>
        /// <param name="width">
        /// The main lengthscale.
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
        public ParticleTrapRight(IMotion motion, double width, double[] startPos, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motion, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            throw new NotImplementedException("Legacy code, untested, update necessary");
            if (startPos.Length != 2)
                throw new ArgumentOutOfRangeException("Spatial dimension does not fit particle definition");

            m_Length = width;
            Aux.TestArithmeticException(width, "Particle width");

            Motion.CharacteristicLength = width;
            Motion.Volume = this.Volume;
            Motion.MomentOfInertia = this.MomentOfInertia;

        }

        [DataMember]
        private readonly double m_Length;

        /// <summary>
        /// The trap is devided into two convex sub particles. Necesarry for the GJK-algorithm in the collision model.
        /// </summary>
        public override int NoOfSubParticles => 2;

        /// <summary>
        /// Area occupied by the particle. 
        /// </summary>
        public override double Volume => (7 * m_Length * m_Length) / 8;

        /// <summary>
        /// Circumference.
        /// </summary>
        public override double Circumference => m_Length * 5;

        /// <summary>
        /// Moment of inertia.
        /// </summary>
        override public double MomentOfInertia => Math.Pow(m_Length, 4) * 0.301627768;

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        protected override double ParticleLevelSetFunction(double[] X, Vector Postion) {
            double r;
            // Falle_Rechts:
            r = Math.Abs(Postion[1] - X[1]);
            r = Math.Max(r, Math.Abs(-X[1] - 0.5 * X[0] + Postion[1] + Postion[0] - m_Length) - Math.Abs(X[1] - Postion[1]));
            r = Math.Max(r, Math.Abs(Postion[0] - X[0] + 0.5 * m_Length));
            r -= 4.5 * m_Length;
            r = -r;
            return r;
        }

        /// <summary>
        /// Returns true if a point is withing the particle.
        /// </summary>
        /// <param name="point">
        /// The point to be tested.
        /// </param>
        /// <param name="tolerance">
        /// tolerance length.
        /// </param>
        protected override bool ParticleContains(Vector point, Vector Position, double tolerance = 0) {
            // only for rectangular cells
            double radiusTolerance = 5 * m_Length + tolerance;
            double distance = point.L2Distance(Position);
            return distance < radiusTolerance;
        }

        /// <summary>
        /// Returns an array with points on the surface of the particle.
        /// </summary>
        /// <param name="hMin">
        /// Minimal cell length. Used to specify the number of surface points.
        /// </param>
        override public MultidimensionalArray GetSurfacePoints(double hMin, double searchAngle, int subParticleID) {
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            int NoOfSurfacePoints = Convert.ToInt32(20 * Circumference / hMin) + 1;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSubParticles, NoOfSurfacePoints, SpatialDim);
            double[] InfinitisemalAngle = GenericBlas.Linspace(-Math.PI / 4, 5 * Math.PI / 4, NoOfSurfacePoints + 1);
            double[] InfinitisemalLength = GenericBlas.Linspace(0, m_Length / 4, NoOfSurfacePoints + 1);
            if (Math.Abs(10 * Circumference / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");

            for (int k = 0; k < NoOfSurfacePoints; k++) {
                SurfacePoints[0, k, 0] = Motion.GetPosition(0)[0] + m_Length / 2 - InfinitisemalLength[k];
                SurfacePoints[0, k, 1] = Motion.GetPosition(0)[1] - m_Length / 2 + 1.5 * SurfacePoints[0, k, 0] + m_Length / 2;
            }

            for (int j = 0; j < NoOfSurfacePoints; j++) {
                SurfacePoints[0, j, 0] = Math.Sign(Math.Cos(InfinitisemalAngle[j])) * m_Length * 7 + Motion.GetPosition(0)[0] + 7 * m_Length / 4;
                SurfacePoints[0, j, 1] = Math.Sign(Math.Sin(InfinitisemalAngle[j])) * m_Length * 7 + Motion.GetPosition(0)[1] + 7 * m_Length / 2;
            }
            for (int j = 0; j < NoOfSurfacePoints; j++) {
                SurfacePoints[1, j, 0] = -Math.Sign(Math.Cos(InfinitisemalAngle[j])) * m_Length * 2.5 + Motion.GetPosition(0)[0];
                SurfacePoints[1, j, 1] = -Math.Sign(Math.Sin(InfinitisemalAngle[j])) * m_Length * 2.5 + Motion.GetPosition(0)[1];
            }
            return SurfacePoints;
        }

        /// <summary>
        /// Returns the legnthscales of a particle.
        /// </summary>
        override public double[] GetLengthScales() {
            return new double[] { m_Length, m_Length };
        }
    }
}

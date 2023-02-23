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

namespace BoSSS.Application.XNSERO_Solver {
    [DataContract]
    [Serializable]
    public class ParticleDisk : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private ParticleDisk() : base() {

        }

        /// <summary>
        /// Constructor for a sphere.
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
        public ParticleDisk(IMotion motion, double radius, double[] startPos, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motion, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            if (startPos.Length != 2)
                throw new ArgumentOutOfRangeException("Spatial dimension does not fit particle definition");

            m_Radius = radius;
            Aux.TestArithmeticException(radius, "Particle radius");

            Motion.CharacteristicLength = radius;
            Motion.Volume = this.Volume;
            Motion.MomentOfInertia = this.MomentOfInertia;

        }

        [DataMember]
        private readonly double m_Radius;

        /// <summary>
        /// Area occupied by the particle.
        /// </summary>
        public override double Volume => Math.PI * m_Radius.Pow2();

        /// <summary>
        /// Circumference. 
        /// </summary>
        public override double Circumference => 2 * Math.PI * m_Radius;

        /// <summary>
        /// Moment of inertia. 
        /// </summary>
        override public double MomentOfInertia => (1 / 2.0) * (Mass * m_Radius.Pow2());

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        protected override double ParticleLevelSetFunction(double[] X, Vector Postion) {
            double x0 = Postion[0];
            double y0 = Postion[1];
            return -(X[0] - x0).Pow2() + -(X[1] - y0).Pow2() + m_Radius.Pow2();
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
            double radiusTolerance = m_Radius + tolerance;
            double distance = point.L2Distance(Motion.GetPosition(0));
            return distance < radiusTolerance;
        }

        /// <summary>
        /// Returns the support point of the particle in the direction specified by a vector.
        /// </summary>
        /// <param name="vector">
        /// A vector. 
        /// </param>
        override public Vector GetSupportPoint(Vector supportVector, Vector Position, Vector Angle, int SubParticleID, double tolerance = 0) {
            double length = Math.Sqrt(supportVector[0].Pow2() + supportVector[1].Pow2());
            double CosT = supportVector[0] / length;
            double SinT = supportVector[1] / length;
            Vector SupportPoint = new Vector(SpatialDim);
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");
            SupportPoint[0] = CosT * (m_Radius + tolerance) + Position[0];
            SupportPoint[1] = SinT * (m_Radius + tolerance) + Position[1];
            if (double.IsNaN(SupportPoint[0]) || double.IsNaN(SupportPoint[1]))
                throw new ArithmeticException("Error trying to calculate point0 Value:  " + SupportPoint[0] + " point1 " + SupportPoint[1]);
            return SupportPoint;
        }

        /// <summary>
        /// Returns the legnthscales of a particle.
        /// </summary>
        override public double[] GetLengthScales() {
            return new double[] { m_Radius, m_Radius };
        }

        public override object Clone() {
            Particle clonedParticle = new ParticleDisk(Motion,
                                                             m_Radius,
                                                             Motion.GetPosition(),
                                                             Motion.GetAngle() * 360 / (2 * Math.PI),
                                                             ActiveStress,
                                                             Motion.GetTranslationalVelocity(),
                                                             Motion.GetRotationalVelocity());
            return clonedParticle;
        }
    }
}


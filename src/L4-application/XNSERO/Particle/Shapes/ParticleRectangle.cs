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

using System;
using System.Runtime.Serialization;
using ilPSP;
using System.Linq;

namespace BoSSS.Application.XNSERO_Solver {
    [DataContract]
    [Serializable]
    public class ParticleRectangle : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private ParticleRectangle() : base() {

        }

        /// <summary>
        /// Constructor for an ellipsoid.
        /// </summary>
        /// <param name="motionInit">
        /// Initializes the motion parameters of the particle (which model to use, whether it is a dry simulation etc.)
        /// </param>
        /// <param name="length">
        /// The length of the horizontal halfaxis.
        /// </param>
        /// <param name="thickness">
        /// The length of the vertical halfaxis.
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
        public ParticleRectangle(IMotion motion, double length, double thickness, double[] startPos, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motion, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            if (startPos.Length != 2)
                throw new ArgumentOutOfRangeException("Spatial dimension does not fit particle definition");

            m_Length = length;
            m_Thickness = thickness;
            Aux.TestArithmeticException(length, "Particle length");
            Aux.TestArithmeticException(thickness, "Particle thickness");

            Motion.CharacteristicLength = GetLengthScales().Max();
            Motion.Volume = this.Volume;
            Motion.MomentOfInertia = this.MomentOfInertia;
        }

        [DataMember]
        private readonly double m_Length;
        [DataMember]
        private readonly double m_Thickness;

        /// <summary>
        /// Circumference of an elliptic particle. Approximated with Ramanujan.
        /// </summary>
        public override double Circumference => 2 * m_Length + 2 * m_Thickness;

        /// <summary>
        /// Moment of inertia of an elliptic particle.
        /// </summary>
        override public double MomentOfInertia => (Mass * (m_Length.Pow2() + m_Thickness.Pow2())) / 12;

        /// <summary>
        /// Area occupied by the particle.
        /// </summary>
        public override double Volume => m_Length * m_Thickness;

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        protected override double ParticleLevelSetFunction(double[] X, Vector Postion) {
            double angle = Motion.GetAngle(0);
            double[] position = Postion;
            double[] tempX = X.CloneAs();
            tempX[0] = X[0] * Math.Cos(angle) + X[1] * Math.Sin(angle);
            tempX[1] = X[0] * Math.Sin(angle) + X[1] * Math.Cos(angle);
            double r = -Math.Max(Math.Abs(tempX[0] - position[0]) - m_Length, Math.Abs(tempX[1] - position[1]) - m_Thickness);
            if (double.IsNaN(r) || double.IsInfinity(r))
                throw new ArithmeticException();
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
            Vector orientation = new Vector(Math.Cos(Motion.GetAngle(0)), -Math.Sin(Motion.GetAngle(0)));
            Vector normalOrientation = new Vector(Math.Sin(Motion.GetAngle(0)), Math.Cos(Motion.GetAngle(0)));
            Vector position = Motion.GetPosition(0);
            double a = m_Length + tolerance;
            double b = m_Thickness + tolerance;
            Vector tempX = new Vector( point * orientation, point * normalOrientation );
            return (Math.Abs(tempX[0] - position[0]) < a && Math.Abs(tempX[1] - position[1]) < b);
        }

        /// <summary>
        /// Returns the support point of the particle in the direction specified by a vector.
        /// </summary>
        /// <param name="vector">
        /// A vector. 
        /// </param>
        override public Vector GetSupportPoint(Vector supportVector, Vector Position, Vector Angle, int SubParticleID, double tolerance = 0) {
            Aux.TestArithmeticException(supportVector, "vector in calc of support point");
            if (supportVector.L2Norm() == 0)
                throw new ArithmeticException("The given vector has no length");

            Vector supportPoint = new Vector(supportVector);
            if (Angle.Dim > 1)
                throw new NotImplementedException("Only 2D support");
            double angle = Angle[0]; // hardcoded 2D
            Vector position = new Vector(Position);
            Vector rotVector = new Vector(supportVector);
            rotVector[0] = supportVector[0] * Math.Cos(angle) - supportVector[1] * Math.Sin(angle);
            rotVector[1] = supportVector[0] * Math.Sin(angle) + supportVector[1] * Math.Cos(angle);
            Vector length = new Vector(position);
            length[0] = m_Length * Math.Cos(angle) - m_Thickness * Math.Sin(angle);
            length[1] = m_Length * Math.Sin(angle) + m_Thickness * Math.Cos(angle);
            for(int d = 0; d < position.Dim; d++) {
                supportPoint[d] = Math.Sign(rotVector[d]) * length[d] + position[d];
            }
            return supportPoint;
        }

        /// <summary>
        /// Returns the legnthscales of a particle.
        /// </summary>
        override public double[] GetLengthScales() {
            return new double[] { m_Length, m_Thickness };
        }
    }
}


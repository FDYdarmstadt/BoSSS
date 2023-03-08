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
    public class ParticleShell : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private ParticleShell() : base() {

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
        public ParticleShell(IMotion motion, double length, double height, double thickness, double[] startPos, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motion, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            throw new NotImplementedException("Legacy code, untested, update necessary");
            if (startPos.Length != 2)
                throw new ArgumentOutOfRangeException("Spatial dimension does not fit particle definition");

            m_Length = length;
            m_Thickness = thickness;
            m_Height = height;
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
        [DataMember]
        private readonly double m_Height;

        /// <summary>
        /// The shell is devided into two convex sub particles. Necessary for the GJK-algorithm in the collision model.
        /// </summary>
        public override int NoOfSubParticles => 3;

        /// <summary>
        /// Circumference of an elliptic particle. Approximated with Ramanujan.
        /// </summary>
        public override double Circumference => 2 * m_Length + 4 * m_Height - 2 * m_Thickness;

        /// <summary>
        /// Moment of inertia of an elliptic particle.
        /// </summary>
        override public double MomentOfInertia => (Mass * (m_Length.Pow2() + m_Height.Pow2()) - Mass * (5 * m_Thickness.Pow2() - 2 * (m_Length + m_Height) * m_Thickness) / ((m_Height - m_Thickness) * (m_Length - 2 * m_Thickness))) / 12;

        /// <summary>
        /// Area occupied by the particle.
        /// </summary>
        public override double Volume => m_Length * m_Height - (m_Height - m_Thickness) * (m_Length - 2 * m_Thickness);

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
            // subparticle _
            double[] subPosition = position.CloneAs();
            subPosition[0] = position[0] - (m_Height - m_Thickness) / 2 * Math.Sin(angle);
            subPosition[1] = position[1] - (m_Height - m_Thickness) / 2 * Math.Cos(angle);
            double subR1 = -Math.Max(Math.Abs(tempX[0] - subPosition[0]) - (m_Length - 2 * m_Thickness) / 2, Math.Abs(tempX[1] - subPosition[1]) - m_Thickness / 2);
            // subparticle | left
            subPosition[0] = position[0] - (m_Length - m_Thickness) / 2 * Math.Cos(angle);
            subPosition[1] = position[1] - (m_Length - m_Thickness) / 2 * Math.Sin(angle);
            double subR2 = -Math.Max(Math.Abs(tempX[0] - subPosition[0]) - m_Thickness / 2, Math.Abs(tempX[1] - subPosition[1]) - m_Height / 2);
            // subparticle | right
            subPosition[0] = position[0] - (-m_Length + m_Thickness) / 2 * Math.Cos(angle);
            subPosition[1] = position[1] - (-m_Length + m_Thickness) / 2 * Math.Sin(angle);
            double subR3 = -Math.Max(Math.Abs(tempX[0] - subPosition[0]) - m_Thickness / 2, Math.Abs(tempX[1] - subPosition[1]) - m_Height / 2);
            // complete
            double r = Math.Max(subR2, subR3);
            r = Math.Max(r, subR1);
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
            double angle = Motion.GetAngle(0);
            Vector position = Motion.GetPosition(0);
            double a = (m_Length - 2 * m_Thickness) / 2 + tolerance;
            double b = m_Height / 2 + tolerance;
            double c = m_Thickness / 2 + tolerance;
            Vector tempX = new Vector(point);
            tempX[0] = point[0] * Math.Cos(angle) - point[1] * Math.Sin(angle);
            tempX[1] = point[0] * Math.Sin(angle) + point[1] * Math.Cos(angle);
            // subparticle _
            Vector subPosition = new Vector(position);
            subPosition[0] = position[0] - (m_Height - m_Thickness) / 2 * Math.Sin(angle);
            subPosition[1] = position[1] - (m_Height - m_Thickness) / 2 * Math.Cos(angle);
            if (Math.Abs(tempX[0] - subPosition[0]) < a && Math.Abs(tempX[1] - subPosition[1]) < c)
                return true;
            // subparticle | left
            subPosition[0] = position[0] - (m_Length - m_Thickness) / 2 * Math.Cos(angle);
            subPosition[1] = position[1] - (m_Length - m_Thickness) / 2 * Math.Sin(angle);
            if (Math.Abs(tempX[0] - subPosition[0]) < c && Math.Abs(tempX[1] - subPosition[1]) < b)
                return true;
            // subparticle | right
            subPosition[0] = position[0] - (-m_Length + m_Thickness) / 2 * Math.Cos(angle);
            subPosition[1] = position[1] - (-m_Length + m_Thickness) / 2 * Math.Sin(angle);
            if (Math.Abs(tempX[0] - subPosition[0]) < c && Math.Abs(tempX[1] - subPosition[1]) < b)
                return true;
            else
               return false;
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

            double[] position = Position;
            if (Angle.Dim > 1)
                throw new NotImplementedException("Only 2D support");
            double angle = Angle[0]; // hardcoded 2D
            double[] subPosition = position.CloneAs();
            double[] length = position.CloneAs();

            switch (SubParticleID) {
                case 1:
                    subPosition[0] = position[0] - (m_Height - m_Thickness) / 2 * Math.Sin(angle);
                    subPosition[1] = position[1] - (m_Height - m_Thickness) / 2 * Math.Cos(angle);
                    length[0] = ((m_Length - 2 * m_Thickness) * Math.Cos(angle) - m_Thickness * Math.Sin(angle)) / 2;
                    length[1] = ((m_Length - 2 * m_Thickness) * Math.Sin(angle) + m_Thickness * Math.Cos(angle)) / 2;
                    break;
                case 2:
                    subPosition[0] = position[0] - (m_Length - m_Thickness) / 2 * Math.Cos(angle);
                    subPosition[1] = position[1] - (m_Length - m_Thickness) / 2 * Math.Sin(angle);
                    length[0] = (m_Thickness * Math.Cos(angle) - m_Height * Math.Sin(angle)) / 2;
                    length[1] = (m_Thickness * Math.Sin(angle) + m_Height * Math.Cos(angle)) / 2;
                    break;
                case 3:
                    subPosition[0] = position[0] - (-m_Length + m_Thickness) / 2 * Math.Cos(angle);
                    subPosition[1] = position[1] - (-m_Length + m_Thickness) / 2 * Math.Sin(angle);
                    length[0] = (m_Thickness * Math.Cos(angle) - m_Height * Math.Sin(angle)) / 2;
                    length[1] = (m_Thickness * Math.Sin(angle) + m_Height * Math.Cos(angle)) / 2;
                    break;
            }

            Vector supportPoint = new Vector(supportVector);
            Vector rotVector = new Vector(supportVector);
            rotVector[0] = supportVector[0] * Math.Cos(angle) - supportVector[1] * Math.Sin(angle);
            rotVector[1] = supportVector[0] * Math.Sin(angle) + supportVector[1] * Math.Cos(angle);
            for(int d = 0; d < position.Length; d++) {
                supportPoint[d] = Math.Sign(rotVector[d]) * length[d] + subPosition[d];
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


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

namespace BoSSS.Application.XNSERO_Solver {
    [DataContract]
    [Serializable]
    public class ParticleSuperEllipsoid : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private ParticleSuperEllipsoid() : base() {

        }

        /// <summary>
        /// Constructor for a superellipsoid.
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
        /// <param name="superEllipsoidExponent">
        /// The exponent of the superellipsoid.
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
        public ParticleSuperEllipsoid(IMotion motion, double length, double thickness, int superEllipsoidExponent, double[] startPos, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motion, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            //throw new NotImplementedException("Legacy code, untested, update necessary");
            if (startPos.Length != 2)
                throw new ArgumentOutOfRangeException("Spatial dimension does not fit particle definition");

            m_Length = length;
            m_Thickness = thickness;
            m_Exponent = superEllipsoidExponent;
            Aux.TestArithmeticException(length, "Particle length");
            Aux.TestArithmeticException(thickness, "Particle thickness");
            Aux.TestArithmeticException(superEllipsoidExponent, "super ellipsoid exponent");

            Motion.CharacteristicLength = GetLengthScales().Max();
            Motion.Volume = this.Volume;
            Motion.MomentOfInertia = this.MomentOfInertia;
        }

        [DataMember]
        private readonly double m_Length;
        [DataMember]
        private readonly double m_Thickness;
        [DataMember]
        private readonly double m_Exponent;

        /// <summary>
        /// Circumference. Extremely rough approximation.
        /// </summary>
        public override double Circumference => (2 * m_Length + 2 * m_Thickness + 2 * Math.PI * m_Thickness) / 2;

        /// <summary>
        /// Area occupied by the particle. 
        /// </summary>
        public override double Volume => 4 * m_Length * m_Thickness;// * (SpecialFunctions.Gamma(1 + 1 / m_Exponent)).Pow2() / SpecialFunctions.Gamma(1 + 2 / m_Exponent);

        /// <summary>
        /// Moment of inertia. 
        /// </summary>
        override public double MomentOfInertia => (1 / 4.0) * Mass * (m_Length * m_Length + m_Thickness * m_Thickness);

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        protected override double ParticleLevelSetFunction(double[] X, Vector Postion) {
            double alpha = -(Motion.GetAngle(0));
            double r;
            r = -Math.Pow(((X[0] - Postion[0]) * Math.Cos(alpha) - (X[1] - Postion[1]) * Math.Sin(alpha)) / m_Length, m_Exponent)
                -Math.Pow(((X[0] - Postion[0]) * Math.Sin(alpha) + (X[1] - Postion[1]) * Math.Cos(alpha)) / m_Thickness, m_Exponent)
                + 1;
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
            Vector orientation = new Vector(Math.Cos(Motion.GetAngle(0)), Math.Sin(Motion.GetAngle(0)));
            Vector normalOrientation = new Vector(-Math.Sin(Motion.GetAngle(0)), Math.Cos(Motion.GetAngle(0)));
            Vector distancePointToPosition = point - Position;
            double a = m_Length + tolerance;
            double b = m_Thickness + tolerance;
            double Superellipsoid = Math.Pow(distancePointToPosition * orientation / a, m_Exponent) + Math.Pow(distancePointToPosition * normalOrientation / b, m_Exponent);
            return Superellipsoid < 1;
        }

        /// <summary>
        /// Returns an array with points on the surface of the particle.
        /// </summary>
        /// <param name="hMin">
        /// Minimal cell length. Used to specify the number of surface points.
        /// </param>
        override public MultidimensionalArray GetSurfacePoints(double dAngle, double searchAngle, int subParticleID) {
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported.");
            double angle = Motion.GetAngle(0);
            int noOfCurrentPointAndNeighbours = 3;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(noOfCurrentPointAndNeighbours, SpatialDim);
            for (int j = 0; j < noOfCurrentPointAndNeighbours; j++) {
                double verticalAxis;
                double horizontalAxis;
                double currentAngle = searchAngle + dAngle * (j - 1);
                if (currentAngle < 0)
                    currentAngle += 2 * Math.PI;
                if (currentAngle > 2 * Math.PI)
                    currentAngle -= 2 * Math.PI;
                if (searchAngle + dAngle * (j - 1) <= Math.PI / 2) {
                    verticalAxis = m_Length * Math.Pow(Math.Abs(Math.Cos(searchAngle + dAngle * (j - 1))), 2 / m_Exponent);
                    horizontalAxis = m_Thickness * Math.Pow(Math.Abs(Math.Sin(searchAngle + dAngle * (j - 1))), 2 / m_Exponent);
                }
                else if (searchAngle + dAngle * (j - 1) > Math.PI / 2 && searchAngle + dAngle * (j - 1) <= Math.PI) {
                    verticalAxis = -m_Length * Math.Pow(Math.Abs(Math.Cos(searchAngle + dAngle * (j - 1))), 2 / m_Exponent);
                    horizontalAxis = m_Thickness * Math.Pow(Math.Abs(Math.Sin(searchAngle + dAngle * (j - 1))), 2 / m_Exponent);
                }
                else if (searchAngle + dAngle * (j - 1) > Math.PI && searchAngle + dAngle * (j - 1) <= 3 * Math.PI / 2) {
                    verticalAxis = -m_Length * Math.Pow(Math.Abs(Math.Cos(searchAngle + dAngle * (j - 1))), 2 / m_Exponent);
                    horizontalAxis = -m_Thickness * Math.Pow(Math.Abs(Math.Sin(searchAngle + dAngle * (j - 1))), 2 / m_Exponent);
                }
                else  {
                    verticalAxis = m_Length * Math.Pow(Math.Abs(Math.Cos(searchAngle + dAngle * (j - 1))), 2 / m_Exponent);
                    horizontalAxis = -m_Thickness * Math.Pow(Math.Abs(Math.Sin(searchAngle + dAngle * (j - 1))), 2 / m_Exponent);
                }
                SurfacePoints[j, 0] = (verticalAxis * Math.Cos(angle) - horizontalAxis * Math.Sin(angle));
                SurfacePoints[j, 1] = (verticalAxis * Math.Sin(angle) + horizontalAxis * Math.Cos(angle));
            }
            //int noOfCurrentPointWithNeighbours = 3;
            //MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(noOfCurrentPointWithNeighbours, SpatialDim);
            //for (int j = 0; j < QuarterSurfacePoints; j++) {
            //    SurfacePoints[0, j, 0] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Length * Math.Cos(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Thickness * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
            //    SurfacePoints[0, j, 1] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Length * Math.Sin(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Thickness * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
            //    SurfacePoints[0, 2 * QuarterSurfacePoints + j - 1, 0] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Length) * Math.Cos(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Thickness * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
            //    SurfacePoints[0, 2 * QuarterSurfacePoints + j - 1, 1] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Length) * Math.Sin(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Thickness * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1]; ;
            //}
            //for (int j = 1; j < QuarterSurfacePoints; j++) {
            //    SurfacePoints[0, 2 * QuarterSurfacePoints - j - 1, 0] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Length) * Math.Cos(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Thickness * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
            //    SurfacePoints[0, 2 * QuarterSurfacePoints - j - 1, 1] = (-(Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Length) * Math.Sin(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Thickness * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
            //    SurfacePoints[0, 4 * QuarterSurfacePoints - j - 2, 0] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Length * Math.Cos(Motion.GetAngle(0)) + Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Thickness * Math.Sin(Motion.GetAngle(0))) + Motion.GetPosition(0)[0];
            //    SurfacePoints[0, 4 * QuarterSurfacePoints - j - 2, 1] = (Math.Pow(Math.Cos(Infinitisemalangle[j]), 2 / m_Exponent) * m_Length * Math.Sin(Motion.GetAngle(0)) - Math.Pow(Math.Sin(Infinitisemalangle[j]), 2 / m_Exponent) * m_Thickness * Math.Cos(Motion.GetAngle(0))) + Motion.GetPosition(0)[1];
            //}
            return SurfacePoints;
        }

        /// <summary>
        /// Returns the legnthscales of a particle.
        /// </summary>
        override public double[] GetLengthScales() {
            return new double[] { m_Length, m_Thickness };
        }
    }
}


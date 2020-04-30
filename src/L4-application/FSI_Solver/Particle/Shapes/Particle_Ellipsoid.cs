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
using ilPSP.Utils;
using FSI_Solver;
using System.Diagnostics;

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_Ellipsoid : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Ellipsoid() : base() {

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
        public Particle_Ellipsoid(ParticleMotionInit motionInit, double length = 4, double thickness = 1, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            m_Length = length;
            m_Thickness = thickness;
            Aux.TestArithmeticException(length, "Particle length");
            Aux.TestArithmeticException(thickness, "Particle thickness");

            Motion.SetParticleMaxLengthscale(GetLengthScales().Max());
            Motion.SetParticleArea(Area);
            Motion.SetParticleMomentOfInertia(MomentOfInertia);
        }

        [DataMember]
        private readonly double m_Length;
        [DataMember]
        private readonly double m_Thickness;

        /// <summary>
        /// Circumference of an elliptic particle. Approximated with Ramanujan.
        /// </summary>
        public override double Circumference => Math.PI * ((m_Length + m_Thickness) + (3 * (m_Length - m_Thickness).Pow2()) / (10 * (m_Length + m_Thickness) + Math.Sqrt(m_Length.Pow2() + 14 * m_Length * m_Thickness + m_Thickness.Pow2())));

        /// <summary>
        /// Moment of inertia of an elliptic particle.
        /// </summary>
        override public double MomentOfInertia => (1 / 4.0) * (Mass_P * (m_Length * m_Length + m_Thickness * m_Thickness));

        /// <summary>
        /// Area occupied by the particle.
        /// </summary>
        public override double Area => m_Length * m_Thickness * Math.PI;

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        public override double LevelSetFunction(double[] X) {
            double angle = -Motion.GetAngle(0);
            double[] position = Motion.GetPosition(0);
            double r = -(((X[0] - position[0]) * Math.Cos(angle) - (X[1] - position[1]) * Math.Sin(angle)) / m_Length).Pow2()
                       - (((X[0] - position[0]) * Math.Sin(angle) + (X[1] - position[1]) * Math.Cos(angle)) / m_Thickness).Pow2()
                       + 1.0;
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
        public override bool Contains(Vector point, double tolerance = 0) {
            Vector orientation = new Vector(Math.Cos(Motion.GetAngle(0)), Math.Sin(Motion.GetAngle(0)));
            Vector normalOrientation = new Vector(-Math.Sin(Motion.GetAngle(0)), Math.Cos(Motion.GetAngle(0)));
            Vector position = Motion.GetPosition(0);
            double a = m_Length + tolerance;
            double b = m_Thickness + tolerance;
            double Ellipse = ((point - position) * orientation).Pow2() / a.Pow2() + ((point - position) * normalOrientation).Pow2() / b.Pow2();
            return Ellipse < 1;
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
            int noOfCurrentPointWithNeighbours = 3;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(noOfCurrentPointWithNeighbours, SpatialDim);
            for (int j = 0; j < noOfCurrentPointWithNeighbours; j++) {
                double verticalAxis = m_Length * Math.Cos(searchAngle + dAngle * (j - 1));
                double horizontalAxis = m_Thickness * Math.Sin(searchAngle + dAngle * (j - 1));
                SurfacePoints[j, 0] = (verticalAxis * Math.Cos(angle) - horizontalAxis * Math.Sin(angle));// + position[0];
                SurfacePoints[j, 1] = (verticalAxis * Math.Sin(angle) + horizontalAxis * Math.Cos(angle));// + position[1];
            }
            return SurfacePoints;
        }

        /// <summary>
        /// Returns the support point of the particle in the direction specified by a vector.
        /// </summary>
        /// <param name="vector">
        /// A vector. 
        /// </param>
        override public Vector GetSupportPoint(Vector supportVector, int SubParticleID) {
            Aux = new FSI_Auxillary();
            Aux.TestArithmeticException(supportVector, "vector in calc of support point");
            if (supportVector.L2Norm() == 0)
                throw new ArithmeticException("The given vector has no length");

            Vector SupportPoint = new Vector(SpatialDim);
            double angle = Motion.GetAngle(0);
            Vector position = new Vector(Motion.GetPosition(0));

            double[,] rotMatrix = new double[2, 2];
            rotMatrix[0, 0] = m_Length * Math.Cos(angle);
            rotMatrix[0, 1] = -m_Thickness * Math.Sin(angle);
            rotMatrix[1, 0] = m_Length * Math.Sin(angle);
            rotMatrix[1, 1] = m_Thickness * Math.Cos(angle);
            double[,] transposeRotMatrix = rotMatrix.CloneAs();
            transposeRotMatrix[0, 1] = rotMatrix[1, 0];
            transposeRotMatrix[1, 0] = rotMatrix[0, 1];

            double[] rotVector = new double[2];
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    rotVector[i] += transposeRotMatrix[i, j] * supportVector[j];
                }
            }
            rotVector.ScaleV(1 / rotVector.L2Norm());

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    SupportPoint[i] += rotMatrix[i, j] * rotVector[j];
                }
                SupportPoint[i] += position[i];
            }
            return SupportPoint;
        }

        /// <summary>
        /// Returns the legnthscales of a particle.
        /// </summary>
        override public double[] GetLengthScales() {
            return new double[] { m_Length, m_Thickness };
        }

        public override object Clone() {
            Particle clonedParticle = new Particle_Ellipsoid(MotionInitializer,
                                                             m_Length,
                                                             m_Thickness,
                                                             Motion.GetPosition(),
                                                             Motion.GetAngle() * 360 / (2 * Math.PI),
                                                             ActiveStress,
                                                             Motion.GetTranslationalVelocity(),
                                                             Motion.GetRotationalVelocity());
            clonedParticle.IsMaster = IsMaster;

            return clonedParticle;
        }
    }
}


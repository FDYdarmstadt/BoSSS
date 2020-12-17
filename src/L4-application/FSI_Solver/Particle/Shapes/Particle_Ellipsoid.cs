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
        /// <param name="halfAxisA">
        /// The length of the horizontal halfaxis.
        /// </param>
        /// <param name="halfAxisB">
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
        public Particle_Ellipsoid(InitializeMotion motionInit, double halfAxisA = 4, double halfAxisB = 1, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            m_Length = halfAxisA;
            m_Thickness = halfAxisB;
            Aux.TestArithmeticException(halfAxisA, "Particle length");
            Aux.TestArithmeticException(halfAxisB, "Particle thickness");

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
            double[] position = Motion.GetPosition(0);
            Vector orientation = Motion.orientationVector;
            double r = -(((X[0] - position[0]) * orientation[0] + (X[1] - position[1]) * orientation[1]) / m_Length).Pow2()
                       -(((X[0] - position[0]) * orientation[1] - (X[1] - position[1]) * orientation[0]) / m_Thickness).Pow2()
                       + 1.0;
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
            Vector orientation = new Vector(Motion.orientationVector);
            Vector position = Motion.GetPosition(0);
            double a = m_Length + tolerance;
            double b = m_Thickness + tolerance;
            double Ellipse = ((point[0] - position[0]) * orientation[0] + (point[1] - position[1]) * orientation[1]).Pow2() / a.Pow2() + ((point[0] - position[0]) * orientation[1] - (point[1] - position[1]) * orientation[0]).Pow2() / b.Pow2();
            return Ellipse < 1;
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
            Vector orientation = Motion.orientationVector;
            Vector position = new Vector(Motion.GetPosition(0));

            double[,] rotMatrix = new double[2, 2];
            rotMatrix[0, 0] = m_Length * orientation[0];
            rotMatrix[0, 1] = -m_Thickness * orientation[1];
            rotMatrix[1, 0] = m_Length * orientation[1];
            rotMatrix[1, 1] = m_Thickness * orientation[0];
            double[,] transposeRotMatrix = rotMatrix.CloneAs();
            transposeRotMatrix[0, 1] = rotMatrix[1, 0];
            transposeRotMatrix[1, 0] = rotMatrix[0, 1];

            double[] rotVector = new double[2];
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    rotVector[i] += transposeRotMatrix[i, j] * supportVector[j];
                }
            }
            rotVector = rotVector.Normalize();

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
            return new Particle_Ellipsoid(MotionInitializer, m_Length, m_Thickness, Motion.GetPosition(), Motion.GetAngle() * 360 / (2 * Math.PI), ActiveStress, Motion.GetTranslationalVelocity(), Motion.GetRotationalVelocity()) {
                IsMaster = IsMaster
            };
        }
    }
}


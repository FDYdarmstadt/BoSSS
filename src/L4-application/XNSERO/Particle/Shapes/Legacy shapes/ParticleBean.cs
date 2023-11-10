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
    public class ParticleBean : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private ParticleBean() : base() {

        }

        /// <summary>
        /// Constructor for a bean. Warning: Legacy code.
        /// </summary>
        /// <param name="motionInit">
        /// Initializes the motion parameters of the particle (which model to use, whether it is a dry simulation etc.)
        /// </param>
        /// <param name="radius">
        /// The main lengthscale of the bean. 
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
        public ParticleBean(IMotion motion, double radius, double[] startPos, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motion, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            throw new NotImplementedException("Legacy code, untested, update necessary");
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
        /// Circumference. Approximated with sphere.
        /// </summary>
        public override double Circumference => 2 * Math.PI * m_Radius;

        /// <summary>
        /// Area occupied by the particle.
        /// </summary>
        public override double Volume => Math.PI * m_Radius * m_Radius;

        /// <summary>
        /// Moment of inertia. Approximated with sphere.
        /// </summary>
        override public double MomentOfInertia => (1 / 2.0) * (Mass * m_Radius * m_Radius);

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        protected override double ParticleLevelSetFunction(double[] X, Vector Postion) {
            double alpha = -Motion.GetAngle(0);
            double[] position = Postion;
            double a = 3.0 * m_Radius.Pow2();
            double b = 1.0 * m_Radius.Pow2();
            return -((((X[0] - position[0]) * Math.Cos(alpha) - (X[1] - position[1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - position[0]) * Math.Sin(alpha) + (X[1] - position[1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - position[0]) * Math.Cos(alpha) - (X[1] - position[1]) * Math.Sin(alpha)).Pow(3) - b * ((X[0] - position[0]) * Math.Sin(alpha) + (X[1] - position[1]) * Math.Cos(alpha)).Pow2());
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
            double alpha = Motion.GetAngle(0);
            Vector position = Motion.GetPosition(0);
            // only for rectangular cells
            double radiusTolerance = 1.0 + tolerance;
            double a = 4.0 * radiusTolerance.Pow2();
            double b = 1.0 * radiusTolerance.Pow2();
            if (-((((point[0] - position[0]) * Math.Cos(alpha) - (point[1] - position[1]) * Math.Sin(alpha)).Pow(2) + ((point[0] - position[0]) * Math.Sin(alpha) + (point[1] - position[1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((point[0] - position[0]) * Math.Cos(alpha) - (point[1] - position[1]) * Math.Sin(alpha)).Pow(3) - b * ((point[0] - position[0]) * Math.Sin(alpha) + (point[1] - position[1]) * Math.Cos(alpha)).Pow2()) > 0) {
                return true;
            }
            return false;
        }

        /// <summary>
        /// Returns the legnthscales of a particle.
        /// </summary>
        override public double[] GetLengthScales() {
            return new double[] { m_Radius, m_Radius };
        }
    }
}


﻿/* =======================================================================
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

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_Sphere : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Sphere() : base() {

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
        public Particle_Sphere(ParticleMotionInit motionInit, double radius, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {

            m_Radius = radius;
            Aux.TestArithmeticException(radius, "Particle radius");

            Motion.GetParticleLengthscale(radius);
            Motion.GetParticleMinimalLengthscale(radius);
            Motion.GetParticleArea(Area);
            Motion.GetParticleMomentOfInertia(MomentOfInertia);

        }

        [DataMember]
        private readonly double m_Radius;

        /// <summary>
        /// Area occupied by the particle.
        /// </summary>
        public override double Area => Math.PI * m_Radius.Pow2();

        /// <summary>
        /// Circumference. 
        /// </summary>
        public override double Circumference => 2 * Math.PI * m_Radius;

        /// <summary>
        /// Moment of inertia. 
        /// </summary>
        override public double MomentOfInertia => (1 / 2.0) * (Mass_P * m_Radius.Pow2());

        /// <summary>
        /// Level set function of the particle.
        /// </summary>
        /// <param name="X">
        /// The current point.
        /// </param>
        public override double LevelSetFunction(double[] X) {
            double x0 = Motion.GetPosition(0)[0];
            double y0 = Motion.GetPosition(0)[1];
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
        public override bool Contains(Vector point, double tolerance = 0) {
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
        override public Vector GetSupportPoint(Vector supportVector, int SubParticleID) {
            double length = Math.Sqrt(supportVector[0].Pow2() + supportVector[1].Pow2());
            double CosT = supportVector[0] / length;
            double SinT = supportVector[1] / length;
            Vector SupportPoint = new Vector(SpatialDim);
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");
            SupportPoint[0] = CosT * m_Radius + Motion.GetPosition(0)[0];
            SupportPoint[1] = SinT * m_Radius + Motion.GetPosition(0)[1];
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
    }
}


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

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_Bean : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Bean() : base() {

        }
        public Particle_Bean(ParticleMotionInit motionInit, double radius, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            radius_P = radius;
            Motion.GetParticleLengthscale(radius);
            Motion.GetParticleArea(Area_P());
            Motion.GetParticleMomentOfInertia(MomentOfInertia_P);
        }

        /// <summary>
        /// Radius of the particle. Not necessary for particles defined by their length and thickness
        /// </summary>
        private readonly double radius_P;

        protected override double Circumference_P {
            get {
                return 2 * Math.PI * radius_P;
            }
        }
        public override double Area_P() {
            // not correct area
            return Math.PI * radius_P * radius_P;
        }
        override public double MomentOfInertia_P {
            get {
                // not correct moment of inertia
                return (1 / 2.0) * (Mass_P * radius_P * radius_P);
            }
        }

        public override double LevelSetFunction(double[] X) {
            double alpha = -(Motion.Angle[0]);
            double a = 3.0 * radius_P.Pow2();
            double b = 1.0 * radius_P.Pow2();
            return -((((X[0] - Motion.position[0][0]) * Math.Cos(alpha) - (X[1] - Motion.position[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - Motion.position[0][0]) * Math.Sin(alpha) + (X[1] - Motion.position[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - Motion.position[0][0]) * Math.Cos(alpha) - (X[1] - Motion.position[0][1]) * Math.Sin(alpha)).Pow(3) - b * ((X[0] - Motion.position[0][0]) * Math.Sin(alpha) + (X[1] - Motion.position[0][1]) * Math.Cos(alpha)).Pow2());
        }

        public override bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false) {
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double radiusTolerance = !WithoutTolerance ? 1.0 + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : 1;
            double a = 4.0 * radiusTolerance.Pow2();
            double b = 1.0 * radiusTolerance.Pow2();
            if (-((((point[0] - Motion.position[0][0]) * Math.Cos(Motion.Angle[0]) - (point[1] - Motion.position[0][1]) * Math.Sin(Motion.Angle[0])).Pow(2) + ((point[0] - Motion.position[0][0]) * Math.Sin(Motion.Angle[0]) + (point[1] - Motion.position[0][1]) * Math.Cos(Motion.Angle[0])).Pow(2)).Pow2() - a * ((point[0] - Motion.position[0][0]) * Math.Cos(Motion.Angle[0]) - (point[1] - Motion.position[0][1]) * Math.Sin(Motion.Angle[0])).Pow(3) - b * ((point[0] - Motion.position[0][0]) * Math.Sin(Motion.Angle[0]) + (point[1] - Motion.position[0][1]) * Math.Cos(Motion.Angle[0])).Pow2()) > 0) {
                return true;
            }
            return false;
        }

        override public double[] GetLengthScales() {
            return new double[] { radius_P, radius_P };
        }
    }
}


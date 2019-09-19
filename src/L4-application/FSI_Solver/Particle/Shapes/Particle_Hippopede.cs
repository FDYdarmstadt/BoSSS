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
using ilPSP.Utils;

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_Hippopede : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Hippopede() : base() {

        }

        public Particle_Hippopede(ParticleMotionInit motionInit, double length, double thickness, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, activeStress, startAngl, startTransVelocity, startRotVelocity) {
            length_P = length;
            thickness_P = thickness;
            Motion.GetParticleLengthscale(GetLengthScales().Max());
            Motion.GetParticleArea(Area_P());
            Motion.GetParticleMomentOfInertia(MomentOfInertia_P);
        }

        /// <summary>
        /// Length of an elliptic particle.
        /// </summary>
        [DataMember]
        public double length_P;

        /// <summary>
        /// Thickness of an elliptic particle.
        /// </summary>
        [DataMember]
        public double thickness_P;

        protected override double Circumference_P {
            get {// not correct circumference
                return Math.PI * length_P * thickness_P;
            }
        }
        public override double Area_P() {
            // not correct area
            return Math.PI * length_P * thickness_P;
        }
        override public double MomentOfInertia_P {
            get {
                // not correct moment of inertia
                return (1 / 2.0) * (Mass_P * length_P * thickness_P);
            }
        }

        public override double LevelSetFunction(double[] X) {
            double a = length_P;
            double b = thickness_P;
            double alpha = -(Motion.angle[0]);
            return -((((X[0] - Motion.position[0][0]) * Math.Cos(alpha) - (X[1] - Motion.position[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - Motion.position[0][0]) * Math.Sin(alpha) + (X[1] - Motion.position[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - Motion.position[0][0]) * Math.Cos(alpha) - (X[1] - Motion.position[0][1]) * Math.Sin(alpha)).Pow2() - b * ((X[0] - Motion.position[0][0]) * Math.Sin(alpha) + (X[1] - Motion.position[0][1]) * Math.Cos(alpha)).Pow2());
        }

        

        public override bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false) {
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double a = !WithoutTolerance ? length_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : length_P;
            double b = !WithoutTolerance ? thickness_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : thickness_P;
            if (-((((point[0] - Motion.position[0][0]) * Math.Cos(Motion.angle[0]) - (point[1] - Motion.position[0][1]) * Math.Sin(Motion.angle[0])).Pow(2) + ((point[0] - Motion.position[0][0]) * Math.Sin(Motion.angle[0]) + (point[1] - Motion.position[0][1]) * Math.Cos(Motion.angle[0])).Pow(2)).Pow2() - a * ((point[0] - Motion.position[0][0]) * Math.Cos(Motion.angle[0]) - (point[1] - Motion.position[0][1]) * Math.Sin(Motion.angle[0])).Pow2() - b * ((point[0] - Motion.position[0][0]) * Math.Sin(Motion.angle[0]) + (point[1] - Motion.position[0][1]) * Math.Cos(Motion.angle[0])).Pow2()) > 0) {
                return true;
            }
            return false;
        }

        public override bool ParticleInternalCell(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false) {
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double a = !WithoutTolerance ? length_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : length_P;
            double b = !WithoutTolerance ? thickness_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : thickness_P;
            if (-((((point[0] - Motion.position[0][0]) * Math.Cos(Motion.angle[0]) - (point[1] - Motion.position[0][1]) * Math.Sin(Motion.angle[0])).Pow(2) + ((point[0] - Motion.position[0][0]) * Math.Sin(Motion.angle[0]) + (point[1] - Motion.position[0][1]) * Math.Cos(Motion.angle[0])).Pow(2)).Pow2() - a * ((point[0] - Motion.position[0][0]) * Math.Cos(Motion.angle[0]) - (point[1] - Motion.position[0][1]) * Math.Sin(Motion.angle[0])).Pow2() - b * ((point[0] - Motion.position[0][0]) * Math.Sin(Motion.angle[0]) + (point[1] - Motion.position[0][1]) * Math.Cos(Motion.angle[0])).Pow2()) > 0) {
                return true;
            }
            return false;
        }

        override public double[] GetLengthScales() {
            return new double[] { length_P, thickness_P };
        }
    }
}


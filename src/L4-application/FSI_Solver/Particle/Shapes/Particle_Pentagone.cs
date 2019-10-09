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
using ilPSP.Utils;

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_Pentagone : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Pentagone() : base() {

        }

        public Particle_Pentagone(ParticleMotionInit motionInit, double width, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            width_P = width;
            Motion.GetParticleLengthscale(width);
            Motion.GetParticleArea(Area_P());
            Motion.GetParticleMomentOfInertia(MomentOfInertia_P);

        }

        /// <summary>
        /// Radius of the particle. Not necessary for particles defined by their length and thickness
        /// </summary>
        private readonly double width_P;

        public override double Area_P() {
            // not correct area
            return (5 * width_P.Pow2()) / 4;
        }
        protected override double Circumference_P {
            get {
                return width_P * 4.41421356;
            }
        }
        override public double MomentOfInertia_P {
            get {
                return Math.Pow(width_P, 4) * 0.2887963;
            }
        }

        public override double LevelSetFunction(double[] X) { // attention: no dependency on angle...
            double r;
            r = Math.Max(X[0] - Motion.GetPosition(0)[0] - width_P, Motion.GetPosition(0)[0] - width_P - X[0]);
            r = Math.Max(r, Motion.GetPosition(0)[1] - 0.5 * width_P - X[1]);
            r = Math.Max(r, Motion.GetPosition(0)[0] - width_P - X[1] - 1.5 * X[0]) + Math.Max(r, X[1] - Motion.GetPosition(0)[1] - 0.5 * width_P);
            r = -r;
            return r;
        }

        public override bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false) {
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double radiusTolerance = !WithoutTolerance ? width_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : 1;
            var distance = point.L2Distance(Motion.GetPosition(0));
            if (distance < (radiusTolerance)) {
                return true;
            }
            return false;
        }

        public override MultidimensionalArray GetSurfacePoints(double hMin) {
            if (spatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            int NoOfSurfacePoints = Convert.ToInt32(20 * Circumference_P / hMin) + 1;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSurfacePoints, spatialDim);
            double[] InfinitisemalAngle = GenericBlas.Linspace(0, 2 * Math.PI, NoOfSurfacePoints + 1);
            if (Math.Abs(10 * Circumference_P / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");

            for (int j = 0; j < NoOfSurfacePoints; j++) {
                SurfacePoints[j, 0] = Math.Cos(InfinitisemalAngle[j]) * width_P + Motion.GetPosition(0)[0];
                SurfacePoints[j, 1] = Math.Sin(InfinitisemalAngle[j]) * width_P + Motion.GetPosition(0)[1];
            }
            return SurfacePoints;
        }

        override public double[] GetLengthScales() {
            return new double[] { width_P, width_P };
        }
    }
}

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
    public class Particle_Ellipsoid : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Ellipsoid() : base() {

        }
        public Particle_Ellipsoid(ParticleMotionInit motionInit, double length = 4, double thickness = 1, double[] startPos = null, double startAngl = 0, double activeStress = 0, double[] startTransVelocity = null, double startRotVelocity = 0) : base(motionInit, startPos, startAngl, activeStress, startTransVelocity, startRotVelocity) {
            length_P = length;
            thickness_P = thickness;
            Motion.GetParticleLengthscale(GetLengthScales().Max());
            Motion.GetParticleArea(Area_P());
            Motion.GetParticleMomentOfInertia(MomentOfInertia_P);
        }

        /// <summary>
        /// Length of an elliptic particle.
        /// </summary>
        private readonly double length_P;

        /// <summary>
        /// Thickness of an elliptic particle.
        /// </summary>
        private readonly double thickness_P;

        public override double Area_P() {
            double a = length_P * thickness_P * Math.PI;
            if (a <= 0.0 || double.IsNaN(a) || double.IsInfinity(a))
                throw new ArithmeticException("Ellipsoid volume/area is " + a);
            return a;
        }

        protected override double Circumference_P {
            get {
                return Math.PI * ((length_P + thickness_P) + (3 * (length_P - thickness_P).Pow2()) / (10 * (length_P + thickness_P) + Math.Sqrt(length_P.Pow2() + 14 * length_P * thickness_P + thickness_P.Pow2())));
            }
        }

        override public double MomentOfInertia_P {
            get {
                return (1 / 4.0) * (Mass_P * (length_P * length_P + thickness_P * thickness_P));
            }
        }

        public override double LevelSetFunction(double[] X) {
            double angle = -Motion.GetAngle(0);
            double[] position = Motion.GetPosition(0);
            double r = -(((X[0] - position[0]) * Math.Cos(angle) - (X[1] - position[1]) * Math.Sin(angle)) / length_P).Pow2()
                       - (((X[0] - position[0]) * Math.Sin(angle) + (X[1] - position[1]) * Math.Cos(angle)) / thickness_P).Pow2()
                       + 1.0;
            if (double.IsNaN(r) || double.IsInfinity(r))
                throw new ArithmeticException();
            return r;
        }

        public override bool Contains(double[] point, double h_min, double h_max = 0, bool WithoutTolerance = false) {
            double angle = Motion.GetAngle(0);
            double[] position = Motion.GetPosition(0);
            // only for rectangular cells
            if (h_max == 0)
                h_max = h_min;
            double radiusTolerance = 1;
            double a = !WithoutTolerance ? length_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : length_P;
            double b = !WithoutTolerance ? thickness_P + Math.Sqrt(h_max.Pow2() + h_min.Pow2()) : thickness_P;
            double Ellipse = ((point[0] - position[0]) * Math.Cos(angle) + (point[1] - position[1]) * Math.Sin(angle)).Pow2() / a.Pow2() + (-(point[0] - position[0]) * Math.Sin(angle) + (point[1] - position[1]) * Math.Cos(angle)).Pow2() / b.Pow2();
            if (Ellipse < radiusTolerance) {
                return true;
            }
            else
                return false;
        }

        override public MultidimensionalArray GetSurfacePoints(double hMin) {
            if (spatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");
            double angle = Motion.GetAngle(0);
            double[] position = Motion.GetPosition(0);
            int NoOfSurfacePoints = Convert.ToInt32(5 * Circumference_P / hMin);
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSubParticles, NoOfSurfacePoints, spatialDim);
            double[] InfinitisemalAngle = GenericBlas.Linspace(0, Math.PI * 2, NoOfSurfacePoints + 1);
            if (Math.Abs(10 * Circumference_P / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");
            for (int j = 0; j < NoOfSurfacePoints; j++) {
                double temp0 = Math.Cos(InfinitisemalAngle[j]) * length_P;
                double temp1 = Math.Sin(InfinitisemalAngle[j]) * thickness_P;
                SurfacePoints[0, j, 0] = (temp0 * Math.Cos(angle) - temp1 * Math.Sin(angle)) + position[0];
                SurfacePoints[0, j, 1] = (temp0 * Math.Sin(angle) + temp1 * Math.Cos(angle)) + position[1];

            }
            return SurfacePoints;
        }

        override public void GetSupportPoint(int SpatialDim, double[] Vector, out double[] SupportPoint) {
            SupportPoint = new double[SpatialDim];
            double angle = Motion.GetAngle(0);
            double[] position = Motion.GetPosition(0);
            if (double.IsNaN(Vector[0]) || double.IsNaN(Vector[1]))
                throw new ArithmeticException("Error trying to calculate VectorVectorVector Value:  " + Vector[0] + " VectorVectorVector " + Vector[1]);
            double[,] B = new double[2, 2];
            B[0, 0] = length_P * Math.Cos(angle);
            if (double.IsNaN(B[0, 0]))
                throw new ArithmeticException("Error trying to calculate Angle Value:  " + B[0, 0] + " length_P " + length_P + " Math.Cos(Angle): " + Math.Cos(angle) + " Angle " + angle);
            B[0, 1] = -thickness_P * Math.Sin(angle);
            B[1, 0] = length_P * Math.Sin(angle);
            B[1, 1] = thickness_P * Math.Cos(angle);
            if (double.IsNaN(angle))
                throw new ArithmeticException("Error trying to calculate Angle Value:  " + angle + " Angle " + angle);
            double[,] BT = B.CloneAs();
            BT[0, 1] = B[1, 0];
            BT[1, 0] = B[0, 1];
            double[] temp = new double[2];
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    temp[i] += BT[i, j] * Vector[j];
                    if (double.IsNaN(temp[0]) || double.IsNaN(temp[1]))
                        throw new ArithmeticException("Error trying to calculate temp Value:  " + temp[i] + "VectorVectorVector Value: " + Vector[j] + " BT[i, j] " + BT[i, j] + " afuduiof" + i + j);
                }
            }
            double BetragTemp = Math.Sqrt(temp[0].Pow2() + temp[1].Pow2());
            if (double.IsNaN(temp[0]) || double.IsNaN(temp[1]))
                throw new ArithmeticException("Error trying to calculate temp Value:  " + temp[0] + " temp " + temp[1] + "VectorVectorVector Value: " + Vector[0] + " VectorVectorVector " + Vector[1]);
            for (int i = 0; i < 2; i++) {
                if (BetragTemp != 0)
                    temp[i] = temp[i] / BetragTemp;
            }
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    SupportPoint[i] += B[i, j] * temp[j];
                    if (double.IsNaN(SupportPoint[i]) || double.IsInfinity(SupportPoint[i]))
                        throw new ArithmeticException("Error trying to calculate SupportPoint Value:  " + SupportPoint[i] + "temp Value: " + temp[j] + " B[i, j] " + B[i, j] + "BetragTemp" + BetragTemp + "ij" + i + j);
                }
                SupportPoint[i] += position[i];
            }
        }

        override public double[] GetLengthScales() {
            return new double[] { length_P, thickness_P };
        }
    }
}


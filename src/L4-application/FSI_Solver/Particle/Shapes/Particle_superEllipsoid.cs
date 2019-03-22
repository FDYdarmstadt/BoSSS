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
using BoSSS.Foundation.XDG;
using ilPSP;
using BoSSS.Foundation.Grid;
using MathNet.Numerics;
using System.Diagnostics;

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_superEllipsoid : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_superEllipsoid() : base() {

        }

        public Particle_superEllipsoid(int Dim, double[] startPos = null, double startAngl = 0) : base(Dim, startPos, startAngl) {


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

        /// <summary>
        /// Exponent of the super ellipsoid. Higher exponent leads to a more "squary" appearance.
        /// </summary>
        [DataMember]
        public int superEllipsoidExponent;

        /// <summary>
        /// %
        /// </summary>
        protected override double AverageDistance {
            get {
                throw new NotImplementedException("todo");
            }
        }

        protected override double Circumference_P {
            get {
                return (2 * length_P + 2 * thickness_P + 2 * Math.PI * thickness_P) / 2;
            }
        }

        protected override double Area_P {
            get {
                return 4 * length_P * thickness_P * (SpecialFunctions.Gamma(1 + 1 / superEllipsoidExponent)).Pow2() / SpecialFunctions.Gamma(1 + 2 / superEllipsoidExponent);
            }

        }
        override public double MomentOfInertia_P {
            get {
                return (1 / 4.0) * Mass_P * (length_P * length_P + thickness_P * thickness_P);
            }
        }
        //override public void UpdateLevelSetFunction()
        //{
        //    double alpha = -(Angle[0]);
        //    Phi_P = delegate (double[] X, double t) {
        //        double r;
        //        r = -Math.Pow(((X[0] - Position[0][0]) * Math.Cos(alpha) - (X[1] - Position[0][1]) * Math.Sin(alpha)) / length_P, superEllipsoidExponent) + -Math.Pow(((X[0] - Position[0][0]) * Math.Sin(alpha) + (X[1] - Position[0][1]) * Math.Cos(alpha)) / thickness_P, 
        //            superEllipsoidExponent) + 1;
        //        if (double.IsNaN(r) || double.IsInfinity(r))
        //            throw new ArithmeticException();
        //        return r;
        //    };
        //}
        public override double Phi_P(double[] X) {
            double alpha = -(Angle[0]);
            double r;
            r = -Math.Pow(
                        ((X[0] - Position[0][0]) * Math.Cos(alpha) - (X[1] - Position[0][1]) * Math.Sin(alpha)) / length_P,
                        superEllipsoidExponent)
                - Math.Pow(
                    ((X[0] - Position[0][0]) * Math.Sin(alpha) + (X[1] - Position[0][1]) * Math.Cos(alpha)) / thickness_P,
                    superEllipsoidExponent)
                + 1;
            if (double.IsNaN(r) || double.IsInfinity(r))
                throw new ArithmeticException();
            return r;
        }
        

        override public CellMask CutCells_P(LevelSetTracker LsTrk) {
            // tolerance is very important
            var radiusTolerance = Math.Min(length_P, thickness_P) + LsTrk.GridDat.Cells.h_minGlobal;// +2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());

            CellMask cellCollection;
            CellMask cells = null;
            double alpha = -(Angle[0]);
            cells = CellMask.GetCellMask(LsTrk.GridDat, 
                X => -(((((X[0] - Position[0][0]) * Math.Cos(alpha) - (X[1] - Position[0][1]) * Math.Sin(alpha)) / length_P).Pow(superEllipsoidExponent) + (((X[0] - Position[0][0]) * Math.Sin(alpha) + (X[1] - Position[0][1]) * Math.Cos(alpha)) / thickness_P).Pow(superEllipsoidExponent)) - radiusTolerance.Pow(superEllipsoidExponent)) > 0);

            CellMask allCutCells = LsTrk.Regions.GetCutCellMask();
            cellCollection = cells.Intersect(allCutCells);
            return cellCollection;
        }
        override public bool Contains(double[] point, LevelSetTracker LsTrk) {
            // only for squared cells
            double radiusTolerance = Math.Max(length_P, thickness_P) + 2.0 * Math.Sqrt(2 * LsTrk.GridDat.Cells.h_minGlobal.Pow2());
            if (-Math.Pow(((point[0] - Position[0][0]) * Math.Cos(Angle[0]) - (point[1] - Position[0][1]) * Math.Sin(Angle[0])) / length_P, superEllipsoidExponent) + -Math.Pow(((point[0] - Position[0][0]) * Math.Sin(Angle[0]) + (point[1] - Position[0][1]) * Math.Cos(Angle[0])) / thickness_P, superEllipsoidExponent) + radiusTolerance.Pow(superEllipsoidExponent) > 0) {
                return true;
            }
            return false;
        }

        override public double ComputeParticleRe(double mu_Fluid) {
            double particleReynolds = 0;
            particleReynolds = Math.Sqrt(TranslationalVelocity[0][0] * TranslationalVelocity[0][0] + TranslationalVelocity[0][1] * TranslationalVelocity[0][1]) * 2 * length_P * particleDensity / mu_Fluid;
            Console.WriteLine("Particle Reynolds number:  " + particleReynolds);
            return particleReynolds;
        }
    }
}


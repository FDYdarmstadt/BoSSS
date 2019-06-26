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
using BoSSS.Foundation;
using System.Linq;
using System.Collections.Generic;
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
        public Particle_Ellipsoid(double[] startPos = null, double startAngl = 0) : base(2, startPos, startAngl) {
            
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
        /// %
        /// </summary>
        protected override double AverageDistance {
            get {
                return 0.5 * (length_P + thickness_P);
            }
        }

        protected override double Area_P {
            get {
                double a = length_P * thickness_P * Math.PI;
                if (a <= 0.0 || double.IsNaN(a) || double.IsInfinity(a))
                    throw new ArithmeticException("Ellipsoid volume/area is " + a);
                return a;
            }
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
        //override public void UpdateLevelSetFunction() {
        //    double alpha = -(Angle[0]);
        //    Phi_P = delegate (double[] X, double t) {
        //        var r = -((((X[0] - Position[0][0]) * Math.Cos(alpha) - (X[1] - Position[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - Position[0][0]) * Math.Sin(alpha) + (X[1] - Position[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + 1.0;
        //        if (double.IsNaN(r) || double.IsInfinity(r))
        //            throw new ArithmeticException();
        //        return r;
        //    };
        //}
        public override double Phi_P(double[] X) {
            double alpha = -(Angle[0]);
            var r = -((((X[0] - Position[0][0]) * Math.Cos(alpha) - (X[1] - Position[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) 
                    -(((X[0] - Position[0][0]) * Math.Sin(alpha) + (X[1] - Position[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) 
                    + 1.0;
            if (double.IsNaN(r) || double.IsInfinity(r))
                throw new ArithmeticException();
            return r;
        }
        override public CellMask CutCells_P(LevelSetTracker LsTrk) {
            // tolerance is very important
            var radiusTolerance = Math.Max(length_P, thickness_P) + LsTrk.GridDat.Cells.h_minGlobal;// +2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());

            CellMask cellCollection;
            CellMask cells = null;
            double alpha = -(Angle[0]);
            cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - Position[0][0]) * Math.Cos(alpha) - (X[1] - Position[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - Position[0][0]) * Math.Sin(alpha) + (X[1] - Position[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0);

            CellMask allCutCells = LsTrk.Regions.GetCutCellMask();
            cellCollection = cells.Intersect(allCutCells);
            return cellCollection;
        }
        public override bool Contains(double[] point, LevelSetTracker LsTrk, bool WithoutTolerance = false) {
            // only for squared cells
            WithoutTolerance = true;
            double radiusTolerance = !WithoutTolerance ? 1.0 + 2.0 * Math.Sqrt(2 * LsTrk.GridDat.Cells.h_minGlobal.Pow2()) : 1;
            double Ellipse = ((point[0] - Position[0][0]) * Math.Cos(Angle[0] + (point[1] - Position[0][1]) * Math.Sin(Angle[0]))).Pow2() / length_P.Pow2() + (-(point[0] - Position[0][0]) * Math.Sin(Angle[0]) + (point[1] - Position[0][1]) * Math.Cos(Angle[0])).Pow2() / thickness_P.Pow2();
            if (Ellipse < radiusTolerance)
            {
                return true;
            }
            else
                return false;
        }

        override public double ComputeParticleRe(double mu_Fluid) {
            double particleReynolds = 0;
            particleReynolds = Math.Sqrt(TranslationalVelocity[0][0] * TranslationalVelocity[0][0] + TranslationalVelocity[0][1] * TranslationalVelocity[0][1]) * 2 * length_P * 1 / mu_Fluid;
            Console.WriteLine("Particle Reynolds number:  " + particleReynolds);
            return particleReynolds;
        }

        override public MultidimensionalArray GetSurfacePoints(LevelSetTracker lsTrk, double[] PositionS, double AngleS)
        {
            int SpatialDim = lsTrk.GridDat.SpatialDimension;
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            double hMin = lsTrk.GridDat.iGeomCells.h_min.Min();
            int NoOfSurfacePoints = Convert.ToInt32(5 * Circumference_P / hMin);
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSubParticles(), NoOfSurfacePoints, SpatialDim);
            double[] InfinitisemalAngle = GenericBlas.Linspace(0, Math.PI * 2, NoOfSurfacePoints + 1);
            if (Math.Abs(10 * Circumference_P / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");
            for (int j = 0; j < NoOfSurfacePoints; j++)
            {
                double temp0 = Math.Cos(InfinitisemalAngle[j]) * length_P;
                double temp1 = Math.Sin(InfinitisemalAngle[j]) * thickness_P;
                SurfacePoints[0, j, 0] = (temp0 * Math.Cos(AngleS) - temp1 * Math.Sin(AngleS)) + PositionS[0];
                SurfacePoints[0, j, 1] = (temp0 * Math.Sin(AngleS) + temp1 * Math.Cos(AngleS)) + PositionS[1];

            }
            return SurfacePoints;
        }

        override public void GetSupportPoint(int SpatialDim, double[] Vector, double[] Position, double Angle, out double[] SupportPoint)
        {
            SupportPoint = new double[SpatialDim];
            if (double.IsNaN(Vector[0]) || double.IsNaN(Vector[1]))
                throw new ArithmeticException("Error trying to calculate VectorVectorVector Value:  " + Vector[0] + " VectorVectorVector " + Vector[1]);
            double[,] B = new double[2, 2];
            B[0, 0] = length_P * Math.Cos(Angle);
            if (double.IsNaN(B[0, 0]))
                throw new ArithmeticException("Error trying to calculate Angle Value:  " + B[0, 0] + " length_P " + length_P + " Math.Cos(Angle): " + Math.Cos(Angle) + " Angle " + Angle);
            B[0, 1] = -thickness_P * Math.Sin(Angle);
            B[1, 0] = length_P * Math.Sin(Angle);
            B[1, 1] = thickness_P * Math.Cos(Angle);
            if (double.IsNaN(Angle))
                throw new ArithmeticException("Error trying to calculate Angle Value:  " + Angle + " Angle " + Angle);
            double[,] BT = B.CloneAs();
            BT[0, 1] = B[1, 0];
            BT[1, 0] = B[0, 1];
            double[] temp = new double[2];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    temp[i] += BT[i, j] * Vector[j];
                    if (double.IsNaN(temp[0]) || double.IsNaN(temp[1]))
                        throw new ArithmeticException("Error trying to calculate temp Value:  " + temp[i] + "VectorVectorVector Value: " + Vector[j] + " BT[i, j] " + BT[i, j] + " afuduiof" + i + j);
                }
            }
            double BetragTemp = Math.Sqrt(temp[0].Pow2() + temp[1].Pow2());
            if (double.IsNaN(temp[0]) || double.IsNaN(temp[1]))
                throw new ArithmeticException("Error trying to calculate temp Value:  " + temp[0] + " temp " + temp[1] + "VectorVectorVector Value: " + Vector[0] + " VectorVectorVector " + Vector[1]);
            for (int i = 0; i < 2; i++)
            {
                if(BetragTemp != 0)
                    temp[i] = temp[i] / BetragTemp;
            }
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    SupportPoint[i] += B[i, j] * temp[j];
                    if (double.IsNaN(SupportPoint[i]) || double.IsInfinity(SupportPoint[i]))
                        throw new ArithmeticException("Error trying to calculate SupportPoint Value:  " + SupportPoint[i] + "temp Value: " + temp[j] + " B[i, j] " + B[i, j] + "BetragTemp" + BetragTemp + "ij" +i +j);
                }
                SupportPoint[i] += Position[i];
            }
        }
        override public double[] GetLengthScales()
        {
            return new double[] { length_P, thickness_P };
        }
    }
}


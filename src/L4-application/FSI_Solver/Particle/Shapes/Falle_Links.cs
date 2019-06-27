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
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Grid.Classic;
using System.Text;
using BoSSS.Platform;
using ilPSP.Tracing;
using System.Diagnostics;
using ilPSP.Utils;
using System.IO;
using System.Collections;
using BoSSS.Platform.Utils.Geom;

namespace BoSSS.Application.FSI_Solver
{
    [DataContract]
    [Serializable]
    public class Particle_Falle_Links : Particle
    {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Falle_Links() : base()
        {

        }

        public Particle_Falle_Links(double[] startPos = null, double startAngl = 0) : base(2, startPos, startAngl)
        {


        }

        /// <summary>
        /// Radius of the particle. Not necessary for particles defined by their length and thickness
        /// </summary>
        [DataMember]
        public double width_P;

        internal override int NoOfSubParticles()
        {
            return 2;
        }

        /// <summary>
        /// %
        /// </summary>
        protected override double AverageDistance
        {
            get
            {
                return width_P;
            }
        }

        protected override double Area_P
        {
            get
            {
                return (7 * width_P * width_P) / 8;
            }
        }
        protected override double Circumference_P
        {
            get
            {
                return width_P * 5;
            }
        }
        override public double MomentOfInertia_P
        {
            get
            {
                return Math.Pow(width_P, 4) * 0.13785958;
            }
        }
        //override public void UpdateLevelSetFunction()
        //{
        //    double alpha = -(Angle[0]);
        //    Phi_P = (X, t) => -(X[0] - Position[0][0]).Pow2() + -(X[1] - Position[0][1]).Pow2() + radius_P.Pow2();
        //}
        public override double Phi_P(double[] X)
        {
            double alpha = -(Angle[0]);
            double r;
            // Rechteck:
            //        r = Math.Max(X[0] - Position[0][0]  - Width_P,  Position[0][0] - Width_P - X[0]);
            //        r = Math.Max(r, X[1] - Position[0][1] - Width_P);
            //        r = Math.Max(r,  Position[0][1] - 0.5*Width_P - X[1]);

            // parallelogram
            //          r = Math.Abs(Position[0][1] + 0.5 * width_P - X[1]);
            //          r = Math.Max(r, Math.Abs(-X[1] - 0.5 * X[0] + Position[0][1] + width_P) + Math.Abs(X[1] - Position[0][1]));
            //          r = Math.Max(r, Math.Abs(Position[0][0] - 0.5 * width_P - X[0]));

            // Falle_Links:
                      r = Math.Abs(Position[0][1] - X[1]);
                      r = Math.Max(r, Math.Abs(-X[1] + 0.5*X[0] + Position[0][1] - Position[0][0] - width_P) - Math.Abs(X[1] - Position[0][1]));
                      r = Math.Max(r, Math.Abs(Position[0][0] - X[0] + 1.5*width_P)); 
                      r = r - 3.5*width_P;

            r = -r;
            return r;
        }

        override public CellMask CutCells_P(LevelSetTracker LsTrk)
        {
            // tolerance is very important
            var radiusTolerance = width_P + LsTrk.GridDat.Cells.h_minGlobal;// +2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());

            CellMask cellCollection;
            CellMask cells = null;
            double alpha = -(Angle[0]);
            cells = CellMask.GetCellMask(LsTrk.GridDat, X => (-(X[0] - Position[0][0]).Pow2() + -(X[1] - Position[0][1]).Pow2() + radiusTolerance.Pow2()) > 0);

            CellMask allCutCells = LsTrk.Regions.GetCutCellMask();
            cellCollection = cells.Intersect(allCutCells);
            return cellCollection;
        }

        public override bool Contains(double[] point, LevelSetTracker LsTrk, bool WithoutTolerance = false) {
           
            // only for squared cells
            double radiusTolerance = width_P + 2.0 * Math.Sqrt(2 * LsTrk.GridDat.Cells.h_minGlobal.Pow2());
            double length_P = 1;
            double thickness_P = 0.2;
            double test = -((((point[0] - Position[0][0]) * Math.Cos(Angle[0]) - (point[1] - Position[0][1]) * Math.Sin(Angle[0])).Pow2()) / length_P.Pow2()) + -(((point[0] - Position[0][0]) * Math.Sin(Angle[0]) + (point[1] - Position[0][1]) * Math.Cos(Angle[0])).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2();
            var distance = point.L2Distance(Position[0]);
            if (distance < (radiusTolerance))
            {
                return true;
            }
            return false;
        }
        override public double ComputeParticleRe(double mu_Fluid)
        {
            double particleReynolds = 0;
            particleReynolds = Math.Sqrt(TranslationalVelocity[0][0] * TranslationalVelocity[0][0] + TranslationalVelocity[0][1] * TranslationalVelocity[0][1]) * 2 * width_P * particleDensity / mu_Fluid;
            Console.WriteLine("Particle Reynolds number:  " + particleReynolds);
            return particleReynolds;
        }

        override public MultidimensionalArray GetSurfacePoints(LevelSetTracker lsTrk, double[] Position, double Angle)
        {
            int SpatialDim = lsTrk.GridDat.SpatialDimension;
            if (SpatialDim != 2)
                throw new NotImplementedException("Only two dimensions are supported at the moment");

            double hMin = lsTrk.GridDat.iGeomCells.h_min.Min();
            int NoOfSurfacePoints = Convert.ToInt32(20 * Circumference_P / hMin) + 1;
            MultidimensionalArray SurfacePoints = MultidimensionalArray.Create(NoOfSubParticles(), NoOfSurfacePoints, 2);
            double[] InfinitisemalAngle = GenericBlas.Linspace(-Math.PI / 4, 5 * Math.PI / 4, NoOfSurfacePoints + 1);
            double[] InfinitisemalLength = GenericBlas.Linspace(0, width_P / 4, NoOfSurfacePoints + 1);
            if (Math.Abs(10 * Circumference_P / hMin + 1) >= int.MaxValue)
                throw new ArithmeticException("Error trying to calculate the number of surface points, overflow");

            for (int k = 0; k < NoOfSurfacePoints; k++)
            {
                SurfacePoints[0, k, 0] = Position[0] - width_P / 2 + InfinitisemalLength[k];
                SurfacePoints[0, k, 1] = Position[1] - width_P / 2 - 1.5 * SurfacePoints[k, 0] + width_P / 2;
            }

            for (int j = 0; j < NoOfSurfacePoints; j++)
            {
                SurfacePoints[0, j, 0] = Math.Sign(Math.Cos(InfinitisemalAngle[j])) * width_P * 7 + Position[0] + 7 * width_P / 4;
                SurfacePoints[0, j, 1] = Math.Sign(Math.Sin(InfinitisemalAngle[j])) * width_P * 7 + Position[1] + 7 * width_P / 2;
            }
            for (int j = 0; j < NoOfSurfacePoints; j++)
            {
                SurfacePoints[1, j, 0] = -Math.Sign(Math.Cos(InfinitisemalAngle[j])) * width_P * 2.5 + Position[0];
                SurfacePoints[1, j, 1] = -Math.Sign(Math.Sin(InfinitisemalAngle[j])) * width_P * 2.5 + Position[1];
            }
            return SurfacePoints;
        }


        override public double[] GetLengthScales()
        {
            return new double[] { width_P, width_P };
        }
    }
}

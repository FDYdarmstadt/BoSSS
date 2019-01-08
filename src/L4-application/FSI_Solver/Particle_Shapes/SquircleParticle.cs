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

namespace BoSSS.Application.FSI_Solver
{
    [DataContract]
    [Serializable]
    public class EllipsoidParticle : Particle
    {
        public EllipsoidParticle(int Dim, int HistoryLength, double[] startPos = null, double startAngl = 0, ParticleShape shape = ParticleShape.spherical) : base(Dim, HistoryLength, startPos, startAngl, shape)
        {
        }
        override public double active_stress_P
        {
            get
            {
                //Approximation formula for circumference according to Ramanujan
                double circumference;
                circumference = Math.PI * ((length_P + thickness_P) + (3 * (length_P - thickness_P).Pow2()) / (10 * (length_P + thickness_P) + Math.Sqrt(length_P.Pow2() + 14 * length_P * thickness_P + thickness_P.Pow2())));
                return 0.5 * circumference * stress_magnitude_P;
            }
        }
        override public double area_P
        {
            get
            {
                return length_P * thickness_P * Math.PI * radius_P;
            }

        }
        override public double MomentOfInertia_P
        {
            get
            {
                return (1 / 4.0) * (mass_P * (length_P * length_P + thickness_P * thickness_P) * radius_P * radius_P);
            }
        }
        override public CellMask cutCells_P(LevelSetTracker LsTrk)
        {
            // tolerance is very important
            var radiusTolerance = radius_P + LsTrk.GridDat.Cells.h_minGlobal;// +2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());

            CellMask cellCollection;
            CellMask cells = null;
            double alpha = -(currentIterAng_P[0]);
            cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0);
            
            CellMask allCutCells = LsTrk.Regions.GetCutCellMask();
            cellCollection = cells.Intersect(allCutCells);
            return cellCollection;
        }
        override public bool Contains(double[] point, LevelSetTracker LsTrk)
        {
            // only for squared cells
            double radiusTolerance = radius_P + 2.0 * Math.Sqrt(2 * LsTrk.GridDat.Cells.h_minGlobal.Pow2());
            if (-((((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow2()) / length_P.Pow2()) + -(((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0)
            {
                return true;
            }     
            return false;
        }
    }

}


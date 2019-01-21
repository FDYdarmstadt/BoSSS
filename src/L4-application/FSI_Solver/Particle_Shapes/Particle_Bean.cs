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
    public class Particle_Bean : Particle
    {
        public Particle_Bean(int Dim, int HistoryLength, double[] startPos = null, double startAngl = 0) : base(Dim, HistoryLength, startPos, startAngl)
        {
            #region Particle history
            // =============================   
            for (int i = 0; i < HistoryLength; i++)
            {
                currentIterPos_P.Add(new double[Dim]);
                currentIterAng_P.Add(new double());
                currentIterVel_P.Add(new double[Dim]);
                currentIterRot_P.Add(new double());
                currentIterForces_P.Add(new double[Dim]);
                temporalForces_P.Add(new double[Dim]);
                currentIterTorque_P.Add(new double());
                temporalTorque_P.Add(new double());
            }
            for (int i = 0; i < 4; i++)
            {
                currentTimePos_P.Add(new double[Dim]);
                currentTimeAng_P.Add(new double());
                currentTimeVel_P.Add(new double[Dim]);
                currentTimeRot_P.Add(new double());
                currentTimeForces_P.Add(new double[Dim]);
                currentTimeTorque_P.Add(new double());
            }
            #endregion

            #region Initial values
            // ============================= 
            if (startPos == null)
            {
                if (Dim == 2)
                {
                    startPos = new double[] { 0.0, 0.0 };
                }
                else
                {
                    startPos = new double[] { 0.0, 0.0, 0.0 };
                }
            }
            currentTimePos_P[0] = startPos;
            currentTimePos_P[1] = startPos;
            //From degree to radiant
            currentTimeAng_P[0] = startAngl * 2 * Math.PI / 360;
            currentTimeAng_P[1] = startAngl * 2 * Math.PI / 360;
            //currentIterVel_P[0][0] = 2e-8;

            UpdateLevelSetFunction();
            #endregion
        }
        override public double active_stress_P
        {
            get
            {
                return 2 * Math.PI * radius_P * stress_magnitude_P;
            }
        }
        override public double Area_P
        {
            get
            {
                // not correct area
                return Math.PI * radius_P * radius_P;
            }
        }
        override public double MomentOfInertia_P
        {
            get
            {
                // not correct moment of inertia
                return (1 / 2.0) * (Mass_P * radius_P * radius_P);
            }
        }
        override public void UpdateLevelSetFunction()
        {
            double a = 3.0 * radius_P.Pow2();
            double b = 1.0 * radius_P.Pow2();
            double alpha = -(currentIterAng_P[0]);
            phi_P = (X, t) => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(3) - b * ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2());
        }
        override public CellMask cutCells_P(LevelSetTracker LsTrk)
        {
            // tolerance is very important
            var radiusTolerance = radius_P + LsTrk.GridDat.Cells.h_minGlobal;// +2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());

            CellMask cellCollection;
            CellMask cells = null;
            double alpha = -(currentIterAng_P[0]);
            double a = 4.0 * radiusTolerance.Pow2();
            double b = 1.0 * radiusTolerance.Pow2();
            cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(2) + ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow(2)).Pow2() - a * ((X[0] - currentIterPos_P[0][0]) * Math.Cos(alpha) - (X[1] - currentIterPos_P[0][1]) * Math.Sin(alpha)).Pow(3) - b * ((X[0] - currentIterPos_P[0][0]) * Math.Sin(alpha) + (X[1] - currentIterPos_P[0][1]) * Math.Cos(alpha)).Pow2()) > 0);

            CellMask allCutCells = LsTrk.Regions.GetCutCellMask();
            cellCollection = cells.Intersect(allCutCells);
            return cellCollection;
        }
        override public bool Contains(double[] point, LevelSetTracker LsTrk)
        {
            // only for squared cells
            double radiusTolerance = radius_P + 2.0 * Math.Sqrt(2 * LsTrk.GridDat.Cells.h_minGlobal.Pow2());
            double a = 4.0 * radiusTolerance.Pow2();
            double b = 1.0 * radiusTolerance.Pow2();
            if (-((((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow(2) + ((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow(2)).Pow2() - a * ((point[0] - currentIterPos_P[0][0]) * Math.Cos(currentIterAng_P[0]) - (point[1] - currentIterPos_P[0][1]) * Math.Sin(currentIterAng_P[0])).Pow(3) - b * ((point[0] - currentIterPos_P[0][0]) * Math.Sin(currentIterAng_P[0]) + (point[1] - currentIterPos_P[0][1]) * Math.Cos(currentIterAng_P[0])).Pow2()) > 0)
            {
                return true;
            }
            return false;
        }
        override public double ComputeParticleRe(double mu_Fluid)
        {
            double particleReynolds = 0;
            particleReynolds = Math.Sqrt(currentIterVel_P[0][0] * currentIterVel_P[0][0] + currentIterVel_P[0][1] * currentIterVel_P[0][1]) * 2 * 4.0 * rho_P / mu_Fluid;
            Console.WriteLine("Particle Reynolds number:  " + particleReynolds);
            return particleReynolds;
        }
    }
}


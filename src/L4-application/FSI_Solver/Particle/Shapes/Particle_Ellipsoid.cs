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

namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class Particle_Ellipsoid : Particle {
        /// <summary>
        /// Empty constructor used during de-serialization
        /// </summary>
        private Particle_Ellipsoid() : base() {

        }
        public Particle_Ellipsoid(int HistoryLength, double[] startPos = null, double startAngl = 0) : base(2, HistoryLength, startPos, startAngl) {
            #region Particle history
            // ============================= 

            int Dim = 2; // spatial dimension of ellipse is always 2

            for (int i = 0; i < HistoryLength; i++) {
                positionAtIteration.Add(new double[Dim]);
                angleAtIteration.Add(new double());
                transVelocityAtIteration.Add(new double[Dim]);
                rotationalVelocityAtIteration.Add(new double());
                hydrodynForcesAtIteration.Add(new double[Dim]);
                hydrodynTorqueAtIteration.Add(new double());
            }
            for (int i = 0; i < 4; i++) {
                positionAtTimestep.Add(new double[Dim]);
                angleAtTimestep.Add(new double());
                transVelocityAtTimestep.Add(new double[Dim]);
                rotationalVelocityAtTimestep.Add(new double());
                hydrodynForcesAtTimestep.Add(new double[Dim]);
                hydrodynTorqueAtTimestep.Add(new double());
            }
            #endregion

            #region Initial values
            // ============================= 
            if (startPos == null) {
                startPos = new double[Dim];
            }
            positionAtTimestep[0] = startPos;
            positionAtTimestep[1] = startPos;
            //From degree to radiant
            angleAtTimestep[0] = startAngl * 2 * Math.PI / 360;
            angleAtTimestep[1] = startAngl * 2 * Math.PI / 360;
            //transVelocityAtIteration[0][0] = 2e-8;

            UpdateLevelSetFunction();
            #endregion
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
        protected override double averageDistance {
            get {
                return 0.5*(length_P + thickness_P);
            }
        }



        override public double Area_P {
            get {
                double a = length_P * thickness_P * Math.PI;
                if (a <= 0.0 || double.IsNaN(a) || double.IsInfinity(a))
                    throw new ArithmeticException("Ellipsoid volume/area is " + a);
                return a;
            }

        }
        public override double Circumference_P {
            get {
                return Math.PI * ((length_P + thickness_P) + (3 * (length_P - thickness_P).Pow2()) / (10 * (length_P + thickness_P) + Math.Sqrt(length_P.Pow2() + 14 * length_P * thickness_P + thickness_P.Pow2())));
            }
        }
        override public double MomentOfInertia_P {
            get {
                return (1 / 4.0) * (Mass_P * (length_P * length_P + thickness_P * thickness_P));
            }
        }
        override public void UpdateLevelSetFunction() {
            double alpha = -(angleAtIteration[0]);
            phi_P = delegate (double[] X, double t) {
                var r = -((((X[0] - positionAtIteration[0][0]) * Math.Cos(alpha) - (X[1] - positionAtIteration[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - positionAtIteration[0][0]) * Math.Sin(alpha) + (X[1] - positionAtIteration[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + 1.0;
                if (double.IsNaN(r) || double.IsInfinity(r))
                    throw new ArithmeticException();
                return r;
            };
        }
        override public CellMask cutCells_P(LevelSetTracker LsTrk) {
            // tolerance is very important
            var radiusTolerance = Math.Max(length_P, thickness_P) + LsTrk.GridDat.Cells.h_minGlobal;// +2.0*Math.Sqrt(2*LsTrk.GridDat.Cells.h_minGlobal.Pow2());

            CellMask cellCollection;
            CellMask cells = null;
            double alpha = -(angleAtIteration[0]);
            cells = CellMask.GetCellMask(LsTrk.GridDat, X => -((((X[0] - positionAtIteration[0][0]) * Math.Cos(alpha) - (X[1] - positionAtIteration[0][1]) * Math.Sin(alpha)).Pow2()) / length_P.Pow2()) + -(((X[0] - positionAtIteration[0][0]) * Math.Sin(alpha) + (X[1] - positionAtIteration[0][1]) * Math.Cos(alpha)).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0);

            CellMask allCutCells = LsTrk.Regions.GetCutCellMask();
            cellCollection = cells.Intersect(allCutCells);
            return cellCollection;
        }
        override public bool Contains(double[] point, LevelSetTracker LsTrk) {
            // only for squared cells
            double radiusTolerance = 1.0 + 2.0 * Math.Sqrt(2 * LsTrk.GridDat.Cells.h_minGlobal.Pow2());
            if (-((((point[0] - positionAtIteration[0][0]) * Math.Cos(angleAtIteration[0]) - (point[1] - positionAtIteration[0][1]) * Math.Sin(angleAtIteration[0])).Pow2()) / length_P.Pow2()) + -(((point[0] - positionAtIteration[0][0]) * Math.Sin(angleAtIteration[0]) + (point[1] - positionAtIteration[0][1]) * Math.Cos(angleAtIteration[0])).Pow2() / thickness_P.Pow2()) + radiusTolerance.Pow2() > 0) {
                return true;
            }
            return false;
        }

        override public double ComputeParticleRe(double mu_Fluid) {
            double particleReynolds = 0;
            particleReynolds = Math.Sqrt(transVelocityAtIteration[0][0] * transVelocityAtIteration[0][0] + transVelocityAtIteration[0][1] * transVelocityAtIteration[0][1]) * 2 * length_P * 1 / mu_Fluid;
            Console.WriteLine("Particle Reynolds number:  " + particleReynolds);
            return particleReynolds;
        }
    }

}


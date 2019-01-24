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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Solution.Utils;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP;

namespace BoSSS.Solution.XheatCommon {


    public class EvaporationAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        /// <param name="_sigma">surface-tension constant</param>
        public EvaporationAtLevelSet(LevelSetTracker LsTrk, double _hVapA, double _Rint, double _Tsat, double _rho, double _sigma) {
            m_LsTrk = LsTrk;
            //this.mEvap = _mEvap;
            //this.hVap = _hVap;
            this.hVapA = _hVapA;
            this.Rint = _Rint;
            this.Tsat = _Tsat;
            this.rho = _rho;
            this.sigma = _sigma;

        }

        //double mEvap;
        //double hVap;

        double hVapA;   // for the identification of the liquid phase
        double Rint;
        double Tsat;
        double rho;     // density of liquid phase 
        double sigma;


        private double ComputeHeatFlux(double T_A, double T_B, double curv, double p_disp) {

            if(hVapA == 0.0)
                return 0.0;

            double pc0 = 0.0; // sigma*curv + p_disp;   // augmented capillary pressure (without nonlinear evaporative masss part)

            double TintMin = 0.0;
            double hVap = 0.0;
            double qEvap = 0.0;
            if(hVapA > 0) {
                hVap = hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rho)));
                if(T_A > TintMin)
                    qEvap = -(T_A - TintMin) / Rint;
            } else if(hVapA < 0) {
                hVap = -hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rho)));
                if(T_B > TintMin)
                    qEvap = (T_B - TintMin) / Rint;
            }

            return qEvap;
        }


        public double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Debug.Assert(cp.ParamsPos[1] == cp.ParamsNeg[1], "curvature must be continuous across interface");
            Debug.Assert(cp.ParamsPos[2] == cp.ParamsNeg[2], "disjoining pressure must be continuous across interface");

            double qEvap = ComputeHeatFlux(cp.ParamsNeg[0], cp.ParamsPos[0], cp.ParamsNeg[1], cp.ParamsNeg[2]);

            //double massFluxEvap = mEvap * hVap;

            double FlxNeg = -0.5 * qEvap;
            double FlxPos = +0.5 * qEvap;

            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            return FlxNeg * vA - FlxPos * vB;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }


        public IList<string> ParameterOrdering {
            get {
                return new string[] { "Temperature0", "Curvature", "DisjoiningPressure" };
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }

    }

}

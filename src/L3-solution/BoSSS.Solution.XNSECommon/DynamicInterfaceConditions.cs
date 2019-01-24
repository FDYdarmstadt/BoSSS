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


namespace BoSSS.Solution.XNSECommon.Operator.DynamicInterfaceConditions {


    public class PrescribedMassFlux : ILevelSetForm {


        LevelSetTracker m_LsTrk;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_d">spatial direction</param>
        /// <param name="_D">spatial dimension</param>
        /// <param name="LsTrk"></param>
        public PrescribedMassFlux(int _d, int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, double _kA, double _kB, double _hVapA, double _Rint, double _Tsat, double _sigma) {
            m_LsTrk = LsTrk;
            if(_d >= _D)
                throw new ArgumentOutOfRangeException();
            this.m_D = _D;
            this.m_d = _d;

            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            //this.M = _M;
            this.kA = _kA;
            this.kB = _kB;
            this.hVapA = _hVapA;
            this.Rint = _Rint;
            //this.TintMin = _TintMin;
            this.Tsat = _Tsat;
            this.sigma = _sigma;
        }

        int m_D;
        int m_d;

        double rhoA;
        double rhoB;

        double kA;
        double kB;
        double hVapA;   // for the identification of the liquid phase
        double Rint;
        //double TintMin;
        double Tsat;
        double sigma;

        //double M;


        //private double ComputeEvaporationMass(double[] GradT_A, double[] GradT_B, double[] n) {

        //    double mEvap = 0.0;

        //    // for testing purposes
        //    double prescribedVolumeFlux = 0.1;
        //    if(hVapA > 0) {
        //        mEvap = -rhoA * prescribedVolumeFlux;
        //    } else {
        //        mEvap = rhoB * prescribedVolumeFlux;
        //    }

        //    // TODO 

        //    return mEvap;
        //}

        private double ComputeEvaporationMass(double T_A, double T_B, double curv, double p_disp) {

            if(hVapA == 0.0)
                return 0.0;

            double pc0 = 0.0; // sigma*curv + p_disp;   // augmented capillary pressure (without nonlinear evaporative masss part)

            double TintMin = 0.0;
            double hVap = 0.0;
            double qEvap = 0.0;
            if(hVapA > 0) {
                hVap = hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rhoA)));
                if(T_A > TintMin)
                    qEvap = -(T_A - TintMin) / Rint;
            } else if(hVapA < 0) {
                hVap = -hVapA;
                TintMin = Tsat * (1 + (pc0 / (hVap * rhoB)));
                if(T_B > TintMin)
                    qEvap = (T_B - TintMin) / Rint;
            }

            return qEvap / hVap;
        }


        public double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Normal = cp.n;

            Debug.Assert(cp.ParamsPos[1] == cp.ParamsNeg[1], "curvature must be continuous across interface");
            Debug.Assert(cp.ParamsPos[2] == cp.ParamsNeg[2], "disjoining pressure must be continuous across interface");

            double M = ComputeEvaporationMass(cp.ParamsNeg[0], cp.ParamsPos[0], cp.ParamsNeg[1], cp.ParamsNeg[2]);
            double massFlux = -M.Pow2() * ((1/rhoA) - (1/rhoB)) * Normal[m_d];

            double p_disp = cp.ParamsNeg[1];

            // augmented capillary pressure
            double acp_jump = 0.0;
            if(!double.IsNaN(p_disp))
                acp_jump = massFlux + p_disp;
            else
                acp_jump = massFlux;


            double FlxNeg = -0.5 * acp_jump;
            double FlxPos = +0.5 * acp_jump;


            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            return FlxNeg * vA - FlxPos * vB;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] {  };
            }
        }


        public IList<string> ParameterOrdering {
            get {
                return new string[] { VariableNames.Temperature, "Curvature", "DisjoiningPressure" }; //{ "GradTempX", "GradTempY", "GradTempZ" }.GetSubVector(0, m_D);
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

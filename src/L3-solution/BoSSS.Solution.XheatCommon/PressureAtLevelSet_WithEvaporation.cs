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
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;
using System.Collections;

namespace BoSSS.Solution.XheatCommon {

    /// <summary>
    /// 
    /// </summary>
    public class GeneralizedPressureFormAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public GeneralizedPressureFormAtLevelSet(int _d, LevelSetTracker LsTrk, double _pSat, double _hVapA) {
            m_d = _d;
            m_LsTrk = LsTrk;

            this.pSat = _pSat;
            this.hVapA = _hVapA;
        }

        int m_d;

        double pSat;
        double hVapA;


        public double LevelSetForm(ref CommonParamsLs inp, double[] pA, double[] pB, double[,] Grad_pA, double[,] Grad_pB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double acc = 0.0;
            if (hVapA > 0.0) {
                acc += (0 - vA) * inp.n[m_d] * pSat;
                acc += (vB - 0) * inp.n[m_d] * pB[0];
            } else {
                acc += (0 - vA) * inp.n[m_d] * pA[0];
                acc += (vB - 0) * inp.n[m_d] * pSat;
            }
            //return (vB - vA) * inp.n[m_d] * pSat;
            return -acc;
        }



        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Pressure };
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
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get { return null; }
        }


    }

}

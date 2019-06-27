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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Solution.XNSECommon.Operator.Pressure {
    
    /// <summary>
    /// 
    /// </summary>
    public class PressureFormAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public PressureFormAtLevelSet(int _d, int _D, LevelSetTracker LsTrk, bool _weighted = false, double _wA = 1.0, double _wB = 1.0) {
            m_d = _d;
            m_D = _D;
            m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentException();

            weighted = _weighted;
            wA = _wA;
            wB = _wB;
        }

        int m_d;
        int m_D;

        bool weighted;
        double wA;
        double wB;

        public double LevelSetForm(ref CommonParamsLs inp, double[] pA, double[] pB, double[,] Grad_pA, double[,] Grad_pB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            if (!weighted) {
                return -(vB - vA) * inp.n[m_d] * 0.5 * (pB[0] + pA[0]);
            } else {
                return -(vB - vA) * inp.n[m_d] * (wA * pB[0] + wB * pA[0]) / (wA + wB);
            }
            
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
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ParameterOrdering {
            get { return null; }
        }
    }

    /// <summary>
    /// 
    /// </summary>
    public class GeneralizedPressureFormAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public GeneralizedPressureFormAtLevelSet(int _d, int _D, LevelSetTracker LsTrk, double _pSat, double _hVapA) {
            m_d = _d;
            m_D = _D;
            m_LsTrk = LsTrk;
            if (_d >= _D)
                throw new ArgumentException();

            this.pSat = _pSat;
            this.hVapA = _hVapA;
        }

        int m_d;
        int m_D;

        double pSat;
        double hVapA;


        public double LevelSetForm(ref CommonParamsLs inp, double[] pA, double[] pB, double[,] Grad_pA, double[,] Grad_pB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            if (hVapA > 0.0)
                return (0 - vA) * inp.n[m_d] * pSat;
            else
                return (vB - 0) * inp.n[m_d] * pSat;

            //return (vB - vA) * inp.n[m_d] * pSat;
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
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ParameterOrdering {
            get { return null; }
        }
    }

}

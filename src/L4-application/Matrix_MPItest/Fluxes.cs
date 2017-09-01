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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;

namespace BoSSS.Application.Matrix_MPItest {

    /// <summary>
    /// fluss fuer du/dx; (Ableitung nach 1. Raumrichtung), bulk-Phase;
    /// </summary>
    class DxFlux : LinearFlux {

        public DxFlux(string varname, double factor) {
            m_varname = varname;
            m_factor = factor;
        }

        string m_varname;
        double m_factor;

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { m_varname };
            }
        }

        protected override double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            return Uin[0]*inp.Normale[0]*m_factor;
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return 0.5*(Uin[0] + Uout[0])*inp.Normale[0]*m_factor;
        }

        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            output[0] = U[0]*m_factor;
        }
    }

    /// <summary>
    /// Fluss fuer du/dx; (Ableitung nach 1. Raumrichtung), common parts for both level-sets;
    /// </summary>
    class LevSetFlx : ILevelSetComponent {

        protected LevelSetTracker m_LsTrk;

        public LevSetFlx(LevelSetTracker _LsTrk, string varname, double factor) {
            m_LsTrk = _LsTrk;
            m_varname = varname;
            m_factor = factor;
        }

        double m_factor;
        string m_varname;

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { m_varname };
            }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public double LevelSetForm(ref CommonParamsLs inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double Flx = 0.5*(U_Pos[0] + U_Neg[0])*inp.n[0];
            return Flx * vA - Flx * vB*m_factor;
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { 
                return m_LsTrk.GetSpeciesId("A"); 
            }
        }

        public SpeciesId NegativeSpecies {
            get { 
                return m_LsTrk.GetSpeciesId("B"); 
            }
        }
    }

    
}

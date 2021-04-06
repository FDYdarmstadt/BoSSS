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

using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;

namespace BoSSS.Solution.RheologyCommon {
    /// <summary>
    /// Volume integral of identity part of constitutive equations.
    /// </summary>
    public class DiffusionInBulk : ConstitutiveEqns_Diffusion, ISpeciesFilter {

        SpeciesId m_spcId;
        string variable;

        /// <summary>
        /// Initialize Diffusion for artificial diffusion
        /// </summary>
        public DiffusionInBulk(double orderOfOne, string _variable, string spcName, SpeciesId spcId) 
            : base(orderOfOne, _variable) {
            this.variable = _variable;
            this.m_spcId = spcId;
            this.ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }
}

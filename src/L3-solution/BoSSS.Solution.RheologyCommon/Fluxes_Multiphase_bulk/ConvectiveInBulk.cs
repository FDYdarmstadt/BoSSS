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

using BoSSS.Foundation.XDG;
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Solution.RheologyCommon {

    /// <summary>
    /// Convective part of constitutive equations in bulk for multiphase.
    /// </summary>
    public class ConvectiveInBulk : ConstitutiveEqns_Convective, ISpeciesFilter {

        /// <summary>
        /// Initialize Convection
        /// </summary>
        public ConvectiveInBulk(int _Component, IncompressibleMultiphaseBoundaryCondMap _BcMap, double Weissenberg, double alpha, string spcName, SpeciesId spcId) : base(_Component, _BcMap, Weissenberg, true, alpha = 1.0) {
            this.ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }
}

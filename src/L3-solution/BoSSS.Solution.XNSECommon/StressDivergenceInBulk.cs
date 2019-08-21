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

namespace BoSSS.Solution.RheologyCommon {
    /// <summary>
    /// Volume integral of identity part of constitutive equations.
    /// </summary>
    public class StressDivergenceInBulk : StressDivergence_Cockburn, ISpeciesFilter {

        SpeciesId m_spcId;
        int Component;
        IncompressibleMultiphaseBoundaryCondMap m_bcMap;
        double m_Reynolds; // Weissenberg number
        double[] Penalty1;
        double Penalty2;

        /// <summary>
        /// Initialize Convection
        /// </summary>
        public StressDivergenceInBulk(int _Component, IncompressibleMultiphaseBoundaryCondMap _BcMap, double _Reynolds, double[] _Penalty1, double _Penalty2, string spcName, SpeciesId spcId) : base(_Component, _BcMap, _Reynolds, _Penalty1, _Penalty2) {
            this.Component = _Component;
            this.m_spcId = spcId;
            this.m_bcMap = _BcMap;
            this.m_Reynolds = _Reynolds;
            this.Penalty1 = _Penalty1;
            this.Penalty2 = _Penalty2;
        }

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }
    }
}

/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer StRoeSTmungsdynamik (chair of fluid dynamics)

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
using XESTSF.Fluxes;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace XESTSF.Fluxes {

    public class RoeSTEnergyFlux_XDG : RoeSTEnergyFlux,ISpeciesFilter {

        private SpeciesId speciesId;

        public RoeSTEnergyFlux_XDG(IBoundaryConditionMap boundaryMap, Material material, string spcName, SpeciesId speciesId, double s_alpha, DGField[] _previous_u, int m_D = 2) : base(boundaryMap, material, s_alpha, _previous_u, m_D) {
            this.speciesId = speciesId;
            this.speciesName = spcName;
            this.ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }
}

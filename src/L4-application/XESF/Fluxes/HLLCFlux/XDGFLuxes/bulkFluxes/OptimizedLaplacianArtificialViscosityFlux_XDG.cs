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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace XDGShock.Fluxes {

    public class OptimizedLaplacianArtificialViscosityFlux_XDG : OptimizedLaplacianArtificialViscosityFlux, ISpeciesFilter {

        private SpeciesId speciesId;

        public OptimizedLaplacianArtificialViscosityFlux_XDG(GridData gridData, string ArgumentVarName, double penaltySafetyFactor, double penaltyFactor, MultidimensionalArray cellLengthScales, string spcName, SpeciesId speciesId) : base(gridData, ArgumentVarName, penaltySafetyFactor, penaltyFactor, cellLengthScales) {
            this.speciesId = speciesId;
            this.ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }
}

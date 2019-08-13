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
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Solution.XNSECommon.Operator.Pressure {

    /// <summary>
    /// pressure gradient, bulk phase, by central differences
    /// </summary>
    public class PressureInSpeciesBulk : PressureGradientLin_d, ISpeciesFilter {

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

        //string m_spcName;
        SpeciesId m_spcId;

        public PressureInSpeciesBulk(int _d, IncompressibleMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId)
            : base(_d, bcMap) {
            base.pressureFunction = bcMap.bndFunction[VariableNames.Pressure + "#" + spcName];
            this.m_bcMap = bcMap;
            //this.m_spcName = spcName;
            this.m_spcId = spcId;
        }

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }

        protected override double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {
            return base.BorderEdgeFlux(ref inp, Uin);
        }

        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {
            base.Flux(ref inp, U, output);
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            return base.InnerEdgeFlux(ref inp, Uin, Uout);
        }
    }
}

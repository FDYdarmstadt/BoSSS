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
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.XNSECommon.Operator.Continuity {

    /// <summary>
    /// velocity jump penalty for the divergence operator, on edges
    /// </summary>
    public class DivergenceInSpeciesBulk_Edge : Divergence_DerivativeSource_Flux, ISpeciesFilter {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="_component">
        /// component of the divergence
        /// </param>
        /// <param name="_bcmap"></param>
        public DivergenceInSpeciesBulk_Edge(int _component, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId, 
            double _rho, double _vorZeichen, bool _RescaleConti)
            : base(_component, _bcmap) {

            rho = _rho;
            m_spcId = spcId;
            //vorZeichen = _vorZeichen;
            this.RescaleConti = _RescaleConti;
            scale = _vorZeichen / ((RescaleConti) ? rho : 1.0);

            int d = base.component;
            base.bndFunction = _bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName];

            this.m_bcmap = _bcmap;
        }

        IncompressibleMultiphaseBoundaryCondMap m_bcmap;

        double rho;

        bool RescaleConti;
        double scale;

        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }



        protected override void InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout, out double FluxInCell, out double FluxOutCell) {
            base.InnerEdgeFlux(ref inp, Uin, Uout, out FluxInCell, out FluxOutCell);
            FluxInCell *= scale;
            FluxOutCell *= scale;
        }

        protected override void BorderEdgeFlux_(ref BoSSS.Foundation.CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            Debug.Assert(Uin.Length == 1);
            Debug.Assert(base.ArgumentOrdering.Count == 1);

            base.BorderEdgeFlux_(ref inp, Uin, out FluxInCell);

            FluxInCell *= scale;
        }



    }

    /// <summary>
    /// volume term for the Divergence / Continuity equation
    /// </summary>
    public class DivergenceInSpeciesBulk_Volume : Divergence_DerivativeSource, ISpeciesFilter {

        public DivergenceInSpeciesBulk_Volume(int _component, int _D, SpeciesId spcId, double _rho, double _vorZeichen, bool _RescaleConti)
            : base(_component, _D) {

            m_spcId = spcId;

            scale = _vorZeichen / ((_RescaleConti) ? _rho : 1.0);
        }

        double scale;

        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        public override double _DerivativeSource(double[] x, double[] Parameters, double[,] GradientU) {
            return base._DerivativeSource(x, Parameters, GradientU) * scale;
        }


    }


}

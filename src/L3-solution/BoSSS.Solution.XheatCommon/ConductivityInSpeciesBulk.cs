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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;


namespace BoSSS.Solution.XheatCommon {

    public class ConductivityInSpeciesBulk : swipConductivity, ISpeciesFilter {


        public ConductivityInSpeciesBulk(double penalty, double sw, ThermalMultiphaseBoundaryCondMap bcMap, int D,
            string spcName, SpeciesId spcId, double _kA, double _kB)
            : base(penalty, D, bcMap) {

            base.m_alpha = sw;
            this.m_bcMap = bcMap;

            this.m_spcId = spcId;

            switch (spcName) {
                case "A": currentk = _kA; complementk = _kB; break;
                case "B": currentk = _kB; complementk = _kA; break;
                default: throw new ArgumentException("Unknown species.");
            }

            double muFactor = Math.Max(currentk, complementk) / currentk;
            base.m_penalty_base = penalty * muFactor;

            base.tempFunction = this.m_bcMap.bndFunction[VariableNames.Temperature + "#" + spcName];
            base.fluxFunction = this.m_bcMap.bndFunction["HeatFlux#" + spcName];

        }


        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        double currentk = double.NaN;
        double complementk = double.NaN;


        ThermalMultiphaseBoundaryCondMap m_bcMap;

        /// <summary>
        /// multiplier for the penalty computation
        /// </summary>
        double m_penalty;


        protected override double Conductivity(double[] Parameters) {
            return currentk;
        }

    }



}

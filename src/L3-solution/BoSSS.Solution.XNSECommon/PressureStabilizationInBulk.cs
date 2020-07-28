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
using BoSSS.Solution.Utils;
using ilPSP;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Pressure stabilization for incompressible flows using an equal-order formulation.
    /// D. A. D. Pietro and A. Ern, Mathematical Aspects of Discontinuous Galerkin Methods. Springer Berlin Heidelberg, 2012.
    /// (Chapter 6.2.4.2)
    /// </summary>
    public class PressureStabilizationInBulk : PressureStabilization, ISpeciesFilter {

        double PressureStabilizationFactor;
        string m_spcName;
        SpeciesId m_spcId;
        //MultidimensionalArray h_max_Edge;

        /// <summary>
        /// Ctor.
        /// </summary>
        public PressureStabilizationInBulk(double PressureStabilizationFactor, double ReynoldsA, double ReynoldsB, string spcName, SpeciesId spcId): base(PressureStabilizationFactor, 0.0) {
            this.PressureStabilizationFactor = PressureStabilizationFactor;
            this.m_spcName = spcName;
            this.m_spcId = spcId;

            switch (spcName) {
                case "A": base.Reynolds = ReynoldsA; break;
                case "B": base.Reynolds = ReynoldsB; break;
                default: throw new ArgumentException("Unknown species.");
            }
        }

        public string validSpeciesId {
            get { return m_spcName; }
        }
    }
}

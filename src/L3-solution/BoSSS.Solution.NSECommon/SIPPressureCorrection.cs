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

using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// SIP discretization of Laplace operator for pressure correction,
    /// with penalty for pressure correction.
    /// </summary>
    public class IP1_Flux_PressureCorrection : SIPLaplace {

        IncompressibleBoundaryCondMap m_BcMap;

        /// <summary>
        /// Ctor.
        /// </summary>
        public IP1_Flux_PressureCorrection(double penalty_base, MultidimensionalArray PenaltyLengthScales, IncompressibleBoundaryCondMap _BcMap)
            : base(penalty_base, PenaltyLengthScales, VariableNames.PresCor) {
            m_BcMap = _BcMap;
        }

        protected override bool IsDirichlet(ref BoSSS.Foundation.CommonParamsBnd inp) {
            IncompressibleBcType edgType = m_BcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Outflow:
                    return true;
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                    return false;
                default:
                    throw new NotImplementedException("unsupported/unknown b.c. - missing implementation;");
            }
        }
    }
}

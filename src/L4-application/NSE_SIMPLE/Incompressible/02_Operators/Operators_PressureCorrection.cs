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
using BoSSS.Solution.NSECommon;

namespace NSE_SIMPLE {

    /// <summary>
    /// SIP discretization for pressure correction with penalty for [[p']].
    /// </summary>
    public class IP1_PressureCorrectionOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>        
        /// <param name="PressureMapping"></param>
        /// <param name="SolverConf"></param>
        public IP1_PressureCorrectionOperator(UnsetteledCoordinateMapping PressureMapping, SolverConfiguration SolverConf)
            : base(PressureMapping, PressureMapping, null, SolverConf, true) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            return (new IP1_Flux_PressureCorrection(SolverConf.PenaltyPressureCorrection, base.GridData.Cells.cj, SolverConf.BcMap)).Operator();
        }
    }
}

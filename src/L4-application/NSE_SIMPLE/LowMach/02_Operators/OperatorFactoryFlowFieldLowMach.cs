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

namespace NSE_SIMPLE {

    /// <summary>
    /// Operators flow field for Low-Mach number flows.
    /// </summary>
    public class OperatorFactoryFlowFieldLowMach : OperatorFactoryFlowFieldVariableDensity {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="VelocityMapping"></param>
        /// <param name="VelocityVectorMapping"></param>
        /// <param name="PressureMapping"></param>
        /// <param name="SolverConf"></param>
        /// <param name="Velocity0"></param>
        /// <param name="Velocity0Mean"></param>
        /// <param name="Temperature0"></param>
        /// <param name="Temperature0Mean"></param>
        /// <param name="eta"></param>
        public OperatorFactoryFlowFieldLowMach(UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping,
            UnsetteledCoordinateMapping PressureMapping,
            SolverConfiguration SolverConf,
            VectorField<SinglePhaseField> Velocity0, VectorField<SinglePhaseField> Velocity0Mean, SinglePhaseField Temperature0, SinglePhaseField Temperature0Mean,
            SinglePhaseField eta)
            : base(VelocityMapping, VelocityVectorMapping, PressureMapping, SolverConf, Velocity0, Velocity0Mean, Temperature0, Temperature0Mean, eta) { }

        protected override SIMPLEOperator GetDivergenceContiOperator(UnsetteledCoordinateMapping PressureMapping, UnsetteledCoordinateMapping VelocityMapping, SinglePhaseField Phi0, SolverConfiguration SolverConf, int SpatComp, SinglePhaseField[] MassFractionsCurrent = null) {
            return (new DivergenceLowMachOperator(PressureMapping, VelocityMapping, Phi0, SolverConf, SpatComp));
        }
    }
}

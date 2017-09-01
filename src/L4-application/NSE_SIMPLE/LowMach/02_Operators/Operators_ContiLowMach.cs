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
using BoSSS.Solution.NSECommon;
using NSE_SIMPLE.LowMach;

namespace NSE_SIMPLE {

    /// <summary>
    /// Divergence operator for Low-Mach continuity equation.
    /// </summary>
    public class DivergenceLowMachOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="PressureMapping"></param>
        /// <param name="VelocityMapping"></param>
        /// <param name="Temperature0"></param>
        /// <param name="SolverConf"></param>
        /// <param name="SpatialComponent"></param>
        public DivergenceLowMachOperator(UnsetteledCoordinateMapping PressureMapping, UnsetteledCoordinateMapping VelocityMapping,
            SinglePhaseField Temperature0, SolverConfiguration SolverConf, int SpatialComponent)
            : base(PressureMapping, VelocityMapping, new SinglePhaseField[] { Temperature0 }, SolverConf, false, SpatialComponent, SpatialComponent,
            MaxUsePerIterMatrix: 3, MaxUsePerIterAffine: 2) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            LowMachSIMPLEControl lowMachControl = SolverConf.Control as LowMachSIMPLEControl;
            return (new Divergence_CentralDifference(SpatialComponent, SolverConf.BcMap, lowMachControl.EoS)).Operator();
        }
    }    
}

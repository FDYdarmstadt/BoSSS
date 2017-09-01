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
using ilPSP.Utils;
using NSE_SIMPLE.LowMach;

namespace NSE_SIMPLE {

    /// <summary>
    /// Convection operator for temperature equation.
    /// </summary>
    public class TemperatureConvectionOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="TemperatureMapping"></param>
        /// <param name="Velocity0"></param>
        /// <param name="Velocity0Mean"></param>
        /// <param name="Temperature0"></param>
        /// <param name="Temperature0Mean"></param>
        /// <param name="SolverConf"></param>
        public TemperatureConvectionOperator(UnsetteledCoordinateMapping TemperatureMapping,
            VectorField<SinglePhaseField> Velocity0, VectorField<SinglePhaseField> Velocity0Mean, SinglePhaseField Temperature0, SinglePhaseField Temperature0Mean,
            SolverConfiguration SolverConf)
            : base(TemperatureMapping, TemperatureMapping,
            ArrayTools.Cat<SinglePhaseField>(Velocity0.ToArray(), Velocity0Mean.ToArray(), Temperature0, Temperature0Mean),
            SolverConf, false) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            LowMachSIMPLEControl lowMachControl = SolverConf.Control as LowMachSIMPLEControl;
            return (new LinearizedScalarConvection(SolverConf.SpatialDimension, SolverConf.BcMap, lowMachControl.EoS)).Operator();
        }
    }

    /// <summary>
    /// Heat conduction for temperature equation.
    /// </summary>
    public class HeatConductionOperator : SIMPLEOperator {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="TemperatureMapping"></param>        
        /// <param name="Temperature0"></param>
        /// <param name="SolverConf"></param>
        public HeatConductionOperator(UnsetteledCoordinateMapping TemperatureMapping, SinglePhaseField Temperature0, SolverConfiguration SolverConf)
            : base(TemperatureMapping, TemperatureMapping, new SinglePhaseField[] { Temperature0 }, SolverConf, false,
            MaxUsePerIterMatrix: 1, MaxUsePerIterAffine: 1) { }

        protected override SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int SpatialComponent, int SpatialDirection) {
            LowMachSIMPLEControl lowMachControl = SolverConf.Control as LowMachSIMPLEControl;
            return (new swipHeatConduction(lowMachControl.Reynolds, lowMachControl.Prandtl, lowMachControl.EoS, SolverConf.PenaltyHeatConduction, 
                base.GridData.Cells.cj,
                SolverConf.BcMap)).Operator();
        }
    }   
}

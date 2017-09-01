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
    /// Collection of operators for temperature equation.
    /// </summary>
    public class OperatorFactoryTemperature {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="TemperatureMapping"></param>
        /// <param name="Velocity0"></param>
        /// <param name="Velocity0Mean"></param>
        /// <param name="Temperature0"></param>
        /// <param name="Temperature0Mean"></param>
        /// <param name="SolverConf"></param>
        public OperatorFactoryTemperature(UnsetteledCoordinateMapping TemperatureMapping, UnsetteledCoordinateMapping PressureMapping,
            VectorField<SinglePhaseField> Velocity0, VectorField<SinglePhaseField> Velocity0Mean, SinglePhaseField Temperature0, SinglePhaseField Temperature0Mean,
            SolverConfiguration SolverConf) {

            this.TemperatureConvection = new TemperatureConvectionOperator(TemperatureMapping, Velocity0, Velocity0Mean, Temperature0, Temperature0Mean, SolverConf);            
            this.HeatConduction = new HeatConductionOperator(TemperatureMapping, Temperature0, SolverConf);            
        }

        /// <summary>
        /// Convection operator for temperature equation.
        /// </summary>
        public SIMPLEOperator TemperatureConvection {
            get;
            private set;
        }

        /// <summary>
        /// Heat conduction for temperature equation.       
        /// </summary>
        public SIMPLEOperator HeatConduction {
            get;
            private set;
        }
    }
}

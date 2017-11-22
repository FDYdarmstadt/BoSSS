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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution;
using CNS.EquationSystem;
using System;

namespace CNS.Solution {

    /// <summary>
    /// Factory for time steppers depending on the configuration options
    /// <see cref="CNSControl.TimeSteppingScheme"/>,
    /// <see cref="CNSControl.ExplicitScheme"/> and
    /// <see cref="CNSControl.ImplicitScheme"/>
    /// </summary>
    public class TimeStepperFactory {

        /// <summary>
        /// Information about the grid.
        /// </summary>
        private GridData gridData;

        /// <summary>
        /// The equation system to be solved
        /// </summary>
        private OperatorFactory equationSystem;

        private CNSControl config;

        private ISpeciesMap speciesMap;

        /// <summary>
        /// Constructs a new factory
        /// </summary>
        /// <param name="config">
        /// Configuration options
        /// </param>
        /// <param name="gridData">
        /// Information about the grid
        /// </param>
        /// <param name="equationSystem">
        /// The equation system to be solved
        /// </param>
        /// <param name="speciesMap"></param>
        public TimeStepperFactory(CNSControl config, GridData gridData, OperatorFactory equationSystem, ISpeciesMap speciesMap) {
            this.config = config;
            this.equationSystem = equationSystem;
            this.gridData = gridData;
            this.speciesMap = speciesMap;
        }

        /// <summary>
        /// Constructs a new time-stepper for the variables represented by
        /// <paramref name="variableMap"/>.
        /// </summary>
        /// <param name="variableMap">
        /// The variables that will be affected by the constructed time stepper
        /// </param>
        /// <param name="parameterMap">Parameter variables</param>
        /// <param name="program"></param>
        /// <returns>
        /// A time stepper for the variables represented by
        /// <paramref name="variableMap"/>.
        /// </returns>
        public ITimeStepper GetTimeStepper<T>(CoordinateMapping variableMap, CoordinateMapping parameterMap, Program<T> program)
            where T : CNSControl, new() {

            return config.ExplicitScheme.Instantiate(
                config,
                equationSystem,
                variableMap,
                parameterMap,
                speciesMap,
                program);
        }
    }
}

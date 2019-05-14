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


namespace CNS.EquationSystem {

    /// <summary>
    /// Implements a flux builder that doesn't do anything
    /// </summary>
    public class NullFluxBuilder : FluxBuilder {

        /// <summary>
        /// The one and only instance of this class
        /// </summary>
        public static readonly NullFluxBuilder Instance = new NullFluxBuilder();

        /// <summary>
        /// Initializes a new instance of the <see cref="NullFluxBuilder"/>
        /// class.
        /// </summary>
        private NullFluxBuilder()
            : base(null, null, null) {
        }

        /// <summary>
        /// Does nothing
        /// </summary>
        /// <param name="op">Irrelevant</param>
        public override void BuildFluxes(Operator op) {
        }
    }
}

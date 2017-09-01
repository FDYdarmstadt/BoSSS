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
using BoSSS.Foundation.Grid;
using CNS.Boundary;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CNS.Diffusion {

    /// <summary>
    /// Dummy Flux for the SIPG Density flux, not needed because of \f$ F_v(U_1)=0\f$ 
    /// </summary>
    public class SIPGDensityFlux : SIPGFlux {

        /// <summary>
        /// Index of the density component, i.e =0
        /// </summary>
        int component;

        /// <summary>
        /// ctor <see cref="SIPGFlux"/>
        /// </summary>
        /// <param name="config"><see cref="SIPGFlux"/></param>
        /// <param name="boundaryMap"><see cref="SIPGFlux"/></param>
        /// <param name="speciesMap"><see cref="SIPGFlux"/></param>
        /// <param name="gridData"><see cref="SIPGFlux"/></param>
        /// <param name="component"><see cref="SIPGFlux"/></param>
        /// <param name="cellMetricFunc"><see cref="SIPGFlux"/></param>
        public SIPGDensityFlux(CNSControl config, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, GridData gridData, int component, Func<MultidimensionalArray> cellMetricFunc)
            : base(config, boundaryMap, speciesMap, gridData, cellMetricFunc) {
            this.component = component;
        }

        /// <summary>
        /// Dummy method, no density update exists, i.e no update of the diffusivity tensor is needed.
        /// See lecture notes Hartmann2008.
        /// </summary>
        /// <param name="state"></param>
        /// <param name="adiabaticWall"></param>
        /// <param name="Tensor"></param>
        /// <param name="cellIndex"></param>
        protected override void UpdateTensorComponent(StateVector state, bool adiabaticWall, double[, ,] Tensor, int cellIndex) {
            // Do nothing
        }

        /// <summary>
        /// Returns the index for the density component, i.e. =0
        /// </summary>
        protected override int Component {
            get {
                return component;
            }
        }
    }
}

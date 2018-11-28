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
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.CompressibleFlowCommon;

namespace CNS.Diffusion {

    /// <summary>
    /// Implements the part of the SIPG Flux specific to the momentum equations
    /// </summary>
    public class SIPGMomentumFlux : SIPGFlux {

        /// <summary>
        /// The index of the momentum component represented by a specific
        /// instance of this class
        /// </summary>
        int component;

        /// <summary>
        /// ctor for the implementation of the SIPG momentum fluxes
        /// </summary>
        /// <param name="config"><see cref="SIPGFlux.config"/></param>
        /// <param name="boundaryMap"><see cref="SIPGFlux.boundaryMap"/></param>
        /// <param name="speciesMap"><see cref="SIPGFlux.speciesMap"/></param>
        /// <param name="gridData"><see cref="SIPGFlux.gridData"/></param>
        /// <param name="component"><see cref="SIPGFlux.Component"/></param>
        /// <param name="cellMetricFunc"><see cref="SIPGFlux"/></param>
        public SIPGMomentumFlux(CNSControl config, BoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, IGridData gridData, int component, Func<MultidimensionalArray> cellMetricFunc)
            : base(config, boundaryMap, speciesMap, gridData, cellMetricFunc) {
            if (component < 1 || component > CNSEnvironment.NumberOfDimensions) {
                throw new ArgumentOutOfRangeException("component");
            }
            this.component = component;
        }

        /// <summary>
        /// <see cref="SIPGFlux.UpdateTensorComponent(StateVector, bool, double[,,], int)"/>
        /// </summary>
        /// <param name="state"></param>
        /// <param name="adiabaticWall"></param>
        /// <param name="Tensor"></param>
        /// <param name="cellIndex"></param>
        protected override void UpdateTensorComponent(StateVector state, bool adiabaticWall, double[, ,] Tensor, int cellIndex) {
            //TO-DO: temperature dependent thermal conductivity/viscosity
            // double ThermalConductivity = state.ThermalConductivity

            double Viscosity = state.GetViscosity(cellIndex);
            double ThermalConductivity = Viscosity;
            double mu_rhoRe = Viscosity / config.ReynoldsNumber / state.Density;

            double v1 = state.Velocity[0];
            double v2 = state.Velocity[1];
            double v3 = state.Velocity[2];

            switch (dimension) {
                case 1:
                    Tensor[0, 0, 0] = -(alphaPlus43) * v1 *mu_rhoRe;
                    Tensor[0, 0, 1] = alphaPlus43 * mu_rhoRe;
                    //Tensor[0, 0, 2] = 0 * eta_rhoRe; 
                    break;
                case 2:
                    switch (component){
                        case 1: // x-momentum
                            //G_xx
                            Tensor[0, 0, 0] = -(alphaPlus43) * v1 * mu_rhoRe;
                            Tensor[0, 0, 1] = alphaPlus43 * mu_rhoRe;
                            //Tensor[0, 0, 2] = 0 * eta_rhoRe;
                            //Tensor[0, 0, 3] = 0 * eta_rhoRe;

                            //G_xy
                            Tensor[0, 1, 0] = -alphaMinus23 * v2 * mu_rhoRe;
                            //Tensor[0, 1, 1] = 0 * eta_rhoRe;
                            Tensor[0, 1, 2] = alphaMinus23 * mu_rhoRe;
                            //Tensor[0, 1, 3] = 0 * eta_rhoRe;

                            // G_yx
                            Tensor[1, 0, 0] = -v2 * mu_rhoRe;
                            //Tensor[1, 0, 1] = 0 * eta_rhoRe;
                            Tensor[1, 0, 2] = 1 * mu_rhoRe;
                            //Tensor[1, 0, 3] = 0 * eta_rhoRe;

                            //G_yy
                            Tensor[1, 1, 0] = -v1 * mu_rhoRe;
                            Tensor[1, 1, 1] = 1 * mu_rhoRe;
                            //Tensor[1, 1, 2] = 0 * eta_rhoRe;
                            //Tensor[1, 1, 3] = 0 * eta_rhoRe;

                            break;
                        case 2: // y-momentum
                            //G_xx
                            Tensor[0, 0, 0] = -v2 * mu_rhoRe;
                            //Tensor[0, 0, 1] = 0 * eta_rhoRe;
                            Tensor[0, 0, 2] = 1 * mu_rhoRe;
                            //Tensor[0, 0, 3] = 0 * eta_rhoRe;

                            //G_xy
                            Tensor[0, 1, 0] = -v2 * mu_rhoRe;
                            Tensor[0, 1, 1] = 1 * mu_rhoRe;
                            //Tensor[0, 1, 2] = 0 * eta_rhoRe;
                            //Tensor[0, 1, 3] = 0 * eta_rhoRe;

                            // G_yx
                            Tensor[1, 0, 0] = -alphaMinus23*v1 * mu_rhoRe;
                            Tensor[1, 0, 1] = alphaMinus23 * mu_rhoRe;
                            //Tensor[1, 0, 2] = 0 * eta_rhoRe;
                            //Tensor[1, 0, 3] = 0 * eta_rhoRe;

                            //G_yy
                            Tensor[1, 1, 0] = -(alphaPlus43) * v2 * mu_rhoRe;
                            //Tensor[1, 1, 1] = 0 * eta_rhoRe;
                            Tensor[1, 1, 2] = alphaPlus43 * mu_rhoRe;
                            //Tensor[1, 1, 3] = 0 * eta_rhoRe;
                            break;
                    }
                    break;
                case 3:
                    switch (component){
                        case 1: // x-momentum
                            //G_xx
                            Tensor[0, 0, 0] = -alphaPlus43 * v1 * mu_rhoRe;
                            Tensor[0, 0, 1] = alphaPlus43 * mu_rhoRe;
                            //Tensor[0, 0, 2] = 0 * eta_rhoRe;
                            //Tensor[0, 0, 3] = 0 * eta_rhoRe;
                            //Tensor[0, 0, 4] = 0 * eta_rhoRe;

                            //G_xy
                            Tensor[0, 1, 0] = -alphaMinus23 * v2 * mu_rhoRe;
                            //Tensor[0, 1, 1] = 0 * eta_rhoRe;
                            Tensor[0, 1, 2] = alphaMinus23 * mu_rhoRe;
                            //Tensor[0, 1, 3] = 0 * eta_rhoRe;
                            //Tensor[0, 1, 4] = 0 * eta_rhoRe;

                            //G_xz
                            Tensor[0, 2, 0] = -alphaMinus23 * v3 * mu_rhoRe;
                            //Tensor[0, 2, 1] = 0 * eta_rhoRe;
                            //Tensor[0, 2, 2] = 0 * eta_rhoRe;
                            Tensor[0, 2, 3] = alphaMinus23 * mu_rhoRe;
                            //Tensor[0, 2, 4] = 0 * eta_rhoRe;


                            // G_yx
                            Tensor[1, 0, 0] = -v2 * mu_rhoRe;
                            //Tensor[1, 0, 1] = 0 * eta_rhoRe;
                            Tensor[1, 0, 2] = 1 * mu_rhoRe;
                            //Tensor[1, 0, 3] = 0 * eta_rhoRe;
                            //Tensor[1, 0, 4] = 0 * eta_rhoRe;

                            //G_yy
                            Tensor[1, 1, 0] = -v1 * mu_rhoRe;
                            Tensor[1, 1, 1] = 1 * mu_rhoRe;
                            //Tensor[1, 1, 2] = 0 * eta_rhoRe;
                            //Tensor[1, 1, 3] = 0 * eta_rhoRe;
                            //Tensor[1, 1, 4] = 0 * eta_rhoRe;

                            //G_yz
                            //Tensor[1, 2, 0] = 0 * eta_rhoRe;
                            //Tensor[1, 2, 1] = 0 * eta_rhoRe;
                            //Tensor[1, 2, 2] = 0 * eta_rhoRe;
                            //Tensor[1, 2, 3] = 0 * eta_rhoRe;
                            //Tensor[1, 2, 4] = 0 * eta_rhoRe;

                            // G_zx
                            Tensor[2, 0, 0] = -v3 * mu_rhoRe;
                            //Tensor[2, 0, 1] = 0 * eta_rhoRe;
                            //Tensor[2, 0, 2] = 0 * eta_rhoRe;
                            Tensor[2, 0, 3] = 1 * mu_rhoRe;
                            //Tensor[2, 0, 4] = 0 * eta_rhoRe;

                            //G_zy
                            //Tensor[2, 1, 0] = 0 * eta_rhoRe;
                            //Tensor[2, 1, 1] = 0 * eta_rhoRe;
                            //Tensor[2, 1, 2] = 0 * eta_rhoRe;
                            //Tensor[2, 1, 3] = 0 * eta_rhoRe;
                            //Tensor[2, 1, 4] = 0 * eta_rhoRe;

                            //G_zz
                            Tensor[2, 2, 0] = -v1 * mu_rhoRe;
                            Tensor[2, 2, 1] = 1 * mu_rhoRe;
                            //Tensor[2, 2, 2] = 0 * eta_rhoRe;
                            //Tensor[2, 2, 3] = 0 * eta_rhoRe;
                            //Tensor[2, 2, 4] = 0 * eta_rhoRe;

                            break;
                        case 2: // y-momentum
                            //G_xx
                            Tensor[0, 0, 0] = -v2 * mu_rhoRe;
                            //Tensor[0, 0, 1] = 0 * eta_rhoRe;
                            Tensor[0, 0, 2] = 1 * mu_rhoRe;
                            //Tensor[0, 0, 3] = 0 * eta_rhoRe;
                            //Tensor[0, 0, 4] = 0 * eta_rhoRe;

                            //G_xy
                            Tensor[0, 1, 0] = -v1 * mu_rhoRe;
                            Tensor[0, 1, 1] = 1 * mu_rhoRe;
                            //Tensor[0, 1, 2] = 0 * eta_rhoRe;
                            //Tensor[0, 1, 3] = 0 * eta_rhoRe;
                            //Tensor[0, 1, 4] = 0 * eta_rhoRe;

                            //G_xz
                            //Tensor[0, 2, 0] = 0 * eta_rhoRe;
                            //Tensor[0, 2, 1] = 0 * eta_rhoRe;
                            //Tensor[0, 2, 2] = 0 * eta_rhoRe;
                            //Tensor[0, 2, 3] = 0 * eta_rhoRe;
                            //Tensor[0, 2, 4] = 0 * eta_rhoRe;

                            // G_yx
                            Tensor[1, 0, 0] = -alphaMinus23*v1 * mu_rhoRe;
                            Tensor[1, 0, 1] = alphaMinus23 * mu_rhoRe;
                            //Tensor[1, 0, 2] = 0 * eta_rhoRe;
                            //Tensor[1, 0, 3] = 0 * eta_rhoRe;
                            //Tensor[1, 0, 4] = 0 * eta_rhoRe;

                            //G_yy
                            Tensor[1, 1, 0] = -alphaPlus43 * v2 * mu_rhoRe;
                            //Tensor[1, 1, 1] = 0 * eta_rhoRe;
                            Tensor[1, 1, 2] = alphaPlus43 * mu_rhoRe;
                            //Tensor[1, 1, 3] = 0 * eta_rhoRe;
                            //Tensor[1, 1, 4] = 0 * eta_rhoRe;

                            // G_yz
                            Tensor[1, 2, 0] = -alphaMinus23 * v3 * mu_rhoRe;
                            //Tensor[1, 2, 1] = 0 * eta_rhoRe;
                            //Tensor[1, 2, 2] = 0 * eta_rhoRe;
                            Tensor[1, 2, 3] = alphaMinus23 * mu_rhoRe;
                            //Tensor[1, 2, 4] = 0 * eta_rhoRe;

                            //G_zx
                            //Tensor[2, 0, 0] = 0 * eta_rhoRe;
                            //Tensor[2, 0, 1] = 0 * eta_rhoRe;
                            //Tensor[2, 0, 2] = 0 * eta_rhoRe;
                            //Tensor[2, 0, 3] = 0 * eta_rhoRe;
                            //Tensor[2, 0, 4] = 0 * eta_rhoRe;

                            //G_zy
                            Tensor[2, 1, 0] = -v3 * mu_rhoRe;
                            //Tensor[2, 1, 1] = 0 * eta_rhoRe;
                            //Tensor[2, 1, 2] = 0 * eta_rhoRe;
                            Tensor[2, 1, 3] = 1 * mu_rhoRe;
                            //Tensor[2, 1, 4] = 0 * eta_rhoRe;

                            //G_zz
                            Tensor[2, 2, 0] = -v2 * mu_rhoRe;
                            //Tensor[2, 2, 1] = 0 * eta_rhoRe;
                            Tensor[2, 2, 2] = 1 * mu_rhoRe;
                            //Tensor[2, 2, 3] = 0 * eta_rhoRe;
                            //Tensor[2, 2, 4] = 0 * eta_rhoRe;
                            break;
                        case 3: // z-momentum
                            //G_xx
                            Tensor[0, 0, 0] = -v3 * mu_rhoRe;
                            //Tensor[0, 0, 1] = 0 * eta_rhoRe;
                            //Tensor[0, 0, 2] = 0 * eta_rhoRe;
                            Tensor[0, 0, 3] = 1 * mu_rhoRe;
                            //Tensor[0, 0, 4] = 0 * eta_rhoRe;

                            //G_xy
                            //Tensor[0, 1, 0] = 0 * eta_rhoRe;
                            //Tensor[0, 1, 1] = 0 * eta_rhoRe;
                            //Tensor[0, 1, 2] = 0 * eta_rhoRe;
                            //Tensor[0, 1, 3] = 0 * eta_rhoRe;
                            //Tensor[0, 1, 4] = 0 * eta_rhoRe;

                            //G_xz
                            Tensor[0, 2, 0] = -v1 * mu_rhoRe;
                            Tensor[0, 2, 1] = 1 * mu_rhoRe;
                            //Tensor[0, 2, 2] = 0 * eta_rhoRe;
                            //Tensor[0, 2, 3] = 0 * eta_rhoRe;
                            //Tensor[0, 2, 4] = 0 * eta_rhoRe;

                            // G_yx
                            //Tensor[1, 0, 0] = 0 * eta_rhoRe;
                            //Tensor[1, 0, 1] = 0 * eta_rhoRe;
                            //Tensor[1, 0, 2] = 0 * eta_rhoRe;
                            //Tensor[1, 0, 3] = 0 * eta_rhoRe;
                            //Tensor[1, 0, 4] = 0 * eta_rhoRe;

                            //G_yy
                            Tensor[1, 1, 0] = -v3 * mu_rhoRe;
                            //Tensor[1, 1, 1] = 0 * eta_rhoRe;
                            //Tensor[1, 1, 2] = 0 * eta_rhoRe;
                            Tensor[1, 1, 3] = 1 * mu_rhoRe;
                            //Tensor[1, 1, 4] = 0 * eta_rhoRe;

                            // G_yz
                            Tensor[1, 2, 0] = -v2 * mu_rhoRe;
                            //Tensor[1, 2, 1] = 0 * eta_rhoRe;
                            Tensor[1, 2, 2] = 1 * mu_rhoRe;
                            //Tensor[1, 2, 3] = 0 * eta_rhoRe;
                            //Tensor[1, 2, 4] = 0 * eta_rhoRe;

                            //G_zx
                            Tensor[2, 0, 0] =  -alphaMinus23 * v1 * mu_rhoRe;
                            Tensor[2, 0, 1] = alphaMinus23 * mu_rhoRe;
                            //Tensor[2, 0, 2] = 0 * eta_rhoRe;
                            //Tensor[2, 0, 3] = 0 * eta_rhoRe;
                            //Tensor[2, 0, 4] = 0 * eta_rhoRe;

                            //G_zy
                            Tensor[2, 1, 0] = -alphaMinus23* v2 * mu_rhoRe;
                            //Tensor[2, 1, 1] = 0 * eta_rhoRe;
                            Tensor[2, 1, 2] = alphaMinus23 * mu_rhoRe;
                            //Tensor[2, 1, 3] = 0 * eta_rhoRe;
                            //Tensor[2, 1, 4] = 0 * eta_rhoRe;

                            //G_zz
                            Tensor[2, 2, 0] = -alphaPlus43* v3 * mu_rhoRe;
                            //Tensor[2, 2, 1] = 0 * eta_rhoRe;
                            //Tensor[2, 2, 2] = 0 * eta_rhoRe;
                            Tensor[2, 2, 3] = alphaPlus43 * mu_rhoRe;
                            //Tensor[2, 2, 4] = 0 * eta_rhoRe;

                            break;
                    }
                    break;
                default:
                    throw new ArgumentOutOfRangeException("dimension");
            }
        }

        /// <summary>
        /// Gives the index of the momentum component
        /// x=1, y=2, z=3
        /// </summary>
        protected override int Component {
            get { return component; }
        }
    }
}

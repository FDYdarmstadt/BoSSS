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

using BoSSS.Foundation.Grid;
using System;
using CNS.Boundary;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace CNS.Diffusion {

    /// <summary>
    /// Implements the part of the SIPG Flux specific to the energy equation
    /// </summary>
    class SIPGEnergyFlux : SIPGFlux {

        /// <summary>
        /// The index of the energy component, i.e. =dimension+1
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
        public SIPGEnergyFlux(CNSControl config, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, GridData gridData, int component, Func<MultidimensionalArray> cellMetricFunc)
            : base(config, boundaryMap, speciesMap, gridData, cellMetricFunc) {
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

            double Viscosity = state.GetViscosity(cellIndex);
            double gamma = config.EquationOfState.HeatCapacityRatio;
            double gamma_Pr = adiabaticWall ? 0 :
                gamma / config.PrandtlNumber;
            double mu_rhoRe = Viscosity / config.ReynoldsNumber / state.Density;
            double MachsScaling = gamma * config.MachNumber * config.MachNumber;

            double v1 = state.Momentum[0] / state.Density;
            double v2 = state.Momentum[1] / state.Density;
            double v3 = state.Momentum[2] / state.Density;
            double E = state.Energy / state.Density;
            double VelocitySquared = v1 * v1 + v2*v2 + v3*v3;

            switch (dimension) {
                case 1:
                    Tensor[0, 0, 0] = -(MachsScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachsScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0, 0, 1] = MachsScaling * (alphaPlus43 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[0, 0, 2] = gamma_Pr * mu_rhoRe;
                    break;
                case 2:
                    // G_xx
                    Tensor[0, 0, 0] = -(MachsScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachsScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0, 0, 1] = MachsScaling * (alphaPlus43 - gamma_Pr)*v1*mu_rhoRe;
                    Tensor[0, 0, 2] = MachsScaling * (1-gamma_Pr)*v2*mu_rhoRe;
                    Tensor[0, 0, 3] = gamma_Pr*mu_rhoRe;

                    // G_xy
                    Tensor[0, 1, 0] = -MachsScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[0, 1, 1] = MachsScaling * v2 * mu_rhoRe;
                    Tensor[0, 1, 2] = MachsScaling * alphaMinus23 * v1 * mu_rhoRe;
                    //Tensor[0, 1, 3] = 0*eta_rhoRe;

                    // G_yx
                    Tensor[1, 0, 0] = -MachsScaling * (alphaPlus43-1.0)*v1*v2*mu_rhoRe;
                    Tensor[1, 0, 1] = MachsScaling * alphaMinus23 * v2*mu_rhoRe;
                    Tensor[1, 0, 2] = MachsScaling * v1 * mu_rhoRe;
                    //Tensor[1, 0, 3] = 0*eta_rhoRe;

                    // G_yy
                    Tensor[1, 1, 0] = -(MachsScaling * (alphaPlus13 * v2 * v2 + VelocitySquared) + gamma_Pr * (E - MachsScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[1, 1, 1] = MachsScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[1, 1, 2] = MachsScaling * (alphaPlus43 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[1, 1, 3] = gamma_Pr * mu_rhoRe;
                    break;

                case 3:
                    // G_xx
                    Tensor[0, 0, 0] = -(MachsScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachsScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0, 0, 1] = MachsScaling * (alphaPlus43 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[0, 0, 2] = MachsScaling * (1 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[0, 0, 3] = MachsScaling * (1 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[0, 0, 4] = gamma_Pr * mu_rhoRe;

                    // G_xy
                    Tensor[0, 1, 0] = -MachsScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[0, 1, 1] = MachsScaling * v2 * mu_rhoRe;
                    Tensor[0, 1, 2] = MachsScaling * alphaMinus23 * v1 * mu_rhoRe;
                    //Tensor[0, 1, 3] = 0 * eta_rhoRe;
                    //Tensor[0, 1, 4] = 0 * eta_rhoRe;

                    // G_xz
                    Tensor[0, 2, 0] = -MachsScaling * alphaPlus13 * v1 * v3 * mu_rhoRe;
                    Tensor[0, 2, 1] = MachsScaling * v3 * mu_rhoRe;
                    //Tensor[0, 2, 2] = 0 * eta_rhoRe;
                    Tensor[0, 2, 3] = MachsScaling * alphaMinus23 * v1 * mu_rhoRe;
                    //Tensor[0, 2, 4] = 0 * eta_rhoRe;


                    // G_yx
                    Tensor[1, 0, 0] = -MachsScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[1, 0, 1] = MachsScaling * alphaMinus23 * v2 * mu_rhoRe;
                    Tensor[1, 0, 2] = MachsScaling * v1 * mu_rhoRe;
                    //Tensor[1, 0, 3] = 0 * eta_rhoRe;
                    //Tensor[1, 0, 4] = 0 * eta_rhoRe;


                    // G_yy
                    Tensor[1, 1, 0] = -(MachsScaling * (alphaPlus13 * v2 * v2 + VelocitySquared) + gamma_Pr * (E - MachsScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[1, 1, 1] = MachsScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[1, 1, 2] = MachsScaling * (alphaPlus43 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[1, 1, 3] = MachsScaling * (1 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[1, 1, 4] = gamma_Pr * mu_rhoRe;

                    // G_yz
                    Tensor[1, 2, 0] = -MachsScaling * alphaPlus13 * v2 * v3 * mu_rhoRe;
                    //Tensor[1, 2, 1] = 0 * eta_rhoRe;
                    Tensor[1, 2, 2] = MachsScaling * v3 * mu_rhoRe;
                    Tensor[1, 2, 3] = MachsScaling * alphaMinus23 * v2 * mu_rhoRe;
                    //Tensor[1, 2, 4] = 0 * eta_rhoRe;

                    // G_zx
                    Tensor[2, 0, 0] = -MachsScaling * alphaPlus13 * v1 * v3 * mu_rhoRe;
                    Tensor[2, 0, 1] = MachsScaling * alphaMinus23 * v3 * mu_rhoRe;
                    //Tensor[2, 0, 2] = 0 * eta_rhoRe;
                    Tensor[2, 0, 3] = MachsScaling * v1 * mu_rhoRe;
                    //Tensor[2, 0, 4] = 0 * eta_rhoRe;

                    // G_zy
                    Tensor[2, 1, 0] = -MachsScaling * alphaPlus13 * v2 * v3 * mu_rhoRe;
                    //Tensor[2, 1, 1] = 0 * eta_rhoRe;
                    Tensor[2, 1, 2] = MachsScaling * alphaMinus23 * v3 * mu_rhoRe;
                    Tensor[2, 1, 3] = MachsScaling * v2 * mu_rhoRe;
                    //Tensor[2, 1, 4] = 0 * eta_rhoRe;

                    // G_zz
                    Tensor[2, 2, 0] = -(MachsScaling * (alphaPlus13 * v3 * v3 + VelocitySquared) + gamma_Pr * (E - MachsScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[2, 2, 1] = MachsScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[2, 2, 2] = MachsScaling * (1 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[2, 2, 3] = MachsScaling * (alphaPlus43 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[2, 2, 4] = gamma_Pr * mu_rhoRe;

                    break;
            }
        }

        /// <summary>
        /// Gives the index of the energy component:
        ///  = dimension+1
        /// </summary>
        protected override int Component {
            get { return component; }
        }
    }
}

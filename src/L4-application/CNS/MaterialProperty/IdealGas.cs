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

using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using System;

namespace CNS.MaterialProperty {

    /// <summary>
    /// Ideal gas law for gases like air at moderate conditions.
    /// </summary>
    [Serializable]
    public class IdealGas : IEquationOfState {

        /// <summary>
        /// Constructs a new ideal gas law
        /// </summary>
        /// <param name="heatCapacityRatio">
        /// The heat capacity ratio
        /// </param>
        public IdealGas(double heatCapacityRatio) {
            this.HeatCapacityRatio = heatCapacityRatio;
        }

        /// <summary>
        /// Constructs an ideal gas law for standard air
        /// </summary>
        public static IdealGas Air {
            get {
                return new IdealGas(1.4);
            }
        }

        /// <summary>
        /// Constructs an ideal gas law for Helium
        /// </summary>
        public static IdealGas Helium {
            get {
                return new IdealGas(1.66);
            }
        }

        #region IEquationOfState Members

        /// <summary>
        /// The heat capacity ratio.
        /// </summary>
        public double HeatCapacityRatio {
            get;
            private set;
        }

        /// <summary>
        /// Calculate the pressure $p$ via
        /// $p = (\kappa - 1) \rho e$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetPressure"/>
        /// </param>
        /// <returns>
        /// \f$ \rho e (\kappa - 1)\f$ 
        /// where
        /// \f$ \rho e\f$  = <see name="StateVector.InnerEnergy"/>
        /// and \f$ \kappa\f$  is the heat capacity ratio
        /// supplied to <see cref="IdealGas.IdealGas"/>.
        /// </returns>
        public double GetPressure(StateVector state) {
            return (HeatCapacityRatio - 1.0) * state.InnerEnergy;
        }

        /// <summary>
        /// Calculates the temperature \f$ T\f$  via
        /// \f$ T = (\kappa - 1) e\f$ 
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetTemperature"/>
        /// </param>
        /// <returns>
        /// \f$ (\kappa - 1) e\f$ 
        /// where
        /// \f$ e\f$  = <see name="StateVector.InnerEnergy"/>
        /// and \f$ \kappa\f$  is the heat capacity ratio.
        /// </returns>
        public double GetTemperature(StateVector state) {
            return (HeatCapacityRatio - 1.0) * state.SpecificInnerEnergy;
        }

        /// <summary>
        /// Calculates the inner energy for an ideal gas.
        /// </summary>
        /// <param name="density"></param>
        /// <param name="pressure"></param>
        /// <returns>
        /// \f$ \frac{p}{(\kappa - 1)}\f$  where
        /// \f$ \kappa\f$  is the heat capacity ratio.
        /// </returns>
        public double GetInnerEnergy(double density, double pressure) {
            return pressure / (HeatCapacityRatio - 1.0);
        }

        /// <summary>
        /// Calculates the local speed of sound \f$a\f$ via
        /// \f$a = \sqrt{\frac{p}{\rho}} / \text{Ma}\f$, where
        /// \f$\text{Ma}\f$ is the reference Mach number
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetSpeedOfSound"/>
        /// </param>
        /// <returns>
        /// \f$ \sqrt{\frac{p}{\rho}} / \text{Ma}\f$ 
        /// where
        /// \f$ \rho\f$  = <see name="StateVector.Density"/>
        /// \f$ p\f$  = <see name="StateVector.Pressure"/>.
        /// </returns>
        public double GetSpeedOfSound(StateVector state) {
            // Equals sqrt{\kappa \frac{p}{\rho}} for Ma = \kappa
            return Math.Sqrt(state.Pressure / state.Density) / state.Material.MachNumber;
        }

        /// <summary>
        /// Calculates the local entropy \f$ S\f$ via
        /// \f$ S = \frac{p}{\rho^\kappa}\f$ 
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetEntropy"/>
        /// </param>
        /// <returns>
        /// \f$ S = \frac{p}{\rho^\kappa}\f$ 
        /// where
        /// \f$ \rho\f$  = <see name="StateVector.Density"/>
        /// \f$ p\f$  = <see name="StateVector.Pressure"/>
        /// and \f$ \kappa\f$  is the heat capacity ratio
        /// supplied to <see cref="IdealGas.IdealGas"/>.
        /// </returns>
        public double GetEntropy(StateVector state) {
            return state.Pressure / Math.Pow(state.Density, HeatCapacityRatio);
        }

        #endregion
    }
}

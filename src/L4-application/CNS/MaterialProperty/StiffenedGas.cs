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
    /// Model equation for stiff fluids (e.g. water) where an ideal gas at a
    /// very pressure is assumed. For details, e.g. see ChangLiou2007. The
    /// generic form of the equation is
    /// \f$ 
    /// p = (\kappa - 1) \rho e - \kappa \pi
    /// \f$ 
    /// with determining constants \f$ \kappa\f$  and
    /// \f$ \pi\f$  (see
    /// <see cref="ReferencePressure"/>).
    /// </summary>
    public class StiffenedGas : IEquationOfState {

        /// <summary>
        /// <see cref="StiffenedGas.StiffenedGas"/>
        /// </summary>
        public readonly double ReferencePressure;

        /// <summary>
        /// Constructs a stiff fluid.
        /// </summary>
        /// <param name="heatCapacityRatio">
        /// The heat capacity ratio.
        /// </param>
        /// <param name="referencePressure">
        /// The reference pressure. A value of zero corresponds to an ideal
        /// gas.
        /// </param>
        public StiffenedGas(double heatCapacityRatio, double referencePressure) {
            this.HeatCapacityRatio = heatCapacityRatio;
            this.ReferencePressure = referencePressure;
        }

        /// <summary>
        /// Instantiates the specific equation of state water.
        /// </summary>
        public static StiffenedGas Water {
            get {
                return new StiffenedGas(7.0, 6.0e8);
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
        /// Calculates the pressure $p$ via
        /// $p = \rho e (\kappa - 1) - \kappa p_\infty$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetPressure"/>
        /// </param>
        /// <returns>
        /// \f$ \rho e (\kappa - 1) - \kappa p_\infty\f$ 
        /// where
        /// \f$ \rho\f$  = <see name="StateVector.Density"/>,
        /// \f$ e\f$  = <see name="StateVector.InnerEnergy"/>
        /// and \f$ \kappa\f$  and
        /// \f$ p_\infty\f$  are the heat capacity ratio
        /// and the reference pressure supplied to
        /// <see cref="StiffenedGas.StiffenedGas"/>, respectively.
        /// </returns>
        public double GetPressure(StateVector state) {
            return (HeatCapacityRatio - 1.0) * state.Density * state.SpecificInnerEnergy - HeatCapacityRatio * ReferencePressure;
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
        /// Calculates the inner energy for a stiffened gas.
        /// </summary>
        /// <param name="density"></param>
        /// <param name="pressure"></param>
        /// <returns>
        /// \f$ \frac{p + \kappa \pi}{(\kappa - 1)}\f$  where
        /// \f$ \kappa\f$  is the heat capacity ratio.
        /// </returns>
        public double GetInnerEnergy(double density, double pressure) {
            return (pressure + HeatCapacityRatio * ReferencePressure) / (HeatCapacityRatio - 1.0);
        }

        /// <summary>
        /// Calculates the speed of sound $a$ via
        /// $a = \sqrt{\kappa \frac{p + \pi}{\rho}}$.
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetSpeedOfSound"/>
        /// </param>
        /// <returns>
        /// \f$ \sqrt{\kappa \frac{p + \pi}{\rho}}\f$ 
        /// where
        /// \f$ \rho\f$  = <see name="StateVector.Density"/>,
        /// \f$ p\f$  = <see name="StateVector.Pressure"/>
        /// and \f$ \kappa\f$  and
        /// \f$ \pi\f$  are the heat capacity ratio
        /// and the reference pressure supplied to
        /// <see cref="StiffenedGas.StiffenedGas"/>, respectively.
        /// </returns>
        public double GetSpeedOfSound(StateVector state) {
            return Math.Sqrt(HeatCapacityRatio * (state.Pressure + ReferencePressure) / state.Density);
        }

        /// <summary>
        /// Calculates the entropy $S$ via
        /// $S = \frac{p + \pi}{\rho^\kappa}$
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetEntropy"/>
        /// </param>
        /// <returns>
        /// \f$ \frac{p + \pi}{\rho^\kappa}\f$ 
        /// where
        /// \f$ \rho\f$  = <see name="StateVector.Density"/>,
        /// \f$ p\f$  = <see name="StateVector.Pressure"/>
        /// and \f$ \kappa\f$  and
        /// \f$ \pi\f$  are the heat capacity ratio
        /// and the reference pressure supplied to
        /// <see cref="StiffenedGas.StiffenedGas"/>, respectively.
        /// </returns>
        /// <remarks>
        /// cf. FarhatEtAl2008
        /// </remarks>
        public double GetEntropy(StateVector state) {
            return (state.Pressure + ReferencePressure) / Math.Pow(state.Density, HeatCapacityRatio);
        }

        #endregion
    }
}

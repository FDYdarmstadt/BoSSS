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
using System.Diagnostics;

namespace CNS.MaterialProperty {

    /// <summary>
    /// An equation of state for modified ideal gases with a lower bound for
    /// the volume of a fluid particle (a.k.a. covolume). The respective
    /// equation of state is sometimes also called Nobel-Abel-EOS. For more
    /// details on the subject, e.g. see Toro2009 or Johnston2005.
    /// </summary>
    /// <remarks>
    /// The covolume (see <see cref="Covolume"/>) is generally denoted by
    /// \f$ b\f$  in formulas.
    /// </remarks>
    [Serializable]
    public class CovolumeGas : IEquationOfState {

        /// <summary>
        /// A lower bound for the volume of a fluid particle; essentially
        /// defines an upper bound for the density given by
        /// \f$ \rho_\text{max} = \frac{1}{b}\f$ .
        /// </summary>
        public readonly double Covolume;

        /// <summary>
        /// Constructs a new covolume gas.
        /// </summary>
        /// <param name="heatCapacityRatio">Heat capacity ratio.</param>
        /// <param name="covolume">
        /// A lower bound for the volume of a fluid particle; essentially
        /// defines an upper bound for the density given by
        /// \f$ \rho_\text{max} = \frac{1}{b}\f$ .
        /// </param>
        public CovolumeGas(double heatCapacityRatio, double covolume) {
            this.HeatCapacityRatio = heatCapacityRatio;
            this.Covolume = covolume;
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
        /// Calculates the pressure.
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetPressure"/>
        /// </param>
        /// <returns>
        /// \f$ 
        /// p = (\kappa - 1) \frac{\rho}{1 - \rho b} e
        /// \f$ 
        /// </returns>
        public double GetPressure(StateVector state) {
            double apparentDensity = state.Density / (1.0 - Covolume * state.Density);
            Debug.Assert(
                apparentDensity > 0.0,
                "Invalid density (might be bigger than 1 / covolume)");

            return (HeatCapacityRatio - 1.0) * apparentDensity * state.SpecificInnerEnergy;
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
        /// Calculates the inner energy for a covolume gas.
        /// </summary>
        /// <param name="density"></param>
        /// <param name="pressure"></param>
        /// <returns>
        /// \f$ \frac{p}{\kappa - 1} (1 - \rho b)\f$  where
        /// \f$ \kappa\f$  is the heat capacity ratio
        /// and \f$ b\f$  is the covolume.
        /// </returns>
        public double GetInnerEnergy(double density, double pressure) {
            return pressure / (HeatCapacityRatio - 1.0) * (1.0 - density * Covolume);
        }


        /// <summary>
        /// Calculates the speed of sound.
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetSpeedOfSound"/>
        /// </param>
        /// <returns>
        /// \f$ 
        /// a = \sqrt{\kappa \frac{p}{\rho (1 - \rho b)}}
        /// \f$ 
        /// </returns>
        public double GetSpeedOfSound(StateVector state) {
            return Math.Sqrt(HeatCapacityRatio * state.Pressure / state.Density / (1.0 - state.Density * Covolume));
        }

        /// <summary>
        /// Calculates the local entropy.
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetEntropy"/>
        /// </param>
        /// <returns>
        /// \f$ 
        /// S = p \left(\frac{1 - \rho b}{\rho}\right)^\gamma
        /// \f$ 
        /// </returns>
        public double GetEntropy(StateVector state) {
            return state.Pressure * Math.Pow(
                (1.0 - state.Density * Covolume) / state.Density,
                HeatCapacityRatio);
        }

        #endregion
    }
}

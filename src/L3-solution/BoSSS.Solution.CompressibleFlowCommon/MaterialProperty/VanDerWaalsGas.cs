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

namespace BoSSS.Solution.CompressibleFlowCommon.MaterialProperty {

    /// <summary>
    /// Equation of state for a <i>real</i> gas described by the van der Waals
    /// law. The equation of state models the influence of attractive forces of
    /// molecules at densities (see <see cref="attraction"/>) as well as a
    /// lower bound for the admissible volume of a fluid particle (a.k.a.
    /// covolume, see <see cref="covolume"/>):
    /// </summary>
    /// <remarks>
    /// The parameter defining the attraction of molecules (see
    /// <see cref="attraction"/>) is denoted by \f$ a\f$ 
    /// in the following formulas. The covolume (see <see cref="covolume"/>) is
    /// denoted by \f$ b\f$ .
    /// </remarks>
    public class VanDerWaalsGas : IEquationOfState {

        /// <summary>
        /// A lower bound for the volume of a fluid particle; essentially
        /// defines an upper bound for the density given by
        /// \f$ \rho_\text{max} = \frac{1}{b}\f$ .
        /// </summary>
        private double covolume;

        /// <summary>
        /// Parameter defining the influence of attractive forces between
        /// molecules at high densities.
        /// </summary>
        private double attraction;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="heatCapacityRatio">Heat capacity ratio</param>
        /// <param name="attraction">
        /// Parameter defining the influence of attractive forces between
        /// molecules at high densities.
        /// </param>
        /// <param name="covolume">
        /// A lower bound for the volume of a fluid particle; essentially
        /// defines an upper bound for the density given by
        /// \f$ \rho_\text{max} = \frac{1}{b}\f$ .
        /// </param>
        public VanDerWaalsGas(double heatCapacityRatio, double attraction, double covolume) {
            this.HeatCapacityRatio = heatCapacityRatio;
            this.covolume = covolume;
            this.attraction = attraction;
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
        /// p = (\kappa - 1) \frac{\rho}{1 - \rho b} e - a \rho^2
        /// \f$ 
        /// </returns>
        public double GetPressure(StateVector state) {
            double apparentDensity = state.Density / (1.0 - covolume * state.Density);
            Debug.Assert(
                apparentDensity > 0.0,
                "Invalid density (might be bigger than 1 / covolume)");

            return (HeatCapacityRatio - 1.0) * apparentDensity * state.SpecificInnerEnergy
                - attraction * state.Density * state.Density;
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <param name="state"></param>
        /// <returns></returns>
        /// <remarks>
        /// Beware: VDW gas implies non-constant heat capacities
        /// </remarks>
        public double GetTemperature(StateVector state) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates the inner energy for a VDW gas.
        /// </summary>
        /// <param name="density"></param>
        /// <param name="pressure"></param>
        /// <returns>
        /// \f$ \frac{p + a \rho^2}{\kappa - 1} (1 - \rho b)\f$ 
        /// where \f$ \kappa\f$  is the heat capacity
        /// ratio, \f$ a\f$  is the attraction coefficient
        /// and \f$ b\f$  is the covolume.
        /// </returns>
        public double GetInnerEnergy(double density, double pressure) {
            return (pressure + attraction * density * density) /
                (HeatCapacityRatio - 1.0) * (1.0 - density * covolume);
        }

        /// <summary>
        /// Calculates the speed of sound.
        /// </summary>
        /// <param name="state">
        /// <see cref="IEquationOfState.GetSpeedOfSound"/>
        /// </param>
        /// <returns>
        /// \f$ 
        /// a = \sqrt{\frac{\kappa p + a \rho^2 (2 b \rho - 1)}{\rho (1 - \rho b)}}
        /// \f$ 
        /// </returns>
        public double GetSpeedOfSound(StateVector state) {
            double p = state.Pressure;
            double rho = state.Density;
            return Math.Sqrt((HeatCapacityRatio * p + attraction * rho * rho * (2.0 * covolume * rho - 1.0)) /
                rho / (1.0 - covolume * rho));
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        /// <param name="state"></param>
        /// <returns></returns>
        public double GetEntropy(StateVector state) {
            throw new NotImplementedException();
        }

        #endregion
    }
}

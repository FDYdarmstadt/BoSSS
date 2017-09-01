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

namespace CNS.MaterialProperty {

    /// <summary>
    /// The dimensionless Sutherland law, see for example Iannelli2012.
    /// </summary>
    [System.Serializable]
    public class SutherlandLaw : IViscosityLaw {

        /// <summary>
        /// Ratio between the material parameter
        /// \f$ S\f$  (in Kelvin) and reference
        /// temperature \f$ T_0\f$  (in Kelvin), i.e.
        /// \f$ c = \frac{S}{T_0}\f$ 
        /// </summary>
        private double temperatureRatio;

        /// <summary>
        /// Constructs Sutherland's law with the given material parameter.
        /// </summary>
        /// <param name="temperatureRatio">
        /// Ratio between the material parameter
        /// \f$ S\f$  (in Kelvin) and reference
        /// temperature \f$ T_0\f$  (in Kelvin), i.e.
        /// \f$ c = \frac{S}{T_0}\f$ 
        /// </param>
        public SutherlandLaw(double temperatureRatio = 110.4 / 273.15) {
            this.temperatureRatio = temperatureRatio;
        }

        #region IViscosityLaw Members

        /// <summary>
        /// Determines the viscosity using Sutherland's law.
        /// </summary>
        /// <param name="temperature">
        /// The dimensionless fluid temperature.
        /// </param>
        /// <param name="cellIndex"></param>
        /// <returns>
        /// \f$ T^{3/2} \frac{1 + c}{T + c}\f$ ,
        /// where \f$ c\f$  is the temperatureRatio
        /// defined in <see cref="SutherlandLaw.SutherlandLaw"/>
        /// </returns>
        public double GetViscosity(double temperature, int cellIndex) {
            return Math.Pow(temperature, 1.5) * (1 + temperatureRatio) / (temperature + temperatureRatio);
        }

        #endregion
    }
}

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
    /// A dimensionless power law for the viscosity, e.g. see Iannelli2012.
    /// </summary>
    [System.Serializable]
    public class PowerLaw : IViscosityLaw {

        /// <summary>
        /// The power to be used in the power expression.
        /// </summary>
        private double power;

        /// <summary>
        /// Constructs a new power law
        /// </summary>
        /// <param name="power">
        /// The power to be used in the power expression
        /// </param>
        public PowerLaw(double power = 0.666) {
            this.power = power;
        }

        #region IViscosityLaw Members

        /// <summary>
        /// Computes the viscosity according to a simple power law
        /// </summary>
        /// <param name="temperature"></param>
        /// <param name="cellIndex"></param>
        /// <returns>
        /// \f$ T^\omega\f$ , where \f$ \omega\f$ 
        /// is an exponent determined from a fit to experiments
        /// </returns>
        public double GetViscosity(double temperature, int cellIndex) {
            return Math.Pow(temperature, power);
        }

        #endregion
    }
}

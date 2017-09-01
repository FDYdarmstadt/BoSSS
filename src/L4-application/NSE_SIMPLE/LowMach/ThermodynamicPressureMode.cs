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

namespace NSE_SIMPLE.LowMach {

    /// <summary>
    /// Option for calculating thermodynamic pressure in Low-Mach approximation, i.e. p0.
    /// </summary>
    public enum ThermodynamicPressureMode {

        /// <summary>
        /// Constant in time.
        /// This option is usually used for open systems,
        /// where p0 is set to the ambient pressure.
        /// </summary>
        Constant,

        /// <summary>
        /// This option can be used for closed systems,
        /// where p0 is determined to conserve mass, i.e.
        /// \f$ p_0 = \frac{m(t=0)}{\int 1/T dV} \f$ .
        /// </summary>
        MassDetermined
    }
}

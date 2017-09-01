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

using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using System;

namespace NSE_SIMPLE.BaseVariableDensity {

    /// <summary>
    /// Provides settings for the variable density SIMPLE algorithm,
    /// i.e. LowMach and multiphase (smooth interface) solvers.
    /// </summary>
    public class VariableDensitySIMPLEControl : SIMPLEControl {

        /// <summary>
        /// Material law for calculating density and viscosity.
        /// </summary>
        [NotNull]
        public MaterialLaw EoS;

        /// <summary>
        /// Froude number.
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public double? Froude;

        /// <summary>
        /// Unity spatial direction vector of gravity.
        /// </summary>
        public double[] GravityDirection;

        /// <summary>
        /// Analytic Solution density.
        /// </summary>
        public Func<double[], double> AnalyticDensity;
        
        /// <summary>
        /// 
        /// </summary>
        public override void Verify() {
            base.Verify();

            if (this.Froude.HasValue) {
                if (GravityDirection == null)
                    throw new ApplicationException("Missing required extended property 'GravityDirection'.");
            }
        }
    }
}

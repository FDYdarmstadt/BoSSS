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
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;


namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Reaction heat source in temperature equation.
    /// </summary>
    public class ReactionHeatSourceLinearized : BoSSS.Solution.Utils.LinearSource {
        string[] m_ArgumentOrdering;
        string[] m_ParameterOrdering;
        double ReactionRate;
        double HeatReleaseFactor;
        double[] ReactionRateConstants;
        double OneOverMolarMass0MolarMass1;
        bool m_speciesTransportOK;
        MaterialLaw EoS;
        double rho;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="HeatReleaseFactor">Heat release computed from the sum of the product of the stoichiometric coefficient, partial heat capacity and molar mass of species alpha for all species. I.e.: sum(alpha = 1.. ns)[v_\alpha cp_alpha M_alpha]. Must be computed locally for non-constant partial heat capacities in later iterations of the code.</param>   
        /// <param name="ReactionRateConstants">0. PreExpFactor/Damköhler number, 1. ActivationTemperature, 2. MassFraction0Exponent, 3. MassFraction1Exponent</param>  
        /// <param name="OneOverMolarMass0MolarMass1"> 1/(M_infty^(a + b -1) * MolarMassFuel^a * MolarMassOxidizer^b). M_infty is the reference for the molar mass steming from non-dimensionalisation of the governing equations.</param>  
        /// <param name="EoS">MaterialLawCombustion</param>  
        public ReactionHeatSourceLinearized(double HeatReleaseFactor, double[] ReactionRateConstants, double OneOverMolarMass0MolarMass1, MaterialLaw EoS, bool speciesTransportOK) {
            m_ArgumentOrdering = new string[] { VariableNames.Temperature };
            m_ParameterOrdering = new string[] { VariableNames.Temperature0, VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 };
            this.HeatReleaseFactor = HeatReleaseFactor;
            this.ReactionRateConstants = ReactionRateConstants;
            this.OneOverMolarMass0MolarMass1 = OneOverMolarMass0MolarMass1;
            this.EoS = EoS;
            this.m_speciesTransportOK = speciesTransportOK;
        }


        /// <summary>
        /// Temperature
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }

        /// <summary>
        /// Temperature, MassFraction0, MassFraction1, MassFraction2, ..., MassFraction_ns at linearization point.
        /// </summary>
        public override IList<string> ParameterOrdering {
            get { return m_ParameterOrdering; }
        }

        /// <summary>
        /// 
        /// </summary>
        protected override double Source(double[] x, double[] parameters, double[] U) {
            rho = EoS.GetDensity(parameters);
            Debug.Assert(!double.IsNaN(rho));
            Debug.Assert(!double.IsInfinity(rho));

            if (m_speciesTransportOK) {              
                ReactionRate = ReactionRateConstants[0] * Math.Exp(-ReactionRateConstants[1] / parameters[0]) * OneOverMolarMass0MolarMass1 * Math.Pow(rho * parameters[1], ReactionRateConstants[2]) * Math.Pow(rho * parameters[2], ReactionRateConstants[3]);

           
            }
            else {
                ReactionRate = 0.0;
            }

            return HeatReleaseFactor * U[0] * ReactionRate;
            
        }
    }
}

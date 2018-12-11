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
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using System.Diagnostics;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Reaction species source in mass transport equation
    /// </summary>
    public class ReactionSpeciesSourceLinearized : BoSSS.Solution.Utils.LinearSource {
        string[] m_ArgumentOrdering;

        double[] StoichiometricCoefficients;
        double[] ReactionRateConstants;
        int SpeciesIndex; //Species index, not to be confused with alpha = SpeciesIndex + 1
        int NumberOfReactants;
        string[] MassFractionNames;
        double OneOverMolarMass0MolarMass1;
        double[] MolarMasses;
        double rho;
        MaterialLaw EoS;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="ReactionRateConstants">constants[0]=PreExpFactor, constants[1]=ActivationTemperature, constants[2]=MassFraction0Exponent, constants[3]=MassFraction1Exponent</param>
        /// <param name="StoichiometricCoefficients"></param>        
        /// <param name="OneOverMolarMass0MolarMass1"> 1/(M_infty^(a + b -1) * MolarMassFuel^a * MolarMassOxidizer^b). M_infty is the reference for the molar mass steming from non-dimensionalisation of the governing equations.</param>  
        /// <param name="MolarMasses">Array of molar masses. 0 Fuel. 1 Oxidizer, 2 to ns products.</param>   
        /// <param name="EoS">MaterialLawCombustion</param>  
        /// <param name="NumerOfReactants">The number of reactants (i.e. ns)</param> 
        /// <param name="SpeciesIndex">Index of the species being balanced. (I.e. 0 for fuel, 1 for oxidizer, 2 for CO2, 3 for water)</param> 
        public ReactionSpeciesSourceLinearized(double[] ReactionRateConstants, double[] StoichiometricCoefficients, double OneOverMolarMass0MolarMass1, double[] MolarMasses, MaterialLaw EoS, int NumberOfReactants, int SpeciesIndex) {
            MassFractionNames = new string[] { VariableNames.MassFraction0, VariableNames.MassFraction1, VariableNames.MassFraction2, VariableNames.MassFraction3 };
            m_ArgumentOrdering = new string[] { MassFractionNames[SpeciesIndex] };
            this.StoichiometricCoefficients = StoichiometricCoefficients;
            this.ReactionRateConstants = ReactionRateConstants;
            this.SpeciesIndex = SpeciesIndex;
            this.NumberOfReactants = NumberOfReactants ;
            this.OneOverMolarMass0MolarMass1 = OneOverMolarMass0MolarMass1;
            this.MolarMasses = MolarMasses;
            this.EoS = EoS;
        }


        /// <summary>
        /// The current argument, if there is one (i.e. when MassFraction0 or MassFraction1 are being balanced)
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {return m_ArgumentOrdering; }
        }

        /// <summary>
        /// Temperature, MassFraction0, MassFraction1, MassFraction 2, MassFraction 3 at the linearization point.
        /// </summary>
        public override IList<string> ParameterOrdering {
            get { return new string[] { VariableNames.Temperature0, VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 }; }
        }



        protected override double Source(double[] x, double[] parameters, double[] U) {

            double ReactionRate = 0.0;
            double ExponentialPart = ReactionRateConstants[0] * Math.Exp(-ReactionRateConstants[1] / parameters[0]);
            rho = EoS.GetDensity(parameters);
            // 0. MassFraction (fuel) balance species source
            if (SpeciesIndex == 0) {
                // M_alpha/(M_1^a * M_2^b) * Da*exp(-Ta/T)*(rho*Y_f)_(k+1)*(rho*Y_f)_(k)^(a-1)*((rho*Y_o)_(k)^b
        
                ReactionRate = ExponentialPart * OneOverMolarMass0MolarMass1 * rho * U[0] * Math.Pow(rho * parameters[1], ReactionRateConstants[2] - 1) * Math.Pow(rho * parameters[2], ReactionRateConstants[3]);
                Debug.Assert(!double.IsNaN(ReactionRate));
                Debug.Assert(!double.IsInfinity(ReactionRate));
            }
            // 1. MassFraction (oxididizer) balance species source
            else if (SpeciesIndex == 1) {
                // M_alpha/(M_1^a * M_2^b) * Da*exp(-Ta/T)*(rho*Y_f)_(k)^(a)*(rho*Y_o)_(k+1)*(rho*Y_o)_(k)^(b-1)
                ReactionRate = ExponentialPart * OneOverMolarMass0MolarMass1 * Math.Pow(rho * parameters[1], ReactionRateConstants[2]) * rho * U[0] * Math.Pow(rho * parameters[2], ReactionRateConstants[3] - 1);
                Debug.Assert(!double.IsNaN(ReactionRate));
                Debug.Assert(!double.IsInfinity(ReactionRate));
            }
            // product balance species source
            else if (SpeciesIndex > 1 && SpeciesIndex < NumberOfReactants) {
                // M_alpha/(M_1^a * M_2^b) * Da*exp(-Ta/T)*(rho*Y_f)_(k)^a*(rho*Y_o)_(k)^b
                ReactionRate = ExponentialPart * OneOverMolarMass0MolarMass1 * Math.Pow(rho * parameters[1], ReactionRateConstants[2]) * Math.Pow(rho * parameters[2], ReactionRateConstants[3]);
                Debug.Assert(!double.IsNaN(ReactionRate));
                Debug.Assert(!double.IsInfinity(ReactionRate));
            }
            else
                throw new System.ArgumentException("Species index cannot be negative or greater than the number of reactants");
            return - MolarMasses[SpeciesIndex] * StoichiometricCoefficients[SpeciesIndex] * ReactionRate;
        }
    }
}

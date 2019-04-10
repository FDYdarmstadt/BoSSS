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
using System.Diagnostics;

namespace BoSSS.Solution.NSECommon {

    //this class is a stub. It needs to be finished. Currently constant values for heat conductivity, diffusivity and partial and mean heat capacity are being used. 

    /// <summary>
    /// Material law for low Mach number flows with combustion number.
    /// </summary>
    public class MaterialLawCombustion : MaterialLawLowMach {

        MaterialParamsMode MatParamsMode;
        double[] MolarMasses;

  
        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="T_ref">Reference temperature - used in Sutherland's law.</param>
        /// <param name="MatParamsMode">The selected material parameter mode.</param>
        /// <param name="MolarMasses">Array of the molar masses of the fuel, oxidizer and products.</param>
        /// <param name="speciesTransportOK"> Used for switching between a density calculated with and without species transport </param>
        public MaterialLawCombustion(double T_ref, MaterialParamsMode MatParamsMode, double[] MolarMasses)
            : base(T_ref, MatParamsMode) {
            this.MatParamsMode = MatParamsMode;
            this.MolarMasses = MolarMasses;


        }

        /// <summary>
        /// Dimensionless ideal gas law for multicomponent flow - returns density as function of
        /// thermodynamic pressure (i.e. p0), temperature, mass fractions and molar masses.
        /// </summary>
        /// <param name="phi">Temperature, MassFractions_SpeciesIndex</param>
        /// <returns>
        /// Density
        /// </returns>
        public override double GetDensity(params double[] phi) {

            if (IsInitialized) {
                double rho;
                double MassFractionsOverMolarFractions;

                if (phi.Length < 4)
                    throw new ArgumentException("Error in density computation. Number of reactants needs to be atleast 3.");

                MassFractionsOverMolarFractions = 0.0;
                for (int n = 1; n < phi.Length; n++) {
                    MassFractionsOverMolarFractions += phi[n] / MolarMasses[n - 1];
                }
                rho = base.ThermodynamicPressure.Current.GetMeanValue(0) / (phi[0] * MassFractionsOverMolarFractions);
                
                Debug.Assert(!(double.IsNaN(rho) || double.IsInfinity(rho)));
               // rho = 1.0;
                return rho;
            }
            else {
                throw new ApplicationException("ThermodynamicPressure is not initialized.");
            }
        }
        
        public override IList<string> ParameterOrdering {
            get {
                return new string[] { VariableNames.Temperature0 , VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0}; 
            }
        }
        //public double GetHeatConductivity(double phi) {
        //    switch (this.MatParamsMode) {
        //        case MaterialParamsMode.Constant:
        //            return 1.0;
        //        case MaterialParamsMode.Sutherland: {
        //                //    throw new NotImplementedException();
        //                return 1.0; // Using a constant value! 
        //            }
        //        case MaterialParamsMode.PowerLaw: {
        //                throw new NotImplementedException();
        //            }
        //        default:
        //            throw new NotImplementedException();
        //    }
        //}

        //public double GetDiffusivity(double phi) {
        //    switch (this.MatParamsMode) {
        //        case MaterialParamsMode.Constant:
        //            return 1.0;
        //        case MaterialParamsMode.Sutherland: {
        //                //    throw new NotImplementedException();
        //                return 1.0; // Using a constant value! 
        //            }
        //        case MaterialParamsMode.PowerLaw: {
        //                throw new NotImplementedException();
        //            }
        //        default:
        //            throw new NotImplementedException();
        //    }
        //}

        //public double GetPartialHeatCapacity(double phi) {
        //    switch (this.MatParamsMode) {
        //        case MaterialParamsMode.Constant:
        //            return 1.0;
        //        case MaterialParamsMode.Sutherland: {
        //                //    throw new NotImplementedException();
        //                return 1.0; // Using a constant value! 
        //            }
        //        case MaterialParamsMode.PowerLaw: {
        //                throw new NotImplementedException();
        //            }
        //        default:
        //            throw new NotImplementedException();
        //    }
        //}

        //public double GetHeatCapacity(double phi) {
        //    return 1.0;
        //}

    }
}

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
using ilPSP.Utils;
using System.Runtime.Serialization;

namespace BoSSS.Solution.NSECommon {

    //this class is a stub. It needs to be finished. Currently constant values for heat conductivity, diffusivity and partial and mean heat capacity are being used. 

    /// <summary>
    /// Material law for low Mach number flows with combustion number.
    /// </summary>
    public class MaterialLawCombustion : MaterialLawLowMach {
        [DataMember]
        MaterialParamsMode MatParamsMode;
        [DataMember]
        double[] MolarMasses;
        [DataMember]
        bool rhoOne;
  
        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="T_ref">Reference temperature - used in Sutherland's law.</param>
        /// <param name="MatParamsMode">The selected material parameter mode.</param>
        /// <param name="MolarMasses">Array of the molar masses of the fuel, oxidizer and products.</param>
        public MaterialLawCombustion(double T_ref, double[] MolarMasses, bool rhoOne, MaterialParamsMode MatParamsMode = MaterialParamsMode.Constant)
            : base(T_ref, MatParamsMode, rhoOne) {
            this.MatParamsMode = MatParamsMode;
            this.MolarMasses = MolarMasses;
            this.rhoOne = rhoOne;
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
                if (!rhoOne) {
                    if (phi.Length < 4)
                        throw new ArgumentException("Error in density computation. Number of reactants needs to be atleast 3.");

                    MassFractionsOverMolarFractions = 0.0;
                    double lastMassFract = 1.0; // The last mass fraction is calculated here, because the other ones are given as arguments and not as parameters.
                    for (int n = 1; n < phi.Length; n++) {
                        lastMassFract -= phi[n]; // Mass fraction calculated as Yn = 1- sum_i^n-1(Y_i);
                    }

                    phi = ArrayTools.Cat(phi, lastMassFract);

                    for (int n = 1; n < phi.Length; n++) {
                        MassFractionsOverMolarFractions += phi[n] / MolarMasses[n - 1];
                    }

                    rho = ThermodynamicPressureValue / (phi[0] * MassFractionsOverMolarFractions);
                    Debug.Assert(!(double.IsNaN(rho) || double.IsInfinity(rho)));
                } else {
                    rho = 1.0;
                }
                return rho;
            } else {
                throw new ApplicationException("ThermodynamicPressure is not initialized.");
            }
        }
        
        public override IList<string> ParameterOrdering {
            get {
                return new string[] { VariableNames.Temperature0 , VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 }; 
            }
        }

        /// <summary>
        /// Calculates heat capacity of air for a given temperature and composition.
        /// Constant pressure heat capacity (c_p), in [KJ/Kg K] 
        /// </summary>
        /// <param name="T">Temperature, in Kelvin</param>
        /// <param name="Y_N2_air">Mass fraction of nitrogen in air</param>
        /// <param name="Y_O2_air">Mass fraction of oxygen in air</param>
        /// <returns>Constant pressure heat capacity (c_p), in [KJ/Kmol K] </returns>
        public double get_cp_air(double T, double Y_N2_air, double Y_O2_air) {
            double cp = 0;
            double R = 8.314; // Gas constant, J/mol K 
            // Cp correlations for O2 and N2. The polynomials used are the NASA Polynomials for cp. http://combustion.berkeley.edu/gri-mech/version30/files30/thermo30.dat
            // For each component two polynomials are given. One for temperatures 300<T<1000 and other for 1000<T<5000
            // N2 Hot
            double[] cp_N2COLD = new double[] { +0.32986769E+01, +0.14082404E-02, -0.39632223E-05, +0.56415153E-08, -0.24448549E-11 };
            double[] cp_O2COLD = new double[] { +0.32129364E+01, +0.11274863E-02, -0.57561505E-06, +0.13138772E-08, -0.87685539E-12 };

            double[] cp_N2HOT = new double[] { +0.29266400E+01, +0.14879767E-02, -0.56847608E-06, +0.10097038E-09, -0.67533513E-14 };
            double[] cp_O2HOT = new double[] { +0.36975782E+01, +0.61351969E-03, -0.12588419E-06, +0.17752815E-10, -0.11364353E-14 };

            if (T <= 1000) {
                double cp_N2 = cp_N2COLD[0] + cp_N2COLD[1] * Math.Pow(T, 1) + cp_N2COLD[2] * Math.Pow(T, 2) + cp_N2COLD[3] * Math.Pow(T, 3) + cp_N2COLD[4] * Math.Pow(T, 4);
                double cp_O2 = cp_O2COLD[0] + cp_O2COLD[1] * Math.Pow(T, 1) + cp_O2COLD[2] * Math.Pow(T, 2) + cp_O2COLD[3] * Math.Pow(T, 3) + cp_O2COLD[4] * Math.Pow(T, 4);
                cp = Y_N2_air * cp_N2 + Y_O2_air * cp_O2;
                if (T < 200) {
                    throw new Exception("Temperature value out of bounds for calculating cp!");
                }
            } else if (T > 1000) {
                double cp_N2 = cp_N2HOT[0] + cp_N2HOT[1] * Math.Pow(T, 1) + cp_N2HOT[2] * Math.Pow(T, 2) + cp_N2HOT[3] * Math.Pow(T, 3) + cp_N2HOT[4] * Math.Pow(T, 4);
                double cp_O2 = cp_O2HOT[0] + cp_O2HOT[1] * Math.Pow(T, 1) + cp_O2HOT[2] * Math.Pow(T, 2) + cp_O2HOT[3] * Math.Pow(T, 3) + cp_O2HOT[4] * Math.Pow(T, 4);
                cp = Y_N2_air * cp_N2 + Y_O2_air * cp_O2;
                if (T > 5000) {
                    throw new Exception("Temperature value out of bounds for calculating cp!");
                }
            } else {
                throw new Exception("äh?");
            }
            double MW_O2 = 32;
            double MW_N2 = 28;
            double avgMW = 1 / (Y_O2_air / MW_O2 + Y_N2_air/MW_N2);
            cp *= R / avgMW; // cp in KJ/Kg K
            return cp; // 
        }



        /// <summary>
        /// Calculates the average molecular weigth of a mixture. 
        /// </summary>
        /// <param name="Mws">Array with molecular weigths</param>
        /// <param name="Ys">Mass Fractions of the mixture</param>
        /// <returns></returns>
        public double getAvgMW(double[] Mws, double[] Ys) {
            double AvgMw;
            double arg = 0;
            for(int i = 0; i < Mws.Length; i++) {
                arg += Ys[i] / Mws[i];
            }
            AvgMw = 1.0 / arg;


            return AvgMw;
        }

    }
}

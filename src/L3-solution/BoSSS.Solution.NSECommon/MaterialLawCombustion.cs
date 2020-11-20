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
using ilPSP;

namespace BoSSS.Solution.NSECommon {

    //this class is a stub. It needs to be finished. Currently constant values for heat conductivity, diffusivity and partial and mean heat capacity are being used. 

    /// <summary>
    /// Material law for low Mach number flows with combustion number.
    /// </summary>
    public class MaterialLawCombustion : MaterialLawLowMach {
        [DataMember] MaterialParamsMode MatParamsMode;
        [DataMember] double[] MolarMasses;
        [DataMember] bool rhoOne;
        [DataMember] bool cpOne;
        [DataMember] public double Q;
        [DataMember] public double TO0;
        [DataMember] public double TF0;
        [DataMember] public double YF0;
        [DataMember] public double YO0;
        [DataMember] public double zst;
        //[DataMember] public double cp;
        [DataMember] public ChemicalConstants CC;
        public double s;
        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="T_ref">Reference temperature - used in Sutherland's law.</param>
        /// <param name="MatParamsMode">The selected material parameter mode.</param>
        /// <param name="rhoOne"></param>
        /// <param name="Q"></param>
        /// <param name="TO0"></param>
        /// <param name="TF0"></param>
        /// <param name="YO0"></param>
        /// <param name="YF0"></param>
        /// <param name="zst"></param>
        /// <param name="CC"></param>
        /// <param name="Prandtl"></param>
        /// <param name="MolarMasses">Array of the molar masses of the fuel, oxidizer and products.</param>
        public MaterialLawCombustion(double T_ref, double[] MolarMasses, MaterialParamsMode MatParamsMode, bool rhoOne, bool _cpOne, double Q, double TO0, double TF0, double YO0, double YF0, double zst, ChemicalConstants CC, double Prandtl)
            : base(T_ref, MatParamsMode, rhoOne, Prandtl) {
            this.MatParamsMode = MatParamsMode;
            this.Prandtl = Prandtl;
            this.MolarMasses = MolarMasses;
            this.rhoOne = rhoOne;
            this.cpOne = _cpOne;
            this.Q = Q;
            this.TO0 = TO0;
            this.TF0 = TF0;
            this.YF0 = YF0;
            this.YO0 = YO0;
            this.zst = zst;
            //this.cp = 1.0;
            this.CC = CC;
            this.s = (CC.nu_O2 * CC.MW_O2) / (CC.nu_CH4 * CC.MW_CH4);

            this.thermoProperties = new ThermodynamicalProperties();
        }

        public ThermodynamicalProperties thermoProperties;

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
            
                if (!rhoOne) {
                    if (phi.Length < 5)
                       throw new ArgumentException("Error in density computation. Number of reactants needs to be atleast 4.");

                    double MassFractionsOverMolarFractions = 0.0;
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
                    Debug.Assert((rho > 0.0));
                } else {
                    rho = 1.0;
                }
                return rho;
            } else {
                throw new ApplicationException("ThermodynamicPressure is not initialized.");
            }
        }
        /// <summary>
        /// Calculation of the heat capacity of the mixture. 
        /// The first element of the argument corresponds to the adimensional temperature
        /// the next N arguments are the mass fractions.
        /// </summary>
        /// <param name="arguments"></param>
        /// <returns></returns>
        public override double GetMixtureHeatCapacity(double[] arguments) {
            double cp;
            if (!cpOne) {

                double lastMassFract = 1.0; // The last mass fraction is calculated here, because the other ones are given as arguments and not as parameters.
                for (int n = 1; n < arguments.Length; n++) {
                    lastMassFract -= arguments[n]; // Mass fraction calculated as Yn = 1- sum_i^n-1(Y_i);
                }

                arguments = ArrayTools.Cat(arguments, lastMassFract);
                double TRef = 300;
                double DimensionalTemperature = arguments[0] * TRef;
                double[] massFractions = arguments.Skip(1).Take(arguments.Length - 1).ToArray();
                string[] names = new string[] { "CH4", "O2", "CO2", "H2O", "N2" };
                cp = thermoProperties.Calculate_Cp_Mixture(massFractions, names, DimensionalTemperature);

            
            } else {
                cp = 1.0;
            }


            return cp;
        }




        public override IList<string> ParameterOrdering {
            get {
                return new string[] { VariableNames.Temperature0 , VariableNames.MassFraction0_0, VariableNames.MassFraction1_0, VariableNames.MassFraction2_0, VariableNames.MassFraction3_0 }; 
            }
        }


        /// <summary>
        /// Calculates local mixture fraction
        /// </summary>     
        /// <returns></returns>
        virtual public double getMixtureFraction(double YF, double YO ) {
            double Z = (s * YF - YO + YO0) / (s * YF0 + YO0);
            return Z;
        }

        /// <summary>
        /// Calculates the global equivalence ratio 
        /// </summary>     
        /// <returns></returns>
        virtual public double getGlobalEquivalenceRatio(double Yf0, double Yox0) {
            Debug.Assert(!(Yf0 < 0));
            Debug.Assert(!(Yox0 < 0));
            double phi = s * Yf0 / Yox0;
            return phi;
        }



        /// <summary>
        /// Calculates the global equivalence ratio 
        /// </summary>     
        /// <returns></returns>
        virtual public double getLocalEquivalenceRatio(double Yf, double Yox) {
            Debug.Assert(!(Yf < -1e3));
            Debug.Assert(!(Yox < -1e3));
            double Z = getMixtureFraction(Yf, Yox);
            double phi = s * (YF0 / YO0) * (Z / (1.0 - Z));
            //Debug.Assert( phi >= -1e-1);
            if (phi.IsNaNorInf()) {
                phi = double.MaxValue;
            }
            return phi;
        }

        /// <summary>
        /// Calculates local activation energy based on the one-Step model from Fernandez Tarrazo, E, A Sanchez, A Linan, and F Williams. “A Simple One-Step Chemistry Model for 
        /// Partially Premixed Hydrocarbon Combustion.” Combustion and Flame 147, no. 1–2 (October 2006): 32–38. https://doi.org/10.1016/j.combustflame.2006.08.001.
        /// </summary>
        /// <param name="Yf"></param>
        /// <param name="Yox"></param>
        /// <returns></returns>
        public double getTa(double Yf, double Yox) {
            double phi = getLocalEquivalenceRatio(Yf, Yox);
            double Ta;
            if (phi <= 0.64) {
                phi = phi < 0 ? 0.0 : phi;
                Ta = 1.0 + 8.250 * (phi - 0.64) * (phi - 0.64);
            } else if (phi > 1.07) {
                //phi = phi < 0 ? 0.0 : phi;
                Ta = 1.0 + 1.443 * (phi - 1.07) * (phi - 1.07);
            } else {
                Ta = 1.0;
            }
            double Ta0 = 15900; 
            Ta *= Ta0;
            return Ta;
        }




        /// <summary>
        /// Calculates local heat release based on the one-Step model from Fernandez Tarrazo, E, A Sanchez, A Linan, and F Williams. “A Simple One-Step Chemistry Model for 
        /// Partially Premixed Hydrocarbon Combustion.” Combustion and Flame 147, no. 1–2 (October 2006): 32–38. https://doi.org/10.1016/j.combustflame.2006.08.001.
        /// </summary>
        /// <param name="Yf"></param>
        /// <param name="Yox"></param>
        /// <returns></returns>
        public double getHeatRelease(double Yf, double Yox) {
            double phi = getLocalEquivalenceRatio(Yf, Yox);
            double q;
            double q0 = 50100; // Heat release per KG fuel, [kJ/kg]
            if (phi <= 1.0) {
                q = 1.0;
            } else {
                q = 1.0 - 0.21 * (phi - 1);
            }
            q *= q0;

         //   q = q < 0 ? 0.0 : q;// Accept only positive or zero value for q.
            return q;
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

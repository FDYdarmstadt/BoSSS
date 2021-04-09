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
using System;
using System.Runtime.Serialization;

namespace BoSSS.Solution.NSECommon {


    /// <summary>
    /// Chemical parameters. 
    /// All parameters are written in SI units (Kg, m, s, KJ, KMol)
    /// </summary>
    [DataContract]
    [Serializable]
    public class ChemicalConstants {
        public ChemicalConstants() {
        }


   
        /// <summary>
        /// Molar mass of Methane, in Kg/Kmol
        /// </summary>
        // Molar Masses
        [DataMember]
        public double MW_CH4 = 16;
        /// <summary>
        /// Molar mass of oxygen, in Kg/Kmol
        /// </summary>
        [DataMember]
        public double MW_O2 = 32;
        /// <summary>
        /// Molar mass of carbon dioxyd , in Kg/Kmol
        /// </summary>
        [DataMember]
        public double MW_CO2 = 44;
        /// <summary>
        /// Molar mass of water, in Kg/Kmol
        /// </summary>
        [DataMember]
        public double MW_H2O = 18;
        /// <summary>
        /// Molar mass of nytrogen, in Kg/Kmol
        /// </summary>
        [DataMember]
        public double MW_N2 = 28;

        /// <summary>
        /// Stoichiometric coefficient in one-step kinetic combustion
        /// </summary>
        [DataMember]
        public double nu_CH4 = -1;

        /// <summary>
        /// Stoichiometric coefficient in one-step kinetic combustion
        /// </summary>
        [DataMember]
        public double nu_O2 = -2;

        /// <summary>
        /// Stoichiometric coefficient in one-step kinetic combustion
        /// </summary>
        [DataMember]
        public double nu_CO2 = 1;

        /// <summary>
        /// Stoichiometric coefficient in one-step kinetic combustion
        /// </summary>
        [DataMember]
        public double nu_H2O = 2;

        /// <summary>
        /// Average molecular weight from air for 23%O_2, 77%N_2
        /// </summary>
        [DataMember]
        public double PM_Air = 28.8288; // kg air/ kMol air  

        /// <summary>
        /// Ideal gas constant,  m3⋅Pa⋅K−1⋅mol−1.
        /// </summary>
        [DataMember]
        public double R_gas = 8.314;// Ideal gas Constant (R), in m3⋅Pa⋅K−1⋅mol−1.

        /// <summary>
        /// Arrhenius preExponential factor of the one-Step combustion of hydrocarbon.
        /// </summary>
        [DataMember]
        public double PreExponentialFactor = 6.9e11; // Pre-exponential factor,  6.9e14 cm3/(mol s) =>  6.9e11m3/(kmol s)
        /// <summary>
        /// Activation temperature of the one-Step combustion of hydrocarbon.
        /// </summary>
        [DataMember]
        public double Ta = 15900; // Activation TEMPERATURE, K

        /// <summary>
        /// Heat release of combustion for a one-step chemistry with (phi &lt; 1), per unit mass
        /// </summary>
        [DataMember]
        public double HeatReleaseMass = 50100; //  KJ/(Kg fuel)

        /// <summary>
        /// Heat release of combustion for a one-step chemistry with (`$phi &lt; 1`$), per molar unit
        /// </summary>
        [DataMember]
        public double HeatReleaseMolar = 802400; //  KJ/(kmol fuel)

        /// <summary>
        /// Calculates the average molecular weigth of a mixture. 
        /// </summary>
        /// <param name="Mws">Array with molecular weigths</param>
        /// <param name="Ys">Mass Fractions of the mixture</param>
        /// <returns></returns>
        public double getAvgMW(double[] Mws, double[] Ys) {
            double AvgMw;
            double arg = 0;
            for (int i = 0; i < Mws.Length; i++) {
                arg += Ys[i] / Mws[i];
            }
            AvgMw = 1.0 / arg;


            return AvgMw;
        }

    }



}

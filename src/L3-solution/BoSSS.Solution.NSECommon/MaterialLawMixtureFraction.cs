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
using System.Runtime.Serialization;
using System.Text;
using BoSSS.Foundation;
using ilPSP;
using NUnit.Framework;

namespace BoSSS.Solution.NSECommon {

    [DataContract]
    [Serializable]
    public class MaterialLawLowMach_MF : MaterialLawCombustion {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="T_ref"> Reference temperature of Sutherland Law </param>
        /// <param name="MolarMasses">Array of molar masses </param>
        /// <param name="MatParamsMode">Material law (constant, sutherland, etc) </param>
        /// <param name="rhoOne">Switch for constant density </param>
        /// <param name="Q"> Adimensionalized heat release </param>
        /// <param name="TO0">Temperature of oxidizer inlet</param>
        /// <param name="TF0">Temperature of fuel inlet</param>
        /// <param name="YO0">Oxigen mass fraction of oxidizer inlet</param>
        /// <param name="YF0">Fuel mass fraction of fuel inlet </param>
        /// <param name="zst">Stoichiometric mixture fraction </param>
        /// <param name="CC"></param>
        /// <param name="Prandtl"></param>
        /// 
        public MaterialLawLowMach_MF(double T_ref, double[] MolarMasses, MaterialParamsMode MatParamsMode, bool rhoOne, double Q, double TO0, double TF0, double YO0, double YF0, double zst, ChemicalConstants CC, double Prandtl) : base(T_ref, MolarMasses,  MatParamsMode,  rhoOne,  Q,  TO0,  TF0,  YO0,  YF0,  zst,  CC, Prandtl) { 
            this.Q = Q;
            this.TO0 = TO0;
            this.TF0 = TF0;
            this.YF0 = YF0;
            this.YO0 = YO0;
            this.zst = zst;
            this.cp = 1.0;
            this.CC = CC;
            this.Prandtl = Prandtl;
            this.MatParamsMode = MatParamsMode;
            this.rhoOne = rhoOne;
            this.s = (CC.nu_O2 * CC.MW_O2) / (CC.nu_CH4 * CC.MW_CH4);

            

        }
        [DataMember] MaterialParamsMode MatParamsMode;
        [DataMember] public bool rhoOne;
        
        //[DataMember] public double Q;
        //[DataMember] public double TO0;
        //[DataMember] public double TF0;
        //[DataMember] public double YF0;
        //[DataMember] public double YO0;
        //[DataMember] public double zst;
        //[DataMember] public double cp;
        //[DataMember] public ChemicalConstants CC;
        //[DataMember] public double s;

        /// <summary>
        /// Calculate density based on the mixture fraction.
        /// </summary>
        /// <param name="Z"></param>
        /// <returns></returns>
        public override double getDensityFromZ(double Z) {
            double res;
          //  Z = repairMixtureFractionValue(Z);
            if (Q > 0) {
                if (!rhoOne) {
                    //Debug.Assert(Z - 1.0 < 1e-4 && Z > -1e-4);
                    double T, Y0, Y1, Y2, Y3, Y4;
                    if (Z >= zst) { // Fuel side
                        T = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * zst * (1 - Z) / (1 - zst);
                        Y0 = YF0 * (Z - zst) / (1 - zst);
                        Y1 = 0;
                        Y2 = -YO0 * (CC.nu_CO2 * CC.MW_CO2) / (CC.nu_O2 * CC.MW_O2) * (1 - Z);
                        Y3 = -YO0 * (CC.nu_H2O * CC.MW_H2O) / (CC.nu_O2 * CC.MW_O2) * (1 - Z);
                        Y4 = (1.0 - YF0) * (1 - Z) + (1.0 - YO0) * Z;
                    } else if (Z < zst) { // Oxydizer side
                        T = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * Z;
                        Y0 = 0;
                        Y1 = YO0 * (1 - Z / zst);
                        Y2 = -YF0 * (CC.nu_CO2 * CC.MW_CO2) / (CC.nu_CH4 * CC.MW_CH4) * Z;
                        Y3 = -YF0 * (CC.nu_H2O * CC.MW_H2O) / (CC.nu_CH4 * CC.MW_CH4) * Z;
                        Y4 = (1.0 - YF0) * (1 - Z) + (1.0 - YO0) * Z;
                    } else {
                        throw new Exception("out of bounds");
                    }
                    Debug.Assert(Math.Abs(1.0 - (Y0 + Y1 + Y2 + Y3 + Y4)) <= 1e-1);
                    double[] densityArguments = new double[] { T, Y0, Y1, Y2, Y3/*, Y4*/ }; // Y4 is calculated internally in the GetDensity method
                    res = base.GetDensity(densityArguments);
                } else {
                    res = 1.0;
                }
            } else {
                res = base.GetDensity(new double[] { 1.0, Z, 1.0 - Z, 0.0, 0.0 });
            }
            return res;
        }

        /// <summary>
        /// Forces the mixture fraction value to be between 0 and 1. 
        /// </summary>
        /// <returns></returns>
        public double repairMixtureFractionValue(double phi) {
            double repairedPhi;
            if (phi < 0.0) {
                repairedPhi = 0.0;
            } else if (phi > 1.0) {
                repairedPhi = 1.0;
            } else {
                repairedPhi = phi;
            }
            return repairedPhi;
        }

        /// <summary>
        ///  The heat conductivity $\lambda. Possibly dependent on the variable <see cref="phi"/> representing the temperature
        ///  Used for the mixture fraction equations
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public override double GetHeatConductivity(double phi) {
            double Temperature = this.getVariableFromZ(phi, VariableNames.Temperature);
            return base.GetViscosity(Temperature);

        }
        /// <summary>
        ///  The diffusivity  D. Possibly dependent on the variable <see cref="phi"/> representing the temperature
        ///  Used for the mixture fraction equations
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public override double GetDiffusivity(double phi) {
            double Temperature = this.getVariableFromZ(phi, VariableNames.Temperature);
            return base.GetViscosity(Temperature);

        }

        /// <summary>
        /// Dimensionless Sutherland's law.
        /// </summary>
        /// <param name="phi">Temperature</param>
        /// <returns>
        /// Dynamic viscosity
        /// </returns>
        public override double GetViscosity(double phi) {
             double Temperature = getVariableFromZ(phi, VariableNames.Temperature); 
            return base.GetViscosity(Temperature);
        }




        /// <summary>
        /// Calculates the global equivalence ratio 
        /// </summary>     
        /// <returns></returns>
        override public double getLocalEquivalenceRatio(double Yf, double Yox) {
            Debug.Assert(!(Yf < 0));
            Debug.Assert(!(Yox < 0));
            double Z = getMixtureFraction(Yf, Yox);
            double phi = s * (YF0 / YO0) * (Z / (1.0 - Z));
            Debug.Assert(phi >= -1e-1);
            if (phi.IsNaNorInf()) {
                phi = double.MaxValue;
            }
            return phi;
            
        }


        /// <summary>
        /// Calculates local mixture fraction
        /// </summary>     
        /// <returns></returns>
        override public double getMixtureFraction(double YF, double YO) {
            double Z = (s * YF - YO + YO0) / (s * YF0 + YO0);
            return Z;
        }

        /// <summary>
        /// Calculates the global equivalence ratio 
        /// </summary>     
        /// <returns></returns>
        override public double getGlobalEquivalenceRatio(double Yf0, double Yox0) {
            Debug.Assert(!(Yf0 < 0));
            Debug.Assert(!(Yox0 < 0));
            double phi = s * Yf0 / Yox0;
            return phi;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="Z"></param>
        /// <param name="id"></param>
        /// <returns></returns>
        public override double getVariableFromZ(double Z, string id) {
            double res;
            Z = repairMixtureFractionValue(Z);
            if (Z >= zst) { // Fuel side
                switch (id) {
                    case VariableNames.Temperature:
                        res = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * zst * (1 - Z) / (1 - zst);
                        break;
                    case VariableNames.MassFraction0:
                        res = YF0 * (Z - zst) / (1 - zst);
                        break;
                    case VariableNames.MassFraction1:
                        res = 0;
                        break;
                    case VariableNames.MassFraction2:
                        res = -YO0 * (CC.nu_CO2 * CC.MW_CO2) / (CC.nu_O2 * CC.MW_O2) * (1 - Z);
                        break;
                    case VariableNames.MassFraction3:
                        res = -YO0 * (CC.nu_H2O * CC.MW_H2O) / (CC.nu_O2 * CC.MW_O2) * (1 - Z);
                        break;
                    case VariableNames.MassFraction4:
                        double YNOxi0 = 1.0 - YO0;
                        double YNFuel0 = 1.0 - YF0;
                        res = YNOxi0 * (1 - Z) + YNFuel0 * Z;
                        break;
                    default:
                        throw new NotImplementedException("Variable " + id + " cannot be derived from mixture Fraction");
                }

            } else if (Z < zst) { // Oxydizer side
                switch (id) {
                    case VariableNames.Temperature:
                        res = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * Z;
                        break;
                    case VariableNames.MassFraction0:
                        res = 0;
                        break;
                    case VariableNames.MassFraction1:
                        res = YO0 * (1 - Z / zst);
                        break;
                    case VariableNames.MassFraction2:
                        res = -YF0 * (CC.nu_CO2 * CC.MW_CO2) / (CC.nu_CH4 * CC.MW_CH4) * Z;
                        break;
                    case VariableNames.MassFraction3:
                        res = -YF0 * (CC.nu_H2O * CC.MW_H2O) / (CC.nu_CH4 * CC.MW_CH4) * Z;
                        break;
                    case VariableNames.MassFraction4:
                        double YNOxi0 = 1.0 - YO0;
                        double YNFuel0 = 1.0 - YF0;
                        res = YNOxi0 * (1 - Z) + YNFuel0 * Z;
                        break;
                    default:
                        throw new NotImplementedException("Variable " + id + " cannot be derived from mixture Fraction");
                }
            } else {
                throw new Exception("out of bounds");
            }
            return res;

        }

        public override IList<string> ParameterOrdering {
            get {
                return new string[] {
                    VariableNames.Temperature0,
                    VariableNames.MassFraction0_0,
                    VariableNames.MassFraction1_0,
                    VariableNames.MassFraction2_0,
                    VariableNames.MassFraction3_0                     
                    //,VariableNames.MassFraction4_0
                };
            }
        }

        //public double Prandtl { get; }
    }
 
}

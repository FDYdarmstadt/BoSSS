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
using System.Runtime.Serialization;

namespace BoSSS.Solution.NSECommon {

    [DataContract]
    [Serializable]
    public class MaterialLawMixtureFractionNew : MaterialLaw_MultipleSpecies {
        [DataMember] public ChemicalConstants CC;

        [DataMember] public double cp;

        //[DataMember] public double Q;

        [DataMember] public bool rhoOne;

        [DataMember] public double s;

        [DataMember] public double TF0;

        [DataMember] public double TO0;

        [DataMember] public double YF0;

        [DataMember] public double YO0;

        [DataMember] public double zst;

        [DataMember] private MaterialParamsMode MatParamsMode;

        /// <summary>
        ///
        /// </summary>
        /// <param name="T_ref_Sutherland"> Reference temperature of Sutherland Law </param>
        /// <param name="MolarMasses">Array of molar masses </param>
        /// <param name="MatParamsMode">Material law (constant, sutherland, etc) </param>
        /// <param name="rhoOne">Switch for constant density </param>
        /// <param name="gasConstant"></param>
        /// <param name="Q"> Adimensionalized heat release , if equal to zero, the unburnt mixture fraction profile will be calculated</param>
        /// <param name="TO0">Temperature of oxidizer inlet</param>
        /// <param name="TF0">Temperature of fuel inlet</param>
        /// <param name="YO0">Oxigen mass fraction of oxidizer inlet</param>
        /// <param name="YF0">Fuel mass fraction of fuel inlet </param>
        /// <param name="zst">Stoichiometric mixture fraction </param>
        /// <param name="CC"></param>
        ///
        /// <param name="Prandtl"></param>
        public MaterialLawMixtureFractionNew(double T_ref_Sutherland, double[] MolarMasses, MaterialParamsMode MatParamsMode, bool rhoOne, double gasConstant, double Q, double TO0, double TF0, double YO0, double YF0, double zst, ChemicalConstants CC, OneStepChemicalModel chemModel, double cpRef, double _smoothingFactor) : base(MolarMasses, MatParamsMode, rhoOne, gasConstant, T_ref_Sutherland, chemModel, _cpRef: 1.0, mycpMode: CpCalculationMode.mixture) {
            base.Q = Q;
            this.CC = CC;

            this.TO0 = TO0;
            this.TF0 = TF0;
            this.YF0 = YF0;
            this.YO0 = YO0;
            this.zst = zst;

            this.cp = cpRef;// GetCpFlameSheet(7.0); // Dimensional value

            this.MatParamsMode = MatParamsMode;


            this.rhoOne = rhoOne;
            this.s = (CC.nu_O2 * CC.MW_O2) / (CC.nu_CH4 * CC.MW_CH4);
            this.SmoothingFactor = _smoothingFactor;
            if(this.SmoothingFactor > 0) {
                Console.WriteLine("The mixture fraction is smoothed with a factor" + SmoothingFactor);
            } else {
                Console.WriteLine("The mixture fraction is not being smoothed");
            }
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

        /// <summary>
        /// Calculate density based on the mixture fraction.
        /// </summary>
        /// <param name="Z"></param>
        /// <returns></returns>
        public override double getDensityFromZ(double Z) {
            double res;
            if (!rhoOne) {
                double T, Y0, Y1, Y2, Y3, Y4;
                T = getVariableFromZ(Z, VariableNames.Temperature);
                Y0 = getVariableFromZ(Z, VariableNames.MassFraction0);
                Y1 = getVariableFromZ(Z, VariableNames.MassFraction1);
                Y2 = getVariableFromZ(Z, VariableNames.MassFraction2);
                Y3 = getVariableFromZ(Z, VariableNames.MassFraction3);
                Y4 = getVariableFromZ(Z, VariableNames.MassFraction4);

                Debug.Assert(Math.Abs(1.0 - (Y0 + Y1 + Y2 + Y3 + Y4)) <= 1e-1);
                double[] densityArguments = new double[] { T, Y0, Y1, Y2, Y3/*, Y4*/ }; // Y4 is calculated internally in the GetDensity method
                res = base.GetDensity(densityArguments);
            } else {
                res = ConstantDensityValue;
            }
            return res;
        }
        /// <summary>
        /////
        ///// </summary>
        //public double ConstantDensityValue { get; set; } = 1.0;
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
        ///
        /// </summary>
        /// <param name="Z"></param>
        /// <param name="id"></param>
        /// <returns></returns>
        public override double getVariableFromZ(double Z, string id) {
            double res;
            Z = repairMixtureFractionValue(Z); // apparently really helps...
            ////////////////////////////////////////////////////////////
            /////////// Irreversible fast chemistry (with chemical reaction)
            ////////////////////////////////////////////////////////////

            bool useSmoothFunctionForMF = SmoothingFactor > 0;

            if (Q > 0) {
                if (!useSmoothFunctionForMF) {
                    if (Z >= zst) { // Fuel side
                        res = getVariableFromZ_FuelSide(Z, id);
                    } else { // Oxydizer side
                        res = getVariableFromZ_AirSide(Z, id);
                    }
                } else {
                    double K = SmoothingFactor;
                    double greater = getVariableFromZ_FuelSide(Z, id);
                    double smaller = getVariableFromZ_AirSide(Z, id);

                    //double bla = Math.Tanh(K * (Z - zst));
                    res = smaller + (1 + Math.Tanh(K * (Z - zst))) / 2 * (greater - smaller);
                }
            } else {
                res = getVariableFromZ_FrozenLimit(Z, id);
            }
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        [DataMember] public double SmoothingFactor = 10;
        /// <summary>
        ///
        /// </summary>
        /// <param name="Z"></param>
        /// <param name="id"></param>
        /// <returns></returns>
        public double getVariableFromZ_FuelSide(double Z, string id) {
            double res;

            switch (id) {
                case VariableNames.Temperature:
                    res = Z * TF0 + (1 - Z) * TO0 + Q * YF0 /*/ cp*/ * zst * (1 - Z) / (1 - zst);
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
            return res;
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="Z"></param>
        /// <param name="id"></param>
        /// <returns></returns>
        public double getVariableFromZ_AirSide(double Z, string id) {
            double res;
            switch (id) {
                case VariableNames.Temperature:
                    res = Z * TF0 + (1 - Z) * TO0 + Q * YF0 /*/ cp */* Z;
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
            return res;
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="Z"></param>
        /// <param name="id"></param>
        /// <returns></returns>
        public double getVariableFromZ_FrozenLimit(double Z, string id) {
            double res;
            ////////////////////////////////////////////////////////////
            /////////// Frozen limit calculation (no chemical reaction)
            ////////////////////////////////////////////////////////////
            switch (id) {
                case VariableNames.Temperature:
                    res = Z * TF0 + (1 - Z) * TO0;
                    break;

                case VariableNames.MassFraction0:
                    res = YF0 * Z;
                    break;

                case VariableNames.MassFraction1:
                    res = YO0 * (1 - Z);
                    break;

                case VariableNames.MassFraction2:
                    res = 0;
                    break;

                case VariableNames.MassFraction3:
                    res = 0;
                    break;

                case VariableNames.MassFraction4:
                    double YNOxi0 = 1.0 - YO0;
                    double YNFuel0 = 1.0 - YF0;
                    res = YNOxi0 * (1 - Z) + YNFuel0 * Z;
                    break;

                default:
                    throw new NotImplementedException("Variable " + id + " cannot be derived from mixture Fraction");
            }
            return res;
        }
        /// <summary>
        /// calculates iteratively temperature at the flame sheet, 
        /// asuming a variable mixture cp formulation.
        /// </summary>
        /// <param name="t_start">Initial estimate of temperature</param>
        /// <returns></returns>
        public double GetCpFlameSheet(double t_start) {
            double zst = this.zst;

            double YF = getVariableFromZ(zst, VariableNames.MassFraction0);
            double YO = getVariableFromZ(zst, VariableNames.MassFraction1);
            double YP1 = getVariableFromZ(zst, VariableNames.MassFraction2);
            double YP2 = getVariableFromZ(zst, VariableNames.MassFraction3);
            double YN = getVariableFromZ(zst, VariableNames.MassFraction4);


            double res = 1000;
            double temperature = t_start;
            while (res > 1e-7) {                
                double cpit = GetMixtureHeatCapacity(new double[] { temperature, YF, YO, YP1, YP2 });
                double T = TO0 + Q * YF0 / cpit * zst;
                res = Math.Abs(temperature - T);
                temperature = T;
            }

            return GetMixtureHeatCapacity(new double[] { temperature, YF, YO, YP1, YP2 });
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
    }
}
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

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Material law for MultiSpecies low Mach number flows
    /// </summary>
    public class MaterialLaw_MultipleSpecies : MaterialLaw {
        [DataMember] private MaterialParamsMode MatParamsMode;
        [DataMember] private CpCalculationMode myCpCalculationsMode;

        [DataMember] public double R;
        [DataMember] public double[] MolarMasses;
        [DataMember] public double T_RefSutherland;

        /// <summary>
        /// Heat of Reaction
        /// </summary>
        [DataMember] public double Q;

        [DataMember]
        protected bool rhoOne;



        [DataMember]
        private double cpRef;

        [DataMember]
        public OneStepChemicalModel m_ChemModel;

        /// <summary>
        /// Ctor.
        /// </summary>
        public MaterialLaw_MultipleSpecies(double[] MolarMasses, MaterialParamsMode MatParamsMode, bool _rhoOne, double gasConstant, double TRefSutherland, OneStepChemicalModel chemModel, double _cpRef, CpCalculationMode mycpMode) {
            this.MatParamsMode = MatParamsMode;
            this.R = gasConstant;
            this.thermoProperties = new ThermodynamicalProperties();
            this.MolarMasses = MolarMasses;
            this.T_RefSutherland = TRefSutherland;
            this.rhoOne = _rhoOne;

            this.cpRef = _cpRef;
            m_ChemModel = chemModel;
            myCpCalculationsMode = mycpMode;
        }

        public ThermodynamicalProperties thermoProperties;

        /// <summary>
        /// true if the ThermodynamicPressure is already initialized
        /// </summary>
        public bool IsInitialized {
            get {
                return (ThermodynamicPressure != -1);
            }
        }

        /// <summary>
        /// Hack to initalize ThermodynamicPressure
        /// </summary>
        /// <param name="ThermodynamicPressure"></param>
        public void Initialize(double ThermodynamicPressure) {
            if (!IsInitialized) {
                this.ThermodynamicPressure = ThermodynamicPressure;
            } else {
                throw new ApplicationException("Initialize() can be called only once.");
            }
        }

        public double ThermodynamicPressure = -1.0;

        /// <summary>
        /// Dimensionless ideal gas law for multicomponent flow - returns density as function of
        /// thermodynamic pressure (i.e. p0), temperature, mass fractions and molar masses.
        /// </summary>
        /// <param name="phi">Temperature, MassFractions_SpeciesIndex</param>
        /// <returns>
        /// Density
        /// </returns>
        public override double GetDensity(params double[] phi) {
            return GetDensity(ThermodynamicPressure, phi);
        }

        /// <summary>
        /// Modes for cp calculation
        /// </summary>
        public enum CpCalculationMode {

            /// <summary>
            /// Use nasa correlation of Nitrogen
            /// Useful for combustion systems where both inflows are diluted in inert.
            /// </summary>
            onlyN2,

            /// <summary>
            /// Use nasa correlations for each component
            /// The mixture cp is calculated by asuming an ideal gas, i.e.
            /// cp_mixture = sum_i(cp,i*Y_i)
            /// </summary>
            mixture,

            /// <summary>
            /// Use a constant value
            /// </summary>
            constant
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="p0"></param>
        /// <param name="phi"></param>
        /// <returns></returns>
        public double GetDensity(double p0, params double[] phi) {
            if (IsInitialized) {
                double rho;

                if (!rhoOne) {
                    if (phi.Length < 2)
                        throw new ArgumentException("Error in density computation. Number of chemical component needs to be atleast 1.");

                    double MassFractionsOverMolarFractions = 0.0;
                    double lastMassFract = 1.0; // The last mass fraction is calculated here, because the other ones are given as arguments and not as parameters.
                    for (int n = 1; n < phi.Length; n++) {
                        lastMassFract -= phi[n]; // Mass fraction calculated as Yn = 1- sum_i^n-1(Y_i);
                    }

                    phi = ArrayTools.Cat(phi, lastMassFract);
                    for (int n = 1; n < phi.Length; n++) {
                        MassFractionsOverMolarFractions += phi[n] / MolarMasses[n - 1];
                    }

                    rho = p0 / (R * phi[0] * MassFractionsOverMolarFractions);

                    Debug.Assert(!(double.IsNaN(rho) || double.IsInfinity(rho)));
                } else {
                    rho = ConstantDensityValue;
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

            double TRef = 300;

            switch (myCpCalculationsMode) {
                case CpCalculationMode.onlyN2:
                    double DimensionalTemperature = arguments[0] * TRef;
                    cp = thermoProperties.getCp("N2", DimensionalTemperature) / this.cpRef; // Using correlation for N2...                                                                                        
                    break;

                case CpCalculationMode.mixture:
                    double DimensionalTemperature2 = arguments[0] * TRef;
                    double lastMassFract = 1.0; // The last mass fraction is calculated here, because the other ones are given as arguments and not as parameters.
                    for (int n = 1; n < arguments.Length; n++) {
                        lastMassFract -= arguments[n]; // Mass fraction calculated as Yn = 1- sum_i^n-1(Y_i);
                    }
                    arguments = ArrayTools.Cat(arguments, lastMassFract);
                    double[] massFractions = arguments.Skip(1).Take(arguments.Length - 1).ToArray();
                    string[] names = new string[] { "CH4", "O2", "CO2", "H2O", "N2" };
                    double cpMixture = thermoProperties.Calculate_Cp_Mixture(massFractions, names, DimensionalTemperature2) ;


                    cp = cpMixture / cpRef;
                    break;

                case CpCalculationMode.constant:
                    cp = ConstantHeatCapacityValue; // => cpRef/cpRef
                    break;

                default:
                    throw new NotImplementedException();
            }

            return cp;
        }

        public override IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

   

        /// <summary>
        ///
        /// </summary>
        public double ConstantDensityValue { get; set; } = 1.0;

        /// <summary>
        ///
        /// </summary>
        public double ConstantViscosityValue { get; set; } = 1.0;


        /// <summary>
        ///
        /// </summary>
        public double ConstantHeatConductivityValue { get; set; } = 1.0;
        /// <summary>
        ///
        /// </summary>
        public double ConstantDiffusivityFactorVal { get; set; } = 1.0;

        /// <summary>
        ///
        /// </summary>
        public double ConstantHeatCapacityValue { get; set; } = 1.0;



        /// <summary>
        /// Dimensionless Sutherland's law.
        /// </summary>
        /// <param name="phi">Temperature</param>
        /// <returns>
        /// Dynamic viscosity
        /// </returns>
        public override double GetViscosity(double phi) {
            phi = Math.Max(0.01, phi);
            double visc = 0; // nondimensional viscosity
            switch (this.MatParamsMode) {
                case MaterialParamsMode.Constant: {
                        visc = ConstantViscosityValue;
                        break;
                    }
                case MaterialParamsMode.Sutherland: {
                        double S = 110.5;
                        visc = Math.Pow(phi, 1.5) * (1 + S / T_RefSutherland) / (phi + S / T_RefSutherland);
                        break;
                    }
                case MaterialParamsMode.PowerLaw: {
                        double exponent = 2.0 / 3.0;//
                        visc = Math.Pow(phi, exponent);
                        break;
                    }
                default:
                    throw new NotImplementedException();
            }
            if (double.IsNaN(visc) || double.IsInfinity(visc) || visc < 0)
                throw new ArithmeticException("Invalid value for viscosity: " + visc);
            return visc;
        }

        /// <summary>
        ///  The heat conductivity $\lambda. Possibly dependent on the variable <see cref="phi"/> representing the temperature
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public virtual double GetHeatConductivity(double phi) {
            phi = Math.Max(0.01, phi);
            double visc = 0; // nondimensional heat conductivity
            switch (this.MatParamsMode) {
                case MaterialParamsMode.Constant: {
                        visc = ConstantHeatConductivityValue;
                        break;
                    }
                case MaterialParamsMode.Sutherland: {
                        double S = 110.5;
                        visc = Math.Pow(phi, 1.5) * (1 + S / T_RefSutherland) / (phi + S / T_RefSutherland);
                        break;
                    }
                case MaterialParamsMode.PowerLaw: {
                        //double exponent = 0.7;
                        double exponent = 2.0 / 3.0;//
                        visc = Math.Pow(phi, exponent);
                        break;
                    }
                default:
                    throw new NotImplementedException();
            }
            if (double.IsNaN(visc) || double.IsInfinity(visc) || visc < 0)
                throw new ArithmeticException("Invalid value for viscosity: " + visc);
            return visc;
        }


        /// <summary>
        /// The mass diffusivity,D  multiplied by rho.
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public virtual double GetDiffusivity(double phi) {
            phi = Math.Max(0.01, phi);
            double visc = 0; // nondimensional heat conductivity
            switch (this.MatParamsMode) {
                case MaterialParamsMode.Constant: {
                        visc = ConstantDiffusivityFactorVal;
                        break;
                    }
                case MaterialParamsMode.Sutherland: {
                        double S = 110.5;
                        visc = Math.Pow(phi, 1.5) * (1 + S / T_RefSutherland) / (phi + S / T_RefSutherland);
                        break;
                    }
                case MaterialParamsMode.PowerLaw: {
                        //double exponent = 0.7;
                        double exponent = 2.0 / 3.0;//
                        visc = Math.Pow(phi, exponent);
                        break;
                    }
                default:
                    throw new NotImplementedException();
            }
            if (double.IsNaN(visc) || double.IsInfinity(visc) || visc < 0)
                throw new ArithmeticException("Invalid value for viscosity: " + visc);
            return visc;
        }


        /// <summary>
        /// Returns thermodynamic pressure as function of inital mass and temperature.
        /// </summary>
        /// <param name="InitialMass"></param>
        /// <param name="Temperature"></param>
        /// <returns></returns>
        public double GetMassDeterminedThermodynamicPressure(double InitialMass, XDGField Temperature) {
            var TemperatureA = Temperature.GetSpeciesShadowField("A");

            var basis = new Basis(TemperatureA.GridDat, 23); // Degree
            var omega = new SinglePhaseField(TemperatureA.Basis); // 1 / T
            //var omega = new XDGField(Temperature.Basis); // 1 / T
            omega.ProjectField(1.0,
                delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes

                    MultidimensionalArray temp = MultidimensionalArray.Create(Len, K);

                    TemperatureA.Evaluate(j0, Len, NS, temp);
                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {
                            result[j, k] = 1 / temp[j, k];
                        }
                    }
                }, new Foundation.Quadrature.CellQuadratureScheme(true, null));

            double p0 = (InitialMass / omega.IntegralOver(null));

            //this.ThermodynamicPressure2.IncreaseHistoryLength(3);
            //this.ThermodynamicPressure2.Current.Clear();
            //this.ThermodynamicPressure2.Current.AccConstant(p0);
            //this.ThermodynamicPressure2.Push();

            //Console.WriteLine("push count : " + ThermodynamicPressure2.PushCount);
            //Console.WriteLine("history length: " + ThermodynamicPressure2.HistoryLength);

            //this.ThermodynamicPressureValue = p0;
            //Console.WriteLine(p0);
            return p0;
        }

        public override double GetMassDeterminedThermodynamicPressure(double InitialMass, SinglePhaseField Temperature) {
            throw new NotImplementedException();
        }

        public override double getVariableFromZ(double Z, string id) {
            throw new NotImplementedException();
        }

        public override double getDensityFromZ(double Z) {
            throw new NotImplementedException();
        }

        public override double GetLambda(double[] VelocityMean, double[] Normal, double ScalarMean) {
            throw new NotImplementedException();
        }

        public override double DiffRho_Temp(double phi) {
            throw new NotImplementedException();
        }

        #region Combustion related utils

        ///// <summary>
        ///// Calculates local mixture fraction
        ///// </summary>
        ///// <returns></returns>
        //virtual public double getMixtureFraction(double YF, double YO ) {
        //    double Z = (s * YF - YO + YO0) / (s * YF0 + YO0);
        //    return Z;
        //}

        ///// <summary>
        ///// Calculates the global equivalence ratio
        ///// </summary>
        ///// <returns></returns>
        //virtual public double getGlobalEquivalenceRatio(double Yf0, double Yox0) {
        //    Debug.Assert(!(Yf0 < 0));
        //    Debug.Assert(!(Yox0 < 0));
        //    double phi = s * Yf0 / Yox0;
        //    return phi;
        //}

        ///// <summary>
        ///// Calculates the global equivalence ratio
        ///// </summary>
        ///// <returns></returns>
        //virtual public double getLocalEquivalenceRatio(double Yf, double Yox) {
        //    Debug.Assert(!(Yf < -1e3));
        //    Debug.Assert(!(Yox < -1e3));
        //    double Z = getMixtureFraction(Yf, Yox);
        //    double phi = s * (YF0 / YO0) * (Z / (1.0 - Z));
        //    //Debug.Assert( phi >= -1e-1);
        //    if (phi.IsNaNorInf()) {
        //        phi = double.MaxValue;
        //    }
        //    return phi;
        //}

        ///// <summary>
        ///// Calculates local activation energy based on the one-Step model from Fernandez Tarrazo, E, A Sanchez, A Linan, and F Williams. “A Simple One-Step Chemistry Model for
        ///// Partially Premixed Hydrocarbon Combustion.” Combustion and Flame 147, no. 1–2 (October 2006): 32–38. https://doi.org/10.1016/j.combustflame.2006.08.001.
        ///// </summary>
        ///// <param name="Yf"></param>
        ///// <param name="Yox"></param>
        ///// <returns></returns>
        //public double getTa(double Yf, double Yox) {
        //    double phi = getLocalEquivalenceRatio(Yf, Yox);
        //    double Ta;
        //    if (phi <= 0.64) {
        //        phi = phi < 0 ? 0.0 : phi;
        //        Ta = 1.0 + 8.250 * (phi - 0.64) * (phi - 0.64);
        //    } else if (phi > 1.07) {
        //        //phi = phi < 0 ? 0.0 : phi;
        //        Ta = 1.0 + 1.443 * (phi - 1.07) * (phi - 1.07);
        //    } else {
        //        Ta = 1.0;
        //    }
        //    double Ta0 = 15900;
        //    Ta *= Ta0;
        //    return Ta;
        //}

        ///// <summary>
        ///// Calculates local heat release based on the one-Step model from Fernandez Tarrazo, E, A Sanchez, A Linan, and F Williams. “A Simple One-Step Chemistry Model for
        ///// Partially Premixed Hydrocarbon Combustion.” Combustion and Flame 147, no. 1–2 (October 2006): 32–38. https://doi.org/10.1016/j.combustflame.2006.08.001.
        ///// </summary>
        ///// <param name="Yf"></param>
        ///// <param name="Yox"></param>
        ///// <returns></returns>
        //public double getHeatRelease(double Yf, double Yox) {
        //    double phi = getLocalEquivalenceRatio(Yf, Yox);
        //    double q;
        //    double q0 = 50100; // Heat release per KG fuel, [kJ/kg]
        //    if (phi <= 1.0) {
        //        q = 1.0;
        //    } else {
        //        q = 1.0 - 0.21 * (phi - 1);
        //    }
        //    q *= q0;

        // //   q = q < 0 ? 0.0 : q;// Accept only positive or zero value for q.
        //    return q;
        //}

        #endregion Combustion related utils
    }

    public class OneStepChemicalModel {

        private OneStepChemicalModel() {
        }

        private bool m_variableParameters;
        private double s;
        private double YO0;
        private double YF0;

        public OneStepChemicalModel(bool variableParameters, double _YO0, double _YF0) {
            m_variableParameters = variableParameters;
            s = 4.0;
            YO0 = _YO0;
            YF0 = _YF0;
        }

        /// <summary>
        /// Calculates local mixture fraction
        /// </summary>
        /// <returns></returns>
        virtual public double getMixtureFraction(double YF, double YO) {
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
            double Ta;
            double Ta0 = 15900;
            if (m_variableParameters) {
                double phi = getLocalEquivalenceRatio(Yf, Yox);

                if (phi <= 0.64) {
                    phi = phi < 0 ? 0.0 : phi;
                    Ta = 1.0 + 8.250 * (phi - 0.64) * (phi - 0.64);
                } else if (phi > 1.07) {
                    //phi = phi < 0 ? 0.0 : phi;
                    Ta = 1.0 + 1.443 * (phi - 1.07) * (phi - 1.07);
                } else {
                    Ta = 1.0;
                }

                Ta *= Ta0;
            } else {
                Ta = Ta0;
            }
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
            double q;
            double q0 = 50100; // Heat release per KG fuel, [kJ/kg]
            if (m_variableParameters) {
                double phi = getLocalEquivalenceRatio(Yf, Yox);

                if (phi <= 1.0) {
                    q = 1.0;
                } else {
                    q = 1.0 - 0.21 * (phi - 1);
                }
                q *= q0;
            } else {
                q = q0;
            }
            //   q = q < 0 ? 0.0 : q;// Accept only positive or zero value for q.
            return q;
        }
    }
}
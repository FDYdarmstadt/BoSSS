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
        /// 
        public MaterialLawLowMach_MF(double T_ref, double[] MolarMasses, MaterialParamsMode MatParamsMode, bool rhoOne, double Q, double TO0, double TF0, double YO0, double YF0, double zst, ChemicalConstants CC, double Prandtl) : base(T_ref, MolarMasses, rhoOne, MatParamsMode, Prandtl) {
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
            
        }
        MaterialParamsMode MatParamsMode;
        public bool rhoOne;
        public double Q;
        public double TO0;
        public double TF0;
        public double YF0;
        public double YO0;
        public double zst;
        public double cp;
        public ChemicalConstants CC;


        /// <summary>
        /// Calculate density based on the mixture fraction.
        /// </summary>
        /// <param name="Z"></param>
        /// <returns></returns>
        public override double getDensityFromZ(double Z) {
            double res;
            Z = repairMixtureFractionValue(Z);
            if (Q > 0) {
                if (!rhoOne) {
                    //Debug.Assert(Z - 1.0 < 1e-4 && Z > -1e-4);
                    double T, Y0, Y1, Y2, Y3, Y4;
                    if (Z >= zst) { // Fuel side
                        T = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * zst * (1 - Z) / (1 - zst);
                        Y0 = YF0 * (Z - zst) / (1 - zst);
                        Y1 = 0;
                        Y2 = -YO0 * (CC.s_CO2 * CC.PM_CO2) / (CC.s_O2 * CC.PM_O2) * (1 - Z);
                        Y3 = -YO0 * (CC.s_H2O * CC.PM_H2O) / (CC.s_O2 * CC.PM_O2) * (1 - Z);
                        Y4 = (1.0 - YF0) * (1 - Z) + (1.0 - YO0) * Z;
                    } else if (Z < zst) { // Oxydizer side
                        T = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * Z;
                        Y0 = 0;
                        Y1 = YO0 * (1 - Z / zst);
                        Y2 = -YF0 * (CC.s_CO2 * CC.PM_CO2) / (CC.s_CH4 * CC.PM_CH4) * Z;
                        Y3 = -YF0 * (CC.s_H2O * CC.PM_H2O) / (CC.s_CH4 * CC.PM_CH4) * Z;
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
            return base.GetHeatConductivity(Temperature);
        }
        /// <summary>
        ///  The diffusivity  D. Possibly dependent on the variable <see cref="phi"/> representing the temperature
        ///  Used for the mixture fraction equations
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public override double GetDiffusivity(double phi) {
            double Temperature = this.getVariableFromZ(phi, VariableNames.Temperature);
            return base.GetDiffusivity(Temperature);

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
                        res = -YO0 * (CC.s_CO2 * CC.PM_CO2) / (CC.s_O2 * CC.PM_O2) * (1 - Z);
                        break;
                    case VariableNames.MassFraction3:
                        res = -YO0 * (CC.s_H2O * CC.PM_H2O) / (CC.s_O2 * CC.PM_O2) * (1 - Z);
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
                        res = -YF0 * (CC.s_CO2 * CC.PM_CO2) / (CC.s_CH4 * CC.PM_CH4) * Z;
                        break;
                    case VariableNames.MassFraction3:
                        res = -YF0 * (CC.s_H2O * CC.PM_H2O) / (CC.s_CH4 * CC.PM_CH4) * Z;
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

        public double Prandtl { get; }
    }
    /// <summary>
    /// Material law for low Mach number flows.
    /// </summary>
    [DataContract]
    [Serializable]
    public class MaterialLawLowMach : MaterialLaw {
        [DataMember]
        public double T_ref;

        [DataMember]
        MaterialParamsMode MatParamsMode;

        [DataMember]
        bool rhoOne;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="T_ref">Reference temperature - used in Sutherland's law.</param>
        /// <param name="MatParamsMode"></param>
        /// <param name="rhoOne"></param>
        /// <param name="Prandtl"></param>
        public MaterialLawLowMach(double T_ref, MaterialParamsMode MatParamsMode, bool rhoOne, double Prandtl)
            : base() {
            this.rhoOne = rhoOne;
            this.Prandtl = Prandtl;
            this.T_ref = T_ref;
            this.MatParamsMode = MatParamsMode;
        }

        /// <summary>
        /// 
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return new string[] {};
            }
        }

        /// <summary>
        /// true if the ThermodynamicPressure is already initialized
        /// </summary>
        protected bool IsInitialized {
            get {
                return ThermodynamicPressure != null;
            }
        }

        public double Prandtl { get; }

        /// <summary>
        /// 
        /// </summary>
        //[NonSerialized]
        [DataMember]
        protected ScalarFieldHistory<SinglePhaseField> ThermodynamicPressure;
        /// <summary>
        /// 
        /// </summary>
        //[NonSerialized]
        [DataMember]
        public double ThermodynamicPressureValue = -1;



        /// <summary>
        /// Hack to initalize ThermodynamicPressure. 
        /// </summary>
        /// <param name="ThermodynamicPressure">
        /// Hack for introducing the value of p0 as a double. has to be changed
        /// </param>
        /// <param name="ThermodynamicPressureValue"></param>
        public void Initialize(ScalarFieldHistory<SinglePhaseField> ThermodynamicPressure, ref double ThermodynamicPressureValue) {
            if (!IsInitialized) {
                this.ThermodynamicPressure = ThermodynamicPressure;
                this.ThermodynamicPressureValue = ThermodynamicPressureValue;
            } else {
                throw new ApplicationException("Initialize() can be called only once.");
            }
        }


        /// <summary>
        /// Hack to initalize ThermodynamicPressure - called by NSE_SIMPLE.VariableSet.Initialize()
        /// </summary>
        /// <param name="ThermodynamicPressure"></param>
        public void Initialize(ScalarFieldHistory<SinglePhaseField> ThermodynamicPressure) {
            if (!IsInitialized) {
                this.ThermodynamicPressure = ThermodynamicPressure;
            } else {
                throw new ApplicationException("Initialize() can be called only once.");
            }
        }

        /// <summary>
        /// Dimensionless ideal gas law - returns density as function of
        /// thermodynamic pressure (i.e. p0) and temperature.
        /// </summary>
        /// <param name="phi">Temperature</param>
        /// <returns>
        /// Density
        /// </returns>
        public override double GetDensity(params double[] phi) {
            if (IsInitialized) {            
                double rho = 1.0;
                if (rhoOne) {
                    rho = 1.0;
                    return rho;
                } else {
                    if (ThermodynamicPressureValue != -1) { // this is a really ugly hack to allow the SIMPLE project to use the p0 DG field. A better solution has to be found                                                    
                        rho = ThermodynamicPressureValue / phi[0];
                    } else {
                        rho = ThermodynamicPressure.Current.GetMeanValue(0) / phi[0];
                    }
                }
                Debug.Assert(!double.IsNaN(rho));
                Debug.Assert(!double.IsInfinity(rho));
                return rho;
            } else {
                throw new ApplicationException("ThermodynamicPressure is not initialized.");
            }
        }


        /// <summary>
        /// Sutherlands law for air. 
        /// </summary>
        /// <param name="T"></param>
        /// <returns>
        /// The viscosity of air at a given temperature in Kg/(m.s)
        /// <see</returns>
        public double getViscosityDim(double T) {
            double S = 110.56;
            double T0 = 273.15; // 
            double viscosity0 = 1.716e-5; //kg/( m s) ==> viscosity at T = 273.15 for air
            double viscosity = viscosity0 * Math.Pow(T / T0, 1.5) * (T0 + S) / (T + S);
            return viscosity;
        }


        /// <summary>
        /// Dimensionless Sutherland's law.
        /// </summary>
        /// <param name="phi">Temperature</param>
        /// <returns>
        /// Dynamic viscosity
        /// </returns>
        public override double GetViscosity(double phi) {
            double visc = 0; // nondimensional viscosity
            switch (this.MatParamsMode) {
                case MaterialParamsMode.Constant: {
                        visc = 1.0;
                        break;
                    }
                case MaterialParamsMode.Sutherland: {
                        double S = 110.56;
                        visc = Math.Pow(phi, 1.5) * (1 + S / T_ref) / (phi + S / T_ref);
                        break;
                    }
                case MaterialParamsMode.PowerLaw: {
                        visc = Math.Pow(phi, 2.0 / 3.0);
                        break;
                    }
                default:
                    throw new NotImplementedException();
            }
            Debug.Assert(!double.IsNaN(visc));
            Debug.Assert(!double.IsInfinity(visc));
            Debug.Assert(visc > 0);
            return visc;
        }

        /// <summary>
        ///  The heat conductivity $\lambda. Possibly dependent on the variable <see cref="phi"/> representing the temperature
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public virtual double GetHeatConductivity(double phi) {

            return GetViscosity(phi) / Prandtl;

        }
        /// <summary>
        /// The mass diffusivity,D  multiplied by rho. 
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public virtual double GetDiffusivity(double phi) {
            return GetViscosity(phi) / Prandtl;

        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public double GetPartialHeatCapacity(double phi) {
            return 1.0;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public double GetHeatCapacity(double phi) {
            double cp = 1.0;
            return cp;
        }


        /// <summary>
        /// Returns the value of the heat capacity ratio Gamma (Cp/Cv)
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public double GetHeatCapacityRatio(double phi) {
            double gamma = 1.4;
            return gamma;
        }

        /// <summary>
        /// Returns thermodynamic pressure as function of inital mass and temperature.
        /// </summary>
        /// <param name="InitialMass"></param>
        /// <param name="Temperature"></param>
        /// <returns></returns>
        public override double GetMassDeterminedThermodynamicPressure(double InitialMass, SinglePhaseField Temperature) {
            SinglePhaseField omega = new SinglePhaseField(Temperature.Basis);
            omega.ProjectField(1.0,
                delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    MultidimensionalArray temp = MultidimensionalArray.Create(Len, K);
                    Temperature.Evaluate(j0, Len, NS, temp);
                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {
                            result[j, k] = 1 / temp[j, k];
                        }
                    }
                }, new Foundation.Quadrature.CellQuadratureScheme(true, null));
            return (InitialMass / omega.IntegralOver(null));
        }

        /// <summary>
        /// Returns value of the parameter Lambda/cp
        /// Can be used for calculating the viscosity if the Prandtl number is known (visc = Pr*(lambda/cp))
        /// and also for calculating (rho*Diff_i = 1/Le*(lambda/cp)) if the lewis number is constant.
        /// </summary>
        /// <param name="Temperature"></param>
        /// <returns></returns>
        public double get_LambdaCp_Term(double Temperature) {
            double TREF = 298;
            double res = 2.58e-5 * (Temperature / TREF); // in Kg/(m.s)
            return res;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="VelocityMean"></param>
        /// <param name="Normal"></param>
        /// <param name="ScalarMean"></param>
        /// <returns></returns>
        public override double GetLambda(double[] VelocityMean, double[] Normal, double ScalarMean) {
            throw new NotImplementedException();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public override double DiffRho_Temp(double phi) {
            throw new NotImplementedException();
        }

        public override double getVariableFromZ(double Z, string id) {
            throw new NotImplementedException();
        }

        public override double getDensityFromZ(double Z) {
            throw new NotImplementedException();
        }
    }
}

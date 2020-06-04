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
        public MaterialLawLowMach(double T_ref, MaterialParamsMode MatParamsMode, bool rhoOne)
            : base() {
            this.rhoOne = rhoOne;
            this.T_ref = T_ref;
            this.MatParamsMode = MatParamsMode;
        }

        /// <summary>
        /// 
        /// </summary>
        public override IList<string> ParameterOrdering {
            get {
                return new string[] { VariableNames.Temperature0/*, VariableNames.Rho */};            
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
                //Debug.Assert(phi[0] > -1* 1e-5); // a small treshold. Temperature shouldnt be negative!
                double rho = 1.0;

                if (rhoOne) {
                    rho = 1;
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
        public double GetHeatConductivity(double phi) {
            switch (this.MatParamsMode) {
                case MaterialParamsMode.Constant:
                    return 1.0;
                case MaterialParamsMode.Sutherland: {
                        double S = 110.5;
                        double viscosity = Math.Pow(phi, 1.5) * (1 + S / T_ref) / (phi + S / T_ref);
                        Debug.Assert(!double.IsNaN(viscosity));
                        Debug.Assert(!double.IsInfinity(viscosity));
                        double cp = 1; // It has always 1 as adimensional value, because is a constant, cp/cpref = 1 
                        double lambda = viscosity * cp / 1.0;//  0.71; // Using the fact that Pr = cp*mu/lambda
                        return lambda;
                    }
                case MaterialParamsMode.PowerLaw: {
                        double viscosity = Math.Pow(phi, 2.0 / 3.0);
                        //double S = 110.5;
                        //double viscosity = Math.Pow(phi, 1.5) * (1 + S / T_ref) / (phi + S / T_ref)*0.71;
                        double lambda = viscosity*0.71;
                        Debug.Assert(!double.IsNaN(lambda));
                        Debug.Assert(!double.IsInfinity(lambda));

                        return lambda; // using viscosity = lambda for Pr = cte...
                    }
                default:
                    throw new NotImplementedException();
            }
        }
        /// <summary>
        /// The mass diffusivity,D  multiplied by rho. 
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public double GetDiffusivity(double phi) {

            switch (this.MatParamsMode) {
                case MaterialParamsMode.Constant:
                    return 1.0;
                case MaterialParamsMode.Sutherland: {
                        //GetHeatConductivity(phi);
                        double S = 110.5;
                        double viscosity = Math.Pow(phi, 1.5) * (1 + S / T_ref) / (phi + S / T_ref);
                        Debug.Assert(!double.IsNaN(viscosity));
                        Debug.Assert(!double.IsInfinity(viscosity));
                        double diff = viscosity /1.0; // Sc = mu / rho*D, and for Lewis = 1, Pr = Sc
                        return diff;  
                    }
                case MaterialParamsMode.PowerLaw: {
                        throw new NotImplementedException();
                    }
                default:
                    throw new NotImplementedException();
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="phi"></param>
        /// <returns></returns>
        public double GetPartialHeatCapacity(double phi) {
            switch (this.MatParamsMode) {
                case MaterialParamsMode.Constant:
                    return 1.0;
                case MaterialParamsMode.Sutherland: {
                        //    throw new NotImplementedException();
                        return 1.0; // Using a constant value! 
                    }
                case MaterialParamsMode.PowerLaw: {
                        throw new NotImplementedException();
                    }
                default:
                    throw new NotImplementedException();
            }
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
 
    }
}

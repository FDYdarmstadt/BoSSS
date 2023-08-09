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
using System.Runtime.Serialization;
using System.Linq;
using System.Text;

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution.XheatCommon {


    [DataContract]
    [Serializable]
    public class ThermalParameters : ICloneable {

        /// <summary>
        /// Include nonlinear terms?
        /// Resp.: \f$ Peclet-number \ll 1 \f$ for vanishing convective term
        /// </summary>
        [DataMember]
        public bool IncludeConvection;

        /// <summary>
        /// density of fluid A
        /// </summary>
        [DataMember]
        public double rho_A;

        /// <summary>
        /// density of fluid B
        /// </summary>
        [DataMember]
        public double rho_B;

        /// <summary>
        /// density of solid C
        /// </summary>
        [DataMember]
        public double rho_C;


        /// <summary>
        /// heat capacity of fluid A
        /// </summary>
        [DataMember]
        public double c_A;

        /// <summary>
        /// heat capacity of fluid B
        /// </summary>
        [DataMember]
        public double c_B;

        /// <summary>
        /// heat capacity of solid C
        /// </summary>
        [DataMember]
        public double c_C;


        /// <summary>
        /// thermal conductivity of fluid A
        /// </summary>
        [DataMember]
        public double k_A;

        /// <summary>
        /// thermal conductivity of fluid B
        /// </summary>
        [DataMember]
        public double k_B;

        /// <summary>
        /// thermal conductivity of solid C
        /// </summary>
        [DataMember]
        public double k_C;

        /// <summary>
        /// thermal expansion coefficient of fluid A
        /// used in Boussionesq approximation of buoyancy
        /// </summary>
        [DataMember]
        public double alpha_A = 0.0;

        /// <summary>
        /// thermal expansion coefficient of fluid B
        /// used in Boussionesq approximation of buoyancy
        /// </summary>
        [DataMember]
        public double alpha_B = 0.0;

        /// <summary>
        /// enthalpy of vaporization, a.k.a. (latent) heat of vaporization. It Is always positive (heat is absorbed by the substance).
        /// </summary>
        [DataMember]
        public double hVap = 0.0;

        /// <summary>
        /// prescribed slip length for Temperature slip
        /// </summary>
        [DataMember]
        public double sliplength = 0.0;

        /// <summary>
        /// enthalpy of vaporization, a.k.a. (latent) heat of vaporization. Is always positive (heat is absorbed by the substance).
        /// The enthalpy of condensation (heat of condensation) is by definition equal to h_vap but with opposite sign (heat is released by the substance).
        /// Therefore the enthalpy of vaporization has to be set according to the liquid phase.
        /// </summary>
        //[DataMember]
        //public double hVap_A = 0.0;

        /// <summary>
        /// enthalpy of vaporization, a.k.a. (latent) heat of vaporization. Is always positive (heat is absorbed by the substance).
        /// The enthalpy of condensation (heat of condensation) is by definition equal to h_vap but with opposite sign (heat is released by the substance).
        /// Therefore the enthalpy of vaporization has to be set according to the liquid phase.
        /// </summary>
        //[DataMember]
        //public double hVap_B = 0.0;

        /// <summary>
        /// saturation temperature, is defined as the temperature of the vapor phase adjacent to the interface
        /// </summary>
        [DataMember]
        public double T_sat = -1.0;

        [DataMember]
        public double p_sat = 0.0;

        /// <summary>
        /// condensation coefficient
        /// </summary>
        [DataMember]
        public double fc = 0.0;

        /// <summary>
        /// individual gas constant
        /// </summary>
        [DataMember]
        public double Rc = 0.0;

        /// <summary>
        /// augmented capillary pressure (for testing purpose)
        /// if negativ, the augmented capillary pressure will be calculated
        /// </summary>
        [DataMember]
        public double pc = -1.0;


        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var cl = new ThermalParameters();
            cl.rho_A = this.rho_A;
            cl.rho_B = this.rho_B;
            cl.rho_C = this.rho_C;

            cl.c_A = this.c_A;
            cl.c_B = this.c_B;
            cl.c_C = this.c_C;

            cl.k_A = this.k_A;
            cl.k_B = this.k_B;
            cl.k_C = this.k_C;

            cl.alpha_A = this.alpha_A;
            cl.alpha_B = this.alpha_B;
            cl.hVap = this.hVap;
            cl.T_sat = this.T_sat;
            cl.p_sat = this.p_sat;
            cl.fc = this.fc;
            cl.Rc = this.Rc;
            cl.pc = this.pc;
            cl.sliplength = this.sliplength;

            return MemberwiseClone();
        }

        /// <summary>
        /// 
        /// </summary>
        public override bool Equals(object obj) {

            var cl = obj as ThermalParameters;
            if(cl == null)
                return false;

            return
                cl.rho_A == this.rho_A &&
            cl.rho_B == this.rho_B &&
            cl.rho_C == this.rho_C &&
            cl.c_A == this.c_A &&
            cl.c_B == this.c_B &&
            cl.c_C == this.c_C &&
            cl.k_A == this.k_A &&
            cl.k_B == this.k_B &&
            cl.k_C == this.k_C &&
            cl.alpha_A == this.alpha_A &&
            cl.alpha_B == this.alpha_B &&
            cl.hVap == this.hVap &&
            cl.T_sat == this.T_sat &&
            cl.p_sat == this.p_sat &&
            cl.fc == this.fc &&
            cl.Rc == this.Rc &&
            cl.pc == this.pc;
            cl.sliplength = this.sliplength;

        }

        /// <summary>
        /// Always the same
        /// </summary>
        public override int GetHashCode() {
            return -1;
        }
    }
}

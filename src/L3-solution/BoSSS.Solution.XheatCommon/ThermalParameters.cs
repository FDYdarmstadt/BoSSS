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
        /// enthalpy of vaporization, a.k.a. (latent) heat of vaporization. Is always positive (heat is absorbed by the substance).
        /// The enthalpy of condensation (heat of condensation) is by definition equal to h_vap but with opposite sign (heat is released by the substance).
        /// Therefore the enthalpy of vaporization has to be set according to the liquid phase.
        /// </summary>
        [DataMember]
        public double hVap_A = 0.0;

        /// <summary>
        /// enthalpy of vaporization, a.k.a. (latent) heat of vaporization. Is always positive (heat is absorbed by the substance).
        /// The enthalpy of condensation (heat of condensation) is by definition equal to h_vap but with opposite sign (heat is released by the substance).
        /// Therefore the enthalpy of vaporization has to be set according to the liquid phase.
        /// </summary>
        [DataMember]
        public double hVap_B = 0.0;

        /// <summary>
        /// prescribed volume flux for testing. 
        /// </summary>
        //[DataMember]
        //public double prescribedVolumeFlux = 0.0;

        /// <summary>
        /// is the interface a material one or is it non-material?
        /// </summary>
        //[DataMember]
        //public bool Material = true;


        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var cl = new ThermalParameters();
            cl.rho_A = this.rho_A;
            cl.rho_B = this.rho_B;
            cl.c_A = this.c_A;
            cl.c_B = this.c_B;
            cl.k_A = this.k_A;
            cl.k_B = this.k_B;
            cl.hVap_A = this.hVap_A;
            cl.hVap_B = this.hVap_B;
            //cl.prescribedVolumeFlux = this.prescribedVolumeFlux;
            //cl.Material = this.Material;
            return cl;
        }


    }
}

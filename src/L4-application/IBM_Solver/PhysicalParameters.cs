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
using System.Runtime.Serialization;
using BoSSS.Solution.Control;
using System.Runtime.InteropServices;

namespace BoSSS.Application.IBM_Solver {

    /// <summary>
    /// Fluid properties (e.g. density & viscosity).
    /// </summary>
    [DataContract]
    [Serializable]
    public class PhysicalParameters : ICloneable {
        
        /// <summary>
        /// Include nonlinear terms?
        /// Resp.: Navier-Stokes vs. Stokes
        /// </summary>
        [DataMember]
        public bool IncludeConvection;

        /// <summary>
        /// Density of the fluid.
        /// </summary>
        [DataMember]
        public double rho_A;

        /// <summary> 
        /// Dynamic viscosity the fluid.
        /// </summary>
        [DataMember]
        public double mu_A;

        /// <summary>
        /// is the interface a material one or is it non-material?
        /// </summary>
        [DataMember]
        public bool Material;

        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var cl = new PhysicalParameters();
            cl.IncludeConvection = this.IncludeConvection;
            cl.rho_A = this.rho_A;
            cl.mu_A = this.mu_A;
            cl.Material = this.Material;
            return cl;
        }
        
    }
}


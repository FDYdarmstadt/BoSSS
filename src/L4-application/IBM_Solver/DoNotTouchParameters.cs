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
using BoSSS.Solution.NSECommon;
using System;
using System.Runtime.Serialization;

namespace BoSSS.Application.IBM_Solver {


    /// <summary>
    /// Advanced settings for the 
    /// </summary>
    [DataContract]
    [Serializable]
    public class DoNotTouchParameters : ICloneable {
        
     
        /// <summary>
        /// Continuity equation: work with div(-) resp. -div(-)
        /// </summary>
        [DataMember]
        public double ContiSign = -1.0;

        /// <summary>
        /// scale continuity equation with one over density
        /// </summary>
        [DataMember]
        public bool RescaleConti = false;

        /// <summary>
        /// stabilization parameter for Local-Lax-Friedrichs flux, phase A
        /// </summary>
        [DataMember]
        public double LFFA = 0.8;

        /// <summary>
        /// stabilization parameter for Local-Lax-Friedrichs flux, phase B
        /// </summary>
        [DataMember]
        public double LFFB = 0.8;

        /// <summary>
        /// Penalty safety factor for the viscous operator.
        /// </summary>
        [DataMember]
        public double PenaltySafety = 4;
        
        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public double CellAgglomerationThreshold = 0.1;
        

        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var cl = new DoNotTouchParameters();
            cl.CellAgglomerationThreshold = this.CellAgglomerationThreshold;
            cl.ContiSign = this.ContiSign;
            cl.RescaleConti = this.RescaleConti;
            cl.LFFA = this.LFFA;
            cl.LFFB = this.LFFB;
            //cl.ViscosityMode = this.ViscosityMode;
            //cl.ViscosityImplementation = this.ViscosityImplementation;
            return cl;
        }
    }
}

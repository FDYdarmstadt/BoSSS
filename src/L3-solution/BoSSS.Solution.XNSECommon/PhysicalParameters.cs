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

namespace BoSSS.Solution.XNSECommon {

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
        /// density of fluid A
        /// </summary>
        [DataMember]
        public double rho_A;

        /// <summary>
        /// density of fluid B
        /// </summary>
        [DataMember]
        public double rho_B;

        /// <summary> dynamic viscosity fluid A </summary>
        [DataMember]
        public double mu_A;

        /// <summary> dynamic viscosity fluid B </summary>
        [DataMember]
        public double mu_B;

        /// <summary>
        /// surface tension
        /// </summary>
        [DataMember]
        public double Sigma;

        /// <summary>
        /// surface shear viscosity
        /// </summary>
        [DataMember]
        public double mu_I = 0.0;

        /// <summary>
        /// surface dilatational viscosity
        /// </summary>
        [DataMember]
        public double lambda_I = 0.0;

        /// <summary>
        /// dissipation coefficient for the effective wall force (fluid A)
        /// </summary>
        [DataMember]
        public double betaS_A = 0.0;

        /// <summary>
        /// dissipation coefficient for the effective wall force (fluid B)
        /// </summary>
        [DataMember]
        public double betaS_B = 0.0;

        /// <summary>
        /// dissipation coefficient for the effective contact line force
        /// </summary>
        [DataMember]
        public double betaL = 0.0;

        /// <summary>
        /// static contact angle
        /// </summary>
        [DataMember]
        public double theta_e = Math.PI / 2.0;

        /// <summary>
        /// prescribed slip length for GNBC
        /// </summary>
        [DataMember]
        public double sliplength = 0.0;

        /// <summary>
        /// is the interface a material one or is it non-material?
        /// </summary>
        [DataMember]
        public bool Material;

        /// <summary>
        /// prescribed volume flux for testing. 
        /// </summary>
        //[DataMember]
        //public double prescribedVolumeFlux = 0.0;

        /// <summary>
        /// Use the artificial surface force (usually only used in manufactured solutions)?
        /// </summary>
        [DataMember]
        public bool useArtificialSurfaceForce = false;

        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var cl = new PhysicalParameters();
            cl.IncludeConvection = this.IncludeConvection;
            cl.rho_A = this.rho_A;
            cl.rho_B = this.rho_B;
            cl.mu_A = this.mu_A;
            cl.mu_B = this.mu_B;
            cl.Sigma = this.Sigma;
            cl.mu_I = this.mu_I;
            cl.lambda_I = this.lambda_I;
            cl.betaS_A = this.betaS_A;
            cl.betaS_B = this.betaS_B;
            cl.betaL = this.betaL;
            cl.theta_e = this.theta_e;
            cl.sliplength = this.sliplength;
            cl.Material = this.Material;
            //cl.prescribedVolumeFlux = this.prescribedVolumeFlux;
            cl.useArtificialSurfaceForce = this.useArtificialSurfaceForce;
            return cl;
        }

        /// <summary>
        /// Constants for a water (species A)/air (species B) pair.
        /// </summary>
        public static PhysicalParameters WaterAir {
            get {
                PhysicalParameters C = new PhysicalParameters();
                C.rho_A = 1000;
                C.rho_B = 1.2;
                C.mu_A = 1.0e-3;
                C.mu_B = 17.1e-6;
                C.Sigma = 72.75e-3;
                return C;
            }
        }

        /// <summary>
        /// Constants for an air (species A)/water (species B) pair.
        /// </summary>
        public static PhysicalParameters AirWater {
            get {
                PhysicalParameters C = new PhysicalParameters();
                C.rho_B = 1000;
                C.rho_A = 1.2;
                C.mu_B = 1.0e-3;
                C.mu_A = 17.1e-6;
                C.Sigma = 72.75e-3;
                return C;
            }
        }


        /// <summary>
        /// All constants set to 1.0.
        /// </summary>
        public static PhysicalParameters AllOne {
            get {
                PhysicalParameters C = new PhysicalParameters();
                C.rho_B = 1.0;
                C.rho_A = 1.0;
                C.mu_B = 1.0;
                C.mu_A = 1.0;
                C.Sigma = 1.0;
                return C;
            }
        }
    }
}


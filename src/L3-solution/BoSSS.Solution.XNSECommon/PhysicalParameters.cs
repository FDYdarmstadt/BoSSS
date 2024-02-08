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
        /// Include Diffusive Term, standard: this is on
        /// </summary>
        [DataMember]
        public bool IncludeDiffusion = true;

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

        /// <summary> Reynolds fluid B for dimensionless </summary>
        [DataMember]
        public double reynolds_B;

        /// <summary> Reynolds fluid A for dimensionless </summary>
        [DataMember]
        public double reynolds_A;

        /// <summary>
        /// surface tension
        /// </summary>
        [DataMember]
        public double Sigma;

        /// <summary>
        /// pressure at the interface for free surface flows 
        /// </summary>
        [DataMember]
        public double pFree = 0.0;

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
        /// coefficient for the surface divergence term in the Boussinesq-Scriven surface model
        /// if negative the value will be computed by (surface dilatational viscosity(lambda_I) - surface shear viscosity(mu_I))
        /// </summary>
        [DataMember]
        public double lambdaI_tilde = -1.0;

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
        /// effective sliplength for the fluid-fluid interface
        /// </summary>
        [DataMember]
        public double slipI = 0.0;

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
        /// Use the artificial surface force (usually only used in manufactured solutions)?
        /// </summary>
        [DataMember]
        public bool useArtificialSurfaceForce = false;

        /// <summary>
        /// clone
        /// </summary>
        public virtual object Clone() {
            var cl = (PhysicalParameters)MemberwiseClone(); // ok for this object, since it contains only value types
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

        public override bool Equals(object obj) {
            var other = obj as PhysicalParameters;
            if (other == null)
                return false;

            return
                this.IncludeConvection == other.IncludeConvection &&
                this.rho_A == other.rho_A &&
                this.rho_B == other.rho_B &&
                this.mu_A == other.mu_A &&
                this.mu_B == other.mu_B &&
                this.mu_I == other.mu_I &&
                this.Sigma == other.Sigma &&
                this.pFree == other.pFree &&
                this.betaS_A == other.betaS_A &&
                this.betaS_B == other.betaS_B &&
                this.betaL == other.betaL &&
                this.theta_e == other.theta_e &&
                this.sliplength == other.sliplength &&
                this.Material == other.Material &&
                this.reynolds_A == other.reynolds_A &&
                this.reynolds_B == other.reynolds_B &&
                this.lambda_I == other.lambda_I &&
                this.lambdaI_tilde == other.lambdaI_tilde &&
                this.useArtificialSurfaceForce == other.useArtificialSurfaceForce;
        }
    }




    [DataContract]
    [Serializable]
    public class PhysicalParametersRheology : PhysicalParameters {

        /// <summary>
        /// viscoelastic dimensionless: ratio between relaxation time and retardation time (fluid A)
        /// </summary>
        [DataMember]
        public double beta_a = 1.0;

        /// <summary>
        /// viscoelastic dimensionless: ratio between relaxation time and retardation time (fluid B)
        /// </summary>
        [DataMember]
        public double beta_b = 1.0;

        /// <summary>
        /// viscoelastic dimensionless: Weissenberg number (fluid A)
        /// </summary>
        [DataMember]
        public double Weissenberg_a = 0.0;

        /// <summary>
        /// viscoelastic dimensionless: Weissenberg number (fluid B)
        /// </summary>
        [DataMember]
        public double Weissenberg_b = 0.0;

        /// <summary>
        /// Giesekus dimensionless: Giesekus factor (fluid A)
        /// </summary>
        [DataMember]
        public double giesekusfactor_a = 0.0;

        /// <summary>
        /// Giesekus dimensionless: Giesekus factor (fluid B)
        /// </summary>
        [DataMember]
        public double giesekusfactor_b = 0.0;

        /// <summary>
        /// Giesekus in fluid A?
        /// </summary>
        [DataMember]
        public bool Giesekus_a = false;

        /// <summary>
        /// Giesekus in fluid B?
        /// </summary>
        [DataMember]
        public bool Giesekus_b = false;


        /// <summary>
        /// clone
        /// </summary>
        public override object Clone() {
            var cl = (PhysicalParameters)MemberwiseClone(); // ok for this object, since it contains only value types
            return cl;
        }
    }
    [DataContract]
    [Serializable]
    public class PhysicalParametersCombustion : PhysicalParameters {

        /// <summary>
        /// Diffusivity factor in fluid A
        /// </summary>
        [DataMember]
        public double rhoD_A = 1.0;

        /// <summary>
        /// Diffusivity factor in fluid B
        /// </summary>
        [DataMember]
        public double rhoD_B = 1.0;

    }



}


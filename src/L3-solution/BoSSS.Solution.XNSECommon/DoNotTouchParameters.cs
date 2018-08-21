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

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// Options of the surface stress tensor \f$ \sigma_{\Gamma} \f$ for the momentum jump condition 
    /// \f[ 
    ///   \llbracket (-p \myMatrix{I} + \mu \myMatrix{D} (\vec{u}) ) \cdot \vec{n} \rrbracket = \divergence{\sigma_{\Gamma}}_{\Gamma},
    /// \f]
    /// </summary>
    public enum SurfaceSressTensor {

        /// <summary>
        /// only the isotropic \f$ \sigma_{\gamma} = \sigma \myMatrix{P}_{\Gamma} \f$
        /// </summary>
        Isotropic,

        /// <summary>
        /// \f$ \mu_I \( \myMatrix{P}_I \grad_I \vec{u} + \( \grad_I \vec{u}\)^T \myMatrix{P}_I \) \f$ (additional dynamic part, resembling the semi-implicit discretization)
        /// </summary>
        SurfaceRateOfDeformation,

        /// <summary>
        /// \f$ \( \lambda_I - \mu_I \) \div_I \vec{u} \myMatrix{P}_I  \f$
        /// </summary>
        SurfaceVelocityDivergence,

        /// <summary>
        /// the Boussinesq-Scriven model, describing a newtonian-like surface with surface shear and dilatational viscosity
        /// </summary>
        FullBoussinesqScriven,

        /// <summary>
        /// resembles <see cref="SurfaceRateOfDeformation"/> with $\mu_I = dt$
        /// </summary>
        SemiImplicit

    }


    /// <summary>
    /// Options for the treatment of the isotropic part of the surface stress tensor 
    /// </summary>
    public enum SurfaceStressTensor_IsotropicMode {
        
        /// <summary>
        /// Curvature is evaluated locally, i.e. the projection of $\divergence{ \nabla \varphi / | \nabla \varphi | $.
        /// onto a DG field is used.
        /// </summary>
        Curvature_Projected,

        /// <summary>
        /// Curvature is evaluated locally at closest points on the level-set, i.e. 
        /// the value of $\divergence{ \nabla \varphi / | \nabla \varphi | $ on the zero-set is extended to the domain.
        /// </summary>
        Curvature_ClosestPoint,


        /// <summary>
        /// A cell-wise mean-curvature is computed from a Laplace-Beltrami-ansatz.
        /// </summary>
        Curvature_LaplaceBeltramiMean,

        /// <summary>
        /// Curvature is computed from a specialized Level-Set as Fourier-series
        /// </summary>
        Curvature_Fourier,
        
        /// <summary>
        /// use a cell-wise Laplace-Beltrami formulation.
        /// </summary>
        LaplaceBeltrami_Local,

        /// <summary>
        /// a Laplace-Beltrami formulation where level-set tangents are averaged at the cell boundary.
        /// </summary>
        LaplaceBeltrami_Flux,

        /// <summary>
        /// a cell-wise Laplace-Beltrami formulation with handling of the contact line
        /// </summary>
        LaplaceBeltrami_ContactLine


    }



    /// <summary>
    /// 
    /// </summary>
    public enum ViscosityMode {
    
        /// <summary>
        /// recommended
        /// </summary>
        Standard,

        /// <summary>
        /// in the special case of \f$ \mu_{\mathfrak{A}} = \mu_{\mathfrak{B}}\f$ ,
        /// this yields a symmetric discretization
        /// </summary>
        /// <remarks>
        /// seem to produce crappy results
        /// </remarks>
        ExplicitTransformation,

        /// <summary>
        /// the full viscous stress tensor is discretized in thw bulk domain, i.e. 
        /// \f[ 
        ///    \mathrm{div} \left( \mu \left( \nabla \vec{u} + \vec{u}^T \right) \right),
        /// \f]
        /// therefore the jump condition can be implemented symmetric.
        /// </summary>
        /// <remarks>
        /// difficult Neumann boundary condition
        /// </remarks>
        FullySymmetric,

        /// <summary>
        /// for test purposes, an implementation that neglects the term
        /// \f[ 
        ///    \mathrm{div} \left( \mu \vec{u}^T \right),
        /// \f]
        /// </summary>
        TransposeTermMissing

    }


    /// <summary>
    /// Advanced settings for the 
    /// </summary>
    [DataContract]
    [Serializable]
    public class DoNotTouchParameters : ICloneable {

        /// <summary>
        /// viscosity discretization
        /// </summary>
        [DataMember]
        public ViscosityMode ViscosityMode = ViscosityMode.FullySymmetric;

        /*
        /// <summary>
        /// viscosity Implementation
        /// H: SIP
        /// SWIP: weighted
        /// </summary>
        [DataMember]
        public ViscosityImplementation ViscosityImplementation = ViscosityImplementation.H;
        */
        
        /// <summary>
        /// Turn the use of ghost penalties on or off, see <br/>
        /// @article{massjung_unfitted_2012,
        ///         title = {An {Unfitted} {Discontinuous} {Galerkin} {Method} {Applied} to {Elliptic} {Interface} {Problems}},
        ///         volume = {50},
        ///         issn = {0036-1429, 1095-7170},
        ///         url = {http://epubs.siam.org/doi/abs/10.1137/090763093},
        ///         doi = {10.1137/090763093},
        ///         language = {en},
        ///         number = {6},
        ///         urldate = {2014-11-03},
        ///         journal = {SIAM Journal on Numerical Analysis},
        ///         author = {Massjung, Ralf},
        ///         month = jan,
        ///         year = {2012},
        ///         pages = {3134--3162}
        /// }
        /// </summary>
        [DataMember]
        public bool UseGhostPenalties = false;

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
        /// Model for the surface stress tensor
        /// </summary>
        [DataMember]
        public SurfaceSressTensor SurfStressTensor = SurfaceSressTensor.Isotropic;

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public bool UseLevelSetStabilization = false;
            
        /// <summary>
        /// implementation variant of the isotropic surface stress
        /// </summary>
        [DataMember]
        public SurfaceStressTensor_IsotropicMode SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux;

        /// <summary>
        /// Expert options regarding the evaluation of the curvature.
        /// </summary>
        [DataMember]
        public CurvatureAlgorithms.FilterConfiguration FilterConfiguration = new CurvatureAlgorithms.FilterConfiguration();

        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var cl = new DoNotTouchParameters() {
                CellAgglomerationThreshold = this.CellAgglomerationThreshold,
                ContiSign = this.ContiSign,
                RescaleConti = this.RescaleConti,
                LFFA = this.LFFA,
                LFFB = this.LFFB,
                PenaltySafety = this.PenaltySafety,
                SurfStressTensor = this.SurfStressTensor,
                SST_isotropicMode = this.SST_isotropicMode,
                UseLevelSetStabilization = this.UseLevelSetStabilization,
                ViscosityMode = this.ViscosityMode,
                //ViscosityImplementation = this.ViscosityImplementation,
                UseGhostPenalties = this.UseGhostPenalties,
                FilterConfiguration = this.FilterConfiguration
            };
            return cl;
        }
    }
}

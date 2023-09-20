﻿/* =======================================================================
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
        SurfaceDivergence,

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
    /// Options for the localization of the slip-region for the generalized navier-slip boundary condition
    /// </summary>
    public enum NavierSlip_Localization {

        /// <summary>
        /// GNBC on the whole boundary
        /// </summary>
        Bulk = 0,

        /// <summary>
        /// GNBC localized to the nearband boundary, elsewhere velocity inlet
        /// </summary>
        Nearband = 1,

        /// <summary>
        /// only the contact line part of GNBC
        /// </summary>
        ContactLine = 2,

        /// <summary>
        /// prescribed prescription of the slip length
        /// </summary>
        Prescribed = 3,

        /// <summary>
        /// write sliplengh for each cell
        /// </summary>
        Everywhere = 4

    }

    /// <summary>
    /// Options for the localization of the slip-region for the generalized navier-slip boundary condition
    /// </summary>
    public enum ThermalSlip_Localization {

        /// <summary>
        /// GNBC on the whole boundary
        /// </summary>
        Bulk = 0,

        /// <summary>
        /// GNBC localized to the nearband boundary, elsewhere velocity inlet
        /// </summary>
        Nearband = 1,

        /// <summary>
        /// only the contact line part of GNBC
        /// </summary>
        ContactLine = 2,

        /// <summary>
        /// prescribed prescription of the slip length
        /// </summary>
        Prescribed = 3,

        /// <summary>
        /// write sliplengh for each cell
        /// </summary>
        Everywhere = 4

    }

    /// <summary>
    /// Options for the boundary type of the immersed boundary
    /// </summary>
    public enum IBM_BoundaryType {

        /// <summary>
        /// No Slip boundary
        /// </summary>
        NoSlip = 0,

        /// <summary>
        /// Navierslip boundary, (with cl handling)
        /// </summary>
        NavierSlip = 1,

        /// <summary>
        /// Freeslip boundary, (with cl handling)
        /// </summary>
        FreeSlip = 2
    }

    /// <summary>
    /// Options for the boundary type of the immersed boundary
    /// </summary>
    public enum IBM_ThermalBoundaryType {

        /// <summary>
        /// No Slip boundary
        /// </summary>
        NoSlip = 0,

        /// <summary>
        /// Slip
        /// </summary>
        ThermalSlip = 1
    }

    /// <summary>
    /// Options for the friction coefficient in GNBC
    /// </summary>
    public enum NavierSlip_SlipLength {

        /// <summary>
        /// prescribed frictions coeff from control file
        /// </summary>
        Prescribed_Beta = 0,

        /// <summary>
        /// prescribed slipLength from control file
        /// </summary>
        Prescribed_SlipLength = 1,

        /// <summary>
        /// friction coeff computed according to h_min
        /// </summary>
        hmin_Grid = 2,

        /// <summary>
        /// friction coeff computed according to h_min / (p + 1)
        /// </summary>
        hmin_DG = 3,

    }


    /// <summary>
    /// 
    /// </summary>
    public enum ViscosityMode {
    
        /// <summary>
        /// the same as <see cref="TransposeTermMissing"/>
        /// </summary>
        Standard,

        
        /// <summary>
        /// the full viscous stress tensor is discretized in the bulk domain, i.e. 
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
        TransposeTermMissing,

        /// <summary>
        /// In viscoelastic case we calculate dimensionless and have the material parameter \beta
        /// such that the dimensionless "Viscosity" is defined as 
        /// \f[ 
        ///     \frac{\beta}{\mathrm{Re}}
        /// \f]
        /// </summary>
        Viscoelastic

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
        */

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
        
        /*
        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public double CellAgglomerationThreshold = 0.1;
        */

        /// <summary>
        /// Model for the surface stress tensor
        /// </summary>
        [DataMember]
        public SurfaceSressTensor SurfStressTensor = SurfaceSressTensor.Isotropic;

        public enum SurfaceTensionForceStabilization {

            None,

            surfaceDeformationRateLocal,

            GradUxGradV,

            surfaceDivergence,

            EdgeDissipation,

        }

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public SurfaceTensionForceStabilization STFstabilization = SurfaceTensionForceStabilization.None;

        /// <summary>
        /// set the the maximum surface tension coefficient value according to given time step restriction
        /// </summary>
        [DataMember]
        public bool SetSurfaceTensionMaxValue = false;

        /// <summary>
        /// switch for free surface flows, thus slip is allowed at the interface
        /// </summary>
        [DataMember]
        public bool freeSurfaceFlow = false;


        /// <summary>
        /// implementation variant of the isotropic surface stress
        /// </summary>
        [DataMember]
        public SurfaceStressTensor_IsotropicMode SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux;

        /// <summary>
        /// in case of a Laplace-Beltrami surface tension computation the curvature is need as parameter
        /// </summary>
        [DataMember]
        public bool CurvatureNeeded = false;

        /// <summary>
        /// Expert options regarding the evaluation of the curvature.
        /// </summary>
        [DataMember]
        public CurvatureAlgorithms.FilterConfiguration FilterConfiguration = CurvatureAlgorithms.FilterConfiguration.NoFilter;


        /// <summary>
        /// See <see cref="NavierSlip_Localization"/>
        /// </summary>
        [DataMember]
        public NavierSlip_Localization GNBC_Localization = NavierSlip_Localization.Bulk;

        /// <summary>
        /// See <see cref="NavierSlip_SlipLength"/>
        /// </summary>
        [DataMember]
        public NavierSlip_SlipLength GNBC_SlipLength = NavierSlip_SlipLength.Prescribed_Beta;

        /// <summary>
        /// See <see cref="NavierSlip_Localization"/>
        /// </summary>
        [DataMember]
        public ThermalSlip_Localization ThermalSlip_Localization = ThermalSlip_Localization.Everywhere;

        /// <summary>
        /// See <see cref="IBM_BoundaryType"/>
        /// </summary>
        [DataMember]
        public IBM_BoundaryType IBM_BoundaryType = IBM_BoundaryType.NoSlip;

        /// <summary>
        /// See <see cref="IBM_ThermalBoundaryType"/>
        /// </summary>
        [DataMember]
        public IBM_ThermalBoundaryType IBM_ThermalBoundaryType = IBM_ThermalBoundaryType.NoSlip;

        //viscoelastic LDG stuff:
        //=========================

        /// <summary>
        /// determines which implementation of objective Term should be used:
        /// =1 only velocity gradient as param
        /// =0 only stress tensor as param
        /// </summary>
        [DataMember]
        public double ObjectiveParam = 1.0;

        /// <summary>
        /// Upwinding factor for convective part in constitutive equation
        /// </summary>
        [DataMember]
        public double alpha = 1.0;

        /// <summary>
        /// Penalty Values LDG (alpha, beta; Lit. Cockburn (2002) Local DG Methods for the Stokes system)
        /// Penalty in Stress Divergence (beta)
        /// </summary>
        [DataMember]
        public double[] Penalty1 = { 0, 0 };

        /// <summary>
        /// Penalty in Constitutive Viscosity (alpha)
        /// </summary>
        [DataMember]
        public double Penalty2 = 1.0;

        /// <summary>
        /// Penalty for pressure/conti (beta)
        /// </summary>
        [DataMember]
        public double[] PresPenalty1 = { 0, 0 };

        /// <summary>
        /// Penalty for pressure (alpha)
        /// </summary>
        [DataMember]
        public double PresPenalty2 = 1.0;

        /// <summary>
        /// penalty for stress in objective term
        /// </summary>
        [DataMember]
        public double StressPenalty = 1.0;

        /// <summary>
        /// double cut cell special handling override <see cref="BoSSS.Foundation.XDG.Quadrature.BruteForceSettingsOverride"/>
        /// </summary>
        [DataMember]
        public bool DoubleCutSpecialQuadrature = false;


        /// <summary>
        /// clone
        /// </summary>
        public object Clone() {
            var cl = (DoNotTouchParameters)MemberwiseClone(); // ok for value type memebers
            return cl;
        }
    }
}

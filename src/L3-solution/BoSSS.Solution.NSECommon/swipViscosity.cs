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
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using BoSSS.Platform;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Foundation;
using ilPSP;

namespace BoSSS.Solution.NSECommon {

    /*
    /// <summary>
    /// A configuration switch to <see cref="swipViscosityBase"/> and decendants, to
    /// switch between different variants of inhomogeneous-diffusion forms
    /// </summary>
    public enum ViscosityImplementation {

        /// <summary>
        /// SWIP-form according to 
        /// <code>
        /// @book{di_pietro_mathematical_2011,
        ///       series = {Math Matiques Et Applications},
        ///       title = {Mathematical Aspects of Discontinuous Galerkin Methods},
        ///       isbn = {9783642229794},
        ///       number = {69},
        ///       publisher = {Springer},
        ///       author = {Di Pietro, Daniele Antonio and Ern, Alexandre},
        ///       year = {2011},
        /// }
        /// </code>
        /// (page 155)
        /// </summary>
        SWIP,

        /// <summary>
        /// an alternate consistent implementation
        /// </summary>
        H
    }
    */

    /// <summary>
    /// A configuration switch to <see cref="swipViscosityBase"/> and decendants, to
    /// switch constant / variable viscosity
    /// and non-dimensionless / dimensionless formulation.
    /// </summary>
    public enum ViscosityOption {

        /// <summary>
        /// Constant viscosity for non-dimesionless formulation.
        /// </summary>
        ConstantViscosity,

        /// <summary>
        /// Constant viscosity for dimesionless formulation.
        /// Reynolds number has to be given.
        /// </summary>
        ConstantViscosityDimensionless,

        /// <summary>
        /// Variable viscosity for non-dimesionless formulation.
        /// </summary>
        VariableViscosity,

        /// <summary>
        /// Variable viscosity for dimesionless formulation.
        /// Reynolds number has to be given.
        /// </summary>
        VariableViscosityDimensionless
    }


    /// <summary>
    /// base class for viscosity terms
    /// </summary>
    public abstract class swipViscosityBase : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="_penalty"></param>
        /// <param name="iComp">
        /// component index
        /// </param>
        /// <param name="D">
        /// spatial dimension.
        /// </param>        
        /// <param name="bcmap"></param>
        /// <param name="_ViscosityMode">
        /// see <see cref="ViscosityOption"/>
        /// </param>
        /// <param name="constantViscosityValue">
        /// Constant value for viscosity. 
        /// Needs to be given for <see cref="ViscosityOption.ConstantViscosity"/>.
        /// </param>
        /// <param name="reynolds">
        /// Reynolds number for dimensionless formulation.
        /// Needs to be given for <see cref="ViscosityOption.ConstantViscosityDimensionless"/> and <see cref="ViscosityOption.VariableViscosityDimensionless"/>.
        /// </param>
        /// <param name="EoS">
        /// Optional material law for calculating the viscosity
        /// as a function of the level-set.
        /// Only available for <see cref="ViscosityOption.VariableViscosity"/> and <see cref="ViscosityOption.VariableViscosityDimensionless"/>.
        /// </param>
        protected swipViscosityBase(
            double _penalty, MultidimensionalArray PenaltyLengthScales, 
            int iComp, int D, IncompressibleBoundaryCondMap bcmap, 
            ViscosityOption _ViscosityMode, double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null, 
            Func<double,int, int, MultidimensionalArray, double> ComputePenalty = null) {
            this.m_penalty = _penalty;
            this.m_ComputePenalty = ComputePenalty;
            this.m_iComp = iComp;
            this.m_D = D;
            this.cj = PenaltyLengthScales;
            velFunction = D.ForLoop(d => bcmap.bndFunction[VariableNames.Velocity_d(d)]);
            EdgeTag2Type = bcmap.EdgeTag2Type;
            this.m_PhysicsMode = bcmap.PhysMode;

            this.m_ViscosityMode = _ViscosityMode;
            switch (_ViscosityMode) {
                case ViscosityOption.ConstantViscosity:
                    if (double.IsNaN(constantViscosityValue))
                        throw new ArgumentException("constantViscosityValue is missing!");
                    this.m_constantViscosityValue = constantViscosityValue;
                    break;
                case ViscosityOption.ConstantViscosityDimensionless:
                    if (double.IsNaN(reynolds))
                        throw new ArgumentException("reynolds number is missing!");
                    this.m_reynolds = reynolds;
                    break;
                case ViscosityOption.VariableViscosity:
                    this.m_EoS = EoS;
                    break;
                case ViscosityOption.VariableViscosityDimensionless:
                    if (double.IsNaN(reynolds))
                        throw new ArgumentException("reynolds number is missing!");
                    this.m_reynolds = reynolds;
                    this.m_EoS = EoS;
                    break;
                default:
                    throw new NotImplementedException();
            }

        }

        /// <summary>
        /// a multiplier which is applied to everything but the penalty terms
        /// </summary>
        protected double m_alpha = 1.0;



        /// <summary>
        /// see <see cref="BoundaryCondMap{BCType}.EdgeTag2Type"/>;
        /// </summary>
        protected IncompressibleBcType[] EdgeTag2Type;

        /// <summary>
        /// component index
        /// </summary>
        protected int m_iComp;

        /// <summary>
        /// Cell-wise length scales for the penalty computation.
        /// </summary>
        MultidimensionalArray cj;

        /// <summary>
        /// multiplier for the penalty computation
        /// </summary>
        double m_penalty;

        /// <summary>
        /// spatial dimension
        /// </summary>
        protected int m_D;

        /// <summary>
        /// Dirichlet boundary values; <br/>
        ///  - 1st index: spatial dimension <br/>
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[][] velFunction;

        /// <summary>
        /// dissipation coefficient for the effective wall force
        /// </summary>
        protected double m_beta = 0.0;


        Func<double,int, int, MultidimensionalArray, double> m_ComputePenalty;

        /// <summary>
        /// see <see cref="ViscosityOption"/>
        /// </summary>
        ViscosityOption m_ViscosityMode;

        /// <summary>
        /// in the case of constant viscosity, the value of the viscosity
        /// </summary>
        double m_constantViscosityValue = double.NaN;

        /// <summary>
        /// Reynolds number for dimensionless formulation.
        /// </summary>
        double m_reynolds = double.NaN;

        /// <summary>
        /// Optional material law for calculating the viscosity
        /// as a function of temperature / level-set.
        /// </summary>
        MaterialLaw m_EoS = null;

        /// <summary>
        /// PhysicsMode needed for ParamterOrdering.
        /// </summary>
        PhysicsMode m_PhysicsMode;

        /// <summary>
        /// the molecular viscosity
        /// </summary>
        virtual protected double Viscosity(double[] Parameters) {
            switch (m_ViscosityMode) {
                case ViscosityOption.ConstantViscosity:
                    return m_constantViscosityValue;
                case ViscosityOption.ConstantViscosityDimensionless:
                    return (1.0 / m_reynolds);
                case ViscosityOption.VariableViscosity:
                    if (m_EoS == null) {
                        return Parameters[0];
                    } else {
                        return m_EoS.GetViscosity(Parameters[0]);
                    }
                case ViscosityOption.VariableViscosityDimensionless:
                    if (m_EoS == null) {
                        return (Parameters[0] / m_reynolds);
                    } else {
                        return (m_EoS.GetViscosity(Parameters[0]) / m_reynolds);
                    }
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        virtual protected double penalty(int jCellIn, int jCellOut) {
            double mu;
            if(m_ComputePenalty != null) {
                mu = m_ComputePenalty(m_penalty, jCellIn, jCellOut, cj);
            } else {
                double cj_in = cj[jCellIn];
                mu = m_penalty * cj_in;
                if(jCellOut >= 0) {
                    double cj_out = cj[jCellOut];
                    mu = Math.Max(mu, m_penalty * cj_out);
                }
            }
            return mu;
        }

        /// <summary>
        /// returns the velocity vector;
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(m_D);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                switch (m_ViscosityMode) {
                    case ViscosityOption.ConstantViscosity:
                    case ViscosityOption.ConstantViscosityDimensionless:
                        return new string[0];
                    case ViscosityOption.VariableViscosity:
                    case ViscosityOption.VariableViscosityDimensionless:
                        if (m_EoS == null) {
                            return new string[] { VariableNames.ViscosityMolecular };
                        } else {
                            switch (m_PhysicsMode) {
                                case PhysicsMode.LowMach:
                                case PhysicsMode.Combustion:
                                    return new string[] { VariableNames.Temperature0 };
                                case PhysicsMode.Multiphase:
                                    return new string[] { VariableNames.LevelSet };
                                case PhysicsMode.Viscoelastic:
                                case PhysicsMode.Incompressible:
                                    throw new ApplicationException("Should not happen.");
                                default:
                                    throw new NotImplementedException();
                            }
                        }
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// very dirty hack to 'inject' an alternate boundary condition value for unit testing,
        /// designed to match <see cref="BoSSS.Application.ipViscosity.TestSolution.U"/>
        /// </summary>
        public Func<int, double[], double> g_Diri_Override;

        /// <summary>
        /// very dirty hack to 'inject' an alternate boundary condition value for unit testing,
        /// designed to match <see cref="BoSSS.Application.ipViscosityTestSolution.dU"/>
        /// </summary>
        public Func<int, double[], int, double> g_Neu_Override;


        /// <summary>
        /// Dirichlet boundary value: the given velocity at the boundary.
        /// </summary>
        protected double g_Diri(double[] X, double time, int EdgeTag, int d) {
            if (this.g_Diri_Override == null) {
                Func<double[], double, double> boundVel = this.velFunction[d][EdgeTag];
                double ret = boundVel(X, time);

                //if (m_UseBoundaryVelocityParameter)
                //    ret += inp.ParameterValuesIn[0];

                return ret;
            } else {
                return g_Diri_Override(d, X);
            }
        }

        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        abstract public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV);

        abstract public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB);

        abstract public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA);
    }

    /// <summary>
    /// \f[ 
    ///   -\operatorname{div} \left( \mu \nabla \vec{u} \right)
    /// \f]
    /// </summary>
    public class swipViscosity_Term1 : swipViscosityBase {

        /// <summary>
        /// ctor; parameter documentation see <see cref="swipViscosityBase.swipViscosityBase"/>.
        /// </summary>
        public swipViscosity_Term1(double _penalty, MultidimensionalArray PenaltyLengthScales, int iComp, int D, IncompressibleBoundaryCondMap bcmap, 
            ViscosityOption _ViscosityMode, double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null, 
            Func<double,int,int,MultidimensionalArray,double> ComputePenalty = null)
            : base(_penalty, PenaltyLengthScales, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS, ComputePenalty) {
            
        }

        //public override IList<string> ArgumentOrdering
        //{
        //    get
        //    {
        //        return new string[] { VariableNames.Velocity_d(m_iComp) };
        //    }
        //}

        public override double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for (int d = 0; d < cpv.D; d++)
                //acc -= GradU[0, d] * GradV[d] * Viscosity(cpv.Parameters) * base.m_alpha;
                acc -= GradU[m_iComp, d] * GradV[d] * Viscosity(cpv.Parameters) * base.m_alpha;
            return -acc;
        }



        public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            double muB = this.Viscosity(inp.Parameters_OUT);


            //switch (base.m_implMode) {
            //    case ViscosityImplementation.H: //
            {
                for (int d = 0; d < inp.D; d++) {
                    //Acc += 0.5 * (muA * _Grad_uA[0, d] + muB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                    //Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
                    Acc += 0.5 * (muA * _Grad_uA[m_iComp, d] + muB * _Grad_uB[m_iComp, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                    Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[m_iComp] - _uB[m_iComp]) * inp.Normale[d];  // symmetry term
                }
                Acc *= base.m_alpha;

                double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
                //Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * muMax; // penalty term
                Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax; // penalty term

                return -Acc;

            }

            //    case ViscosityImplementation.SWIP: //
            //    {
            //            for (int d = 0; d < inp.D; d++) {
            //                //Acc += (muB * muA * _Grad_uA[0, d] + muA * muB * _Grad_uB[0, d]) / (muA + muB) * (_vA - _vB) * inp.Normale[d];  // consistency term
            //                //Acc += (muB * muA * _Grad_vA[d] + muA * muB * _Grad_vB[d]) / (muA + muB) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            //                Acc += (muB * muA * _Grad_uA[m_iComp, d] + muA * muB * _Grad_uB[m_iComp, d]) / (muA + muB) * (_vA - _vB) * inp.Normale[d];  // consistency term
            //                Acc += (muB * muA * _Grad_vA[d] + muA * muB * _Grad_vB[d]) / (muA + muB) * (_uA[m_iComp] - _uB[m_iComp]) * inp.Normale[d];  // symmetry term
            //            }
            //            Acc *= base.m_alpha;
            //            //Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * (2 * muA * muB) / (muA + muB); // penalty term
            //            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * (2 * muA * muB) / (muA + muB); // penalty term

            //            return -Acc;
            //    }
            //    default:
            //    throw new NotImplementedException();
            //}
        }


        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag) {
            if (base.g_Neu_Override == null) {
                return 0.0;
            } else {
                double Acc = 0;
                for (int i = 0; i < base.m_D; i++) {
                    Acc += N[i] * g_Neu_Override(base.m_iComp, X, i);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++

                    double g_D = base.g_Diri(inp.X, inp.time, inp.EdgeTag, m_iComp);
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: {
                    for (int d = 0; d < inp.D; d++) {
                        double nd = inp.Normale[d];
                        //Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                        //Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        Acc += (muA * _Grad_uA[m_iComp, d]) * (_vA) * nd;
                        Acc += (muA * _Grad_vA[d]) * (_uA[m_iComp] - g_D) * nd;
                    }
                    Acc *= base.m_alpha;

                    //Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                    Acc -= muA * (_uA[m_iComp] - g_D) * (_vA - 0) * pnlty;
                    break;
                    //       }
                    //    default:
                    //        throw new NotImplementedException();

                    //}
                    //break;
                }
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.NavierSlip_Linear: {
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: {

                    int D = inp.D;
                    double g_D;

                    for (int dN = 0; dN < D; dN++) {
                        g_D = base.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);

                        for (int dD = 0; dD < D; dD++) {
                            // consistency
                            Acc += muA * (inp.Normale[dN] * _Grad_uA[dN, dD] * inp.Normale[dD]) * (_vA * inp.Normale[m_iComp]) * base.m_alpha;
                            // symmetry
                            Acc += muA * (inp.Normale[m_iComp] * _Grad_vA[dD] * inp.Normale[dD]) * (_uA[dN] - g_D) * inp.Normale[dN] * base.m_alpha;
                        }

                        // penalty
                        Acc -= muA * ((_uA[dN] - g_D) * inp.Normale[dN]) * ((_vA - 0) * inp.Normale[m_iComp]) * pnlty;
                    }

                    if (edgType == IncompressibleBcType.NavierSlip_Linear && m_beta > 0.0) {

                        double[,] P = new double[D, D];
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                double nn = inp.Normale[d1] * inp.Normale[d2];
                                if (d1 == d2) {
                                    P[d1, d2] = 1 - nn;
                                } else {
                                    P[d1, d2] = -nn;
                                }
                            }
                        }

                        // tangential dissipation force term
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                g_D = base.g_Diri(inp.X, inp.time, inp.EdgeTag, d2);
                                Acc -= (m_beta * P[d1, d2] * (_uA[d2] - g_D)) * (P[d1, m_iComp] * _vA) * base.m_alpha;
                            }
                        }

                    }

                    //break;
                    //    }
                    //default:
                    //    throw new NotImplementedException();
                    //}
                    break;
                }
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {
                    // Atmospheric outlet/pressure outflow: hom. Neumann
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    double g_N = g_Neu(inp.X, inp.Normale, inp.EdgeTag);
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: {
                    Acc += muA * g_N * _vA * base.m_alpha;
                    //    break;


                    //    default:
                    //        throw new NotImplementedException();

                    //}

                    break;
                }
                case IncompressibleBcType.Pressure_Dirichlet: {
                    // Dirichlet boundary condition for pressure.
                    // Inner values of velocity gradient are taken, i.e.
                    // no boundary condition for the velocity (resp. velocity gradient) is imposed.                        
                    //switch (base.m_implMode) {
                    //    case ViscosityImplementation.H:
                    //    case ViscosityImplementation.SWIP: {
                    for (int d = 0; d < inp.D; d++) {
                        //Acc += (muA * _Grad_uA[0, d]) * (_vA) * inp.Normale[d];
                        Acc += (muA * _Grad_uA[m_iComp, d]) * (_vA) * inp.Normale[d];
                    }
                    Acc *= base.m_alpha;
                    //            break;
                    //        }
                    //    default:
                    //        throw new NotImplementedException();
                    //}
                    break;
                }
                default:
                    throw new NotImplementedException();
            }

            return -Acc;
        }
    }



    /// <summary>
    /// Solver mode for Swip2 and Swip3 terms.
    /// </summary>
    public enum ViscositySolverMode {
        /// <summary>
        /// Coupled solving for all velocity components.
        /// </summary>
        FullyCoupled,

        /// <summary>
        /// Segregated solving of single velocity components.
        /// </summary>
        Segregated
    }

    /// <summary>
    /// \f[ 
    ///   - \operatorname{div} \left( \mu (\partial_d \vec{u}) \right)
    /// \f]
    /// </summary>
    public class swipViscosity_Term2 : swipViscosityBase {

        private ViscositySolverMode ViscSolverMode;

        /// <summary>
        /// ctor; parameter documentation see <see cref="swipViscosityBase.swipViscosityBase"/>.
        /// </summary>
        public swipViscosity_Term2(double _penalty, MultidimensionalArray PenaltyLengthScales, int iComp, int D, IncompressibleBoundaryCondMap bcmap,
            ViscosityOption _ViscosityMode, ViscositySolverMode ViscSolverMode = ViscositySolverMode.FullyCoupled,
            double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null)
            : base(_penalty, PenaltyLengthScales, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS) {

            this.ViscSolverMode = ViscSolverMode;
        }

        public override double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for (int d = 0; d < cpv.D; d++)
                // we want to:
                //    sum(  \partial_{m_iComp} u_d  * \partial_{d} v, d=0..D-1)
                acc += GradU[d, base.m_iComp] * GradV[d] * Viscosity(cpv.Parameters) * base.m_alpha;
            return acc;
        }


        public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            double muB = this.Viscosity(inp.Parameters_OUT);


            //switch (base.m_implMode) {
            //    case ViscosityImplementation.H: {
            for (int i = 0; i < inp.D; i++) {
                // consistency term
                Acc += 0.5 * (muA * _Grad_uA[i, m_iComp] + muB * _Grad_uB[i, m_iComp]) * (_vA - _vB) * inp.Normale[i];
                // symmetry term
                switch (ViscSolverMode) {
                    case ViscositySolverMode.FullyCoupled:
                    Acc += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[i] - _uB[i]) * inp.Normale[m_iComp];
                    break;
                    case ViscositySolverMode.Segregated:
                    if (i == m_iComp)
                        Acc += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[i] - _uB[i]) * inp.Normale[m_iComp];
                    break;
                    default:
                    throw new NotImplementedException();
                }
            }
            Acc *= base.m_alpha;

            // penalty term
            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax;

            return -Acc;
            //}

            //    case ViscosityImplementation.SWIP: {
            //            for (int i = 0; i < inp.D; i++) {
            //                Acc += (muB * muA * _Grad_uA[i, m_iComp] + muA * muB * _Grad_uB[i, m_iComp]) / (muA + muB) * (_vA - _vB) * inp.Normale[i];  // consistency term
            //                Acc += (muB * muA * _Grad_vA[i] + muA * muB * _Grad_vB[i]) / (muA + muB) * (_uA[i] - _uB[i]) * inp.Normale[m_iComp];  // symmetry term
            //            }
            //            Acc *= base.m_alpha;

            //            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * 2 * muA * muB / (muA + muB); // penalty term

            //            return -Acc;
            //        }
            //    default:
            //        throw new NotImplementedException();
            //}
        }


        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag) {
            if (base.g_Neu_Override == null) {
                //return 0.0;

                throw new NotSupportedException("Neumann BC. for the \\/U^T -- term is problematic!");

            } else {
                double Acc = 0;
                for (int i = 0; i < base.m_D; i++) {
                    Acc += N[i] * g_Neu_Override(i, X, base.m_iComp);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                        // inhom. Dirichlet b.c.
                        // +++++++++++++++++++++
                        double g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp);
                        //switch (base.m_implMode) {
                        //    case ViscosityImplementation.H:
                        //    case ViscosityImplementation.SWIP: {
                                    for (int i = 0; i < inp.D; i++) {
                                        // consistency
                                        Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normale[i];
                                        // symmetry
                                        switch (ViscSolverMode) {
                                            case ViscositySolverMode.FullyCoupled:
                                                Acc += (muA * _Grad_vA[i]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normale[m_iComp];
                                                break;
                                            case ViscositySolverMode.Segregated:
                                                if (i == m_iComp)
                                                    Acc += (muA * _Grad_vA[i]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normale[m_iComp];
                                                break;
                                            default:
                                                throw new NotImplementedException();
                                        }
                                    }
                                    Acc *= base.m_alpha;

                                    // penalty
                                    Acc -= muA * (_uA[m_iComp] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp)) * (_vA - 0) * pnlty;
                        //            break;
                        //        }
                        //    default:
                        //        throw new NotImplementedException();

                        //}
                        break;
                    }
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.NavierSlip_Linear: {
                        //switch (m_implMode) {
                        //    case ViscosityImplementation.H:
                        //    case ViscosityImplementation.SWIP: {

                                    int D = inp.D;
                                    double g_D;

                                    for (int dN = 0; dN < D; dN++) {
                                        for (int dD = 0; dD < D; dD++) {
                                            // consistency
                                            Acc += muA * (inp.Normale[dN] * _Grad_uA[dD, dN] * inp.Normale[dD]) * (_vA * inp.Normale[m_iComp]) * base.m_alpha;
                                            // symmetry
                                            switch (ViscSolverMode) {
                                                case ViscositySolverMode.FullyCoupled:
                                                    g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dD);
                                                    Acc += muA * (inp.Normale[dN] * _Grad_vA[dN] * inp.Normale[m_iComp]) * (_uA[dD] - g_D) * inp.Normale[dD] * base.m_alpha;
                                                    break;
                                                case ViscositySolverMode.Segregated:
                                                default:
                                                    throw new NotImplementedException();
                                            }
                                        }
                                        g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);
                                        // penalty
                                        Acc -= muA * ((_uA[dN] - g_D) * inp.Normale[dN]) * ((_vA - 0) * inp.Normale[m_iComp]) * pnlty;
                                    }

                        //            break;
                        //        }
                        //    default:
                        //        throw new NotImplementedException();
                        //}
                        break;
                }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {
                        //switch (base.m_implMode) {
                        //    case ViscosityImplementation.H:
                        //    case ViscosityImplementation.SWIP: {
                                    if (base.g_Neu_Override == null) {
                                        // Inner values of velocity gradient are taken, i.e.
                                        // no boundary condition for the velocity (resp. velocity gradient) is imposed.
                                        for (int i = 0; i < inp.D; i++) {
                                            Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normale[i];
                                        }
                                    } else {
                                        double g_N = g_Neu(inp.X, inp.Normale, inp.EdgeTag);
                                        Acc += muA * g_N * _vA;
                                    }
                                    Acc *= base.m_alpha;
                                    //break;
                        //        }
                        //    default:
                        //        throw new NotImplementedException();
                        //}
                        break;
                    }
                default:
                    throw new NotSupportedException();
            }

            return -Acc;
        }
    }


    /// <summary>
    /// \f[ 
    ///   \frac{2}{3} \operatorname{div} \left( \mu \myMatrix{I} \operatorname{div} ( \vec{u} )  \right)
    /// \f]
    /// </summary>
    public class swipViscosity_Term3 : swipViscosityBase {

        private ViscositySolverMode ViscSolverMode;

        /// <summary>
        /// ctor; parameter documentation see <see cref="swipViscosityBase.swipViscosityBase"/>.
        /// </summary>
        public swipViscosity_Term3(double _penalty, MultidimensionalArray PenaltyLengthScales, int iComp, int D, IncompressibleBoundaryCondMap bcmap,
            ViscosityOption _ViscosityMode, ViscositySolverMode ViscSolverMode = ViscositySolverMode.FullyCoupled,
            double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null)
            : base(_penalty, PenaltyLengthScales, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS) {

            this.ViscSolverMode = ViscSolverMode;
        }

        public override double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for (int d = 0; d < cpv.D; d++)
                acc -= GradU[d, d] * GradV[base.m_iComp] * Viscosity(cpv.Parameters) * base.m_alpha;
            return acc * (2.0 / 3.0);
        }


        public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.jCellIn, inp.jCellOut);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            double muB = this.Viscosity(inp.Parameters_OUT);


            //switch (base.m_implMode) {
            //    case ViscosityImplementation.H: {
                        for (int i = 0; i < inp.D; i++) {
                            // consistency term
                            Acc += 0.5 * (muA * _Grad_uA[i, i] + muB * _Grad_uB[i, i]) * (_vA - _vB) * inp.Normale[m_iComp];
                            // symmetry term
                            switch (ViscSolverMode) {
                                case ViscositySolverMode.FullyCoupled:
                                    Acc += 0.5 * (muA * _Grad_vA[m_iComp] + muB * _Grad_vB[m_iComp]) * (_uA[i] - _uB[i]) * inp.Normale[i];
                                    break;
                                case ViscositySolverMode.Segregated:
                                    if (i == m_iComp)
                                        Acc += 0.5 * (muA * _Grad_vA[m_iComp] + muB * _Grad_vB[m_iComp]) * (_uA[i] - _uB[i]) * inp.Normale[i];
                                    break;
                                default:
                                    throw new NotImplementedException();
                            }
                        }
                        Acc *= base.m_alpha;

                        // penalty term
                        double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
                        Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax;

                        return Acc * (2.0 / 3.0);
            //        }

            //    case ViscosityImplementation.SWIP: {
            //            for (int i = 0; i < inp.D; i++) {
            //                Acc += (muB * muA * _Grad_uA[i, i] + muA * muB * _Grad_uB[i, i]) / (muA + muB) * (_vA - _vB) * inp.Normale[m_iComp];  // consistency term
            //                Acc += (muB * muA * _Grad_vA[m_iComp] + muA * muB * _Grad_vB[m_iComp]) / (muA + muB) * (_uA[i] - _uB[i]) * inp.Normale[i];  // symmetry term
            //            }
            //            Acc *= base.m_alpha;

            //            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * 2 * muA * muB / (muA + muB); // penalty term

            //            return Acc * (2.0 / 3.0);
            //        }
            //    default:
            //        throw new NotImplementedException();
            //}
        }



        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag) {
            if (base.g_Neu_Override == null) {
                //return 0.0;

                throw new NotSupportedException("Neumann BC. for the \\/U^T -- term is problematic!");

            } else {
                double Acc = 0;
                for (int i = 0; i < base.m_D; i++) {
                    Acc += N[m_iComp] * g_Neu_Override(i, X, i);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                        // inhom. Dirichlet b.c.
                        // +++++++++++++++++++++
                        double g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp);
                        //switch (base.m_implMode) {
                        //    case ViscosityImplementation.H:
                        //    case ViscosityImplementation.SWIP: {
                                    for (int i = 0; i < inp.D; i++) {
                                        // consistency
                                        Acc += (muA * _Grad_uA[i, i]) * (_vA) * inp.Normale[m_iComp];
                                        // symmetry
                                        switch (ViscSolverMode) {
                                            case ViscositySolverMode.FullyCoupled:
                                                Acc += (muA * _Grad_vA[m_iComp]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normale[i];
                                                break;
                                            case ViscositySolverMode.Segregated:
                                                if (i == m_iComp)
                                                    Acc += (muA * _Grad_vA[m_iComp]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normale[i];
                                                break;
                                            default:
                                                throw new NotImplementedException();
                                        }
                                    }
                                    Acc *= base.m_alpha;

                                    // penalty
                                    Acc -= muA * (_uA[m_iComp] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp)) * (_vA - 0) * pnlty;
                        //            break;
                        //        }

                        //    default:
                        //        throw new NotImplementedException();

                        //}
                        break;
                    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {
                        //switch (base.m_implMode) {
                        //    case ViscosityImplementation.H:
                        //    case ViscosityImplementation.SWIP: {
                                    if (base.g_Neu_Override == null) {
                                        // Inner values of velocity gradient are taken, i.e.
                                        // no boundary condition for the velocity (resp. velocity gradient) is imposed.
                                        for (int i = 0; i < inp.D; i++) {
                                            Acc += (muA * _Grad_uA[i, i]) * (_vA) * inp.Normale[m_iComp];
                                        }
                                    } else {
                                        double g_N = g_Neu(inp.X, inp.Normale, inp.EdgeTag);
                                        Acc += muA * g_N * _vA;
                                    }
                                    Acc *= base.m_alpha;
                        //            break;
                        //        }
                        //    default:
                        //        throw new NotImplementedException();
                        //}
                        break;
                    }
                default:
                    throw new NotSupportedException();
            }

            return Acc * (2.0 / 3.0);
        }
    }
}

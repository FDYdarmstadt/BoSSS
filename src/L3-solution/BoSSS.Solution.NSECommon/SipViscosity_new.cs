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
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// base class for viscosity terms
    /// </summary>
    public abstract class SipViscosityNewBase : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, IEquationComponentCoefficient, IEquationComponentChecking, ISupportsJacobianComponent {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="_penaltyBase"></param>
        /// <param name="iComp">
        /// component index
        /// </param>
        /// <param name="D">
        /// spatial dimension.
        /// </param>
        /// <param name="ViscosityTerms"></param>        
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
        /// Flag for activating/deactivating the optimized fluxes
        /// <param name="ignoreVectorized"></param>
        protected SipViscosityNewBase(
                                    double _penaltyBase,
                                    int iComp, int D, ViscosityTermsSwitch ViscosityTerms, IncompressibleBoundaryCondMap bcmap, /*TermSwitches ViscosityTerms,*/
                                    ViscosityOption _ViscosityMode, double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null, bool ignoreVectorized = false) {
            //Func<double, int, int, MultidimensionalArray, double> ComputePenalty = null) {
            this.m_penalty_base = _penaltyBase;
            //this.m_ComputePenalty = ComputePenalty;
            this.m_iComp = iComp;
            this.m_D = D;
            //this.cj = PenaltyLengthScales;
            velFunction = D.ForLoop(d => bcmap.bndFunction[VariableNames.Velocity_d(d)]);
            EdgeTag2Type = bcmap.EdgeTag2Type;
            this.m_PhysicsMode = bcmap.PhysMode;
            this.m_ignoreVectorized = ignoreVectorized;
            this.m_ViscosityMode = _ViscosityMode;
            m_ViscositySwitches = ViscosityTerms;

            m_useArgumentsForViscosity = true;
        //m_ignoreVectorized = true;
            switch (_ViscosityMode) {
                case ViscosityOption.ConstantViscosity:
                    if (double.IsNaN(constantViscosityValue))
                        throw new ArgumentException("constantViscosityValue is missing!");
                    this.m_constantViscosityValue = constantViscosityValue;
                    break;
                case ViscosityOption.ConstantViscosityDimensionless:
                    if (double.IsNaN(reynolds))
                        throw new ArgumentException("Reynolds number is missing!");
                    this.m_reynolds = reynolds;
                    break;
                case ViscosityOption.VariableViscosity:
                    this.m_EoS = EoS;
                    break;
                case ViscosityOption.VariableViscosityDimensionless:
                    if (double.IsNaN(reynolds))
                        throw new ArgumentException("Reynolds number is missing!");
                    this.m_reynolds = reynolds;
                    this.m_EoS = EoS;
                    break;
                default:
                    throw new NotImplementedException();

            }

            switch (bcmap.PhysMode) {
                case PhysicsMode.Incompressible: 
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(m_D), new string[] { });
                    break;
                case PhysicsMode.MixtureFraction:
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.MixtureFraction);
                ScalarFunction = bcmap.bndFunction[VariableNames.MixtureFraction]; // 
                break;
                case PhysicsMode.LowMach:
                case PhysicsMode.Combustion:
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.Temperature); // Viscosity is only dependent on Temperature.
                ScalarFunction = bcmap.bndFunction[VariableNames.Temperature];
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
        /// spatial dimension
        /// </summary>
        protected int m_D;

        /// <summary>
        /// Dirichlet boundary values; 
        ///  - 1st index: spatial dimension 
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[][] velFunction;


        /// <summary>
        /// Dirichlet boundary values of temperature; 
        ///  - 1st index: spatial dimension 
        ///  - 2nd index: edge tag
        /// </summary>
        protected Func<double[], double, double>[] ScalarFunction;


        /// <summary>
        /// dissipation coefficient for the effective wall force
        /// </summary>
        protected double m_beta = 0.0;

        /// <summary>
        /// slip-length for the navier-slip BC
        /// </summary>
        protected MultidimensionalArray Lslip;


        //Func<double, int, int, MultidimensionalArray, double> m_ComputePenalty;

        /// <summary>
        /// see <see cref="ViscosityOption"/>
        /// </summary>
        protected ViscosityOption m_ViscosityMode;

        /// <summary>
        /// in the case of constant viscosity, the value of the viscosity
        /// </summary>
        double m_constantViscosityValue = double.NaN;

        /// <summary>
        /// Reynolds number for dimensionless formulation.
        /// </summary>
        protected double m_reynolds = double.NaN;

        /// <summary>
        /// Optional material law for calculating the viscosity
        /// as a function of temperature / level-set.
        /// </summary>
        MaterialLaw m_EoS = null;

        /// <summary>
        /// Optional switch for ignoring the vectorized versions of fluxes. May be used for testing 
        /// </summary>
        bool m_ignoreVectorized = false;

        /// <summary>
        /// Really ugly switch used to decide which values are to be used in the arguments of viscosity, either arguments or parameters.
        /// </summary>
       protected bool m_useArgumentsForViscosity = true;

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
                        if (Parameters.Length > 1) {
                            return m_EoS.GetViscosity(Parameters);
                        } else {
                            return m_EoS.GetViscosity(Parameters[0]);
                        }
                    }
                case ViscosityOption.VariableViscosityDimensionless:
                    if (m_EoS == null) {
                        return (Parameters[0] / m_reynolds);
                    } else {
                        if (Parameters.Length > 1) {
                            return (m_EoS.GetViscosity(Parameters) / m_reynolds);
                        } else {
                            return (m_EoS.GetViscosity(Parameters[0]) / m_reynolds);
                        }
                     
                    }
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            switch (m_ViscosityMode) {
                case ViscosityOption.ConstantViscosity:
                case ViscosityOption.ConstantViscosityDimensionless:
                    return new IEquationComponent[] { this };

                case ViscosityOption.VariableViscosity:
                case ViscosityOption.VariableViscosityDimensionless:
                    var DivergenceDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
                    var DivergenceDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
                    return new IEquationComponent[] { DivergenceDerivEdg, DivergenceDerivVol };

                default:
                    throw new NotImplementedException();
            }
        }



        /// <summary>
        /// Cell-wise length scales for the penalty computation.
        /// </summary>
        MultidimensionalArray cj;

        /// <summary>
        /// base multiplier for the penalty computation
        /// </summary>
        protected double m_penalty_base;





        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        /// <param name="cs"></param>
        /// <param name="DomainDGdeg"></param>
        /// <param name="TestDGdeg"></param>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            m_D = cs.GrdDat.SpatialDimension;
            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;

            if (cs.UserDefinedValues.Keys.Contains("SlipLengths"))
                Lslip = (MultidimensionalArray)cs.UserDefinedValues["SlipLengths"];
            // Set the Reynolds number to a user defined value contained in the CoefficientSet cs
            // Useful in case that the Reynolds number changes during a simulation...
            if (cs.UserDefinedValues.Keys.Contains("Reynolds"))
                m_reynolds = (double)cs.UserDefinedValues["Reynolds"];

            if (cs.UserDefinedValues.Keys.Contains("VelocityMultiplier"))
                VelocityMultiplier = (double)cs.UserDefinedValues["VelocityMultiplier"];
        }
        /// <summary>
        /// 
        /// </summary>
        public double VelocityMultiplier = 1.0;

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected double penalty(IGridData g, int jCellIn, int jCellOut, int iEdge) {

            double penaltySizeFactor_A = 1.0 / cj[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(m_penalty));
            Debug.Assert(!double.IsInfinity(m_penalty));

            return penaltySizeFactor * m_penalty * m_penalty_base;
        }
        /// <summary>
        /// A list with all the arguments of this flux
        /// </summary>
        string[] m_ArgumentOrdering;

        /// <summary>
        /// returns the velocity vector;
        /// </summary>
        public virtual IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
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
                            case PhysicsMode.Incompressible:
                            case PhysicsMode.MixtureFraction:
                                    return new string[] { };
                                case PhysicsMode.LowMach:
                                case PhysicsMode.Combustion:
                                    return new string[] { /*VariableNames.Temperature0*/ };
                                case PhysicsMode.Multiphase:
                                    return new string[] { VariableNames.LevelSet };
                                case PhysicsMode.RANS:
                                    return new string[] { "k0", "omega0" }; // TODO generalize this for other turbulence models
                                case PhysicsMode.Viscoelastic:
                          
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
        /// Velocity Dirichlet boundary value: the given velocity at the boundary.
        /// </summary>
        protected double g_Diri(double[] X, double time, int EdgeTag, int d) {
            if (this.g_Diri_Override == null) {
                return this.velFunction[d][EdgeTag](X, time) * VelocityMultiplier;
            } else {
                return g_Diri_Override(d, X);
            }
        }



        /// <summary>
        /// Temperature Dirichlet boundary value: the given temperature at the boundary.
        /// </summary>
        protected double T_Diri(double[] X, double time, int EdgeTag, int d) {        
            return this.ScalarFunction[EdgeTag](X, time);
        }


        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV/* | TermActivationFlags.V | TermActivationFlags.GradV*/);
            }
        }

        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV;

            }
        }

        bool IEquationComponentChecking.IgnoreVectorizedImplementation {
            get {
                return m_ignoreVectorized;
            }
        }




        abstract public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV);

        abstract public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB);

        abstract public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA);



        public ViscosityTermsSwitch m_ViscositySwitches;
        public ViscosityTermsSwitch ViscosityTermsSwitches {
            get {
                return m_ViscositySwitches;
            }
        }




        public ViscositySolverMode ViscSolverMode;


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


    }

    /// <summary>
    /// Implements the operator 
    /// \f[ 
    ///   -\operatorname{div} \left( \mu \nabla \vec{u} \right) +   - \operatorname{div} \left( \mu (\partial_d \vec{u}) \right) +\frac{2}{3} \operatorname{div} \left( \mu \myMatrix{I} \operatorname{div} ( \vec{u} )  \right)
    /// \f]
    /// for a variable \mu value
    /// </summary>
    public class SipViscosity_Variable : SipViscosityNewBase, INonlinVolumeForm_GradV,
        INonlinEdgeForm_GradV,
        INonlinEdgeForm_V {

        /// <summary>
        /// ctor; parameter documentation see <see cref="SipViscosityNewBase.SipViscosityNewBase"/>.
        /// </summary>
        public SipViscosity_Variable(double _penalty, int iComp, int D, ViscosityTermsSwitch myviscosityTerms, IncompressibleBoundaryCondMap bcmap,
                                   ViscosityOption _ViscosityMode, double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null, bool ignoreVectorized = false)
            //Func<double, int, int, MultidimensionalArray, double> ComputePenalty = null)
            : base(_penalty, iComp, D, myviscosityTerms, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS, ignoreVectorized) {
            this.ViscSolverMode = ViscositySolverMode.FullyCoupled;
        }

        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu_1(double[] X, double[] N, int EdgeTag) {
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

        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu_2(double[] X, double[] N, int EdgeTag) {
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
        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu_3(double[] X, double[] N, int EdgeTag) {
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

        public double BoundaryEdgeForm_1(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;
            double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);//, inp.GridDat.Cells.cj);
                                                                                     //double muA = this.Viscosity(inp.Parameters_IN);


            int NOargs = _uA.Length;
            double[] viscArgsIN = m_useArgumentsForViscosity ? _uA.GetSubVector(inp.D, NOargs - inp.D) : inp.Parameters_IN;
            double muA = this.Viscosity(viscArgsIN);


            //double muA = this.Viscosity(_uA.GetSubVector(inp.D, 1));
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++

                    double g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, m_iComp);

                    for (int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        Acc += (muA * _Grad_uA[m_iComp, d]) * (_vA) * nd;
                        Acc += (muA * _Grad_vA[d]) * (_uA[m_iComp] - g_D) * nd;
                    }
                    Acc *= base.m_alpha;

                    double temperatureOut = this.T_Diri(inp.X, inp.time, inp.EdgeTag, m_iComp);
                    double[] viscArgsOut = new double[] { temperatureOut };
                    double muB = this.Viscosity(viscArgsOut);

                    // penalty term          
                    double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;

                    Acc -= muMax * (_uA[m_iComp] - g_D) * (_vA - 0) * pnlty;
                    break;
                }
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry: {

                        int D = inp.D;
                        double g_D;

                        for (int dN = 0; dN < D; dN++) {

                            g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);

                            for (int dD = 0; dD < D; dD++) {
                                // consistency
                                Acc += muA * (inp.Normal[dN] * _Grad_uA[dN, dD] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                                // symmetry
                                Acc += muA * (inp.Normal[m_iComp] * _Grad_vA[dD] * inp.Normal[dD]) * (_uA[dN] - g_D) * inp.Normal[dN] * base.m_alpha;
                            }

                            // penalty
                            Acc -= muA * ((_uA[dN] - g_D) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[m_iComp]) * pnlty;
                        }
                        break;
                    }
                case IncompressibleBcType.NavierSlip_Linear: {

                        double ls = Lslip[inp.jCellIn];
                        if (ls == 0.0)
                            goto case IncompressibleBcType.Velocity_Inlet;

                        if (ls > 0)
                            m_beta = muA / ls;


                        int D = inp.D;
                        double g_D;

                        for (int dN = 0; dN < D; dN++) {
                            g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);

                            for (int dD = 0; dD < D; dD++) {
                                // consistency
                                Acc += muA * (inp.Normal[dN] * _Grad_uA[dN, dD] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                                // symmetry
                                Acc += muA * (inp.Normal[m_iComp] * _Grad_vA[dD] * inp.Normal[dD]) * (_uA[dN] - g_D) * inp.Normal[dN] * base.m_alpha;
                            }

                            // penalty
                            Acc -= muA * ((_uA[dN] - g_D) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[m_iComp]) * pnlty;
                        }


                        double[,] P = new double[D, D];
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                double nn = inp.Normal[d1] * inp.Normal[d2];
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
                                g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, d2);
                                Acc -= (m_beta * P[d1, d2] * (_uA[d2] - g_D)) * (P[d1, m_iComp] * _vA) * base.m_alpha;
                            }
                        }

                        break;
                    }
                //case IncompressibleBcType.NavierSlip_localized: {

                //    double ls = Lslip[inp.jCellIn];
                //    if(ls > 0.0) {
                //        m_beta = muA / ls;
                //        goto case IncompressibleBcType.NavierSlip_Linear;
                //    } else {
                //        goto case IncompressibleBcType.Velocity_Inlet;
                //    }
                //}
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet: {
                        // Atmospheric outlet/pressure outflow: hom. Neumann
                        // +++++++++++++++++++++++++++++++++++++++++++++++++
                        double g_N = g_Neu_1(inp.X, inp.Normal, inp.EdgeTag);

                        Acc += muA * g_N * _vA * base.m_alpha;

                        break;
                    }
                case IncompressibleBcType.Pressure_Dirichlet: {
                        // Dirichlet boundary condition for pressure.
                        // Inner values of velocity gradient are taken, i.e.
                        // no boundary condition for the velocity (resp. velocity gradient) is imposed.                        

                        for (int d = 0; d < inp.D; d++) {
                            //Acc += (muA * _Grad_uA[0, d]) * (_vA) * inp.Normale[d];
                            Acc += (muA * _Grad_uA[m_iComp, d]) * (_vA) * inp.Normal[d];
                        }
                        Acc *= base.m_alpha;

                        break;
                    }
                default:
                    throw new NotImplementedException();
            }

            return -Acc;
        }

        public double BoundaryEdgeForm_2(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);
            int NOargs = _uA.Length;
            double[] viscArgsIN = m_useArgumentsForViscosity ? _uA.GetSubVector(inp.D, NOargs - inp.D) : inp.Parameters_IN;
            double muA = this.Viscosity(viscArgsIN);

            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                        // inhom. Dirichlet b.c.
                        // +++++++++++++++++++++


                        for (int i = 0; i < inp.D; i++) {
                            // consistency
                            Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
                            // symmetry
                            switch (ViscSolverMode) {
                                case ViscositySolverMode.FullyCoupled:
                                    Acc += (muA * _Grad_vA[i]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[m_iComp];
                                    break;
                                case ViscositySolverMode.Segregated:
                                    if (i == m_iComp)
                                        Acc += (muA * _Grad_vA[i]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[m_iComp];
                                    break;
                                default:
                                    throw new NotImplementedException();
                            }
                        }
                        Acc *= base.m_alpha;

                    double temperatureOut = this.T_Diri(inp.X, inp.time, inp.EdgeTag, m_iComp);
                    double[] viscArgsOut = new double[] { temperatureOut };
                    double muB = this.Viscosity(viscArgsOut);

                    // penalty term          
                    double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;


                    // penalty
                    Acc -= muMax * (_uA[m_iComp] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp)) * (_vA - 0) * pnlty;

                        break;
                    }
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry: {

                        int D = inp.D;
                        double g_D;

                        for (int dN = 0; dN < D; dN++) {
                            for (int dD = 0; dD < D; dD++) {
                                // consistency
                                Acc += muA * (inp.Normal[dN] * _Grad_uA[dD, dN] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                                // symmetry
                                switch (ViscSolverMode) {
                                    case ViscositySolverMode.FullyCoupled:
                                        g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dD);
                                        Acc += muA * (inp.Normal[dN] * _Grad_vA[dN] * inp.Normal[m_iComp]) * (_uA[dD] - g_D) * inp.Normal[dD] * base.m_alpha;
                                        break;
                                    case ViscositySolverMode.Segregated:
                                    default:
                                        throw new NotImplementedException();
                                }
                            }
                            g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);
                            // penalty
                            Acc -= muA * ((_uA[dN] - g_D) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[m_iComp]) * pnlty;
                        }

                        break;
                    }
                case IncompressibleBcType.NavierSlip_Linear: {

                        double ls = Lslip[inp.jCellIn];
                        if (ls == 0.0)
                            goto case IncompressibleBcType.Velocity_Inlet;
                        else
                            goto case IncompressibleBcType.FreeSlip;


                    }
                //case IncompressibleBcType.NavierSlip_localized: {

                //        double ls = Lslip[inp.jCellIn];
                //        if(ls > 0.0) {
                //            m_beta = muA / ls;
                //            goto case IncompressibleBcType.NavierSlip_Linear;
                //        } else {
                //            goto case IncompressibleBcType.Velocity_Inlet;
                //        }
                //    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet: {

                        if (base.g_Neu_Override == null) {
                            // Inner values of velocity gradient are taken, i.e.
                            // no boundary condition for the velocity (resp. velocity gradient) is imposed.
                            for (int i = 0; i < inp.D; i++) {
                                Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
                            }
                        } else {
                            double g_N = g_Neu_2(inp.X, inp.Normal, inp.EdgeTag);
                            Acc += muA * g_N * _vA;
                        }
                        Acc *= base.m_alpha;

                        break;
                    }
                default:
                    throw new NotSupportedException("Not Implemented boundary condition");
            }

            return -Acc;
        }

        public double BoundaryEdgeForm_3(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;
            double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);//, inp.GridDat.Cells.cj);
            int NOargs = _uA.Length;
            double[] viscArgsIN = m_useArgumentsForViscosity ? _uA.GetSubVector(inp.D, NOargs - inp.D) : inp.Parameters_IN;
            double muA = this.Viscosity(viscArgsIN);

            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                        // inhom. Dirichlet b.c.
                        // +++++++++++++++++++++                      

                        for (int i = 0; i < inp.D; i++) {
                            // consistency
                            Acc += (muA * _Grad_uA[i, i]) * (_vA) * inp.Normal[m_iComp];
                            // symmetry
                            switch (ViscSolverMode) {
                                case ViscositySolverMode.FullyCoupled:
                                    Acc += (muA * _Grad_vA[m_iComp]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[i];
                                    break;
                                case ViscositySolverMode.Segregated:
                                    if (i == m_iComp)
                                        Acc += (muA * _Grad_vA[m_iComp]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[i];
                                    break;
                                default:
                                    throw new NotImplementedException();
                            }
                        }
                        Acc *= base.m_alpha;


                    double temperatureOut = this.T_Diri(inp.X, inp.time, inp.EdgeTag, m_iComp);
                    double[] viscArgsOut = new double[] { temperatureOut };
                    double muB = this.Viscosity(viscArgsOut);

                    // penalty term          
                    double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
                    // penalty
                    Acc -= muMax * (_uA[m_iComp] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp)) * (_vA - 0) * pnlty;

                        break;
                    }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.ScalarDirichlet_PressureOutlet: {

                        if (base.g_Neu_Override == null) {
                            // Inner values of velocity gradient are taken, i.e.
                            // no boundary condition for the velocity (resp. velocity gradient) is imposed.
                            for (int i = 0; i < inp.D; i++) {
                                Acc += (muA * _Grad_uA[i, i]) * (_vA) * inp.Normal[m_iComp];
                            }
                        } else {
                            double g_N = g_Neu_3(inp.X, inp.Normal, inp.EdgeTag);
                            Acc += muA * g_N * _vA;
                        }
                        Acc *= base.m_alpha;

                        break;
                    }
                default:
                    throw new NotSupportedException();
            }

            return Acc * (2.0 / 3.0);
        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double acc1 = 0;
            double acc2 = 0;
            double acc3 = 0;
            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_u) != 0) {
                acc1 = BoundaryEdgeForm_1(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
            }

            // Add SWIPViscosity2 term
            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_uT) != 0) {
                acc2 = BoundaryEdgeForm_2(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
            }
            // Add SWIPViscosity3 term
            if ((m_ViscositySwitches & ViscosityTermsSwitch.divU) != 0) {
                acc3 = BoundaryEdgeForm_3(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
            }

            return (acc1 + acc2 + acc3);
        }


        public override double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double acc1 = 0;
            double acc2 = 0;
            double acc3 = 0;

            double pnlty = this.penalty(inp.GridDat, inp.jCellIn, inp.jCellOut, inp.iEdge);
            //double muA = this.Viscosity(inp.Parameters_IN);
            //double muB = this.Viscosity(inp.Parameters_OUT);

            int NOargsIN = _uA.Length;
            int NOargsOut = _uB.Length;
            double[] viscArgsIN = m_useArgumentsForViscosity ? _uA.GetSubVector(inp.D, NOargsIN - inp.D) : inp.Parameters_IN;
            double[] viscArgsOUT = m_useArgumentsForViscosity ? _uB.GetSubVector(inp.D, NOargsOut - inp.D) : inp.Parameters_IN;
            double muA = this.Viscosity(viscArgsIN);
            double muB = this.Viscosity(viscArgsOUT);

            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;

            // Add SWIPViscosity1 term

            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_u) != 0) {
                for (int i = 0; i < inp.D; i++) {
                    acc1 += 0.5 * (muA * _Grad_uA[m_iComp, i] + muB * _Grad_uB[m_iComp, i]) * (_vA - _vB) * inp.Normal[i];  // consistency term  
                    acc1 += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[m_iComp] - _uB[m_iComp]) * inp.Normal[i];  // symmetry term
                }
                acc1 *= base.m_alpha;
                acc1 -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax; // penalty term          
            }


            // Add SWIPViscosity2 term
            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_uT) != 0) {
                for (int i = 0; i < inp.D; i++) {
                    // consistency term
                    acc2 += 0.5 * (muA * _Grad_uA[i, m_iComp] + muB * _Grad_uB[i, m_iComp]) * (_vA - _vB) * inp.Normal[i];
                    // symmetry term
                    switch (ViscSolverMode) {
                        case ViscositySolverMode.FullyCoupled:
                            acc2 += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[i] - _uB[i]) * inp.Normal[m_iComp];
                            break;
                        case ViscositySolverMode.Segregated:
                            if (i == m_iComp)
                                acc2 += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[i] - _uB[i]) * inp.Normal[m_iComp];
                            break;
                        default:
                            throw new NotImplementedException();
                    }
                }
                acc2 *= base.m_alpha;
                // penalty term
                acc2 -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax;
            }

            // Add SWIPViscosity3 term
            if ((m_ViscositySwitches & ViscosityTermsSwitch.divU) != 0) {

                for (int i = 0; i < inp.D; i++) {
                    // consistency term
                    acc3 += 0.5 * (muA * _Grad_uA[i, i] + muB * _Grad_uB[i, i]) * (_vA - _vB) * inp.Normal[m_iComp];
                    // symmetry term
                    switch (ViscSolverMode) {
                        case ViscositySolverMode.FullyCoupled:
                            acc3 += 0.5 * (muA * _Grad_vA[m_iComp] + muB * _Grad_vB[m_iComp]) * (_uA[i] - _uB[i]) * inp.Normal[i];
                            break;
                        case ViscositySolverMode.Segregated:
                            if (i == m_iComp)
                                acc3 += 0.5 * (muA * _Grad_vA[m_iComp] + muB * _Grad_vB[m_iComp]) * (_uA[i] - _uB[i]) * inp.Normal[i];
                            break;
                        default:
                            throw new NotImplementedException();
                    }
                }

                acc3 *= base.m_alpha;

                // penalty term
                acc3 -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax;
                acc3 *= (2.0 / 3.0);

            }


            return (-acc1 + -acc2 + acc3);
        }


        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc1 = 0;
            double acc2 = 0;
            double acc3 = 0;
            int NOargs = U.Length;
            double[] viscArgs = m_useArgumentsForViscosity ? U.GetSubVector(cpv.D, NOargs - cpv.D) : cpv.Parameters;
            double visc = this.Viscosity(viscArgs);


            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_u) != 0) {
                for (int d = 0; d < cpv.D; d++)
                    acc1 -= GradU[m_iComp, d] * GradV[d];
                acc1 *= visc * base.m_alpha;
            }


            // Add SWIPViscosity2 term
            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_uT) != 0) {
                for (int d = 0; d < cpv.D; d++)
                    acc2 -= GradU[d, base.m_iComp] * GradV[d];
                acc2 *= visc * base.m_alpha;

            }
            // Add SWIPViscosity3 term
            if ((m_ViscositySwitches & ViscosityTermsSwitch.divU) != 0) {
                for (int d = 0; d < cpv.D; d++)
                    acc3 -= GradU[d, d] * GradV[base.m_iComp] * visc * base.m_alpha;
                acc3 *= (2.0 / 3.0);
            }

            double acc = (-acc1 - acc2 + acc3);
            return acc;
        }




        void INonlinVolumeForm_GradV.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) {

            int NumofCells = prm.Len;
            int NumOfNodes = f.GetLength(1); // no of nodes per cell
            Debug.Assert(f.GetLength(0) == NumofCells);
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            int _NOArgumentsForViscosity = this.ArgumentOrdering.Count - m_D; // Dont include the 2  or 3 velocity components
            double[] Parameters = m_useArgumentsForViscosity ? new double[_NOArgumentsForViscosity]: new double[_NOParams];

            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_u) != 0) {

                for (int cell = 0; cell < NumofCells; cell++) { // loop over cells...
                    for (int node = 0; node < NumOfNodes; node++) { // loop over nodes... 

                        //for (int np = 0; np < _NOParams; np++) {
                        //    Parameters[np] = prm.ParameterVars[np][cell, node];
                        //}
                        if (m_useArgumentsForViscosity) {
                            for (int np = 0; np < _NOArgumentsForViscosity; np++) {
                                Parameters[np] = U[m_D+np][cell,node];
                            }
                        } else {
                            for (int np = 0; np < _NOParams; np++) {
                                Parameters[np] = prm.ParameterVars[np][cell, node];
                            }
                        }


                        double viscosity = Viscosity(Parameters) * base.m_alpha;
                        Debug.Assert(!double.IsNaN(viscosity));
                        Debug.Assert(!double.IsInfinity(viscosity));

                        for (int d = 0; d < prm.GridDat.SpatialDimension; d++) {
                            // Add SWIPViscosity1 term
                            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_u) != 0) {
                                f[cell, node, d] += viscosity * GradU[m_iComp][cell, node, d];
                            }
                            // Add SWIPViscosity2 term
                            if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_uT) != 0) {
                                f[cell, node, d] += viscosity * GradU[d][cell, node, m_iComp];
                            }
                            // Add SWIPViscosity3 term
                            if ((m_ViscositySwitches & ViscosityTermsSwitch.divU) != 0) {
                                f[cell, node, m_iComp] -= viscosity * GradU[d][cell, node, d] * (2.0 / 3.0);
                            }
                        }
                    }
                }
            }

        }

        void INonlinInnerEdgeForm_GradV.NonlinInternalEdge_GradV(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout, MultidimensionalArray fIN, MultidimensionalArray fOT) { // OK :)
            int NumOfEdges = efp.Len;
            Debug.Assert(fIN.GetLength(0) == NumOfEdges);
            Debug.Assert(fOT.GetLength(0) == NumOfEdges);
            int NumOfNodes = fIN.GetLength(1); // no of nodes per cell
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            int _NOArgumentsForViscosity = this.ArgumentOrdering.Count - m_D; // Dont include the 2  or 3 velocity components
            double[] ParametersIN = m_useArgumentsForViscosity ? new double[_NOArgumentsForViscosity] : new double[_NOParams];
            double[] ParametersOT = m_useArgumentsForViscosity ? new double[_NOArgumentsForViscosity] : new double[_NOParams];


            for (int edges = 0; edges < NumOfEdges; edges++) { // loop over edges...
                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                    double uJump = 0.5 * (Uin[m_iComp][edges, node] - Uout[m_iComp][edges, node]);

                    //for (int np = 0; np < _NOParams; np++) {
                    //    ParametersIN[np] = efp.ParameterVars_IN[np][edges, node];
                    //    ParametersOT[np] = efp.ParameterVars_OUT[np][edges, node];

                    //}

                    if (m_useArgumentsForViscosity) {
                        for (int np = 0; np < _NOArgumentsForViscosity; np++) {
                            ParametersIN[np] = Uin[m_D + np][edges, node];
                            ParametersOT[np] = Uout[m_D + np][edges, node];
                        }
                    } else {
                        for (int np = 0; np < _NOParams; np++) {
                            ParametersIN[np] = efp.ParameterVars_IN[np][edges, node];
                            ParametersOT[np] = efp.ParameterVars_OUT[np][edges, node];
                        }
                    }
                
                    double viscosityIN = Viscosity(ParametersIN);
                    double viscosityOT = Viscosity(ParametersOT);
                    Debug.Assert(!double.IsNaN(viscosityIN));
                    Debug.Assert(!double.IsNaN(viscosityOT));
                    Debug.Assert(!double.IsInfinity(viscosityIN));
                    Debug.Assert(!double.IsInfinity(viscosityOT));



                    double fluxIn = viscosityIN * uJump;
                    double fluxOut = viscosityOT * uJump;


                    if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_u) != 0) {
                        for (int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                            double n = efp.Normals[edges, node, d];

                            fIN[edges, node, d] -= fluxIn * n;
                            fOT[edges, node, d] -= fluxOut * n;
                            fIN[edges, node, d] *= base.m_alpha;
                            fOT[edges, node, d] *= base.m_alpha;
                        }
                    }


                    // Add SWIPViscosity2 term
                    if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_uT) != 0) {
                        double n = efp.Normals[edges, node, m_iComp];
                        switch (ViscSolverMode) {
                            case ViscositySolverMode.FullyCoupled:
                                for (int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                                    fIN[edges, node, d] -= 0.5 * viscosityIN * (Uin[d][edges, node] - Uout[d][edges, node]) * n * base.m_alpha;
                                    fOT[edges, node, d] -= 0.5 * viscosityOT * (Uin[d][edges, node] - Uout[d][edges, node]) * n * base.m_alpha;
                                }
                                break;
                            case ViscositySolverMode.Segregated:
                                fIN[edges, node, m_iComp] -= 0.5 * viscosityIN * (Uin[m_iComp][edges, node] - Uout[m_iComp][edges, node]) * n * base.m_alpha;
                                fOT[edges, node, m_iComp] -= 0.5 * viscosityOT * (Uin[m_iComp][edges, node] - Uout[m_iComp][edges, node]) * n * base.m_alpha;
                                break;
                            default:
                                throw new NotImplementedException();

                        }
                    }
                    // Add SWIPViscosity3 term
                    if ((m_ViscositySwitches & ViscosityTermsSwitch.divU) != 0) {
                        double n;
                        switch (ViscSolverMode) {
                            case ViscositySolverMode.FullyCoupled:
                                for (int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                                    n = efp.Normals[edges, node, d];
                                    fIN[edges, node, m_iComp] -= 0.5 * viscosityIN * (Uin[d][edges, node] - Uout[d][edges, node]) * n * base.m_alpha * (-2.0 / 3.0);
                                    fOT[edges, node, m_iComp] -= 0.5 * viscosityOT * (Uin[d][edges, node] - Uout[d][edges, node]) * n * base.m_alpha * (-2.0 / 3.0);
                                }
                                break;
                            case ViscositySolverMode.Segregated:
                                n = efp.Normals[edges, node, m_iComp];
                                fIN[edges, node, m_iComp] -= 0.5 * viscosityIN * (Uin[m_iComp][edges, node] - Uout[m_iComp][edges, node]) * n * base.m_alpha * (-2.0 / 3.0);
                                fOT[edges, node, m_iComp] -= 0.5 * viscosityOT * (Uin[m_iComp][edges, node] - Uout[m_iComp][edges, node]) * n * base.m_alpha * (-2.0 / 3.0);
                                break;
                            default:
                                throw new NotImplementedException();

                        }
                    }

                }
            }




        }

        void INonlinBoundaryEdgeForm_GradV.NonlinBoundaryEdge_GradV(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin, MultidimensionalArray f) {
            int L = efp.Len;
            Debug.Assert(f.GetLength(0) == L);
            int K = f.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            Debug.Assert(_NOParams == efp.ParameterVars_IN.Length);
            int _NOargs = this.ArgumentOrdering.Count;
            Debug.Assert(_NOargs == Uin.Length);
            Debug.Assert(_NOargs == GradUin.Length);

            CommonParamsBnd cpv;
            cpv.GridDat = efp.GridDat;
            cpv.Parameters_IN = new double[_NOParams];
            cpv.Normal = new Vector(D); ;
            cpv.X = new Vector(D);
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            byte[] EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for (int l = 0; l < L; l++) { // loop over edges ...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for (int k = 0; k < K; k++) { // loop over nodes...

                    for (int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        Debug.Assert(!double.IsNaN(cpv.Parameters_IN[np]));
                        Debug.Assert(!double.IsInfinity(cpv.Parameters_IN[np]));
                    }

                    for (int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for (int na = 0; na < _NOargs; na++) {
                        if (Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if (GradUin[na] != null) {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                            }
                        }
                    }

                    for (int d = 0; d < D; d++) {
                        _GradV_in[d] = 1;
                        f[l, k, d] += this.BoundaryEdgeForm(ref cpv, _U_in, _GradU_in, _V_in, _GradV_in);
                        _GradV_in[d] = 0;
                    }
                }
            }
        }

        void INonlinInnerEdgeForm_V.NonlinInternalEdge_V(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout, MultidimensionalArray fin, MultidimensionalArray fot) {

            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fot.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per edge

            for (int edge = 0; edge < NumOfCells; edge++) { // loop over edges...
                int iEdge = efp.e0 + edge;

                int jCellIn = efp.GridDat.iGeomEdges.CellIndices[iEdge, 0];
                int jCellOut = efp.GridDat.iGeomEdges.CellIndices[iEdge, 1];
                double pnlty = penalty(efp.GridDat, jCellIn, jCellOut, iEdge);

                int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
                int _NOArgumentsForViscosity = this.ArgumentOrdering.Count - m_D; // Dont include the 2  or 3 velocity components


                double[] ParametersIN = m_useArgumentsForViscosity ? new double[_NOArgumentsForViscosity] : new double[_NOParams];
                double[] ParametersOT = m_useArgumentsForViscosity ? new double[_NOArgumentsForViscosity] : new double[_NOParams];

                for (int node = 0; node < NumOfNodes; node++) {
                    // SIPG Flux Loops
                    if (m_useArgumentsForViscosity) {
                        for (int np = 0; np < _NOArgumentsForViscosity; np++) {
                            ParametersIN[np] = Uin[m_D + np][edge, node];
                            ParametersOT[np] = Uout[m_D + np][edge, node];
                        }
                    } else {
                        for (int np = 0; np < _NOParams; np++) {
                            ParametersIN[np] = efp.ParameterVars_IN[np][edge, node];
                            ParametersOT[np] = efp.ParameterVars_OUT[np][edge, node];
                        }
                    }

                    double viscosityIN = Viscosity(ParametersIN);
                    double viscosityOT = Viscosity(ParametersOT);
                    Debug.Assert(!double.IsNaN(viscosityIN));
                    Debug.Assert(!double.IsNaN(viscosityOT));
                    Debug.Assert(!double.IsInfinity(viscosityIN));
                    Debug.Assert(!double.IsInfinity(viscosityOT));


                    if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_u) != 0) {
                        double flux = 0.0;
                        for (int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                            double n = efp.Normals[edge, node, d];
                            flux -= 0.5 * (viscosityIN * GradUin[m_iComp][edge, node, d] + viscosityOT * GradUout[m_iComp][edge, node, d]) * n * base.m_alpha;
                        }
                        flux += Math.Max(viscosityIN, viscosityOT) * (Uin[m_iComp][edge, node] - Uout[m_iComp][edge, node]) * pnlty;

                        fin[edge, node] += flux;
                        fot[edge, node] -= flux;

                    }


                    // Add SWIPViscosity2 term
                    if ((m_ViscositySwitches & ViscosityTermsSwitch.grad_uT) != 0) {
                        double flux = 0.0;
                        for (int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                            double n = efp.Normals[edge, node, d];
                            flux -= 0.5 * (viscosityIN * GradUin[d][edge, node, m_iComp] + viscosityOT * GradUout[d][edge, node, m_iComp]) * n * base.m_alpha;
                        }
                        flux += Math.Max(viscosityIN, viscosityOT) * (Uin[m_iComp][edge, node] - Uout[m_iComp][edge, node]) * pnlty;

                        fin[edge, node] += flux;
                        fot[edge, node] -= flux;

                    }
                    // Add SWIPViscosity3 term
                    if ((m_ViscositySwitches & ViscosityTermsSwitch.divU) != 0) {
                        double flux = 0.0;
                        for (int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                            double n = efp.Normals[edge, node, m_iComp];
                            flux -= 0.5 * (viscosityIN * GradUin[d][edge, node, d] + viscosityOT * GradUout[d][edge, node, d]) * n * base.m_alpha;
                        }
                        flux += Math.Max(viscosityIN, viscosityOT) * (Uin[m_iComp][edge, node] - Uout[m_iComp][edge, node]) * pnlty;

                        fin[edge, node] += flux * (-2.0 / 3.0);
                        fot[edge, node] -= flux * (-2.0 / 3.0);
                    }


                }
            }




        }

        void INonlinBoundaryEdgeForm_V.NonlinBoundaryEdge_V(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin, MultidimensionalArray f) {
            int L = efp.Len;
            Debug.Assert(f.GetLength(0) == L);
            int K = f.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;
            int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            Debug.Assert(_NOParams == efp.ParameterVars_IN.Length);
            int _NOargs = this.ArgumentOrdering.Count;
            Debug.Assert(_NOargs == Uin.Length);
            Debug.Assert(_NOargs == GradUin.Length);

            CommonParamsBnd cpv;
            cpv.GridDat = efp.GridDat;
            cpv.Parameters_IN = new double[_NOParams];
            cpv.Normal = new double[D];
            cpv.X = new double[D];
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            byte[] EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;



            for (int l = 0; l < L; l++) { // loop over edges...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for (int k = 0; k < K; k++) { // loop over nodes...

                    for (int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        Debug.Assert(!double.IsNaN(cpv.Parameters_IN[np]));
                        Debug.Assert(!double.IsInfinity(cpv.Parameters_IN[np]));
                    }

                    for (int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for (int na = 0; na < _NOargs; na++) {
                        if (Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if (GradUin[na] != null) {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for (int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                            }
                        }
                    }

                    _V_in = 1;
                    f[l, k] += this.BoundaryEdgeForm(ref cpv, _U_in, _GradU_in, _V_in, _GradV_in);
                }
            }

        }


    }

    /// <summary>
    /// A configuration switch to <see cref="swipViscosityBase"/> and decendants, to
    /// switch constant / variable viscosity
    /// and non-dimensionless / dimensionless formulation.
    /// </summary>
    public enum ViscosityOption2 {

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
    /// 
    /// </summary>
    [Flags]
    public enum ViscosityTermsSwitch {
        /// <summary>
        /// Nothing
        /// </summary>
        None = 0,

        /// <summary>
        ///  SWIP 1
        ///  </summary>
        grad_u = 0x1,

        /// <summary>
        ///  SWIP 2
        /// </summary>
        grad_uT = 0x2,

        /// <summary>
        /// SWIP 3
        /// </summary>
        divU = 0x4
    }


}

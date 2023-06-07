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

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// A configuration switch to <see cref="SipViscosityBase"/> and decendants, to
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
    public abstract class SipViscosityBase : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, IEquationComponentCoefficient, IEquationComponentChecking, ISupportsJacobianComponent, IDGdegreeConstraint {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="_penaltyBase">
        /// base/safety factor for the penalty factor; should **not** depend on geometry or DG degree, these contributions are factored in by <see cref="penalty(IGridData, int, int, int)"/>
        /// </param>
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
        protected SipViscosityBase(double _penaltyBase,
                                    int iComp, int D, IncompressibleBoundaryCondMap bcmap,
                                    ViscosityOption _ViscosityMode, double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null//, bool ignoreVectorized = false
                                                                                                                                                                    ) {
            this.m_penalty_base = _penaltyBase;
            if(m_penalty_base.IsNaNorInf() || m_penalty_base < 0)
                throw new ArgumentOutOfRangeException();
            this.m_iComp = iComp;
            this.m_D = D;
            velFunction = D.ForLoop(d => bcmap.bndFunction[VariableNames.Velocity_d(d)]);
            EdgeTag2Type = bcmap.EdgeTag2Type;
            this.m_PhysicsMode = bcmap.PhysMode;
            //this.m_ignoreVectorized = ignoreVectorized;
            this.m_ViscosityMode = _ViscosityMode;
            switch(_ViscosityMode) {
                case ViscosityOption.ConstantViscosity:
                    if(double.IsNaN(constantViscosityValue))
                        throw new ArgumentException("constantViscosityValue is missing!");
                    this.m_constantViscosityValue = constantViscosityValue;
                    break;
                case ViscosityOption.ConstantViscosityDimensionless:
                    if(double.IsNaN(reynolds))
                        throw new ArgumentException("Reynolds number is missing!");
                    this.m_reynolds = reynolds;
                    break;
                case ViscosityOption.VariableViscosity:
                    this.m_EoS = EoS;
                    break;
                case ViscosityOption.VariableViscosityDimensionless:
                    if(double.IsNaN(reynolds))
                        throw new ArgumentException("Reynolds number is missing!");
                    this.m_reynolds = reynolds;
                    this.m_EoS = EoS;
                    break;
                default:
                    throw new NotImplementedException();
            }

        }

        /// <summary>
        /// checks that we have at least DG degree 1
        /// </summary>
        public bool IsValidDomainDegreeCombination(int[] DomainDegreesPerVariable, int CodomainDegree) {
            if (DomainDegreesPerVariable.Min() < 1)
                return false;
            if (CodomainDegree < 1)
                return false;

            return true;
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
        ViscosityOption m_ViscosityMode;

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
        
        ///// <summary>
        ///// Optional switch for ignoring the vectorized versions of fluxes. May be used for testing 
        ///// </summary>
        //bool m_ignoreVectorized = false;

        /// <summary>
        /// PhysicsMode needed for ParamterOrdering.
        /// </summary>
        PhysicsMode m_PhysicsMode;

        /// <summary>
        /// Viscosity term coefficient - in dependence of the equation this could be one-over-Reynolds, molecular viscosity, etc.;
        /// Dependence on <paramref name="Velocity"/> and <paramref name="VelocityGrad"/> allows non-linear viscosity models.
        /// </summary>
        virtual protected double Viscosity(double[] Parameters, double[] Velocity, double[,] VelocityGrad) {
            switch(m_ViscosityMode) {
                case ViscosityOption.ConstantViscosity:
                return m_constantViscosityValue;
                case ViscosityOption.ConstantViscosityDimensionless:
                return (1.0 / m_reynolds);
                case ViscosityOption.VariableViscosity:
                if(m_EoS == null) {
                    return Parameters[0];
                } else {
                    if(Parameters.Length > 1) {
                        return m_EoS.GetViscosity(Parameters);
                    } else {
                        return m_EoS.GetViscosity(Parameters[0]);
                    }
                }
                case ViscosityOption.VariableViscosityDimensionless:
                if(m_EoS == null) {
                    return (Parameters[0] / m_reynolds);
                } else {
                    return (m_EoS.GetViscosity(Parameters[0]) / m_reynolds);
                }
                default:
                throw new NotImplementedException();
            }

        }

        // store instances of object used in 'Viscosity' -- theoretically, this should produce less stress on the garbage collector
        double[] m_Parameters;
        double[] m_Velocity;
        double[,] m_VelocityGrad;

        protected double Viscosity(int i, int n, int D, MultidimensionalArray[] Parameters, MultidimensionalArray[] Velocity, MultidimensionalArray[] VelocityGrad) {
            int noOfParams = Parameters != null ? Parameters.Length : 0;
            if(m_Parameters == null)
                m_Parameters = new double[noOfParams];

            //int NoOfVelocities1 = Velocity != null ? Velocity.Length : 0;
            //int NoOfVelocities2 = VelocityGrad != null ? VelocityGrad.Length : 0;
            
            if(m_Velocity == null)
                m_Velocity = new double[D];

            //int D = VelocityGrad != null ? VelocityGrad[0].GetLength(2) : 0;
            if(m_VelocityGrad == null) 
                m_VelocityGrad = new double[D, D];

            
            for(int f = 0; f < noOfParams; f++) {
                m_Parameters[f] = Parameters[f][i, n];
            }
            for(int f = 0; f < D; f++) {
                if(Velocity[f] != null)
                    m_Velocity[f] = Velocity[f][i, n];
                else
                    m_Velocity[f] = 0.0;
                if(VelocityGrad[f] != null) {
                    for(int d = 0; d < D; d++) {
                        m_VelocityGrad[f, d] = VelocityGrad[f][i, n, d];
                    }
                } else {
                    for(int d = 0; d < D; d++) {
                        m_VelocityGrad[f, d] = 0.0;
                    }
                }
            }

            return Viscosity(m_Parameters, m_Velocity, m_VelocityGrad);
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
                    throw new NotImplementedException("Nonlinear dependence - todo.");
                   
                    //return new IEquationComponent[] { this };


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
        private double m_penalty_base;


        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            m_D = cs.GrdDat.SpatialDimension;
            double _D = m_D;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;
            if(cs.UserDefinedValues != null)
            {
                if (cs.UserDefinedValues.Keys.Contains("SlipLengths"))
                    Lslip = (MultidimensionalArray)cs.UserDefinedValues["SlipLengths"];
                // Set the Reynolds number to a user defined value contained in the CoefficientSet cs
                // Useful in case that the Reynolds number changes during a simulation...
                if (cs.UserDefinedValues.Keys.Contains("Reynolds"))
                    m_reynolds = (double)cs.UserDefinedValues["Reynolds"];
            }
        }


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
            Debug.Assert(!double.IsInfinity(m_penalty_base));
            Debug.Assert(!double.IsInfinity(m_penalty));

            double µ = penaltySizeFactor * m_penalty * m_penalty_base;
            if (µ <= 0 || µ.IsNaNorInf()) {
                string errStr = ($"Inf/NaN/Non-Positive in penalty comp: {µ}; (m_penalty = {m_penalty}, m_penalty_base = {m_penalty_base}, jCellIn = {jCellIn}, jCellOut = {jCellOut}, cj_in = {cj[jCellIn]}, cj_out = {(jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0)}, penaltySizeFactor_A = {penaltySizeFactor_A}, penaltySizeFactor_B = {penaltySizeFactor_B})");
                throw new ArithmeticException(errStr);
            }
            return µ;
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
                switch(m_ViscosityMode) {
                    case ViscosityOption.ConstantViscosity:
                    case ViscosityOption.ConstantViscosityDimensionless:
                        return new string[0];
                    case ViscosityOption.VariableViscosity:
                    case ViscosityOption.VariableViscosityDimensionless:
                        if(m_EoS == null) {
                            return new string[] { VariableNames.ViscosityMolecular };
                        } else {
                            switch(m_PhysicsMode) {
                                case PhysicsMode.LowMach:
                                case PhysicsMode.Combustion:
                                    return new string[] { VariableNames.Temperature0 };
                                case PhysicsMode.Multiphase:
                                    return new string[] { VariableNames.LevelSet };
                                case PhysicsMode.RANS:
                                    return new string[] {"k0", "omega0"}; // TODO generalize this for other turbulence models
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
        /// very dirty hack to 'inject' an alternate boundary condition value for unit testing
        /// </summary>
        public Func<int, double[], double> g_Diri_Override;

        /// <summary>
        /// very dirty hack to 'inject' an alternate boundary condition value for unit testing
        /// </summary>
        public Func<int, double[], int, double> g_Neu_Override;

        /// <summary>
        /// Velocity gradient part for stress boundary condition
        /// - 1st index: spatial dimension 
        /// - 2nd index: edge tag
        /// </summary>
        public Func<double[], double, double>[][] g_Neu_GradU;


        /// <summary>
        /// Dirichlet boundary value: the given velocity at the boundary.
        /// </summary>
        protected virtual double g_Diri(double[] X, double time, int EdgeTag, int d) {
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

        virtual public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV);
            }
        }

        virtual public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        virtual public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        bool IEquationComponentChecking.IgnoreVectorizedImplementation {
            get {
                return false;
            }
        }

        abstract public double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV);

        abstract public double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB);

        abstract public double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA);


    }

    /// <summary>
    /// \f[ 
    ///   -\mathrm{div} \left( \mu \nabla \vec{u} \right)
    /// \f]
    /// </summary>
    public class SipViscosity_GradU : SipViscosityBase,
        INonlinVolumeForm_GradV,
        INonlinEdgeForm_GradV,
        INonlinEdgeForm_V  //
    {

        /// <summary>
        /// ctor; parameter documentation see <see cref="SipViscosityBase.SipViscosityBase"/>.
        /// </summary>
        public SipViscosity_GradU(double _penalty_safety, int iComp, int D, IncompressibleBoundaryCondMap bcmap,
                                   ViscosityOption _ViscosityMode, double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null)
            : base(_penalty_safety, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS) {
            
            //g_Neu_GradU = D.ForLoop(d => bcmap.bndFunction[VariableNames.Velocity_GradientVector(D).GetRow(iComp)[d]]);
        }

      
        public override double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                //acc -= GradU[0, d] * GradV[d] * Viscosity(cpv.Parameters) * base.m_alpha;
                acc -= GradU[m_iComp, d] * GradV[d];

            if(acc != 0.0)
                acc *= Viscosity(cpv.Parameters, U, GradU) * base.m_alpha;

            return -acc;
        }



        public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;
            double pnlty = this.penalty(inp.GridDat, inp.jCellIn, inp.jCellOut, inp.iEdge);
            double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
            double muB = this.Viscosity(inp.Parameters_OUT, _uB, _Grad_uB);


            for(int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (muA * _Grad_uA[m_iComp, d] + muB * _Grad_uB[m_iComp, d]) * (_vA - _vB) * inp.Normal[d];  // consistency term  
                Acc += 0.5 * (muA * _Grad_vA[d] + muB * _Grad_vB[d]) * (_uA[m_iComp] - _uB[m_iComp]) * inp.Normal[d];  // symmetry term
            }
            Acc *= base.m_alpha;

            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            //Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * muMax; // penalty term
            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax; // penalty term 
            return -Acc;
        }


        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag, double t = 0.0) {
            if(base.g_Neu_Override == null) {
                double Acc = 0.0;
                if (!g_Neu_GradU.IsNullOrEmpty()) {
                    for (int i = 0; i < base.m_D; i++) {
                        Acc += N[i] * g_Neu_GradU[i][EdgeTag](X, t);
                    }
                }
                return Acc;
            } else {
                double Acc = 0.0;
                for(int i = 0; i < base.m_D; i++) {
                    Acc += N[i] * g_Neu_Override(base.m_iComp, X, i);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;
            double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);//, inp.GridDat.Cells.cj);
            
            double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++

                    double g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, m_iComp);

                    for(int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        //Acc += (muA * _Grad_uA[0, d]) * (_vA) * nd;
                        //Acc += (muA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;
                        Acc += (muA * _Grad_uA[m_iComp, d]) * (_vA) * nd;
                        Acc += (muA * _Grad_vA[d]) * (_uA[m_iComp] - g_D) * nd;
                    }
                    Acc *= base.m_alpha;

                    //Acc -= muA * (_uA[0] - g_D) * (_vA - 0) * pnlty;
                    Acc -= muA * (_uA[m_iComp] - g_D) * (_vA - 0) * pnlty;
                    break;
                }
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry: {

                    int D = inp.D;
                    double g_D;

                    for(int dN = 0; dN < D; dN++) {

                        g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);

                        for(int dD = 0; dD < D; dD++) {
                            // consistency
                            Acc += muA * (inp.Normal[dN] * _Grad_uA[dN, dD] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                            // symmetry
                            Acc += muA * (inp.Normal[m_iComp] * _Grad_vA[dD] * inp.Normal[dD]) * (_uA[dN] - g_D) * inp.Normal[dN] * base.m_alpha;
                        }

                        // penalty
                        Acc -= muA * ((_uA[dN] - g_D) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[m_iComp]) * pnlty;
                        //Acc = 0;
                    }
                    break;
                }
                case IncompressibleBcType.NavierSlip_Linear: {

                    double ls = Lslip[inp.jCellIn];
                    if(ls == 0.0)
                        goto case IncompressibleBcType.Velocity_Inlet;

                    if(ls > 0)
                        m_beta = muA / ls;


                    int D = inp.D;
                    double g_D;

                    for(int dN = 0; dN < D; dN++) {
                        g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);

                        for(int dD = 0; dD < D; dD++) {
                            // consistency
                            Acc += muA * (inp.Normal[dN] * _Grad_uA[dN, dD] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                            // symmetry
                            Acc += muA * (inp.Normal[m_iComp] * _Grad_vA[dD] * inp.Normal[dD]) * (_uA[dN] - g_D) * inp.Normal[dN] * base.m_alpha;
                        }

                        // penalty
                        Acc -= muA * ((_uA[dN] - g_D) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[m_iComp]) * pnlty;
                    }


                    double[,] P = new double[D, D];
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            double nn = inp.Normal[d1] * inp.Normal[d2];
                            if(d1 == d2) {
                                P[d1, d2] = 1 - nn;
                            } else {
                                P[d1, d2] = -nn;
                            }
                        }
                    }

                    // tangential dissipation force term
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
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
                case IncompressibleBcType.Pressure_Outlet: {
                    // Atmospheric outlet/pressure outflow: hom. Neumann
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    double g_N = g_Neu(inp.X, inp.Normal, inp.EdgeTag, inp.time);

                    //for (int d = 0; d < inp.D; d++) {
                        //Acc += (muA * _Grad_uA[m_iComp, d]) * (_vA) * inp.Normal[d];
                        //Acc += (muA * g_N[d]) * (_vA) * inp.Normal[d];
                    //}

                    Acc += muA * g_N * _vA * base.m_alpha;

                    break;
                }
                case IncompressibleBcType.Pressure_Dirichlet: {
                    // Dirichlet boundary condition for pressure.
                    // Inner values of velocity gradient are taken, i.e.
                    // no boundary condition for the velocity (resp. velocity gradient) is imposed.                        

                    for(int d = 0; d < inp.D; d++) {
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

        
        
        void INonlinVolumeForm_GradV.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) {
            int NumofCells = prm.Len;
            int NumOfNodes = f.GetLength(1); // no of nodes per cell
            Debug.Assert(f.GetLength(0) == NumofCells);
            //int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            //double[] Parameters = new double[_NOParams];

            int D = prm.GridDat.SpatialDimension;
            for(int cell = 0; cell < NumofCells; cell++) { // loop over cells...
                for(int node = 0; node < NumOfNodes; node++) { // loop over nodes... 
                    //for(int np = 0; np < _NOParams; np++) {
                    //    Parameters[np] = prm.ParameterVars[np][cell, node];
                    //}
                    double viscosity = Viscosity(cell, node, D, prm.ParameterVars, U, GradU) * base.m_alpha;
                    Debug.Assert(!double.IsNaN(viscosity));
                    Debug.Assert(!double.IsInfinity(viscosity));

                    for(int d = 0; d < prm.GridDat.SpatialDimension; d++) {
                        f[cell, node, d] += viscosity * GradU[m_iComp][cell, node, d];
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
            //double[] ParametersIN = new double[_NOParams];
            //double[] ParametersOT = new double[_NOParams];

            //int NumOfArguments = 1;
            //Debug.Assert(NumOfArguments == Uin.Length);
            //Debug.Assert(NumOfArguments == Uout.Length);
            //Debug.Assert(NumOfArguments == GradUin.Length);
            //Debug.Assert(NumOfArguments == GradUout.Length);

            int D = efp.GridDat.SpatialDimension;

            for(int edges = 0; edges < NumOfEdges; edges++) { // loop over edges...
                for(int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                    double uJump = 0.5 * (Uin[m_iComp][edges, node] - Uout[m_iComp][edges, node]);

                    //for(int np = 0; np < _NOParams; np++) {
                    //    ParametersIN[np] = efp.ParameterVars_IN[np][edges, node];
                    //    ParametersOT[np] = efp.ParameterVars_OUT[np][edges, node];
                    //}
                    double viscosityIN = Viscosity(edges, node, D, efp.ParameterVars_IN, Uin, GradUin);
                    double viscosityOT = Viscosity(edges, node, D, efp.ParameterVars_OUT, Uout, GradUout);


                    Debug.Assert(!double.IsNaN(viscosityIN));
                    Debug.Assert(!double.IsNaN(viscosityOT));
                    Debug.Assert(!double.IsInfinity(viscosityIN));
                    Debug.Assert(!double.IsInfinity(viscosityOT));



                    double fluxIn = viscosityIN * uJump;
                    double fluxOut = viscosityOT * uJump;
                    for(int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                        double n = efp.Normals[edges, node, d];
                        fIN[edges, node, d] -= fluxIn * n;
                        fOT[edges, node, d] -= fluxOut * n;
                        fIN[edges, node, d] *= base.m_alpha;
                        fOT[edges, node, d] *= base.m_alpha;
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
            cpv.Normal = new Vector(D);;
            cpv.X = new Vector(D);
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            byte[] EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for(int l = 0; l < L; l++) { // loop over edges ...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for(int k = 0; k < K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        Debug.Assert(!double.IsNaN(cpv.Parameters_IN[np]));
                        Debug.Assert(!double.IsInfinity(cpv.Parameters_IN[np]));
                    }

                    for(int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for(int na = 0; na < _NOargs; na++) {
                        if(Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if(GradUin[na] != null) {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                            }
                        }
                    }

                    for(int d = 0; d < D; d++) {
                        _GradV_in[d] = 1;
                        f[l, k, d] += this.BoundaryEdgeForm(ref cpv, _U_in, _GradU_in, _V_in, _GradV_in);
                        _GradV_in[d] = 0;
                    }
                }
            } 
        }
        //*/
        
        void INonlinInnerEdgeForm_V.NonlinInternalEdge_V(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout, MultidimensionalArray fin, MultidimensionalArray fot) {
            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fot.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;

            for(int i = 0; i < NumOfCells; i++) { // loop over cells...
                int iEdge = efp.e0 + i;

                int jCellIn = efp.GridDat.iGeomEdges.CellIndices[iEdge, 0];
                int jCellOut = efp.GridDat.iGeomEdges.CellIndices[iEdge, 1];
                double pnlty = penalty(efp.GridDat, jCellIn, jCellOut, iEdge);

                //int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
                //double[] ParametersIN = new double[_NOParams];
                //double[] ParametersOT = new double[_NOParams];
                for(int node = 0; node < NumOfNodes; node++) {
                    // SIPG Flux Loops
                    //for(int np = 0; np < _NOParams; np++) {
                    //    ParametersIN[np] = efp.ParameterVars_IN[np][cell, node];
                    //    ParametersOT[np] = efp.ParameterVars_OUT[np][cell, node];
                    //}
                    double viscosityIN = Viscosity(i, node, D, efp.ParameterVars_IN, Uin, GradUin);
                    double viscosityOT = Viscosity(i, node, D, efp.ParameterVars_OUT, Uout, GradUout);

                    Debug.Assert(!double.IsNaN(viscosityIN));
                    Debug.Assert(!double.IsNaN(viscosityOT));
                    Debug.Assert(!double.IsInfinity(viscosityIN));
                    Debug.Assert(!double.IsInfinity(viscosityOT));

                    double flux = 0.0;
                    for(int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                        double n = efp.Normals[i, node, d];
                        flux -= 0.5 * (viscosityIN * GradUin[m_iComp][i, node, d] + viscosityOT * GradUout[m_iComp][i, node, d]) * n * base.m_alpha;
                    }
                    flux += Math.Max(viscosityIN, viscosityOT) * (Uin[m_iComp][i, node] - Uout[m_iComp][i, node]) * pnlty;
                
                    fin[i, node] += flux;
                    fot[i, node] -= flux;
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

 

            for(int l = 0; l < L; l++) { // loop over edges...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for(int k = 0; k < K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                        Debug.Assert(!double.IsNaN(cpv.Parameters_IN[np]));
                        Debug.Assert(!double.IsInfinity(cpv.Parameters_IN[np]));
                    }

                    for(int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for(int na = 0; na < _NOargs; na++) {
                        if(Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if(GradUin[na] != null) {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                            }
                        }
                    }

                    _V_in = 1;
                    f[l, k] += this.BoundaryEdgeForm(ref cpv, _U_in, _GradU_in, _V_in, _GradV_in);
                }
            }
        }
        //*/
    }


    /// <summary>
    /// \f[ 
    ///   - \mathrm{div} \left( \mu (\partial_d \vec{u})^T \right)
    /// \f]
    /// </summary>
    public class SipViscosity_GradUtransp : SipViscosityBase, INonlinVolumeForm_GradV,
        INonlinEdgeForm_GradV,
        INonlinEdgeForm_V {

        //private ViscositySolverMode ViscSolverMode;

        /// <summary>
        /// ctor; parameter documentation see <see cref="SipViscosityBase.SipViscosityBase"/>.
        /// </summary>
        public SipViscosity_GradUtransp(double _penalty, int iComp, int D, IncompressibleBoundaryCondMap bcmap,
                                   ViscosityOption _ViscosityMode, /*ViscositySolverMode ViscSolverMode = ViscositySolverMode.FullyCoupled,*/
                                   double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null)
            : base(_penalty, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS) {
            //this.ViscSolverMode = ViscSolverMode;
        }

        public override double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            double visc = Viscosity(cpv.Parameters, U, GradU);
            for(int d = 0; d < cpv.D; d++)
                // we want to:
                //    sum(  \partial_{m_iComp} u_d  * \partial_{d} v, d=0..D-1)
                acc += GradU[d, base.m_iComp] * GradV[d] * visc * base.m_alpha;
            return acc;
        }


        public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.GridDat, inp.jCellIn, inp.jCellOut, inp.iEdge);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
            double muB = this.Viscosity(inp.Parameters_OUT, _uB, _Grad_uB);


            for(int i = 0; i < inp.D; i++) {
                // consistency term
                Acc += 0.5 * (muA * _Grad_uA[i, m_iComp] + muB * _Grad_uB[i, m_iComp]) * (_vA - _vB) * inp.Normal[i];
                // symmetry term
                //switch(ViscSolverMode) {
                //    case ViscositySolverMode.FullyCoupled:
                        Acc += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[i] - _uB[i]) * inp.Normal[m_iComp];
                //        break;
                //    case ViscositySolverMode.Segregated:
                //        if(i == m_iComp)
                //            Acc += 0.5 * (muA * _Grad_vA[i] + muB * _Grad_vB[i]) * (_uA[i] - _uB[i]) * inp.Normal[m_iComp];
                //        break;
                //    default:
                //        throw new NotImplementedException();
                //}
            }
            Acc *= base.m_alpha;

            // penalty term
            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax;

            return -Acc;
        }


        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag, double t = 0.0) {
            if(base.g_Neu_Override == null) {
                //return 0.0;
               
                //double Acc = 0;
                //for (int i = 0; i < base.m_D; i++) {
                //    Acc += N[i] * g_Neu_GradU[i][EdgeTag](X, t);
                //}
                //return Acc;
                
                throw new NotSupportedException("Neumann BC. for the \\/U^T -- term is problematic!");

            } else {
                double Acc = 0;
                for(int i = 0; i < base.m_D; i++) {
                    Acc += N[i] * g_Neu_Override(i, X, base.m_iComp);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++
                     

                    for(int i = 0; i < inp.D; i++) {
                        // consistency
                        Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
                        // symmetry
                        //switch(ViscSolverMode) {
                        //    case ViscositySolverMode.FullyCoupled:
                                Acc += (muA * _Grad_vA[i]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[m_iComp];
                        //        break;
                        //    case ViscositySolverMode.Segregated:
                        //        if(i == m_iComp)
                        //            Acc += (muA * _Grad_vA[i]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[m_iComp];
                        //        break;
                        //    default:
                        //        throw new NotImplementedException();
                        //}
                    }
                    Acc *= base.m_alpha;

                    // penalty
                    Acc -= muA * (_uA[m_iComp] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp)) * (_vA - 0) * pnlty;

                    break;
                }
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry: {

                    int D = inp.D;
                    double g_D;

                    for(int dN = 0; dN < D; dN++) {
                        for(int dD = 0; dD < D; dD++) {
                            // consistency
                            Acc += muA * (inp.Normal[dN] * _Grad_uA[dD, dN] * inp.Normal[dD]) * (_vA * inp.Normal[m_iComp]) * base.m_alpha;
                            // symmetry
                            //switch(ViscSolverMode) {
                            //    case ViscositySolverMode.FullyCoupled:
                                    g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dD);
                                    Acc += muA * (inp.Normal[dN] * _Grad_vA[dN] * inp.Normal[m_iComp]) * (_uA[dD] - g_D) * inp.Normal[dD] * base.m_alpha;
                            //        break;
                            //    case ViscositySolverMode.Segregated:
                            //    default:
                            //        throw new NotImplementedException();
                            //}
                        }
                        g_D = this.g_Diri(inp.X, inp.time, inp.EdgeTag, dN);
                        // penalty
                        Acc -= muA * ((_uA[dN] - g_D) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[m_iComp]) * pnlty;
                    }

                    break;
                }
                case IncompressibleBcType.NavierSlip_Linear: {

                    double ls = Lslip[inp.jCellIn];
                    if(ls == 0.0)
                        goto case IncompressibleBcType.Velocity_Inlet;
                    else
                        goto case IncompressibleBcType.FreeSlip;


                }
                //case IncompressibleBcType.InnerValues: {

                //        for (int i = 0; i < inp.D; i++) {
                //            Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
                //        }
                //        Acc *= base.m_alpha;

                //        break;
                //    }
                case IncompressibleBcType.Pressure_Dirichlet:                    
                    // Inner values of velocity gradient are taken, i.e.
                    // no boundary condition for the velocity (resp. velocity gradient) is imposed.                        
                    for (int i = 0; i < inp.D; i++) {
                        Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
                    }                    
                    Acc *= base.m_alpha;
                    break;
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {
                    //if (!base.g_Neu_GradU.IsNullOrEmpty()) {
                    //    double g_N = g_Neu(inp.X, inp.Normal, inp.EdgeTag);
                    //    Acc += muA * g_N * _vA;
                    //} else
                    if (base.g_Neu_Override == null) {
                        // !!!!! for now B.C. is only imposed via the GradU Term explicitly !!!!!
                        // Inner values of velocity gradient are taken, i.e.
                        // no boundary condition for the velocity (resp. velocity gradient) is imposed.
                        for (int i = 0; i < inp.D; i++) {
                            Acc += (muA * _Grad_uA[i, m_iComp]) * (_vA) * inp.Normal[i];
                        }
                    } else {
                        double g_N = g_Neu(inp.X, inp.Normal, inp.EdgeTag);
                        Acc += muA * g_N * _vA;
                    }
                    Acc *= base.m_alpha;
                    break;
                }
                default:
                    throw new NotSupportedException();
            }

            return -Acc;
        }

        void INonlinVolumeForm_GradV.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) { //OK:)

            int NumofCells = prm.Len;
            int NumOfNodes = f.GetLength(1); // no of nodes per cell
            Debug.Assert(f.GetLength(0) == NumofCells);
            //int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            //double[] Parameters = new double[_NOParams];
            int D = prm.GridDat.SpatialDimension;

            for(int cell = 0; cell < NumofCells; cell++) { // loop over cells...
                for(int node = 0; node < NumOfNodes; node++) { // loop over nodes... 

                    //for(int np = 0; np < _NOParams; np++) {
                    //    Parameters[np] = prm.ParameterVars[np][cell, node];
                    //}
                    double viscosity = Viscosity(cell, node, D, prm.ParameterVars, U, GradU) * base.m_alpha;

                    for(int d = 0; d < prm.GridDat.SpatialDimension; d++) {
                        f[cell, node, d] += viscosity * GradU[d][cell, node, m_iComp];
                    }
                }
            }
        }

        void INonlinInnerEdgeForm_GradV.NonlinInternalEdge_GradV(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout, MultidimensionalArray fIN, MultidimensionalArray fOT) { // OK :)
            int NumOfEdges = efp.Len;
            Debug.Assert(fIN.GetLength(0) == NumOfEdges);
            Debug.Assert(fOT.GetLength(0) == NumOfEdges);
            int NumOfNodes = fIN.GetLength(1); // no of nodes per cell
            //int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            //double[] ParametersIN = new double[_NOParams];
            //double[] ParametersOT = new double[_NOParams];
            int D = efp.GridDat.SpatialDimension;


            for(int edges = 0; edges < NumOfEdges; edges++) { // loop over edges...
                for(int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    //for(int np = 0; np < _NOParams; np++) {
                    //    ParametersIN[np] = efp.ParameterVars_IN[np][edges, node];
                    //    ParametersOT[np] = efp.ParameterVars_OUT[np][edges, node];
                    //}
                    double viscosityIN = Viscosity(edges, node, D, efp.ParameterVars_IN, Uin, GradUin);
                    double viscosityOT = Viscosity(edges, node, D, efp.ParameterVars_OUT, Uout, GradUout);

                    double n = efp.Normals[edges, node, m_iComp];
                    //switch(ViscSolverMode) {
                    //    case ViscositySolverMode.FullyCoupled:
                            for(int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                                fIN[edges, node, d] -= 0.5 * viscosityIN * (Uin[d][edges, node] - Uout[d][edges, node]) * n * base.m_alpha;
                                fOT[edges, node, d] -= 0.5 * viscosityOT * (Uin[d][edges, node] - Uout[d][edges, node]) * n * base.m_alpha;
                            }
                    //        break;
                    //    case ViscositySolverMode.Segregated:
                    //        fIN[edges, node, m_iComp] -= 0.5 * viscosityIN * (Uin[m_iComp][edges, node] - Uout[m_iComp][edges, node]) * n * base.m_alpha;
                    //        fOT[edges, node, m_iComp] -= 0.5 * viscosityOT * (Uin[m_iComp][edges, node] - Uout[m_iComp][edges, node]) * n * base.m_alpha;
                    //        break;
                    //    default:
                    //        throw new NotImplementedException();
                    //
                    //}
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
            cpv.Normal = new double[D];
            cpv.X = new double[D];
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            byte[] EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for(int l = 0; l < L; l++) { // loop over edges ...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for(int k = 0; k < K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for(int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for(int na = 0; na < _NOargs; na++) {
                        if(Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if(GradUin[na] != null) {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                            }
                        }
                    }

                    for(int d = 0; d < D; d++) {
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
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;

            for(int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;

                int jCellIn = efp.GridDat.iGeomEdges.CellIndices[iEdge, 0];
                int jCellOut = efp.GridDat.iGeomEdges.CellIndices[iEdge, 1];
                double pnlty = penalty(efp.GridDat, jCellIn, jCellOut, iEdge);

                int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
                //double[] ParametersIN = new double[_NOParams];
                //double[] ParametersOT = new double[_NOParams];
                for(int node = 0; node < NumOfNodes; node++) {
                    // SIPG Flux Loops
                    //for(int np = 0; np < _NOParams; np++) {
                    //    ParametersIN[np] = efp.ParameterVars_IN[np][cell, node];
                    //    ParametersOT[np] = efp.ParameterVars_OUT[np][cell, node];
                    //}
                    double viscosityIN = Viscosity(cell, node, D, efp.ParameterVars_IN, Uin, GradUin);
                    double viscosityOT = Viscosity(cell, node, D, efp.ParameterVars_OUT, Uout, GradUout);


                    double flux = 0.0;
                    for(int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                        double n = efp.Normals[cell, node, d];
                        flux -= 0.5 * (viscosityIN * GradUin[d][cell, node, m_iComp] + viscosityOT * GradUout[d][cell, node, m_iComp]) * n * base.m_alpha;
                    }
                    flux += Math.Max(viscosityIN, viscosityOT) * (Uin[m_iComp][cell, node] - Uout[m_iComp][cell, node]) * pnlty;

                    fin[cell, node] += flux ;
                    fot[cell, node] -= flux ;
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



            for(int l = 0; l < L; l++) { // loop over edges...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for(int k = 0; k < K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for(int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for(int na = 0; na < _NOargs; na++) {
                        if(Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if(GradUin[na] != null) {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for(int d = 0; d < D; d++) {
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
    /// \f[ 
    ///   \frac{2}{3} \mathrm{div} \left( \mu \myMatrix{I} \mathrm{div} ( \vec{u} )  \right)
    /// \f]
    /// </summary>
    public class SipViscosity_divU : SipViscosityBase, INonlinVolumeForm_GradV,
        INonlinEdgeForm_GradV,
        INonlinEdgeForm_V {

        //private ViscositySolverMode ViscSolverMode;

        /// <summary>
        /// ctor; parameter documentation see <see cref="SipViscosityBase.SipViscosityBase"/>.
        /// </summary>
        public SipViscosity_divU(double _penalty, int iComp, int D, IncompressibleBoundaryCondMap bcmap,
                                   ViscosityOption _ViscosityMode/*, ViscositySolverMode ViscSolverMode = ViscositySolverMode.FullyCoupled*/,
                                   double constantViscosityValue = double.NaN, double reynolds = double.NaN, MaterialLaw EoS = null)
            : base(_penalty, iComp, D, bcmap, _ViscosityMode, constantViscosityValue, reynolds, EoS) {

            //this.ViscSolverMode = ViscSolverMode;
        }

        public override double VolumeForm(ref Foundation.CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double visc = Viscosity(cpv.Parameters, U, GradU);
            double acc = 0;
            for(int d = 0; d < cpv.D; d++)
                acc -= GradU[d, d] * GradV[base.m_iComp] * visc * base.m_alpha;
            return acc * (2.0 / 3.0);
        }


        public override double InnerEdgeForm(ref Foundation.CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = this.penalty(inp.GridDat, inp.jCellIn, inp.jCellOut, inp.iEdge);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
            double muB = this.Viscosity(inp.Parameters_OUT, _uB, _Grad_uB);

            for(int i = 0; i < inp.D; i++) {
                // consistency term
                Acc += 0.5 * (muA * _Grad_uA[i, i] + muB * _Grad_uB[i, i]) * (_vA - _vB) * inp.Normal[m_iComp];
                // symmetry term
                Acc += 0.5 * (muA * _Grad_vA[m_iComp] + muB * _Grad_vB[m_iComp]) * (_uA[i] - _uB[i]) * inp.Normal[i];

            }
            Acc *= base.m_alpha;

            // penalty term
            double muMax = (Math.Abs(muA) > Math.Abs(muB)) ? muA : muB;
            Acc -= (_uA[m_iComp] - _uB[m_iComp]) * (_vA - _vB) * pnlty * muMax;

            return Acc * (2.0 / 3.0);
        }



        /// <summary>
        /// Neumann boundary value;
        /// </summary>
        double g_Neu(double[] X, double[] N, int EdgeTag) {
            if(base.g_Neu_Override == null) {
                //return 0.0;

                throw new NotSupportedException("Neumann BC. for the \\/U^T -- term is problematic!");

            } else {
                double Acc = 0;
                for(int i = 0; i < base.m_D; i++) {
                    Acc += N[m_iComp] * g_Neu_Override(i, X, i);
                }
                return Acc;
            }
        }


        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;
            double pnlty = 2 * this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);//, inp.GridDat.Cells.cj);
            double muA = this.Viscosity(inp.Parameters_IN, _uA, _Grad_uA);
            IncompressibleBcType edgType = base.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case IncompressibleBcType.Velocity_Inlet:
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann: {
                    // inhom. Dirichlet b.c.
                    // +++++++++++++++++++++                      

                    for(int i = 0; i < inp.D; i++) {
                        // consistency
                        Acc += (muA * _Grad_uA[i, i]) * (_vA) * inp.Normal[m_iComp];
                        // symmetry
                        Acc += (muA * _Grad_vA[m_iComp]) * (_uA[i] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, i)) * inp.Normal[i];
                    }
                    Acc *= base.m_alpha;

                    // penalty
                    Acc -= muA * (_uA[m_iComp] - this.g_Diri(inp.X, inp.time, inp.EdgeTag, base.m_iComp)) * (_vA - 0) * pnlty;

                    break;
                }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {

                    if(base.g_Neu_Override == null) {
                        // Inner values of velocity gradient are taken, i.e.
                        // no boundary condition for the velocity (resp. velocity gradient) is imposed.
                        for(int i = 0; i < inp.D; i++) {
                            Acc += (muA * _Grad_uA[i, i]) * (_vA) * inp.Normal[m_iComp];
                        }
                    } else {
                        double g_N = g_Neu(inp.X, inp.Normal, inp.EdgeTag);
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

        void INonlinVolumeForm_GradV.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) {

            int NumofCells = prm.Len;
            int NumOfNodes = f.GetLength(1); // no of nodes per cell
            Debug.Assert(f.GetLength(0) == NumofCells);
            int D = prm.GridDat.SpatialDimension;

            //int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
            //double[] Parameters = new double[_NOParams];
            for(int cell = 0; cell < NumofCells; cell++) { // loop over cells...
                for(int node = 0; node < NumOfNodes; node++) { // loop over nodes... 
                    //for(int np = 0; np < _NOParams; np++) {
                    //    Parameters[np] = prm.ParameterVars[np][cell, node];
                    //}
                    double viscosity = Viscosity(cell, node, D, prm.ParameterVars, U, GradU);
                    for(int d = 0; d < prm.GridDat.SpatialDimension; d++) {
                        f[cell, node, m_iComp] -= viscosity * GradU[d][cell, node, d] * (2.0 / 3.0) * base.m_alpha;
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
            //double[] ParametersIN = new double[_NOParams];
            //double[] ParametersOT = new double[_NOParams];
            int D = efp.GridDat.SpatialDimension;

            for(int edges = 0; edges < NumOfEdges; edges++) { // loop over edges...
                for(int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    //for(int np = 0; np < _NOParams; np++) {
                    //    ParametersIN[np] = efp.ParameterVars_IN[np][edges, node];
                    //    ParametersOT[np] = efp.ParameterVars_OUT[np][edges, node];
                    //}

                    double viscosityIN = Viscosity(edges, node, D, efp.ParameterVars_IN, Uin, GradUin);
                    double viscosityOT = Viscosity(edges, node, D, efp.ParameterVars_OUT, Uout, GradUout);

                    double n;
                    for(int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                        n = efp.Normals[edges, node, d];
                        fIN[edges, node, m_iComp] -= 0.5 * viscosityIN * (Uin[d][edges, node] - Uout[d][edges, node]) * n * base.m_alpha * (-2.0 / 3.0);
                        fOT[edges, node, m_iComp] -= 0.5 * viscosityOT * (Uin[d][edges, node] - Uout[d][edges, node]) * n * base.m_alpha * (-2.0 / 3.0);
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
            cpv.Normal = new double[D];
            cpv.X = new double[D];
            cpv.time = efp.time;

            double[] _GradV_in = new double[D];
            double[,] _GradU_in = new double[_NOargs, D];
            double[] _U_in = new double[_NOargs];
            double _V_in = 0.0;
            byte[] EdgeTags = efp.GridDat.iGeomEdges.EdgeTags;

            for(int l = 0; l < L; l++) { // loop over edges ...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for(int k = 0; k < K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for(int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for(int na = 0; na < _NOargs; na++) {
                        if(Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if(GradUin[na] != null) {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = 0;
                            }
                        }
                    }

                    for(int d = 0; d < D; d++) {
                        _GradV_in[d] = 1;
                        f[l, k, d] += this.BoundaryEdgeForm(ref cpv, _U_in, _GradU_in, _V_in, _GradV_in);
                        _GradV_in[d] = 0;
                    }
                }
            }
        }


        void INonlinInnerEdgeForm_V.NonlinInternalEdge_V(ref EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout, MultidimensionalArray fin, MultidimensionalArray fot) {
            int NumOfEdges = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfEdges);
            Debug.Assert(fot.GetLength(0) == NumOfEdges);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell
            int D = efp.GridDat.SpatialDimension;

            for(int i = 0; i < NumOfEdges; i++) { // loop over cells...
                int iEdge = efp.e0 + i;

                int jCellIn = efp.GridDat.iGeomEdges.CellIndices[iEdge, 0];
                int jCellOut = efp.GridDat.iGeomEdges.CellIndices[iEdge, 1];
                double pnlty = penalty(efp.GridDat, jCellIn, jCellOut, iEdge);

                int _NOParams = this.ParameterOrdering == null ? 0 : this.ParameterOrdering.Count;
                //double[] ParametersIN = new double[_NOParams];
                //double[] ParametersOT = new double[_NOParams];
                for(int node = 0; node < NumOfNodes; node++) {
                    // SIPG Flux Loops
                    //for(int np = 0; np < _NOParams; np++) {
                    //    ParametersIN[np] = efp.ParameterVars_IN[np][i, node];
                    //    ParametersOT[np] = efp.ParameterVars_OUT[np][i, node];
                    //}
                    double viscosityIN = Viscosity(i, node, D, efp.ParameterVars_IN, Uin, GradUin);
                    double viscosityOT = Viscosity(i, node, D, efp.ParameterVars_OUT, Uout, GradUout);


                    double flux = 0.0;
                    for(int d = 0; d < efp.GridDat.SpatialDimension; d++) {
                        double n = efp.Normals[i, node, m_iComp];
                        flux -= 0.5 * (viscosityIN * GradUin[d][i, node, d] + viscosityOT * GradUout[d][i, node, d]) * n * base.m_alpha;
                    }
                    flux += Math.Max(viscosityIN, viscosityOT) * (Uin[m_iComp][i, node] - Uout[m_iComp][i, node]) * pnlty;

                    fin[i, node] += flux * (-2.0 / 3.0);
                    fot[i, node] -= flux * (-2.0 / 3.0);
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



            for(int l = 0; l < L; l++) { // loop over edges...
                cpv.iEdge = efp.e0 + l;
                cpv.EdgeTag = EdgeTags[cpv.iEdge];

                for(int k = 0; k < K; k++) { // loop over nodes...

                    for(int np = 0; np < _NOParams; np++) {
                        cpv.Parameters_IN[np] = efp.ParameterVars_IN[np][l, k];
                    }

                    for(int d = 0; d < D; d++) {
                        cpv.Normal[d] = efp.Normals[l, k, d];
                        cpv.X[d] = efp.Nodes[l, k, d];
                    }

                    for(int na = 0; na < _NOargs; na++) {
                        if(Uin[na] != null) {
                            _U_in[na] = Uin[na][l, k];
                        } else {
                            _U_in[na] = 0;
                        }
                        if(GradUin[na] != null) {
                            for(int d = 0; d < D; d++) {
                                _GradU_in[na, d] = GradUin[na][l, k, d];
                            }
                        } else {
                            for(int d = 0; d < D; d++) {
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
}

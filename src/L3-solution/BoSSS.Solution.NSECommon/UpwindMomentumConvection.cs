using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.NSECommon {
    
    /// <summary>
    /// Upwind-based convection operator for the momentum equation;
    /// nonlinear implementation, but supports <see cref="DifferentialOperator.GetJacobiOperator"/>.
    /// </summary>
    public class UpwindMomentumConvection : IVolumeForm, IEdgeForm, ISupportsJacobianComponent {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;


        IncompressibleBoundaryCondMap m_bcmap;

        /// <summary>
        /// Component index of the momentum equation.
        /// </summary>
        protected int m_component;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] velFunction;

        /// <summary>
        /// Fluid density
        /// </summary>
        double m_rho;

        /// <summary>
        /// Ctor for common part of incompressible and low Mach number flows.
        /// </summary>
        public UpwindMomentumConvection(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, double __rho) {
            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;
            m_component = _component;
            m_rho = __rho;

            velFunction = new Func<double[], double, double>[m_bcmap.MaxEdgeTagNo, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                velFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d)], d);
        }


        /// <summary>
        /// flux at the boundary
        /// </summary>
        double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            var _Uin = new Vector(Uin);
            Debug.Assert(inp.D == _Uin.Dim);

            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NavierSlip_Linear:
                case IncompressibleBcType.Velocity_Inlet: {

                    var _Uot = new Vector(inp.D);
                    for (int d = 0; d < inp.D; d++)
                        _Uot[d] = velFunction[inp.EdgeTag, d](inp.X, inp.time);

                    return (_Uot * inp.Normal) * _Uot[m_component] * m_rho;
                }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {


                    return (_Uin * inp.Normal) * _Uin[m_component] * m_rho;
                }
                default:
                throw new NotImplementedException("Boundary condition not implemented!");
            }
        }

        /// <summary>
        /// Nonlinear upwind flux
        /// </summary>
        double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            var _Uin = new Vector(Uin);
            var _Uot = new Vector(Uout);

            Vector Umean =  (_Uin + _Uot) * 0.5;
            if (Umean * inp.Normal > 0) {
                r = _Uin * inp.Normal * _Uin[m_component];
            } else {
                r = _Uot * inp.Normal * _Uot[m_component];
            }

            return r * m_rho;
        }

        /// <summary>
        /// returns
        /// \f[ 
        ///   \vec{v} \cdot u_d,
        /// \f]
        /// where \f$ \vec{v}\f$  is the linearization point.
        /// For variable density the result is multiplied by \f$ \rho\f$ .
        /// </summary>
        void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            var _U = new Vector(U);
            Debug.Assert(inp.D == _U.Dim);

            output.SetV(_U * _U[m_component] * m_rho);
        }

        /// <summary>
        /// 
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double[] R = new double[m_SpatialDimension];
            Flux(ref cpv, U, R);
            return -R.InnerProd(GradV);
        }

        /// <summary>
        /// 
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return InnerEdgeFlux(ref inp, _uIN, _uOUT) * (_vIN - _vOUT);
        }

        /// <summary>
        ///
        /// </summary>
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return BorderEdgeFlux(ref inp, _uA) * _vA;
        }

        /// <summary>
        /// 
        /// </summary>
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var ConvDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            //var ConvDerivEdg = new OldEdgeFormDifferentiator(this);
            var ConvDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { ConvDerivEdg, ConvDerivVol };
        }

        /// <summary>
        /// the velocity vector
        /// </summary>
        public IList<string> ArgumentOrdering => VariableNames.VelocityVector(m_SpatialDimension);

        /// <summary>
        /// No Parameters
        /// </summary>
        public IList<string> ParameterOrdering => null;

        /// <summary>
        /// 
        /// </summary>
        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        /// <summary>
        /// 
        /// </summary>
        public bool IgnoreVectorizedImplementation => false;

        /// <summary>
        /// 
        /// </summary>
        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        /// <summary>
        /// 
        /// </summary>
        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;
    }




    /// <summary>
    /// Upwind-based convection operator for the momentum equation;
    /// nonlinear implementation, but supports <see cref="DifferentialOperator.GetJacobiOperator"/>.
    /// </summary>
    public class LocalLaxFriedrichsConvection : IVolumeForm, IEdgeForm, ISupportsJacobianComponent, ISpeciesFilter {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        protected int m_SpatialDimension;

        /// <summary>
        /// 
        /// </summary>
        IncompressibleBoundaryCondMap m_bcmap;

        /// <summary>
        /// Component index of the momentum equation.
        /// </summary>
        protected int m_component;

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] velFunction;

        /// <summary>
        /// Fluid density
        /// </summary>
        double m_rho;
        
        /// <summary>
        /// LLF value for convective operators without density,
        /// i.e. incompressible momentum equation and
        /// Level-Set advection for multiphase flows.
        /// </summary>
        static public double GetLambda(Vector VelocityMean, Vector Normal) {
            Debug.Assert(VelocityMean.Dim == Normal.Dim, "Mismatch in dimensions!");
            double V_n = VelocityMean * Normal;
            double Lambda = V_n*2;
            return Math.Abs(Lambda);
        }

        /// <summary>
        /// Ctor for common part of incompressible and low Mach number flows.
        /// </summary>
        public LocalLaxFriedrichsConvection(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, double __rho, string Species) {
            if(SpatDim < 2 || SpatDim > 3)
                throw new ArgumentException("unknown spatial dimension");
            if(_component < 0 || _component >= SpatDim)
                throw new ArgumentException("component index out of range");
            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;
            m_component = _component;
            m_rho = __rho;
            this.ValidSpecies = Species;

            velFunction = new Func<double[], double, double>[m_bcmap.MaxEdgeTagNo, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                velFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d)], d);
        }


        /// <summary>
        /// flux at the boundary
        /// </summary>
        double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            var _Uin = new Vector(Uin);
            Debug.Assert(inp.D == _Uin.Dim);

            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NavierSlip_Linear:
                case IncompressibleBcType.Velocity_Inlet: {

                    var _Uot = new Vector(inp.D);
                    for (int d = 0; d < inp.D; d++)
                        _Uot[d] = velFunction[inp.EdgeTag, d](inp.X, inp.time);

                    return (_Uot * inp.Normal) * _Uot[m_component] * m_rho;
                }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {


                    return (_Uin * inp.Normal) * _Uin[m_component] * m_rho;
                }
                default:
                throw new NotImplementedException("Boundary condition not implemented!");
            }
        }

        /// <summary>
        /// Nonlinear Local Lax-Friedrichs flux
        /// </summary>
        double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            var _Uin = new Vector(Uin);
            var _Uot = new Vector(Uout);

            Vector Umean =  (_Uin + _Uot) * 0.5;
            r += (Umean * inp.Normal) * Umean[m_component];

            double LambdaIn = GetLambda(_Uin, inp.Normal);
            double LambdaOt = GetLambda(_Uot, inp.Normal);
            double Lambda = Math.Max(LambdaIn, LambdaOt);
            double uJump = Uin[m_component] - Uout[m_component];
            r += Lambda * uJump*0.5;

            return r * m_rho;
        }

        /// <summary>
        /// returns
        /// \f[ 
        ///   \vec{v} \cdot u_d,
        /// \f]
        /// where \f$ \vec{v}\f$  is the linearization point.
        /// For variable density the result is multiplied by \f$ \rho\f$ .
        /// </summary>
        void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            var _U = new Vector(U);
            Debug.Assert(inp.D == _U.Dim);

            output.SetV(_U * _U[m_component] * m_rho);
        }

        /// <summary>
        /// 
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double[] R = new double[m_SpatialDimension];
            Flux(ref cpv, U, R);
            return -R.InnerProd(GradV);
        }

        /// <summary>
        /// 
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return InnerEdgeFlux(ref inp, _uIN, _uOUT) * (_vIN - _vOUT);
        }

        /// <summary>
        ///
        /// </summary>
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return BorderEdgeFlux(ref inp, _uA) * _vA;
        }

        /// <summary>
        /// 
        /// </summary>
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            var ConvDerivEdg = new EdgeFormDifferentiator(this, SpatialDimension);
            var ConvDerivVol = new VolumeFormDifferentiator(this, SpatialDimension);
            return new IEquationComponent[] { ConvDerivEdg, ConvDerivVol };
        }

        /// <summary>
        /// the velocity vector
        /// </summary>
        public IList<string> ArgumentOrdering => VariableNames.VelocityVector(m_SpatialDimension);

        /// <summary>
        /// No Parameters
        /// </summary>
        public IList<string> ParameterOrdering => null;

        /// <summary>
        /// 
        /// </summary>
        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        /// <summary>
        /// 
        /// </summary>
        public bool IgnoreVectorizedImplementation => false;

        /// <summary>
        /// 
        /// </summary>
        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        /// <summary>
        /// 
        /// </summary>
        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        /// <summary>
        /// 
        /// </summary>
        public string ValidSpecies {
            get;
            private set;
        }
    }
}

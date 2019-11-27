using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IBM_Solver {
 
    /// <summary>
    /// Upwind-based convection operator
    /// </summary>
    public class UpwindConvection : IVolumeForm, IEdgeForm, ISupportsJacobianComponent {

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
        public UpwindConvection(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component, double __rho) {
            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;
            m_component = _component;
            m_rho = __rho;

            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
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

            Vector Umean = (_Uin + _Uot) * 0.5;
            if(Umean*inp.Normal >= 0) {
                r = _Uin * inp.Normal * _Uin[m_component];
            } else {
                r = _Uot * inp.Normal * _Uot[m_component];
            }

            return r*m_rho;
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

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double[] R = new double[m_SpatialDimension];
            Flux(ref cpv, U, R);
            return -R.InnerProd(GradV);
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return InnerEdgeFlux(ref inp, _uIN, _uOUT) * (_vIN - _vOUT);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return BorderEdgeFlux(ref inp, _uA) * _vA;
        }

        public IEquationComponent[] GetJacobianComponents() {
            var ConvDerivEdg = new EdgeFormDifferentiator(this);
            var ConvDerivVol = new VolumeFormDifferentiator(this);
            return new IEquationComponent[] { ConvDerivEdg, ConvDerivVol };
        }

        public IList<string> ArgumentOrdering => VariableNames.VelocityVector(m_SpatialDimension);

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public bool IgnoreVectorizedImplementation => false;

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;
    }


    public class ConvectionAtIB : ILevelSetForm, ISupportsJacobianComponent {
        public ConvectionAtIB(LevelSetTracker LsTrk, int _d, int _D, double fluidDensity, bool UseMovingMesh) {
            m_LsTrk = LsTrk;
            m_D = _D;
            m_d = _d;
            fDensity = fluidDensity;
            m_UseMovingMesh = UseMovingMesh;
        }

        int m_D;
        int m_d;
        double fDensity;
        bool m_UseMovingMesh;
        LevelSetTracker m_LsTrk;

        // Use Fluxes as in Bulk Convection
        LinearizedConvection NegFlux;


        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
           
        public double LevelSetForm(ref CommonParams cp, 
            double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            Debug.Assert(m_D == cp.D);

            if (m_UseMovingMesh) {
                return 0.0;
            } else {
                Vector Uin = new Vector(U_Neg);
                return (Uin * cp.Normal) * Uin[m_d] * fDensity * v_Neg;
            }

        }

        public IEquationComponent[] GetJacobianComponents() {
            return new IEquationComponent[] { new NotImplemntedClass(this) };
        }

        class NotImplemntedClass : ILevelSetForm {
            public NotImplemntedClass(ConvectionAtIB __owner) {
                m_owner = __owner;
            }
            ConvectionAtIB m_owner;

            public int LevelSetIndex => m_owner.LevelSetIndex;

            public SpeciesId PositiveSpecies => m_owner.PositiveSpecies;

            public SpeciesId NegativeSpecies => m_owner.NegativeSpecies;

            public TermActivationFlags LevelSetTerms => m_owner.LevelSetTerms;

            public IList<string> ArgumentOrdering => m_owner.ArgumentOrdering;

            public IList<string> ParameterOrdering => m_owner.ParameterOrdering;

            public double LevelSetForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
                throw new NotImplementedException();
            }
        }


    }

}

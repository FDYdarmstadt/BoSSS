using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
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

namespace IBM_Solver {
    class VolumeFormDifferentiator : IVolumeForm {

        public VolumeFormDifferentiator(IVolumeForm vf) {
            m_VolForm = vf;
            m_eps = Math.Sqrt(BLAS.MachineEps);
            m_ParameterOrdering = ArrayTools.Cat(vf.ArgumentOrdering.Select(fn => fn + "_lin"), vf.ParameterOrdering);
        }

        double m_eps;

        IVolumeForm m_VolForm;

        public TermActivationFlags VolTerms => m_VolForm.VolTerms;

        public IList<string> ArgumentOrdering => m_VolForm.ArgumentOrdering;
        
        string[] m_ParameterOrdering;

        public IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }

        public bool IgnoreVectorizedImplementation => false;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double ret = 0.0;
            int GAMMA = m_VolForm.ArgumentOrdering.Count;


            double eps = m_eps;
            double[] Utmp = new double[GAMMA];
            Utmp.SetSubVector(cpv.Parameters, 0, GAMMA);
            
            int NoOfParams = m_VolForm.ParameterOrdering != null ? m_VolForm.ParameterOrdering.Count : 0;
            double[] OrgParams = new double[NoOfParams];
            OrgParams.SetSubVector(cpv.Parameters, GAMMA, NoOfParams);
            CommonParamsVol clonedParams = cpv;
            Debug.Assert(object.ReferenceEquals(cpv, clonedParams) == false);
            clonedParams.Parameters = OrgParams;


            //double[] dU = new double[GAMMA];
            double f0 = m_VolForm.VolumeForm(ref clonedParams, Utmp, null, V, GradV);

            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (   ((m_VolForm.VolTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0) 
                    && (U[iVar] != 0.0)) { // perf. opt.
                    // add perturbation
                    double bkup = Utmp[iVar];
                    double delta = Math.Abs(Utmp[iVar]) * eps;
                    Utmp[iVar] += delta;

                    // flux eval
                    double f1 = m_VolForm.VolumeForm(ref clonedParams, Utmp, null, V, GradV);

                    // restore un-perturbed state
                    Utmp[iVar] = bkup;

                    // compute finite difference
                    double dU_iVar = (f1 - f0) / delta;

                    // inner product
                    ret += dU_iVar * U[iVar];
                }

                if (   ((m_VolForm.VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxGradV)) != 0) 
                    && (GradU.GetRow(iVar).L2NormPow2() != 0.0)) {

                    throw new NotImplementedException("todo");
                }
            }
            //ret += GenericBlas.InnerProd(dU, U);




            return ret;
        }
    }

    class EdgeFormDifferentiator : IEdgeForm {

        public EdgeFormDifferentiator(IEdgeForm ef) {
            m_EdgForm = ef;
            m_eps = Math.Sqrt(BLAS.MachineEps);
            m_ParameterOrdering = ArrayTools.Cat(ef.ArgumentOrdering.Select(fn => fn + "_lin"), ef.ParameterOrdering);
        }

        double m_eps;

        IEdgeForm m_EdgForm;

        public TermActivationFlags InnerEdgeTerms => m_EdgForm.InnerEdgeTerms;

        public TermActivationFlags BoundaryEdgeTerms => m_EdgForm.BoundaryEdgeTerms;

        public IList<string> ArgumentOrdering => m_EdgForm.ArgumentOrdering;
        
        string[] m_ParameterOrdering;

        public IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }

        public bool IgnoreVectorizedImplementation => false;

        public double InnerEdgeForm(ref CommonParams inp, double[] U_IN, double[] U_OT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double ret = 0.0;
            int GAMMA = m_EdgForm.ArgumentOrdering.Count;

            double eps = m_eps;
            
            int NoOfParams = m_EdgForm.ParameterOrdering != null ? m_EdgForm.ParameterOrdering.Count : 0;
            double[] OrgParams_IN = new double[NoOfParams];
            double[] OrgParams_OT = new double[NoOfParams];
            OrgParams_IN.SetSubVector(inp.Parameters_IN, GAMMA, NoOfParams);
            OrgParams_OT.SetSubVector(inp.Parameters_OUT, GAMMA, NoOfParams);

            CommonParams clonedParams = inp;
            Debug.Assert(object.ReferenceEquals(inp, clonedParams) == false);
            clonedParams.Parameters_IN = OrgParams_IN;
            clonedParams.Parameters_OUT = OrgParams_OT;

            double[] U_IN_temp = new double[GAMMA];
            double[] U_OT_temp = new double[GAMMA];
            U_IN_temp.SetSubVector(inp.Parameters_IN, 0, GAMMA);
            U_OT_temp.SetSubVector(inp.Parameters_OUT, 0, GAMMA);


            //double[] dU = new double[GAMMA];
            double f0 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, null, null, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_EdgForm.InnerEdgeTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)) {
                    if (U_IN[iVar] != 0.0) {
                        // add perturbation
                        double bkup = U_IN_temp[iVar];
                        double delta = Math.Abs(U_IN_temp[iVar]) * eps;
                        U_IN_temp[iVar] += delta;

                        // flux eval
                        double f1 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, null, null, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);

                        // restore un-perturbed state
                        U_IN_temp[iVar] = bkup;

                        // compute finite difference
                        double dU_iVar = (f1 - f0) / delta;

                        // inner product
                        ret += dU_iVar * U_IN[iVar];
                    }
                    
                    if(U_OT[iVar] != 0.0) {
                        // add perturbation
                        double bkup = U_OT_temp[iVar];
                        double delta = Math.Abs(U_OT_temp[iVar]) * eps;
                        U_OT_temp[iVar] += delta;

                        // flux eval
                        double f1 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, null, null, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);

                        // restore un-perturbed state
                        U_OT_temp[iVar] = bkup;

                        // compute finite difference
                        double dU_iVar = (f1 - f0) / delta;

                        // inner product
                        ret += dU_iVar * U_OT[iVar];
                    }
                }

            }



            return ret;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotImplementedException();
        }
    }

    /*
    /// <summary>
    /// Upwind-based convection operator
    /// </summary>
    public class UpwindConvection : IVolumeForm, IEdgeForm {

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
        /// Ctor for common part of incompressible and low Mach number flows.
        /// </summary>
        /// <param name="SpatDim"></param>
        /// <param name="_bcmap"></param>
        /// <param name="_component"></param>
        public UpwindConvection(int SpatDim, IncompressibleBoundaryCondMap _bcmap, int _component) {
            m_SpatialDimension = SpatDim;
            m_bcmap = _bcmap;
            m_component = _component;

            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < SpatDim; d++)
                velFunction.SetColumn(m_bcmap.bndFunction[VariableNames.Velocity_d(d)], d);
        }




        /// <summary>
        /// flux at the boundary
        /// </summary>
        double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {

            IncompressibleBcType edgeType = m_bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.NoSlipNeumann:
                case IncompressibleBcType.FreeSlip:
                case IncompressibleBcType.SlipSymmetry:
                case IncompressibleBcType.NavierSlip_Linear:
                case IncompressibleBcType.Velocity_Inlet: {

                


                    // Setup params
                    // ============
                    Foundation.CommonParams inp2;
                    inp2.GridDat = inp.GridDat;
                    inp2.Normale = inp.Normale;
                    inp2.iEdge = inp.iEdge;
                    inp2.Parameters_IN = inp.Parameters_IN;
                    inp2.X = inp.X;
                    inp2.time = inp.time;

                    // Dirichlet value for velocity
                    // ============================
                    double Uout = velFunction[inp.EdgeTag, m_component](inp.X, inp.time);

                    // Specify Parameters_OUT
                    // ======================
                    inp2.Parameters_OUT = new double[inp.Parameters_IN.Length];

                    // Outer values for Velocity and VelocityMean
                    for (int j = 0; j < m_SpatialDimension; j++) {

                        inp2.Parameters_OUT[j] = velFunction[inp.EdgeTag, j](inp.X, inp.time);

                        // Velocity0MeanVectorOut is set to zero, i.e. always LambdaIn is used.
                        inp2.Parameters_OUT[m_SpatialDimension + j] = 0.0;
                    }

                     r = InnerEdgeFlux(ref inp2, Uin, new double[] { Uout });


                    return r;
                }
                case IncompressibleBcType.Pressure_Dirichlet:
                case IncompressibleBcType.Outflow:
                case IncompressibleBcType.Pressure_Outlet: {
                    double r = 0.0;
                    

                    

                    return r;
                }
                default:
                throw new NotImplementedException("Boundary condition not implemented!");
            }
        }

        /// <summary>
        /// bla bla bla
        /// </summary>
        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            double r = 0.0;

            // Calculate central part
            // ======================

            double rhoIn = 1.0;
            double rhoOut = 1.0;


            // 2 * {u_i * u_j} * n_j,
            // resp. 2 * {rho * u_i * u_j} * n_j for variable density
            r += rhoIn * Uin[0] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
            r += rhoOut * Uout[0] * (inp.Parameters_OUT[0] * inp.Normale[0] + inp.Parameters_OUT[1] * inp.Normale[1]);
            if (m_SpatialDimension == 3) {
                r += rhoIn * Uin[0] * inp.Parameters_IN[2] * inp.Normale[2] + rhoOut * Uout[0] * inp.Parameters_OUT[2] * inp.Normale[2];
            }

            // Calculate dissipative part
            // ==========================

            double[] VelocityMeanIn = new double[m_SpatialDimension];
            double[] VelocityMeanOut = new double[m_SpatialDimension];
            for (int d = 0; d < m_SpatialDimension; d++) {
                VelocityMeanIn[d] = inp.Parameters_IN[m_SpatialDimension + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[m_SpatialDimension + d];
            }


            return r;
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
            output[0] = U[0] * inp.Parameters[0];
            output[1] = U[0] * inp.Parameters[1];
            if (m_SpatialDimension == 3) {
                output[2] = U[0] * inp.Parameters[2];
            }

            if (m_bcmap.PhysMode == PhysicsMode.LowMach || m_bcmap.PhysMode == PhysicsMode.Multiphase) {

                double rho = EoS.GetDensity(inp.Parameters[2 * m_SpatialDimension]);
                for (int d = 0; d < m_SpatialDimension; d++)
                    output[d] *= rho;
            }

            if (m_bcmap.PhysMode == PhysicsMode.Combustion) {
                double[] args = new double[NumberOfReactants + 1];
                for (int n = 0; n < NumberOfReactants + 1; n++) {
                    args[n] = inp.Parameters[2 * m_SpatialDimension + n];
                }
                double rho = EoS.GetDensity(args);
                for (int d = 0; d < m_SpatialDimension; d++)
                    output[d] *= rho;
            }


        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double[] R = new double[m_SpatialDimension];
            Flux(ref cpv, U, R);
            return R.InnerProd(GradV);
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return InnerEdgeFlux(ref inp, _uIN, _uOUT) * (_vIN - _vOUT);
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return BorderEdgeFlux(ref inp, _uA) * _vA;
        }

        public IList<string> ArgumentOrdering => VariableNames.VelocityVector(m_SpatialDimension);

        public IList<string> ParameterOrdering => null;

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public bool IgnoreVectorizedImplementation => false;

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;
    }

*/
}

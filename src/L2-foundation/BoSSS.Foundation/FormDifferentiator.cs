using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation {

    abstract public class FormDifferentiatorCommon : IEquationComponent, IEquationComponentChecking {


        /// <summary>
        /// %
        /// </summary>
        protected FormDifferentiatorCommon(IEquationComponent vf, TermActivationFlags __Terms, int SpatialDimension) {
            m_OrgForm = vf;
            m_eps = Math.Sqrt(BLAS.MachineEps);
            m_SpatialDimension = SpatialDimension;

            Terms = __Terms;
            if ((Terms & (TermActivationFlags.GradUxV | TermActivationFlags.UxV)) != 0)
                Terms |= TermActivationFlags.V;
            if ((Terms & (TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV)) != 0)
                Terms |= TermActivationFlags.GradV;

            // Parameters for form derivative
            // ==============================

            ParamUreq = (Terms & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0;
            ParamGradUreq = (Terms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0;

            m_ParameterOrdering = new string[0];

            // original parameters
            OffsetOrgParams = m_ParameterOrdering.Length;
            if (vf.ParameterOrdering != null && vf.ParameterOrdering.Count > 0) {
                m_ParameterOrdering = m_ParameterOrdering.Cat(vf.ParameterOrdering);
            }

            // U in linearization point
            OffsetUparams = m_ParameterOrdering.Length;
            if (ParamUreq && vf.ArgumentOrdering != null && vf.ArgumentOrdering.Count > 0)
                m_ParameterOrdering = m_ParameterOrdering.Cat(vf.ArgumentOrdering.Select(fn => fn + "_lin"));

            // derivatives in linearization point
            OffsetGradUparams = m_ParameterOrdering.Length;
            if (ParamGradUreq && vf.ArgumentOrdering != null && vf.ArgumentOrdering.Count > 0) {
                foreach (string arg in vf.ArgumentOrdering) {
                    m_ParameterOrdering = m_ParameterOrdering.Cat(SpatialDimension.ForLoop(d => arg + "_lin_d[" + d + "]"));
                }
            }
        }

        /// <summary>
        /// required terms for the3 linearization 
        /// </summary>
        protected TermActivationFlags Terms;

        /// <summary>
        /// 
        /// </summary>
        protected int m_SpatialDimension;

        /// <summary>
        /// Differentiation w.r.t. Trial Variable required.
        /// </summary>
        protected bool ParamUreq;

        /// <summary>
        /// Differentiation w.r.t. Trial Variable Gradient required.
        /// </summary>
        protected bool ParamGradUreq;

        /// <summary>
        /// Offset into parameters to access Trial variable (called 'U') value at linearization point 
        /// </summary>
        protected int OffsetUparams;


        /// <summary>
        /// Offset into parameters to access Trial variable gradient (called 'U') value at linearization point 
        /// </summary>
        protected int OffsetGradUparams;

        /// <summary>
        /// Offset into parameters to access parameters of original form
        /// </summary>
        protected int OffsetOrgParams;
        
        /// <summary>
        /// Relative finite difference length
        /// </summary>
        protected double m_eps;

        IEquationComponent m_OrgForm;

       
        /// <summary>
        /// %
        /// </summary>
        public IList<string> ArgumentOrdering => m_OrgForm.ArgumentOrdering;

        string[] m_ParameterOrdering;

        /// <summary>
        /// %
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public bool IgnoreVectorizedImplementation => false;

        /// <summary>
        /// Extract parameters for original form.
        /// </summary>
        protected void GetOrgParams(double[] Parameters, out double[] OrgParams) {
            int NoOfParams = m_OrgForm.ParameterOrdering != null ? m_OrgForm.ParameterOrdering.Count : 0;
            OrgParams = new double[NoOfParams];
            OrgParams.SetSubVector(Parameters, OffsetOrgParams, NoOfParams);
        }

        /// <summary>
        /// 
        /// </summary>
        protected double GetTmpTrialVals(double[] Parameters, out double[] Utmp, out double[,] GradUtmp) {
            int GAMMA = m_OrgForm.ArgumentOrdering.Count;
            int D = m_SpatialDimension;
            double eps = m_eps;

            double delta = 0;
            Utmp = new double[GAMMA];
            GradUtmp = new double[GAMMA, D];
            for (int iVar = 0; iVar < GAMMA; iVar++) {

                Utmp[iVar] = Parameters[iVar + OffsetUparams];
                delta += Utmp[iVar].Pow2();
                for (int d = 0; d < D; d++) {
                    GradUtmp[iVar, d] = Parameters[iVar * D + d + OffsetGradUparams];
                    delta += GradUtmp[iVar, d].Pow2();
                }
            }

            delta = Math.Min(eps, Math.Sqrt(delta) * eps);
            Debug.Assert(delta > 0);
            return delta;
        }
    }



    /// <summary>
    /// Differentiation of a volume form, used e.g.to obtain a Jacobian of an operator, see <see cref="SpatialOperator.GetJacobiOperator"/>.
    /// </summary>
    public class VolumeFormDifferentiator : FormDifferentiatorCommon, IVolumeForm {

        IVolumeForm m_VolForm;

        /// <summary>
        /// %
        /// </summary>
        public VolumeFormDifferentiator(IVolumeForm vf, int SpatialDimension) : base(vf, vf.VolTerms, SpatialDimension) {
            m_VolForm = vf;
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return base.Terms;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double ret = 0.0;
            int GAMMA = m_VolForm.ArgumentOrdering.Count;
            int D = cpv.D;
            Debug.Assert(D == m_SpatialDimension, "Spatial Dimension Mismatch.");

            double eps = m_eps;

            //Utmp.SetSubVector(cpv.Parameters, 0, GAMMA);

            double delta = GetTmpTrialVals(cpv.Parameters, out var Utmp, out var GradUtmp);
            
            CommonParamsVol clonedParams = cpv;
            Debug.Assert(object.ReferenceEquals(cpv, clonedParams) == false);
            GetOrgParams(cpv.Parameters, out clonedParams.Parameters);

            //double[] dU = new double[GAMMA];
            double f0 = m_VolForm.VolumeForm(ref clonedParams, Utmp, GradUtmp, V, GradV);
            ret += f0; // affine contribution - contains V and GradV contribution

            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_VolForm.VolTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)
                    ) {
                    //&& (U[iVar] != 0.0)) { // perf. opt.
                    // add perturbation
                    double bkup = Utmp[iVar];
                    Debug.Assert(delta > 0);
                    Utmp[iVar] += delta;

                    // flux eval
                    double f1 = m_VolForm.VolumeForm(ref clonedParams, Utmp, GradUtmp, V, GradV);
                    if (double.IsInfinity(f1))
                        throw new ArithmeticException();
                    if (double.IsNaN(f1))
                        throw new ArithmeticException();

                    // restore un-perturbed state
                    Utmp[iVar] = bkup;

                    // compute finite difference
                    double dU_iVar = (f1 - f0) / delta;
                    if (double.IsInfinity(dU_iVar))
                        throw new ArithmeticException();
                    if (double.IsNaN(dU_iVar))
                        throw new ArithmeticException();

                    // inner product
                    ret += dU_iVar * U[iVar];
                    ret -= dU_iVar * Utmp[iVar]; // subtract affine contribution
                }

                if (((m_VolForm.VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxGradV)) != 0)
                    //&& (GradU.GetRow(iVar).L2NormPow2() != 0.0)) {
                    ) {

                    for (int d = 0; d < D; d++) {

                        double bkup = GradUtmp[iVar, d];
                        GradUtmp[iVar, d] += delta;

                        // flux eval
                        double f1 = m_VolForm.VolumeForm(ref clonedParams, Utmp, GradUtmp, V, GradV);
                        if (double.IsInfinity(f1))
                            throw new ArithmeticException();
                        if (double.IsNaN(f1))
                            throw new ArithmeticException();

                        // restore un-perturbed state
                        GradUtmp[iVar, d] = bkup;

                        // compute finite difference
                        double dU_iVar = (f1 - f0) / delta;
                        if (double.IsInfinity(dU_iVar))
                            throw new ArithmeticException();
                        if (double.IsNaN(dU_iVar))
                            throw new ArithmeticException();

                        // inner product
                        ret += dU_iVar * GradU[iVar, d];
                        ret -= dU_iVar * GradUtmp[iVar, d]; // subtract affine contribution
                    }
                }
            }

            return ret;
        }
    }


    /// <summary>
    /// Differentiation of a edge form, used e.g.to obtain a Jacobian of an operator, see <see cref="SpatialOperator.GetJacobiOperator"/>.
    /// </summary>
    public class EdgeFormDifferentiator : FormDifferentiatorCommon, IEdgeForm {

        IEdgeForm m_EdgForm;

        /// <summary>
        /// ctor
        /// </summary>
        public EdgeFormDifferentiator(IEdgeForm ef, int SpatialDimension) : 
            base(ef, ef.InnerEdgeTerms | ef.BoundaryEdgeTerms, SpatialDimension) //
        {
            m_EdgForm = ef;
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags InnerEdgeTerms {
            get {
                return base.Terms;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return base.Terms;
            }
        }
                

        public double InnerEdgeForm(ref CommonParams inp, double[] U_IN, double[] U_OT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double ret = 0.0;
            int GAMMA = m_EdgForm.ArgumentOrdering.Count;


            CommonParams clonedParams = inp;
            Debug.Assert(object.ReferenceEquals(inp, clonedParams) == false);
            GetOrgParams(inp.Parameters_IN, out clonedParams.Parameters_IN);
            GetOrgParams(inp.Parameters_OUT, out clonedParams.Parameters_OUT);


            double deltaIn = GetTmpTrialVals(inp.Parameters_IN, out var U_IN_temp, out var GradU_IN_temp);
            double deltaOt = GetTmpTrialVals(inp.Parameters_IN, out var U_OT_temp, out var GradU_OT_temp);
            double delta = Math.Max(deltaIn, deltaOt);


            //double[] dU = new double[GAMMA];
            double f0 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, null, null, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
            ret += f0; // affine contribution - contains V and GradV contribution
            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_EdgForm.InnerEdgeTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)) {
                    //if (U_IN[iVar] != 0.0) {
                    {
                        // add perturbation
                        double bkup = U_IN_temp[iVar];
                        U_IN_temp[iVar] += delta;

                        // flux eval
                        double f1 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, null, null, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
                        if (double.IsInfinity(f1))
                            throw new ArithmeticException();
                        if (double.IsNaN(f1))
                            throw new ArithmeticException();

                        // restore un-perturbed state
                        U_IN_temp[iVar] = bkup;

                        // compute finite difference
                        double dU_iVar = (f1 - f0) / delta;
                        if (double.IsInfinity(dU_iVar))
                            throw new ArithmeticException();
                        if (double.IsNaN(dU_iVar))
                            throw new ArithmeticException();

                        // inner product
                        ret += dU_iVar * U_IN[iVar];
                        ret -= dU_iVar * U_IN_temp[iVar]; // subtract affine contribution
                    }

                    //if(U_OT[iVar] != 0.0) {
                    {
                        // add perturbation
                        double bkup = U_OT_temp[iVar];
                        U_OT_temp[iVar] += delta;

                        // flux eval
                        double f1 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, null, null, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
                        if (double.IsInfinity(f1))
                            throw new ArithmeticException();
                        if (double.IsNaN(f1))
                            throw new ArithmeticException();

                        // restore un-perturbed state
                        U_OT_temp[iVar] = bkup;

                        // compute finite difference
                        double dU_iVar = (f1 - f0) / delta;
                        if (double.IsInfinity(dU_iVar))
                            throw new ArithmeticException();
                        if (double.IsNaN(dU_iVar))
                            throw new ArithmeticException();

                        // inner product
                        ret += dU_iVar * U_OT[iVar];
                        ret -= dU_iVar * U_OT_temp[iVar]; // subtract affine contribution
                    }
                }

            }

            return ret;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] U_IN, double[,] _Grad_uA, double _vIN, double[] _Grad_vIN) {
            double ret = 0.0;
            int GAMMA = m_EdgForm.ArgumentOrdering.Count;

            double eps = m_eps;

            int NoOfParams = m_EdgForm.ParameterOrdering != null ? m_EdgForm.ParameterOrdering.Count : 0;
            double[] OrgParams_IN = new double[NoOfParams];
            OrgParams_IN.SetSubVector(inp.Parameters_IN, GAMMA, NoOfParams);

            CommonParamsBnd clonedParams = inp;
            Debug.Assert(object.ReferenceEquals(inp, clonedParams) == false);
            clonedParams.Parameters_IN = OrgParams_IN;

            double[] U_IN_temp = new double[GAMMA];
            U_IN_temp.SetSubVector(inp.Parameters_IN, 0, GAMMA);

            //double[] dU = new double[GAMMA];
            double f0 = m_EdgForm.BoundaryEdgeForm(ref clonedParams, U_IN_temp, null, _vIN, _Grad_vIN);
            ret += f0; // affine contribution - contains V and GradV contribution
            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_EdgForm.InnerEdgeTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)) {
                    //if (U_IN[iVar] != 0.0) {
                    {
                        // add perturbation
                        double bkup = U_IN_temp[iVar];
                        double delta = Math.Abs(U_IN[iVar]);//, Math.Abs(U_INtm[iVar]) * eps;
                        delta = eps;
                        Debug.Assert(delta > 0);
                        U_IN_temp[iVar] += delta;

                        // flux eval
                        double f1 = m_EdgForm.BoundaryEdgeForm(ref clonedParams, U_IN_temp, null, _vIN, _Grad_vIN);
                        if (double.IsInfinity(f1))
                            throw new ArithmeticException();
                        if (double.IsNaN(f1))
                            throw new ArithmeticException();

                        // restore un-perturbed state
                        U_IN_temp[iVar] = bkup;

                        // compute finite difference
                        double dU_iVar = (f1 - f0) / delta;
                        if (double.IsInfinity(dU_iVar))
                            throw new ArithmeticException();
                        if (double.IsNaN(dU_iVar))
                            throw new ArithmeticException();

                        // inner product
                        ret += dU_iVar * U_IN[iVar];
                        ret -= dU_iVar * U_IN_temp[iVar]; // subtract affine contribution
                    }
                }
            }

            return ret;
        }
    }


}

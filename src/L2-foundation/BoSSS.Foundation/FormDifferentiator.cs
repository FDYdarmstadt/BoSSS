using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation {

    /// <summary>
    /// Differentiation of a volume form, used e.g.to obtain a Jacobian of an operator, see <see cref="SpatialOperator.GetJacobiOperator"/>.
    /// </summary>
    public class VolumeFormDifferentiator : IVolumeForm {

        /// <summary>
        /// %
        /// </summary>
        public VolumeFormDifferentiator(IVolumeForm vf, int SpatialDimension) {
            m_VolForm = vf;
            m_eps = Math.Sqrt(BLAS.MachineEps);

            VolTerms = vf.VolTerms;
            if ((VolTerms & (TermActivationFlags.GradUxV | TermActivationFlags.UxV)) != 0)
                VolTerms |= TermActivationFlags.V;
            if ((VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV)) != 0)
                VolTerms |= TermActivationFlags.GradV;

            ParamUreq = (vf.VolTerms & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0;
            ParamGradUreq = (vf.VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0;
            
            m_ParameterOrdering = new string[0];

            if (ParamUreq && vf.ArgumentOrdering != null && vf.ArgumentOrdering.Count > 0)
                m_ParameterOrdering = m_ParameterOrdering.Cat(vf.ArgumentOrdering.Select(fn => fn + "_lin"));

            OffsetGradUparams = m_ParameterOrdering.Length;
            if (ParamGradUreq && vf.ArgumentOrdering != null && vf.ArgumentOrdering.Count > 0) {
                foreach (string arg in vf.ArgumentOrdering) {
                    m_ParameterOrdering = m_ParameterOrdering.Cat(SpatialDimension.ForLoop(d => arg + "_lin_d[" + d + "]"));
                }
            }

            OffsetOrgParams = m_ParameterOrdering.Length;
            if (vf.ParameterOrdering != null && vf.ParameterOrdering.Count > 0) {
                m_ParameterOrdering = m_ParameterOrdering.Cat(vf.ParameterOrdering);
            }

        }

        bool ParamUreq;

        bool ParamGradUreq;
        int OffsetGradUparams;

        int OffsetOrgParams;

        double m_eps;

        IVolumeForm m_VolForm;

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags VolTerms {
            get;
            private set;
        }

        /// <summary>
        /// %
        /// </summary>
        public IList<string> ArgumentOrdering => m_VolForm.ArgumentOrdering;

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
        /// %
        /// </summary>
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
            ret += f0; // affine contribution - contains V and GradV contribution

            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_VolForm.VolTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)
                    ) {
                    //&& (U[iVar] != 0.0)) { // perf. opt.
                    // add perturbation
                    double bkup = Utmp[iVar];
                    double delta = Math.Abs(U[iVar]) * eps;
                    delta = eps;
                    Debug.Assert(delta > 0);
                    Utmp[iVar] += delta;

                    // flux eval
                    double f1 = m_VolForm.VolumeForm(ref clonedParams, Utmp, null, V, GradV);
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
                    && (GradU.GetRow(iVar).L2NormPow2() != 0.0)) {

                    throw new NotImplementedException("todo");
                }
            }

            return ret;
        }
    }


    /// <summary>
    /// Differentiation of a edge form, used e.g.to obtain a Jacobian of an operator, see <see cref="SpatialOperator.GetJacobiOperator"/>.
    /// </summary>
    public class EdgeFormDifferentiator : IEdgeForm {

        /// <summary>
        /// ctor
        /// </summary>
        public EdgeFormDifferentiator(IEdgeForm ef, int SpatialDimension) {
            m_EdgForm = ef;
            m_eps = Math.Sqrt(BLAS.MachineEps);

            InnerEdgeTerms = ef.InnerEdgeTerms;
            if ((InnerEdgeTerms & (TermActivationFlags.GradUxV | TermActivationFlags.UxV)) != 0)
                InnerEdgeTerms |= TermActivationFlags.V;
            if ((InnerEdgeTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV)) != 0)
                InnerEdgeTerms |= TermActivationFlags.GradV;

            BoundaryEdgeTerms = ef.BoundaryEdgeTerms;
            if ((BoundaryEdgeTerms & (TermActivationFlags.GradUxV | TermActivationFlags.UxV)) != 0)
                BoundaryEdgeTerms |= TermActivationFlags.V;
            if ((BoundaryEdgeTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV)) != 0)
                BoundaryEdgeTerms |= TermActivationFlags.GradV;

            ParamUreq = ((ef.InnerEdgeTerms | ef.BoundaryEdgeTerms) & (TermActivationFlags.UxGradV | TermActivationFlags.UxV)) != 0;
            ParamGradUreq = ((ef.InnerEdgeTerms | ef.BoundaryEdgeTerms) & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0;
            
            m_ParameterOrdering = new string[0];

            if (ParamUreq && ef.ArgumentOrdering != null && ef.ArgumentOrdering.Count > 0)
                m_ParameterOrdering = m_ParameterOrdering.Cat(ef.ArgumentOrdering.Select(fn => fn + "_lin"));

            OffsetGradUparams = m_ParameterOrdering.Length;
            if (ParamGradUreq && ef.ArgumentOrdering != null && ef.ArgumentOrdering.Count > 0) {
                foreach (string arg in ef.ArgumentOrdering) {
                    m_ParameterOrdering = m_ParameterOrdering.Cat(SpatialDimension.ForLoop(d => arg + "_lin_d[" + d + "]"));
                }
            }

            OffsetOrgParams = m_ParameterOrdering.Length;
            if (ef.ParameterOrdering != null && ef.ParameterOrdering.Count > 0) {
                m_ParameterOrdering = m_ParameterOrdering.Cat(ef.ParameterOrdering);
            }

        }

        bool ParamUreq;

        bool ParamGradUreq;
        int OffsetGradUparams;

        int OffsetOrgParams;

        double m_eps;

        IEdgeForm m_EdgForm;

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags InnerEdgeTerms {
            get;
            private set;
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags BoundaryEdgeTerms {
            get;
            private set;
        }

        /// <summary>
        /// %
        /// </summary>
        public IList<string> ArgumentOrdering => m_EdgForm.ArgumentOrdering;

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
            ret += f0; // affine contribution - contains V and GradV contribution
            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_EdgForm.InnerEdgeTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)) {
                    //if (U_IN[iVar] != 0.0) {
                    {
                        // add perturbation
                        double bkup = U_IN_temp[iVar];
                        double delta = Math.Abs(U_IN[iVar]) * eps;
                        delta = eps;
                        Debug.Assert(delta > 0);
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
                        double delta = Math.Abs(U_OT[iVar]) * eps;
                        delta = eps;
                        Debug.Assert(delta > 0);
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

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
                if (ParamUreq) {
                    Utmp[iVar] = Parameters[iVar + OffsetUparams];
                    delta += Utmp[iVar].Pow2();
                }
                if (ParamGradUreq) {
                    for (int d = 0; d < D; d++) {
                        GradUtmp[iVar, d] = Parameters[iVar * D + d + OffsetGradUparams];
                        delta += GradUtmp[iVar, d].Pow2();
                    }
                }
            }

            delta = Math.Max(eps, Math.Sqrt(delta) * eps);
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

            double delta = GetTmpTrialVals(cpv.Parameters, out var Utmp, out var GradUtmp);
            
            CommonParamsVol clonedParams = cpv;
            Debug.Assert(object.ReferenceEquals(cpv, clonedParams) == false);
            GetOrgParams(cpv.Parameters, out clonedParams.Parameters);

            double f0 = m_VolForm.VolumeForm(ref clonedParams, Utmp, GradUtmp, V, GradV);
            ret += f0; // affine contribution - contains V and GradV contribution

            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_VolForm.VolTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)
                    ) {
                    //&& (U[iVar] != 0.0)) { // perf. opt.
                    ret += Diff(ref Utmp[iVar], U[iVar], ref clonedParams, Utmp, GradUtmp, V, GradV, delta, f0);
                }

                if (((m_VolForm.VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxGradV)) != 0)
                    //&& (GradU.GetRow(iVar).L2NormPow2() != 0.0)) {
                    ) {

                    for (int d = 0; d < D; d++) {
                        ret += Diff(ref GradUtmp[iVar, d], GradU[iVar, d], ref clonedParams, Utmp, GradUtmp, V, GradV, delta, f0);
                    }
                }
            }

            return ret;
        }

        private double Diff(ref double PertubVar, double Var,
            ref CommonParamsVol clonedParams, 
            double[] Utmp, double[,] GradUtmp, double V, double[] GradV, 
            double delta, double f0) {
            
            // add perturbation
            double bkup = PertubVar;
            PertubVar += delta;
            Debug.Assert(delta > 0);

            // flux eval
            double f1 = m_VolForm.VolumeForm(ref clonedParams, Utmp, GradUtmp, V, GradV);
#if DEBUG
            if (double.IsInfinity(f1))
                throw new ArithmeticException();
            if (double.IsNaN(f1))
                throw new ArithmeticException();
#endif            

            // compute finite difference
            double dU_iVar = (f1 - f0) / delta;
            if (double.IsInfinity(dU_iVar))
                throw new ArithmeticException();
            if (double.IsNaN(dU_iVar))
                throw new ArithmeticException();

            // restore un-perturbed state
            PertubVar = bkup;

            // inner product
            double ret = 0;
            ret += dU_iVar * Var;
            ret -= dU_iVar * PertubVar; // subtract affine contribution
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

        //static public int SetDir = 0;
        static public int Dir = 0;


        public double InnerEdgeForm(ref CommonParams inp, double[] U_IN, double[] U_OT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double ret = 0.0;
            int GAMMA = m_EdgForm.ArgumentOrdering.Count;
            int D = inp.D;
            Debug.Assert(D == m_SpatialDimension, "spatial dimension mismatch");
            
            CommonParams clonedParams = inp;
            Debug.Assert(object.ReferenceEquals(inp, clonedParams) == false);
            GetOrgParams(inp.Parameters_IN, out clonedParams.Parameters_IN);
            GetOrgParams(inp.Parameters_OUT, out clonedParams.Parameters_OUT);
            
            double deltaIn = GetTmpTrialVals(inp.Parameters_IN, out var U_IN_temp, out var GradU_IN_temp);
            double deltaOt = GetTmpTrialVals(inp.Parameters_OUT, out var U_OT_temp, out var GradU_OT_temp);
            double delta = Math.Max(deltaIn, deltaOt);


            //SetDir = 0;
            
            double f0 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
            
            ret += f0; // affine contribution - contains V and GradV contribution

            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_EdgForm.InnerEdgeTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)) {
                    //if (U_IN[iVar] != 0.0) {
                    {
                        //if (U_IN[iVar] != 0.0 && (_vIN != 0 || _vOUT != 0))
                        //    Console.Write("");

                        ret += Diff(ref U_IN_temp[iVar], U_IN[iVar], 
                            ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT, 
                            delta, f0);

                       
                    }

                    //if(U_OT[iVar] != 0.0) {
                    {
                        //if (U_OT[iVar] != 0.0 && (_vIN != 0 || _vOUT != 0))
                        //    Console.Write("");

                        ret += Diff(ref U_OT_temp[iVar], U_OT[iVar], 
                            ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT, 
                            delta, f0);

                       
                    }
                }

                if (((m_EdgForm.InnerEdgeTerms & (TermActivationFlags.GradUxV | TermActivationFlags.GradUxGradV)) != 0)) {
                    for (int d = 0; d < D; d++) {
                        //if (U_IN[iVar] != 0.0) {
                        {
                            ret += Diff(ref GradU_IN_temp[iVar, d], _Grad_uIN[iVar, d], 
                                ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT, 
                                delta, f0);
                        }

                        //if(U_OT[iVar] != 0.0) {
                        {
                            ret += Diff(ref GradU_OT_temp[iVar, d], _Grad_uOUT[iVar, d], 
                                ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT, 
                                delta, f0);
                        }
                    }
                }

            }

            //Dir = 0;
            return ret;
        }

        private double Diff(
            ref double PertubVar, double Var,
            ref CommonParams clonedParams, 
            double[] U_IN_temp, double[] U_OT_temp, double[,] GradU_IN_temp, double[,] GradU_OT_temp,
            double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT, double delta, double f0) {

            // add perturbation
            double bkup = PertubVar;
            PertubVar += delta;

            // flux eval
            double f1 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
#if DEBUG
            if (double.IsInfinity(f1))
                throw new ArithmeticException();
            if (double.IsNaN(f1))
                throw new ArithmeticException();
#endif

            // compute finite difference
            double dU_iVar = (f1 - f0) / delta;
            if (double.IsInfinity(dU_iVar))
                throw new ArithmeticException();
            if (double.IsNaN(dU_iVar))
                throw new ArithmeticException();

            // restore un-perturbed state
            PertubVar = bkup;

            // inner product
            double ret = 0;
            ret += dU_iVar * Var;
            ret -= dU_iVar * PertubVar; // subtract affine contribution
            return ret;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] U_IN, double[,] _Grad_uIn, double _vIN, double[] _Grad_vIN) {
            double ret = 0.0;
            int GAMMA = m_EdgForm.ArgumentOrdering.Count;
            int D = inp.D;
            Debug.Assert(D == m_SpatialDimension, "spatial dimension mismatch");

            CommonParamsBnd clonedParams = inp;
            Debug.Assert(object.ReferenceEquals(inp, clonedParams) == false);
            GetOrgParams(inp.Parameters_IN, out clonedParams.Parameters_IN);
            
            double delta = GetTmpTrialVals(inp.Parameters_IN, out var U_IN_temp, out var GradU_IN_temp);

            double f0 = m_EdgForm.BoundaryEdgeForm(ref clonedParams, U_IN_temp, GradU_IN_temp, _vIN, _Grad_vIN);
            ret += f0; // affine contribution - contains V and GradV contribution

            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_EdgForm.InnerEdgeTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)) {
                    //if (U_IN[iVar] != 0.0) {
                    {
                        ret += DiffBnd(ref U_IN_temp[iVar], U_IN[iVar], 
                            ref clonedParams, U_IN_temp, GradU_IN_temp, _vIN, _Grad_vIN, delta, f0);
                    }
                }

                if (((m_EdgForm.InnerEdgeTerms & (TermActivationFlags.GradUxV | TermActivationFlags.GradUxGradV)) != 0)) {
                    for (int d = 0; d < D; d++) {
                        //if (U_IN[iVar] != 0.0) {
                        {
                            ret += DiffBnd(ref GradU_IN_temp[iVar, d], _Grad_uIn[iVar, d], 
                                ref clonedParams, U_IN_temp, GradU_IN_temp, _vIN, _Grad_vIN, delta, f0);
                        }
                    }
                }
            }

            return ret;
        }

        private double DiffBnd(ref double PertubVar, double Var,
            ref CommonParamsBnd clonedParams, double[] U_IN_temp, double[,] GradU_IN_temp, double _vIN, double[] _Grad_vIN, double delta, double f0) {

            // add perturbation
            double bkup = PertubVar;
            PertubVar += delta;
            Debug.Assert(delta > 0);

            // flux eval
            double f1 = m_EdgForm.BoundaryEdgeForm(ref clonedParams, U_IN_temp, GradU_IN_temp, _vIN, _Grad_vIN);
#if DEBUG
            if (double.IsInfinity(f1))
                throw new ArithmeticException();
            if (double.IsNaN(f1))
                throw new ArithmeticException();
#endif

            // compute finite difference
            double dU_iVar = (f1 - f0) / delta;
            if (double.IsInfinity(dU_iVar))
                throw new ArithmeticException();
            if (double.IsNaN(dU_iVar))
                throw new ArithmeticException();
            
            // restore un-perturbed state
            PertubVar = bkup;

            // inner product
            double ret = 0;
            ret += dU_iVar * Var;
            ret -= dU_iVar * PertubVar; // subtract affine contribution
            return ret;

        }
    }
}

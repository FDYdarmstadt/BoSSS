using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Base-Functionality for utilities for computation of Spatial Operator Jacobians, (<see cref="DifferentialOperator.GetJacobiOperator"/>,
    /// i.e. approximate differentiation of equation components.
    /// </summary>
    abstract public class FormDifferentiatorCommon : IEquationComponent, IEquationComponentChecking, IEquationComponentCoefficient, IParameterHandling, ISpeciesFilter, IEquationComponentSpeciesNotification, ILevelSetEquationComponentCoefficient, IMultitreadSafety {


        /// <summary>
        /// %
        /// </summary>
        protected FormDifferentiatorCommon(IEquationComponent vf, TermActivationFlags __Terms, int SpatialDimension) {
            m_OrgForm = vf;
            m_eps = Math.Sqrt(BLAS.MachineEps);
            // DIAGNOSTIC (coordinator-requested FD-step-size sensitivity check, this round):
            // allow overriding the FD differentiation step size to test empirically whether a
            // specific M-entry's small discrepancy against an independent reconstruction scales
            // with the FD step (truncation-error signature) or stays fixed (pointing elsewhere).
            {
                string epsOverride = System.Environment.GetEnvironmentVariable("DIAG_FD_EPS_OVERRIDE");
                if (epsOverride != null && double.TryParse(epsOverride, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double epsVal))
                    m_eps = epsVal;
            }
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
                NoOfOrgParams = vf.ParameterOrdering.Count;
            } else {
                NoOfOrgParams = 0;
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
        /// required terms for the linearization 
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
        /// number of parameters in original component
        /// </summary>
        protected int NoOfOrgParams;

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

        /// <summary>
        /// passes the coefficients to original form
        /// </summary>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if (m_OrgForm is IEquationComponentCoefficient eqc) {
                eqc.CoefficientUpdate(cs, DomainDGdeg, TestDGdeg);
            }
        }


        /// <summary>
        /// redirects to parameter update of original form; Jacobi parameter are updated in <see cref="JacobianParamUpdate"/>
        /// </summary>
        public void MyParameterUpdate(DGField[] Arguments, DGField[] Parameters) {

            if (m_OrgForm is IParameterHandling ph) {
                var orgParameters = Parameters.GetSubVector(OffsetOrgParams, NoOfOrgParams);
                Parameters.SetAll(default(DGField));

                ph.MyParameterUpdate(Arguments, orgParameters);
                Parameters.SetSubVector(orgParameters, OffsetOrgParams, NoOfOrgParams); // restore what we got from original component
            } else {
                // nop
                Parameters.SetAll(default(DGField)); // very secret signal that the parameter was not really updated
            }
        }

        /// <summary>
        /// redirects to parameter update of original form; Jacobi parameter are updated in <see cref="JacobianParamUpdate"/>
        /// </summary>
        public DGField[] MyParameterAlloc(DGField[] Arguments) {
            if (m_OrgForm is IParameterHandling ph) {
                var orgRet = ph.MyParameterAlloc(Arguments);
                var ret = new DGField[this.m_ParameterOrdering.Length];
                ret.SetSubVector(orgRet, OffsetOrgParams, NoOfOrgParams);
                return ret;
            } else {
                var ParamNames = this.ParameterOrdering;
                return ParamNames != null ? new DGField[ParamNames.Count] : new DGField[0];
            }
        }

        /// <summary>
        /// <see cref="IEquationComponentSpeciesNotification.SetParameter"/>
        /// </summary>
        public void SetParameter(string speciesName) {
            if(m_OrgForm is IEquationComponentSpeciesNotification ecsn) {
                ecsn.SetParameter(speciesName);
            }
        }

        virtual public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            if (m_OrgForm is ILevelSetEquationComponentCoefficient eqc) {
                eqc.CoefficientUpdate(csA,csB,DomainDGdeg,TestDGdeg);
            }
        }

        /// <summary>
        /// <see cref="ISpeciesFilter.ValidSpecies"/>
        /// </summary>
        public string ValidSpecies {
            get {
                if (m_OrgForm is ISpeciesFilter sf) {
                    return sf.ValidSpecies;
                } else {
                    return null;
                }

            }
        }


        public override string ToString() {
            return "XDG-" + this.GetType().Name + ": " + m_OrgForm.ToString() + "(" + m_OrgForm.GetType().FullName + ")";
        }

        bool m_IsMultithreadSafe;

        public bool IsMultithreadSafe {
            get {
                SetupMultithread();
                return m_IsMultithreadSafe;
            }
        }

        bool m_SetupMultithread_executed = false;

        void SetupMultithread() {
            if (m_SetupMultithread_executed == true)
                return;
            m_SetupMultithread_executed = true;

            if (this.m_OrgForm is IMultitreadSafety ms) {
                if (ms.IsMultithreadSafe) {
                    // nothing to do
                    m_IsMultithreadSafe = true;
                } else {
                    var clone = ms.CloneForThread();
                    if (clone != null) {
                        m_IsMultithreadSafe = true;
                    } else {
                        m_IsMultithreadSafe = false;
                    }
                }

                m_SetupMultithread_executed = true;
            } else {
                // per default we just assume nothing goes wrong.
                m_IsMultithreadSafe = true;
            }
        }

        protected abstract FormDifferentiatorCommon Clone();

        public IEquationComponent CloneForThread() {
            SetupMultithread();
            if((m_OrgForm as IMultitreadSafety)?.CloneForThread() != null) 
                return Clone();
            else
                return null;
        }

        public object GetPadlock() {
            SetupMultithread();
            if (this.m_OrgForm is IMultitreadSafety ms)
                return ms.GetPadlock();
            else
                return this.m_OrgForm;
        }
    }



    /// <summary>
    /// Differentiation of a volume form, used to obtain a Jacobian of an operator, see <see cref="DifferentialOperator.GetJacobiOperator"/>.
    /// </summary>
    public class VolumeFormDifferentiator : FormDifferentiatorCommon, IVolumeForm {

        IVolumeForm m_VolForm;

        // DIAGNOSTIC (output cap added retroactively, coordinator-requested housekeeping):
        // DIAG_EULER_VOL_CHECK has no built-in per-cell/edge scoping in some call patterns and
        // was found to generate on the order of 1e6+ lines (hundreds of MB) over a full
        // multi-iteration solve when left unbounded, which previously filled the machine's temp
        // drive. Caps total output the same way this session's newer diagnostics do.
        private static int s_eulerVolCheckCount = 0;

        protected override FormDifferentiatorCommon Clone() {
            var VolFormClone = (m_VolForm as IMultitreadSafety)?.CloneForThread() as IVolumeForm;
            if (VolFormClone != null)
                return new VolumeFormDifferentiator(VolFormClone, m_SpatialDimension);
            else
                return null;
        }


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

            // DIAGNOSTIC (coordinator-requested root-cause confirmation, final round): dump the
            // "_lin_d[d]" SNAPSHOT gradient (GradUtmp, sourced from JacobianParamUpdate's
            // L2-projected "Derivative()" parameter field, NOT the exact per-node invJac-
            // transformed gradient) for a target cell, to directly compare against an
            // independently-computed EXACT gradient at the same physical point -- tests the
            // hypothesis that this snapshot is a degree-preserving L2 projection of the true
            // (rational, non-polynomial on curved cells) gradient, hence inexact specifically in
            // curvature-coupled directions, which would explain why cross-variable (e.g.
            // Velocity's) derivatives that MULTIPLY by this snapshot come out slightly wrong on
            // curved cells while same-variable (Temperature's own) derivatives -- which only test
            // sensitivity to a perturbation, not the snapshot's absolute value -- stay exact.
            if (System.Environment.GetEnvironmentVariable("DIAG_LINSNAPSHOT_DUMP") == "1") {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_LINSNAPSHOT_CELL");
                int targetCell = targetStr != null ? int.Parse(targetStr) : -1;
                if (targetCell < 0 || cpv.jCell == targetCell) {
                    Console.Error.WriteLine(string.Format(
                        "DIAG_LINSNAPSHOT: form={0} jCell={1} X=({2:G17},{3:G17},{4:G17}) GradUtmp[0,:]=[{5:G17},{6:G17},{7:G17}]",
                        m_VolForm.GetType().Name, cpv.jCell, cpv.Xglobal[0], cpv.Xglobal[1], cpv.Xglobal[2],
                        GradUtmp[0, 0], D > 1 ? GradUtmp[0, 1] : double.NaN, D > 2 ? GradUtmp[0, 2] : double.NaN));
                }
            }

            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables

                if (((m_VolForm.VolTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)
                    ) {
                    //&& (U[iVar] != 0.0)) { // perf. opt.
                    ret += Diff(ref Utmp[iVar], U[iVar], ref clonedParams, Utmp, GradUtmp, V, GradV, delta, f0, -1, iVar);
                }

                if (((m_VolForm.VolTerms & (TermActivationFlags.GradUxGradV | TermActivationFlags.GradUxV)) != 0)
                    //&& (GradU.GetRow(iVar).L2NormPow2() != 0.0)) {
                    ) {

                    for (int d = 0; d < D; d++) {
                        ret += Diff(ref GradUtmp[iVar, d], GradU[iVar, d], ref clonedParams, Utmp, GradUtmp, V, GradV, delta, f0, d, iVar);
                    }
                }
            }

            // DIAGNOSTIC (env-gated, zero cost unless set): same Euler-homogeneity identity check
            // as EdgeFormDifferentiator's DIAG_EULER_CHECK, but for the VOLUME differentiator --
            // the edge-form check came back clean on both flat and curved faces, but isolating
            // VolumeFormDifferentiator alone (DIVFORM_ISOLATION_NOEDGE=1) reproduces the FULL
            // curved-cell aggregate residual, so the defect must be here instead.
            if (System.Environment.GetEnvironmentVariable("DIAG_EULER_VOL_CHECK") == "1" && s_eulerVolCheckCount < 400) {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_EULER_VOL_CHECK_CELL");
                int targetCell = targetStr != null ? int.Parse(targetStr) : -1;
                if (targetCell < 0 || cpv.jCell == targetCell) {
                    System.Threading.Interlocked.Increment(ref s_eulerVolCheckCount);
                    bool probeIsZero = U.All(v => v == 0.0);
                    Console.WriteLine(string.Format(
                        "DIAG_EULER_VOL: form={11} jCell={0} X=({1:0.########},{2:0.########},{3:0.########}) GradV=({4:0.########},{5:0.########},{6:0.########}) probeZero={7} U=[{8}] V={9:0.####} f0={10:0.##########} ret={12:0.##########} ret+f0={13:0.####E-00}",
                        cpv.jCell, cpv.Xglobal[0], cpv.Xglobal[1], cpv.Xglobal[2],
                        GradV[0], GradV[1], GradV[2],
                        probeIsZero,
                        string.Join(",", U.Select(v => v.ToString("0.####"))),
                        V, f0, m_VolForm.ToString(), ret, ret + f0));
                }
            }

            return ret;
        }

        private double Diff(ref double PertubVar, double Var,
            ref CommonParamsVol clonedParams,
            double[] Utmp, double[,] GradUtmp, double V, double[] GradV,
            double delta, double f0,
            int d, int iVar) {

            // add perturbation
            double bkup = PertubVar;
            PertubVar += delta;
            Debug.Assert(delta > 0);

            // flux eval
            double f1 = m_VolForm.VolumeForm(ref clonedParams, Utmp, GradUtmp, V, GradV);

            /*
            double[] f1x = new double[100];
            var cclonedParams = clonedParams;
            void ParaFlux (int iThread) {
                f1x[iThread] = m_VolForm.VolumeForm(ref cclonedParams, Utmp, GradUtmp, V, GradV);
            }

            Parallel.For(0, f1x.Length, ParaFlux);

            double totErr = f1x.Select(f1x_i => (f1x_i - f1).Pow2()).Sum();
            if(totErr != 0) {
                throw new ArithmeticException();
            }
            */

#if DEBUG
            if (double.IsInfinity(f1))
                throw new ArithmeticException();
            if (double.IsNaN(f1))
                throw new ArithmeticException();
#endif            

            // compute finite difference
            double dU_iVar = (f1 - f0) / delta;
            if (double.IsInfinity(dU_iVar))
                throw new ArithmeticException("Got INF while differentiation of volume part of " + m_VolForm.GetType().Name);
            if (double.IsNaN(dU_iVar))
                throw new ArithmeticException("Got NAN while differentiation of volume part of " + m_VolForm.GetType().Name);

            // restore un-perturbed state
            PertubVar = bkup;

            // test
            /*
            if(V == 1 && GradV.L2Norm() == 0) {
                double u = Utmp[1];
                double v = Utmp[2];
                double Tx = GradUtmp[0, 0];
                double Ty = GradUtmp[0, 1];

                double soll;
                if(d < 0 && iVar == 1) {
                    soll = Tx;
                } else if(d < 0 && iVar == 2) {
                    soll = Ty;
                } else if(d == 0 && iVar == 0) {
                    soll = u;
                } else if(d == 1 && iVar == 0) {
                    soll = v;
                } else {
                    soll = 0;
                }

                Debug.Assert((dU_iVar - soll).Abs() < delta*1000);

            }
            */



            // inner product
            double ret = 0;
            ret += dU_iVar * Var;
            ret -= dU_iVar * PertubVar; // subtract affine contribution
            return ret;
        }


        public override string ToString() {
            return base.ToString();
        }
    }


    /// <summary>
    /// Differentiation of a edge form, used e.g.to obtain a Jacobian of an operator, see <see cref="DifferentialOperator.GetJacobiOperator"/>.
    /// </summary>
    public class EdgeFormDifferentiator : FormDifferentiatorCommon, IEdgeForm {

        IEdgeForm m_EdgForm;

        // DIAGNOSTIC (coordinator-requested FD-differentiation metric-consistency check, this
        // round): output cap shared across InnerEdgeForm's INIT/AFTER-F0 prints and Diff's
        // BEFORE-F1/AFTER-F1 prints, so DIAG_METRIC_CONSISTENCY=1 can't flood the log across a
        // full multi-iteration Newton solve -- only the first few call sequences are needed to
        // answer whether Normal/X ever change within one InnerEdgeForm invocation.
        private static int s_metricConsistencyCount = 0;

        // DIAGNOSTIC (output cap added retroactively, coordinator-requested housekeeping):
        // DIAG_EULER_CHECK (both the InnerEdgeForm and BoundaryEdgeForm sites below) has no
        // built-in cap and was found to generate on the order of 1e6+ lines (hundreds of MB to
        // ~1GB) over a full multi-iteration solve when scoped only by cell (not edge), which
        // previously filled the machine's temp drive. Shared across both call sites, capped the
        // same way this session's newer diagnostics are.
        private static int s_eulerCheckCount = 0;

        protected override FormDifferentiatorCommon Clone() {
            var EdgFormClone = (m_EdgForm as IMultitreadSafety)?.CloneForThread() as IEdgeForm;
            if (EdgFormClone != null)
                return new EdgeFormDifferentiator(EdgFormClone, m_SpatialDimension);
            else
                return null;
        }



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
            int D = inp.D;
            Debug.Assert(D == m_SpatialDimension, "spatial dimension mismatch");

            CommonParams clonedParams = inp;
            Debug.Assert(object.ReferenceEquals(inp, clonedParams) == false);
            GetOrgParams(inp.Parameters_IN, out clonedParams.Parameters_IN);
            GetOrgParams(inp.Parameters_OUT, out clonedParams.Parameters_OUT);

            double deltaIn = GetTmpTrialVals(inp.Parameters_IN, out var U_IN_temp, out var GradU_IN_temp);
            double deltaOt = GetTmpTrialVals(inp.Parameters_OUT, out var U_OT_temp, out var GradU_OT_temp);
            double delta = Math.Max(deltaIn, deltaOt);

            // DIAGNOSTIC (temporary, env-var gated, zero cost unless set): dump the raw Utmp/GradUtmp
            // actually read from the "_lin"/"_lin_d" parameter buffers for a specific target edge, to
            // directly check whether the DIAGONAL (cr=0,cc=0 / cr=1,cc=1) block's inputs match the
            // exact analytic linearization-point gradient at this quadrature point -- independent of
            // any periodic cross-block rotation question.
            if (System.Environment.GetEnvironmentVariable("DIAG_DIAG_BLOCK") == "1") {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_DIAG_BLOCK_EDGE");
                if (targetStr != null && int.TryParse(targetStr, out int targetEdge) && inp.iEdge == targetEdge) {
                    Console.WriteLine(string.Format(
                        "DIAG_DIAGBLK: dispatch={14} edge={0} jIn={1} jOt={2} X=({3:0.######},{4:0.######},{5:0.######}) U_IN_temp=[{6}] U_OT_temp=[{7}] GradU_IN_temp[1,:]=[{8:0.########},{9:0.########},{10:0.########}] GradU_OT_temp[1,:]=[{11:0.########},{12:0.########},{13:0.########}]",
                        inp.iEdge, inp.jCellIn, inp.jCellOut, inp.X[0], inp.X[1], inp.X[2],
                        string.Join(",", U_IN_temp.Select(v => v.ToString("0.########"))),
                        string.Join(",", U_OT_temp.Select(v => v.ToString("0.########"))),
                        GradU_IN_temp.GetLength(0) > 1 ? GradU_IN_temp[1, 0] : double.NaN,
                        GradU_IN_temp.GetLength(0) > 1 ? GradU_IN_temp[1, 1] : double.NaN,
                        GradU_IN_temp.GetLength(0) > 1 ? GradU_IN_temp[1, 2] : double.NaN,
                        GradU_OT_temp.GetLength(0) > 1 ? GradU_OT_temp[1, 0] : double.NaN,
                        GradU_OT_temp.GetLength(0) > 1 ? GradU_OT_temp[1, 1] : double.NaN,
                        GradU_OT_temp.GetLength(0) > 1 ? GradU_OT_temp[1, 2] : double.NaN,
                        BoSSS.Foundation.Quadrature.Linear.DiagDispatchTag.CurrentTermType ?? "?"));
                }
            }

            //SetDir = 0;

            // DIAGNOSTIC (coordinator-requested FD-differentiation metric-consistency check,
            // this round): geometry (Normal, X) is a property of the quadrature NODE, not of
            // the perturbed trial DOF value -- mathematically, clonedParams.Normal/.X must be
            // bit-identical across the f0 call below and every f1 call inside Diff() for this
            // same InnerEdgeForm invocation, since only U_IN_temp/U_OT_temp/GradU_*_temp
            // entries are perturbed, never clonedParams itself. But clonedParams is passed
            // `ref` into m_EdgForm.InnerEdgeForm, so a real physics form COULD in principle
            // (bug or not) write back into it. Verify empirically rather than assuming.
            bool diagMetricConsistency = System.Environment.GetEnvironmentVariable("DIAG_METRIC_CONSISTENCY") == "1"
                && System.Environment.GetEnvironmentVariable("DIAG_METRIC_CONSISTENCY_CELL") is string mcCellStr
                && int.TryParse(mcCellStr, out int mcCell)
                && (inp.jCellIn == mcCell || inp.jCellOut == mcCell)
                && s_metricConsistencyCount < 400;
            if (diagMetricConsistency) {
                System.Threading.Interlocked.Increment(ref s_metricConsistencyCount);
                Console.WriteLine(string.Format(
                    "DIAG_METRIC_CONSISTENCY[INIT]: edge={0} jIn={1} jOt={2} form={7} Normal=({3:G17},{4:G17},{5:G17}) X=({6:G17})",
                    inp.iEdge, inp.jCellIn, inp.jCellOut, clonedParams.Normal[0], clonedParams.Normal[1], clonedParams.Normal[2],
                    string.Join(",", Enumerable.Range(0, clonedParams.X.Dim).Select(d => clonedParams.X[d].ToString("G17"))),
                    m_EdgForm.GetType().Name));
            }

            double f0 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);

            if (diagMetricConsistency) {
                System.Threading.Interlocked.Increment(ref s_metricConsistencyCount);
                Console.WriteLine(string.Format(
                    "DIAG_METRIC_CONSISTENCY[AFTER-F0]: edge={0} jIn={1} jOt={2} Normal=({3:G17},{4:G17},{5:G17}) X=({6:G17}) f0={7:G17}",
                    inp.iEdge, inp.jCellIn, inp.jCellOut, clonedParams.Normal[0], clonedParams.Normal[1], clonedParams.Normal[2],
                    string.Join(",", Enumerable.Range(0, clonedParams.X.Dim).Select(d => clonedParams.X[d].ToString("G17"))),
                    f0));
            }

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

            //if(debug) {
            //    double retAlt = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN, U_OT, _Grad_uIN, _Grad_uOUT, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
            //    if(retAlt != 0.0)
            //        Console.Write("");
            //    Debug.Assert((ret - retAlt).Abs() <= 1e-10);
            //}

            // DIAGNOSTIC (temporary, env-var gated): dump the actual probe values (_vIN/_vOUT,
            // the test-function slots LinearFormVectorizer's GetSourceProbe/GetCoeffWithSourceProbe
            // manipulate) and the final RETURNED "ret" value for this specific edge -- to see
            // exactly what EdgeFormDifferentiator.InnerEdgeForm hands back to the affine ("b")
            // assembly's V[l,k,cr] slots, and whether it matches the expected
            // f0-sum(dU*Utmp) affine identity given the (already-confirmed-correct) Utmp/GradUtmp.
            if (System.Environment.GetEnvironmentVariable("DIAG_DIAG_BLOCK") == "1") {
                string targetStr2 = System.Environment.GetEnvironmentVariable("DIAG_DIAG_BLOCK_EDGE");
                if (targetStr2 != null && int.TryParse(targetStr2, out int targetEdge2) && inp.iEdge == targetEdge2) {
                    Console.WriteLine(string.Format(
                        "DIAG_DIAGBLK_RET: dispatch={0} edge={1} jIn={2} jOt={3} X=({4:0.######},{5:0.######},{6:0.######}) _vIN={7:0.########} _vOUT={8:0.########} f0={9:0.########} ret={10:0.########}",
                        BoSSS.Foundation.Quadrature.Linear.DiagDispatchTag.CurrentTermType ?? "?",
                        inp.iEdge, inp.jCellIn, inp.jCellOut, inp.X[0], inp.X[1], inp.X[2],
                        _vIN, _vOUT, f0, ret));
                }
            }

            // DIAGNOSTIC (env-gated, zero cost unless set): direct Euler-homogeneity identity
            // check -- for a degree-2-homogeneous two-sided form (r homogeneous in (Uin,Uout)
            // jointly, at fixed Normal/v), Euler's theorem predicts ret(U_IN=U_OT=0) == -f0
            // exactly. Dumps the OUTER probe U_IN/U_OT (as received by THIS differentiator call,
            // distinct from the internal Utmp/GradUtmp used to compute f0), inp.Normal, f0, ret,
            // and ret+f0 (should be ~0 for an all-zero probe) for a specific target edge, to
            // directly see whether the identity itself breaks on curved cells, or holds per-call
            // and the defect is elsewhere (quadrature summation/weighting across nodes).
            if (System.Environment.GetEnvironmentVariable("DIAG_EULER_CHECK") == "1" && s_eulerCheckCount < 400) {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_EULER_CHECK_EDGE");
                string targetCellStr = System.Environment.GetEnvironmentVariable("DIAG_EULER_CHECK_CELL");
                bool cellMatch = targetCellStr == null || (int.TryParse(targetCellStr, out int targetCell0) && (inp.jCellIn == targetCell0 || inp.jCellOut == targetCell0));
                if (cellMatch && (targetStr == null || (int.TryParse(targetStr, out int targetEdge) && inp.iEdge == targetEdge))) {
                    System.Threading.Interlocked.Increment(ref s_eulerCheckCount);
                    bool probeIsZero = U_IN.All(v => v == 0.0) && U_OT.All(v => v == 0.0);
                    bool isPeriodic = false;
                    try {
                        var tags = inp.GridDat.iGeomEdges.EdgeTags;
                        if (inp.iEdge >= 0 && inp.iEdge < tags.Length)
                            isPeriodic = tags[inp.iEdge] >= BoSSS.Foundation.Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG;
                    } catch { }
                    // DIAGNOSTIC (coordinator-requested kink-proximity classification, this
                    // round): compute the actual normal-velocity quantities (vINxN/vOUTxN) that
                    // drive this class of form's Math.Min/Max Roe-type kink, using U_IN_temp/
                    // U_OT_temp (== Utmp at the point f0/ret were just computed -- restored to
                    // their unperturbed values by Diff() after every probe, so still valid here).
                    // velOffset=1 for HeatEq's [Temperature,VelX,VelY,VelZ] ArgumentOrdering,
                    // velOffset=0 for momentum's [VelX,VelY,VelZ]-only ArgumentOrdering.
                    double vINxN_diag = double.NaN, vOUTxN_diag = double.NaN;
                    try {
                        string formName = m_EdgForm.GetType().Name;
                        int velOffset = formName.Contains("HeatConvection") ? 1 : (formName.Contains("Convection") ? 0 : -1);
                        if (velOffset >= 0 && velOffset + D <= U_IN_temp.Length && velOffset + D <= U_OT_temp.Length) {
                            vINxN_diag = 0; vOUTxN_diag = 0;
                            for (int dd = 0; dd < D; dd++) {
                                vINxN_diag += U_IN_temp[velOffset + dd] * clonedParams.Normal[dd];
                                vOUTxN_diag += U_OT_temp[velOffset + dd] * clonedParams.Normal[dd];
                            }
                        }
                    } catch { }
                    Console.WriteLine(string.Format(
                        "DIAG_EULER: form={17} edge={0} jIn={1} jOt={2} periodic={18} X=({3:0.########},{4:0.########},{5:0.########}) N=({6:0.########},{7:0.########},{8:0.########}) probeZero={9} U_IN=[{10}] U_OT=[{11}] vIN={12:0.####} vOT={13:0.####} f0={14:0.##########} ret={15:0.##########} ret+f0={16:0.####E-00} vINxN={19:G6} vOUTxN={20:G6}",
                        inp.iEdge, inp.jCellIn, inp.jCellOut,
                        inp.X[0], inp.X[1], inp.X[2],
                        inp.Normal[0], inp.Normal[1], inp.Normal[2],
                        probeIsZero,
                        string.Join(",", U_IN.Select(v => v.ToString("0.####"))),
                        string.Join(",", U_OT.Select(v => v.ToString("0.####"))),
                        _vIN, _vOUT, f0, ret, ret + f0, m_EdgForm.ToString(), isPeriodic,
                        vINxN_diag, vOUTxN_diag));
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

            // DIAGNOSTIC (coordinator-requested FD-differentiation metric-consistency check,
            // this round): see the matching INIT/AFTER-F0 prints in InnerEdgeForm above -- same
            // filter/cap, applied here to every perturbed (f1) call.
            bool diagMetricConsistencyD = System.Environment.GetEnvironmentVariable("DIAG_METRIC_CONSISTENCY") == "1"
                && System.Environment.GetEnvironmentVariable("DIAG_METRIC_CONSISTENCY_CELL") is string mcCellStrD
                && int.TryParse(mcCellStrD, out int mcCellD)
                && (clonedParams.jCellIn == mcCellD || clonedParams.jCellOut == mcCellD)
                && s_metricConsistencyCount < 400;
            if (diagMetricConsistencyD) {
                System.Threading.Interlocked.Increment(ref s_metricConsistencyCount);
                var xLocalBefore = clonedParams.X;
                Console.WriteLine(string.Format(
                    "DIAG_METRIC_CONSISTENCY[BEFORE-F1]: edge={0} jIn={1} jOt={2} Normal=({3:G17},{4:G17},{5:G17}) X=({6:G17})",
                    clonedParams.iEdge, clonedParams.jCellIn, clonedParams.jCellOut,
                    clonedParams.Normal[0], clonedParams.Normal[1], clonedParams.Normal[2],
                    string.Join(",", Enumerable.Range(0, xLocalBefore.Dim).Select(d => xLocalBefore[d].ToString("G17")))));
            }

            // flux eval
            double f1 = m_EdgForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);

            if (diagMetricConsistencyD) {
                System.Threading.Interlocked.Increment(ref s_metricConsistencyCount);
                var xLocalAfter = clonedParams.X;
                Console.WriteLine(string.Format(
                    "DIAG_METRIC_CONSISTENCY[AFTER-F1]: edge={0} jIn={1} jOt={2} Normal=({3:G17},{4:G17},{5:G17}) X=({6:G17}) f1={7:G17}",
                    clonedParams.iEdge, clonedParams.jCellIn, clonedParams.jCellOut,
                    clonedParams.Normal[0], clonedParams.Normal[1], clonedParams.Normal[2],
                    string.Join(",", Enumerable.Range(0, xLocalAfter.Dim).Select(d => xLocalAfter[d].ToString("G17"))),
                    f1));
            }
#if DEBUG
            if (double.IsInfinity(f1))
                throw new ArithmeticException();
            if (double.IsNaN(f1))
                throw new ArithmeticException();
#endif

            // compute finite difference
            double dU_iVar = (f1 - f0) / delta;
            if (double.IsInfinity(dU_iVar))
                throw new ArithmeticException("Got INF while differentiation of inner edge part of " + m_EdgForm.GetType().Name);
            if (double.IsNaN(dU_iVar))
                throw new ArithmeticException("Got NAN while differentiation of inner edge part of " + m_EdgForm.GetType().Name);

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
                if (((m_EdgForm.BoundaryEdgeTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)) {
                    //if (U_IN[iVar] != 0.0) {
                    {
                        ret += DiffBnd(ref U_IN_temp[iVar], U_IN[iVar],
                            ref clonedParams, U_IN_temp, GradU_IN_temp, _vIN, _Grad_vIN, delta, f0);
                    }
                }

                if (((m_EdgForm.BoundaryEdgeTerms & (TermActivationFlags.GradUxV | TermActivationFlags.GradUxGradV)) != 0)) {
                    for (int d = 0; d < D; d++) {
                        //if (U_IN[iVar] != 0.0) {
                        {
                            ret += DiffBnd(ref GradU_IN_temp[iVar, d], _Grad_uIn[iVar, d],
                                ref clonedParams, U_IN_temp, GradU_IN_temp, _vIN, _Grad_vIN, delta, f0);
                        }
                    }
                }
            }

            //if(debug) {
            //    double retAlt = m_EdgForm.BoundaryEdgeForm(ref clonedParams, U_IN, _Grad_uIn, _vIN, _Grad_vIN);
            //    if(retAlt != 0.0)
            //        Console.Write("");
            //    Debug.Assert((ret - retAlt).Abs() <= 1e-10);
            //}

            // DIAGNOSTIC (coordinator-requested check, this round): same
            // Euler-homogeneity/self-consistency identity check as
            // InnerEdgeForm's own DIAG_EULER_CHECK above, but for
            // BoundaryEdgeForm -- the newly-discovered contributor
            // (convection's ConstantTemperature wall term) that has never
            // been checked with this rigor. Dumps the outer probe U_IN (as
            // received by this differentiator call), f0, ret, and ret+f0
            // (should be ~0 for an all-zero probe, mirroring the interior
            // form's already-established ret+f0~0 identity) for a specific
            // target edge.
            if (System.Environment.GetEnvironmentVariable("DIAG_EULER_CHECK") == "1" && s_eulerCheckCount < 400) {
                string targetStr = System.Environment.GetEnvironmentVariable("DIAG_EULER_CHECK_EDGE");
                string targetCellStr = System.Environment.GetEnvironmentVariable("DIAG_EULER_CHECK_CELL");
                bool cellMatch = targetCellStr == null || (int.TryParse(targetCellStr, out int targetCell0) && inp.jCellIn == targetCell0);
                if (cellMatch && (targetStr == null || (int.TryParse(targetStr, out int targetEdge) && inp.iEdge == targetEdge))) {
                    System.Threading.Interlocked.Increment(ref s_eulerCheckCount);
                    bool probeIsZero = U_IN.All(v => v == 0.0);
                    Console.WriteLine(string.Format(
                        "DIAG_EULER_BND: form={0} edge={1} jCellIn={2} EdgeTag={3} X=({4:0.########},{5:0.########},{6:0.########}) N=({7:0.########},{8:0.########},{9:0.########}) probeZero={10} U_IN=[{11}] vIN={12:0.####} f0={13:0.##########} ret={14:0.##########} ret+f0={15:0.####E-00}",
                        m_EdgForm.GetType().ToString(), inp.iEdge, inp.jCellIn, inp.EdgeTag,
                        inp.X[0], inp.X[1], inp.X[2],
                        inp.Normal[0], inp.Normal[1], inp.Normal[2],
                        probeIsZero,
                        string.Join(",", U_IN.Select(v => v.ToString("0.####"))),
                        _vIN, f0, ret, ret + f0));
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
                throw new ArithmeticException("Got INF while differentiation of boundary edge part of " + m_EdgForm.GetType().Name);
            if (double.IsNaN(dU_iVar))
                throw new ArithmeticException("Got NAN while differentiation of boundary edge part of " + m_EdgForm.GetType().Name);

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
    /// Differentiation of an edge form, used e.g.to obtain a Jacobian of an operator, see <see cref="DifferentialOperator.GetJacobiOperator"/>.
    /// In principal the same as <see cref="EdgeFormDifferentiator"/>, but using <see cref="ILevelSetForm"/>
    /// </summary>
    public class LevelSetFormDifferentiator : FormDifferentiatorCommon, ILevelSetForm, ILevelSetEquationComponentCoefficient {

        ILevelSetForm m_LvsetForm;

        protected override FormDifferentiatorCommon Clone() {
            var LvsetFormClone = (m_LvsetForm as IMultitreadSafety)?.CloneForThread() as ILevelSetForm;
            if (LvsetFormClone != null)
                return new LevelSetFormDifferentiator(LvsetFormClone, m_SpatialDimension);
            else
                return null;
        }

        /// <summary>
        /// ctor
        /// </summary>
        public LevelSetFormDifferentiator(ILevelSetForm ef, int SpatialDimension) :
            base(ef, ef.LevelSetTerms, SpatialDimension) //
        {
            m_LvsetForm = ef;
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags LevelSetTerms {
            get {
                return base.Terms;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public int LevelSetIndex {
            get {
                return m_LvsetForm.LevelSetIndex;
            }
        }
        /// <summary>
        /// %
        /// </summary>
        public string PositiveSpecies {
            get {
                return m_LvsetForm.PositiveSpecies;
            }
        }
        /// <summary>
        /// %
        /// </summary>
        public string NegativeSpecies {
            get {
                return m_LvsetForm.NegativeSpecies;
            }
        }


        //static public int SetDir = 0;
        static public int Dir = 0;


        public double InnerEdgeForm(ref CommonParams inp, double[] U_IN, double[] U_OT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double ret = 0.0;
            int GAMMA = m_LvsetForm.ArgumentOrdering.Count;
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

            double f0 = m_LvsetForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);

            ret += f0; // affine contribution - contains V and GradV contribution

            for (int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (((m_LvsetForm.LevelSetTerms & (TermActivationFlags.UxV | TermActivationFlags.UxGradV)) != 0)) {
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

                if (((m_LvsetForm.LevelSetTerms & (TermActivationFlags.GradUxV | TermActivationFlags.GradUxGradV)) != 0)) {
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

        /// <summary>
        /// passes the coefficients to original form
        /// </summary>
        override public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            if (m_LvsetForm is ILevelSetEquationComponentCoefficient eqc) {
                eqc.CoefficientUpdate(csA, csB, DomainDGdeg, TestDGdeg);
            }
            base.CoefficientUpdate(csA,csB, DomainDGdeg, TestDGdeg);
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
            double f1 = m_LvsetForm.InnerEdgeForm(ref clonedParams, U_IN_temp, U_OT_temp, GradU_IN_temp, GradU_OT_temp, _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
#if DEBUG
            if (double.IsInfinity(f1))
                throw new ArithmeticException();
            if (double.IsNaN(f1))
                throw new ArithmeticException();
#endif

            // compute finite difference
            double dU_iVar = (f1 - f0) / delta;
            if (double.IsInfinity(dU_iVar))
                throw new ArithmeticException("Got INF while differentiation of inner edge part of " + m_LvsetForm.GetType().Name);
            if (double.IsNaN(dU_iVar))
                throw new ArithmeticException("Got NAN while differentiation of inner edge part of " + m_LvsetForm.GetType().Name);

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

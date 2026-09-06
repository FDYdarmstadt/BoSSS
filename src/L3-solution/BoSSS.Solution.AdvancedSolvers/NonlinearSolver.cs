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
using BoSSS.Foundation;
using ilPSP.LinSolvers;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using ilPSP;
using ilPSP.Utils;
using System.Diagnostics;
using System.CodeDom;
using BoSSS.Solution.Tecplot;
using ilPSP.Tracing;
using Microsoft.CodeAnalysis.CSharp.Syntax;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// Evaluation or linearization/matrix assembly of the operator;
    /// this delegate is, so-to-say, the connection between the used <see cref="IDifferentialOperator"/> and its evaluation/linearization,
    /// which can be used to build a <see cref="MultigridOperator"/>
    /// </summary>
    /// <param name="Matrix"></param>
    /// <param name="Affine"></param>
    /// <param name="Linearization">
    /// - false: operator evaluation
    /// - true: linearization
    /// </param>
    /// <param name="MassMatrix"></param>
    /// <param name="CurrentState">
    /// Current state of the solution
    /// </param>
    /// <param name="OberFrickelHack">
    /// the original operator that somehow produced the matrix; yes, this API is convoluted piece-of-shit
    /// </param>
    /// <remarks>
    /// As a recipe, this function does:
    /// 1. compute linearization/evaluation at the current state , i.e. 
    /// 1.1 in the case of `Linearization == true`,
    ///     compute the operator matrix (<see cref="IDifferentialOperator.LinearizationHint"/>, <see cref="IDifferentialOperator.GetMatrixBuilder"/>, <see cref="IDifferentialOperator.GetJacobiOperator"/>, <see cref="IDifferentialOperator.GetFDJacobianBuilder"/>)
    /// 1.2 in the case of `Linearization == false`,
    ///     evaluate the operator (<see cref="IDifferentialOperator.GetEvaluatorEx"/>)
    /// 2. add timestepping terms (<see cref="IDifferentialOperator.TemporalOperator"/>), also depending on the actual timestepping scheme (e.g. BDF or Runge-Kutta)
    /// 3. Compute the mass matrix
    /// 4. perform the agglomeration (<see cref="LevelSetTracker.GetAgglomerator"/>
    /// </remarks>
    public delegate void OperatorEvalOrLin(out BlockMsrMatrix Matrix, out double[] Affine, out BlockMsrMatrix MassMatrix, DGField[] CurrentState, bool Linearization, out IDifferentialOperator OberFrickelHack);


    /// <summary>
    /// base-class for non-linear solvers
    /// </summary>
    public abstract class NonlinearSolver {

        /// <summary>
        /// ctor
        /// </summary>
        public NonlinearSolver(OperatorEvalOrLin __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq, MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig) {
            m_AssembleMatrix = __AssembleMatrix;
            m_AggBasisSeq = __AggBasisSeq.ToArray();
            m_MultigridOperatorConfig = __MultigridOperatorConfig;
        }

        /// <summary>
        /// Evaluation and linearization of PDE to solve
        /// </summary>
        protected OperatorEvalOrLin m_AssembleMatrix;


        IDifferentialOperator m_AbstractOperator;

        /// <summary>
        /// spatial operator provided by <see cref="m_AssembleMatrix"/>;
        /// </summary>
        protected IDifferentialOperator AbstractOperator {
            get {
                return m_AbstractOperator;
            }
            set {
                if (m_AbstractOperator == null)
                    m_AbstractOperator = value;
                if (!object.ReferenceEquals(m_AbstractOperator, value)) {
                    throw new ApplicationException("Hey buddy, you should not change the operator during the solver life-time");
                }
            }
        }


        /// <summary>
        /// Multigrid basis
        /// - 1st index: Multigrid level
        /// - 2nd index: variable index
        /// </summary>
        protected AggregationGridBasis[][] m_AggBasisSeq;

        /// <summary>
        /// required for construction of <see cref="CurrentLin"/>
        /// </summary>
        protected MultigridOperator.ChangeOfBasisConfig[][] m_MultigridOperatorConfig;

        /// <summary>
        /// Called at every iteration; the arguments are 
        ///  - iteration index 
        ///  - current solution 
        ///  - current residual 
        ///  - current multigrid operator
        /// </summary>
        public event Action<int, double[], double[], MultigridOperator> IterationCallback;

        /// <summary>
        /// Triggers <see cref="IterationCallback"/>.
        /// </summary>
        protected void OnIterationCallback(int iterIndex, double[] currentSol, double[] currentRes, MultigridOperator Mgop) {
            if (IterationCallback != null)
                IterationCallback(iterIndex, currentSol, currentRes, Mgop);
        }



        /// <summary>
        /// Helper routine for the initial phase of <see cref="SolverDriver{S}"/>
        /// </summary>
        /// <param name="X">initial guess</param>
        /// <param name="RHS"></param>
        /// <param name="Sol1">
        /// The initial solution, transformed to the aggregation multigrid basis, see <see cref="m_AggBasisSeq"/>.
        /// </param>
        /// <param name="Res1">
        /// The residual of the initial solution
        /// </param>
        protected void Init<S>(CoordinateVector X, S RHS, out double[] Sol1, out double[] Res1)
            where S : IList<double> //
        {
            this.ProblemMapping = X.Mapping;

            int Lraw = X.Mapping.LocalLength;  // length of Solution/RHS in original space
            if (RHS != null) {
                if (RHS.Count != Lraw)
                    throw new ArgumentException();
                this.RHSRaw = RHS.ToArray();
            } else {
                this.RHSRaw = null;
            }
            this.UpdateLinearization(X.Mapping.Fields, 1.0);

            int Ltrf = this.CurrentLin.Mapping.LocalLength;


            // set initial guess (input) as first approximation to the solution ...
            Sol1 = new double[Ltrf];
            this.CurrentLin.TransformSolInto(X, Sol1);

            // ... and evaluate its residual
            Res1 = new double[Ltrf];
            this.EvalLinearizedResidual(Sol1, ref Res1);
        }


        /// <summary>
        /// Template for implementation of the solver routine.
        /// </summary>
        /// <typeparam name="S"></typeparam>
        /// <param name="X">
        /// On entry, an initial guess to the linear system.
        /// On exit, hopefully the solution to the nonlinear system.
        /// </param>
        /// <param name="RHS">
        /// If not equal null, must be passed to <see cref="RHSRaw"/>.
        /// </param>
        /// <returns>
        /// - true: solver algorithm successfully converged
        /// - false: something went wrong
        /// </returns>
        abstract public bool SolverDriver<S>(CoordinateVector X, S RHS)
           where S : IList<double>;

        /// <summary>
        /// Number of nonlinear iterations in last call to <see cref="SolverDriver{S}(CoordinateVector, S)"/>
        /// </summary>
        public int NoOfNonlinearIter {
            get;
            protected set;
        }


        /// <summary>
        /// Preconditioner/solver configuration for the linearized problem
        /// </summary>
        public ISolverFactory PrecondConfig;





        /// <summary>
        /// Current linearization of the nonlinear operator: the linearized
        /// system is given as 
        /// <see cref="CurrentLin"/>*X = <see cref="LinearizationRHS"/>.
        /// </summary>
        protected MultigridOperator CurrentLin;

        /// <summary>
        /// number of DOFs in linearization of last iteration
        /// </summary>
        public long EssentialDOFs {
            get {
                return CurrentLin?.Mapping.TotalLength ?? 0;
            }
        }


        /// <summary>
        /// Optional RHS to the nonlinear system, 
        /// </summary>
        protected double[] RHSRaw;

        /// <summary>
        /// Current affine part of the nonlinear operator.
        /// </summary>
        protected double[] LinearizationRHS;



        /// <summary>
        /// Evaluation of the nonlinear operator.
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="CurrentState">
        /// Current state of DG fields
        /// </param>
        /// <param name="Output"></param>
        /// <param name="HomotopyValue">
        /// <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
        /// </param>
        /// <param name="ApplyRef">
        /// apply additional modification due to free-mean-value fixing (aka. pressure reference point), <see cref="MultigridOperator.FreeMeanValue"/>
        /// </param>        
        protected void EvaluateOperator(double alpha, IEnumerable<DGField> CurrentState, double[] Output, double HomotopyValue, bool ApplyRef = false) {
            using (var tr = new FuncTrace()) {
                if (alpha != 1.0)
                    throw new NotSupportedException("some moron has removed this");

                SetHomotopyValue(HomotopyValue);

                // the real call:
                // DIAGNOSTIC (coordinator-requested bounded check, this round): tag which of the
                // two assembly paths (F-path raw eval vs M-path linearization, below) is
                // currently in flight via an env var, so deeply-nested equation-component code
                // (e.g. LinearizedConvection.cs's InnerEdgeFlux_impl, called identically by both)
                // can label its own DIAG_FLUX_RAW/DIAG_FLUX_LAMBDA prints accordingly. An env var
                // (not a cross-assembly static field) is used deliberately: NSECommon does not
                // reference AdvancedSolvers. No effect on solver behavior otherwise.
                System.Environment.SetEnvironmentVariable("DIAG_CURRENT_ASSEMBLY_PATH", $"F-PATH-RAWEVAL@{EvaluationCounter}");
                // PERFORMANCE-INVESTIGATION DIAGNOSTIC (coordinator-requested, this round,
                // temporary): wall-clock timing for this real per-iteration F-path (residual)
                // assembly call, gated by DIAG_ASSEMBLY_TIMING=1. Zero extra cost when unset.
                bool diagAsmTiming = System.Environment.GetEnvironmentVariable("DIAG_ASSEMBLY_TIMING") == "1";
                var diagAsmSwF = diagAsmTiming ? System.Diagnostics.Stopwatch.StartNew() : null;
                this.m_AssembleMatrix(out BlockMsrMatrix DummyMtx, out double[] OpEvalRaw, out BlockMsrMatrix MassMtxRaw, CurrentState.ToArray(), false, out var abstractOp);
                if (diagAsmTiming) {
                    diagAsmSwF.Stop();
                    Console.Error.WriteLine($"DIAG_ASSEMBLY_TIMING: tag=F-PATH evalCount={EvaluationCounter} elapsedMs={diagAsmSwF.Elapsed.TotalMilliseconds:0.###}");
                }
                System.Environment.SetEnvironmentVariable("DIAG_CURRENT_ASSEMBLY_PATH", "unset");
                if (DummyMtx != null)
                    // only evaluation ==> OpMatrix must be null
                    throw new ApplicationException($"The provided {typeof(OperatorEvalOrLin).Name} is not correctly implemented.");
                this.AbstractOperator = abstractOp;

               
#if DEBUG
                const int TEST_INTERVALL = 10;
#else
                const int TEST_INTERVALL = 1000;
#endif
                if (EvaluationCounter % TEST_INTERVALL == 0) { // do the following, expensive check only for every TEST_INTERVALL-th evaluation.

                    // Comparison of linearization and evaluation:
                    // -------------------------------------------
                    //
                    // Note that, in BoSSS, currently the Linearization of f(u) around u0 is defines as
                    //     f(u) ≈ M(u0)*u + b(u0),
                    // instead of the typical Taylor series representation f(u) ≈ f(u0) + ∂f(u0)*(u-u0).
                    // The relation between the BoSSS-representation and the Taylor series is
                    //     M(u0) = ∂f(u0),
                    //     b(u0) = f(u0) - ∂f(u0)*u0.
                    // Therefore, we check that
                    //     M(u0)*u0 + b(u0) = f(u0).
                    //
                
                    System.Environment.SetEnvironmentVariable("DIAG_CURRENT_ASSEMBLY_PATH", $"M-PATH-JACOBIAN@{EvaluationCounter}");
                    this.m_AssembleMatrix(out BlockMsrMatrix LinMtx, out double[] OpAffine, out _, CurrentState.ToArray(), true, out _);
                    System.Environment.SetEnvironmentVariable("DIAG_CURRENT_ASSEMBLY_PATH", "unset");
                    var Check = OpAffine.CloneAs();
                    LinMtx.SpMV(1.0, new CoordinateVector(CurrentState), 1.0, Check);

                    var err = Check.CloneAs();
                    err.AccV(-1.0, OpEvalRaw);
                    double l2_err = err.MPI_L2Norm();
                    double comp = Math.Sqrt(Math.Max(OpEvalRaw.MPI_L2Norm(), Check.MPI_L2Norm()) * BLAS.MachineEps + BLAS.MachineEps);

                    // DIAGNOSTIC (coordinator-requested mesh-refinement/dilution investigation):
                    // dump the worst per-DOF Check/OpEvalRaw/err for every codomain variable,
                    // RESTRICTED to DOFs on periodic-boundary-touching cells, UNCONDITIONALLY at
                    // every TEST_INTERVALL-th evaluation -- i.e. regardless of whether the
                    // aggregate "Mismatch" check below actually trips. This exists specifically to
                    // separate two hypotheses for why finer production grids stop printing the
                    // "Mismatch" warning: (a) the absolute per-DOF discrepancy at a periodic cell
                    // genuinely shrinks with mesh refinement, or (b) it stays the same absolute
                    // size while `comp`/the overall norm grows, diluting it below the l2_err>comp
                    // aggregate threshold. Reading the raw per-DOF numbers directly, at every
                    // resolution, sidesteps needing the aggregate check to fire at all.
                    if (System.Environment.GetEnvironmentVariable("DIAG_PERIODIC_SCAN") == "1") {
                        try {
                            var worstByField = new Dictionary<string, (int idx, double abs, int jCell, int n)>();
                            for (int i = 0; i < err.Length; i++) {
                                this.ProblemMapping.LocalFieldCoordinateIndex(i, out int ifld, out int jCell, out int n);

                                BoSSS.Foundation.Grid.IGridData gridDat = (ifld >= 0 && ifld < this.ProblemMapping.BasisS.Count)
                                    ? this.ProblemMapping.BasisS[ifld].GridDat
                                    : null;
                                bool isPeriodicCell = false;
                                if (gridDat != null) {
                                    var E2C = gridDat.iGeomEdges.CellIndices;
                                    var tags = gridDat.iGeomEdges.EdgeTags;
                                    int nEdges = tags.Length;
                                    for (int e = 0; e < nEdges; e++) {
                                        if ((E2C[e, 0] == jCell || E2C[e, 1] == jCell)
                                            && tags[e] >= BoSSS.Foundation.Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG) {
                                            isPeriodicCell = true;
                                            break;
                                        }
                                    }
                                }
                                if (!isPeriodicCell)
                                    continue;

                                string eqName = (this.AbstractOperator != null
                                    && this.AbstractOperator.CodomainVar != null
                                    && ifld >= 0 && ifld < this.AbstractOperator.CodomainVar.Count)
                                    ? this.AbstractOperator.CodomainVar[ifld]
                                    : $"field#{ifld}";

                                double a = Math.Abs(err[i]);
                                if (!worstByField.TryGetValue(eqName, out var cur) || a > cur.abs) {
                                    worstByField[eqName] = (i, a, jCell, n);
                                }
                            }

                            Console.Error.WriteLine(
                                $"DIAG_PERIODIC_SCAN: evalCount={EvaluationCounter} l2_err={l2_err:G17} comp={comp:G17} tripped={(l2_err > comp)} " +
                                $"||OpEvalRaw||={OpEvalRaw.MPI_L2Norm():G17} ||Check||={Check.MPI_L2Norm():G17} ||err||={l2_err:G17}");
                            foreach (var kv in worstByField) {
                                string eqName = kv.Key;
                                var (idx, absv, jCell, n) = kv.Value;
                                Console.Error.WriteLine(
                                    $"DIAG_PERIODIC_SCAN:   field={eqName} worstPeriodicIdx={idx} jCell={jCell} node={n} " +
                                    $"Check={Check[idx]:G17} OpEvalRaw={OpEvalRaw[idx]:G17} err={err[idx]:G17} |err|={absv:G17}");
                                // DIAGNOSTIC (coordinator-requested cheap-reproduction-scan, this
                                // round): also report THIS field's own worst-periodic-DOF row nnz
                                // -- DIAG_MROW above only ever dumps the GLOBAL worst DOF (which,
                                // at reduced DG degree, may belong to a different field), so this
                                // gives a per-field nnz reading independent of which field happens
                                // to be globally worst.
                                if (System.Environment.GetEnvironmentVariable("DIAG_PERIODIC_MROW") == "1") {
                                    try {
                                        long[] pcols = null;
                                        double[] pvals = null;
                                        int pnnz = LinMtx.GetRow(idx, ref pcols, ref pvals);
                                        int realNnz = 0;
                                        for (int k = 0; k < pnnz; k++) if (pvals[k] != 0.0) realNnz++;
                                        Console.Error.WriteLine(
                                            $"DIAG_PERIODIC_SCAN:   -> field={eqName} ownWorstRow nnz={realNnz} (raw entries={pnnz})");
                                    } catch (Exception mrowEx3) {
                                        Console.Error.WriteLine($"DIAG_PERIODIC_SCAN:   -> field={eqName} ownWorstRow nnz FAILED: {mrowEx3.Message}");
                                    }
                                }
                            }
                            // DIAGNOSTIC (coordinator-requested degree-1 HeatEq reconstruction,
                            // this round): full-row dump for an EXPLICIT target index -- DIAG_MROW
                            // above only ever fires for the GLOBAL worst DOF (at reduced degree
                            // this is MomentumY, not HeatEq), so this lets us get HeatEq's own
                            // worst row's full column list even when it's not globally worst.
                            // Piggybacks on this block's existing every-10th-call cadence, so it's
                            // inherently bounded/cheap -- no separate cap needed.
                            string targetRowStr = System.Environment.GetEnvironmentVariable("DIAG_TARGET_MROW");
                            if (targetRowStr != null && int.TryParse(targetRowStr, out int targetRowIdx)) {
                                try {
                                    long[] tcols = null;
                                    double[] tvals = null;
                                    int tnnz = LinMtx.GetRow(targetRowIdx, ref tcols, ref tvals);
                                    var u0vec = new CoordinateVector(CurrentState);
                                    var sb2 = new System.Text.StringBuilder();
                                    sb2.Append($"DIAG_TARGET_MROW: targetIdx={targetRowIdx} evalCount={EvaluationCounter} Check={Check[targetRowIdx]:G17} OpEvalRaw={OpEvalRaw[targetRowIdx]:G17} OpAffine={OpAffine[targetRowIdx]:G17} rawEntries={tnnz}\n");
                                    for (int k = 0; k < tnnz; k++) {
                                        long c = tcols[k];
                                        double v = tvals[k];
                                        if (v == 0.0) continue;
                                        string colDesc = "?";
                                        double u0val = double.NaN;
                                        try {
                                            this.ProblemMapping.LocalFieldCoordinateIndex((int)c, out int cifld, out int cjCell, out int cn);
                                            string cName = (this.AbstractOperator != null && this.AbstractOperator.DomainVar != null
                                                && cifld >= 0 && cifld < this.AbstractOperator.DomainVar.Count)
                                                ? this.AbstractOperator.DomainVar[cifld] : $"field#{cifld}";
                                            colDesc = $"{cName} cell={cjCell} mode={cn}";
                                        } catch { }
                                        try { u0val = u0vec[(int)c]; } catch { }
                                        sb2.Append($"       col={c} [{colDesc}] M={v:G17} u0={u0val:G17}\n");
                                    }
                                    Console.Error.Write(sb2.ToString());

                                    // DIAGNOSTIC (coordinator-requested degree-1 HeatEq
                                    // reconstruction, this round): frozen-state ablation for the
                                    // TARGET row specifically -- same "assemble twice, ablate one
                                    // env var, restore, no re-solve" pattern as DIAG_FIXED_ABLATION
                                    // above, but applied to targetRowIdx's own M-row instead of
                                    // the global worstIdx. This isolates a term's contribution to
                                    // targetRowIdx's columns AT THE EXACT SAME u0 as the unablated
                                    // dump above, avoiding the "two separate re-solves converge to
                                    // slightly different u0" noise that a real (re-solved) ablation
                                    // comparison would introduce for u0-dependent (nonlinear) terms.
                                    string ablVar = System.Environment.GetEnvironmentVariable("DIAG_TARGET_MROW_ABLATE_VAR");
                                    string ablVal = System.Environment.GetEnvironmentVariable("DIAG_TARGET_MROW_ABLATE_VAL");
                                    if (ablVar != null && ablVal != null) {
                                        string savedAblVal = System.Environment.GetEnvironmentVariable(ablVar);
                                        try {
                                            System.Environment.SetEnvironmentVariable(ablVar, ablVal);
                                            this.m_AssembleMatrix(out BlockMsrMatrix LinMtx_ta, out double[] OpAffine_ta, out _, CurrentState.ToArray(), true, out _);
                                            long[] ta_cols = null;
                                            double[] ta_vals = null;
                                            int ta_nnz = LinMtx_ta.GetRow(targetRowIdx, ref ta_cols, ref ta_vals);
                                            var sb3 = new System.Text.StringBuilder();
                                            sb3.Append($"DIAG_TARGET_MROW_ABLATED: targetIdx={targetRowIdx} evalCount={EvaluationCounter} ablVar={ablVar} ablVal={ablVal} OpAffine_abl={OpAffine_ta[targetRowIdx]:G17} rawEntries={ta_nnz}\n");
                                            for (int k = 0; k < ta_nnz; k++) {
                                                long c = ta_cols[k];
                                                double v = ta_vals[k];
                                                if (v == 0.0) continue;
                                                string colDesc = "?";
                                                try {
                                                    this.ProblemMapping.LocalFieldCoordinateIndex((int)c, out int cifld, out int cjCell, out int cn);
                                                    string cName = (this.AbstractOperator != null && this.AbstractOperator.DomainVar != null
                                                        && cifld >= 0 && cifld < this.AbstractOperator.DomainVar.Count)
                                                        ? this.AbstractOperator.DomainVar[cifld] : $"field#{cifld}";
                                                    colDesc = $"{cName} cell={cjCell} mode={cn}";
                                                } catch { }
                                                sb3.Append($"       col={c} [{colDesc}] M_abl={v:G17}\n");
                                            }
                                            Console.Error.Write(sb3.ToString());
                                        } catch (Exception ablEx) {
                                            Console.Error.WriteLine($"DIAG_TARGET_MROW_ABLATED: (threw) {ablEx.Message}");
                                        } finally {
                                            System.Environment.SetEnvironmentVariable(ablVar, savedAblVal);
                                        }
                                    }
                                } catch (Exception tgtEx) {
                                    Console.Error.WriteLine($"DIAG_TARGET_MROW: (threw) {tgtEx.Message}");
                                }
                            }
                        } catch (Exception diagEx2) {
                            Console.Error.WriteLine($"    -> (DIAG_PERIODIC_SCAN itself failed, ignoring: {diagEx2.Message})");
                        }
                    }

                    if (l2_err > comp) {
                        //tr.Error($"EvaluateOperator: Mismatch between operator linearization and evaluation: Operator matrix-Jacobian distance: {l2_err}, relative: {l2_err / comp} (comparison value {comp})");
                        Console.Error.WriteLine($"Mismatch between operator linearization and evaluation: absolute dist.: {l2_err}, relative: {l2_err/comp} (comparison denominator {comp})");
                        //throw new ArithmeticException($"Mismatch between operator linearization and evaluation: Operator matrix-Jacobian distance: { l2_err }, relative: { l2_err/comp} (comparison value { comp})");

                        // DIAGNOSTIC-ONLY (never throws, never alters solver state): identify
                        // WHERE (on this MPI rank) the worst contribution to the mismatch above
                        // lives -- which equation/variable, which cell, and whether that cell
                        // touches a periodic boundary and/or is an XDG cut-cell. Added to
                        // distinguish "this is the known periodic-rotation convection-form defect"
                        // from "this is something else entirely" (e.g. XDG cut-cell mechanics,
                        // wall BCs) on real production runs. Wrapped in try/catch so a failure in
                        // this reporting code can never destabilize the actual solve.
                        try {
                            int worstIdx = -1;
                            double worstAbs = -1;
                            for (int i = 0; i < err.Length; i++) {
                                double a = Math.Abs(err[i]);
                                if (a > worstAbs) {
                                    worstAbs = a;
                                    worstIdx = i;
                                }
                            }
                            if (worstIdx >= 0) {
                                this.ProblemMapping.LocalFieldCoordinateIndex(worstIdx, out int ifld, out int jCell, out int n);

                                string eqName = (this.AbstractOperator != null
                                    && this.AbstractOperator.CodomainVar != null
                                    && ifld >= 0 && ifld < this.AbstractOperator.CodomainVar.Count)
                                    ? this.AbstractOperator.CodomainVar[ifld]
                                    : $"field#{ifld}";

                                BoSSS.Foundation.Grid.IGridData gridDat = (ifld >= 0 && ifld < this.ProblemMapping.BasisS.Count)
                                    ? this.ProblemMapping.BasisS[ifld].GridDat
                                    : null;

                                bool isPeriodicCell = false;
                                if (gridDat != null) {
                                    var E2C = gridDat.iGeomEdges.CellIndices;
                                    var tags = gridDat.iGeomEdges.EdgeTags;
                                    int nEdges = tags.Length;
                                    for (int e = 0; e < nEdges; e++) {
                                        if ((E2C[e, 0] == jCell || E2C[e, 1] == jCell)
                                            && tags[e] >= BoSSS.Foundation.Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG) {
                                            isPeriodicCell = true;
                                            break;
                                        }
                                    }
                                }

                                bool isXdg = ifld >= 0 && ifld < this.ProblemMapping.BasisS.Count
                                    && this.ProblemMapping.BasisS[ifld] is XDGBasis;
                                bool? isCutCell = null;
                                if (isXdg) {
                                    var xb = (XDGBasis)this.ProblemMapping.BasisS[ifld];
                                    if (xb.Tracker != null)
                                        isCutCell = xb.Tracker.Regions.GetCutCellMask().Contains(jCell);
                                }

                                Console.Error.WriteLine(
                                    $"    -> worst mismatch on MPI rank {ilPSP.Environment.MPIEnv.MPI_Rank}: |err|={worstAbs:0.####E-00} at equation '{eqName}', local cell={jCell}, DG node={n}; "
                                    + $"touchesPeriodicBoundary={isPeriodicCell}; isXdgField={isXdg}"
                                    + (isXdg ? $", isCutCell={(isCutCell.HasValue ? isCutCell.Value.ToString() : "unknown (no Tracker)")}" : ""));

                                // DIAGNOSTIC (coordinator-requested full from-scratch DOF
                                // reconstruction, this round): print the exact two assembled
                                // values being compared (Check[worstIdx] = M(u0)*u0+b(u0), and
                                // OpEvalRaw[worstIdx] = f(u0)) at FULL precision, so an
                                // independently-computed from-scratch summation can be diffed
                                // directly against BOTH of the real code's own values, not just
                                // their difference.
                                Console.Error.WriteLine(
                                    $"    -> DIAG_DOFVALUE: worstIdx={worstIdx} Check[worstIdx]={Check[worstIdx]:G17} OpEvalRaw[worstIdx]={OpEvalRaw[worstIdx]:G17} err={err[worstIdx]:G17}");

                                // DIAGNOSTIC (coordinator-requested independent Jacobian-row
                                // reconstruction, this round): dump the FULL nonzero row of
                                // LinMtx (M) for worstIdx -- every occupied column index, its
                                // decoded (field,cell,mode), and its actual matrix value -- so an
                                // independent finite-difference reconstruction can be compared
                                // column-by-column against the real code's own M entries, rather
                                // than guessing which columns matter.
                                try {
                                    long[] cols = null;
                                    double[] vals = null;
                                    int nnz = LinMtx.GetRow(worstIdx, ref cols, ref vals);
                                    var sb = new System.Text.StringBuilder();
                                    sb.Append($"    -> DIAG_MROW: worstIdx={worstIdx} nnz={nnz}\n");
                                    for (int k = 0; k < nnz; k++) {
                                        long c = cols[k];
                                        double v = vals[k];
                                        if (v == 0.0) continue;
                                        string colDesc = "?";
                                        try {
                                            this.ProblemMapping.LocalFieldCoordinateIndex((int)c, out int cifld, out int cjCell, out int cn);
                                            string cName = (this.AbstractOperator != null && this.AbstractOperator.DomainVar != null
                                                && cifld >= 0 && cifld < this.AbstractOperator.DomainVar.Count)
                                                ? this.AbstractOperator.DomainVar[cifld] : $"field#{cifld}";
                                            colDesc = $"{cName} cell={cjCell} mode={cn}";
                                        } catch { }
                                        sb.Append($"       col={c} [{colDesc}] M={v:G17}\n");
                                    }
                                    Console.Error.Write(sb.ToString());
                                } catch (Exception mrowEx) {
                                    Console.Error.WriteLine($"    -> DIAG_MROW: (threw) {mrowEx.Message}");
                                }

                        // DIAGNOSTIC (coordinator-requested disambiguation, this round):
                        // evaluate a "counterfactual ablated" Check/OpEvalRaw at this EXACT
                        // SAME CurrentState (u0) -- NOT by re-solving, just two extra one-shot
                        // evaluations with a target edge's conduction (or other) contribution
                        // temporarily toggled off via the SAME env-var-gated ablation
                        // mechanism used for the real (re-solved) ablation tests, then
                        // restored immediately after. This isolates whether an ablation
                        // finding reflects a genuine defect AT THIS baseline state, vs. an
                        // artifact of re-solving to a different converged state.
                        if (System.Environment.GetEnvironmentVariable("DIAG_FIXED_ABLATION") == "1") {
                            string targetEdgeVar = System.Environment.GetEnvironmentVariable("DIAG_FIXED_ABLATION_EDGEVAR") ?? "XNSFE_REPRO_DISABLE_EDGE_COND";
                            string targetEdgeVal = System.Environment.GetEnvironmentVariable("DIAG_FIXED_ABLATION_EDGEVAL") ?? "35";
                            string savedVal = System.Environment.GetEnvironmentVariable(targetEdgeVar);
                            try {
                                System.Environment.SetEnvironmentVariable(targetEdgeVar, targetEdgeVal);

                                this.m_AssembleMatrix(out BlockMsrMatrix DummyMtx2, out double[] OpEvalRaw_abl, out BlockMsrMatrix _, CurrentState.ToArray(), false, out var _abstractOp2);
                                this.m_AssembleMatrix(out BlockMsrMatrix LinMtx_abl, out double[] OpAffine_abl, out _, CurrentState.ToArray(), true, out _);
                                var Check_abl = OpAffine_abl.CloneAs();
                                LinMtx_abl.SpMV(1.0, new CoordinateVector(CurrentState), 1.0, Check_abl);

                                var err_abl = Check_abl.CloneAs();
                                err_abl.AccV(-1.0, OpEvalRaw_abl);

                                int worstIdx_abl = -1;
                                double worstAbs_abl = -1;
                                for (int i = 0; i < err_abl.Length; i++) {
                                    double a = Math.Abs(err_abl[i]);
                                    if (a > worstAbs_abl) { worstAbs_abl = a; worstIdx_abl = i; }
                                }
                                string worstDesc_abl = "?";
                                try {
                                    this.ProblemMapping.LocalFieldCoordinateIndex(worstIdx_abl, out int aifld, out int ajCell, out int an);
                                    string aName = (this.AbstractOperator != null && this.AbstractOperator.CodomainVar != null
                                        && aifld >= 0 && aifld < this.AbstractOperator.CodomainVar.Count)
                                        ? this.AbstractOperator.CodomainVar[aifld] : $"field#{aifld}";
                                    worstDesc_abl = $"{aName} cell={ajCell} mode={an}";
                                } catch { }

                                Console.Error.WriteLine(
                                    $"    -> DIAG_FIXED_ABLATION: ({targetEdgeVar}={targetEdgeVal}, SAME u0, NO re-solve) " +
                                    $"at ORIGINAL worstIdx={worstIdx}: Check_abl={Check_abl[worstIdx]:G17} OpEvalRaw_abl={OpEvalRaw_abl[worstIdx]:G17} err_abl={err_abl[worstIdx]:G17}  " +
                                    $"(unablated at same idx/u0: Check={Check[worstIdx]:G17} OpEvalRaw={OpEvalRaw[worstIdx]:G17} err={err[worstIdx]:G17})");
                                Console.Error.WriteLine(
                                    $"    -> DIAG_FIXED_ABLATION: full-array worst under ablation: worstIdx_abl={worstIdx_abl} [{worstDesc_abl}] |err_abl|={worstAbs_abl:G17}  " +
                                    $"(vs unablated full-array worst |err|={worstAbs:G17} at idx={worstIdx})");
                            } catch (Exception fixedAblEx) {
                                Console.Error.WriteLine($"    -> DIAG_FIXED_ABLATION: (threw) {fixedAblEx.Message}");
                            } finally {
                                System.Environment.SetEnvironmentVariable(targetEdgeVar, savedVal);
                            }
                        }
                            }
                        } catch (Exception diagEx) {
                            Console.Error.WriteLine($"    -> (mismatch-location diagnostic itself failed, ignoring: {diagEx.Message})");
                        }
				    } else {
                        //Console.WriteLine($"MATCH between operator linearization and evaluation: absolute dist.:  {l2_err}, relative: {l2_err / comp} (comparison value {comp})");
                    }
                    //*/
                }
                EvaluationCounter++;

                CurrentLin.TransformRhsInto(OpEvalRaw, Output, ApplyRef);
            }
        }

        int EvaluationCounter = 0;


        /// <summary>
        /// Update of the homotopy value 
        /// </summary>
        /// <param name="HomotopyValue">
        /// must be between 0 and 1 (both including), <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
        /// </param>
        protected void SetHomotopyValue(double HomotopyValue) {
            if (HomotopyValue < 0)
                throw new ArgumentOutOfRangeException();
            if (HomotopyValue > 1)
                throw new ArgumentOutOfRangeException();


            if (AbstractOperator == null && HomotopyValue != 1) {
                // do one call just to attain the spatial operator;
                // not very efficient, but ok for the moment

                var dummyState = this.CurrentLin.BaseGridProblemMapping.BasisS.Select(delegate (Basis b) {
                    DGField r;
                    if (b is XDGBasis xb)
                        r = new XDGField(xb);
                    else
                        r = new SinglePhaseField(b);
                    return r;
                }).ToArray();

                this.m_AssembleMatrix(out var OpMtxRaw, out _, out _, dummyState, false, out var _Dummy);
                Debug.Assert(OpMtxRaw == null); // only evaluation ==> OpMatrix must be null
                this.AbstractOperator = _Dummy;
            }

            if (AbstractOperator == null && HomotopyValue != 1.0) {
                throw new NotSupportedException("unable to attain operator for homotopy update");
            }

            if (AbstractOperator != null && HomotopyValue != AbstractOperator.CurrentHomotopyValue)
                // set homotopy value
                AbstractOperator.CurrentHomotopyValue = HomotopyValue;
        }



        /// <summary>
        /// Updating the <see cref="CurrentLin"/> -- operator;
        /// </summary>
        /// <param name="CurrentState">input, linearization point</param>
        /// <param name="U0">output, linearization point, after external update, transformed back</param>
        /// <param name="HomotopyValue">
        /// <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
        /// </param>
        protected void Update(IEnumerable<DGField> CurrentState, double[] U0, double HomotopyValue) {

            this.UpdateLinearization(CurrentState, HomotopyValue);
            CurrentLin.TransformSolInto(new CoordinateVector(CurrentState), U0);
        }


        /// <summary>
        /// Updating the <see cref="CurrentLin"/> -- operator;
        /// </summary>
        /// <param name="CurrentState">linearization point</param>
        /// <param name="HomotopyValue">
        /// <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
        /// </param>
        protected void UpdateLinearization(IEnumerable<DGField> CurrentState, double HomotopyValue) {
            if (!(this.ProblemMapping.BasisS.Count == CurrentState.Count()))
                throw new ArgumentException("mismatch in number of fields.");

            SetHomotopyValue(HomotopyValue);

            // the real call:
            // PERFORMANCE-INVESTIGATION DIAGNOSTIC (coordinator-requested, this round,
            // temporary): wall-clock timing for this real per-iteration M-path (Jacobian
            // linearization) assembly call, gated by DIAG_ASSEMBLY_TIMING=1. Zero extra cost
            // when unset. `LinearizationCounter` is not usable as a call tag here -- it only
            // ever increments inside the pre-existing, permanently-commented-out block below.
            bool diagAsmTimingM = System.Environment.GetEnvironmentVariable("DIAG_ASSEMBLY_TIMING") == "1";
            var diagAsmSwM = diagAsmTimingM ? System.Diagnostics.Stopwatch.StartNew() : null;
            this.m_AssembleMatrix(out BlockMsrMatrix OpMtxRaw, out double[] OpAffineRaw, out BlockMsrMatrix MassMtxRaw, CurrentState.ToArray(),
                true, // `true` signals that we want a linearization
                out IDifferentialOperator abstractOperator);
            if (diagAsmTimingM) {
                diagAsmSwM.Stop();
                s_asmTimingMCount++;
                Console.Error.WriteLine($"DIAG_ASSEMBLY_TIMING: tag=M-PATH call#={s_asmTimingMCount} elapsedMs={diagAsmSwM.Elapsed.TotalMilliseconds:0.###}");
            }
            AbstractOperator = abstractOperator;
            
            
#if DEBUG
            const int TEST_INTERVALL = 10;
#else
            const int TEST_INTERVALL = 1000;
#endif
            if (LinearizationCounter % TEST_INTERVALL == 0) {
                // do the following, expensive check only for every TEST_INTERVALL-th evaluation.
                // Comparison of linearization and evaluation:
                // -------------------------------------------
                //
                // Note that, in BoSSS, currently the Linearization of f(u) around u0 is defines as
                //     f(u) ≈ M(u0)*u + b(u0),
                // instead of the typical Taylor series representation f(u) ≈ f(u0) + ∂f(u0)*(u-u0).
                // The relation between the BoSSS-representation and the Taylor series is
                //     M(u0) = ∂f(u0),
                //     b(u0) = f(u0) - ∂f(u0)*u0.
                // Therefore, we check that
                //     M(u0)*u0 + b(u0) = f(u0).
                //

                /*
                // the evaluation
                this.m_AssembleMatrix(out BlockMsrMatrix DummyMtx, out double[] OpEvalRaw, out _, CurrentState.ToArray(), 
                    false, // `false` signals tat we want an evaluation
                    out _);
                if (DummyMtx != null)
                    // only evaluation ==> OpMatrix must be null
                    throw new ApplicationException($"The provided {typeof(OperatorEvalOrLin).Name} is not correctly implemented.");
                
                var Check = OpAffineRaw.CloneAs();
                OpMtxRaw.SpMV(1.0, new CoordinateVector(CurrentState), 1.0, Check);

                var err = Check.CloneAs();
                err.AccV(-1.0, OpEvalRaw);
                double l2_err = err.MPI_L2Norm();
                double comp = Math.Sqrt(Math.Max(OpEvalRaw.MPI_L2Norm(), Check.MPI_L2Norm()) * BLAS.MachineEps + BLAS.MachineEps);

                if (l2_err > comp) {
                    Console.Error.WriteLine($"UpdateLinearization: Mismatch between operator linearization and evaluation: Operator matrix-Jacobian distance: {l2_err}, relative: {l2_err/comp} (comparison value {comp})");
                    //throw new ArithmeticException($"Mismatch between operator linearization and evaluation: Operator matrix-Jacobian distance: { l2_err }, relative: { l2_err/comp} (comparison value { comp})");
                }
                

                LinearizationCounter++;
                
                {
                    this.m_AssembleMatrix(out BlockMsrMatrix OpMtxRaw2, out double[] OpAffineRaw2, out BlockMsrMatrix _, CurrentState.ToArray(),
                        true, // `true` signals that we want a linearization
                        out _);

                    OpMtxRaw2.Acc(-1.0, OpMtxRaw);
                    OpAffineRaw2.AccV(-1.0, OpAffineRaw);

                    Console.WriteLine($"Linearization change: {OpMtxRaw2.InfNorm()}, {OpAffineRaw2.MPI_L2Norm()}");
                }
                */

            }


            // setup of the multigrid operator
            CurrentLin = new MultigridOperator(this.m_AggBasisSeq, this.ProblemMapping,
                OpMtxRaw.CloneAs(), MassMtxRaw,
                this.m_MultigridOperatorConfig,
                AbstractOperator);
            OpAffineRaw = OpAffineRaw.CloneAs();

            // RHS of the linearization
            if (this.RHSRaw != null)
                OpAffineRaw.AccV(-1.0, this.RHSRaw);
            if (LinearizationRHS == null || LinearizationRHS.Length != this.CurrentLin.Mapping.LocalLength)
                LinearizationRHS = new double[this.CurrentLin.Mapping.LocalLength];
            else
                LinearizationRHS.ClearEntries();
            CurrentLin.TransformRhsInto(OpAffineRaw, this.LinearizationRHS, true);
            this.LinearizationRHS.ScaleV(-1.0);

            //writing to text files
            //{
            //    Console.WriteLine("Writing to text file");
            //    OpMtxRaw.SaveToTextFileSparse("OpMtxRaw-" + LinearizationCounter);
            //    MassMtxRaw.SaveToTextFileSparse("MassMtxRaw-" + LinearizationCounter);
            //    CurrentLin.LeftChangeOfBasis.SaveToTextFileSparse("LeftChangeOfBasis-" + LinearizationCounter);
            //    CurrentLin.RightChangeOfBasis.SaveToTextFileSparse("RightChangeOfBasis-" + LinearizationCounter);
            //    CurrentLin.OperatorMatrix.SaveToTextFileSparse("OperatorMatrix-" + LinearizationCounter);
            //    this.LinearizationRHS.SaveToTextFile("RHS-" + LinearizationCounter);
            //    LinearizationCounter++;
            //}

        }

        int LinearizationCounter = 0;

        // PERFORMANCE-INVESTIGATION DIAGNOSTIC (coordinator-requested, this round, temporary):
        // call counter for DIAG_ASSEMBLY_TIMING's M-PATH tag, see UpdateLinearization.
        static int s_asmTimingMCount = 0;


        /// <summary>
        /// This method tests that, 
        /// if any entry in <see cref="IDifferentialOperator.FreeMeanValue"/> is true, 
        /// a change in the mean/average value of the respective variable must **not** have any effect on the residual.
        /// </summary>
        public bool TestFreeMeanValue(CoordinateVector SolutionVec, double HomotopyValue) {
            using (var ft = new FuncTrace()) {
                ft.InfoToConsole = false;
                ft.Info("Checking for free mean value.");

                int L = this.CurrentLin.Mapping.LocalLength;

                if (CurrentLin.FreeMeanValue.Any()) {

                    double[] ResidualBeforMeanCor = new double[L];
                    double[] ResidualAfterMeanCor = new double[L];

                    EvaluateOperator(1, SolutionVec.Mapping.Fields, ResidualBeforMeanCor, HomotopyValue);
                    double SolNormA = SolutionVec.MPI_L2Norm();


                    DGField[] flds = SolutionVec.Mapping.Fields.ToArray();
                    bool[] FreeMeanValue = CurrentLin.FreeMeanValue;
                    if (flds.Length != FreeMeanValue.Length)
                        throw new ApplicationException();

                    const double arbitrary_distortion_value = 2000.1234;

                    for (int iFld = 0; iFld < flds.Length; iFld++) {
                        if (FreeMeanValue[iFld]) {
                            flds[iFld].AccConstant(arbitrary_distortion_value);
                        }
                    }

                    EvaluateOperator(1, SolutionVec.Mapping.Fields, ResidualAfterMeanCor, HomotopyValue);
                    double SolNormB = SolutionVec.L2Norm();

                    for (int iFld = 0; iFld < flds.Length; iFld++) {
                        if (FreeMeanValue[iFld]) {
                            flds[iFld].AccConstant(-arbitrary_distortion_value);
                        }
                    }

                    double Dist = ResidualBeforMeanCor.MPI_L2Dist(ResidualAfterMeanCor);
                    double ResNormA = ResidualBeforMeanCor.MPI_L2Norm();
                    double ResNormB = ResidualAfterMeanCor.MPI_L2Norm();

                    /*
                    // test code begin ....
                    double[] __ResidualDifference = ResidualAfterMeanCor.CloneAs();
                    __ResidualDifference.AccV(-1.0, ResidualBeforMeanCor);
                    DGField[] __ResidualDifferenceDg = this.CurrentLin.ProlongateRhsToDg(__ResidualDifference, "residualDifference");
                    double[] __ResidualDifferenceNorms = __ResidualDifferenceDg.Select(dgf => (dgf is XDGField xdgf) ? xdgf.L2NormAllSpecies() : dgf.L2Norm()).ToArray();

                    var __LevelSets = __ResidualDifferenceDg.Select(dfg => dfg as XDGField).FirstOrDefault()?.Basis?.Tracker?.LevelSets ?? null;
                    if(__LevelSets != null) {
                        foreach(var ls in __LevelSets) {
                            if(ls is DGField dgLs)
                                dgLs.AddToArray(ref __ResidualDifferenceDg);
                        }
                    }

                    ft.Info($"`FreeMeanValue` check: Original residual: {ResNormA}; residual after distortion {ResNormB}; distance is {Dist}, per component: {__ResidualDifferenceNorms.ToConcatString()}");
                    //Tecplot.Tecplot.PlotFields(__ResidualDifferenceDg, "ResiDiff", 0, 2);
                        
                    // test code end...
                    */

                    double RefVal = Math.Max(Math.Max(BLAS.MachineEps.Sqrt(), Math.Max(ResNormA, ResNormB) * 1e-7), Math.Abs(SolNormA - SolNormB) * 1e-7);
                    if (Dist > RefVal) {
                        double[] ResidualDifference = ResidualAfterMeanCor.CloneAs();
                        ResidualDifference.AccV(-1.0, ResidualBeforMeanCor);
                        DGField[] ResidualDifferenceDg = this.CurrentLin.ProlongateRhsToDg(ResidualDifference, "residualDifference");
                        double[] ResidualDifferenceNorms = ResidualDifferenceDg.Select(dgf => (dgf is XDGField xdgf) ? xdgf.L2NormAllSpecies() : dgf.L2Norm()).ToArray();

                        var LevelSets = ResidualDifferenceDg.Select(dfg => dfg as XDGField).FirstOrDefault()?.Basis?.Tracker?.LevelSets ?? null;
                        if(LevelSets != null) {
                            foreach(var ls in LevelSets) {
                                if(ls is DGField dgLs)
                                    dgLs.AddToArray(ref ResidualDifferenceDg);
                            }
                        }

                        ft.Error($"Something seems wrong with `FreeMeanValue`: drastic change of operator residual; Original residual: {ResNormA}; residual after distortion {ResNormB}; distance is {Dist}, per component: {ResidualDifferenceNorms.ToConcatString()}, reference value {RefVal}");
                        //Console.Error.WriteLine($"Something seems wrong with `FreeMeanValue`: drastic change of operator residual; Original residual: {ResNormA}; residual after distortion {ResNormB}; distance is {Dist}, per component: {ResidualDifferenceNorms.ToConcatString()}, reference value {RefVal}");
                        Tecplot.Tecplot.PlotFields(ResidualDifferenceDg, "ResidualDifference33", 0, 2);
                        //throw new ArithmeticException($"Something seems wrong with `FreeMeanValue`: drastic change of operator residual; Original residual: {ResNormA}; residual after distortion {ResNormB}; distance is {Dist}, reference value {RefVal}");
                        return false;
                    } else {
                        return true;
                    }
                }
            
                return true;
            }
        }


        /// <summary>
        /// DG Coordinate mapping for the base DG/XDG space, without any multigrid fuzz
        /// </summary>
        protected UnsetteledCoordinateMapping ProblemMapping {
            get;
            set;
        }

        /// <summary>
        /// Residual of linearized system, i.e.
        /// <paramref name="out_Resi"/> := RHS - M*<paramref name="in_U"/>
        /// </summary>
        protected void EvalLinearizedResidual(double[] in_U, ref double[] out_Resi) {
            if (this.LinearizationRHS.Length != in_U.Length)
                throw new ApplicationException("internal error");
            if (out_Resi.Length != in_U.Length)
                out_Resi = new double[in_U.Length];

            double L2_in_U = in_U.L2Norm();
            double L2_RHS = this.LinearizationRHS.L2Norm();

            out_Resi.SetV(this.LinearizationRHS, 1.0);
            CurrentLin.OperatorMatrix.SpMV(-1.0, in_U, 1.0, out_Resi);

            double L2_out_Resi = out_Resi.L2Norm();
        }

        /// <summary>
        /// Inner product with respect to the current mass matrix.
        /// 
        /// Note: this is the canonical inner product of the underlying DG-space, since 
        /// for a DG/XDG field represented in an arbitrary basis $\phi_{j}$ one verifies that
        /// \[
        ///         (u, v) = 
        ///     \int_\Omega 
        ///         \left( \sum_{j} \phi_{j} \tilde{u}_{j} \right) 
        ///         \left( \sum_{l} \phi_{l} \tilde{v}_{l} \right) 
        ///     dV = 
        ///     \sum_{j l} \tilde{u}_{j} \tilde{v}_{l} ( \phi_{j}, \phi_{l} )
        ///     =
        ///       \tilde{u}^T M \tilde{v},
        /// \]
        /// where $M$ denotes the mass matrix ($M_{j l} = ( \phi_ { j}, \phi_ { l} )$).
        /// </summary>
        protected double InnerProduct<T1, T2>(T1 vecA, T1 vecB)
            where T1 : IList<double>
            where T2 : IList<double> //
        {
            int L = this.CurrentLin.Mapping.LocalLength;

            if (vecA.Count != L)
                throw new ArgumentException("wrong length.");
            if (vecB.Count != L)
                throw new ArgumentException("wrong length.");


            var MaMa = this.CurrentLin.MassMatrix;

            if (MaMa == null) {
                return vecA.MPI_InnerProd(vecB);
            } else {
                double[] tmp = new double[L];
                MaMa.SpMV(1.0, vecB, 0.0, tmp);
                return vecA.MPI_InnerProd(tmp);
            }
        }

        /// <summary>
        /// Norm induced by <see cref="InnerProduct{T1, T2}(T1, T1)"/>
        /// </summary>
        protected double Norm<T1>(T1 vec)
            where T1 : IList<double> //
        {
            return InnerProduct<T1,T1>(vec, vec).Sqrt();

            /*
            int L0 = CurrentLin.BaseGridProblemMapping.LocalLength;
            double[] tmp = new double[L0];
            //CurrentLin.TransformSolFrom(tmp, vec);
            CurrentLin.TransformRhsFrom(tmp, vec);

            var XFields = CurrentLin.BaseGridProblemMapping.BasisS.Select(bs => new XDGField(bs as XDGBasis)).ToArray();
            var Coordv = new CoordinateVector(XFields);
            Coordv.SetV(tmp);

            double norm_acc = 0;
            for (int i = 0; i < XFields.Length; i++) {
                norm_acc += XFields[i].L2NormSpecies("A").Pow2();
                norm_acc += XFields[i].L2NormSpecies("B").Pow2();
            }

            //Console.WriteLine("New art of norm: " + norm_acc.Sqrt());
            return norm_acc.Sqrt();
            */
        }
    }
}

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
using System.Threading.Tasks;
using BoSSS.Foundation;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Connectors.Matlab;
using ilPSP.Tracing;
using System.IO;
using System.Diagnostics;
using BoSSS.Foundation.XDG;
using NUnit.Framework;

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// Newton solver with various Globalization (i.e. enlargement of convergence radius) options.
    /// </summary>
    /// <remarks>
    /// Further reading:
    /// - Pawlowski et. al., 2006, Globalization Techniques for Newton–Krylov Methods and Applications to the Fully Coupled Solution of the Navier–Stokes Equations, SIAM Review, Vol. 48, No. 4, pp 700-721.
    /// - Pawlowski et. al., 2008, Inexact Newton Dogleg Methods, SIAM Journal on Numerical Analysis, Vol. 46, No. 4, pp 2112-2132.
    /// - Kelley, 2003, Solving Nonlinear Equations with Newton’s Method. Fundamentals of Algorithms. Society for Industrial and Applied Mathematics, 2003. https://doi.org/10.1137/1.9780898718898.
    /// - Kikker, Kummer, Oberlack, 2021, A fully coupled high‐order discontinuous Galerkin solver for viscoelastic fluid flow, IJNMF, No. 6, Vol. 93, 2021, 1736--1758.
    /// </remarks>
    public class Newton : NonlinearSolver {
        /// <summary>
        /// ctor
        /// </summary>
        public Newton(OperatorEvalOrLin __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq, MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig) :
            base(__AssembleMatrix, __AggBasisSeq, __MultigridOperatorConfig) //
        {
        }

        /// <summary>
        /// Maximum number of Newton iterations
        /// </summary>
        public int MaxIter = 50;

        /// <summary>
        /// Minimum number of Newton iterations
        /// </summary>
        public int MinIter = 1;


        /// <summary>
        /// Number of iterations, where Jacobi is not updated. Also known as constant newton method. Default 1, means regular newton.
        /// </summary>
        public int constant_newton_it = 1;

        /// <summary>
        /// If, for a specific homotopy parameter value, Newton does not converges successfully,
        /// (within this number of iterations) a roll-back to the last solution is done and the step with is reduced
        /// </summary>
        public int HomotopyStepLongFail = 20;

        /// <summary>
        /// Convergence criterion for nonlinear iteration;
        /// Setting it to exactly zero should enforce the nonlinear solver to iterate until no improvement on the termination criterion
        /// can be achieved; then, it terminates; this is described in
        /// Kikker, Kummer, Oberlack:
        /// A fully coupled high‐order discontinuous {Galerkin} solver for viscoelastic fluid flow,
        /// IJNMF, No. 6, Vol. 93, 2021, 1736--1758,
        /// section 3.3.1.
        /// </summary>
        public double ConvCrit {
            get {
                return m_ConvCrit;
            }
            set {
                m_ConvCrit = value;
            }
        }

        double m_ConvCrit = 0.0;

        /// <summary>
        /// Maximum number of step-length iterations
        /// </summary>
        public double maxStep = 30;


        /// <summary>
        /// Options for (approximate) solution to the linearizes system
        /// </summary>
        public enum ApproxInvJacobianOptions {

            /// <summary>
            /// Using a matrix-free GMRES implementation, <see cref="Newton.MatrixFreeGMRES"/>
            /// </summary>
            MatrixFreeGMRES = 1,

            /// <summary>
            /// Using the solver <see cref="NonlinearSolver.PrecondConfig"/> for computing Newton corrections
            /// </summary>
            ExternalSolver = 2
        }

        /// <summary>
        /// Options for (approximate) solution to the linearized system
        /// </summary>
        public ApproxInvJacobianOptions ApproxJac = ApproxInvJacobianOptions.ExternalSolver;

        /// <summary>
        /// Options for Globalization, i.e. means to ensure convergence of Newton iterations
        /// when the initial guess is far away from the solution.
        /// </summary>
        public enum GlobalizationOption {

            /// <summary>
            /// Dogleg method according to
            /// Pawlowski et. al., 2006, Globalization Techniques for Newton–Krylov Methods and Applications to the Fully Coupled Solution of the Navier–Stokes Equations, SIAM Review, Vol. 48, No. 4, pp 700-721.
            /// </summary>
            Dogleg = 1,

            /// <summary>
            /// Parabolic line search according to
            /// Kelley, C., Solving Nonlinear Equations with Newton’s Method. Fundamentals of Algorithms. Society for Industrial and Applied Mathematics, 2003. https://doi.org/10.1137/1.9780898718898.
            /// </summary>
            LineSearch = 2,

            /// <summary>
            /// Globalization turned off
            /// </summary>
            None = 0
        }

        /// <summary>
        /// Options for Globalization, i.e. means to ensure convergence of Newton iterations
        /// when the initial guess is far away from the solution.
        /// </summary>
        public GlobalizationOption Globalization = GlobalizationOption.Dogleg;

    
        /// <summary>
        /// Switch the use of the Homotopy-Path (<see cref="IDifferentialOperator.HomotopyUpdate"/>) on/off
        /// </summary>
        public bool UseHomotopy {
            get {
                if(this.AbstractOperator == null)
                    return false;

                return this.AbstractOperator.HomotopyUpdate.Count > 0;
            }
        }


        /// <summary>
        /// Main solver routine
        /// </summary>
        public override bool SolverDriver<S>(CoordinateVector SolutionVec, S RHS) {

            (var gnSuccess, int NoOfIter) = GlobalizedNewton(SolutionVec, RHS); // note: we have to run the default branch first, before we can query 'UseHomotopy' 
            base.NoOfNonlinearIter = NoOfIter;
            if (gnSuccess == false && UseHomotopy == true) {
                (gnSuccess, NoOfIter) = HomotopyNewton(SolutionVec, RHS);
                base.NoOfNonlinearIter = NoOfIter;
                return gnSuccess;
            } else {
                return gnSuccess;
            }
        }


        /// <summary>
        /// Main solver routine
        /// </summary>
        public (bool success, int NoOfIter) GlobalizedNewton<S>(CoordinateVector SolutionVec, S RHS) where S : IList<double> {
            using(var tr = new FuncTrace()) {

                bool success = false;

                // Initialization
                // =============

                double[] CurSol, // "current (approximate) solution", i.e.
                    CurRes; // residual associated with 'CurSol'


                using(new BlockTrace("Slv Init", tr)) {
                    base.Init(SolutionVec, RHS, out CurSol, out CurRes);
                };

                this.CurrentLin.TransformSolFrom(SolutionVec, CurSol);
                EvaluateOperator(1, SolutionVec.Mapping.ToArray(), CurRes, 1.0);
                if(UseHomotopy) { // after the first operator eval, we can access the 'base.AbstractOperator'
                    // don't run this branch - use the Homotopy branch
                    return (false, 0);
                }

                // intial residual evaluation
                double norm_CurRes = CurRes.MPI_L2Norm();
                OnIterationCallback(0, CurSol.CloneAs(), CurRes.CloneAs(), this.CurrentLin);
                double fnorminit = norm_CurRes;
                List<double> normHistory = new List<double>();
                normHistory.Add(norm_CurRes);


                // Main Loop
                // =========

                int itc = 0;
                double TrustRegionDelta = -1; // only used for dogleg (aka Trust-Region) method
                using(new BlockTrace("Slv Iter", tr)) {
                    //while((norm_CurRes > ConvCrit * fnorminit + ConvCrit
                    //    && itc < MaxIter)
                    //    || itc < MinIter) {
                    while(true) {
                        if(CheckTermination(ref success, norm_CurRes, fnorminit, normHistory, itc))
                            break;

                        itc++;
                        NewtonStep(SolutionVec, itc, CurSol, CurRes, 1.0, ref norm_CurRes, ref TrustRegionDelta);
                        normHistory.Add(norm_CurRes);
                    }
                }


                return (success, itc);

            }
        }

        bool CheckTermination(ref bool success, double norm_CurRes, double fnorminit, List<double> normHistory, int itc) {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = false;
                tr.Info($"Checking termination criterion: Iter {itc}, Current residual norm {norm_CurRes}, Initial Residual Norm {fnorminit}");
                
                 if(ConvCrit > 0) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // **if** the convergence criterion is set by the user, we should aim for it
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    tr.Info("NEWTON: going for convergence criterion " + ConvCrit);


                    if(itc < MinIter) {
                        tr.Info("below minimal number of iterations");
                        return false;
                    }

                    if(norm_CurRes <= ConvCrit * fnorminit + ConvCrit) {
                        success = true;
                        tr.Info("Newton iterations SUCCESS!");
                        return true;
                    }

                    if(itc >= MaxIter) {
                        // run out of iterations
                        return true;
                    }

                    return false;
                } else {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // ... otherwise, we do the best we can and iterate until no further improvement can be made
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    tr.Info($"NEWTON: trying to converge as far as possible (ConvCrit is {ConvCrit})");
                    bool terminateLoop = false;



                    //
                    // See:
                    // Kikker, Kummer, Oberlack:
                    // A fully coupled high‐order discontinuous {Galerkin} solver for viscoelastic fluid flow,
                    // IJNMF, No. 6, Vol. 93, 2021, 1736--1758,
                    // section 3.3.1
                    //

                    const double _MinConvCrit = 1e-6;

                    double LastAverageNormReduction() {
                        const int N = 3; // look at latest 3 residuals
                        if(N < 1)
                            throw new ArgumentException();

                        if(normHistory.Count - N < 0)
                            return 1e100; // ignore if we have not at least 'N'  residuals so far.

                        int i0 = normHistory.Count - N;

                        // take the (minimum) skyline to be immune against residual oscillations...
                        double[] NormHistorySkyline = new double[normHistory.Count];
                        NormHistorySkyline[0] = normHistory[0];
                        for(int i = 1; i < normHistory.Count; i++) {
                            NormHistorySkyline[i] = Math.Min(NormHistorySkyline[i - 1], normHistory[i]);
                        }


                        double Avg = 0;
                        double Count = 0;
                        for(int i = i0; i < normHistory.Count - 1; i++) { // look at the last 'N' residual norms...
                            double ResNormReductionFactor = NormHistorySkyline[i] / Math.Max(NormHistorySkyline[i + 1], double.Epsilon);
                            if(ResNormReductionFactor < 1)
                                ResNormReductionFactor = 1; // should never happen anyway due to skylining...
                            Count = Count + 1;
                            Avg = Avg + ResNormReductionFactor;
                        }
                        return Avg / Count;
                    }

                    if(itc >= MinIter) {
                        // only terminate if we reached the minimum number of iterations
                        if(norm_CurRes <= _MinConvCrit * fnorminit + _MinConvCrit) {
                            // reached minimum convergence criterion

                            double ALNR = LastAverageNormReduction();
                            tr.Info($"minimal build-in convergence criterion reached ({_MinConvCrit * fnorminit + _MinConvCrit}), LastAverageNormReduction is {ALNR}, {norm_CurRes} ...");
#if TEST
                            ALNR = 0;
#endif
                            // continue until the solver stalls
                            if(ALNR <= 1.5) {
                                tr.Info("... no sufficient further progress - terminating.");
                                success = true;
                                terminateLoop = true;
                            } else {
                                tr.Info("... but continue as long as we make progress;");
                                //success = true;
                                //terminateLoop = true;
                            }
                        } else {
                            tr.Info($"minimal build-in convergence criterion NOT reached (current residual is {norm_CurRes}, limit is {_MinConvCrit * fnorminit + _MinConvCrit}) yet.");
                        }

                        if(itc >= MaxIter) {
                            // run out of iterations
                            tr.Info($"Maximum number of iterations reached ({MaxIter}) - terminating.");
                            terminateLoop = true;
                        }
                    } else {
                        tr.Info($"**NOT** terminating now, minimal number of iterations (iteration {itc}, minimum is {MinIter}) not reached yet; CurRes = {norm_CurRes}, threshold is {_MinConvCrit * fnorminit + _MinConvCrit}, initial norm was {fnorminit}");
                    }

                    return terminateLoop;
                }
            }
        }
        
     
        /// <summary>
        /// Main solver routine with homotopy;
        /// </summary>
        public (bool success, int NoOfIter) HomotopyNewton<S>(CoordinateVector SolutionVec, S RHS) where S : IList<double> {

            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = false;
                // Initialization
                // =============

                bool success = false;

                double[] CurSol, // "current (approximate) solution", i.e.
                    CurRes; // residual associated with 'CurSol'
                
                                
                using(new BlockTrace("Slv Init", tr)) {
                    base.Init(SolutionVec, RHS, out CurSol, out CurRes);
                };

                this.CurrentLin.TransformSolFrom(SolutionVec, CurSol); // CurSol -> SolutionVec
                EvaluateOperator(1, SolutionVec.Fields, CurRes, 1.0); // seems redundant, but don't dare to remove!

                double HomotopyParameter = 0.0;
                EvaluateOperator(1, SolutionVec.Fields, CurRes, HomotopyParameter);


                // intial residual evaluation
                double norm_CurRes = CurRes.MPI_L2Norm();
                OnIterationCallback(0, CurSol.CloneAs(), CurRes.CloneAs(), this.CurrentLin);
                double fnorminit = norm_CurRes;
                List<double> normHistory = new List<double>();
                if(HomotopyParameter >= 1.0)
                    normHistory.Add(norm_CurRes);

                // Homotopy Solution extrapolation
                // ===============================

                var AcceptedHomoSolutions = new List<(double HomoValue, double[] HomoSolution, double norm_Res)>();

                // Lagrange polynomial (aka. nodal interpolation) over the last 'N' accepted solutions
                double Lagrange(int N, int n, double alpha) {
                    if(N > AcceptedHomoSolutions.Count)
                        throw new ArgumentOutOfRangeException();
                    if(n < 0 || n >= N)
                        throw new ArgumentOutOfRangeException();

                    int M0 = AcceptedHomoSolutions.Count - N;
                    double polyVal = 1.0;

                    double alpha_n = AcceptedHomoSolutions[n + M0].HomoValue;
                    for(int i = 0; i < N; i++) {
                        if(i != n) {
                            double alpha_i = AcceptedHomoSolutions[i + M0].HomoValue;
                            polyVal *= (alpha - alpha_i) / (alpha_n - alpha_i);
                        }
                    }

                    return polyVal;
                }

                // extrapolation from accepted solutions to new homotopy parameter values
                void SolutionExtrapolation(double newHomoParam, int N) {
                    if(N > AcceptedHomoSolutions.Count || N < 0)
                        throw new ArgumentOutOfRangeException();

                    SolutionVec.Clear();
                    int M0 = AcceptedHomoSolutions.Count - N;
                    for(int n = 0; n < N; n++) {
                        if(Math.Abs(Lagrange(N, n, AcceptedHomoSolutions[n + M0].HomoValue) - 1.0) > 1.0e-10)
                            throw new ArithmeticException("Error in Lagrange");
                        for(int i = 0; i < N; i++) {
                            if(i != n) {
                                if(Math.Abs(Lagrange(N, i, AcceptedHomoSolutions[n + M0].HomoValue)) > 1.0e-10)
                                    throw new ArithmeticException("Error in Lagrange (2)");
                            }
                        }

                        double sss = Lagrange(N, n, newHomoParam);
                        SolutionVec.AccV(sss, AcceptedHomoSolutions[n + M0].HomoSolution);
                    }
                }

                // hard-coded constants
                // ===============================
                const double HomotopyStepAcceptedFactor = 100.0; // a solution for a specific homotopy value is accepted,
                //                                                 if the residual is below the convergence criterion, scaled by this factor

                double allowedIncrease = 1e6; // initial value for acceptable residual increase when the homotopy parameter is increased.

                // preventing step width being to long:
                // ---------------------------------------
             
                //const int HomotopyStepLongFail = 20; // If, for a specific homotopy parameter value, Newton does not converges successfully,
                //   (within this number of iterations) a roll-back to the last solution is done and the step with is reduced
                const double StepWidthReductionFactor = 0.2; // Reduction factor (must be smaller 1.0) for the homotopy step (if homotopy step failed)

                // preventing step width being to short:
                // ---------------------------------------
             
                const int HomotopyStepShortFail = 3; // If, for a specific homotopy parameter value, Newton converges successfully
                //       with less than this number of iterations, the step width should be increased.
                const double StepWidthIncreaseFactor = 8; // Increase factor (must be bigger then the inverse of the reduction factor) for the homotopy step


                // main Loop
                // ===============================
                
                int IterCounter = 0;
                int HomoStepCounter = 0; // no of iterations spent in current Homotopy step
                double TrustRegionDelta = -1; // only used for dogleg (aka Trust-Region) method
                using(new BlockTrace("Slv Iter", tr)) {

                    bool ResidualCrit(double fac = 1.0) {
                        double __ConvCrit = 1.0e-6;
                        double thresh = __ConvCrit * fnorminit * fac + __ConvCrit * fac;
                        bool ret = norm_CurRes <= thresh;
                        tr.Info($"Pre-final convergence in Homotopy loop: {ret} (threshold = {thresh}, cur Res = {norm_CurRes}, init Res = {fnorminit}, ConvCrit = {ConvCrit}, fac = {fac})");
                        return ret;
                    }


                    //while(!ResidualCrit() || IterCounter < MinIter || HomotopyParameter < 1.0) {
                    while(true) {
                        //if(IterCounter > MaxIter) {
                        //    Console.WriteLine($"Failed Newton: maximum number of iterations {MaxIter} exceeded.");
                        //    break;
                        //}


                        if(HomotopyParameter >= 1.0) {
                            if(CheckTermination(ref success, norm_CurRes, fnorminit, normHistory, IterCounter))
                                break;
                        } else {
                            if(IterCounter >= MinIter) {
                                if(IterCounter >= MaxIter) {
                                    tr.Error($"Failed Newton: maximum number of iterations {MaxIter} exceeded.");
                                    break;
                                }
                            }
                        }

                        // test 

                        bool GoodForHomoIncrease = ResidualCrit(HomotopyStepAcceptedFactor) && HomotopyParameter < 1.0;
                        bool BadHomo = (HomoStepCounter >= HomotopyStepLongFail
                                        //|| TrustRegionDelta < 5e-6 // Don't use this! It may cause a wrong/unintended fail -- and reduction of homotopy parameter -- at the end of the iteration, when we are almost converged.
                        ) && (AcceptedHomoSolutions.Count > 0);

                        if (TrustRegionDelta < 5e-6) {
                            tr.Info("Homotopy decrease because of low trust region");
                        }

                        if(!GoodForHomoIncrease && BadHomo) {
                            allowedIncrease *= StepWidthReductionFactor;
                            
                            if(allowedIncrease < 1) {
                                tr.Error($"Failed Newton: unable to raise Homotopy parameter without loosing convergence.");
                                break;
                            }

                            var ll = AcceptedHomoSolutions.Last();
                            HomotopyParameter = ll.HomoValue;
                            SolutionVec.SetV(ll.HomoSolution);
                            this.UpdateLinearization(SolutionVec.Fields, HomotopyParameter);
                            
                            this.CurrentLin.TransformSolInto(SolutionVec, CurSol);
                            this.EvaluateOperator(1, SolutionVec.Fields, CurRes, HomotopyParameter);
                            norm_CurRes = CurRes.MPI_L2Norm();

                            //Console.WriteLine($"Homo Rollback: homotopy value {Homotopy}, old Res norm {ll.norm_Res}, now {norm_CurRes}");
                        }

                        if(GoodForHomoIncrease) {
                            // +++++++++++++++
                            // accept solution
                            // +++++++++++++++
                            AcceptedHomoSolutions.Add((HomotopyParameter, SolutionVec.ToArray(), norm_CurRes));
                            while(AcceptedHomoSolutions.Count > 8)
                                AcceptedHomoSolutions.RemoveAt(0);
                        }


                        if(GoodForHomoIncrease || BadHomo) {
                            // +++++++++++++++++++++++++
                            // Update the homotopy value
                            // +++++++++++++++++++++++++

                            Debug.Assert(AcceptedHomoSolutions.Count > 0);

                            double targetResNorm = norm_CurRes * allowedIncrease;

                            tr.Info($"Trying to raise homotopy value (target residual is {targetResNorm}) ...");

                            int i;
                            int LastN = -1;

                            //int PlotCount = 0;

                            for(i = 0; i < 100; i++) {
                                double tryHomotopyValue;
                                if(i == 0)
                                    tryHomotopyValue = 1.0; // nail it exactly, avoid any round-off-error
                                else
                                    tryHomotopyValue = Math.Pow(2.0, -i) * (1.0 - HomotopyParameter) + HomotopyParameter;

                                LastN = -1;
                                for(int N = AcceptedHomoSolutions.Count; N > 0; N--) {
                                    SolutionExtrapolation(tryHomotopyValue, N);
                                    EvaluateOperator(1, SolutionVec.Fields, CurRes, tryHomotopyValue);
                                    double norm_TryValue = CurRes.MPI_L2Norm();

                                    //Console.WriteLine($"({i},{N}) resNorm: {norm_TryValue}");

                                    if(norm_TryValue < targetResNorm) {
                                        HomotopyParameter = tryHomotopyValue;

                                        if(N > 1) {
                                            this.UpdateLinearization(SolutionVec.Fields, HomotopyParameter);
                            
                                            this.CurrentLin.TransformSolInto(SolutionVec, CurSol);
                                            this.EvaluateOperator(1, SolutionVec.Fields, CurRes, HomotopyParameter);
                                            norm_CurRes = CurRes.MPI_L2Norm();
                                            norm_TryValue = norm_CurRes;
                                        }


                                        norm_CurRes = norm_TryValue;
                                        LastN = N;
                                        break;
                                    }
                                }

                                if(LastN > 0)
                                    break;
                            }

                            tr.Info($"Raised to {HomotopyParameter}, (new residual {norm_CurRes}, tried {i+1} times, used {LastN}-th order extrapolation.");

                            if(HomoStepCounter < HomotopyStepShortFail)
                                allowedIncrease *= StepWidthIncreaseFactor;
                            HomoStepCounter = 0;
                        }


                        // perform a Newton step
                        // ---------------------
                        IterCounter++;
                        HomoStepCounter++;
                        NewtonStep(SolutionVec, IterCounter, CurSol, CurRes, HomotopyParameter, ref norm_CurRes, ref TrustRegionDelta);
                        if(HomotopyParameter >= 1.0)
                            normHistory.Add(norm_CurRes);
                    }
                }



                return (success, IterCounter);
            }
        }

        private void NewtonStep(CoordinateVector SolutionVec, int itc, double[] CurSol, double[] CurRes, double HomotopyValue, ref double norm_CurRes, ref double TrustRegionDelta) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = false;
                // computation of Newton step
                // --------------------------

                double[] step = new double[CurSol.Length];

                // How should the inverse of the Jacobian be approximated?
                if (ApproxJac == ApproxInvJacobianOptions.MatrixFreeGMRES) {
                    // ++++++++++++++++++++++++++
                    // Option: Matrix-Free GMRES
                    // ++++++++++++++++++++++++++


                    using (var mtxFreeSlv = new MatrixFreeGMRES() { owner = this, HomotopyValue = HomotopyValue }) {
                        double thresh = norm_CurRes * 1e-5;
                        mtxFreeSlv.TerminationCriterion = (iter, R0_l2, R_l2) => {
                            return ((R_l2 > thresh) && (iter < 100), R_l2 <= thresh);
                        };

                        step = mtxFreeSlv.Krylov(SolutionVec, CurSol, CurRes, out double errstep);
                        step.ScaleV(-1);

                    }

                } else if (ApproxJac == ApproxInvJacobianOptions.ExternalSolver) {
                    // +++++++++++++++++++++++++++++
                    // Option: use 'external' solver
                    // +++++++++++++++++++++++++++++

                    using (var solver = this.PrecondConfig.CreateInstance(CurrentLin)) {

                        step.ClearEntries();
                        solver.ResetStat();

                        if (solver is IProgrammableTermination pt) {
                            // iterative solver with programmable termination is used - so use it

                            double thresh = norm_CurRes * 1e-5;
                            tr.Info($"Inexact Newton: setting convergence threshold to {thresh:0.##E-00}");
                            pt.TerminationCriterion = (iter, R0_l2, R_l2) => {
                                return ((R_l2 > thresh) && (iter < 100), R_l2 < thresh);
                            };
                        }

                        //dgREs = CurrentLin.ProlongateRhsToDg(CurRes, "Rhs_");
                        //Console.WriteLine("RHS in ref cell: " + dgREs[2].GetMeanValue(CurrentLin.ReferenceCell_local));
                        solver.Solve(step, CurRes);
                        step.ScaleV(-1);
                    }
                } else {
                    throw new NotImplementedException($"approximation option {ApproxJac} for the Jacobian seems not to be existent.");
                }

                
                // globalization
                // -------------
                double[] OldSolClone;
                if (base.AbstractOperator.SolverSafeguard != null) {
                    OldSolClone = SolutionVec.ToArray();
                } else {
                    OldSolClone = null;
                }

                //var DgOldSol = CurrentLin.ProlongateSolToDg(CurSol, "OldSol_");
                //var DgStep = CurrentLin.ProlongateSolToDg(step, "Step_");
                //DGField pressure = SolutionVec.Mapping.Fields[2];
                //Console.WriteLine("Mean value before correction: " + pressure.GetMeanValue(CurrentLin.ReferenceCell_local));

                tr.Info("Using Globalization: " + Globalization);
                switch (Globalization) {
                    case GlobalizationOption.Dogleg:
                        DogLeg(SolutionVec, CurSol, CurRes, step, HomotopyValue, itc, ref TrustRegionDelta);
                        break;

                    case GlobalizationOption.LineSearch:
                        LineSearch(SolutionVec, CurSol, CurRes, step, HomotopyValue);
                        break;

                    case GlobalizationOption.None:
                        ZeroGlobalization(SolutionVec, CurSol, CurRes, step, HomotopyValue);
                        break;

                    default:
                        throw new NotImplementedException();
                }

                if (base.AbstractOperator.SolverSafeguard != null) {
                    var newSol = SolutionVec.Fields.ToArray();
                    var oldSol = newSol.Select(f => f.CloneAs()).ToArray();
                    var oldSolVec = new CoordinateVector(oldSol);
                    oldSolVec.SetV(OldSolClone, 1.0);

                    base.AbstractOperator.SolverSafeguard(oldSol, newSol);
                }

                // fix the pressure
                // ----------------
                if (CurrentLin.FreeMeanValue.Any()) {
                    if (itc == 0 || itc % 5 == 0) // execute this expensive test not to often.
                        base.TestFreeMeanValue(SolutionVec, HomotopyValue);

                    DGField[] flds = SolutionVec.Mapping.Fields.ToArray();
                    bool[] FreeMeanValue = CurrentLin.FreeMeanValue;
                    if (flds.Length != FreeMeanValue.Length)
                        throw new ApplicationException();

                    int RefCellLocal = CurrentLin.ReferenceCell_local;

                    double[] MeanValues = new double[flds.Length];
                    for (int iFld = 0; iFld < flds.Length; iFld++) {

                        if (FreeMeanValue[iFld]) {
                            //double mean = flds[iFld].GetMeanValueTotal(null);

                            if (RefCellLocal >= 0) {
                                MeanValues[iFld] = flds[iFld].GetMeanValue(RefCellLocal);

                            }

                            //flds[iFld].AccConstant(-mean);
                        }
                    }
                    MeanValues = MeanValues.MPISum();
                    for (int iFld = 0; iFld < flds.Length; iFld++) {

                        if (FreeMeanValue[iFld]) {
                            tr.Info("Mean value in reference cell: " + MeanValues[iFld]);
                            flds[iFld].AccConstant(-MeanValues[iFld]);
                            //Console.WriteLine("Mean value after correction: " + flds[iFld].GetMeanValue(RefCellLocal));
                        }
                    }
                }

                // update linearization
                // --------------------
                if (itc % constant_newton_it == 0) {

                    base.Update(SolutionVec.Mapping.Fields, CurSol, HomotopyValue);

                    if (constant_newton_it != 1) {
                        Console.WriteLine("Jacobian is updated: it {0}", itc);
                    }
                }

                //// plotting during Newton iterations:  
                //var DgSolution = CurrentLin.ProlongateSolToDg(CurSol, "Sol_");
                //Tecplot.Tecplot.PlotFields(DgSolution.Cat(DgOldSol, DgStep, dgREs), "DuringNewton-" + itc, itc, 3);

                // residual evaluation & callback
                // ------------------------------
                EvaluateOperator(1, SolutionVec.Mapping.Fields, CurRes, HomotopyValue, ApplyRef: true);
                norm_CurRes = CurRes.MPI_L2Norm();

                //double norm_CurResX = base.Norm(CurRes);
                //tr.Info($"Norm BEFORE going into Iter Callback: {norm_CurResX:0.####E-00}");

                OnIterationCallback(itc, CurSol.CloneAs(), CurRes.CloneAs(), this.CurrentLin);
            }
        }
        /// <summary>
        /// Newton Globalization via parabolic line search
        /// </summary>
        /// <param name="SolutionVec">
        /// output: updated solution in original DG coordinates
        /// </param>
        /// <param name="CurSol">
        /// input: current solution in the preconditioned DG coordinates
        /// </param>
        /// <param name="CurRes">
        /// input: residual for <paramref name="CurSol"/>
        /// </param>
        /// <param name="step">
        /// input: Newton step
        /// </param>
        /// <param name="HomotopyValue">
        /// current parameter for the homotopy curve
        /// </param>
        /// <returns>
        /// Updated Solution
        /// </returns>
        private void LineSearch(CoordinateVector SolutionVec, double[] CurSol, double[] CurRes, double[] step, double HomotopyValue) {
            using (var tr = new FuncTrace()) {
                //tr.InfoToConsole = true;
                double[] TempSol;
                double[] TempRes;

                double alpha = 1E-4;
                double sigma1 = 0.5;

                // Start line search
                double lambda = 1;
                double lamm = 1;
                double lamc = lambda;
                int iarm = 0;
                TempSol = CurSol.CloneAs();
                TempSol.AccV(lambda, step);
                this.CurrentLin.TransformSolFrom(SolutionVec, TempSol);

                TempRes = new double[TempSol.Length];
                EvaluateOperator(1, SolutionVec.Mapping.Fields, TempRes, HomotopyValue);

                double nf0 = base.Norm(CurRes);  //.L2NormPow2().MPISum().Sqrt();
                double nft = base.Norm(TempRes); //.L2NormPow2().MPISum().Sqrt();
                double ff0 = nf0 * nf0; // residual norm of current solution ^2
                double ffc = nft * nft; //
                double ffm = nft * nft; //

                tr.Info($"LineSearch: Current Solution Residual: {nf0:0.####E-00}");
                tr.Info($"LineSearch: Full Newton Step Residual: {nft:0.####E-00}");

                // Control of the the step size
                while (nft >= (1 - alpha * lambda) * nf0 && iarm < maxStep) {

                    // Line search starts here
                    if (iarm == 0)
                        lambda = sigma1 * lambda;
                    else
                        lambda = parab3p(lamc, lamm, ff0, ffc, ffm); // ff0: curent sol, ffc: most recent reduction, ffm: previous reduction

                    // Update x;
                    TempSol = CurSol.CloneAs();
                    TempSol.AccV(lambda, step);
                    lamm = lamc;
                    lamc = lambda;

                    this.CurrentLin.TransformSolFrom(SolutionVec, TempSol);
                    EvaluateOperator(1, SolutionVec.Mapping.Fields, TempRes, HomotopyValue);
                    
                    nft = base.Norm(TempRes); //.L2NormPow2().MPISum().Sqrt();
                    ffm = ffc;
                    ffc = nft * nft;
                    iarm++;

                    tr.Info("LineSearch: Reduced     Step Residual: " + nft + " (lambda = " + lambda + ")");
                    
                }
                

                // transform solution back to 'original domain'
                // to perform the linearization at the new point...
                // (and for Level-Set-Updates ...)
                this.CurrentLin.TransformSolFrom(SolutionVec, TempSol);
            }
        }


        /// <summary>
        /// Zero Globalization
        /// </summary>
        /// <param name="SolutionVec">
        /// output: updated solution in original DG coordinates
        /// </param>
        /// <param name="CurSol">
        /// input: current solution in the preconditioned DG coordinates
        /// </param>
        /// <param name="CurRes">
        /// input: residual for <paramref name="CurSol"/>
        /// </param>
        /// <param name="step">
        /// input: Newton step
        /// </param>
        /// <param name="HomotopyValue">
        /// current paramter for the homotopy curve
        /// </param>
        /// <returns>
        /// Updated Solution
        /// </returns>
        private void ZeroGlobalization(CoordinateVector SolutionVec, double[] CurSol, double[] CurRes, double[] step, double HomotopyValue) {

            double[] TempSol;


  
            TempSol = CurSol.CloneAs();
            TempSol.AccV(1.0, step);
            this.CurrentLin.TransformSolFrom(SolutionVec, TempSol);

           


        }

        /// <summary>
        /// Newton Globalization via Dogleg Method (a trust region approach)
        /// </summary>
        /// <param name="SolutionVec">
        /// output: updated solution in original DG coordinates
        /// </param>
        /// <param name="CurSol">
        /// input: current solution in the preconditioned DG coordinates
        /// </param>
        /// <param name="CurRes">
        /// input: residual for <paramref name="CurSol"/>
        /// </param>
        /// <param name="stepIN">
        /// input: (inexact) Newton step
        /// </param>
        /// <param name="NewtonIterCnt">
        /// newton iteration counter, starts at 1
        /// </param>
        /// <param name="TrustRegionDelta">
        /// Estimated in <paramref name="NewtonIterCnt"/>==1;
        /// value must be stored externally for later iterations.
        /// </param>
        /// <param name="HomotopyValue">
        /// <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
        /// </param>
        /// <remarks>
        /// See:
        /// - Pawlowski et. al., 2006, Globalization Techniques for Newton–Krylov Methods and Applications to the Fully Coupled Solution of the Navier–Stokes Equations, SIAM Review, Vol. 48, No. 4, pp 700-721.
        /// - Pawlowski et. al., 2008, Inexact Newton Dogleg Methods, SIAM Journal on Numerical Analysis, Vol. 46, No. 4, pp 2112-2132.
        /// </remarks>
        private void DogLeg(CoordinateVector SolutionVec, double[] CurSol, double[] CurRes, double[] stepIN, double HomotopyValue, int NewtonIterCnt, ref double TrustRegionDelta) {

            // algorithm constants taken from [Pawlowski et. al. 2006]
            // =======================================================
            const double t = 1.0e-4;
            const double delta_min = 1e-6;
            const double delta_max = 1e10;

            // initial estimate of trust region width in first iteration
            // =========================================================
            if(NewtonIterCnt < 1)
                throw new ArgumentException();
            if(NewtonIterCnt == 1) {
                double norm_step = stepIN.MPI_L2Norm();
                if(norm_step < delta_min)
                    TrustRegionDelta = 2 * delta_min;
                else
                    TrustRegionDelta = norm_step;

                TrustRegionDelta = Math.Min(delta_max, TrustRegionDelta);
            }


            // TODO change later: always use a very small initial trust region!
            // TrustRegionDelta = TrustRegionDelta / 10;

            if(TrustRegionDelta < delta_min || TrustRegionDelta > delta_max)
                throw new ArithmeticException("trust region width out of allowed range");

            // compute Cauchy point
            // ====================
            double[] stepCP;
            {
                // step 1: calculate direction of Cauchy point / direction of steepest decent
                double[] dk = new double[CurSol.Length];
                BlockMsrMatrix JacTransp = this.CurrentLin.OperatorMatrix.Transpose();
                JacTransp.SpMV(-1.0, CurRes, 0.0, dk); // eventually, replace this by a SpMV - transpose


                // step 2: computing Cauchy point
                double[] Mdk = new double[dk.Length];
                this.CurrentLin.OperatorMatrix.SpMV(1.0, dk, 0.0, Mdk);
                double[] a = (new[] { CurRes.InnerProd(Mdk), Mdk.L2NormPow2() }).MPISum();

                double lambda = -a[0] / Math.Max(double.Epsilon*100, a[1]);

                stepCP = dk;
                stepCP.ScaleV(lambda);
                dk = null; // invalidate

                // test:
                double[] ResAtCP = CurRes.CloneAs();
                this.CurrentLin.OperatorMatrix.SpMV(1.0, stepCP, 1.0, ResAtCP);

                double ResNormAtCP = ResAtCP.MPI_L2Norm();
                double ResNormAtX0 = CurRes.MPI_L2Norm();

                Assert.LessOrEqual(ResNormAtCP, ResNormAtX0, "Something wrong in calculation of Cauchy Point -- residual of linear model increased.");
            }

            // using only Cauchy-Point
            //var _NewSol = CurSol.CloneAs();
            //_NewSol.AccV(1.0, stepCP);
            //this.CurrentLin.TransformSolFrom(SolutionVec, _NewSol);
            //return;

            // find point on Dogleg curve, within the trust region
            // ===================================================
            double[] step = new double[stepIN.Length];
            double l2_stepCP = stepCP.MPI_L2Norm();
            double l2_stepIN = stepIN.MPI_L2Norm();
            double[] NewSol;
            void PointOnDogleg(double _TrustRegionDelta) {
                if(l2_stepIN <= _TrustRegionDelta) {
                    // use Newton Step
                    //Console.WriteLine($"       -------- using Newton step (delta = {_TrustRegionDelta})");
                    step.SetV(stepIN);
                } else {
                    if(l2_stepCP < _TrustRegionDelta) {
                        // interpolate between Cauchy-point and Newton-step
                        Console.WriteLine($"Info: Newton solver - between Cauchy point and Newton step (delta = {_TrustRegionDelta})");
                        Debug.Assert(l2_stepCP * 0.99999 <= _TrustRegionDelta); // Cauchy Point is INSIDE   trust region
                        Debug.Assert(l2_stepIN * 1.00001 >= _TrustRegionDelta); // Newton Step  is OUTSIDE  trust region
                        Debug.Assert(l2_stepCP <= l2_stepIN * 1.00001); // consequently, the Newton Step must be larger than the Cauchy point

                        //double tau = (l2_stepCP - _TrustRegionDelta) / (l2_stepCP - l2_stepIN); // nominator and denominator must be negative!
                        //                                                                        //double tau2 = (l2_stepCP + TrustRegionDelta) / (l2_stepCP - l2_stepIN); // other solution of quadratic problem, probably negative
                        double tau;
                        {
                            double A = l2_stepCP, B = l2_stepIN, C = stepCP.MPI_ddot(stepIN);
                            tau = (A.Pow2() - C + Math.Sqrt((A.Pow2() + B.Pow2() - 2 * C) * _TrustRegionDelta.Pow2() - A.Pow2() * B.Pow2() + C.Pow2()))
                                / (A.Pow2() + B.Pow2() - 2 * C);
                            if(!(tau >= -0.00001 && tau <= 1.00001) || tau.IsNaN() || tau.IsInfinity())
                                throw new ArithmeticException();
                        }
                        // do interpolation
                        step.SetV(stepCP, (1 - tau));
                        step.AccV(tau, stepIN);
                        Debug.Assert(Math.Abs((step.MPI_L2Norm() / _TrustRegionDelta) - 1.0) <= 1.0e-3, "interpolation step went wrong");
                    } else {
                        // use reduced Cauchy point
                        Console.WriteLine($"Info: Newton solver - using Cauchy point (delta = {_TrustRegionDelta})");
                        Debug.Assert(l2_stepCP * 1.00001 >= _TrustRegionDelta); // Cauchy Point is outside trust region
                        step.SetV(stepCP, alpha: (_TrustRegionDelta / l2_stepCP));
                        Debug.Assert(Math.Abs((step.MPI_L2Norm() / _TrustRegionDelta) - 1.0) <= 1.0e-3, "scaling step went wrong");
                    }
                }

                NewSol = CurSol.CloneAs();
                NewSol.AccV(1.0, step);
            }

            PointOnDogleg(TrustRegionDelta);

            // check and adapt trust region
            // ============================
            double l2_CurRes = CurRes.MPI_L2Norm();

            double[] temp = new double[CurRes.Length];

            // predicted residual reduction
            double pred() {

                temp.SetV(CurRes);
                this.CurrentLin.OperatorMatrix.SpMV(1.0, step, 1.0, temp);

                return l2_CurRes - temp.MPI_L2Norm();
            }

            // actual residual reduction
            double ared() {
                this.CurrentLin.TransformSolFrom(SolutionVec, NewSol);
                base.EvaluateOperator(1.0, SolutionVec.Fields, temp, HomotopyValue);

                return l2_CurRes - temp.MPI_L2Norm();
            }

            // trust region adaptation loop
            double last_ared = ared();
            double last_pred = pred();
            while(last_ared < t * last_pred) {
                double newTrustRegionDelta = TrustRegionDelta * 0.5;
                if(newTrustRegionDelta <= delta_min)
                    break;

                PointOnDogleg(TrustRegionDelta);

                TrustRegionDelta = Math.Max(delta_min, newTrustRegionDelta);


                last_ared = ared();
                last_pred = pred();
            }

            // final trust region update
            // =========================
            {
                // taken from [Pawlovski et.al. 2006], section 2.4;
                // originally from J. E. Dennis, Jr. and R. B. Schnabel. Numerical Methods for Unconstrained Optimization and Nonlinear Equations. Series in Automatic Computation. Prentice-Hall, Englewood Cliffs, NJ, 1983.

                const double rho_s = 0.1;
                const double rho_e = 0.75;
                const double beta_s = 0.25;
                const double beta_e = 4.0;

                if(last_ared / last_pred < rho_s && l2_stepIN < TrustRegionDelta) {
                    TrustRegionDelta = Math.Max(l2_stepIN, delta_min);
                } else if(last_ared / last_pred < rho_s) {
                    TrustRegionDelta = Math.Max(beta_s * TrustRegionDelta, delta_min); // shrinking
                } else if(last_ared / last_pred > rho_e) {
                    TrustRegionDelta = Math.Min(beta_e * TrustRegionDelta, delta_max); // enhancing
                }
            }



            // return updated solution
            // =======================

            // remark: 'SolutionVec' should be up to date due to the last call to 'ared()'

            //this.CurrentLin.TransformSolFrom(SolutionVec, NewSol);
            return;
        }


        /// <summary>
        /// container class for all matrix-free GMRES routines
        /// </summary>
        class MatrixFreeGMRES : ISolverSmootherTemplate, IProgrammableTermination {

            /// <summary>
            /// ctor
            /// </summary>
            public MatrixFreeGMRES() {
                this.TerminationCriterion = DefaultTermination;
            }


            /// <summary>
            /// %
            /// </summary>
            public Newton owner;

            /// <summary>
            /// <see cref="IDifferentialOperator.CurrentHomotopyValue"/>
            /// </summary>
            internal double HomotopyValue;

            /*
            /// <summary>
            /// Maximum number of GMRES(m) restarts
            /// </summary>
            public int restart_limit = 1000;
            */

            /// <summary>
            /// Maximum dimension of the krylov subspace. Equals m in GMRES(m)
            /// </summary>
            public int maxKrylovDim = 200;

            public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion { 
                get;
                set;
            }

            static private (bool bNotTerminate, bool bSuccess) DefaultTermination(int iter, double R0_l2, double R_l2) {
                if(iter > 100)
                    return (false, false); // fail

                if(R_l2 < R0_l2 * 10e-8 + 10e-8)
                    return (false, true); // success

                return (true, false); // keep running
            }


            int ThisRunFirstIter;
            double rho0; // residual in first run

            bool Termination(double R_l2) {
                bool ret = TerminationCriterion(this.ThisLevelIterations - ThisRunFirstIter, rho0, R_l2).bNotTerminate;
                if(ret)
                    Converged = true;
                return ret;
            }

            /// <summary>
            /// <see cref="ISolverSmootherTemplate.IterationsInNested"/>
            /// </summary>
            public int IterationsInNested {
                get {
                    if(this.PrecondInstance != null) {
                        return this.PrecondInstance.ThisLevelIterations;
                    } else {
                        return 0;
                    }
                }
            }

            /// <summary>
            /// <see cref="ISolverSmootherTemplate.ThisLevelIterations"/>
            /// </summary>
            public int ThisLevelIterations {
                get;
                private set;
            }

            /// <summary>
            /// true if the last solver run actually worked, which is not the case very often
            /// </summary>
            public bool Converged {
                get;
                private set;
            }

            /*
            /// <summary>
            /// Convergence for Krylov and GMRES iterations
            /// </summary>
            public double GMRESConvCrit = 1e-6;
            */

            /// <summary>
            /// Preconditioned GMRES, using <see cref="NonlinearSolver.PrecondConfig"/> as a preconditioner
            /// </summary>
            /// <param name="SolutionVec">Current Point</param>
            /// <param name="f0">Function at current point</param>
            /// <param name="xinit">initial iterate</param>
            /// <param name="errstep">error of step</param>
            /// <param name="currentX"></param>
            /// <returns></returns>
            double[] Solve(CoordinateVector SolutionVec, double[] currentX, double[] f0, double[] xinit, out double errstep) {
                using(var tr = new FuncTrace()) {
                    int n = f0.Length;

                    int reorth = 1; // Orthogonalization method -> 1: Brown/Hindmarsh condition, 3: Always reorthogonalize

                    // RHS of the linear equation system
                    double[] b = new double[n];
                    b.AccV(1, f0);

                    double[] x = new double[n];
                    double[] r = new double[n];

                    int Nloc = owner.CurrentLin.OperatorMatrix.RowPartitioning.LocalLength;
                    long Ntot = owner.CurrentLin.OperatorMatrix.RowPartitioning.TotalLength;

                    r = b;

                    //Initial solution
                    if(xinit.L2Norm() != 0) {
                        x = xinit.CloneAs();
                        r.AccV(-1, dirder(SolutionVec, currentX, x, f0));
                    }
                    // Precond = null;
                    if(this.PrecondInstance != null) {
                        var temp2 = r.CloneAs();
                        r.ClearEntries();
                        //this.OpMtxRaw.InvertBlocks(OnlyDiagonal: false, Subblocks: false).SpMV(1, temp2, 0, r);
                        this.PrecondInstance.Solve(r, temp2);
                    }

                    int m = maxKrylovDim;
                    double[][] V = (m + 1).ForLoop(i => new double[Nloc]); //   V(1:n,1:m+1) = zeros(n,m);
                    MultidimensionalArray H = MultidimensionalArray.Create(m + 1, m + 1); //   H(1:m,1:m) = zeros(m,m);
                    double[] c = new double[m + 1];
                    double[] s = new double[m + 1];
                    double[] y;
                    double rho = r.L2NormPow2().MPISum().Sqrt();
                    errstep = rho;
                    double[] g = new double[m + 1];
                    g[0] = rho;
                    rho0 = rho;

                    Console.WriteLine("Error NewtonGMRES:   " + rho);

                    // Termination of entry
                    //if(rho < GMRESConvCrit)
                    if (Termination(rho))
                        return currentX.ToArray();

                    V[0].SetV(r, alpha: (1.0 / rho));
                    double beta = rho;
                    int k = 1;

                    while((!Termination(rho)) && k <= m) {
                        V[k].SetV(dirder(SolutionVec, currentX, V[k - 1], f0));
                        //CurrentLin.OperatorMatrix.SpMV(1.0, V[k-1], 0.0, temp3);
                        // Call directional derivative
                        //V[k].SetV(f0);

                        ThisLevelIterations++;

                        if(this.PrecondInstance != null) {
                            var temp3 = V[k].CloneAs();
                            V[k].ClearEntries();
                            //this.OpMtxRaw.InvertBlocks(false,false).SpMV(1, temp3, 0, V[k]);
                            this.PrecondInstance.Solve(V[k], temp3);
                        }

                        double normav = V[k].MPI_L2Norm();

                        // Modified Gram-Schmidt
                        for(int j = 1; j <= k; j++) {
                            H[j - 1, k - 1] = GenericBlas.InnerProd(V[k], V[j - 1]).MPISum();
                            V[k].AccV(-H[j - 1, k - 1], V[j - 1]);
                        }
                        H[k, k - 1] = V[k].L2NormPow2().MPISum().Sqrt();
                        double normav2 = H[k, k - 1];


                        // Reorthogonalize ?
                        if((reorth == 1 && Math.Round(normav + 0.001 * normav2, 3) == Math.Round(normav, 3)) || reorth == 3) {
                            for(int j = 1; j <= k; j++) {
                                double hr = GenericBlas.InnerProd(V[k], V[j - 1]).MPISum();
                                H[j - 1, k - 1] = H[j - 1, k - 1] + hr;
                                V[k].AccV(-hr, V[j - 1]);
                            }
                            H[k, k - 1] = V[k].L2NormPow2().MPISum().Sqrt();
                        }

                        // Watch out for happy breakdown
                        if(H[k, k - 1] != 0)
                            V[k].ScaleV(1 / H[k, k - 1]);



                        // Form and store the information for the new Givens rotation
                        //if (k > 1) {
                        //    // for (int i = 1; i <= k; i++) {
                        //    H.SetColumn(k - 1, givapp(c.GetSubVector(0, k - 1), s.GetSubVector(0, k - 1), H.GetColumn(k - 1), k - 1));
                        //    //}
                        //}

                        // Givens rotation from SoftGMRES
                        double temp;
                        for(int l = 1; l <= k - 1; l++) {
                            // apply Givens rotation, H is Hessenbergmatrix
                            temp = c[l - 1] * H[l - 1, k - 1] + s[l - 1] * H[l + 1 - 1, k - 1];
                            H[l + 1 - 1, k - 1] = -s[l - 1] * H[l - 1, k - 1] + c[l - 1] * H[l + 1 - 1, k - 1];
                            H[l - 1, k - 1] = temp;
                        }
                        //	 [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
                        rotmat(out c[k - 1], out s[k - 1], H[k - 1, k - 1], H[k + 1 - 1, k - 1]);
                        temp = c[k - 1] * g[k - 1]; //                       % approximate residual norm
                        H[k - 1, k - 1] = c[k - 1] * H[k - 1, k - 1] + s[k - 1] * H[k + 1 - 1, k - 1];
                        H[k + 1 - 1, k - 1] = 0.0;


                        // Don't divide by zero if solution has  been found
                        var nu = (H[k - 1, k - 1].Pow2() + H[k, k - 1].Pow2()).Sqrt();
                        if(nu != 0) {
                            //c[k - 1] = H[k - 1, k - 1] / nu;
                            //s[k - 1] = H[k, k - 1] / nu;
                            //H[k - 1, k - 1] = c[k - 1] * H[k - 1, k - 1] - s[k - 1] * H[k, k - 1];
                            //H[k, k - 1] = 0;

                            // givapp for g
                            g[k + 1 - 1] = -s[k - 1] * g[k - 1];
                            g[k - 1] = temp;

                            //var w1 = c[k - 1] * g[k - 1] - s[k - 1] * g[k];
                            //var w2 = s[k - 1] * g[k - 1] + c[k - 1] * g[k];
                            //g[k - 1] = w1;
                            //g[k] = w2;
                        }

                        rho = Math.Abs(g[k]);

                        Console.WriteLine("Error NewtonGMRES:   " + rho);

                        k++;

                    }

                    Console.WriteLine("GMRES completed after:   " + k + "steps");

                    k--;



                    // update approximation and exit
                    //y = H(1:i,1:i) \ g(1:i);
                    y = new double[k];
                    H.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { k - 1, k - 1 })
                        .Solve(y, g.GetSubVector(0, k));

                    int totalIter = k;

                    // x = x + V(:,1:i)*y;
                    for(int ii = 0; ii < k; ii++) {
                        x.AccV(y[ii], V[ii]);
                    }

                    // update approximation and exit
                    //using (StreamWriter writer = new StreamWriter(m_SessionPath + "//GMRES_Stats.txt", true)) {
                    //    writer.WriteLine("");
                    //}

                    errstep = rho;

                    return x;
                }
            }

            


            /// <summary>
            /// Driver routine
            /// </summary>
            public double[] Krylov(CoordinateVector SolutionVec, double[] currentX, double[] f0, out double errstep) {
                //this.m_AssembleMatrix(out OpMtxRaw, out OpAffineRaw, out MassMtxRaw, SolutionVec.Mapping.Fields.ToArray());
                
                ThisRunFirstIter = ThisLevelIterations;
                Converged = false;
                rho0 = double.NaN;


                double[] step = Solve(SolutionVec, currentX, f0, new double[currentX.Length], out errstep);
                int kinn = 0;
                Console.WriteLine("Error Krylov:   " + errstep);




                //while(kinn < restart_limit && errstep > GMRESConvCrit) {
                while(!Termination(errstep)) { 
                    kinn++;

                    step = Solve(SolutionVec, currentX, f0, step, out errstep);

                    Console.WriteLine("Error Krylov:   " + errstep);
                }

                return step;

            }

            /// <summary>
            /// Finite difference directional derivative Approximate f'(x) w
            /// C.T.Kelley, April 1, 2003
            /// This code comes with no guarantee or warranty of any kind.
            /// </summary>
            /// <param name="SolutionVec">Solution point</param>
            /// <param name="w">Direction</param>
            /// <param name="f0">f0, usually has been calculated earlier</param>
            /// <param name="linearization">True if the Operator should be linearized and evaluated afterwards</param>
            /// <param name="currentX">linearization point</param>
            /// <returns></returns>
            double[] dirder(CoordinateVector SolutionVec, double[] currentX, double[] w, double[] f0, bool linearization = false) {
                using(var tr = new FuncTrace()) {
                    double epsnew = 1E-7;

                    int n = SolutionVec.Length;
                    double[] fx = new double[f0.Length];

                    // Scale the step
                    if(w.L2NormPow2().MPISum().Sqrt() == 0) {
                        fx.Clear();
                        return fx;
                    }

                    var normw = w.L2NormPow2().MPISum().Sqrt();

                    double xs = GenericBlas.InnerProd(currentX, w).MPISum() / normw;

                    if(xs != 0) {
                        epsnew = epsnew * Math.Max(Math.Abs(xs), 1) * Math.Sign(xs);
                    }
                    epsnew = epsnew / w.L2NormPow2().MPISum().Sqrt();

                    var del = currentX.CloneAs();

                    del.AccV(epsnew, w);

                    double[] temp = new double[SolutionVec.Length];

                    temp.CopyEntries(SolutionVec);

                    owner.CurrentLin.TransformSolFrom(SolutionVec, del);

                    // Just evaluate linearized operator
                    //var OpAffineRaw = this.LinearizationRHS.CloneAs();
                    //this.CurrentLin.OperatorMatrix.SpMV(1.0, new CoordinateVector(SolutionVec.Mapping.Fields.ToArray()), 1.0, OpAffineRaw);
                    //CurrentLin.TransformRhsInto(OpAffineRaw, fx);
                    if(linearization == false) {
                        owner.EvaluateOperator(1.0, SolutionVec.Mapping.Fields, fx, this.HomotopyValue);
                    }
                    //else {
                    //    this.m_AssembleMatrix(out OpMtxRaw, out OpAffineRaw, out MassMtxRaw, SolutionVec.Mapping.Fields.ToArray(), true);
                    //    OpMtxRaw.SpMV(1.0, new CoordinateVector(SolutionVec.Mapping.Fields.ToArray()), 1.0, OpAffineRaw);
                    //    CurrentLin.TransformRhsInto(OpAffineRaw, fx);
                    //}

                    SolutionVec.CopyEntries(temp);

                    // (f1 - f0) / epsnew
                    fx.AccV(1, f0);
                    fx.ScaleV(1 / epsnew);

                    return fx;

                }

            }

            /// <summary>
            /// Apply a sequence of k Givens rotations, used within gmres codes.
            /// C.T.Kelley, April 1, 2003
            /// This code comes with no guarantee or warranty of any kind.
            /// </summary>
            /// <param name="c"></param>
            /// <param name="s"></param>
            /// <param name="vin"></param>
            /// <param name="k"></param>
            /// <returns></returns>
            static double[] givapp(double[] c, double[] s, double[] vin, int k) {
                double[] vrot = vin;
                double w1, w2;

                for(int i = 1; i < k; i++) {
                    w1 = c[i - 1] * vrot[i - 1] - s[i - 1] * vrot[i];
                    w2 = s[i - 1] * vrot[i - 1] + c[i - 1] * vrot[i];
                    vrot[i - 1] = w1;
                    vrot[i] = w2;
                }
                return vrot;
            }

            /// <summary>
            /// Compute the Givens rotation matrix parameters for a and b.
            /// </summary>
            static void rotmat(out double c, out double s, double a, double b) {
                double temp;
                if(b == 0.0) {
                    c = 1.0;
                    s = 0.0;
                } else if(Math.Abs(b) > Math.Abs(a)) {
                    temp = a / b;
                    s = 1.0 / Math.Sqrt(1.0 + temp * temp);
                    c = temp * s;
                } else {
                    temp = b / a;
                    c = 1.0 / Math.Sqrt(1.0 + temp * temp);
                    s = temp * c;
                }
            }

            /// <summary>
            /// init
            /// </summary>
            public void Init(MultigridOperator op) {
                var baseMapping = op.BaseGridProblemMapping;
                DGField[] ff = new DGField[baseMapping.BasisS.Count];
                for(int i = 0; i < ff.Length; i++) {
                    var b = baseMapping.BasisS[i];
                    if(b is XDGBasis xb) {
                        ff[i] = new XDGField(xb);
                    } else {
                        ff[i] = new SinglePhaseField(b);
                    }
                }
                OrgDG = new CoordinateVector(ff);

                if(PrecondInstance != null) {
                    PrecondInstance.Dispose();
                    PrecondInstance = null;
                }

                if(owner.PrecondConfig != null)
                    PrecondInstance = owner.PrecondConfig.CreateInstance(op);
            }

            /// <summary>
            /// Instance of <see cref="NonlinearSolver.PrecondConfig"/>
            /// </summary>
            ISolverSmootherTemplate PrecondInstance;

            CoordinateVector OrgDG; 

            /// <summary>
            /// Solution routine as defined by interface
            /// </summary>
            public void Solve<U, V>(U X, V B)
                where U : IList<double>
                where V : IList<double> //
            {
                double[] sol = this.Krylov(OrgDG, X.ToArray(), B.ToArray(), out double errstep);
                X.SetV(sol, 1.0);
            }

            /// <summary>
            /// <see cref="ISolverSmootherTemplate.ResetStat"/>
            /// </summary>
            public void ResetStat() {
                ThisLevelIterations = 0;
                Converged = false;
                this.ThisRunFirstIter = 0;
            }

            public object Clone() {
                throw new NotImplementedException();
            }

            public void Dispose() {
                if(PrecondInstance != null)
                    PrecondInstance.Dispose();
                PrecondInstance = null;
            }

            public long UsedMemory() {
                throw new NotImplementedException();
            }
        }



        /// <summary>
        /// Apply three-point safeguarded parabolic model for a line search.
        /// C.T.Kelley, April 1, 2003
        /// This code comes with no guarantee or warranty of any kind.
        /// function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
        /// input:
        ///        lambdac = current steplength
        ///        lambdam = previous steplength
        ///        ff0 = value of \| F(x_c) \|^2
        ///        ffc = value of \| F(x_c + \lambdac d) \|^2
        ///        ffm = value of \| F(x_c + \lambdam d) \|^2
        ///
        /// output:
        /// lambdap = new value of lambda given parabolic model
        ///
        /// internal parameters:
        /// sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
        /// </summary>
        /// <param name="lambdac"></param>
        /// <param name="lambdam"></param>
        /// <param name="ff0"></param>
        /// <param name="ffc"></param>
        /// <param name="ffm"></param>
        /// <returns></returns>
        static double parab3p(double lambdac, double lambdam, double ff0, double ffc, double ffm) {
            double sigma0 = 0.1;
            double sigma1 = 0.5;

            double c2 = lambdam * (ffc - ff0) - lambdac * (ffm - ff0);
            if(c2 >= 0)
                return sigma1 * lambdac;
            double c1 = lambdac * lambdac * (ffm - ff0) - lambdam * lambdam * (ffc - ff0);
            double lambdap = -c1 * 0.5 / c2;
            if(lambdap < sigma0 * lambdac) lambdap = sigma0 * lambdac;
            if(lambdap > sigma1 * lambdac) lambdap = sigma1 * lambdac;

            return lambdap;
        }
    }
}

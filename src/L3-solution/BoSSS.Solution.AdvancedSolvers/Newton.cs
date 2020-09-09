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
using System.Security.Cryptography.X509Certificates;

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// Implementation based on presudocode from Kelley, C.
    /// Solving Nonlinear Equations with Newton’s Method. Fundamentals of Algorithms.
    /// Society for Industrial and Applied Mathematics, 2003. https://doi.org/10.1137/1.9780898718898.
    /// </summary>
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
        /// Maximum number of GMRES(m) restarts
        /// </summary>
        public int restart_limit = 1000;

        /// <summary>
        /// Number of iterations, where Jacobi is not updated. Also known as constant newton method. Default 1, means regular newton.
        /// </summary>
        public int constant_newton_it = 1;

        /// <summary>
        /// Maximum dimension of the krylov subspace. Equals m in GMRES(m)
        /// </summary>
        public int maxKrylovDim = 30;

        /// <summary>
        /// Convergence criterion for nonlinear iteration
        /// </summary>
        public double ConvCrit = 1e-6;

        /// <summary>
        /// Maximum number of step-length iterations
        /// </summary>
        public double maxStep = 30;

        /// <summary>
        /// Convergence for Krylov and GMRES iterations
        /// </summary>
        public double GMRESConvCrit = 1e-6;

        /// <summary>
        /// Options for (approximate) solution to the linearizes system
        /// </summary>
        public enum ApproxInvJacobianOptions {

            /// <summary>
            /// Using a matrix-free GMRES implementation, optionally employing <see cref="Newton.linsolver"/> (typically not matrix-free) as a preconditioner.
            /// </summary>
            MatrixFreeGMRES = 1,

            /// <summary>
            /// Using the solver <see cref="Newton.linsolver"/> for computing Newton corrections
            /// </summary>
            ExternalSolver = 2
        }

        /// <summary>
        /// Options for (approximate) solution to the linearizes system
        /// </summary>
        public ApproxInvJacobianOptions ApproxJac = ApproxInvJacobianOptions.MatrixFreeGMRES;

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
            LineSearch = 2
        }

        /// <summary>
        /// Options for Globalization, i.e. means to ensure convergence of Newton iterations
        /// when the initial guess is far away from the solution.
        /// </summary>
        public GlobalizationOption Globalization = GlobalizationOption.LineSearch;

        public string m_SessionPath;

        /// <summary>
        /// Optional external solver for computation of Newton corrections
        /// </summary>
        public ISolverSmootherTemplate linsolver;

        public bool UsePresRefPoint;

        /// <summary>
        /// Prints the step reduction factor
        /// </summary>
        public bool printLambda = false;


        /// <summary>
        /// Main solver routine
        /// </summary>
        public override void SolverDriver<S>(CoordinateVector SolutionVec, S RHS) {

            using (var tr = new FuncTrace()) {


                int itc = 0;
                double[] CurSol, // "current (approximate) solution", i.e.
                    CurRes; // residual associated with 'CurSol'
                //TempRes;
                double rat;


                // Eval_F0

                using (new BlockTrace("Slv Init", tr)) {
                    base.Init(SolutionVec, RHS, out CurSol, out CurRes);
                };

                this.CurrentLin.TransformSolFrom(SolutionVec, CurSol);
                EvaluateOperator(1, SolutionVec.Mapping.ToArray(), CurRes);

                Console.WriteLine("Residual base.init:   " + CurRes.L2NormPow2().MPISum().Sqrt());
                //base.EvalResidual(x, ref f0);

                // fnorm
                double fnorm = CurRes.L2NormPow2().MPISum().Sqrt();
                double fNormo = 1;
                double errstep;
                double[] step = new double[CurSol.Length];
                double TrustRegionDelta = -1; // only used for dogleg (aka Trust-Region) method
                //double[] stepOld = new double[CurSol.Length];
                //BlockMsrMatrix CurrentJac;
                //bool secondCriteriumConverged = false;
                OnIterationCallback(itc, CurSol.CloneAs(), CurRes.CloneAs(), this.CurrentLin);
                double fnorminit = fnorm;
                using (new BlockTrace("Slv Iter", tr)) {
                    while ((fnorm > ConvCrit * fnorminit + ConvCrit
                        //&& secondCriteriumConverged == false
                        && itc < MaxIter)
                        || itc < MinIter) {
                        rat = fnorm / fNormo;
                        //if (Math.Abs(fNormo - fnorm) < 1e-12)
                        //    break;
                        fNormo = fnorm;
                        itc++;

                        // computation of Newton step
                        // --------------------------


                        // How should the inverse of the Jacobian be approximated?
                        if (ApproxJac == ApproxInvJacobianOptions.MatrixFreeGMRES) {
                            // ++++++++++++++++++++++++++
                            // Option: Matrix-Free GMRES
                            // ++++++++++++++++++++++++++

                            if (Precond != null) {
                                Precond.Init(CurrentLin);
                            }
                            //base.EvalResidual(x, ref f0);
                            step = Krylov(SolutionVec, CurSol, CurRes, out errstep);
                            step.ScaleV(-1);
                        } else if (ApproxJac == ApproxInvJacobianOptions.ExternalSolver) {
                            // +++++++++++++++++++++++++++++
                            // Option: use 'external' solver
                            // +++++++++++++++++++++++++++++



                            var solver = this.linsolver;
                            solver.Init(CurrentLin);
                            step.ClearEntries();
                            var check = CurRes.CloneAs();
                            solver.ResetStat();

                            if (solver is IProgrammableTermination pt) {
                                // iterative solver with programmable termination is used -

                                double f0_L2 = CurRes.MPI_L2Norm();
                                double thresh = f0_L2 * 1e-5;
                                Console.WriteLine($"Inexact Newton: setting convergence threshold to {thresh:0.##E-00}");
                                pt.TerminationCriterion = (iter, R0_l2, R_l2) => {
                                    return (R_l2 > thresh) && (iter < 100);
                                };


                            }

                            solver.Solve(step, CurRes);
                            step.ScaleV(-1);

                        } else {
                            throw new NotImplementedException($"approximation option {ApproxJac} for the Jacobian seems not to be existent.");
                        }

                        // globalization
                        // -------------
                        switch (Globalization) {
                            case GlobalizationOption.Dogleg:
                                DogLeg(SolutionVec, CurSol, CurRes, step, itc, ref TrustRegionDelta);
                                break;

                            case GlobalizationOption.LineSearch:
                                LineSearch(SolutionVec, CurSol, CurRes, step);
                                break;

                            default:
                                throw new NotImplementedException();
                        }

                        // fix the pressure
                        // ----------------
                        if (UsePresRefPoint == false) {

                            if (SolutionVec.Mapping.Fields[2] is XDGField Xpres) {
                                DGField presSpA = Xpres.GetSpeciesShadowField("A");
                                DGField presSpB = Xpres.GetSpeciesShadowField("B");
                                var meanpres = presSpB.GetMeanValueTotal(null);
                                presSpA.AccConstant(-1.0 * meanpres);
                                presSpB.AccConstant(-1.0 * meanpres);
                            } else {
                                DGField pres = SolutionVec.Mapping.Fields[2];
                                var meanpres = pres.GetMeanValueTotal(null);
                                pres.AccConstant(-1.0 * meanpres);
                            }
                        }



                        // update linearization
                        // --------------------
                        if (itc % constant_newton_it == 0) {
                            //base.UpdateLinearization(SolutionVec.Mapping.Fields);
                            base.Update(SolutionVec.Mapping.Fields, ref CurSol);
                            if (constant_newton_it != 1) {
                                Console.WriteLine("Jacobian is updated: it {0}", itc);
                            }
                        }

                        // residual evaluation & callback
                        // ------------------------------
                        EvaluateOperator(1, SolutionVec.Mapping.Fields, CurRes);
                        fnorm = CurRes.MPI_L2Norm();


                        OnIterationCallback(itc, CurSol.CloneAs(), CurRes.CloneAs(), this.CurrentLin);

                        #region second criterium
                        // Just testing. According to "Pawlowski et al. - 2006 - Globalization Techniques for Newton–Krylov Methods"
                        // this criterium is useful to "ensure that even finer physical details of the flow and are resolved"


                        //double[] WMat_s = new double[x.Length];
                        //double psi_r = 1e-3;
                        //double psi_a = 1e-8;
                        //double[] truestep = step;
                        //truestep.ScaleV(lambda);
                        //for(int i = 0; i < WMat_s.Length; i++) {
                        //    WMat_s[i] = 1 / (psi_r * x[i] + psi_a) * truestep[i];
                        //}
                        //WMat_s.CheckForNanOrInfV();

                        //double secondCriterium = WMat_s.L2Norm()/x.Length;

                        //if(secondCriterium < 1) {
                        //    secondCriteriumConverged = true;
                        //}
                        ////Console.WriteLine("Norm Of the second criterium {0}", secondCriterium);
                        ////if((fnorm < ConvCrit * fnorminit * 0 + ConvCrit))
                        ////    Console.WriteLine("Criterium 1 fulfilled");
                        ////if((secondCriteriumConverged == true))
                        ////    Console.WriteLine("Criterium 2 fulfilled");
                        ////if(itc > MaxIter)
                        ////    Console.WriteLine("Criterium 3 fulfilled");

                        #endregion

                    }
                }


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
        /// <returns>
        /// Updated Solution
        /// </returns>
        private void LineSearch(CoordinateVector SolutionVec, double[] CurSol, double[] CurRes, double[] step) {

            double[] TempSol;
            double[] TempRes;

            double alpha = 1E-4;
            double sigma1 = 0.5;


            // Start line search
            double lambda = 1;
            double lamm = 1;
            double lamc = lambda;
            double iarm = 0;
            TempSol = CurSol.CloneAs();
            TempSol.AccV(lambda, step);
            this.CurrentLin.TransformSolFrom(SolutionVec, TempSol);

            TempRes = new double[TempSol.Length];
            EvaluateOperator(1, SolutionVec.Mapping.Fields, TempRes);

            double nft = TempRes.L2NormPow2().MPISum().Sqrt();
            double nf0 = CurRes.L2NormPow2().MPISum().Sqrt();
            double ff0 = nf0 * nf0; // residual norm of current solution ^2
            double ffc = nft * nft; //
            double ffm = nft * nft; //


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
                EvaluateOperator(1, SolutionVec.Mapping.Fields, TempRes);

                nft = TempRes.L2NormPow2().MPISum().Sqrt();
                ffm = ffc;
                ffc = nft * nft;
                iarm++;

                if (printLambda)
                    Console.WriteLine("    Residuum:  " + nft + " lambda = " + lambda);

            }
            // transform solution back to 'original domain'
            // to perform the linearization at the new point...
            // (and for Level-Set-Updates ...)
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
        /// <remarks>
        /// See:
        /// - Pawlowski et. al., 2006, Globalization Techniques for Newton–Krylov Methods and Applications to the Fully Coupled Solution of the Navier–Stokes Equations, SIAM Review, Vol. 48, No. 4, pp 700-721.
        /// - Pawlowski et. al., 2008, Inexact Newton Dogleg Methods, SIAM Journal on Numerical Analysis, Vol. 46, No. 4, pp 2112-2132.
        /// </remarks>
        private void DogLeg(CoordinateVector SolutionVec, double[] CurSol, double[] CurRes, double[] stepIN, int NewtonIterCnt, ref double TrustRegionDelta) {

            // algorithm constants taken from [Pawlowski et. al. 2006]
            // =======================================================
            const double t = 1.0e-4;
            const double delta_min = 1e-6;
            const double delta_max = 1e10;

            // initial estimate of trust region width in first iteration
            // =========================================================
            if (NewtonIterCnt < 1)
                throw new ArgumentException();
            if (NewtonIterCnt == 1) {
                double norm_step = stepIN.MPI_L2Norm();
                if (norm_step < delta_min)
                    TrustRegionDelta = 2 * delta_min;
                else
                    TrustRegionDelta = norm_step;

                TrustRegionDelta = Math.Min(delta_max, TrustRegionDelta);
            }


            // TODO change later: always use a very small initial trust region!
            // TrustRegionDelta = TrustRegionDelta / 10;

            if (TrustRegionDelta < delta_min || TrustRegionDelta > delta_max)
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

                double lambda = -a[0] / a[1];

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
                if (l2_stepIN <= _TrustRegionDelta) {
                    // use Newton Step
                    Console.WriteLine($"       -------- using Newton step (delta = {_TrustRegionDelta})");
                    step.SetV(stepIN);
                } else {
                    if (l2_stepCP < _TrustRegionDelta) {
                        // interpolate between Cauchy-point and Newton-step
                        Console.WriteLine($"       -------- between Cauchy point and Newton step (delta = {_TrustRegionDelta})");
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
                            if (!(tau >= -0.00001 && tau <= 1.00001) || tau.IsNaN() || tau.IsInfinity())
                                throw new ArithmeticException();
                        }
                        // do interpolation
                        step.SetV(stepCP, (1 - tau));
                        step.AccV(tau, stepIN);
                        Debug.Assert(Math.Abs((step.MPI_L2Norm() / _TrustRegionDelta) - 1.0) <= 1.0e-3, "interpolation step went wrong");
                    } else {
                        // use reduced Cauchy point
                        Console.WriteLine($"       -------- using Cauchy point (delta = {_TrustRegionDelta})");
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
                base.EvaluateOperator(1.0, SolutionVec.Fields, temp);

                return l2_CurRes - temp.MPI_L2Norm();
            }

            // trust region adaptation loop
            double last_ared = ared();
            double last_pred = pred();
            while (last_ared < t * last_pred) {
                double newTrustRegionDelta = TrustRegionDelta * 0.5;
                if (newTrustRegionDelta <= delta_min)
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

                if (last_ared / last_pred < rho_s && l2_stepIN < TrustRegionDelta) {
                    TrustRegionDelta = Math.Max(l2_stepIN, delta_min);
                } else if (last_ared / last_pred < rho_s) {
                    TrustRegionDelta = Math.Max(beta_s * TrustRegionDelta, delta_min); // shrinking
                } else if (last_ared / last_pred > rho_e) {
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
        /// Preconditioned GMRES, using <see cref="NonlinearSolver.Precond"/> as a preconditioner
        /// </summary>
        /// <param name="SolutionVec">Current Point</param>
        /// <param name="f0">Function at current point</param>
        /// <param name="xinit">initial iterate</param>
        /// <param name="errstep">error of step</param>
        /// <param name="currentX"></param>
        /// <returns></returns>
        double[] MatrixFreeGMRES(CoordinateVector SolutionVec, double[] currentX, double[] f0, double[] xinit, out double errstep) {
            using (var tr = new FuncTrace()) {
                int n = f0.Length;

                int reorth = 1; // Orthogonalization method -> 1: Brown/Hindmarsh condition, 3: Always reorthogonalize

                // RHS of the linear equation system
                double[] b = new double[n];
                b.AccV(1, f0);

                double[] x = new double[n];
                double[] r = new double[n];

                int Nloc = base.CurrentLin.OperatorMatrix.RowPartitioning.LocalLength;
                int Ntot = base.CurrentLin.OperatorMatrix.RowPartitioning.TotalLength;

                r = b;

                //Initial solution
                if (xinit.L2Norm() != 0) {
                    x = xinit.CloneAs();
                    r.AccV(-1, dirder(SolutionVec, currentX, x, f0));
                }
                // Precond = null;
                if (Precond != null) {
                    var temp2 = r.CloneAs();
                    r.ClearEntries();
                    //this.OpMtxRaw.InvertBlocks(OnlyDiagonal: false, Subblocks: false).SpMV(1, temp2, 0, r);
                    Precond.Solve(r, temp2);
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

                Console.WriteLine("Error NewtonGMRES:   " + rho);

                // Termination of entry
                if (rho < GMRESConvCrit)
                    return SolutionVec.ToArray();

                V[0].SetV(r, alpha: (1.0 / rho));
                double beta = rho;
                int k = 1;

                while ((rho > GMRESConvCrit) && k <= m) {
                    V[k].SetV(dirder(SolutionVec, currentX, V[k - 1], f0));
                    //CurrentLin.OperatorMatrix.SpMV(1.0, V[k-1], 0.0, temp3);
                    // Call directional derivative
                    //V[k].SetV(f0);

                    if (Precond != null) {
                        var temp3 = V[k].CloneAs();
                        V[k].ClearEntries();
                        //this.OpMtxRaw.InvertBlocks(false,false).SpMV(1, temp3, 0, V[k]);
                        Precond.Solve(V[k], temp3);
                    }

                    double normav = V[k].L2NormPow2().MPISum().Sqrt();

                    // Modified Gram-Schmidt
                    for (int j = 1; j <= k; j++) {
                        H[j - 1, k - 1] = GenericBlas.InnerProd(V[k], V[j - 1]).MPISum();
                        V[k].AccV(-H[j - 1, k - 1], V[j - 1]);
                    }
                    H[k, k - 1] = V[k].L2NormPow2().MPISum().Sqrt();
                    double normav2 = H[k, k - 1];


                    // Reorthogonalize ?
                    if ((reorth == 1 && Math.Round(normav + 0.001 * normav2, 3) == Math.Round(normav, 3)) || reorth == 3) {
                        for (int j = 1; j <= k; j++) {
                            double hr = GenericBlas.InnerProd(V[k], V[j - 1]).MPISum();
                            H[j - 1, k - 1] = H[j - 1, k - 1] + hr;
                            V[k].AccV(-hr, V[j - 1]);
                        }
                        H[k, k - 1] = V[k].L2NormPow2().MPISum().Sqrt();
                    }

                    // Watch out for happy breakdown
                    if (H[k, k - 1] != 0)
                        V[k].ScaleV(1 / H[k, k - 1]);



                    // Form and store the information for the new Givens rotation
                    //if (k > 1) {
                    //    // for (int i = 1; i <= k; i++) {
                    //    H.SetColumn(k - 1, givapp(c.GetSubVector(0, k - 1), s.GetSubVector(0, k - 1), H.GetColumn(k - 1), k - 1));
                    //    //}
                    //}

                    // Givens rotation from SoftGMRES
                    double temp;
                    for (int l = 1; l <= k - 1; l++) {
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
                    if (nu != 0) {
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
                for (int ii = 0; ii < k; ii++) {
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
        double[] Krylov(CoordinateVector SolutionVec, double[] currentX, double[] f0, out double errstep) {
            //this.m_AssembleMatrix(out OpMtxRaw, out OpAffineRaw, out MassMtxRaw, SolutionVec.Mapping.Fields.ToArray());
            double[] step = MatrixFreeGMRES(SolutionVec, currentX, f0, new double[currentX.Length], out errstep);
            int kinn = 0;
            Console.WriteLine("Error Krylov:   " + errstep);

            while (kinn < restart_limit && errstep > GMRESConvCrit) {
                kinn++;

                step = MatrixFreeGMRES(SolutionVec, currentX, f0, step, out errstep);

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
        /// <returns></returns>
        public double[] dirder(CoordinateVector SolutionVec, double[] currentX, double[] w, double[] f0, bool linearization = false) {
            using (var tr = new FuncTrace()) {
                double epsnew = 1E-7;

                int n = SolutionVec.Length;
                double[] fx = new double[f0.Length];

                // Scale the step
                if (w.L2NormPow2().MPISum().Sqrt() == 0) {
                    fx.Clear();
                    return fx;
                }

                var normw = w.L2NormPow2().MPISum().Sqrt();

                double xs = GenericBlas.InnerProd(currentX, w).MPISum() / normw;

                if (xs != 0) {
                    epsnew = epsnew * Math.Max(Math.Abs(xs), 1) * Math.Sign(xs);
                }
                epsnew = epsnew / w.L2NormPow2().MPISum().Sqrt();

                var del = currentX.CloneAs();

                del.AccV(epsnew, w);

                double[] temp = new double[SolutionVec.Length];

                temp.CopyEntries(SolutionVec);

                this.CurrentLin.TransformSolFrom(SolutionVec, del);

                // Just evaluate linearized operator
                //var OpAffineRaw = this.LinearizationRHS.CloneAs();
                //this.CurrentLin.OperatorMatrix.SpMV(1.0, new CoordinateVector(SolutionVec.Mapping.Fields.ToArray()), 1.0, OpAffineRaw);
                //CurrentLin.TransformRhsInto(OpAffineRaw, fx);
                if (linearization == false) {
                    EvaluateOperator(1.0, SolutionVec.Mapping.Fields, fx);
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

            for (int i = 1; i < k; i++) {
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
            if (b == 0.0) {
                c = 1.0;
                s = 0.0;
            } else if (Math.Abs(b) > Math.Abs(a)) {
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
            if (c2 >= 0)
                return sigma1 * lambdac;
            double c1 = lambdac * lambdac * (ffm - ff0) - lambdam * lambdam * (ffc - ff0);
            double lambdap = -c1 * 0.5 / c2;
            if (lambdap < sigma0 * lambdac) lambdap = sigma0 * lambdac;
            if (lambdap > sigma1 * lambdac) lambdap = sigma1 * lambdac;

            return lambdap;
        }


    }

}

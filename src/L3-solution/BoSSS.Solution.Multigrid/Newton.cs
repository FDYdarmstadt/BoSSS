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

namespace BoSSS.Solution.Multigrid {
    /// <summary>
    /// Implementation based on presudocode from Kelley, C. 
    /// Solving Nonlinear Equations with Newton’s Method. Fundamentals of Algorithms. 
    /// Society for Industrial and Applied Mathematics, 2003. https://doi.org/10.1137/1.9780898718898.
    /// </summary>
    public class Newton : NonlinearSolver {
        public Newton(AssembleMatrixDel __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq, MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig) :
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
        public int restart_limit = 10;


        /// <summary>
        /// Maximum dimension of the krylov subspace. Equals m in GMRES(m)
        /// </summary>
        public int maxKrylovDim = 30;

        /// <summary>
        /// Convergence criterium for nonlinear iteration
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

        public CoordinateVector m_SolutionVec;

        public enum ApproxInvJacobianOptions { GMRES = 1, DirectSolver = 2, DirectSolverHybrid = 3, DirectSolverOpMatrix =4 }

        public ApproxInvJacobianOptions ApproxJac = ApproxInvJacobianOptions.DirectSolverOpMatrix;

        public MsrMatrix currentPrecMatrix = null;

        public string m_SessionPath;


        //bool solveVelocity = true;

        //double VelocitySolver_ConvergenceCriterion = 1e-5;

        //double StressSolver_ConvergenceCriterion = 1e-5;

        public override void SolverDriver<S>(CoordinateVector SolutionVec, S RHS) {

            using (var tr = new FuncTrace()) {

                m_SolutionVec = SolutionVec;

                int itc;
                itc = 0;
                double[] x, xt, xOld, f0, deltaX, ft;
                double rat;
                double alpha = 1E-4;
                //double sigma0 = 0.1;
                double sigma1 = 0.5;
                //double maxarm = 20;
                //double gamma = 0.9;

                // Eval_F0 
                using (new BlockTrace("Slv Init", tr)) {
                    base.Init(SolutionVec, RHS, out x, out f0);
                };

                Console.WriteLine("Residual base.init:   " + f0.L2NormPow2().MPISum().Sqrt());


                deltaX = new double[x.Length];
                xt = new double[x.Length];
                ft = new double[x.Length];


                this.CurrentLin.TransformSolFrom(SolutionVec, x);
                base.EvalResidual(x, ref f0);


                // fnorm
                double fnorm = f0.L2NormPow2().MPISum().Sqrt();
                double fNormo = 1;
                double errstep;
                double[] step = new double[x.Length];
                double[] stepOld = new double[x.Length];
                MsrMatrix CurrentJac;

                Console.WriteLine("Start residuum for nonlinear iteration:  " + fnorm);

                OnIterationCallback(itc, x.CloneAs(), f0.CloneAs(), this.CurrentLin);

                //int[] Velocity_idx = SolutionVec.Mapping.GetSubvectorIndices(false, 0, 1, 2);
                //int[] Stresses_idx = SolutionVec.Mapping.GetSubvectorIndices(false, 3, 4, 5);

                //int[] Velocity_fields = new int[] { 0, 1, 2 };
                //int[] Stress_fields = new int[] { 3, 4, 5 };

                //int NoCoupledIterations = 1;

                using (new BlockTrace("Slv Iter", tr)) {
                    while (fnorm > ConvCrit && itc < MaxIter) {
                        rat = fnorm / fNormo;
                        fNormo = fnorm;
                        itc++;

                        // How should the inverse of the Jacobian be approximated?
                        if (ApproxJac == ApproxInvJacobianOptions.GMRES) {
                            CurrentJac = new MsrMatrix(x.Length);
                            if (Precond != null) {
                                Precond.Init(CurrentLin);
                            }
                            step = Krylov(SolutionVec, x, f0, out errstep);
                        } else if (ApproxJac == ApproxInvJacobianOptions.DirectSolver) {
                            CurrentJac = diffjac(SolutionVec, x, f0);
                            CurrentJac.SaveToTextFileSparse("Jacobi");
                            CurrentLin.OperatorMatrix.SaveToTextFileSparse("OpMatrix");
                            var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
                            solver.DefineMatrix(CurrentJac);
                            step.ClearEntries();
                            solver.Solve(step, f0);

                        }
                        else if (ApproxJac == ApproxInvJacobianOptions.DirectSolverHybrid) {
                            //EXPERIMENTAL_____________________________________________________________________
                            MultidimensionalArray OpMatrixMatl = MultidimensionalArray.Create(x.Length, x.Length);
                            CurrentJac = diffjac(SolutionVec, x, f0);
                            //Console.WriteLine("Calling MATLAB/Octave...");
                            using (BatchmodeConnector bmc = new BatchmodeConnector())
                            {
                                bmc.PutSparseMatrix(CurrentJac, "Jacobi");
                                bmc.PutSparseMatrix(CurrentLin.OperatorMatrix, "OpMatrix");
                                bmc.Cmd("Jacobi(abs(Jacobi) < 10^-6)=0; dim = length(OpMatrix);");
                                bmc.Cmd("dim = length(OpMatrix);");
                                bmc.Cmd(@"for i=1:dim

                                if (i >= 16) && (i <= 33) && (i + 6 <= 33) && (i + 12 <= 33)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end

                                if (i >= 49) && (i <= 66) && (i + 6 <= 66) && (i + 12 <= 66)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end

                                if (i >= 82) && (i <= 99) && (i + 6 <= 99) && (i + 12 <= 99)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end

                                if (i >= 115) && (i <= 132) && (i + 6 <= 132) && (i + 12 <= 132)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end

                                if (i >= 148) && (i <= 165) && (i + 6 <= 165) && (i + 12 <= 165)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end 

                                if (i >= 181) && (i <= 198) && (i + 6 <= 198) && (i + 12 <= 198)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end  

                                if (i >= 214) && (i <= 231) && (i + 6 <= 231) && (i + 12 <= 231)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 247) && (i <= 264) && (i + 6 <= 264) && (i + 12 <= 264)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 280) && (i <= 297) && (i + 6 <= 297) && (i + 12 <= 297)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 313) && (i <= 330) && (i + 6 <= 330) && (i + 12 <= 330)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 346) && (i <= 363) && (i + 6 <= 363) && (i + 12 <= 363)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 379) && (i <= 396) && (i + 6 <= 396) && (i + 12 <= 396)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 412) && (i <= 429) && (i + 6 <= 429) && (i + 12 <= 429)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 445) && (i <= 462) && (i + 6 <= 462) && (i + 12 <= 462)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 478) && (i <= 495) && (i + 6 <= 495) && (i + 12 <= 495)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                if (i >= 511) && (i <= 528) && (i + 6 <= 528) && (i + 12 <= 528)
                                    OpMatrix(i, i) = Jacobi(i, i);
                                    OpMatrix(i, i + 6) = Jacobi(i, i + 6);
                                    OpMatrix(i + 6, i) = Jacobi(i + 6, i);
                                    OpMatrix(i + 12, i + 6) = Jacobi(i + 12, i + 6);
                                    OpMatrix(i + 6, i + 12) = Jacobi(i + 6, i + 12);
                                end
                                end");
                                bmc.Cmd("OpMatrix = full(OpMatrix)");
                                bmc.GetMatrix(OpMatrixMatl, "OpMatrix");
                                bmc.Execute(false);
                            }

                            MsrMatrix OpMatrix = OpMatrixMatl.ToMsrMatrix();
                            var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();

                            //Console.WriteLine("USING HIGH EXPERIMENTAL OPMATRIX WITH JAC! only for p=1, GridLevel=2");
                            solver.DefineMatrix(OpMatrix);
                            //______________________________________________________________________________________________________

                            step.ClearEntries();
                            solver.Solve(step, f0);


                        } else if (ApproxJac == ApproxInvJacobianOptions.DirectSolverOpMatrix) {
                            CurrentJac = new MsrMatrix(x.Length);
                            CurrentJac = CurrentLin.OperatorMatrix.ToMsrMatrix();
                            var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
                            solver.DefineMatrix(CurrentJac);
                            step.ClearEntries();
                            solver.Solve(step, f0);


                        } else {
                            throw new NotImplementedException("Your approximation option for the jacobian seems not to be existent.");
                        }

                        //if (itc > NoCoupledIterations)
                        //{
                        //    if (solveVelocity)
                        //    {
                        //        Console.WriteLine("stress correction = 0");
                        //        foreach (int idx in Stresses_idx)
                        //        {
                        //            step[idx] = 0.0;
                        //        }
                        //    }
                        //    else
                        //    {
                        //        Console.WriteLine("velocity correction = 0");
                        //        foreach (int idx in Velocity_idx)
                        //        {
                        //            step[idx] = 0.0;
                        //        }
                        //    }
                        //}

                        // Start line search
                        xOld = x;
                        double lambda = 1;
                        double lamm = 1;
                        double lamc = lambda;
                        double iarm = 0;
                        xt = x.CloneAs();
                        xt.AccV(lambda, step);
                        this.CurrentLin.TransformSolFrom(SolutionVec, xt);
                        EvaluateOperator(1, SolutionVec.Mapping.Fields, ft);
                        var nft = ft.L2NormPow2().MPISum().Sqrt(); var nf0 = f0.L2NormPow2().MPISum().Sqrt(); var ff0 = nf0 * nf0; var ffc = nft * nft; var ffm = nft * nft;

                        // Control of the the step size
                        while (nft >= (1 - alpha * lambda) * nf0 && iarm < maxStep) {

                            // Line search starts here

                            if (iarm == 0)
                                lambda = sigma1 * lambda;
                            else
                                lambda = parab3p(lamc, lamm, ff0, ffc, ffm);

                            // Update x;
                            xt = x.CloneAs();
                            xt.AccV(lambda, step);
                            lamm = lamc;
                            lamc = lambda;

                            this.CurrentLin.TransformSolFrom(SolutionVec, xt);

                            EvaluateOperator(1, SolutionVec.Mapping.Fields, ft);

                            nft = ft.L2NormPow2().MPISum().Sqrt();
                            ffm = ffc;
                            ffc = nft * nft;
                            iarm++;

                            Console.WriteLine("Step size:  " + lambda + "with Residuum:  " + nft);
                        }
                        // transform solution back to 'original domain'
                        // to perform the linearization at the new point...
                        // (and for Level-Set-Updates ...)
                        this.CurrentLin.TransformSolFrom(SolutionVec, xt);

                        // update linearization
                        base.Update(SolutionVec.Mapping.Fields, ref xt);

                        // residual evaluation & callback
                        base.EvalResidual(xt, ref ft);

                        // EvaluateOperator(1, SolutionVec.Mapping.Fields, ft);

                        //base.Init(SolutionVec, RHS, out x, out f0);

                        fnorm = ft.L2NormPow2().MPISum().Sqrt();

                        x = xt;
                        f0 = ft.CloneAs();

                        //if (itc > NoCoupledIterations)
                        //{

                        //    double coupledL2Res = 0.0;
                        //    if (solveVelocity)
                        //    {
                        //        foreach (int idx in Velocity_idx)
                        //        {
                        //            coupledL2Res += f0[idx].Pow2();
                        //        }
                        //    }
                        //    else
                        //    {
                        //        foreach (int idx in Stresses_idx)
                        //        {
                        //            coupledL2Res += f0[idx].Pow2();
                        //        }
                        //    }
                        //    coupledL2Res = coupledL2Res.Sqrt();

                        //    Console.WriteLine("coupled residual = {0}", coupledL2Res);

                        //    if (solveVelocity && coupledL2Res < this.VelocitySolver_ConvergenceCriterion)
                        //    {
                        //        Console.WriteLine("SolveVelocity = false");
                        //        this.solveVelocity = false;
                        //    }
                        //    else if (!solveVelocity && coupledL2Res < this.StressSolver_ConvergenceCriterion)
                        //    {
                        //        Console.WriteLine("SolveVelocity = true");
                        //        this.solveVelocity = true;
                        //    }
                        //}

                        OnIterationCallback(itc, x.CloneAs(), f0.CloneAs(), this.CurrentLin);


                    }
                }

                SolutionVec = m_SolutionVec;

            }
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="SolutionVec">Current Point</param>
        /// <param name="f0">Function at current point</param>
        /// <param name="xinit">initial iterate</param>
        /// <param name="errstep">error of step</param>
        /// <returns></returns>
        public double[] GMRES(CoordinateVector SolutionVec, double[] currentX, double[] f0, double[] xinit, out double errstep) {
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

                    //using (StreamWriter writer = new StreamWriter(m_SessionPath + "//GMRES_Stats.txt", true))
                    //{
                    //writer.WriteLine(k + "   " + rho);
                    //}


                    //Console.WriteLine("Error NewtonGMRES:   " + rho );

                    k++;

                }


                k--;

                //Console.WriteLine("GMRES completed after:   " + k + "steps");

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

        public double[] Krylov(CoordinateVector SolutionVec, double[] currentX, double[] f0, out double errstep) {
            //this.m_AssembleMatrix(out OpMtxRaw, out OpAffineRaw, out MassMtxRaw, SolutionVec.Mapping.Fields.ToArray());
            double[] step = GMRES(SolutionVec, currentX, f0, new double[currentX.Length], out errstep);
            int kinn = 0;
            Console.WriteLine("Error Krylov:   " + errstep);

            while (kinn < restart_limit && errstep > GMRESConvCrit) {
                kinn++;

                step = GMRES(SolutionVec, currentX, f0, step, out errstep);

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
        /// <param name="f0">f0, usally has been calculated earlier</param>
        /// <returns></returns>
        public double[] dirder(CoordinateVector SolutionVec, double[] currentX, double[] w, double[] f0) {
            using (var tr = new FuncTrace()) {
                double epsnew = 1E-7;

                int n = SolutionVec.Length;
                double[] fx = new double[f0.Length];

                // Scale the step
                if (w.L2Norm().MPISum() == 0) {
                    fx.Clear();
                    return fx;
                }

                var normw = w.L2NormPow2().MPISum().Sqrt();

                double xs = GenericBlas.InnerProd(currentX, w).MPISum() / normw;

                if (xs != 0) {
                    epsnew = epsnew * Math.Max(Math.Abs(xs), 1) * Math.Sign(xs);
                }
                epsnew = epsnew / w.L2Norm().MPISum();

                var del = currentX.CloneAs();

                del.AccV(epsnew, w);

                double[] temp = new double[SolutionVec.Length];

                temp.CopyEntries(SolutionVec);

                this.CurrentLin.TransformSolFrom(SolutionVec, del);

                // Just evaluate linearized operator
                //var OpAffineRaw = this.LinearizationRHS.CloneAs();
                //this.CurrentLin.OperatorMatrix.SpMV(1.0, new CoordinateVector(SolutionVec.Mapping.Fields.ToArray()), 1.0, OpAffineRaw);
                //CurrentLin.TransformRhsInto(OpAffineRaw, fx);

                //EvaluateOperator(1.0, SolutionVec.Mapping.Fields, fx);

                this.m_AssembleMatrix(out OpMtxRaw, out OpAffineRaw, out MassMtxRaw, SolutionVec.Mapping.Fields.ToArray());

                OpMtxRaw.SpMV(1.0, new CoordinateVector(SolutionVec.Mapping.Fields.ToArray()), 1.0, OpAffineRaw);

                CurrentLin.TransformRhsInto(OpAffineRaw, fx);

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
        public double[] givapp(double[] c, double[] s, double[] vin, int k) {
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
        public double parab3p(double lambdac, double lambdam, double ff0, double ffc, double ffm) {
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


        BlockMsrMatrix OpMtxRaw, MassMtxRaw;
        double[] OpAffineRaw;

        /// <summary>
        /// Computes a forward difference jacobian and returns the dense jacobian
        /// </summary>
        /// <param name="SolutionVec"></param>
        /// <param name="currentX"></param>
        /// <param name="f0"></param>
        /// <returns></returns>
        public MsrMatrix diffjac(CoordinateVector SolutionVec, double[] currentX, double[] f0) {
            int n = currentX.Length;
            MsrMatrix jac = new MsrMatrix(n);

            var temp = new double[n];

            for (int i = 0; i < n; i++) {
                var zz = new double[n];
                zz[i] = 1;
                temp = dirder(SolutionVec, currentX, zz, f0);
                for (int j = 0; j < n; j++) {
                    jac[j, i] = temp[j];
                }
            }

            return jac;
        }

        /// <summary>
        /// Evaluation of the nonlinear operator.
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="CurrentState">
        /// Current state of DG fields
        /// </param>
        /// <param name="beta">
        /// Pre-scaling of <paramref name="Output"/>.
        /// </param>
        /// <param name="Output"></param>
        //public override void EvaluateOperator(double alpha, IEnumerable<DGField> CurrentState, double[] Output) {
        //    BlockMsrMatrix OpMtxRaw, MassMtxRaw;
        //    double[] OpAffineRaw;
        //    this.m_AssembleMatrix(out OpMtxRaw, out OpAffineRaw, out MassMtxRaw, CurrentState.ToArray());

        //    OpMtxRaw.SpMV(alpha, new CoordinateVector(CurrentState.ToArray()), 1.0, OpAffineRaw);

        //    CurrentLin.TransformRhsInto(OpAffineRaw, Output);

        //    // Inverse of current PrexMatrix
        //    if (currentPrecMatrix != null) {
        //        var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
        //        solver.DefineMatrix(currentPrecMatrix);
        //        var temp = Output.CloneAs();
        //        Output.ClearEntries();
        //        solver.Solve(Output, temp);
        //    }

        //}

        /// <summary>
        /// Gathers RHS and solution vector on process 0 for more than one MPI process.
        /// </summary>
        /// <param name="__b">local part of rhs</param>
        /// <param name="__x">local part of solution/initial guess</param>
        /// <param name="_x">gathered rhs</param>
        /// <param name="_b">gathered solution vectors</param>
        private void GatherOnProc0(double[] __x, double[] __b, out double[] _x, out double[] _b) {
            int size, rank;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);

            if (size > 1) {
                // gather rhs on processor 0
                // +++++++++++++++++++++++++

                if (rank == 0) {
                    _x = new double[OpMtxRaw.RowPartitioning.TotalLength];
                    _b = new double[OpMtxRaw.RowPartitioning.TotalLength];

                    Array.Copy(__b, 0, _b, 0, OpMtxRaw.RowPartitioning.LocalLength);

                    unsafe {
                        fixed (double* pb = &_b[0]) {
                            for (int rcv_rank = 1; rcv_rank < size; rcv_rank++) {
                                MPI_Status stat;
                                csMPI.Raw.Recv((IntPtr)(pb + OpMtxRaw.RowPartitioning.GetI0Offest(rcv_rank)), OpMtxRaw.RowPartitioning.GetLocalLength(rcv_rank), csMPI.Raw._DATATYPE.DOUBLE, rcv_rank, 342346 + rcv_rank, csMPI.Raw._COMM.WORLD, out stat);
                            }
                        }
                    }

                } else {
                    // send my part to P0
                    unsafe {
                        fixed (double* pb = &__b[0]) {
                            csMPI.Raw.Send((IntPtr)pb, OpMtxRaw.RowPartitioning.LocalLength, csMPI.Raw._DATATYPE.DOUBLE, 0, 342346 + rank, csMPI.Raw._COMM.WORLD);
                        }
                    }

                    _x = null;
                    _b = null;
                }
            } else {
                _x = __x;
                _b = __b;
            }

        }

        /// <summary>
        /// Scatters solution vector from process 0 to other MPI processors.
        /// </summary>
        /// <param name="__x">
        /// input; long vector on proc 0
        /// </param>
        /// <param name="_x">
        /// output;
        /// </param>
        private void ScatterFromProc0(double[] __x, double[] _x) {
            int size, rank;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);

            if (size > 1) {
                // distribute solution to other processors
                // +++++++++++++++++++++++++++++++++++++++

                if (rank == 0) {
                    Array.Copy(_x, 0, __x, 0, OpMtxRaw.RowPartitioning.LocalLength);

                    unsafe {
                        fixed (double* px = &_x[0]) {
                            for (int targ_rank = 1; targ_rank < size; targ_rank++) {
                                csMPI.Raw.Send((IntPtr)(px + OpMtxRaw.RowPartitioning.GetI0Offest(targ_rank)), OpMtxRaw.RowPartitioning.GetLocalLength(targ_rank), csMPI.Raw._DATATYPE.DOUBLE, targ_rank, 4444 + targ_rank, csMPI.Raw._COMM.WORLD);
                            }
                        }
                    }

                } else {
                    unsafe {
                        fixed (double* px = &__x[0]) {
                            MPI_Status _st;
                            csMPI.Raw.Recv((IntPtr)px, OpMtxRaw.RowPartitioning.LocalLength, csMPI.Raw._DATATYPE.DOUBLE, 0, 4444 + rank, csMPI.Raw._COMM.WORLD, out _st);
                        }
                    }
                }
            } else {
                if (!object.ReferenceEquals(__x, _x)) {
                    int L = __x.Length;
                    if (__x.Length != _x.Length)
                        throw new ApplicationException("internal error.");
                    Array.Copy(_x, __x, L);
                }

            }
        }

    }


}

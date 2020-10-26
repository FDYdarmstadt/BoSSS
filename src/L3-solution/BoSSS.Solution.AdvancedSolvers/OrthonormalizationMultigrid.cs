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
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using MPI.Wrappers;
using System.Numerics;
using System.Diagnostics;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Statistic;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// 
    /// </summary>
    public class OrthonormalizationMultigrid : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination {

        /// <summary>
        /// ctor
        /// </summary>
        public OrthonormalizationMultigrid() {
            TerminationCriterion = DefaultTermination;
        }

        private bool DefaultTermination(int iter, double R0_l2, double R_l2) {
            if (iter > 100)
                return false;

            if (R_l2 < R0_l2 * 10e-8 + 10e-8)
                return false;

            return true;
        }

        /// <summary>
        /// ~
        /// </summary>
        public Func<int, double, double, bool> TerminationCriterion {
            get;
            set;
        }


        /// <summary>
        /// The matrix at this level.
        /// </summary>
        public BlockMsrMatrix OpMatrix;

        /// <summary>
        /// passed in <see cref="Init"/>
        /// </summary>
        MultigridOperator m_MgOperator;


        bool CoarseOnLovwerLevel = true;

        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(MultigridOperator op) {
            using (var tr = new FuncTrace()) {
                this.m_MgOperator = op;
                var Mtx = op.OperatorMatrix;
                var MgMap = op.Mapping;
                if(op.LevelIndex == 0)
                    viz = new MGViz(op);

                if (!Mtx.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!Mtx.ColPartition.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Column partitioning mismatch.");

                MxxHistory.Clear();
                SolHistory.Clear();

                
                // set operator
                // ============
                this.OpMatrix = Mtx;


                // initiate coarser level
                // ======================
                if (this.CoarserLevelSolver == null) {
                    //throw new NotSupportedException("Missing coarse level solver.");
                    Console.WriteLine("OrthonormalizationMultigrid: running without coarse solver.");
                } else {
                    if (op.CoarserLevel != null) {
                        this.CoarserLevelSolver.Init(op.CoarserLevel);
                        CoarseOnLovwerLevel = true;
                    } else {
                        Console.WriteLine("OrthonormalizationMultigrid: running coarse solver on same level.");
                        this.CoarserLevelSolver.Init(op);
                        CoarseOnLovwerLevel = false;
                    }
                }

                // init smoother
                // =============
                if (PreSmoother != null)
                    PreSmoother.Init(op);
                if (PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                    PostSmoother.Init(op);
            }
        }


#pragma warning disable 0649
        MGViz viz;
#pragma warning restore 0649

        /// <summary>
        /// coarse-level correction; can be defined either
        /// - on this level (then the coarse solver may perform its of prolongation/restriction), or
        /// - on coarser level, then prolongation/restriction is handled by this solver.
        /// </summary>
        public ISolverSmootherTemplate CoarserLevelSolver;
        
        /// <summary>
        /// high frequency solver before coarse grid correction
        /// </summary>
        public ISolverSmootherTemplate PreSmoother;

        /// <summary>
        /// to be removed.
        /// </summary>
        public ISolverSmootherTemplate DebugSmoother;
        
        /// <summary>
        /// high frequency solver before coarse grid correction
        /// </summary>
        public ISolverSmootherTemplate PostSmoother;
        
        
        
        public int m_omega = 1;
        
        public int MaxKrylovDim = int.MaxValue;

        //public bool SpectralAnalysis;
        /// <summary>
        /// 
        /// </summary>
        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        /// <summary>
        /// computes the residual on this level.
        /// </summary>
        /// <param name="B">input: RHS of the system</param>
        /// <param name="X">input: solution approximation</param>
        /// <param name="Res">output: on exit <paramref name="B"/> - <see cref="OpMatrix"/>*<paramref name="X"/></param>
        public void Residual(double[] Res, double[] X, double[] B) {
            Debug.Assert(Res.Length == m_MgOperator.Mapping.LocalLength);
            Debug.Assert(X.Length == m_MgOperator.Mapping.LocalLength);
            Debug.Assert(B.Length == m_MgOperator.Mapping.LocalLength);
            int L = Res.Length;
            Array.Copy(B, Res, L);
            OpMatrix.SpMV(-1.0, X, 1.0, Res);
            //Res.AccV(1.0, B);


            /*
            int L = Res.Count;

            var OpMatrix_coarse = m_MgOperator.CoarserLevel.OperatorMatrix;
            int Lc = OpMatrix_coarse._RowPartitioning.LocalLength;

            double[] Xc = new double[Lc];
            double[] Bc = new double[Lc];
            m_MgOperator.CoarserLevel.Restrict(X, Xc);
            m_MgOperator.CoarserLevel.Restrict(B, Bc);

            double[] Res_coarse = Bc; Bc = null;
            OpMatrix_coarse.SpMV(-1.0, Xc, 1.0, Res_coarse);

            m_MgOperator.CoarserLevel.Restrict(Res, Res_coarse);

            double[] Res1 = new double[L];// Res.ToArray();
            m_MgOperator.CoarserLevel.Prolongate(-1.0, Res1, 0.0, Res_coarse);

            double alpha = ParallelBlas.MPI_ddot(Res, Res1) / Res1.L2NormPow2().MPISum();



            double[] RemRes = Res.ToArray();
            RemRes.AccV(-alpha, Res1);

            Res1.ScaleV(alpha);


            int iLevel = m_MgOperator.LevelIndex;
            Console.Write("    ");
            for(int iLv = 0; iLv <= iLevel; iLv++) {
                Console.Write("  ");
            }

            Console.WriteLine(iLevel + " " + Res.L2Norm().ToString("0.####E-00") + "   " + RemRes.L2Norm().ToString("0.####E-00") + "   " + Res1.L2Norm().ToString("0.####E-00") + "   " + GenericBlas.InnerProd(RemRes, Res1));

            //viz.PlotVectors(new double[][] { Res.ToArray(), Res1, Res_coarse, RemRes }, new[] { "Residual", "ProloResi", "CarseResi", "RestResi" });
            */
        }


        /// <summary>
        /// solution guesses from smoothers
        /// </summary>
        List<double[]> SolHistory = new List<double[]>();
        
        /// <summary>
        /// - orthonormal system of matrix-vector products;
        /// - the i-th entry is equal to  <see cref="OpMatrix"/>*<see cref="SolHistory"/>[i]
        /// </summary>
        List<double[]> MxxHistory = new List<double[]>();

        
        void AddSol(ref double[] X) {
            using (new FuncTrace()) {
                AddSolCore(ref X);

                /*
                var map = m_MgOperator.Mapping;
                int NoOfVar = map.NoOfVariables;

                for (int iVar = 0; iVar < NoOfVar; iVar++) {
                    var subi = new SubBlockSelector(m_MgOperator.Mapping);
                    subi.VariableSelector(iVar);
                    var bobi = new BlockMask(subi);
                    var bi = bobi.LocalMask;

                    var X_iVar = new double[X.Length];
                    X_iVar.AccV(1.0, X, bi, bi);

                    AddSolCore(ref X_iVar);
                }
                */

                // split solution in high and low modes
                /*
                var Xlo = X.CloneAs();
                var Xhi = X.CloneAs();

                int J = m_MgOperator.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int Np = m_MgOperator.Mapping.AggBasis[0].MaximalLength;
                if (J * Np != X.Length)
                    throw new NotSupportedException("experimental stuff failed");

                int Ncut = 4;
                for(int j = 0; j < J; j++) {
                    for(int n = 0; n < Np; n++) {
                        int l = m_MgOperator.Mapping.LocalUniqueIndex(0, j, n);
                        if (n < Ncut)
                            Xhi[l] = 0;
                        else
                            Xlo[l] = 0;
                    }
                }

                __AddSol(ref Xlo);
                __AddSol(ref Xhi);
                //*/
            }
        }
        

        void AddSolCore(ref double[] X) {
            Debug.Assert(SolHistory.Count == MxxHistory.Count);
            Debug.Assert(X.Length == OpMatrix._RowPartitioning.LocalLength);
            int L = X.Length;
            int KrylovDim = SolHistory.Count;

            double[] Mxx = new double[L];
            OpMatrix.SpMV(1.0, X, 0.0, Mxx);

            //double NormInitial = Mxx.MPI_L2Norm();

            for (int jj = 0; jj < 2; jj++) { // re-orthogonalisation, loop-limit to 2; See book of Saad, p 162, section 6.3.2
                for (int i = 0; i < KrylovDim; i++) {
                    Debug.Assert(!object.ReferenceEquals(Mxx, MxxHistory[i]));
                    double beta = BLAS.ddot(L, Mxx, 1, MxxHistory[i], 1).MPISum();
                    BLAS.daxpy(L, -beta, SolHistory[i], 1, X, 1);
                    BLAS.daxpy(L, -beta, MxxHistory[i], 1, Mxx, 1);
                }

                //double NormAfter = Mxx.MPI_L2Norm();
                //Console.WriteLine("   orthonormalization norm reduction: " + (NormAfter/NormInitial));

                double gamma = 1.0 / Mxx.MPI_L2Norm();
                //double gamma = 1.0 / BLAS.dnrm2(L, Mxx, 1).Pow2().MPISum().Sqrt();
                BLAS.dscal(L, gamma, Mxx, 1);
                BLAS.dscal(L, gamma, X, 1);
            }


            SolHistory.Add(X);
            MxxHistory.Add(Mxx);
            X = null;
        }

        double __MinimizeResidual(double[] outX, double[] Sol0, double[] Res0, double[] outRes, bool diagnosis = false) {
            using (new FuncTrace()) {
                Debug.Assert(SolHistory.Count == MxxHistory.Count);
                Debug.Assert(outX.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(Sol0.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(Res0.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(outRes.Length == m_MgOperator.Mapping.LocalLength);
                
                int KrylovDim = SolHistory.Count;
                int L = outX.Length;

                /*
                double[] alpha = new double[KrylovDim];
                for (int i = 0; i < KrylovDim; i++) {
                    alpha[i] = BLAS.ddot(L, MxxHistory[i], 1, Res0, 1);
                }
                alpha = alpha.MPISum();
                for(int i = 0; i < KrylovDim; i++) {
                    if(alpha[i].IsNaNorInf())
                        throw new ArithmeticException("Numerical breakdown in residual minimization.");
                }


                Array.Copy(Sol0, outX, L);
                Array.Copy(Res0, outRes, L);
                for (int i = 0; i < KrylovDim; i++) {
                    //outX.AccV(alpha[i], SolHistory[i]);
                    //outRes.AccV(-alpha[i], MxxHistory[i]);
                    BLAS.daxpy(L, alpha[i], SolHistory[i], 1, outX, 1);
                    BLAS.daxpy(L, -alpha[i], MxxHistory[i], 1, outRes, 1);
                }*/

                int i = KrylovDim - 1;
                double alpha_i = BLAS.ddot(L, MxxHistory[i], 1, Res0, 1).MPISum();
                BLAS.daxpy(L, alpha_i, SolHistory[i], 1, outX, 1);
                BLAS.daxpy(L, -alpha_i, MxxHistory[i], 1, outRes, 1);




                double ResNorm =  BLAS.dnrm2(L, outRes, 1).Pow2().MPISum().Sqrt();

               

                return ResNorm;

                /* we cannot do the following 
                // since the 'MxxHistory' vectors form an orthonormal system,
                // the L2-norm is the L2-Norm of the 'alpha'-coordinates (Parceval's equality)
                return alpha.L2Norm();            
                */
            }
        }


        double MinimizeResidual(double[] outX, double[] Sol0, double[] Res0, double[] outRes, bool diagnosis = false) {
            using (new FuncTrace()) {
                Debug.Assert(SolHistory.Count == MxxHistory.Count);
                Debug.Assert(outX.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(Sol0.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(Res0.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(outRes.Length == m_MgOperator.Mapping.LocalLength);
                
                int KrylovDim = SolHistory.Count;
                int L = outX.Length;

                
                double[] alpha = new double[KrylovDim];
                for (int i = 0; i < KrylovDim; i++) {
                    alpha[i] = BLAS.ddot(L, MxxHistory[i], 1, Res0, 1);
                }
                alpha = alpha.MPISum();
                for(int i = 0; i < KrylovDim; i++) {
                    if(alpha[i].IsNaNorInf())
                        throw new ArithmeticException("Numerical breakdown in residual minimization.");
                }


                Array.Copy(Sol0, outX, L);
                Array.Copy(Res0, outRes, L);
                for (int i = 0; i < KrylovDim; i++) {
                    //outX.AccV(alpha[i], SolHistory[i]);
                    //outRes.AccV(-alpha[i], MxxHistory[i]);
                    BLAS.daxpy(L, alpha[i], SolHistory[i], 1, outX, 1);
                    BLAS.daxpy(L, -alpha[i], MxxHistory[i], 1, outRes, 1);
                }



                double ResNorm = outRes.MPI_L2Norm();


                Sol0.SetV(outX);
                Res0.SetV(outRes);


                if(KrylovDim > MaxKrylovDim) {
                    int leastSignificantVec = alpha.Select(x => x.Abs()).IndexOfMin();
                    MxxHistory.RemoveAt(leastSignificantVec);
                    SolHistory.RemoveAt(leastSignificantVec);
                }

               

                return ResNorm;

                /* we cannot do the following 
                // since the 'MxxHistory' vectors form an orthonormal system,
                // the L2-norm is the L2-Norm of the 'alpha'-coordinates (Parceval's equality)
                return alpha.L2Norm();            
                */
            }
        }




        ISparseSolver PlottiesFullsolver = null;


        /// <summary>
        /// Debug plotting / only used for development
        /// </summary>
        /// <param name="RawCorr">correction applied</param>
        /// <param name="rl">current residual (for <paramref name="_xl"/>)</param>
        /// <param name="_xl">current solution</param>
        /// <param name="_xl_prev">solution before correction</param>
        /// <param name="B">rhs of the system</param>
        void PlottyMcPlot(double[] rl, double[] _xl, double[] _xl_prev, double[] RawCorr, double[] B) {
            return;
            if (viz == null)
                return;

            if (this.m_MgOperator.LevelIndex > 0)
                return;


            if (_xl_prev == null)
                _xl_prev = new double[rl.Length];
            if (RawCorr == null)
                RawCorr = new double[rl.Length];

            double[] rlcc = rl.CloneAs();
            rlcc.Normalize();

            if (PlottiesFullsolver == null) {
                PlottiesFullsolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver() {
                    CacheFactorization = true
                };
                PlottiesFullsolver.DefineMatrix(m_MgOperator.OperatorMatrix);
            }

            double[] optCorr = new double[rlcc.Length];
            PlottiesFullsolver.Solve(optCorr, rl);

            double[] OrthoCorr = _xl.CloneAs();
            PlottiesFullsolver.Solve(OrthoCorr, B);

            double[] deltaCorr = OrthoCorr.CloneAs();
            deltaCorr.AccV(-1.0, optCorr);


            string[] Names = new[] { "Solution",
                                     "LastCorrection",
                                     "ExactSol",
                                     "ExactCorrection",
                                     "DeltaCorrection",
                                     "Residual",
                                     "NormalizedResidual" };
            this.viz.PlotVectors(new double[][] {
                               _xl.ToArray(), //            current solution
                                RawCorr, //                 last precond result (correction)
                                OrthoCorr, //         
                                optCorr, //                 optimal correction
                                deltaCorr, //               delta of precond and optimal correction
                                rl.ToArray(), //            current residual
                                rlcc //                     normed residual
                            }, Names);

        }

        /// <summary>
        /// the multigrid iterations for a linear problem
        /// </summary>
        /// <param name="_xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="_B">the right-hand-side of the problem</param>
        public void Solve<U, V>(U _xl, V _B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (new FuncTrace()) {
                double[] B, X;
                if (_B is double[])
                    B = _B as double[];
                else
                    B = _B.ToArray();
                if (_xl is double[])
                    X = _xl as double[];
                else
                    X = _xl.ToArray();

                // clear history of coarse solvers
                MxxHistory.Clear();
                SolHistory.Clear();


                //// in case of spectral analysis
                //if (this.m_MgOperator.LevelIndex == 0 && SpectralAnalysis)
                //{
                //    // Set RHS to zero and introduce random intitial guess respectively error
                //    Console.WriteLine("Performing Spectral Analysis, inserting initial error ...");
                //    B.Clear();
                //    X.Clear();
                //    var rand = new Random();
                //    X = Enumerable.Repeat(0, X.Length).Select(i => rand.NextDouble() * 2 - 1).ToArray();                    
                //}

                int L = X.Length;
                int Lc;
                if (this.CoarserLevelSolver != null && CoarseOnLovwerLevel)
                    Lc = m_MgOperator.CoarserLevel.Mapping.LocalLength;
                else
                    Lc = -1;

                double[] rl = new double[L];
                double[] rlc = Lc > 0 ? new double[Lc] : null;


                double[] Sol0 = X.CloneAs();
                double[] Res0 = new double[L];
                Residual(Res0, Sol0, B);
                Array.Copy(Res0, rl, L);


                /*
                ISparseSolver Fullsolver = null;
                string[] Names = new[] { "Solution", "LastCorrection", "ExactCorrection", "DeltaCorrection", "Residual", "NormalizedResidual" };
                if (this.m_MgOperator.LevelIndex == 0 && viz != null) {
                    double[] rlcc = rl.CloneAs();
                    rlcc.Normalize();

                    Fullsolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver() {
                        CacheFactorization = true
                    };
                    Fullsolver.DefineMatrix(m_MgOperator.OperatorMatrix);

                    double[] deltaCorr = new double[rlcc.Length];
                    double[] optCorr = new double[rlcc.Length];
                    double[] dummy = new double[rlcc.Length];

                    this.viz.PlotVectors(new double[][] {
                               _xl.ToArray(), //     current solution
                                dummy, //            last precond result (correction)
                                optCorr, //          optimal correction
                                deltaCorr, //        delta of precond and optimal correction
                                rl.ToArray(), //     current residual
                                rlcc //              normed residual
                            }, Names);
                }
                */

                //PlottyMcPlot(rl, X, null, null, B);
                double[] Xprev = null, Corr = null;
                if (PlottiesFullsolver != null) {
                    Xprev = X.CloneAs();
                    Corr = new double[Xprev.Length];
                }

                double iter0_resNorm = Res0.MPI_L2Norm();
                double resNorm = iter0_resNorm;
                this.IterationCallback?.Invoke(0, Sol0, Res0, this.m_MgOperator);

                //var tmpX = new double[L];

                for (int iIter = 1; true; iIter++) {
                    if (!TerminationCriterion(iIter, iter0_resNorm, resNorm)) {
                        Converged = true;
                        break;
                    }

                    // pre-smoother
                    // ------------

                    {

                        // Test: residual is already computed by 'MinimizeResidual' form previous iter; 
#if DEBUG
                        {
                            double[] rTest = new double[rl.Length];
                            Residual(rTest, X, B); // Residual on this level; 
                            double resDist = rTest.MPI_L2Dist(rl);
                            double resNormTst = rl.MPI_L2Norm();
                            if(resDist > resNormTst * 10e-5)
                                Console.WriteLine($"Residual vector (before pre-smoother) is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
                            //Debug.Assert(resDist <= resNormTst * 10e-5, $"Residual vector is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
                        }
#endif

                        // compute correction
                        double[] PreCorr = new double[L];
                        var oldRl = rl.CloneAs();
                        PreSmoother.Solve(PreCorr, rl); // Vorglättung
                        if (Corr != null) // only for plotting/debugging
                            Corr.SetV(PreCorr);

                        {
                            var checkPreCor = new double[L];
                            DebugSmoother.Solve(PreCorr, oldRl);

                            var diffCorr = PreCorr.CloneAs(); diffCorr.AccV(-1.0, checkPreCor);
                            this.viz.PlotVectors(new double[][] {
                                PreCorr, checkPreCor, diffCorr
                            }, new[] { "PMGon", "noPMG", "diff" });
                        }



                        //tmpX.SetV(X);
                        //tmpX.AccV(1.0, PreCorr);
                        //SpecAnalysisSample(iIter, PreCorr, "smooth1");

                        // orthonormalization and residual minimization
                        AddSol(ref PreCorr);
                        if (Xprev != null)
                            Xprev.SetV(X);
                        resNorm = MinimizeResidual(X, Sol0, Res0, rl);

                        SpecAnalysisSample(iIter, X, "ortho1");

                        if (!TerminationCriterion(iIter, iter0_resNorm, resNorm)) {
                            Converged = true;
                            break;
                        }
                    }

                    PlottyMcPlot(rl, X, Xprev, Corr, B);


                    // coarse grid correction
                    // ----------------------
                    // Test: Residual on this level / already computed by 'MinimizeResidual' above
#if DEBUG
                    {
                        double[] rTest = new double[rl.Length];
                        Residual(rTest, X, B); // Residual on this level; 
                                               // Test also fails if convergence criterium is to strict because then machine accuracy is reached
                        double resDist = rTest.MPI_L2Dist(rl);
                        double resNormTst = rl.MPI_L2Norm();
                        if(resDist > resNormTst * 10e-5)
                            Console.WriteLine($"Residual vector (after pre-smoother/before coarse-correction) is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
                        //Debug.Assert(resDist <= resNormTst * 10e-5, $"Residual vector is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
                    }
#endif

                    for (int i = 0; i < m_omega; i++) {
                        if (this.CoarserLevelSolver != null) {

                            double[] vl = new double[L];
                            if (CoarseOnLovwerLevel) {
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // coarse grid solver defined on COARSER MESH LEVEL:
                                // this solver must perform restriction and prolongation
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                // restriction of residual
                                this.m_MgOperator.CoarserLevel.Restrict(rl, rlc);

                                // Berechnung der Grobgitterkorrektur
                                double[] vlc = new double[Lc];
                                this.CoarserLevelSolver.Solve(vlc, rlc);

                                // Prolongation der Grobgitterkorrektur
                                this.m_MgOperator.CoarserLevel.Prolongate(1.0, vl, 1.0, vlc);


                            } else {
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // coarse grid solver defined on the SAME MESH LEVEL:
                                // performs (probably) its own restriction/prolongation
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                // Berechnung der Grobgitterkorrektur
                                this.CoarserLevelSolver.Solve(vl, rl);
                            }


                            // record correction for Debug-Plotting 
                            if (Corr != null)
                                Corr.SetV(vl);

                            //tmpX.SetV(X);
                            //tmpX.AccV(1.0, vl);
                            //SpecAnalysisSample(iIter, tmpX, "cgc");

                            // orthonormalization and residual minimization
                            AddSol(ref vl);
                            if (Xprev != null)
                                Xprev.SetV(X);
                            resNorm = MinimizeResidual(X, Sol0, Res0, rl);

                            SpecAnalysisSample(iIter, X, "ortho2");

                            if (!TerminationCriterion(iIter, iter0_resNorm, resNorm)) {
                                Converged = true;
                                break;
                            }

                        }

                        PlottyMcPlot(rl, X, Xprev, Corr, B);

                        // post-smoother
                        // -------------

                        for (int g = 0; g < 2; g++) { // doppelt hält besser
                                                      // Test: Residual on this level / already computed by 'MinimizeResidual' above
#if DEBUG
                            {
                                double[] rTest = new double[rl.Length];
                                Residual(rTest, X, B); // Residual on this level; 
                                double resDist = rTest.MPI_L2Dist(rl);
                                double resNormTst = rl.MPI_L2Norm();
                                if(resDist > resNormTst * 10e-5)
                                    Console.WriteLine($"Residual vector (before post-smoother run #{g}) is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
                                //Debug.Assert(resDist <= resNormTst * 10e-5, $"Residual vector is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
                            }
#endif

                            // compute correction
                            double[] PostCorr = new double[L];
                            PostSmoother.Solve(PostCorr, rl); // Vorglättung
                            if (Corr != null)
                                Corr.SetV(PostCorr);

                            //tmpX.SetV(X);
                            //tmpX.AccV(1.0, PreCorr);
                            //SpecAnalysisSample(iIter, tmpX, "smooth2_" + g);

                            // orthonormalization and residual minimization
                            AddSol(ref PostCorr);
                            if (Xprev != null)
                                Xprev.SetV(X);
                            resNorm = MinimizeResidual(X, Sol0, Res0, rl);

                            SpecAnalysisSample(iIter, X, "ortho3_" + g);

                            if (!TerminationCriterion(iIter, iter0_resNorm, resNorm)) {
                                Converged = true;
                                break;
                            }
                        }

                        // iteration callback
                        // ------------------

                        this.ThisLevelIterations++;

                        IterationCallback?.Invoke(iIter, X, rl, this.m_MgOperator);

                    }


                    // solution copy
                    // =============
                    if (!ReferenceEquals(_xl, X)) {
                        _xl.SetV(X);
                    }
                }
            }
        }

        /// <summary>
        /// the multigrid iterations for a linear problem (experimental version)
        /// </summary>
        /// <param name="_xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="_B">the right-hand-side of the problem</param>
        public void Solve__experimento<U, V>(U _xl, V _B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (new FuncTrace()) {
                double[] B, X;
                if (_B is double[])
                    B = _B as double[];
                else
                    B = _B.ToArray();
                if (_xl is double[])
                    X = _xl as double[];
                else
                    X = _xl.ToArray();

                //// clear history, makes a small difference on coarse levels, which one is better?
                //MxxHistory.Clear();
                //SolHistory.Clear();


                //// in case of spectral analysis
                //if (this.m_MgOperator.LevelIndex == 0 && SpectralAnalysis) {
                //    // Set RHS to zero and introduce random intitial guess respectively error
                //    Console.WriteLine("Performing Spectral Analysis, inserting initial error ...");
                //    B.Clear();
                //    X.Clear();
                //    var rand = new Random();
                //    X = Enumerable.Repeat(0, X.Length).Select(i => rand.NextDouble() * 2 - 1).ToArray();
                //}

                int L = X.Length;
                int Lc;
                if (this.CoarserLevelSolver != null && CoarseOnLovwerLevel)
                    Lc = m_MgOperator.CoarserLevel.Mapping.LocalLength;
                else
                    Lc = -1;

                //double[] rl = new double[L];
                double[] rlc = Lc > 0 ? new double[Lc] : null;

                // Residual of initial solution guess
                double[] Res = new double[L];
                Residual(Res, X, B);


                // solution loop
                double iter0_resNorm = Res.MPI_L2Norm();
                double resNorm = iter0_resNorm;
                this.IterationCallback?.Invoke(0, X, Res, this.m_MgOperator);
                for (int iIter = 1; true; iIter++) {
                    if (!TerminationCriterion(iIter, iter0_resNorm, resNorm)) {
                        Converged = true;
                        break;
                    }


                    // coarse grid correction
                    // ----------------------

                    if (this.CoarserLevelSolver != null) {
                        double[] vl = new double[L];
                        if (CoarseOnLovwerLevel) {
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // coarse grid solver defined on COARSER MESH LEVEL:
                            // this solver must perform restriction and prolongation
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            // restriction of residual
                            this.m_MgOperator.CoarserLevel.Restrict(Res, rlc);

                            // Berechnung der Grobgitterkorrektur
                            double[] vlc = new double[Lc];
                            for (int i = 0; i < m_omega; i++)
                                this.CoarserLevelSolver.Solve(vlc, rlc);

                            // Prolongation der Grobgitterkorrektur
                            this.m_MgOperator.CoarserLevel.Prolongate(1.0, vl, 1.0, vlc);


                        } else {
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            // coarse grid solver defined on the SAME MESH LEVEL:
                            // performs (probably) its own restriction/prolongation
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

                            // Berechnung der Grobgitterkorrektur
                            this.CoarserLevelSolver.Solve(vl, Res);
                        }

                        // correct solution & update residual
                        var Xprev = X.CloneAs();
                        BLAS.daxpy(L, 1.0, vl, 1, X, 1); // add correction: X = X + vl
                        Residual(Res, X, B);

                        IterationCallback?.Invoke(iIter * 2 - 1, X, Res, this.m_MgOperator);

                        PlottyMcPlot(Res, X, Xprev, vl, B);


                        // check termination:
                        if (!TerminationCriterion(iIter, iter0_resNorm, resNorm)) {
                            Converged = true;
                            break;
                        }
                    }

                    // smoother
                    // ---------

                    for (int g = 0; g < 2; g++) {

                        // compute correction
                        double[] PreCorr = new double[L];
                        PreSmoother.Solve(PreCorr, Res); // 

                        // orthonormalization and residual minimization
                        double[] RawCorr = PreCorr.CloneAs();
                        AddSol(ref PreCorr);
                        double[] newX = new double[L];
                        double[] newRes = new double[L];
                        resNorm = MinimizeResidual(newX, X, Res, newRes);

                        PlottyMcPlot(newRes, newX, X, RawCorr, B);

                        X.SetV(newX);
                        Res.SetV(newRes);

                        if (!TerminationCriterion(iIter, iter0_resNorm, resNorm)) {
                            Converged = true;
                            break;
                        }


                    }

                    // iteration callback
                    // ------------------

                    this.ThisLevelIterations++;

                    IterationCallback?.Invoke(iIter * 2, X, Res, this.m_MgOperator);

                }


                // solution copy
                // =============
                //IterationCallback?.Invoke(iIter + 1, X, rl, this.m_MgOperator);
                if (!ReferenceEquals(_xl, X)) {
                    _xl.SetV(X);
                }
            }
        }

        private double[] cloneofX;

        /// <summary>
        /// RealX is left unchanged, no worries
        /// </summary>
        /// <param name="iter"></param>
        /// <param name="RealX"></param>
        /// <param name="name"></param>
        private void SpecAnalysisSample(int iter, double[] RealX, string name) {
            if (iter % 5 != 0 && iter!=1)
                return;
            if (cloneofX == null) cloneofX = new double[RealX.Length];
            cloneofX.SetV(RealX);
            if (this.m_MgOperator.LevelIndex == 0)
                ExtractSamples?.Invoke(iter, cloneofX, name);
        }

        /// <summary>
        /// gefrickel, could be integrated into IterationCallback
        /// without individual string of course
        /// </summary>
        public Action<int, double[],string> ExtractSamples {
            get;
            set;
        }




        /// <summary>
        /// Currently not used
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Sol0"></param>
        /// <param name="rl"></param>
        /// <param name="Res0"></param>
        /// <param name="B"></param>
        /// <param name="L"></param>
        private void Restart(double[] X, double[] Sol0, double[] rl, double[] Res0, double[] B, int L) {
            MxxHistory.Clear();
            SolHistory.Clear();
            Array.Copy(X, Sol0, L);
            Residual(rl, X, B);
            Array.Copy(rl, Res0, L);
            Console.WriteLine("      restarted with residual = " + rl.L2Norm() + " on MG level " + m_MgOperator.LevelIndex);
        }



        /// <summary>
        /// ~
        /// </summary>
        public int IterationsInNested {
            get {
                int iter = 0;

                if (PreSmoother != null)
                    iter += this.PreSmoother.IterationsInNested + this.PreSmoother.ThisLevelIterations;

                if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                    iter += this.PostSmoother.IterationsInNested + this.PostSmoother.ThisLevelIterations;

                iter += this.CoarserLevelSolver.IterationsInNested + this.CoarserLevelSolver.ThisLevelIterations;

                return iter;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public int ThisLevelIterations {
            get;
            private set;
        }

        /// <summary>
        /// %
        /// </summary>
        public bool Converged {
            get;
            private set;
        }
        
        /// <summary>
        /// %
        /// </summary>
        public void ResetStat() {
            this.Converged = false;
            this.ThisLevelIterations = 0;
            if (this.PreSmoother != null)
                this.PreSmoother.ResetStat();
            if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                this.PostSmoother.ResetStat();
            if (this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.ResetStat();
        }

        /// <summary>
        /// %
        /// </summary>
        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }
    }
}

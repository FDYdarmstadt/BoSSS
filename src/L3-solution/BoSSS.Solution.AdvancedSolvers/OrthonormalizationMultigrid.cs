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
using System.IO;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// 
    /// </summary>
    public class OrthonormalizationMultigrid : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination {

        ///// <summary>
        /////  Individual configuration of <see cref="OrthonormalizationMultigrid"/>
        ///// </summary>
        //public class myConfig : ConfigBase {
        //    /// <summary>
        //    /// config cctor of <see cref="OrthonormalizationMultigrid"/>
        //    /// </summary>
        //    public myConfig() {
        public int MaxKrylovDim = int.MaxValue;
        //    }
            /// <summary>
            /// No restriction and prolongation on coarsest grid
            /// </summary>
            public bool CoarseOnLovwerLevel = true;
            
            /// <summary>
            /// W-cycle
            /// </summary>
            public int m_omega = 1;

        //    /// <summary>
        //    /// ~
        //    /// </summary>
        //    /// <returns></returns>
        //    public override ISolverSmootherTemplate GetInstance() {
        //        return new OrthonormalizationMultigrid();
        //    }
        //}

        //myConfig m_config;
        
        ///// <summary>
        ///// ~
        ///// </summary>
        //public ConfigBase Config {
        //    get {
        //        if (m_config == null)
        //            m_config = new myConfig();
        //        return m_config;
        //    }
        //}

        /// <summary>
        /// ctor
        /// </summary>
        public OrthonormalizationMultigrid() {
            TerminationCriterion = DefaultTermination;
        }

        static private (bool bNotTerminate, bool bSuccess) DefaultTermination(int iter, double R0_l2, double R_l2) {
            if(iter > 100)
                return (false, false); // fail

            if(R_l2 < R0_l2 * 10e-8 + 10e-8)
                return (false, true); // success

            return (true, false); // keep running
        }


        /// <summary>
        /// ~
        /// </summary>
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
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
                TrackMemory(1);
                if (!Mtx.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!Mtx.ColPartition.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Column partitioning mismatch.");

                MxxHistory.Clear();
                SolHistory.Clear();
                Alphas.Clear();
                

                // set operator
                // ============
                this.OpMatrix = Mtx;

#if TEST
                Console.WriteLine("level: {0} cells: {1} degree: {2}", op.LevelIndex, op.Mapping.LocalNoOfBlocks, op.Degrees[0]);
#endif
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
#if TEST
                try {
                    op.OperatorMatrix.CheckForNanOrInfM();
                } catch (Exception ex) {
                    Console.WriteLine("Arithmetic exception in OperatorMatrix");
                    Console.WriteLine(ex.Message);
                    op.OperatorMatrix.SaveToTextFileSparse("A");
                }
                
#endif
            }
        }


        private void TrackMemory(int pos) {
#if TEST
            if (m_MgOperator.LevelIndex != 0) return;
            long memWork = 0, memPrivate = 0, memGC = 0;
            Process myself = Process.GetCurrentProcess();
            {
                try {
                    memWork = myself.WorkingSet64 / (1024 * 1024);
                    memPrivate = myself.PrivateMemorySize64 / (1024 * 1024);
                    memGC = GC.GetTotalMemory(false) / (1024 * 1024);
                    memWork = memWork.MPISum();
                    memPrivate = memPrivate.MPISum();
                    memGC = memGC.MPISum();
                } catch (Exception e) {
                    Console.WriteLine(e.Message);
                }
            }
            if (m_MgOperator.Mapping.MpiRank == 0) {
                var strw = m_MTracker;
                if (m_MTracker == null) {
                    m_MTracker = new StreamWriter("MemoryTrack", true);
                    strw = m_MTracker;
                    strw.Write("pos\t");
                    strw.Write("workingset\t");
                    strw.Write("private\t");
                    strw.Write("GC\t\n");
                }
                strw.Write(pos + "\t");
                strw.Write(memWork + "\t");
                strw.Write(memPrivate + "\t");
                strw.Write(memGC + "\t\n");
                strw.Flush();
            }
#endif
        }

        StreamWriter m_MTracker = null;

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

        ///// <summary>
        ///// to be removed.
        ///// </summary>
        //public ISolverSmootherTemplate DebugSmoother;
        
        /// <summary>
        /// high frequency solver before coarse grid correction
        /// </summary>
        public ISolverSmootherTemplate PostSmoother;
        
        
        

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
            using (new FuncTrace()) {
                Debug.Assert(Res.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(X.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(B.Length == m_MgOperator.Mapping.LocalLength);
                int L = Res.Length;
                //Res.SetV(B);
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


        /// <summary>
        /// scaling factors which were applied to <see cref="SolHistory"/> to approximate the solution
        /// </summary>
        List<(double,double,int)> Alphas = new List<(double,double,int)>();

        void AddSol(ref double[] X) {
            using(var ft = new FuncTrace()) {

                double FillXwithRandom(double[] __X, double[] __Mxx) {
                    double __NormMxx;
                    ft.Error("Solution norm is exactly 0.0; trying with a random vector instead to recover.");
                    __X.FillRandom();

                    __Mxx.ClearEntries();
                    OpMatrix.SpMV(1.0, __X, 0.0, __Mxx);
                    __NormMxx = __Mxx.MPI_L2Norm();

                    if(__NormMxx == 0)
                        throw new ArithmeticException("Numerical breakdown: norm of matrix * solution after using RANDOM NUMBERS (!) is " + __NormMxx);
                    if(__NormMxx.IsNaNorInf())
                        throw new ArithmeticException("Numerical breakdown: norm of matrix * solution after using RANDOM NUMBERS (!) is " + __NormMxx);
                    return __NormMxx;
                }



                Debug.Assert(SolHistory.Count == MxxHistory.Count);
                Debug.Assert(X.Length == OpMatrix._RowPartitioning.LocalLength);
                int L = X.Length;
                int KrylovDim = SolHistory.Count;


                double[] Mxx = new double[L];
                OpMatrix.SpMV(1.0, X, 0.0, Mxx);
                double NormMxx = Mxx.MPI_L2Norm();
                if(NormMxx.IsNaNorInf())
                    throw new ArithmeticException("Numerical breakdown: norm of matrix * solution is " + NormMxx);

                if(NormMxx == 0) {
                    // solution is really strange: try with a random vector.
                    NormMxx = FillXwithRandom(X, Mxx);
                }

                // scale Mxx to norm 1, should allow better stability when compared to the Krylov vectors who have all norm 1;
                Mxx.ScaleV(1.0 / NormMxx);
                X.ScaleV(1.0 / NormMxx);

                for (int jj = 0; jj < 10; jj++) { // re-orthogonalisation, loop-limit to 2; See book of Saad, p 162, section 6.3.2

                    for (int i = 0; i < KrylovDim; i++) {
                        Debug.Assert(!object.ReferenceEquals(Mxx, MxxHistory[i]));
                        double beta = BLAS.ddot(L, Mxx, 1, MxxHistory[i], 1).MPISum();
                        BLAS.daxpy(L, -beta, SolHistory[i], 1, X, 1);
                        BLAS.daxpy(L, -beta, MxxHistory[i], 1, Mxx, 1);
                    }

                    //double NormAfter = Mxx.MPI_L2Norm();
                    //Console.WriteLine("   orthonormalization norm reduction: " + (NormAfter/NormInitial));

                    double NewMxxNorm = Mxx.MPI_L2Norm();

                    if (NewMxxNorm <= 1E-5) {
                        using (new BlockTrace("re-orthonormalization", ft)) {
                            // a lot of canceling out has occurred.
                            // do another loop, to ensure ortho-normality:
                            //
                            // Mxx and X lost more than 5 Magnitudes, so we might have lost a couple of digits of accuracy

                            // Note: although |Mxx| might be small, |X| can be quite large (and probably random).
                            // To ensure stability, we must start over with a re-scaled X!
                            // We have to re-scale what is remaining of X:
                            double Xnorm = X.MPI_L2Norm();
                            Console.WriteLine("severe cancellation may have occurred, attempting re-orthonormalization");
                            ft.Info("severe cancellation may have occurred, norm after orthogonalization is " + NewMxxNorm + "; norm of X: " + Xnorm + " Doing Re-orthonormalization (" + (jj + 1) + ")");

                            if (Xnorm < 1e-200) {
                                // prohibits div by 0, if we got zero solution  
                                NormMxx = FillXwithRandom(X, Mxx);
                            } else {
                                X.ScaleV(1.0 / Xnorm);
                                Mxx.ClearEntries();
                                OpMatrix.SpMV(1.0, X, 0.0, Mxx);
                                NormMxx = Mxx.MPI_L2Norm();

                                if (NormMxx == 0)
                                    // solution is really strange: try with a random vector.
                                    NormMxx = FillXwithRandom(X, Mxx);
                            }

                            double gamma = 1 / NewMxxNorm;
                            BLAS.dscal(L, gamma, Mxx, 1);
                            BLAS.dscal(L, gamma, X, 1);
                        }
                    } else {
                        double gamma = 1 / NewMxxNorm;
                        BLAS.dscal(L, gamma, Mxx, 1);
                        BLAS.dscal(L, gamma, X, 1);
                        break;
                    }
                }


                SolHistory.Add(X);
                MxxHistory.Add(Mxx);
                X = null;
            }
        }

        /// <summary>
        /// optimized version of <see cref="MinimizeResidual"/>
        /// </summary>
        /// <param name="outX">
        /// - input: the solution build upon all vectors in the Krylov space **except the last one**
        /// - output: on exit, updated by the contribution from the most recent Krylov basis
        /// </param>
        /// <param name="Res0">
        /// residual with respect to the initial guess
        /// </param>
        /// <param name="Sol0">
        /// initial guess; not used
        /// </param>
        /// <param name="outRes">
        /// - input: the residual with respect to the input value of <paramref name="outX"/>
        /// - output: on exit, the residual with respect to the output value of <paramref name="outX"/>
        /// </param>
        /// <param name="id"></param>
        double MinimizeResidual(double[] outX, double[] Sol0, double[] Res0, double[] outRes, int id) {
            using(var ft = new FuncTrace()) {
                Debug.Assert(SolHistory.Count == MxxHistory.Count);
                Debug.Assert(outX.Length == m_MgOperator.Mapping.LocalLength);
                //Debug.Assert(Sol0.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(Res0.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(outRes.Length == m_MgOperator.Mapping.LocalLength);



                int KrylovDim = SolHistory.Count;
                int L = outX.Length;
                
                if(Alphas.Count != KrylovDim -1) {
                    throw new ApplicationException();
                }

                // note: this implementation is based on the assumption, that the 
                // factors alpha_0 ... alpha_(i-1) are equal to the previous call;
                // therefore, it is only necessary to apply the contributions from the most recent Krylov vector

                var oldResiNorm = outRes.MPI_L2Norm();

                int i = KrylovDim - 1;
                double alpha_i = BLAS.ddot(L, MxxHistory[i], 1, Res0, 1).MPISum();
                BLAS.daxpy(L, alpha_i, SolHistory[i], 1, outX, 1); // accumulate solution correction...
                BLAS.daxpy(L, -alpha_i, MxxHistory[i], 1, outRes, 1); // ...and its effect on the residual

                double xNorm = outX.MPI_L2Norm();
                ft.Info("|x|: " + xNorm + ", last alpha was " + alpha_i);


                /*
                // -------------------
                var rjCheck = Res0.CloneAs();
                OpMatrix.SpMV(-1.0, outX, 1.0, rjCheck);
                double dist = rjCheck.L2Dist(outRes);

                
                var OrthoCheck = MultidimensionalArray.Create(i + 1, i + 1);
                for(int k = 0; k <= i; k++) {
                    for(int w = 0; w <= i; w++) {
                        OrthoCheck[k, w] = MxxHistory[k].InnerProd(MxxHistory[w]).MPISum();
                    }
                }
                OrthoCheck.AccEye(-1.0);
                var orthoErr = OrthoCheck.InfNorm();

                double newResiNorm = outRes.L2Norm();
                Debug.Assert(dist <= Math.Max(outRes.L2Norm(), rjCheck.L2Norm()) * 1e-10);
                Debug.Assert(newResiNorm < oldResiNorm * 1.0001, "residual should shrink");
                Debug.Assert(orthoErr < 1.0e-7);
                // ------------------------
                */


                double ResNorm = outRes.MPI_L2Norm();
                Alphas.Add((alpha_i, oldResiNorm/ResNorm, id));
                return ResNorm;
            }
        }

        /*
        /// <summary>
        /// original residual minimization, which re-computes all contributions
        /// </summary>
        double MinimizeResidual(double[] outX, double[] Sol0, double[] Res0, double[] outRes, int id) {
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
                    BLAS.daxpy(L, alpha[i], SolHistory[i], 1, outX, 1);
                    BLAS.daxpy(L, -alpha[i], MxxHistory[i], 1, outRes, 1);
                }

                this.Alphas.Clear();
                for(int i = 0; i < alpha.Length; i++)
                    this.Alphas.Add((alpha[i], double.NaN, -2));

                double ResNorm = outRes.MPI_L2Norm();

                if(KrylovDim > MaxKrylovDim) {
                    int leastSignificantVec = alpha.Select(x => x.Abs()).IndexOfMin();
                    MxxHistory.RemoveAt(leastSignificantVec);
                    SolHistory.RemoveAt(leastSignificantVec);
                }

               

                return ResNorm;
            }
        }

        */

        private void CatchThisExclamationmark(double[] Vector, string name, int position) {
            double L2norm = Vector.L2Norm();
            Console.WriteLine("Norm of {0} at {1}: {2}", name, position, L2norm);
            try {
                Vector.CheckForNanOrInfV();
            } catch (Exception ex) {
                Console.WriteLine("Arithmetic Exception in {0}, at {1}:",name,position);
                Console.WriteLine(ex.Message);
                Vector.SaveToTextFile(name+"_"+position);
            }
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
            //Solve__experimento(_xl, _B);
            //return;


            using(var f = new FuncTrace()) {
                if(this.m_MgOperator.LevelIndex == 0) {
                    GC.Collect();
                }



                double[] B, X;
                if(_B is double[])
                    B = _B as double[];
                else
                    B = _B.ToArray();
                if(_xl is double[])
                    X = _xl as double[];
                else
                    X = _xl.ToArray();

                // clear history of coarse solvers
                MxxHistory.Clear();
                SolHistory.Clear();
                Alphas.Clear();


                int L = X.Length;
                int Lc;
                if(this.CoarserLevelSolver != null && CoarseOnLovwerLevel)
                    Lc = m_MgOperator.CoarserLevel.Mapping.LocalLength;
                else
                    Lc = -1;

                double[] Res = new double[L];
                double[] ResCoarse = Lc > 0 ? new double[Lc] : null;


                double[] Sol0 = X.CloneAs();
                double[] Res0 = new double[L];
                Residual(Res0, Sol0, B);
                Array.Copy(Res0, Res, L);

#if TEST
                int pos=0;
                CatchThisExclamationmark(Res,"Res", pos);
                CatchThisExclamationmark(Sol0, "Sol", pos);
                pos++;
#endif

                double iter0_resNorm = Res0.MPI_L2Norm();
                double resNorm = iter0_resNorm;
                this.IterationCallback?.Invoke(0, Sol0, Res0, this.m_MgOperator);

                
                for(int iIter = 1; true; iIter++) {
                    var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if(!termState.bNotTerminate) {
                        Converged = termState.bSuccess;
                        break;
                    }

                    TrackMemory(2);

                    // pre-smoother
                    // ------------

                    {
                        VerivyCurrentResidual(X, B, Res, iIter);


                        // compute correction
                        //var oldRl = rl.CloneAs();

                        if(PreSmoother != null) {
                            double[] PreCorr = new double[L];
                            PreSmoother.Solve(PreCorr, Res); // Vorglättung

                            // orthonormalization and residual minimization
                            AddSol(ref PreCorr);
                            resNorm = MinimizeResidual(X, Sol0, Res0, Res, 1);

                        }
#if TEST
                        CatchThisExclamationmark(Res, "Res", pos);
                        CatchThisExclamationmark(X, "Sol", pos);
                        pos++;
#endif

                        //SpecAnalysisSample(iIter, X, "ortho1");
                        var termState2 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                        if(!termState2.bNotTerminate) {
                            Converged = termState2.bSuccess;
                            break;
                        }
                        
                    }

                    //PlottyMcPlot(rl, X, Xprev, Corr, B);


                    // coarse grid correction
                    // ----------------------
                    // Test: Residual on this level / already computed by 'MinimizeResidual' above
                    VerivyCurrentResidual(X, B, Res, iIter);


                    for(int i = 0; i < m_omega; i++) {
                        if(this.CoarserLevelSolver != null) {

                            double[] vl = new double[L];
                            if(CoarseOnLovwerLevel) {
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // coarse grid solver defined on COARSER MESH LEVEL:
                                // this solver must perform restriction and prolongation
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                using (new BlockTrace("Restriction",f)) {
                                    // restriction of residual
                                    this.m_MgOperator.CoarserLevel.Restrict(Res, ResCoarse);
                                }
                                // Berechnung der Grobgitterkorrektur
                                double[] vlc = new double[Lc];
                                this.CoarserLevelSolver.Solve(vlc, ResCoarse);
                                using (new BlockTrace("Prolongation", f)) {
                                    // Prolongation der Grobgitterkorrektur
                                    this.m_MgOperator.CoarserLevel.Prolongate(1.0, vl, 1.0, vlc);
                                }

                            } else {
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // coarse grid solver defined on the SAME MESH LEVEL:
                                // performs (probably) its own restriction/prolongation
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                // Berechnung der Grobgitterkorrektur
                                this.CoarserLevelSolver.Solve(vl, Res);
                            }



                            // orthonormalization and residual minimization
                            AddSol(ref vl);
                            resNorm = MinimizeResidual(X, Sol0, Res0, Res, 2);


#if TEST
                            CatchThisExclamationmark(Res, "Res", pos);
                            CatchThisExclamationmark(X, "Sol", pos);
                            pos++;
#endif
                            //SpecAnalysisSample(iIter, X, "ortho2");



                        }
                    } // end of coarse-solver loop

                    var termState3 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if(!termState3.bNotTerminate) {
                        Converged = termState3.bSuccess;
                        break;
                    }



                    // post-smoother
                    // -------------

                    for(int g = 0; g < 2; g++) { // doppelt hält besser
                                                 // Test: Residual on this level / already computed by 'MinimizeResidual' above
                        
                        VerivyCurrentResidual(X, B, Res, iIter);

                        
                        if (PostSmoother != null) {
                            double[] PostCorr = new double[L];
                            PostSmoother.Solve(PostCorr, Res); // compute correction (Nachglättung)
                            AddSol(ref PostCorr); //orthonormalization and residual minimization
                            resNorm = MinimizeResidual(X, Sol0, Res0, Res, 3 + g);
                        }

#if TEST
                        CatchThisExclamationmark(Res, "Res", pos);
                        CatchThisExclamationmark(X, "Sol", pos);
                        pos++;
#endif

                        //SpecAnalysisSample(iIter, X, "ortho3_" + g);
                    } // end of post-smoother loop

                    // iteration callback
                    // ------------------
                    this.ThisLevelIterations++;
                    IterationCallback?.Invoke(iIter, X, Res, this.m_MgOperator);

                    SpecAnalysisSample(iIter, X, "_");

                    var termState4 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if(!termState4.bNotTerminate) {
                        Converged = termState4.bSuccess;
                        break;
                    }

                } // end of solver iterations

                // solution copy
                // =============
                if(!ReferenceEquals(_xl, X)) {
                    _xl.SetV(X);
                }

            } // end of functrace
        }

        private void VerivyCurrentResidual(double[] X, double[] B, double[] Res, int iter) {
            using (new FuncTrace()) {
#if DEBUG
            {
#else
                if (iter % 20 == 0 && iter > 1) {
#endif
                    double[] rTest = new double[Res.Length];
                    Residual(rTest, X, B); // Residual on this level; 
                                           // Test also fails if convergence criterium is to strict because then machine accuracy is reached
                    double resDist = rTest.MPI_L2Dist(Res);
                    //Console.WriteLine("verified Residual: " + resDist);
                    double resNormTst = Res.MPI_L2Norm();
                    double XnormTest = X.MPI_L2Norm();
                    if (resDist > resNormTst * 10e-5 + XnormTest * 1e-5)
                        throw new ArithmeticException($"Residual vector (after pre-smoother/before coarse-correction) is not up-to-date: distance is {resDist}, reference value {resNormTst}");
                    //Debug.Assert(resDist <= resNormTst * 10e-5, $"Residual vector is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
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
            
            if(ExtractSamples == null || this.m_MgOperator.LevelIndex > 0)
                return; // zero overhead


            if (cloneofX == null) 
                cloneofX = new double[RealX.Length];
            cloneofX.SetV(RealX);
            ExtractSamples.Invoke(iter, cloneofX, name);
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

        bool m_verbose = true;

        public long UsedMemory() {
            long Memory = 0;
            Memory += MemoryOfMultigrid();
            Memory += MemoryOfSmoother();
            return Memory;
        }

        public long MemoryOfSmoother() {
            long Memory = 0;
            if (this.CoarserLevelSolver is OrthonormalizationMultigrid)
                Memory += (this.CoarserLevelSolver as OrthonormalizationMultigrid).MemoryOfSmoother();
            if(PreSmoother != null) Memory += PreSmoother.UsedMemory();
            if (PostSmoother != null) Memory += PostSmoother.UsedMemory();
            return Memory;
        }

        public long MemoryOfMultigrid() {
            long Memory = 0;
            if (this.CoarserLevelSolver is OrthonormalizationMultigrid)
                Memory += (this.CoarserLevelSolver as OrthonormalizationMultigrid).MemoryOfMultigrid();
            int SizeSol = this.SolHistory.Count() * this.SolHistory[0].Length * sizeof(double);
            int SizeMxx = this.MxxHistory.Count() * this.MxxHistory[0].Length * sizeof(double);
            int SizeAlpha = this.Alphas.Count() * (sizeof(double)*2+sizeof(int));
            Memory += (SizeSol + SizeMxx + SizeAlpha);
            return Memory;
        }

        public void Dispose() {
            TrackMemory(3);
            //if (m_MTracker != null) m_MTracker.Dispose();
            if (m_verbose && this.m_MgOperator.LevelIndex == 0) {
                Console.WriteLine($"OrthoMG - total memory: {UsedMemory()} MB");
                Console.WriteLine($"OrthoMG - internal memory: {MemoryOfMultigrid()} MB");
                Console.WriteLine($"OrthoMG - smoother memory: {MemoryOfSmoother()} MB");
            }
            if(this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.Dispose();
            this.SolHistory.Clear();
            this.MxxHistory.Clear();
            this.Alphas.Clear();
            this.SolHistory = null;
            this.MxxHistory = null;
            this.Alphas = null;

            if (PreSmoother != null) this.PreSmoother.Dispose();
            if (PostSmoother != null) this.PostSmoother.Dispose();
            this.PreSmoother = null;
            this.PostSmoother = null;
            this.CoarserLevelSolver = null;
        }
    }
}

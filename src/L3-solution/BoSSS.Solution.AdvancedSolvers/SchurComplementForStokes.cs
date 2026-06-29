using BoSSS.Foundation.Quadrature;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.monkey;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework.Constraints;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// A global Schur complement solver for Stokes problems.
    /// </summary>
    public class SchurComplementForStokes : ISolverSmootherTemplate, IProgrammableTermination {


        public SchurComplementForStokes() { 
            this.TerminationCriterion = m_config.DefaultTermination;
        }


        [Serializable]
        public class Config : IterativeSolverConfig {

            public Config() {
                base.WriteDebugInfoLevel = 3;
            }


            override public string Name => "Schur Complement";

            override public string Shortname => "SchurComp";

            [DataMember]
            public ISolverFactory ConvectionDiffusionSolver = new AdvancedSolvers.OrthoMGSchwarzConfig() {
#if DEBUG
                CoarseKickIn = 100, // for debugging purposes, enforce a later switch to direct solver, so that we have a non-trivial multigrid depth on coarse meshes.
                NoOfMultigridLevels = 2
#endif
            };

            override public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                var instance = new SchurComplementForStokes();
                instance.Init(level);
                return instance;
            }

            [DataMember]
            public int NumberOfIterationsPerSmoother = 4;


            /// <summary>
            /// - if set to 1, a this performs a V-cycle
            /// - if set to 2 or higher, a W-cycle, resp. WW-cycle, WWW-cycle, etc.
            /// </summary>
            [DataMember]
            public int NoOfCoarseSolverSweeps = 1;

            /// <summary>
            /// 
            /// </summary>
            public int NoOfPostSmootherSweeps = 2;
        }

        Config m_config = new Config();


        public int IterationsInNested {
            get;
            private set;
        }


        public int ThisLevelIterations {
            get;
            private set;
        }

        public bool Converged {
            get;
            private set;
        }
        public Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get;
            set;
        }

        public object Clone() {
            throw new NotImplementedException();
        }

        public void Dispose() {

        }

        public void Init(MultigridOperator op) {
            int D = op.GridData.SpatialDimension;
            int NoOfVars = op.Mapping.NoOfVariables;
            if(NoOfVars != D + 1) {
                throw new NotSupportedException($"expecting {D + 1} variables, got {NoOfVars}");
            }

            // create the sub-multigrid solver
            int[] VelVarIndices = D.ForLoop(d => d);
            var cdMultigridOp = op.GetSubOperator(VelVarIndices);

            // we init the convection-diffusion solver on the top level;
            //  -> this should create a sequence of (recursive) multigrid solvers
            //  -> we pass this sequence along the sequence for the Schur-Complement solver
            var __ConvectionDiffusionSolver = this.m_config.ConvectionDiffusionSolver.CreateInstance(cdMultigridOp);

            InitCoarseRecursive(op, __ConvectionDiffusionSolver, cdMultigridOp);
        }

        MPI_Comm m_comm;

        void InitCoarseRecursive(MultigridOperator op, ISolverSmootherTemplate opSolver, MultigridOperator cdMultigridOp) {
            this.ConvectionDiffusionSolver = opSolver;
            int D = op.GridData.SpatialDimension;
            m_comm = op.OperatorMatrix.MPI_Comm;
            m_iLevel = op.LevelIndex;
            this.PostSmootherArithmeticExceptionCount = 0;
            this.CoarseArithmeticExceptionCount = 0;
            FullOperator = op;

            int[] VelVarIndices = D.ForLoop(d => d);
            int[] PrsVarIndices = new int[] { D };
            VelocityDofIndices = op.Mapping.GetSubvectorIndices(VelVarIndices);
            PressureDofIndices = op.Mapping.GetSubvectorIndices(PrsVarIndices);
            FullSystemPartitioning = op.OperatorMatrix.RowPartitioning;

            // extract the sub-matrices
            PressureGradient = op.OperatorMatrix.GetSubMatrix(VelocityDofIndices, PressureDofIndices);
            VelocityDivergence = op.OperatorMatrix.GetSubMatrix(PressureDofIndices, VelocityDofIndices);
            M22Matrix = op.OperatorMatrix.GetSubMatrix(PressureDofIndices, PressureDofIndices);
            SchurPartitioning = VelocityDivergence.RowPartitioning;

            // init smoother
            m_Smoother = new MatrixFreeChebyshevPreconditioner(this);
            m_Smoother.Init(PressureDofIndices.Length);

            // orthomin on this level
            ortho = new CoreOrthonormalizationProcedure(this.EvaluateSchurComplementMatrix, VelocityDivergence.RowPartitioning);

            // recursion
            var coarseSolver = GetCoarseConvectionDiffusionSolver();
            if(coarseSolver != null) {
                m_Coarse = new SchurComplementForStokes();
                m_Coarse.InitCoarseRecursive(op.CoarserLevel, coarseSolver, cdMultigridOp.CoarserLevel);
            }
        }

        int m_iLevel;

        SchurComplementForStokes m_Coarse;

        MatrixFreeSmoother m_Smoother;

        CoreOrthonormalizationProcedure ortho;

        MultigridOperator FullOperator;

        void SchurRestrict<T1, T2>(T1 IN_fine, T2 OUT_coarse)
            where T1 : IList<double>
            where T2 : IList<double> //
        {
            int LFine = FullOperator.Mapping.LocalLength;
            int LCoarse = FullOperator.CoarserLevel.Mapping.LocalLength;
            double[] IN_fine_full = new double[LFine];
            IN_fine_full.AccVi64(1.0, IN_fine, this.PressureDofIndices, default(long[]), -FullOperator.Mapping.i0, 0);

            double[] OUT_coarse_full = new double[LCoarse];
            FullOperator.CoarserLevel.Restrict(IN_fine_full, OUT_coarse_full);

            OUT_coarse.ClearEntries();
            OUT_coarse.AccVi64(1.0, OUT_coarse_full, default(long[]), m_Coarse.PressureDofIndices, 0, -m_Coarse.FullOperator.Mapping.i0);
        }

        void SchurProlongate<T1, T2>(double alpha, T1 OUT_fine, double beta, T2 IN_coarse)
            where T1 : IList<double>
            where T2 : IList<double> //
        {
            int LFine = FullOperator.Mapping.LocalLength;
            int LCoarse = FullOperator.CoarserLevel.Mapping.LocalLength;
            double[] IN_coarse_full = new double[LCoarse];
            IN_coarse_full.AccVi64(1.0, IN_coarse, m_Coarse.PressureDofIndices, default(long[]), -m_Coarse.FullOperator.Mapping.i0, 0);

            double[] OUT_fine_full = new double[LFine];
            FullOperator.CoarserLevel.Prolongate(1.0, OUT_fine_full, 0.0, IN_coarse_full);

            OUT_fine.ScaleV(beta);
            OUT_fine.AccVi64(alpha, OUT_fine_full, default(long[]), m_Coarse.PressureDofIndices, 0, -m_Coarse.FullOperator.Mapping.i0);
        }


        ISolverSmootherTemplate GetCoarseConvectionDiffusionSolver() {
            var current = this.ConvectionDiffusionSolver;
            if(current is OrthonormalizationMultigrid omg) {
                return omg.CoarserLevelSolver;
            } else {

                return null;
            }
        }

        IPartitioning FullSystemPartitioning;

        IPartitioning SchurPartitioning;

        ISolverSmootherTemplate ConvectionDiffusionSolver;

        long[] VelocityDofIndices;

        long[] PressureDofIndices;


        BlockMsrMatrix PressureGradient;

        BlockMsrMatrix VelocityDivergence;

        /// <summary>
        /// lower right block of the system matrix;
        /// typically, this is zero; otherwise it might contain stabilization terms.
        /// </summary>
        BlockMsrMatrix M22Matrix;


        int PostSmootherArithmeticExceptionCount = 0;
        int CoarseArithmeticExceptionCount = 0;

        void EvaluateSchurComplementMatrix<U, V>(double alpha, V Input, double beta, U Result)
            where U : IList<double>
            where V : IList<double> // 
        {
            EvaluateSchurComplementMatrix<U, V>(alpha, Input, beta, Result, false);
        }

        void EvaluateSchurComplementMatrix<U, V>(double alpha, V Input, double beta, U Result, bool lowRes = false)
            where U : IList<double>
            where V : IList<double> // 
        {
            int L = VelocityDivergence.RowPartitioning.LocalLength;
            if(PressureGradient.ColPartitioning.LocalLength != L)
                throw new ApplicationException("mismatch - number of columns, pressure gradient matrix");
            if(M22Matrix != null) {
                if(M22Matrix.RowPartitioning.LocalLength != L)
                    throw new ApplicationException("mismatch - number of rows, M22 matrix");
                if(M22Matrix.ColPartitioning.LocalLength != L)
                    throw new ApplicationException("mismatch - number of columns, M22 matrix");
            }

            if(Result.Count != L)
                throw new ArgumentException("mismatch, length of input vector");
            if(Input.Count != L)
                throw new ArgumentException("mismatch, length of input vector");

            double[] ResultXbeta = null;
            if(beta != 0.0)
                ResultXbeta = Result.ToArray();

            int LV = VelocityDivergence.ColPartitioning.LocalLength;
            double[] PresGrad = new double[LV];
            double[] CDsol = new double[LV];

            // evaluate pressure gradient
            PressureGradient.SpMV(1.0, Input, 0.0, PresGrad);

            // Solve Conv/Diff System
            Func<int, double, double, (bool bNotTerminate, bool bSuccess)> orgTermCrit = null;
            if(lowRes && ConvectionDiffusionSolver is IProgrammableTermination spt) {
                orgTermCrit = spt.TerminationCriterion;
                double MachineEpsilon = GenericBlas.MachineEps;
                spt.TerminationCriterion = delegate (int i, double r0, double r) {
                    var ret = ((i <= 1 || r > r0 * 0.1) && (r0 > MachineEpsilon || r > MachineEpsilon) && (i <= 10), true);
                    // do at least one iteration until r < 0.1 r0 or max. 10 iterations, with safeguard for tiny residuals
                    return ret;
                };
            }
            try {
                ConvectionDiffusionSolver.Solve(CDsol, PresGrad);
            } finally {
                if(orgTermCrit != null)
                    (ConvectionDiffusionSolver as IProgrammableTermination).TerminationCriterion = orgTermCrit;
            }

            // take velocity divergence
            VelocityDivergence.SpMV(-1.0, CDsol, 0.0, Result);

            // add M22 contribution
            if(M22Matrix != null) {
                M22Matrix.SpMV(1.0, Input, 1.0, Result);
            }

            if(alpha != 1.0)
                Result.ScaleV(alpha);
            if(beta != 0.0)
                Result.AccV(beta, ResultXbeta);
        }


        public void ResetStat() {
            IterationsInNested = 0;
            ThisLevelIterations = 0;
            Converged = false;
        }

        Stopwatch ThisLevelTime = new Stopwatch();
        Stopwatch CrseLevelTime = new Stopwatch();

        public void Solve<U, V>(U _xl, V _B)
            where U : IList<double>
            where V : IList<double> //
        {
            using(new FuncTrace()) {
                var part = FullSystemPartitioning;
                long i0 = part.i0;
                int[] loc_VelocityDofIndices = VelocityDofIndices.Select(glb_idx => checked((int)(glb_idx - i0))).ToArray();
                int[] loc_PressureDofIndices = PressureDofIndices.Select(glb_idx => checked((int)(glb_idx - i0))).ToArray();
                int Lvel = loc_VelocityDofIndices.Length;
                int Lprs = loc_PressureDofIndices.Length;

                var B1 = _B.GetSubVector(loc_VelocityDofIndices);
                var B2 = _B.GetSubVector(loc_PressureDofIndices);

                // build the RHS of the Schur-Complement equation
                var schurRHS = B2; B2 = null;
                var invM11b1 = new double[Lvel];
                this.ConvectionDiffusionSolver.Solve(invM11b1, B1);
                this.VelocityDivergence.SpMV(-1.0, invM11b1, 1.0, schurRHS);

                // solve for the pressure -- the main part of the work
                var pressureSol = new double[Lprs];
                SolveSchurCompl(pressureSol, schurRHS);

                // solve for velocity
                var velocitySol = new double[Lvel];
                var momRhs = B1; B1 = null;
                this.PressureGradient.SpMV(-1.0, pressureSol, 1.0, momRhs);
                this.ConvectionDiffusionSolver.Solve(velocitySol, momRhs);

                // combine velocity and pressure solution
                _xl.SetSubVector(velocitySol, loc_VelocityDofIndices);
                _xl.SetSubVector(pressureSol, loc_PressureDofIndices);
            }
        }

        void SchurResidual(double[] Res, double[] X, double[] B) {
            using(new FuncTrace()) {
                Debug.Assert(Res.Length == SchurPartitioning.LocalLength);
                Debug.Assert(X.Length == SchurPartitioning.LocalLength);
                Debug.Assert(B.Length == SchurPartitioning.LocalLength);
                int L = Res.Length;
                //Res.SetV(B);
                Array.Copy(B, Res, L);
                EvaluateSchurComplementMatrix(-1.0, X, 1.0, Res);
                //Res.AccV(1.0, B);
            }
        }

        public void SolveSchurCompl(double[] X, double[] B) {
            using(var f = new FuncTrace()) {
                ThisLevelTime.Start();

                int L = X.Length;
                double[] Res = new double[L];

                double[] ResCoarse;
                int Lc;
                if(this.m_Coarse != null) {
                    Lc = m_Coarse.SchurPartitioning.LocalLength;
                    ResCoarse = new double[Lc];
                } else {
                    Lc = -1;
                    ResCoarse = null;
                }

                double[] Sol0 = X.CloneAs();
                double[] Res0 = new double[L];
                SchurResidual(Res0, Sol0, B);
                Array.Copy(Res0, Res, L);

                if(ortho.Norm(Res0) <= 0) {
                    double normB = ortho.Norm(B);
                    double normX = ortho.Norm(X);
                    f.Error($" **** Residual is 0.0: |X| = {normB}; |B| = {normB}; |Res0| = {ortho.Norm(Res0)}");


                }

                double iter0_resNorm = ortho.Norm(Res0);
                double resNorm = iter0_resNorm;

                // clear history of coarse solvers
                ortho.Clear();
                bool bIterate = true;
                bool skipPreSmoothInFollowingIters = false;

                WriteDebug(f, 0, resNorm, null);

                MatrixFreeSmoother PreSmoother = this.m_Smoother;
                MatrixFreeSmoother PostSmoother = this.m_Smoother;

                double MachineEpsilon = GenericBlas.MachineEps;

                int iIter;
                for(iIter = 1; bIterate; iIter++) {

                    WriteDebug(f, iIter, resNorm, "initial start");

                    var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if(!termState.bNotTerminate) {
                        Converged = termState.bSuccess;
                        f.Info($"termState bSuccess? {Converged}");
                        break;
                    } else {

                    }



                    // pre-smoother
                    // ------------

                    // Be aware that we have two options for skipping the pre-smoother: one for the following iterations (determined at post smoother), the second for skipping in general from config.
                    if(PreSmoother != null && !skipPreSmoothInFollowingIters) {

                        double[] PreCorr = new double[L];
                        PreSmoother.Solve(PreCorr, Res); // Vorglättung

                        // orthonormalization and residual minimization
                        resNorm = ortho.AddSolAndMinimizeResidual(ref PreCorr, X, Sol0, Res0, Res, "presmoothL" + m_iLevel + "It" + iIter);

                        WriteDebug(f, iIter, resNorm, " pre-smoother applied");

                        //SpecAnalysisSample(iIter, X, "ortho1");
                        var termState2 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                        if(!termState2.bNotTerminate) {
                            Converged = termState2.bSuccess;
                            f.Info($"termState2 bSuccess? {Converged}");
                            break;
                        }

                        skipPreSmoothInFollowingIters = false;
                    }



                    // coarse grid correction
                    // ----------------------
                    if(this.m_Coarse != null) {
                        CrseLevelTime.Start();
                        double resNorm_b4Coarse = resNorm;
                        for(int i = 0; i < m_config.NoOfCoarseSolverSweeps; i++) {


                            // Berechnung der Grobgitterkorrektur
                            SchurRestrict(Res, ResCoarse);
                            double[] vlc = new double[Lc];
                            this.m_Coarse.TerminationCriterion = delegate (int i, double r0, double r) {
                                var ret = ((i <= 1 || r > r0 * 0.1) && (r0 > MachineEpsilon || r > MachineEpsilon) && (i <= 20), true);
                                // do at least one iteration until r < 0.1 r0 or max. 20 iterations, with safeguard for tiny residuals
                                return ret;
                            };

                            this.m_Coarse.SolveSchurCompl(vlc, ResCoarse);

                            // Prolongation der Grobgitterkorrektur
                            double[] vl = new double[L];
                            SchurProlongate(1.0, vl, 1.0, vlc);

                            // orthonormalization and residual minimization
                            if(vl.ContainsNanOrInf().MPIOr() && CoarseArithmeticExceptionCount < 20) {
                                CoarseArithmeticExceptionCount++;
                                f.Error("Coarse solver failed " + CoarseArithmeticExceptionCount);
                                vl = null;
                            } else {
                                resNorm = ortho.AddSolAndMinimizeResidual(ref vl, X, Sol0, Res0, Res, "coarsecorL" + m_iLevel + "It" + iIter);
                            }


                        } // end of coarse-solver loop

                        WriteDebug(f, iIter, resNorm, "coarse-solver applied");

                        var termState3 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                        if(!termState3.bNotTerminate) {
                            Converged = termState3.bSuccess;
                            f.Info($"termState3 bSuccess? {Converged}");
                            break;
                        }
                        CrseLevelTime.Stop();
                    }



                    // post-smoother
                    // -------------
                    if(PostSmoother != null) {
                        bool termPost = false;


                        for(int g = 0; g < m_config.NoOfPostSmootherSweeps; g++) { // Test: Residual on this level / already computed by 'MinimizeResidual' above


                            double[] PostCorr = new double[L];

                            bool fail = false;
                            try {
                                PostSmoother.Solve(PostCorr, Res); // compute correction (Nachglättung)
                            } catch(ArithmeticException ae) {
                                f.Error($"Smoother fail on Rank {FullOperator.DgMapping.MpiRank}: " + ae.ToString());
                                fail = true;
                            }
                            fail = fail.MPIOr();

                            if((PostCorr.ContainsNanOrInf() || fail).MPIOr()) {
                                PostSmootherArithmeticExceptionCount++;
                                f.Error("Post-Smoother produces NAN/INF at sweep " + g + " for " + PostSmootherArithmeticExceptionCount + "-th time.");

                                if(Res.ContainsNanOrInf()) {
                                    Console.WriteLine("... so does RHS");
                                    throw new ArithmeticException("Post-Smoother produces NAN/INF at sweep " + g + " RHS already corrupted.");
                                } else {
                                    Console.WriteLine("... although RHS is regular");
                                }


                            } else {
                                resNorm = ortho.AddSolAndMinimizeResidual(ref PostCorr, X, Sol0, Res0, Res, "pstsmthL" + m_iLevel + "It" + iIter + "-sw" + g);
                                WriteDebug(f, iIter, resNorm, "post-smoother applied");

                                var termState4 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                                if(!termState4.bNotTerminate) {
                                    Converged = termState4.bSuccess;
                                    f.Info($"termState4 bSuccess? {Converged}");
                                    termPost = true;
                                    break;
                                }

                                if(ortho.CancellationTriggered) {
                                    skipPreSmoothInFollowingIters = true; // most of the time, the pre-smoother does nothing different than the post-smoother; So, if we cancel post-smoothing there is no need to do pre-smoothing in the next loop.

                                    break;
                                }
                            }
                        }



                        if(termPost)
                            break;

                    } // end of post-smoother loop

                    // iteration callback
                    // ------------------
                    this.ThisLevelIterations++;

                }// end of solver iterations

                WriteDebug(f, iIter, resNorm, "final");

                ThisLevelTime.Stop();
            } // end of functrace
        }


        /// <summary>
        /// Debug-Output should always go through the <see cref="FuncTrace"/>;
        /// If one would like to see it on screen, turn it on by setting `tr.InfoToConsole = true`.
        /// </summary>
        void WriteDebug(FuncTrace tr, int iter, double res, string text) {
            if(this.m_config.WriteDebugInfoLevel < 2)
                return;

            int iLevel = m_iLevel;

            if(iLevel >= 0)
                tr.InfoToConsole = true; // set to `true` if you want to see output
            tr.Info($"{string.Concat(Enumerable.Repeat("-", iLevel))} SchurComplementForStokes, current level={iLevel}, iteration={iter} {(text != null ? " - " + text : "")} and res norm: {res}");

            return;
        }

        public long UsedMemory() {
            long Memory = 0;

            if(m_iLevel != 0)
                Memory += this.ConvectionDiffusionSolver?.UsedMemory() ?? 0;
            Memory += this.ortho?.GetMemUsage() ?? 0;

            Memory += PressureGradient?.UsedMemory ?? 0;
            Memory += VelocityDivergence?.UsedMemory ?? 0;
            Memory += M22Matrix?.UsedMemory ?? 0;

            return Memory;
        }

        abstract public class MatrixFreeSmoother {
            protected readonly SchurComplementForStokes m_Owner;
            protected int NoOfPressureDOFs;

            public MatrixFreeSmoother(
                SchurComplementForStokes __Owner) {
                m_Owner = __Owner;
            }

            public virtual void Init(int __NoOfPressureDoFs) {
                NoOfPressureDOFs = __NoOfPressureDoFs;
            }


            /// <summary>
            /// Computes an approximation <paramref name="z"/> to the solution `zEx` of the system `A * zEx = ` <paramref name="r"/>.
            /// </summary>
            /// <param name="z">
            /// Output, approximate solution 
            /// </param>
            /// <param name="r">
            /// Input, RHS of the system
            /// </param>            
            abstract public void Solve(
                IList<double> z,
                IList<double> r);

        }


        public class MatrixFreeChebyshevPreconditioner : MatrixFreeSmoother {

            public MatrixFreeChebyshevPreconditioner(SchurComplementForStokes __Owner) : base(__Owner) {
                m_degree = __Owner.m_config.NumberOfIterationsPerSmoother;
            }

            int m_degree = 4;

            public override void Init(int __NoOfPressureDoFs) {
                using(var f = new FuncTrace()) {
                    f.InfoToConsole = m_Owner.m_config.WriteDebugInfoLevel >= 1;
                    base.Init(__NoOfPressureDoFs);

                    double lambdaMaxEstimate = EstimateLargestEigenvalue(iterations: 20);
                    f.Info("MatrixFreeChebyshevPreconditioner, largest eigenvalue " + lambdaMaxEstimate);

                    (m_lambdaMin, m_lambdaMax) =
                        MatrixFreeChebyshevPreconditioner.ChooseSmoothingInterval(
                            lambdaMaxEstimate,
                            lowerFraction: 0.1,
                            safetyFactor: 1.2
                        );
                }
            }

            double m_lambdaMin, m_lambdaMax;

            /// <summary>
            /// Computes an approximation <paramref name="z"/> to the solution `zEx` of the system `A * zEx = ` <paramref name="r"/>.
            /// </summary>
            /// <param name="z">
            /// Output, approximate solution 
            /// </param>
            /// <param name="r">
            /// Input, RHS of the system
            /// </param>
            public override void Solve(IList<double> z, IList<double> r) {
                Solve(z, r, degree: m_degree, m_lambdaMin, m_lambdaMax);
            }


            /// <summary>
            /// Compuztes an approximation <paramref name="z"/> to the solution `zEx` of the system `A * zEx = ` <paramref name="r"/>.
            /// </summary>
            /// <param name="z">
            /// Output, approximate solution 
            /// </param>
            /// <param name="r">
            /// Input, RHS of the system
            /// </param>
            /// <param name="degree"></param>
            /// <param name="lambdaMin"></param>
            /// <param name="lambdaMax"></param>
            public void Solve(
                    IList<double> z,
                    IList<double> r,
                    int degree,
                    double lambdaMin,
                    double lambdaMax) {
                z.ClearEntries();

                double[] residual = new double[NoOfPressureDOFs];
                double[] Sz = new double[NoOfPressureDOFs];

                residual.SetV(r);

                // Chebyshev roots in [lambdaMin, lambdaMax]
                for(int k = 1; k <= degree; ++k) {
                    double theta = Math.PI * (2.0 * k - 1.0) / (2.0 * degree);

                    double lambdaK =
                        0.5 * (lambdaMax + lambdaMin)
                        + 0.5 * (lambdaMax - lambdaMin) * Math.Cos(theta);

                    double omega = 1.0 / lambdaK;

                    // z <- z + omega * residual
                    z.AccV(omega, residual);

                    // residual <- r - Sz
                    Sz.ClearEntries();
                    m_Owner.EvaluateSchurComplementMatrix(1.0, z, 0.0, Sz, true);

                    for(int i = 0; i < NoOfPressureDOFs; ++i)
                        residual[i] = r[i] - Sz[i];
                }
            }

            // Estimates largest eigenvalue using power iteration.
            // Assumes S is approximately symmetric positive definite on the relevant subspace.
            public double EstimateLargestEigenvalue(
                int iterations = 20,
                int seed = 1) {
                var comm = m_Owner.m_comm;
                Random rng = new Random(seed);

                double[] x = new double[NoOfPressureDOFs];
                double[] y = new double[NoOfPressureDOFs];

                for(int i = 0; i < NoOfPressureDOFs; ++i)
                    x[i] = rng.NextDouble() - 0.5;

                x.MPI_Normalize();

                double lambda = 0.0;

                for(int it = 0; it < iterations; ++it) {
                    y.ClearEntries();
                    m_Owner.EvaluateSchurComplementMatrix(1.0, x, 0.0, y, true);

                    lambda = x.MPI_InnerProd(y, m_Owner.m_comm) / x.MPI_L2NormPow2(m_Owner.m_comm);

                    y.MPI_Normalize(m_Owner.m_comm);
                    x.SetV(y);
                }

               return Math.Abs(lambda);
            }

            // Chooses a Chebyshev smoothing interval.
            // This is NOT the true spectral interval; it targets high-frequency modes.
            public static (double lambdaMin, double lambdaMax) ChooseSmoothingInterval(
                double estimatedLambdaMax,
                double lowerFraction = 0.1,
                double safetyFactor = 1.2) {
                double lambdaMax = safetyFactor * estimatedLambdaMax;
                double lambdaMin = lowerFraction * lambdaMax;

                if(lambdaMin <= 0.0)
                    throw new ArgumentException("lambdaMin must be positive.");

                return (lambdaMin, lambdaMax);
            }
        }

    }
}

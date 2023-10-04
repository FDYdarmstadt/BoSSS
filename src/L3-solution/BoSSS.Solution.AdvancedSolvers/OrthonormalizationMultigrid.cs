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
using System.Runtime.Serialization;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using System.Net.NetworkInformation;
using BoSSS.Solution.Tecplot;
using System.Xml.Linq;
using NUnit.Common;
using ilPSP.LinSolvers.monkey;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// 
    /// </summary>
    /// <remarks>
    /// First described in: 
    ///  Vinsome, P.K.W.: 
	///  Orthomin, an Iterative Method for Solving Sparse Sets of Simultaneous Linear Equations,
    ///  SPE Symposium on Numerical Simulation of Reservoir Performance, 
    ///  doi: 10.2118/5729-MS, Los Angeles, California, 1976
    /// </remarks>
    class CoreOrthonormalizationProcedure : ICloneable {


        public CoreOrthonormalizationProcedure(ISparseMatrix __OpMatrix) {
            MxxHistory = new List<double[]>();
            SolHistory = new List<double[]>();
            Alphas = new List<(double, double, string)>();
            OpMatrix = __OpMatrix;
        }

        public void Clear() {
            MxxHistory.Clear();
            SolHistory.Clear();
            Alphas.Clear();
        }


        public long GetMemUsage() {
            int SizeSol = SolHistory != null && SolHistory.Count() > 0 ? this.SolHistory.Count() * this.SolHistory[0].Length * sizeof(double) : 0;
            int SizeMxx = MxxHistory != null && MxxHistory.Count() > 0 ? this.MxxHistory.Count() * this.MxxHistory[0].Length * sizeof(double) : 0;
            int SizeAlpha = Alphas != null && Alphas.Count() > 0 ? this.Alphas.Count() * (sizeof(double) * 2 + sizeof(int)) : 0;
            return SizeAlpha + SizeMxx + SizeSol;
        }


        /// <summary>
        /// 
        /// </summary>
        public ISparseMatrix OpMatrix {
            get;
        }


        /// <summary>
        /// solution guesses from smoothers
        /// </summary>
        public List<double[]> SolHistory = new List<double[]>();

        /// <summary>
        /// - orthonormal system of matrix-vector products;
        /// - the i-th entry is equal to  <see cref="OpMatrix"/>*<see cref="SolHistory"/>[i]
        /// </summary>
        public List<double[]> MxxHistory = new List<double[]>();


        /// <summary>
        /// scaling factors which were applied to <see cref="SolHistory"/> to approximate the solution
        /// - 1st item: the scaling applied to the respective solutions provided through <see cref="AddSol"/>
        /// - 2nd: relative residual reduction (typically greater or equal to 1)
        /// - 3rd: user id data
        /// </summary>
        public List<(double alpha_i, double RelResReduction, string id)> Alphas = new List<(double, double, string)>();


        /// <summary>
        /// 
        /// </summary>
        public double FullMinimization(double[] outX, double[] Sol0, double[] Res0, double[] outRes) {
            //var ids = Alphas.Select(ttt => ttt.id);
            Alphas.Clear();
            double LastResi = -1;
            while(Alphas.Count < SolHistory.Count) {
                LastResi = MinimizeResidual(outX, Sol0, Res0, outRes, "lulu");
                //ids = ids.Skip(1);
            }
            return LastResi;
        }


        /// <summary>
        /// Residual minimization in the space spanned by <see cref="SolHistory"/>/<see cref="MxxHistory"/>
        /// </summary>
        /// <param name="outX">
        /// - input: the solution build upon all vectors in the Krylov space **except the last one**
        /// - output: on exit, updated by the contribution from the most recent Krylov basis
        /// </param>
        /// <param name="Res0">
        /// residual with respect to the initial guess
        /// </param>
        /// <param name="Sol0">
        /// initial guess; required for checking <paramref name="outRes"/>,
        /// resp. for re-computing <paramref name="outRes"/> in the case of severe round-of errors.
        /// </param>
        /// <param name="outRes">
        /// - input: the residual with respect to the input value of <paramref name="outX"/>
        /// - output: on exit, the residual with respect to the output value of <paramref name="outX"/>
        /// </param>
        /// <param name="id">
        /// name/information for debugging and analysis
        /// </param>
        /// <param name="plot"></param>
        private double MinimizeResidual(double[] outX, double[] Sol0, double[] Res0, double[] outRes, string id, bool plot = false) {
            using (var ft = new FuncTrace()) {
                Debug.Assert(SolHistory.Count == MxxHistory.Count);
                Debug.Assert(outX.Length == OpMatrix.ColPartition.LocalLength);
                //Debug.Assert(Sol0.Length == m_MgOperator.Mapping.LocalLength);
                Debug.Assert(Res0.Length == OpMatrix.RowPartitioning.LocalLength);
                Debug.Assert(outRes.Length == OpMatrix.RowPartitioning.LocalLength);



                int KrylovDim = SolHistory.Count;
                int L = outX.Length;

                //if (Alphas.Count != KrylovDim - 1) {
                //    throw new ApplicationException();
                //}

                // note: this implementation is based on the assumption, that the 
                // factors alpha_0 ... alpha_(i-1) are equal to the previous call;
                // therefore, it is only necessary to apply the contributions from the most recent Krylov vector

                var oldResiNorm = Norm(outRes);

                int i = KrylovDim - 1;
                double alpha_i = InnerProduct(MxxHistory[i], Res0);
                BLAS.daxpy(L, alpha_i, SolHistory[i], 1, outX, 1); // accumulate solution correction...
                BLAS.daxpy(L, -alpha_i, MxxHistory[i], 1, outRes, 1); // ...and its effect on the residual



                if (m_UsedMoreThanOneOrthoCycle || KrylovDim % RecalcResiPeriod == 0) { // 40 should be a decent compromise between performance and stability
                  
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // once in a while, re-compute the residual to get rid of round-off errors which might accumulate
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


                    outRes.SetV(Res0);
                    Sol0.AccV(-1.0, outX);
                    this.OpMatrix.SpMV(1.0, Sol0, 1.0, outRes);
                    Sol0.AccV(+1.0, outX);
                }


                double ResNorm = Norm(outRes);
                double tol = Math.Max(ResNorm, oldResiNorm) * 0.2;
                tol = Math.Max(tol, 1.0e-7);
                if (ResNorm > oldResiNorm + tol)
                    throw new ArithmeticException($"residual increase (L = {L}): old norm: {oldResiNorm}, new norm {ResNorm}, reduction factor: {ResNorm / oldResiNorm}");


                if (plot)
                    Console.WriteLine($"reduction factor: {ResNorm / oldResiNorm}: old norm: {oldResiNorm}, new norm {ResNorm}");




                Alphas.Add((alpha_i, oldResiNorm / ResNorm, id));
                return ResNorm;
            }
        }

        public int RecalcResiPeriod = 1;
 

        bool m_UsedMoreThanOneOrthoCycle = false;

        /// <summary>
        /// 
        /// </summary>
        internal (double UsablePart, bool CancellationTriggered) AddSol(ref double[] X, string name) {
            using (var ft = new FuncTrace()) {

                m_UsedMoreThanOneOrthoCycle = false;

                double FillXwithRandom(double[] __X, double[] __Mxx, int iSweep) {
                    double __NormMxx;
                    ft.Error(name + ": Solution norm is exactly 0.0 -- OR -- re-orthonormalization failed; trying with a random vector instead to recover; sweep " + iSweep);
                    __X.FillRandom();

                    __Mxx.ClearEntries();
                    OpMatrix.SpMV(1.0, __X, 0.0, __Mxx);
                    __NormMxx = Norm(__Mxx);

                    if (__NormMxx == 0)
                        throw new ArithmeticException(name + ": Numerical breakdown: norm of matrix * solution after using RANDOM NUMBERS (!) is " + __NormMxx);
                    if (__NormMxx.IsNaNorInf())
                        throw new ArithmeticException(name + ": Numerical breakdown: norm of matrix * solution after using RANDOM NUMBERS (!) is " + __NormMxx);
                    return __NormMxx;
                }

                Debug.Assert(SolHistory.Count == MxxHistory.Count);
                Debug.Assert(X.Length == OpMatrix.RowPartitioning.LocalLength);
                int L = X.Length;
                int KrylovDim = SolHistory.Count;


                double[] Mxx = new double[L];
                OpMatrix.SpMV(1.0, X, 0.0, Mxx);
                double NormMxx = Norm(Mxx);
                if (NormMxx.IsNaNorInf())
                    throw new ArithmeticException(name + ": Numerical breakdown: norm of matrix * solution is " + NormMxx);

                if (NormMxx == 0) {
                    // solution is really strange: try with a random vector.
                    double Xnorm = Norm(X);
                    ft.Error("Problem on Entry into Orthonormalization: NormMxx = " + NormMxx + ", Xnorm = " + Xnorm);
                    NormMxx = FillXwithRandom(X, Mxx, 0);
                }

                // scale Mxx to norm 1, should allow better stability when compared to the Krylov vectors who have all norm 1;
                BLAS.dscal(L, 1.0 / NormMxx, Mxx, 1);
                BLAS.dscal(L, 1.0 / NormMxx, X, 1);

                const int MaxOrtho = 10;
                double MxxRemainderAfterOrthoNorm = 0;
                bool _CancellationTriggered = false;
                for (int jj = 0; jj <= MaxOrtho; jj++) { // re-orthogonalisation, loop-limit to 10; See also book of Saad, p 156, section 6.3.2

                    m_UsedMoreThanOneOrthoCycle = jj > 0;

                    for (int i = 0; i < KrylovDim; i++) {
                        Debug.Assert(!object.ReferenceEquals(Mxx, MxxHistory[i]));
                        double beta = InnerProduct(Mxx, MxxHistory[i]);
                        BLAS.daxpy(L, -beta, SolHistory[i], 1, X, 1);
                        BLAS.daxpy(L, -beta, MxxHistory[i], 1, Mxx, 1);
                    }

                    //double NormAfter = Norm(Mxx);
                    //Console.WriteLine("   orthonormalization norm reduction: " + (NormAfter/NormInitial));

                    double NewMxxNorm = Norm(Mxx);
                    if (jj == 0)
                        MxxRemainderAfterOrthoNorm = NewMxxNorm;

                    //{
                    //    long _L_mpi = ((long)X.Length).MPISum(comm);
                    //    double _Xnorm = Norm(X);
                    //    Console.Out.WriteLine("Orthonormalization (R" + OpMatrix.RowPartitioning.MpiRank + "): " + name + "." + jj +  ", norm after orthogonalization is " + NewMxxNorm + "; norm of X: " + _Xnorm + "; L = " + _L_mpi);
                    //}

                    if (NewMxxNorm <= 1E-5) {
                        using (new BlockTrace("re-orthonormalization", ft)) {
                            _CancellationTriggered = true;
                            if (jj < MaxOrtho) {

                                // a lot of canceling out has occurred; |Mxx| dropped by several magnitudes.
                                // do another loop, to ensure ortho-normality:
                                // Mxx and X lost more than 5 Magnitudes, so we might have lost a couple of digits of accuracy

                                // Note: although |Mxx| might be small, |X| can be quite large (and probably random).
                                // To ensure stability, we must start over with a re-scaled X!
                                // We have to re-scale what is remaining of X:
                                double Xnorm = Norm(X);
                                long L_mpi = ((long)X.Length).MPISum(comm);
                                Console.Out.WriteLine("Orthonormalization (Rank" + OpMatrix.RowPartitioning.MpiRank + "): Severe cancellation after " + name + ", norm after ortho is " + NewMxxNorm + "; norm of X: " + Xnorm + " Doing Re-orthonormalization (" + (jj + 1) + "); L = " + L_mpi);
                                ft.Info("Orthonormalization: Severe cancellation after " + name + ", norm after ortho is " + NewMxxNorm + "; norm of X: " + Xnorm + " Doing Re-orthonormalization (" + (jj + 1) + "); L = " + L_mpi);

                                if ((Xnorm < 1e-200) || (jj == MaxOrtho - 1)) {
                                    // prohibits div by 0, if we got zero solution  
                                    ft.Error("Xnorm = " + Xnorm + " (almost zero solution branch)");
                                    NormMxx = FillXwithRandom(X, Mxx, jj + 1);
                                } else {
                                    X.ScaleV(1.0 / Xnorm);
                                    Mxx.ClearEntries();
                                    OpMatrix.SpMV(1.0, X, 0.0, Mxx);
                                    NormMxx = Norm(Mxx);

                                    if (NormMxx == 0) {
                                        // solution is really strange: try with a random vector.
                                        ft.Error("Xnorm = " + Xnorm + " (exactly zero solution branch)");
                                        NormMxx = FillXwithRandom(X, Mxx, jj + 1);
                                    }
                                }

                                double gamma = 1 / NewMxxNorm;
                                BLAS.dscal(L, gamma, Mxx, 1);
                                BLAS.dscal(L, gamma, X, 1);
                            } else {
                                throw new ArithmeticException(name + ": numerical breakdown");
                            }
                        }
                    } else {
                        double gamma = 1 / NewMxxNorm;
                        BLAS.dscal(L, gamma, Mxx, 1);
                        BLAS.dscal(L, gamma, X, 1);
                        break; // normally, we should terminate after 1 orthonormalization cycle
                    }
                }

                SolHistory.Add(X);
                MxxHistory.Add(Mxx);
                X = null;

                return (MxxRemainderAfterOrthoNorm, _CancellationTriggered);
            }
        }

        MPI_Comm comm {
            get {
                return OpMatrix.MPI_Comm;
            }
        }

        public double[] Diagonal {
            get;
            set;
        }


        /// <summary>
        /// Inner Product use for orthonormalization (not necessarily the standard l2-product)
        /// </summary>
        public double InnerProduct(double[] A, double[] B) {
            if (Diagonal == null) {
                Debug.Assert(A.Length == B.Length);
                return BLAS.ddot(A.Length, A, 1, B, 1).MPISum(comm);
            } else {
                int L = Diagonal.Length; 
                Debug.Assert(A.Length == L);
                Debug.Assert(B.Length == L);

                double acc = 0;
                for(int i = 0; i < L; i++) {
                    if (Diagonal[i] <= 0)
                        throw new ArithmeticException("fucked up diagonal");
                    acc += A[i] * B[i] * Diagonal[i];
                }
                return acc.MPISum(comm);
            }
        }

        /// <summary>
        /// norm induced by <see cref="InnerProduct(double[], double[])"/>
        /// </summary>
        public double Norm(double[] A) {
            if (Diagonal == null) {
                return A.MPI_L2Norm(comm);
            } else {
                return InnerProduct(A, A).Sqrt();
            }
        }


        /// <summary>
        /// Residual minimization in the space spanned by <see cref="SolHistory"/>/<see cref="MxxHistory"/>
        /// </summary>
        /// <param name="outX">
        /// - input: the solution build upon all vectors in the Krylov space **except the last one**
        /// - output: on exit, updated by the contribution from the most recent Krylov basis
        /// </param>
        /// <param name="Res0">
        /// residual with respect to the initial guess
        /// </param>
        /// <param name="Sol0">
        /// initial guess; required for checking <paramref name="outRes"/>,
        /// resp. for re-computing <paramref name="outRes"/> in the case of severe round-of errors.
        /// </param>
        /// <param name="outRes">
        /// - input: the residual with respect to the input value of <paramref name="outX"/>
        /// - output: on exit, the residual with respect to the output value of <paramref name="outX"/>
        /// </param>
        /// <param name="name">
        /// identifier 
        /// </param>
        /// <param name="xGuess">
        /// guess for the solution (typically, preconditioner output); will be modified during the orthonormalization
        /// </param>
        /// <param name="print"></param>
        public double AddSolAndMinimizeResidual(ref double[] xGuess, double[] outX, double[] Sol0, double[] Res0, double[] outRes, string name, bool print = false) {            
            var (UsablePart, _CancellationTriggered) = AddSol(ref xGuess, name);
            if (print) {
                Console.WriteLine($" ................ {name} Usable Part {UsablePart}");
            }
            double newResNorm = MinimizeResidual(outX, Sol0, Res0, outRes, name, print);
            this.CancellationTriggered = _CancellationTriggered;

            //if(name != null)
            //    Console.WriteLine($"   {name} \t\t\t{UsablePart:0.####e-00}\t{Alphas.Last().RelResReduction:0.####e-00}\t{newResNorm:0.####e-00}");


            return newResNorm;
        }

        public object Clone() {
            var R = new CoreOrthonormalizationProcedure(this.OpMatrix);
            R.Alphas.AddRange(this.Alphas);
            R.CancellationTriggered = this.CancellationTriggered;
            R.Diagonal = this.Diagonal?.CloneAs();
            R.MxxHistory.AddRange(this.MxxHistory.Select(a => a.CloneAs())); 
            R.SolHistory.AddRange(this.SolHistory.Select(a => a.CloneAs()));
            R.RecalcResiPeriod = this.RecalcResiPeriod;
            R.m_UsedMoreThanOneOrthoCycle = this.m_UsedMoreThanOneOrthoCycle;

            return R;
        }

        public bool CancellationTriggered;
    }
 
    /// <summary>
    /// A recursive multigrid method, where the convergence (i.e. non-divergence) of each mesh level is guaranteed by an
    /// orthonormalization approach, similar to flexible GMRES, 
    /// as described in:
    ///   BoSSS: A package for multigrid extended discontinuous Galerkin methods; Kummer, Florian; Weber, Jens; Smuda, Martin; Computers &amp; Mathematics with Applications 81
    ///   see https://www.sciencedirect.com/science/article/abs/pii/S0898122120301917?via%3Dihub
    ///   see also https://tubiblio.ulb.tu-darmstadt.de/121465/
    ///   
    /// One instance of this class represents only one level of the multigrid method; coarser 
    /// levels must be added by configuring the <see cref="OrthonormalizationMultigrid.CoarserLevelSolver"/> member.
    /// </summary>
    /// <remarks>
    /// Residual minimization through orthonormalization (ORTHOMIN) is first described in:
    ///  Vinsome, P.K.W.: 
    ///  Orthomin, an Iterative Method for Solving Sparse Sets of Simultaneous Linear Equations,
    ///  SPE Symposium on Numerical Simulation of Reservoir Performance, 
    ///  doi: 10.2118/5729-MS, Los Angeles, California, 1976
    /// </remarks>    
    public class OrthonormalizationMultigrid : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination, ISubsystemSolver {


        /// <summary>
        /// Individual configuration of <see cref="OrthonormalizationMultigrid"/>
        /// </summary>
        [DataContract]
        [Serializable]
        public class Config : IterativeSolverConfig {


            /// <summary>
            /// Ctor
            /// </summary>
            public Config() {
            }



            /// <summary>
            /// config cctor of <see cref="OrthonormalizationMultigrid"/>
            /// </summary>
            public int MaxKrylovDim = int.MaxValue;

            /// <summary>
            /// - True: the default value: <see cref="OrthonormalizationMultigrid.CoarserLevelSolver"/> is initialized and solved on coarser level
            /// - false: <see cref="OrthonormalizationMultigrid.CoarserLevelSolver"/> is initialized on the same level, but it may perform tis own restriction
            /// </summary>
            [DataMember]
            public bool CoarseOnLovwerLevel = true;

            /// <summary>
            /// - if set to 1, a this performs a V-cycle
            /// - if set to 2 or higher, a W-cycle, resp. WW-cycle, WWW-cycle, etc.
            /// </summary>
            [DataMember]
            public int m_omega = 1;

            /// <summary>
            /// 
            /// </summary>
            override public string Name {
                get { return "OrthonormalizationMultigrid"; }
            }

            /// <summary>
            /// 
            /// </summary>
            override public string Shortname {
                get { return "OrthMG"; }
            }

            /// <summary>
            /// factory
            /// </summary>
            override public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                var instance = new OrthonormalizationMultigrid();
                instance.myConfig = this;
                instance.Init(level);
                instance.TerminationCriterion = base.DefaultTermination;
                return instance;
            }

            /// <summary>
            /// 
            /// </summary>
            public int NoOfPostSmootherSweeps = 2;
        }


        Config myConfig;

        /// <summary>
        /// Solver configuration
        /// </summary>
        public Config config {
            get {
                return myConfig;
            }
        }

        static int counter_objectId = 1;

        public int objectId;

        /// <summary>
        /// ctor
        /// </summary>
        public OrthonormalizationMultigrid() {
            myConfig = new Config();
            TerminationCriterion = myConfig.DefaultTermination;
            objectId = counter_objectId;
            counter_objectId++;
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
        public BlockMsrMatrix OpMatrix => m_OpMapPair.OperatorMatrix;

        /// <summary>
        /// passed in <see cref="InitImpl"/>
        /// </summary>
        IOperatorMappingPair m_OpMapPair;


        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(IOperatorMappingPair op) {
            InitImpl(op);
        }

        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(MultigridOperator op) {
            InitImpl(op);

        }


        CoreOrthonormalizationProcedure ortho;

        void InitImpl(IOperatorMappingPair op) {
            using (var tr = new FuncTrace()) {
                if (object.ReferenceEquals(op, m_OpMapPair))
                    return; // already initialized
                else
                    this.Dispose();


                this.m_OpMapPair = op;
                var Mtx = op.OperatorMatrix;
                var MgMap = op.DgMapping;
                
                if (!Mtx.RowPartitioning.EqualsPartition(MgMap))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!Mtx.ColPartition.EqualsPartition(MgMap))
                    throw new ArgumentException("Column partitioning mismatch.");

              

                ortho = new CoreOrthonormalizationProcedure(Mtx);


                // set operator
                // ============

                // initiate coarser level
                // ======================
                if (this.CoarserLevelSolver == null) {
                    //throw new NotSupportedException("Missing coarse level solver.");
                    tr.Info("OrthonormalizationMultigrid: running without coarse solver.");
                } else {
                    if (op is MultigridOperator mgOp) {
                        if (myConfig.CoarseOnLovwerLevel && mgOp.CoarserLevel != null) {
                            this.CoarserLevelSolver.Init(mgOp.CoarserLevel);
                        } else {
                            tr.Info("OrthonormalizationMultigrid: running coarse solver on same level.");
                            this.CoarserLevelSolver.Init(mgOp);
                        }
                    } else {
                        if(myConfig.CoarseOnLovwerLevel == false && this.CoarserLevelSolver is ISubsystemSolver ssCoarse) {
                            ssCoarse.Init(op);
                        } else {
                            throw new NotSupportedException($"Unable to initialize coarse-level-solver if operator is not a {typeof(MultigridOperator)}");
                        }
                    }
                }

                // init smoother
                // =============
                if (PreSmoother != null) {
                    if (PreSmoother is ISubsystemSolver ssPreSmother) {
                        ssPreSmother.Init(op);
                    } else {
                        if (op is MultigridOperator mgOp) {
                            PreSmoother.Init(mgOp);
                        } else {
                            throw new NotSupportedException($"Unable to initialize pre-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                        }
                    }
                }
                if (PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother)) {
                    //PostSmoother.Init(op);
                    if (PostSmoother is ISubsystemSolver ssPostSmother) {
                        ssPostSmother.Init(op);
                    } else {
                        if (op is MultigridOperator mgOp) {
                            PostSmoother.Init(mgOp);
                        } else {
                            throw new NotSupportedException($"Unable to initialize post-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                        }
                    }
                }
                
            }
        }

        bool m_AdditionalPostSmoothersInitialized = false;

        /// <summary>
        /// Deferred initialization of the <see cref="AdditionalPostSmoothers"/>;
        /// these are fall-back-options, so on a normal run, we don't want to waste time and resources for their initialization.
        /// </summary>
        /// <exception cref="NotSupportedException"></exception>
        void InitAdditionalPostSmooters() {
            using (new FuncTrace()) {
                if (AdditionalPostSmoothers != null) {
                    var op = this.m_OpMapPair;
                    foreach (var ps in AdditionalPostSmoothers) {
                        if (ps is ISubsystemSolver ssPostSmother) {
                            ssPostSmother.Init(op);
                        } else {
                            if (op is MultigridOperator mgOp) {
                                ps.Init(mgOp);
                            } else {
                                throw new NotSupportedException($"Unable to initialize post-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                            }
                        }
                    }

                    m_AdditionalPostSmoothersInitialized = true;
                }
            }
        }

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
        /// high frequency solver after coarse grid correction
        /// </summary>
        public ISolverSmootherTemplate PostSmoother;

        /// <summary>
        /// high frequency solver after coarse grid correction
        /// </summary>
        public ISolverSmootherTemplate[] AdditionalPostSmoothers;

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
                Debug.Assert(Res.Length == m_OpMapPair.DgMapping.LocalLength);
                Debug.Assert(X.Length == m_OpMapPair.DgMapping.LocalLength);
                Debug.Assert(B.Length == m_OpMapPair.DgMapping.LocalLength);
                int L = Res.Length;
                //Res.SetV(B);
                Array.Copy(B, Res, L);
                OpMatrix.SpMV(-1.0, X, 1.0, Res);
                //Res.AccV(1.0, B);


            }
        }

        static double Inner(MultigridOperator mgOp, double[] a, double[] b) {
            int L = mgOp.Mapping.LocalLength;
            if (a.Length != L)
                throw new ArgumentException();
            if(b.Length != L)
                throw new ArgumentException();

            if(mgOp.MassMatrix == null) {
                return a.MPI_InnerProd(b, mgOp.OperatorMatrix.MPI_Comm);
            } else {
                double[] Mb = new double[L];
                mgOp.MassMatrix.SpMV(1.0, b, 0.0, Mb);
                return a.MPI_InnerProd(Mb, mgOp.MassMatrix.MPI_Comm);
            }


        }


        static public BlockMsrMatrix AgglomMassMatrix;

        /// <summary>
        /// the multigrid iterations for a linear problem
        /// </summary>
        /// <param name="_xl">on input, the initial guess; on exit, the result of the multigrid iteration</param>
        /// <param name="_B">the right-hand-side of the problem</param>
        public void Solve<U, V>(U _xl, V _B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (var f = new FuncTrace()) {
                ThisLevelTime.Start();
                double[] B, X;
                if (_B is double[])
                    B = _B as double[];
                else
                    B = _B.ToArray();
                if (_xl is double[])
                    X = _xl as double[];
                else
                    X = _xl.ToArray();

                int L = X.Length;
                double[] Res = new double[L];

                double[] ResCoarse;
                int Lc;
                if (this.CoarserLevelSolver != null && myConfig.CoarseOnLovwerLevel) {
                    Lc = ((MultigridOperator)m_OpMapPair).CoarserLevel.Mapping.LocalLength;
                    ResCoarse = new double[Lc];
                } else {
                    Lc = -1;
                    ResCoarse = null;
                }



                double[] Sol0 = X.CloneAs();
                double[] Res0 = new double[L];
                Residual(Res0, Sol0, B);
                Array.Copy(Res0, Res, L);

                if(ortho.Norm(Res0) <= 0) {
                    double normB = ortho.Norm(B);
                    double normX = ortho.Norm(X);
                    f.Error($" **** Residual is 0.0: |X| = {normB}; |B| = {normB}; |Res0| = {ortho.Norm(Res0)}");


                }

                int iLevel = ((m_OpMapPair as MultigridOperator)?.LevelIndex ?? -1);

                double iter0_resNorm = ortho.Norm(Res0);
                double resNorm = iter0_resNorm;
                this.IterationCallback?.Invoke(0, Sol0, Res0, this.m_OpMapPair as MultigridOperator);
                                
                // clear history of coarse solvers
                ortho.Clear();
                bool bIterate = true;

                bool skipPreSmooth = false;


                ISolverSmootherTemplate[] allSmooters;
                int iPostSmooter;
                {
                    allSmooters = new ISolverSmootherTemplate[(PostSmoother != null ? 1 : 0) + (AdditionalPostSmoothers?.Length ?? 0)];
                    if (PostSmoother != null)
                        allSmooters[allSmooters.Length - 1] = PostSmoother;
                    if (AdditionalPostSmoothers != null)
                        Array.Copy(AdditionalPostSmoothers, allSmooters, AdditionalPostSmoothers.Length);
                    iPostSmooter = allSmooters.Length - 1;
                }


                int iIter;
                for (iIter = 1; bIterate; iIter++) {
                    var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if (!termState.bNotTerminate) {
                        Converged = termState.bSuccess;
                        break;
                    } else {

                    }

                    // pre-smoother
                    // ------------

                    {
                        if (PreSmoother != null && skipPreSmooth == false) {
                            VerivyCurrentResidual(X, B, Res, iIter);

                            double[] PreCorr = new double[L];
                            PreSmoother.Solve(PreCorr, Res); // Vorglättung

                            // orthonormalization and residual minimization
                            resNorm = ortho.AddSolAndMinimizeResidual(ref PreCorr, X, Sol0, Res0, Res, "presmoothL" + iLevel);

                            //SpecAnalysisSample(iIter, X, "ortho1");
                            var termState2 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                            if (!termState2.bNotTerminate) {
                                Converged = termState2.bSuccess;
                                break;
                            }

                            skipPreSmooth = false;
                        }
                    }

                    // coarse grid correction
                    // ----------------------
                    CrseLevelTime.Start();
                    // Test: Residual on this level / already computed by 'MinimizeResidual' above
                    VerivyCurrentResidual(X, B, Res, iIter);

                    double resNorm_b4Coarse = resNorm;
                    for (int i = 0; i < myConfig.m_omega; i++) {
                        if (this.CoarserLevelSolver != null) {

                            double[] vl = new double[L];
                            if (myConfig.CoarseOnLovwerLevel) {

                                var _MgOperator = m_OpMapPair as MultigridOperator;

                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // coarse grid solver defined on COARSER MESH LEVEL:
                                // this solver must perform restriction and prolongation
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                using (new BlockTrace("Restriction", f)) {
                                    // restriction of residual
                                    _MgOperator.CoarserLevel.Restrict(Res, ResCoarse);
                                }
                                // Berechnung der Grobgitterkorrektur
                                double[] vlc = new double[Lc];
                                this.CoarserLevelSolver.Solve(vlc, ResCoarse);
                                using (new BlockTrace("Prolongation", f)) {
                                    // Prolongation der Grobgitterkorrektur
                                    _MgOperator.CoarserLevel.Prolongate(1.0, vl, 1.0, vlc);
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
                            resNorm = ortho.AddSolAndMinimizeResidual(ref vl, X, Sol0, Res0, Res, "coarsecorL" + iLevel);

                        }
                    } // end of coarse-solver loop

                    var termState3 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if (!termState3.bNotTerminate) {
                        Converged = termState3.bSuccess;
                        break;
                    }
                    CrseLevelTime.Stop();

                    // post-smoother
                    // -------------
                    if (PostSmoother != null || AdditionalPostSmoothers != null) {
                        bool termPost = false;

                        

                        for (int g = 0; g < config.NoOfPostSmootherSweeps; g++) { // Test: Residual on this level / already computed by 'MinimizeResidual' above

                            ISolverSmootherTemplate _PostSmoother = allSmooters[iPostSmooter];// g % allSmooters.Length];
                            if (iPostSmooter > 0)
                                InitAdditionalPostSmooters();
                            VerivyCurrentResidual(X, B, Res, iIter); // 
                            double[] PostCorr = new double[L];

                            _PostSmoother.Solve(PostCorr, Res); // compute correction (Nachglättung)
                            

                            if (PostCorr.ContainsForNanOrInfV()) {
                                Console.Error.WriteLine("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g);

                                if (Res.ContainsForNanOrInfV())
                                    Console.WriteLine("... so does RHS");
                                else
                                    Console.WriteLine("... although RHS is regular");

                                //var viz = new MGViz(m_OpMapPair as MultigridOperator);
                                //var __RHS = viz.ProlongateRhsToDg(Res, "RES");
                                //var __SOL = viz.ProlongateRhsToDg(PostCorr, "SOL");
                                //Tecplot.Tecplot.PlotFields(__SOL.Cat(__RHS), "iilufail", 0.0, 0);
                                throw new ArithmeticException("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g);


                            } else {
                                resNorm = ortho.AddSolAndMinimizeResidual(ref PostCorr, X, Sol0, Res0, Res, "pstsmthL" +  iLevel + "-sw" + g);



                                var termState4 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                                if (!termState4.bNotTerminate) {
                                    Converged = termState4.bSuccess;
                                    termPost = true;
                                    break;
                                }

                                if (ortho.CancellationTriggered) {
                                    skipPreSmooth = true; // most of the time, the pre-smoother does nothing different than the post-smoother; So, if we cancel post-smoothing there is no need to do pre-smoothing in the next loop.

                                    iPostSmooter++;
                                    if(iPostSmooter >= allSmooters.Length)
                                        iPostSmooter = 0;

                                    
                                   

                                    break;
                                }
                            }
                        }
                        

                        if (termPost)
                            break;

                    } // end of post-smoother loop


                    
                    // iteration callback
                    // ------------------
                    this.ThisLevelIterations++;
                    IterationCallback?.Invoke(iIter, X, Res, this.m_OpMapPair as MultigridOperator);

                } // end of solver iterations

                IterationCallback?.Invoke(iIter, X, Res, this.m_OpMapPair as MultigridOperator);


                // solution copy
                // =============
                if (!ReferenceEquals(_xl, X)) {
                    _xl.SetV(X);
                }
                ThisLevelTime.Stop();
            } // end of functrace
        }


      
        /// <summary>
        /// For performance optimization, the <see cref="OrthonormalizationMultigrid"/>
        /// assumes that <see cref="PreSmoother"/> and <see cref="PostSmoother"/>
        /// update the residual on exit.
        /// </summary>
        private void VerivyCurrentResidual(double[] X, double[] B, double[] Res, int iter) {
            using (var tr = new FuncTrace()) {
                /*
#if DEBUG
            {
#else
                if (iter % 20 == 0 && iter > 1) {
#endif
                    double[] rTest = new double[Res.Length];
                    Residual(rTest, X, B); // Residual on this level; 
                                           // Test also fails if convergence criterium is to strict because then machine accuracy is reached
                                           //Console.WriteLine("verified Residual: " + resDist);
                    var ss = (new double[] { rTest.L2DistPow2(Res), Res.L2NormPow2(), X.L2NormPow2() }).MPISum(m_MgOperator.DgMapping.MPI_Comm);
                    double resDist = ss[0].Sqrt();
                    double resNormTst = ss[1].Sqrt();
                    double XnormTest = ss[2].Sqrt();
                    tr.Info($"Residual vector check iter {iter}: distance is {resDist}, reference value {resNormTst}");
                    if (resDist > resNormTst * 10e-5 + XnormTest * 1e-5)
                        throw new ArithmeticException($"Residual vector (after pre-smoother/before coarse-correction) is not up-to-date: distance is {resDist}, reference value {resNormTst}");
                    //Debug.Assert(resDist <= resNormTst * 10e-5, $"Residual vector is not up-to-date: distance is {resDist}, reference value ${resNormTst}");
                }

                */
            }
        }

      
        Stopwatch ThisLevelTime = new Stopwatch();
        Stopwatch CrseLevelTime = new Stopwatch();


        /// <summary>
        /// ~
        /// </summary>
        public int IterationsInNested {
            get {
                int iter = 0;

                iter += (this.PreSmoother?.IterationsInNested ?? 0) + (this.PreSmoother?.ThisLevelIterations ?? 0);

                if (this.PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother))
                    iter += this.PostSmoother.IterationsInNested + this.PostSmoother.ThisLevelIterations;

                iter += (this.CoarserLevelSolver?.IterationsInNested ?? 0) + (this.CoarserLevelSolver?.ThisLevelIterations ?? 0);

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
            if (PreSmoother != null) Memory += PreSmoother.UsedMemory();
            if (PostSmoother != null) Memory += PostSmoother.UsedMemory();
            return Memory;
        }

        public long MemoryOfMultigrid() {
            long Memory = 0;
            if (this.CoarserLevelSolver is OrthonormalizationMultigrid)
                Memory += (this.CoarserLevelSolver as OrthonormalizationMultigrid).MemoryOfMultigrid();

            Memory += ortho?.GetMemUsage() ?? 0;
            return Memory;
        }

        public void Dispose() {
            //if (m_MTracker != null) m_MTracker.Dispose();
            if (m_verbose && m_OpMapPair != null && m_OpMapPair is MultigridOperator _mgop) {
                int lv = _mgop.LevelIndex;
                Console.WriteLine($"OrthoMG lv {lv} - total memory: {UsedMemory() / (1024 * 1024)} MB");
                Console.WriteLine($"OrthoMG lv {lv} - internal memory: {MemoryOfMultigrid() / (1024 * 1024)} MB");
                Console.WriteLine($"OrthoMG lv {lv} - smoother memory: {MemoryOfSmoother() / (1024 * 1024)} MB");

                Console.WriteLine($"OrthoMG lv {lv} - total runtime: {ThisLevelTime.Elapsed.TotalSeconds} sec");
                Console.WriteLine($"OrthoMG lv {lv} - coarse runtime: {CrseLevelTime.Elapsed.TotalSeconds:F1} sec ({100*CrseLevelTime.Elapsed.TotalSeconds/ThisLevelTime.Elapsed.TotalSeconds:F1})");
            }
            if (this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.Dispose();
            //this.CoarserLevelSolver = null; // don't delete - we need this again for the next init

            ortho?.Clear();
            ortho = null;
            this.m_OpMapPair = null;

            this.PreSmoother?.Dispose();
            this.PostSmoother?.Dispose();
            //this.PreSmoother = null; // don't delete - we need this again for the next init
            //this.PostSmoother = null;  // don't delete - we need this again for the next init

            if(AdditionalPostSmoothers != null && m_AdditionalPostSmoothersInitialized) {
                foreach(var aps in AdditionalPostSmoothers)
                    aps?.Dispose();
                m_AdditionalPostSmoothersInitialized = false;
            }
        }


    }
}

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
using System.Runtime.InteropServices;
using System.Threading;
using ilPSP.Kraypis;
using static System.Reflection.Metadata.BlobBuilder;
using ilPSP.LinSolvers.PARDISO;
using System.Threading.Tasks;
using Newtonsoft.Json.Linq;

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
    class CoreOrthonormalizationProcedureTP : ICloneable {


        public CoreOrthonormalizationProcedureTP(ISparseMatrix __OpMatrix) {
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
                    ft.Info($"reduction factor: {ResNorm / oldResiNorm}: old norm: {oldResiNorm}, new norm {ResNorm}");




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
        public double AddSolAndMinimizeResidual(ref double[] xGuess, double[] outX, double[] Sol0, double[] Res0, double[] outRes, string name) {            
            var (_, _CancellationTriggered) = AddSol(ref xGuess, name);
            
            double newResNorm = MinimizeResidual(outX, Sol0, Res0, outRes, name, false);
            this.CancellationTriggered = _CancellationTriggered;

            //if(name != null)
            //    Console.WriteLine($"   {name} \t\t\t{UsablePart:0.####e-00}\t{Alphas.Last().RelResReduction:0.####e-00}\t{newResNorm:0.####e-00}");


            return newResNorm;
        }

        public object Clone() {
            var R = new CoreOrthonormalizationProcedureTP(this.OpMatrix);
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
    public class TaskParallelOrthoMG : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination, ISubsystemSolver {


		/// <summary>
		/// Individual configuration of <see cref="TaskParallelOrthoMG"/>
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
            public bool CoarseOnLowerLevel = true;

            /// <summary>
            /// - if set to 1, a this performs a V-cycle
            /// - if set to 2 or higher, a W-cycle, resp. WW-cycle, WWW-cycle, etc.
            /// </summary>
            [DataMember]
            public int m_omega = 1;

            /// <summary> skip the pre-smoother </summary>
            [DataMember]
            public bool SkipPreSmoother = false;

            /// <summary> Pre-smoother and coarse grid correction do not work sequential (i.e., residual from presmoother is not supplied to coarse grid solver) </summary>
            [DataMember]
            public bool NonSerialPreSmoother = false;

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
                var instance = new TaskParallelOrthoMG();
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
        public TaskParallelOrthoMG() {
            myConfig = new Config();
            TerminationCriterion = myConfig.DefaultTermination;
            worldCommRank = GetRank(csMPI.Raw._COMM.WORLD);
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
            if (opCommRestrictionOperator == null || opCommProlongationOperator == null)
                throw new ArgumentException("Restriction and prolongation operator not set. Are you sure that the operator is called correctly?");

			InitImpl(op);
        }

        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(MultigridOperator op) {
			opCommRestrictionOperator = op.GetRestrictionOperator;
			opCommProlongationOperator = op.GetPrologonationOperator;
			opLeftChangeOfBasisMatrix = op.LeftChangeOfBasis;
			opRightChangeOfBasisMatrix = op.RightChangeOfBasis;

            if (op.OperatorMatrix.MPI_Comm != csMPI.Raw._COMM.WORLD)
                throw new Exception("Task parallel OrthoMG (finest level) should be initiated with an operator in world communicator");

			InitImpl(op);
        }

        public void CommitCoarseningOperators(BlockMsrMatrix restrictionOperator, BlockMsrMatrix prolongationOperator) {
            if (myTask == parallelTask.Smoother)
                return;

            return;
			// nothing to do
		}
		// ? a stack of restriction operators for coarse levels ?

		BlockMsrMatrix opCommRestrictionOperator = null;
		BlockMsrMatrix opCommProlongationOperator = null;
		BlockMsrMatrix opLeftChangeOfBasisMatrix = null;
		BlockMsrMatrix opRightChangeOfBasisMatrix = null;
		BlockMsrMatrix subCommSmootherOpMatrix = null;
		BlockMsrMatrix subCommCoarseOpMatrix = null;
		BlockMsrMatrix smootherPermutation = null;
		BlockMsrMatrix coarsePermutation = null;

		List<(long, long)> columnMappingOpToSmoother = new List<(long, long)>(); //cell-based indexes for the mapping from the operator to the smoother matrix (source idx, target idx)
		List<(long, long)> columnMappingOpToCoarse = new List<(long, long)>();   // cell-based indexes for the mapping from the operator to the coarse matrix (source idx, target idx)

		CoreOrthonormalizationProcedureTP ortho;

		private int? m_NoOfCoarseProcs;
		private int? m_NoOfSmootherProcs;

		public int NoOfCoarseProcs {
			get => m_NoOfCoarseProcs ?? throw new InvalidOperationException("Value not yet assigned.");
			set {
				if (m_NoOfCoarseProcs.HasValue)
					throw new InvalidOperationException("Value already assigned and cannot be changed.");
				m_NoOfCoarseProcs = value;
			}
		}

		public int NoOfSmootherProcs {
			get => m_NoOfSmootherProcs ?? throw new InvalidOperationException("Value not yet assigned.");
			set {
				if (m_NoOfSmootherProcs.HasValue)
					throw new InvalidOperationException("Value already assigned and cannot be changed.");
				m_NoOfSmootherProcs = value;
			}
		}


		int GetNumberOfCoarseBlocks() {
            return opCommSize / 4;
		}

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

                ortho = new CoreOrthonormalizationProcedureTP(Mtx);

                int NoOfCoarseBlocks = GetNumberOfCoarseBlocks();
                SplitCommunicator(NoOfCoarseBlocks);
                CommitCoarseningOperators(opCommRestrictionOperator, opCommProlongationOperator);
                CreateCoarseOpMatrix();

				// set operator
				// ============

				// initiate coarser level
				// ======================
				InitCoarse();

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

		ICoordinateMapping MgMap => m_OpMapPair.DgMapping;

        void CreateCoarseOpMatrix() {
			using (var tr = new FuncTrace()) {
				//if (opCommRestrictionOperator == null || opCommProlongationOperator == null)
				//	throw new ArgumentException("Restriction and prolongation operator not set. Are you sure that the operator is called correctly?");

                //mapping
                ChangeCommunicators();
				// 


				// create coarse level
				// ===================
				//m_OpMapPair = new OperatorMappingPair(opCommRestrictionOperator, opCommProlongationOperator, opLeftChangeOfBasisMatrix, opRightChangeOfBasisMatrix);
			}
		}

		void ChangeCommunicators() {
            if (MgMap is MultigridMapping mapping) {
				//smoother (should be called by all procs in opComm)
                var (distForSmoother, columnMappingOpToSmoother) = DistributeMapping(mapping, NoOfSmootherProcs, 0);
				var GlobDOFidx = CellIndexToDOFs(mapping, distForSmoother);
				var Permutation = GetPermutationMatrix(mapping, GlobDOFidx); //this is technically a permutation matrix
				var LocalizedMatrix = BlockMsrMatrix.Multiply(Permutation.Matrix, OpMatrix); // transfer matrix to the target procs (changed the position of rows not columns)
				LocalizedMatrix.ChangeColumnIndices(columnMappingOpToSmoother); //switches column indices with respect to the mapping
				LocalizedMatrix.ChangeCommPattern(subComm, LocalizedMatrix._RowPartitioning); // calculates new external indices with respec to the new partitioning (still opComm)
				subCommSmootherOpMatrix = LocalizedMatrix.GetSubMatrix(Permutation.RowIndices, Permutation.RowIndices, subComm); // finally create a new matrix with the new Comm
				smootherPermutation = Permutation.Matrix;

				// coarse (should be called by all procs in opComm)
				var (distForCoarse, columnMappingOpToCoarse) = DistributeMapping(mapping, NoOfCoarseProcs, NoOfSmootherProcs); // offset is set according to remaning elements = (0 + NoOfSmootherProcs) 
				GlobDOFidx = CellIndexToDOFs(mapping, distForCoarse);
				Permutation = GetPermutationMatrix(mapping, GlobDOFidx); //this is technically a permutation matrix
				LocalizedMatrix = BlockMsrMatrix.Multiply(Permutation.Matrix, OpMatrix); // transfer matrix to the target procs (changed the position of rows not columns)
				LocalizedMatrix.ChangeColumnIndices(columnMappingOpToCoarse); //switches column indices with respect to the mapping
				LocalizedMatrix.ChangeCommPattern(subComm, LocalizedMatrix._RowPartitioning); // calculates new external indices with respec to the new partitioning (still opComm)
				subCommCoarseOpMatrix = LocalizedMatrix.GetSubMatrix(Permutation.RowIndices, Permutation.RowIndices, subComm); // finally create a new matrix with the new Comm
				coarsePermutation = Permutation.Matrix;

#if DEBUG
				TestMatrices(distForSmoother, GlobDOFidx, mapping, smootherPermutation, subCommSmootherOpMatrix, "s_", verbose);
				TestMatrices(distForCoarse, GlobDOFidx, mapping, coarsePermutation, subCommCoarseOpMatrix, "c_", verbose);
#endif
			} else {
                throw new NotImplementedException();
            }
		}

		void TestMatrices(IEnumerable<int> dist, (long i0Cell, int lenCell)[] GlobDOFidx, MultigridMapping mapping, BlockMsrMatrix PermutationMtx, BlockMsrMatrix subMtx, string tag = "t_", bool WriteToFiles = false) {
			//Create a test vector
			Random rnd = new Random(44); //seed is 44
			var TestVector = new double[mapping.LocalLength];
			for (int i = 0; i < TestVector.Length; i++) {
				TestVector[i] = rnd.NextDouble();
			}

			// first test SpMV with opCommMatrix
			var opResult = new double[mapping.LocalLength];
			OpMatrix.SpMV(1.0, TestVector, 0.0, opResult);
				
			// test vector to subComm
			var TestVectorSub = new double[PermutationMtx._RowPartitioning.LocalLength];
			PermutationMtx.SpMV(1.0, TestVector, 0.0, TestVectorSub);

			// result at the subComm
			var subOpResult = new double[PermutationMtx._RowPartitioning.LocalLength];
			subMtx.SpMV(1.0, TestVectorSub, 0.0, subOpResult);

			// map back the subComm result to the opComm
			var TransposeRedist = PermutationMtx.Transpose();
			var subOpResultBackToOp = new double[TransposeRedist._RowPartitioning.LocalLength];
			TransposeRedist.SpMV(1.0, subOpResult, 0.0, subOpResultBackToOp);

			if (WriteToFiles) {
				dist.SaveToTextFileDebug($"{tag}distForSmoother.txt");
				GlobDOFidx.SaveToTextFileDebug($"{tag}GlobDOFidx.txt");

				OpMatrix.SaveToTextFileSparseDebug($"{tag}OpMatrix.txt"); //parallel
				OpMatrix.SaveToTextFileSparse($"{tag}OpMatrix.txt"); /// combine and write as a single file

				TestVector.SaveToTextFileDebug($"{tag}TestVector", ".txt");
				opResult.SaveToTextFileDebug($"{tag}opResult", ".txt");

				subMtx.SaveToTextFileSparseDebug($"{tag}localBlock.txt");
				subMtx.SaveToTextFileSparse($"{tag}localBlock.txt");
				
				PermutationMtx.SaveToTextFileSparseDebug($"{tag}Redist.txt");
				PermutationMtx.SaveToTextFileSparse($"{tag}Redist.txt");

				subOpResult.SaveToTextFileDebug($"{tag}subOpResult", ".txt");
				subOpResultBackToOp.SaveToTextFileDebug($"{tag}subOpResultBackToOp", ".txt");
			}

			for (int i = 0; i < subOpResultBackToOp.Length; i++) {
				if (Math.Abs(opResult[i] - subOpResultBackToOp[i]) > Math.Pow(10, -12))
					throw new Exception("Something odd with redistribution matrix, the solutions do not hold");
			}
		}

		/// <summary>
		/// Creates the DOF tuples for the local cells and exchange them with the other processors
		/// </summary>
		/// <param name="map">current multigrid mapping at this level</param>
		/// <param name="targetProcs">target processor for each cell (index: local cell index, value: proc rank on opComm)</param>
		/// <returns>global DOF data distributed w.r.t. <paramref name="targetProcs"/>.</returns>
		(long i0Cell, int lenCell)[] CellIndexToDOFs(MultigridMapping map, int[] targetProcs) {
			IGridData g = map.AggGrid;
			int NoOfCells = targetProcs.Length;
			int J = g.iLogicalCells.NoOfLocalUpdatedCells;
			Debug.Assert(J == NoOfCells); // 

			var myList = new (long i0Cell, int lenCell)[J]; //for cell currently on this proc (may be distributed to another or not)

			// Create matrix entry information from cell indices
			for (int jCellLoc = 0; jCellLoc < J; jCellLoc++) {
                myList[jCellLoc] = (map.GetCellI0(jCellLoc), map.GetLength(jCellLoc));
			}

			myList.SaveToTextFileDebug("t_listLocalOfDOfs");
			//exchange the data with processors depending on the list targetProcs
			var rcvData = ExchangeDOFData(myList, targetProcs); 

			return rcvData;
		}

		/// <summary>
		/// Exchanges the DOF data with the other processors
		/// </summary>
		/// <param name="localList">DOF data for the current local cells</param>
		/// <param name="targetProcs">target processor for each cell (index: local cell index, value: proc rank on opComm)</param>
		/// <returns>exchanged DOF data with respect to <paramref name="targetProcs"/>.</returns>
		(long i0Cell, int lenCell)[] ExchangeDOFData((long i0Cell, int lenCell)[] localList, int[] targetProcs) {
			// Make a list for each processor
			var myList = new List<(long i0Cell, int lenCell)>();
			List<(long i0Cell, int lenCell)>[] commBuffers = opCommSize.ForLoop(i => new List<(long i0Cell, int lenCell)>());
            for (int i = 0; i < targetProcs.Length; i++) {
                int proc = targetProcs[i];
                commBuffers[proc].Add(localList[i]);
			}

			// Create send data
			Dictionary<int, (long i0Cell, int lenCell)[]> commBuffersGlobal = new Dictionary<int, (long i0Cell, int lenCell)[]>();
            for (int p = 0; p < opCommSize; p++) { //loop over other processsors
				if (p == opCommRank)
					continue;

                commBuffersGlobal.Add(p, commBuffers[p].ToArray());
				commBuffers[p].SaveToTextFileDebug($"t_SendListfrom{opCommRank}to{p}", "txt");
			}

			// Exchange the data
			IDictionary<int, (long i0Cell, int lenCell)[]> rcvData = SerialisationMessenger.ExchangeData(commBuffersGlobal);
			foreach (var data in rcvData) {
                myList.AddRange(data.Value);
				data.Value.SaveToTextFileDebug($"t_recListfrom{rcvData.Keys}to{opCommRank}", "txt");
			}

			myList.AddRange(commBuffers[opCommRank]); //add the indices that is already on the processor and was designated to be on this processor
			(long i0Cell, int lenCell)[] filteredSorted = myList.Where(x => x.lenCell > 0)
	                                                            .OrderBy(x => x.i0Cell)
	                                                            .ToArray();

			return filteredSorted;
		}

		private (BlockMsrMatrix Matrix, List<long> RowIndices, List<long> ColIndices) GetPermutationMatrix(MultigridMapping mapping, (long i0Global, int CellLen)[] globDOFsData) {
			using (new FuncTrace()) {

				BlockPartitioning TargetPartitioning = GetPartitioning(globDOFsData);
				BlockMsrMatrix PermutationMatrix = new BlockMsrMatrix(TargetPartitioning, mapping);
                List<long> RowIndices = new List<long>(); 
                List<long> ColIndices= new List<long>();
				{
					int cnt = 0;
						if (globDOFsData != null && globDOFsData.Length > 0) {
							int NoCells = globDOFsData.Length;

							for (int j = 0; j < NoCells; j++) {
								//int Len = part.lenCell[j];
								int Len = globDOFsData[j].CellLen;
								long i0Row = TargetPartitioning.i0 + cnt;
								long i0Col = globDOFsData[j].i0Global;   //part.i0Cell[j];

								PermutationMatrix.AccBlock(i0Row, i0Col, 1.0, MultidimensionalArray.CreateEye(Len));

								for (int k = 0; k < Len; k++) {
									RowIndices.Add(i0Row + k);
									ColIndices.Add(i0Col + k);
									//columnMappingOpToSmoother.Add((i0Col + k, i0Row + k));
							}

								cnt += Len;
							}
						}


					
				}

#if DEBUG
				{
						if (RowIndices != null) {

							int L = RowIndices.Count;
							for (int l = 0; l < L; l++) {
								if (l > 0) {
									Debug.Assert(RowIndices[l] > RowIndices[l - 1], "Error, Row indexing is not strictly increasing for some reason.");
									Debug.Assert(ColIndices[l] > ColIndices[l - 1], "Error, Column indexing is not strictly increasing for some reason.");
								}
							}

						}

					
				}
#endif


				return (PermutationMatrix, RowIndices, ColIndices);
			}
		}

		private int GetRank(MPI_Comm op_comm) {
			csMPI.Raw.Comm_Rank(op_comm, out int m_rank);
			return m_rank;
		}

		MPI_Comm myComm => subComm != null ? subComm : opComm;

		MPI_Comm subComm;
		int subCommRank; 
		int subCommSize; 

		MPI_Comm opComm => m_OpMapPair.OperatorMatrix.MPI_Comm;
		int opCommRank => m_OpMapPair.OperatorMatrix._RowPartitioning.MpiRank;
		int opCommSize => m_OpMapPair.OperatorMatrix._RowPartitioning.MpiSize;

		readonly int worldCommRank;

		int[] GetNoOfSpeciesList(MultigridMapping Map) {

			int J = Map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
			int[] NoOfSpecies = new int[J];

			XdgAggregationBasis xb = (XdgAggregationBasis)(Map.AggBasis.FirstOrDefault(b => b is XdgAggregationBasis));

			if (xb != null) {
				for (int jCell = 0; jCell < J; jCell++) {
					NoOfSpecies[jCell] = xb.GetNoOfSpecies(jCell);
				}
			} else {
				NoOfSpecies.SetAll(1);
			}
            
			// MPI gather on rank 0
			int MPIsz = opCommSize;
			int[] rcvCount = MPIsz.ForLoop(r => Map.AggGrid.CellPartitioning.GetLocalLength(r));
			return NoOfSpecies.MPIGatherv(rcvCount, 0, opComm);
		}

		/// <summary>
		/// Distributes the cells to new number of procs on opComm
		/// </summary>
		/// <param name="MGMapping">Mapping used for cells</param>
		/// <param name="NoOfParts">Number of target partitions</param>
		/// <param name="CoreOffset">Offset for MPI ranks</param>
		/// <returns></returns>
		(int[] partLoc, IList<(long, long)> columnMapping) DistributeMapping(MultigridMapping MGMapping, int NoOfParts, int CoreOffset) {
			using (new FuncTrace()) {
				var columnMapping = new List<(long sourceIdx, long targetIdx)>();
				int[] partLoc;
				if (NoOfParts == 1) {
					columnMapping = Enumerable.Range(0, (int)MGMapping.AggGrid.CellPartitioning.TotalLength).Select(i => ((long)i, (long)i)).ToList();
					partLoc = Enumerable.Repeat(0, MGMapping.AggGrid.CellPartitioning.LocalLength).ToArray();
					return (partLoc, columnMapping);
				}


				(int[] xadj, int[] adjncy) = GetCurrentAggGridGraphForMetis(MGMapping);
				int[] NoOfSpecies = GetNoOfSpeciesList(MGMapping);

				int MPIrnk = opCommRank;
				int[] part;
				if (MPIrnk == 0) { //call metis on only the rank 0 (opComm)
					int ncon = 1;
					int edgecut = 0;
					int[] options = new int[METIS.METIS_NOPTIONS];
					METIS.SETDEFAULTOPTIONS(options);

					options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = 1; // 
					options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 10; // This is the default refinement iterations
					options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 30; // Maximum imbalance of 3 percent (this is the default kway clustering)
					options[(int)METIS.OptionCodes.METIS_OPTION_NUMBERING] = 0;

					int J = xadj.Length - 1;
					part = new int[J];
					Debug.Assert(xadj.Where(idx => idx > adjncy.Length).Count() == 0);
					Debug.Assert(adjncy.Where(j => j >= J).Count() == 0);

					METIS.PARTGRAPHKWAY(
							ref J, ref ncon,
							xadj,
							adjncy.ToArray(),
							NoOfSpecies,
							null,
							null,
							ref NoOfParts,
							null,
							null,
							options,
							ref edgecut,
							part);
					if (CoreOffset != 0)
						part = part.Select(i => i + CoreOffset).ToArray();
				} else {
					part = null;
				}

				//broadcast to every processor. (this is necessary to create column mapping)
				var partGlob = part.MPIBroadcast(0,opComm); 

				int[] i0part = new int[NoOfParts]; //new i0partitioning for each part
				for (int p = 1; p < NoOfParts; p++)
					i0part[p] = i0part[p-1] + partGlob.Where(b => b == p-1).Count();

				// create column mapping
				for (int b = 0; b < partGlob.Length; b++) {
					int targetRank = partGlob[b];
					columnMapping.Add((b, i0part[targetRank]));
					i0part[targetRank]++; //iterate the index for the next cell
				}

				partLoc = partGlob.GetSubVector((int)MGMapping.AggGrid.CellPartitioning.i0, MGMapping.AggGrid.CellPartitioning.LocalLength);
				return (partLoc, columnMapping);
			}
		}

		(int[] xadj, int[] adj) GetCurrentAggGridGraphForMetis(MultigridMapping Map) {
			using (new FuncTrace()) {
                var aggGrid = Map.AggGrid;
				int Jloc = aggGrid.iLogicalCells.NoOfLocalUpdatedCells;
				int Jglb = checked((int)aggGrid.CellPartitioning.TotalLength);
				long[] GlidxExt = aggGrid.iParallel.GlobalIndicesExternalCells;

				var comm = opComm;
				int MPIrnk = opCommRank;
				int MPIsiz = opCommSize;

				int[][] Neighbors = aggGrid.iLogicalCells.CellNeighbours;
				int L = Neighbors.Sum(ll => ll.Length);
				int[] adj = new int[L]; // METIS input; neighbor vertices
				int[] xadj = new int[Jloc + (MPIrnk == MPIsiz - 1 ? 1 : 0)]; // METIS input: offset into `adj`
				int cnt = 0;
				int j0 = checked((int)aggGrid.CellPartitioning.i0);
				for (int j = 0; j < Jloc; j++) {
					xadj[j] = cnt;
					int[] Neighbors_j = Neighbors[j];
					for (int iN = 0; iN < Neighbors_j.Length; iN++) {
						int jNeigh = Neighbors_j[iN];
						if (jNeigh < Jloc) {
							// local cell
							adj[cnt] = jNeigh + j0;
						} else {
							adj[cnt] = checked((int)GlidxExt[jNeigh - Jloc]);
						}
						cnt++;
					}
				}
				Debug.Assert(cnt == L);
				if (MPIrnk == MPIsiz - 1)
					xadj[Jloc] = cnt;
				Debug.Assert(xadj[0] == 0);

				// convert `xadj` to global indices
				int[] locLengths = L.MPIAllGather(comm);
				Debug.Assert(locLengths.Length == MPIsiz);
				int[] xadjOffsets = new int[MPIsiz];
				for (int i = 1; i < locLengths.Length; i++) {
					xadjOffsets[i] = xadjOffsets[i - 1] + locLengths[i - 1];
				}
				for (int j = 0; j < xadj.Length; j++)
					xadj[j] += xadjOffsets[MPIrnk];

				// gather `xadj` and `adj` on Rank 0
				int[] xadjGl;
				{
					int[] rcvCount = MPIsiz.ForLoop(r => aggGrid.CellPartitioning.GetLocalLength(r));
					rcvCount[rcvCount.Length - 1] += 1; // the final length which is added on the last processor
					xadjGl = xadj.MPIGatherv(rcvCount, 0, comm);
				}

				int[] adjGl = adj.MPIGatherv(locLengths, 0, comm);
				if (MPIrnk == 0)
					Debug.Assert(xadjGl.Last() == adjGl.Length);

				// return
				return (xadjGl, adjGl);
			}
		}

		//int[] DistributeRows(BlockMsrMatrix M, int NoOfParts) {
		//	using (new FuncTrace()) {
		//		(int[] xadj, int[] adjncy) = GetCSRforRows(M);
		//		int[] NoOfSpecies = GetNoOfSpeciesList(M.);

		//		int MPIrnk = op.Mapping.MpiRank;
		//		int[] part;
		//		if (MPIrnk == 0) {
		//			int ncon = 1;
		//			int edgecut = 0;
		//			int[] options = new int[METIS.METIS_NOPTIONS];
		//			METIS.SETDEFAULTOPTIONS(options);

		//			options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = 1; // 
		//			options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 10; // This is the default refinement iterations
		//			options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 30; // Maximum imbalance of 3 percent (this is the default kway clustering)
		//			options[(int)METIS.OptionCodes.METIS_OPTION_NUMBERING] = 0;

		//			int J = xadj.Length - 1;
		//			part = new int[J];
		//			Debug.Assert(xadj.Where(idx => idx > adjncy.Length).Count() == 0);
		//			Debug.Assert(adjncy.Where(j => j >= J).Count() == 0);
	
		//			METIS.PARTGRAPHKWAY(
		//					ref J, ref ncon,
		//					xadj,
		//					adjncy.ToArray(),
		//					NoOfSpecies,
		//					null,
		//					null,
		//					ref NoOfParts,
		//					null,
		//					null,
		//					options,
		//					ref edgecut,
		//					part);
		//		} else {
		//			part = null;
		//		}

		//		// scatter back to processors
		//		int MPIsz = op.Mapping.MpiSize;
		//		int[] sendCount = MPIsz.ForLoop(r => op.GridData.CellPartitioning.GetLocalLength(r));
		//		var partLoc = part.MPIScatterv(sendCount, 0, op.OperatorMatrix.MPI_Comm);
		//		return partLoc;
		//	}
		//}

        (int[] xadj, int[] adjncy) GetCSRforRows(BlockMsrMatrix M) {
			int n = (int)M.NoOfRows; // or M.N if defined that way
			int[] xadj = new int[n + 1];
			List<int> adjncyList = new List<int>();

			xadj[0] = 0;
			for (int i = 0; i < n; i++) {
                // Assuming M is stored in a row-wise structure.
                // Get the column indices of nonzero entries for row i.
                // Exclude the diagonal if desired.
                long[] NonZeroColumns = null;
                M.GetOccupiedColumnIndices(i, ref NonZeroColumns);
				foreach (long col in NonZeroColumns) {
					if (col != i) // optionally exclude the diagonal
					{
						adjncyList.Add((int)col);
					}
				}
				xadj[i + 1] = adjncyList.Count;
			}

			int[] adjncy = adjncyList.ToArray();

            return (xadj, adjncy);
		}

		/// <summary>
		/// Creates a new communicator for the coarse level solver.
		/// </summary>
		/// <param name="NoOfCoarseProcs"></param>
		void SplitCommunicator(int NumberOfCoarseProcs) {
			//The zero-th rank is the one assigned for fine mesh for the ease of operations
			NoOfCoarseProcs = Math.Max(1, NumberOfCoarseProcs);
			NoOfSmootherProcs = opCommSize - NoOfCoarseProcs;

			if (NoOfCoarseProcs < 1 || NoOfSmootherProcs < 1)
                throw new ArgumentOutOfRangeException($"Number of coarse processors and smoother processors must be greater than 0. Coarse: {NoOfCoarseProcs}, Smoother: {NoOfSmootherProcs}");

			if (opCommRank < NoOfSmootherProcs) {
				myTask = parallelTask.Smoother;
				csMPI.Raw.CommSplit(opComm, 0, opCommRank, out subComm);
			} else {
				myTask = parallelTask.Coarse;
				csMPI.Raw.CommSplit(opComm, 1, opCommRank, out subComm);
			}

			csMPI.Raw.Comm_Rank(subComm, out subCommRank);
			csMPI.Raw.Comm_Size(subComm, out subCommSize);

            if(verbose)
			    Console.WriteLine($"The proc with worldRank-{worldCommRank} with opRank{opCommRank} is assigned to {subCommRank}of{subCommSize} new com{subComm.m1}");

            CreateMpiRankMappings();
		}

        void CreateMpiRankMappings() {
            // Each process sends its rank as its value.
            int totalCount = opCommSize; //  number of processes in the communicator

			// Prepare managed buffers.
			int[] sendBuffer = new int[] { subCommRank }; //workaround for fixed statement
			CommToSubCommMapping = new int[totalCount];

			// Pin the buffers to obtain their addresses.
			GCHandle sendHandle = GCHandle.Alloc(sendBuffer, GCHandleType.Pinned);
			GCHandle recvHandle = GCHandle.Alloc(CommToSubCommMapping, GCHandleType.Pinned);

			IntPtr sendPtr = sendHandle.AddrOfPinnedObject();
			IntPtr recvPtr = recvHandle.AddrOfPinnedObject();

			// Perform the Allgather call.
			csMPI.Raw.Allgather(sendPtr, 1, csMPI.Raw._DATATYPE.INT, recvPtr, 1, csMPI.Raw._DATATYPE.INT, opComm);

			// Display gathered values on each process.
			Console.WriteLine($"Process {opCommRank} gathered:");
			for (int i = 0; i < totalCount; i++) {
				Console.WriteLine($"Rank {i}: {CommToSubCommMapping[i]}");
			}
            WorldToSubCommMapping = CommToSubCommMapping;// ??
		}

        bool verbose = true;

        int[] CommToSubCommMapping;

		int[] WorldToSubCommMapping;

        parallelTask myTask = parallelTask.All;

		(int[] xadj, int[] adj) GetGridGraphForMetis(MultigridOperator op) {
			using (new FuncTrace()) {
				int Jloc = op.GridData.iLogicalCells.NoOfLocalUpdatedCells;
				int Jglb = checked((int)op.GridData.CellPartitioning.TotalLength);
				long[] GlidxExt = op.GridData.iParallel.GlobalIndicesExternalCells;

				var comm = op.OperatorMatrix.MPI_Comm;
				int MPIrnk = op.OperatorMatrix.RowPartitioning.MpiRank;
				int MPIsiz = op.OperatorMatrix.RowPartitioning.MpiSize;

				int[][] Neighbors = op.GridData.iLogicalCells.CellNeighbours;
				int L = Neighbors.Sum(ll => ll.Length);
				int[] adj = new int[L]; // METIS input; neighbor vertices
				int[] xadj = new int[Jloc + (MPIrnk == MPIsiz - 1 ? 1 : 0)]; // METIS input: offset into `adj`
				int cnt = 0;
				int j0 = checked((int)op.GridData.CellPartitioning.i0);
				for (int j = 0; j < Jloc; j++) {
					xadj[j] = cnt;
					int[] Neighbors_j = Neighbors[j];
					for (int iN = 0; iN < Neighbors_j.Length; iN++) {
						int jNeigh = Neighbors_j[iN];
						if (jNeigh < Jloc) {
							// local cell
							adj[cnt] = jNeigh + j0;
						} else {
							adj[cnt] = checked((int)GlidxExt[jNeigh - Jloc]);
						}
						cnt++;
					}
				}
				Debug.Assert(cnt == L);
				if (MPIrnk == MPIsiz - 1)
					xadj[Jloc] = cnt;
				Debug.Assert(xadj[0] == 0);

				// convert `xadj` to global indices
				int[] locLengths = L.MPIAllGather(comm);
				Debug.Assert(locLengths.Length == MPIsiz);
				int[] xadjOffsets = new int[MPIsiz];
				for (int i = 1; i < locLengths.Length; i++) {
					xadjOffsets[i] = xadjOffsets[i - 1] + locLengths[i - 1];
				}
				for (int j = 0; j < xadj.Length; j++)
					xadj[j] += xadjOffsets[MPIrnk];

				// gather `xadj` and `adj` on Rank 0
				int[] xadjGl;
				{
					int[] rcvCount = MPIsiz.ForLoop(r => op.GridData.CellPartitioning.GetLocalLength(r));
					rcvCount[rcvCount.Length - 1] += 1; // the final length which is added on the last processor
					xadjGl = xadj.MPIGatherv(rcvCount, 0, comm);
				}

				int[] adjGl = adj.MPIGatherv(locLengths, 0, comm);
				if (MPIrnk == 0)
					Debug.Assert(xadjGl.Last() == adjGl.Length);

				// return
				return (xadjGl, adjGl);
			}
		}

		BlockPartitioning GetPartitioning((long i0Global, int CellLen)[] DOFs) {
			using (new FuncTrace()) {
				var LnCell = new List<int>();
				var i0Cell = new List<long>();

				int cnt = 0;

				if (DOFs != null && DOFs.Length > 0) {
					int NoCells = DOFs.Length;// part.jGlobal.Length;
					for (int j = 0; j < NoCells; j++) {
						int Len = DOFs[j].CellLen; // part.lenCell[j];

						i0Cell.Add(cnt);
						LnCell.Add(Len);
						cnt += Len;
					}
				}

				int LocalLength = cnt;
				var partitioning = new BlockPartitioning(LocalLength, i0Cell, LnCell, opComm, i0isLocal: true);
				return partitioning;
			}
		}

		private void InitCoarse() {
            using (var tr = new FuncTrace()) {
                var op = this.m_OpMapPair;
                if (this.CoarserLevelSolver == null) {
                    //throw new NotSupportedException("Missing coarse level solver.");
                    tr.Info("OrthonormalizationMultigrid: running without coarse solver.");
                } else {
                    if (op is MultigridOperator mgOp) { // the first call
                        if (myConfig.CoarseOnLowerLevel && mgOp.CoarserLevel != null) {
                            this.CoarserLevelSolver.Init(mgOp.CoarserLevel);
                        } else {
                            tr.Info("OrthonormalizationMultigrid: running coarse solver on same level.");
                            this.CoarserLevelSolver.Init(mgOp);
                        }
                    } else {
                        if (myConfig.CoarseOnLowerLevel == false && this.CoarserLevelSolver is ISubsystemSolver ssCoarse) {
                            ssCoarse.Init(op);
                        } else {
                            throw new NotSupportedException($"Unable to initialize coarse-level-solver if operator is not a {typeof(MultigridOperator)}");
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
                //Res.SetV(partGlob);
                Array.Copy(B, Res, L);
                OpMatrix.SpMV(-1.0, X, 1.0, Res);
                //Res.AccV(1.0, partGlob);


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

        enum parallelTask { 
            All,
            Smoother,
            Coarse
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
			//Debugger.Launch();
			using (var f = new FuncTrace()) {
                f.StdoutOnAllRanks();

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
                if (this.CoarserLevelSolver != null && myConfig.CoarseOnLowerLevel) {
                    Lc = ((MultigridOperator)m_OpMapPair).CoarserLevel.Mapping.LocalLength;
                    ResCoarse = new double[Lc];
                } else {
                    Lc = -1;
                    ResCoarse = null;
                }

				parallelTask myTask = parallelTask.All;
				MPI_Comm newComm = m_OpMapPair.DgMapping.MPI_Comm;

				//if (PreSmoother is SchwarzForCoarseMesh schwarz) {
				//	Console.WriteLine($"Tot number of Schwarz blocks {schwarz.config.NoOfBlocks}");
				//	if (schwarz.config.NoOfBlocks < m_OpMapPair.DgMapping.MpiSize) {
				//		Console.WriteLine("A possible speed up");
				//		Console.WriteLine($"Rank-{m_OpMapPair.DgMapping.MpiRank} old com{m_OpMapPair.DgMapping.MPI_Comm.m1}");
    //                    schwarz.CommitInitialXandB(X, partGlob);
                        
    //                    if (m_OpMapPair.DgMapping.MpiRank < schwarz.config.NoOfBlocks) {
    //                        myTask = parallelTask.Smoother;
    //                        csMPI.Raw.CommSplit(m_OpMapPair.DgMapping.MPI_Comm, 0, m_OpMapPair.DgMapping.MpiRank, out newComm);
    //                    } else {
    //                        myTask = parallelTask.Coarse;
    //                        csMPI.Raw.CommSplit(m_OpMapPair.DgMapping.MPI_Comm, 1, m_OpMapPair.DgMapping.MpiRank, out newComm);
    //                    }
    //                    subComm = newComm;
				//		var newOP = ReInitImpl(m_OpMapPair, newComm);
				//		Console.WriteLine($"Rank-{m_OpMapPair.DgMapping.MpiRank} new com{newComm.m1}");

    //                }
				//}

				csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

				bool done = false;
				int SmootherGroupLeader = 0;
				int CoarseGroupLeader = m_OpMapPair.DgMapping.MpiSize - 1;
                done.MPIAnd(newComm);
				double[] Sol0 = X.CloneAs();
                double[] Res0 = new double[L];
                Residual(Res0, Sol0, B);
                Array.Copy(Res0, Res, L);

                if(ortho.Norm(Res0) <= 0) {
                    double normB = ortho.Norm(B);
                    double normX = ortho.Norm(X);
                    f.Error($" **** Residual is 0.0: |X| = {normB}; |partGlob| = {normB}; |Res0| = {ortho.Norm(Res0)}");


                }

                int iLevel = ((m_OpMapPair as MultigridOperator)?.LevelIndex ?? -1);

                void WriteDebug(int iter, double res, string text) {
                    if (iLevel >= 0)
					    Console.WriteLine($"{string.Concat(Enumerable.Repeat("-", iLevel))} OrthoMG, current level={iLevel}, iteration={iter} {(text != null ? " - " + text : "")} and res norm: {res}");
					
                    return;
                }

                double iter0_resNorm = ortho.Norm(Res0);
                double resNorm = iter0_resNorm;
                this.IterationCallback?.Invoke(0, Sol0, Res0, this.m_OpMapPair as MultigridOperator);
                                
                // clear history of coarse solvers
                ortho.Clear();
                bool bIterate = true;



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
				WriteDebug(0, resNorm, $"NonSerialPreSmoother for iterative solver is {(config.NonSerialPreSmoother && !config.SkipPreSmoother ? "activated" : "deactivated")}");
				
                

				int iIter;
                for (iIter = 1; bIterate; iIter++) {
                    WriteDebug(iIter, resNorm, "initial start");
					csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

					var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if (!termState.bNotTerminate) {
                        Converged = termState.bSuccess;
                        break;
                    } else {

                    }


					// pre-smoother
					// ------------
					if (myTask == parallelTask.All || myTask == parallelTask.Smoother) {
						VerivyCurrentResidual(X, B, Res, iIter);
						double[] PreCorr = new double[L];

						if (m_OpMapPair.DgMapping.MpiRank == SmootherGroupLeader) {
							byte completionSignal = 0;
							MPI_Request req;
							unsafe {
								byte* pSignal = &completionSignal;
								IntPtr ptr = (IntPtr)pSignal;
								csMPI.Raw.Irecv(ptr, 1,
												csMPI.Raw._DATATYPE.BYTE,
												CoarseGroupLeader, 0,
												m_OpMapPair.DgMapping.MPI_Comm, out req);

							}
							int k = 1;
							while (!done) {
								PreSmoother.Solve(PreCorr, Res); // Vorglättung
								Thread.Sleep(1000);
								Console.WriteLine($"{k}-second waiting for completionSignal={completionSignal}");
								csMPI.Raw.Test(ref req, out done, out MPI_Status status);
								done.MPIOr(newComm);
								k++;
							}
							Console.WriteLine($"waited {k}-second for completionSignal={completionSignal}");

						} else {
							while (!done) {
								PreSmoother.Solve(PreCorr, Res); // Vorglättung
								Thread.Sleep(1000);
								done.MPIOr(newComm);
							}

						}

                        if (PreSmoother is SchwarzForCoarseMesh schwarz2) { 
							schwarz2.PropagateBackToOpComm(PreCorr);
						}
						// orthonormalization and residual minimization
						resNorm = ortho.AddSolAndMinimizeResidual(ref PreCorr, X, Sol0, Res0, Res, "presmoothL" + iLevel);

						WriteDebug(iIter, resNorm, " pre-smoother applied");

						//SpecAnalysisSample(iIter, X, "ortho1");
						var termState2 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
						if (!termState2.bNotTerminate) {
							Converged = termState2.bSuccess;
							break;
						}

					}

					

				

                    // coarse grid correction
                    // ----------------------
                    if (myTask == parallelTask.All || myTask == parallelTask.Coarse) {

						if (m_OpMapPair.DgMapping.MpiRank == CoarseGroupLeader) {
							Console.WriteLine("Waiting for 20 seconds...");
							Thread.Sleep(20000);  // 20 seconds delay
							Console.WriteLine("Done waiting.");

							//while (!done) {
							//	Thread.Sleep(1000);
							//	done.MPIOr(newComm);
							//}

							byte completionSignal = 1;
							unsafe {
								byte* pSignal = &completionSignal;
								IntPtr ptr = (IntPtr)pSignal;
								csMPI.Raw.Send(ptr, 1, csMPI.Raw._DATATYPE.BYTE, SmootherGroupLeader, 0, m_OpMapPair.DgMapping.MPI_Comm);
							}
							Console.WriteLine("Sent signal");
						} else {
							Thread.Sleep(19000);
							done = true;

						}
						Thread.Sleep(2000000);

						CrseLevelTime.Start();
                        // Test: Residual on this level / already computed by 'MinimizeResidual' above
                        VerivyCurrentResidual(X, B, Res, iIter);

                        double resNorm_b4Coarse = resNorm;
                        for (int i = 0; i < myConfig.m_omega; i++) {
                            if (this.CoarserLevelSolver != null) {

                                double[] vl = new double[L];
                                if (myConfig.CoarseOnLowerLevel) {

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
                                if (vl.ContainsNanOrInf().MPIOr() && CoarseArithmeticExceptionCount < 20) {
                                    CoarseArithmeticExceptionCount++;
                                    Console.WriteLine("Coarse solver failed " + CoarseArithmeticExceptionCount);
                                    if (CoarseArithmeticExceptionCount == 1) {
                                        BlockMsrMatrix coarseMtx = myConfig.CoarseOnLowerLevel ? (m_OpMapPair as MultigridOperator).CoarserLevel.OperatorMatrix : m_OpMapPair.OperatorMatrix;
                                        coarseMtx.SaveToTextFileSparse("FailedCoarseMatrix.txt");
                                    }
                                    vl = null;
                                    this.CoarserLevelSolver.Dispose();
                                    this.InitCoarse();
                                } else {
                                    resNorm = ortho.AddSolAndMinimizeResidual(ref vl, X, Sol0, Res0, Res, "coarsecorL" + iLevel);
                                }

                            }
                        } // end of coarse-solver loop
                        WriteDebug(iIter, resNorm, "coarse-solver applied");

                        var termState3 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                        if (!termState3.bNotTerminate) {
                            Converged = termState3.bSuccess;
                            break;
                        }
                        CrseLevelTime.Stop();                       
					}

                    if (myTask == parallelTask.All || myTask == parallelTask.Smoother) {

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

                                bool fail = false;
                                try {
                                    _PostSmoother.Solve(PostCorr, Res); // compute correction (Nachglättung)
                                } catch (ArithmeticException ae) {
                                    f.Error($"Smoother fail on Rank {m_OpMapPair.DgMapping.MpiRank}: " + ae.ToString());
                                    fail = true;
                                }
                                fail = fail.MPIOr();

                                if ((PostCorr.ContainsNanOrInf() || fail).MPIOr()) {
                                    PostSmootherArithmeticExceptionCount++;
                                    Console.Error.WriteLine("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g + " for " + PostSmootherArithmeticExceptionCount + "-th time.");

                                    if (Res.ContainsNanOrInf()) {
                                        Console.WriteLine("... so does RHS");
                                        throw new ArithmeticException("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g + " RHS already corrupted.");
                                    } else
                                        Console.WriteLine("... although RHS is regular");

                                    //var viz = new MGViz(m_OpMapPair as MultigridOperator);
                                    //var __RHS = viz.ProlongateRhsToDg(Res, "RES");
                                    //var __SOL = viz.ProlongateRhsToDg(PostCorr, "SOL");
                                    //Tecplot.Tecplot.PlotFields(__SOL.Cat(__RHS), "iilufail", 0.0, 0);
                                    //

                                    if (PostSmootherArithmeticExceptionCount < 20) {
                                        PostCorr = null;
                                        _PostSmoother.Dispose();
                                        if (PostSmoother is ISubsystemSolver ssPostSmother) {
                                            ssPostSmother.Init(this.m_OpMapPair);
                                        } else {
                                            if (this.m_OpMapPair is MultigridOperator mgOp) {
                                                PostSmoother.Init(mgOp);
                                            } else {
                                                throw new NotSupportedException($"Unable to initialize post-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
                                            }
                                        }
                                    } else {
                                        throw new ArithmeticException("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g);
                                    }

                                } else {
                                    resNorm = ortho.AddSolAndMinimizeResidual(ref PostCorr, X, Sol0, Res0, Res, "pstsmthL" + iLevel + "-sw" + g);
                                    WriteDebug(iIter, resNorm, "post-smoother applied");



                                    var termState4 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                                    if (!termState4.bNotTerminate) {
                                        Converged = termState4.bSuccess;
                                        termPost = true;
                                        break;
                                    }

                                    if (ortho.CancellationTriggered) {
                                        //skipPreSmoothInFollowingIters = true; // most of the time, the pre-smoother does nothing different than the post-smoother; So, if we cancel post-smoothing there is no need to do pre-smoothing in the next loop.

                                        iPostSmooter++;
                                        if (iPostSmooter >= allSmooters.Length)
                                            iPostSmooter = 0;




                                        break;
                                    }
                                }
                            }


                            if (termPost)
                                break;

                        } // end of post-smoother loop
                    }

                    
                    // iteration callback
                    // ------------------
                    this.ThisLevelIterations++;
                    IterationCallback?.Invoke(iIter, X, Res, this.m_OpMapPair as MultigridOperator);

                } // end of solver iterations

                IterationCallback?.Invoke(iIter, X, Res, this.m_OpMapPair as MultigridOperator);
				WriteDebug(iIter, resNorm, "final");


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
                    Residual(rTest, X, partGlob); // Residual on this level; 
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

        int PostSmootherArithmeticExceptionCount = 0;
        int CoarseArithmeticExceptionCount = 0;


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

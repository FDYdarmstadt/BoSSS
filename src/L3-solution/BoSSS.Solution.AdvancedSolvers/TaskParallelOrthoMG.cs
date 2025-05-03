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
using System.Reflection.Emit;
using System.Drawing;
using NUnit.Framework.Internal;
using BoSSS.Foundation.IO;
using System.ComponentModel.Design;

namespace BoSSS.Solution.AdvancedSolvers {

	enum TpTaskType {
		All,
		Smoother,
		Coarse
	}

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

	// this is designed to hold information getting from the WORLD
	public class TaskParallelMGOperator : IOperatorMappingPair {

		/// <summary>
		/// An ad-hoc information storing class for World-Level MG operator
		/// Constructor (responsible for world to sub communicators)
		/// should be defined on world communicator and then will be used on sub communicators with the solver
		/// Basically this defines the operator matrix and the prolongation matrix and transfers them to the responsible processors without changing the communicator
		/// Technically each level is responsible for the prolongation matrix to the finer level and the operator matrix of the current level
		/// Every level works with three different communicators: level from finer level (coarse comm of the finer), 
		/// and the smoother and the coarse for this level (still on World but the processors are designated)
		/// In each level first processors are are for the smoother and the rest is for the coarse,
		/// leading that 0-th rank is always the smoother with the finest and last rank is for the coarse on the coarsest level
		/// This can be called from the main program or from a finer level of another solver, so the finest level might have or not a prolongation op.
		/// </summary>
		/// <param name="OperatorMatrix"></param>
		/// <param name="RestrictionMatrix"></param>
		/// <param name="ProlongationMatrix"></param>
		/// <param name="Mapping"></param>
		/// <param name="size"></param>
		public TaskParallelMGOperator(BlockMsrMatrix OperatorMatrix, BlockMsrMatrix ProlongationMatrix, MultigridMapping Mapping, int size = 0, IOperatorMappingPair finerLevel = null) { 
            m_OpMtx = OperatorMatrix;
			m_ProlMtx = ProlongationMatrix;
			m_MultigridMapping = Mapping;
			NoOfThisProcs = size == 0 ? m_OpMtx.RowPartitioning.MpiSize : size;
			FinerLevel = finerLevel;
			CheckMPICorrectness();

			NoOfSmootherProcs = Math.Max(1, NoOfThisProcs / 4);
			(m_xadj,m_adj) = GetCurrentAggGridGraphForMetis(m_MultigridMapping);
            m_NoOfSpecies = GetNoOfSpeciesList(m_MultigridMapping);
			(ThisCellI0s, ThisNewCellMapping) = DistributeMapping(m_MultigridMapping, NoOfThisProcs);
			(SmootherCellI0s, SmootherNewCellMapping) = DistributeMapping(m_MultigridMapping, NoOfSmootherProcs);
			(CoarseCellI0s, CoarseNewCellMapping) = DistributeMapping(m_MultigridMapping, NoOfCoarseProcs);
            CellToDOFdata = CellIndexToDOFs(m_MultigridMapping);
			CalculateWorldToSubDistribution();
			//SplitCommunicator();

			// To Do:
			// 1. Check if the change of basis is correct
			// 2. Check if the prolongation matrix is correct
			// 3. Check if the operator matrix is correct (write a test function)
		}

		public IBlockPartitioning DefaultPartition => m_MultigridMapping;
		public IBlockPartitioning FinerLevelCoarsePartitiong => FinerLevel is TaskParallelMGOperator fine ? fine.CoarseTargetPartitioning : DefaultPartition;
		public IList<(long Source, long Target)> FinerLevelCoarseCellMapping => FinerLevel is TaskParallelMGOperator fine ? fine.CoarseNewCellMapping : null;

		public IBlockPartitioning ThisTargetPartitioning;
		public IBlockPartitioning CoarseTargetPartitioning;
		public IBlockPartitioning SmootherTargetPartitioning;

		void CheckMPICorrectness() {
			Console.WriteLine($"TaskParallelMGOperator: {m_MultigridMapping.MPI_Comm.m1} rank {m_MultigridMapping.MpiRank} of {m_MultigridMapping.MpiSize} with {m_MultigridMapping.MpiRank} of {m_MultigridMapping.MpiSize} and {NoOfThisProcs} size and {myTask}");
			if (m_MultigridMapping.MPI_Comm != csMPI.Raw._COMM.WORLD)
				throw new ArgumentException("The mg mapping must be defined on the WORLD communicator.");

			//if (!thisLevelWorldRankings.Contains(worldCommRank))
			//	//throw new ArgumentOutOfRangeException($"The world rank {worldCommRank} is not in the range of this level. The range is {thisLevelWorldRankings[0]}-{thisLevelWorldRankings[thisLevelWorldRankings.Length - 1]}." +
			//		$" The smoothers are {smoothLevelWorldRankings[0]}-{smoothLevelWorldRankings[smoothLevelWorldRankings.Length - 1]} and the coarse are {coarseLevelWorldRankings[0]}-{coarseLevelWorldRankings[coarseLevelWorldRankings.Length - 1]}.");
		}


		int[] thisLevelWorldRankings => Enumerable.Range(worldMPIOffset, NoOfThisProcs).ToArray();
		int[] smoothLevelWorldRankings => Enumerable.Range(worldMPIOffset, NoOfSmootherProcs).ToArray();
		int[] coarseLevelWorldRankings => Enumerable.Range(worldMPIOffset + NoOfSmootherProcs, NoOfCoarseProcs).ToArray();
		//int[] ThisCommToSubCommMapping;

		public MPI_Comm currentComm => m_MultigridMapping.MPI_Comm;

		public int worldMPIOffset => worldCommSize - NoOfThisProcs; //the offset of this operator matrix in the world communicator, the mpi rank of the first processor in this level

		public int worldCommRank => m_MultigridMapping.MpiRank;
		public int worldCommSize => m_MultigridMapping.MpiSize;
		//public MPI_Comm thisComm => FinerLevel == null ? m_MultigridMapping.MPI_Comm : FinerLevel is TaskParallelMGOperator fine ? fine.subComm : FinerLevel.DgMapping.MPI_Comm;
		//public int thisCommRank => FinerLevel is TaskParallelMGOperator fine ? fine.subCommRank : FinerLevel.DgMapping.MpiRank;
		//public MPI_Comm subComm;
		//int subCommRank;
		//int subCommSize;


		TpTaskType myTask = TpTaskType.All; // this is used to determine the task of the current processor (smoother or coarse) in the sub communicator	

		//bool verbose = true;

		internal (long i0Cell, int lenCell)[] localBlocksForThisLevel;
		internal (long i0Cell, int lenCell)[] localBlocksForSmoother;
		internal (long i0Cell, int lenCell)[] localBlocksForCoarse;

		void CalculateWorldToSubDistribution() {
			localBlocksForThisLevel = GetLocalDistribution(CellToDOFdata, ThisNewCellMapping, ThisCellI0s, worldMPIOffset, NoOfThisProcs);
			localBlocksForSmoother = GetLocalDistribution(CellToDOFdata, SmootherNewCellMapping, SmootherCellI0s, worldMPIOffset, NoOfSmootherProcs);
			localBlocksForCoarse = GetLocalDistribution(CellToDOFdata, CoarseNewCellMapping, CoarseCellI0s, worldMPIOffset + NoOfSmootherProcs, NoOfCoarseProcs);

			ThisTargetPartitioning = GetPartitioning(localBlocksForThisLevel, currentComm);
			SmootherTargetPartitioning = GetPartitioning(localBlocksForSmoother, currentComm);
			CoarseTargetPartitioning = GetPartitioning(localBlocksForCoarse, currentComm);
			
			m_OpMtx = ChangeOpPartitioning(ThisTargetPartitioning, ThisNewCellMapping, "Op"); //this level of 
			m_OpMtx_smoother = ChangeOpPartitioning(SmootherTargetPartitioning, SmootherNewCellMapping, "OpS"); //this level of solving is only for smoother so op should partitioned for smoother (reserved for optimization)
			m_ProlMtx = ChangeProlPartitioning(ThisTargetPartitioning, ThisNewCellMapping, FinerLevelCoarsePartitiong, FinerLevelCoarseCellMapping, "Pro"); //same processors differernt level of mg operator (different number of cells for row and column)
		}

		BlockMsrMatrix ChangeOpPartitioning(IBlockPartitioning targetPartitioning, IList<(long Source, long Target)> newBlockIndices = null, string tag = "Op") {
			if (verbose) { 
			m_OpMtx.SaveToTextFileSparseDebug($"lvl{level}_{tag}_oldOp.txt");
			m_OpMtx.SaveToTextFileSparse($"lvl{level}_{tag}_oldOp.txt");
			}
			BlockMsrMatrix ret;

			if (m_OpMtx.RowPartitioning == targetPartitioning && m_OpMtx.ColPartition == targetPartitioning && newBlockIndices == null) //this can happen at the first level of tp
				ret = m_OpMtx.CloneAs();
			else
				ret = m_OpMtx.CloneAs().ChangePartitioning(targetPartitioning, targetPartitioning, newBlockIndices, newBlockIndices);

			if (verbose) {
				ret.SaveToTextFileSparseDebug($"lvl{level}_{tag}_newOp.txt");
				ret.SaveToTextFileSparse($"lvl{level}_{tag}_newOp.txt");
			}

			return ret;
		}

		BlockMsrMatrix ChangeProlPartitioning(IBlockPartitioning thisPartitioning, IList<(long Source, long Target)> thisNewBlockIndices, IBlockPartitioning finerPartitioning, IList<(long Source, long Target)> finerNewBlockIndices, string tag = "c") {
			if (m_ProlMtx == null)
				return null;

			if (verbose) {
				m_ProlMtx.SaveToTextFileSparse($"lvl{level}_{tag}_oldPro.txt");
				m_ProlMtx.SaveToTextFileSparseDebug($"lvl{level}_{tag}_oldPro.txt");
			}
			var ret = m_ProlMtx.ChangePartitioning(finerPartitioning, thisPartitioning, finerNewBlockIndices, thisNewBlockIndices);

			if (verbose) { 
			ret.SaveToTextFileSparse($"lvl{level}_{tag}_newPro.txt");
			ret.SaveToTextFileSparseDebug($"lvl{level}_{tag}_newPro.txt");
			}

			return ret;
		}

		bool verbose = true;

		BlockPartitioning GetPartitioning((long i0Global, int CellLen)[] DOFs, MPI_Comm comm) {
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
				var partitioning = new BlockPartitioning(LocalLength, i0Cell, LnCell, comm, i0isLocal: true);
				return partitioning;
			}
		}

		(long i0Cell, int lenCell)[] GetLocalDistribution((long i0Cell, int lenCell)[] globalDOFs, List<(long Source, long Target)> cellColumnMapping, long[] targeti0s, int procOffset, int procSize) {
			int newRank = worldCommRank - procOffset;
			if (newRank < 0 || newRank >= procSize)
				return new (long i0Cell, int lenCell)[0];

			var ret = new List<(long i0Cell, int lenCell)>();
			var newBlocki0 = (int)targeti0s[newRank];
			var newBlockiE = (int)targeti0s[newRank + 1];

			// if cell is designated to be on this processor, add it to the list
			for (long iCell = 0; iCell < cellColumnMapping.Count; iCell++) {
				var iTargetCell = cellColumnMapping[(int)iCell].Target; //designated cell index
				var iSourceCell = cellColumnMapping[(int)iCell].Source;

				if (iTargetCell >= newBlocki0 && iTargetCell < newBlockiE) { // if designated cell index falls into the range of this processor
					var sourceBlock = globalDOFs[iSourceCell]; // get the current block index (source)
					ret.Add(sourceBlock);
				}
			}
			ret.Sort((x, y) => x.i0Cell.CompareTo(y.i0Cell)); // sort the list according to the cell index
			Debug.Assert(newBlockiE - newBlocki0 == ret.Count);
			return ret.ToArray();
		}

		public int level => FinerLevel is TaskParallelMGOperator fine ? fine.level + 1 : 0;
        public int NoOfThisProcs; 
        public int NoOfSmootherProcs;
        public int NoOfCoarseProcs => NoOfThisProcs - NoOfSmootherProcs;

		public (long i0Cell, int lenCell)[] CellToDOFdata; //global data (all of them) and global indices
		public List<(long source, long target)> ThisNewCellMapping; //global data (all of them) and global indices
		public List<(long source, long target)> CoarseNewCellMapping; //global data (all of them) and global indices
		public List<(long source, long target)> SmootherNewCellMapping; //global data (all of them) and global indices

		public long[] ThisCellI0s;
		public long[] CoarseCellI0s;
        public long[] SmootherCellI0s;

		public int[] m_xadj;
		public int[] m_adj;
        public int[] m_NoOfSpecies;

		//public long celli0;
		//public int localCellLength;
		//public int totalCellLength;

		public IOperatorMappingPair FinerLevel;
		public TaskParallelMGOperator CoarserLevel;

		public BlockMsrMatrix m_OpMtx;
		public BlockMsrMatrix m_OpMtx_smoother;
		public BlockMsrMatrix m_ProlMtx;
		public MultigridMapping m_MultigridMapping;

		//public BlockMsrMatrix old_OpMtx;
		//public BlockMsrMatrix old_ProlMtx;

		public BlockMsrMatrix OperatorMatrix => m_OpMtx;
		public BlockMsrMatrix ProlongationMatrix => m_ProlMtx;
		public ICoordinateMapping DgMapping => m_MultigridMapping;

		(int[] xadj, int[] adj) GetCurrentAggGridGraphForMetis(MultigridMapping Map) {
			using (new FuncTrace()) {
				var aggGrid = Map.AggGrid;
				int Jloc = aggGrid.iLogicalCells.NoOfLocalUpdatedCells;
				int Jglb = checked((int)aggGrid.CellPartitioning.TotalLength);
				long[] GlidxExt = aggGrid.iParallel.GlobalIndicesExternalCells;

				var comm = Map.MPI_Comm;
				int MPIrnk = Map.MpiRank;
				int MPIsiz = Map.MpiSize;

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

				// convert `m_xadj` to global indices
				int[] locLengths = L.MPIAllGather(comm);
				Debug.Assert(locLengths.Length == MPIsiz);
				int[] xadjOffsets = new int[MPIsiz];
				for (int i = 1; i < locLengths.Length; i++) {
					xadjOffsets[i] = xadjOffsets[i - 1] + locLengths[i - 1];
				}
				for (int j = 0; j < xadj.Length; j++)
					xadj[j] += xadjOffsets[MPIrnk];

				// gather `m_xadj` and `adj` on Rank 0
				int[] xadjGl;
				{
					int[] rcvCount = MPIsiz.ForLoop(r => aggGrid.CellPartitioning.GetLocalLength(r));
					rcvCount[rcvCount.Length - 1] += 1; // the final length which is added on the last processor
					xadjGl = xadj.MPIGatherv(rcvCount, worldMPIOffset, comm); // store at the first responsible processor of this level
				}

				int[] adjGl = adj.MPIGatherv(locLengths, worldMPIOffset, comm); // store at the first responsible processor of this level
				if (MPIrnk == worldMPIOffset)
					Debug.Assert(xadjGl.Last() == adjGl.Length);

				// return
				return (xadjGl, adjGl);
			}
		}

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
			int MPIsz = Map.MpiSize;
			int[] rcvCount = MPIsz.ForLoop(r => Map.AggGrid.CellPartitioning.GetLocalLength(r));
			return NoOfSpecies.MPIGatherv(rcvCount, worldMPIOffset, Map.MPI_Comm);
		}

		public (long[] cellI0s, List<(long source, long target)> colMapping) DistributeMapping(MultigridMapping MGMapping, int NoOfParts) {
			using (new FuncTrace()) {
				var columnMapping = new List<(long sourceIdx, long targetIdx)>();
				if (NoOfParts == 1) {
					columnMapping = Enumerable.Range(0, (int)MGMapping.AggGrid.CellPartitioning.TotalLength).Select(i => ((long)i, (long)i)).ToList();
					return (new long[] { 0, (int)MGMapping.AggGrid.CellPartitioning.TotalLength }, columnMapping);
				}

				int MPIrnk = MGMapping.MpiRank;
				int[] part;
				if (MPIrnk == worldMPIOffset) { //call metis on only the rank 0 (opComm)
					int ncon = 1;
					int edgecut = 0;
					int[] options = new int[METIS.METIS_NOPTIONS];
					METIS.SETDEFAULTOPTIONS(options);

					options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = 1; // 
					options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 10; // This is the default refinement iterations
					options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 30; // Maximum imbalance of 3 percent (this is the default kway clustering)
					options[(int)METIS.OptionCodes.METIS_OPTION_NUMBERING] = 0;

					int J = m_xadj.Length - 1;
					part = new int[J];
					Debug.Assert(m_xadj.Where(idx => idx > m_adj.Length).Count() == 0);
					Debug.Assert(m_adj.Where(j => j >= J).Count() == 0);

					METIS.PARTGRAPHKWAY(
							ref J, ref ncon,
							m_xadj,
							m_adj.ToArray(),
							m_NoOfSpecies,
							null,
							null,
							ref NoOfParts,
							null,
							null,
							options,
							ref edgecut,
							part);
				} else {
					part = null;
				}

				//broadcast to every processor. (this is necessary to create column mapping)
				var partGlob = part.MPIBroadcast(worldMPIOffset, MGMapping.MPI_Comm);

				int[] i0part = new int[NoOfParts+1]; //new i0partitioning for each part
				for (int p = 1; p < NoOfParts+1; p++)
					i0part[p] = i0part[p - 1] + partGlob.Where(b => b == p - 1).Count();

                long[] cell_i0s = i0part.Select(i0 => (long)i0).ToArray();

				// create column mapping
				for (int b = 0; b < partGlob.Length; b++) {
					int targetRank = partGlob[b];
					columnMapping.Add((b, i0part[targetRank]));
					i0part[targetRank]++; //iterate the index for the next cell
				}

				return (cell_i0s, columnMapping);
			}
		}

		(long i0Cell, int lenCell)[] CellIndexToDOFs(MultigridMapping map) {
			int J = map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
			var myList = new (long i0Cell, int lenCell)[J]; //for cell currently on this proc (may be distributed to another or not)

			// Create matrix entry information from cell indices
			for (int jCellLoc = 0; jCellLoc < J; jCellLoc++) {
				myList[jCellLoc] = (map.GetCellI0(jCellLoc), map.GetLength(jCellLoc));
			}
			(long i0Cell, int lenCell)[][] globList = myList.MPIAllGatherO(map.MPI_Comm);
			var flatGlobArray = globList.SelectMany(x => x).ToArray();

			return flatGlobArray;
		}

	}

	public class StandAloneOperatorMappingPairWithGridData : IOperatorMappingPair {
		BlockMsrMatrix m_OpMtx;
		PseudoCoordinateMapping m_Mapping;
		public int[] m_xadj;
		public int[] m_adj;
		public int[] m_NoOfSpecies;
		public (long i0Cell, int lenCell)[] CellToDOFdata;

		/// <summary>
		/// ctor for the stand alone DG operator mapping pair
		/// Assumption: block indices are equilavent to cell indices
		/// </summary>
		/// <param name="OperatorMtx"></param>
		/// <param name="cellIndexMapping">old cell indicis to new cell indices for grid data (cells often get re-distributed) </param>
		/// <param name="xadj"></param>
		/// <param name="adj"></param>
		/// <param name="NoOfSpecies">a weight information</param>
		public StandAloneOperatorMappingPairWithGridData(BlockMsrMatrix OperatorMtx, List<(long source, long target)> cellIndexMapping, int[] xadj, int[] adj, int[] NoOfSpecies, bool GridDataOn = true) {
			m_OpMtx = OperatorMtx;

			Debug.Assert(cellIndexMapping.Count == OperatorMtx._RowPartitioning.TotalNoOfBlocks);
			m_Mapping = new PseudoCoordinateMapping(OperatorMtx._RowPartitioning);
			// if the master/leader proc, get the grid/csr data and convert to the new cell indices
			if (OperatorMtx._RowPartitioning.MpiRank == 0 && GridDataOn) {
				m_xadj = new int[xadj.Length];
				m_adj = new int[adj.Length];
				m_NoOfSpecies = new int[NoOfSpecies.Length];
				RemapCSRAndSpieces(xadj, adj, NoOfSpecies, cellIndexMapping, out m_xadj, out m_adj, out m_NoOfSpecies);
			}

			if (GridDataOn)
				CellToDOFdata = CellIndexToDOFs(m_Mapping);
		}

		/// <summary>
		/// A very ad-hoc and ugly way to deal with DOF data
		/// </summary>
		/// <param name="map"></param>
		/// <returns></returns>
		(long i0Cell, int lenCell)[] CellIndexToDOFs(IBlockPartitioning map) {
			int J = map.LocalNoOfBlocks;
			var myList = new (long i0Cell, int lenCell)[J]; //for cell currently on this proc (may be distributed to another or not)

			// Create matrix entry information from cell indices
			for (int jCellLoc = 0; jCellLoc < J; jCellLoc++) {
				myList[jCellLoc] = (map.GetBlockI0(jCellLoc), map.GetBlockLen(jCellLoc));
			}
			(long i0Cell, int lenCell)[][] globList = myList.MPIAllGatherO(map.MPI_Comm);
			var flatGlobArray = globList.SelectMany(x => x).ToArray();

			return flatGlobArray;
		}

		/// <summary>
		/// Remaps the CSR and values of the matrix to the new cell indices
		/// </summary>
		/// <param name="xadj"></param>
		/// <param name="adj"></param>
		/// <param name="NoOfSpecies"></param>
		/// <param name="cellIndexMapping"></param>
		/// <param name="newXadj"></param>
		/// <param name="newAdj"></param>
		/// <param name="newNoOfSpecies"></param>
		void RemapCSRAndSpieces(	int[] xadj, int[] adj, int[] NoOfSpecies,
								List<(long source, long target)> cellIndexMapping,
								out int[] newXadj, 	out int[] newAdj, out int[] newNoOfSpecies) {
			int N = xadj.Length - 1;
			// ensure mapping size matches
			Debug.Assert(cellIndexMapping.Count == N, "cellIndexMapping must cover every row");

			// 1) build old→new permutation array
			var perm = new int[N];
			foreach (var (src, tgt) in cellIndexMapping) {
				Debug.Assert(src >= 0 && src < N, $"src {src} out of range");
				Debug.Assert(tgt >= 0 && tgt < N, $"tgt {tgt} out of range");
				perm[src] = (int)tgt;
			}

#if DEBUG
			// ensure perm is a true permutation
			var seen = new bool[N];
			for (int i = 0; i < N; i++) {
				Debug.Assert(!seen[perm[i]], $"duplicate target {perm[i]}");
				seen[perm[i]] = true;
			}
#endif

			// 2) invert it: newRow→oldRow
			var invPerm = new int[N];
			for (int oldRow = 0; oldRow < N; oldRow++) {
				invPerm[perm[oldRow]] = oldRow;
			}

#if DEBUG
			// sanity‐check invertibility
			for (int i = 0; i < N; i++) {
				Debug.Assert(invPerm[perm[i]] == i, "invPerm not inverse of perm");
			}
#endif
			// 3) rebuild CSR
			var adjList = new List<int>(adj.Length);
			newXadj = new int[N + 1];
			for (int newRow = 0; newRow < N; newRow++) {
				newXadj[newRow] = adjList.Count;
				int oldRow = invPerm[newRow];
				Debug.Assert(oldRow >= 0 && oldRow < N, "oldRow out of range");
				for (int k = xadj[oldRow]; k < xadj[oldRow + 1]; k++) {
					int oldNb = adj[k];
					Debug.Assert(oldNb >= 0 && oldNb < N, $"neighbor {oldNb} out of range");
					adjList.Add(perm[oldNb]);
				}
			}
			newXadj[N] = adjList.Count;
			Debug.Assert(newXadj[0] == 0);
			Debug.Assert(newXadj[N] == adj.Length, "total edges count changed");
			newAdj = adjList.ToArray();
			Debug.Assert(newAdj.Length == adj.Length);

			// 4) remap NoOfSpecies
			newNoOfSpecies = new int[N];
			foreach (var (src, tgt) in cellIndexMapping) {
				newNoOfSpecies[tgt] = NoOfSpecies[src];
			}
			Debug.Assert(newNoOfSpecies.Length == N);
		}

		public BlockMsrMatrix OperatorMatrix => m_OpMtx;

		public ICoordinateMapping DgMapping => m_Mapping;

		class PseudoCoordinateMapping : BlockPartitioning, ICoordinateMapping {

			public PseudoCoordinateMapping(IBlockPartitioning blockPartitioning) : base(blockPartitioning) {
				// this is used for the stand alone MG solver
			}	

			public int NoOfVariables => throw new NotImplementedException();

			public int[] DgDegree => throw new NotImplementedException();

			public int NoOfExternalCells => throw new NotImplementedException();

			public int SpatialDimension => throw new NotImplementedException();

			public int NoOfLocalUpdatedCells => throw new NotImplementedException();

			public int LocalCellCount => throw new NotImplementedException();

			public SpeciesId[] UsedSpecies => throw new NotImplementedException();

			public int GetLength(int jLoc) {
				throw new NotImplementedException();
			}

			public int GetNoOfSpecies(int jCell) {
				throw new NotImplementedException();
			}

			public int GetSpeciesIndex(int jCell, SpeciesId SId) {
				throw new NotImplementedException();
			}

			public long GlobalUniqueIndex(int ifld, int jCell, int jSpec, int n) {
				throw new NotImplementedException();
			}

			public long GlobalUniqueIndex(int ifld, int jCell, int n) {
				throw new NotImplementedException();
			}

			public bool IsXDGvariable(int iVar) {
				throw new NotImplementedException();
			}

			public int LocalUniqueIndex(int ifld, int jCell, int iSpec, int n) {
				throw new NotImplementedException();
			}

			public int LocalUniqueIndex(int ifld, int jCell, int n) {
				throw new NotImplementedException();
			}
		}
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
		public BlockMsrMatrix OpMatrix => thisCommOpMatrix; // m_OpMapPair.OperatorMatrix;
		public BlockMsrMatrix WorldCommOpMatrix => m_OpMapPair.OperatorMatrix;

		TaskParallelMGOperator TpMapping => m_OpMapPair as TaskParallelMGOperator;

		/// <summary>
		/// passed in <see cref="InitImpl"/>
		/// </summary>
		IOperatorMappingPair m_OpMapPair;

        /// <summary>
        /// defines the problem matrix
        /// </summary>
        public void Init(IOperatorMappingPair op) {
            if (op is TaskParallelMGOperator TP)
                InitImpl(TP);
            else if (op is MultigridOperator MG)
                Init(MG);
            else 
				throw new NotImplementedException("TaskParallelOrthoMG: Init(IOperatorMappingPair) not implemented yet");
        }

        /// <summary>
        /// entry point for the TaskParallelOrthoMG (finest level)
        /// </summary>
        public void Init(MultigridOperator op) {
			Debugger.Launch();

			if (op.OperatorMatrix.MPI_Comm != csMPI.Raw._COMM.WORLD)
				throw new Exception("Task parallel OrthoMG (finest level) should be initiated with an operator in world communicator");

			this.FinerLevelCoarseComm = csMPI.Raw._COMM.WORLD;
			this.FinerLevelCoarseCommRank = op.Mapping.MpiRank;

			int WorldSize = op.Mapping.MpiSize;
			var thisTP = new TaskParallelMGOperator(op.OperatorMatrix, op.GetPrologonationOperator, op.Mapping, WorldSize);
			TaskParallelMGOperator finerTP = thisTP;
			if (verbose) {
				op.OperatorMatrix.SaveToTextFileSparseDebug($"OperatorMatrix_0.txt");
				op.OperatorMatrix.SaveToTextFileSparse($"OperatorMatrix_0.txt");
			}

			int level = 1;
			for (MultigridOperator op_lv = op.CoarserLevel; op_lv != null; op_lv = op_lv.CoarserLevel) {
				if (verbose) {
					op_lv.OperatorMatrix.SaveToTextFileSparseDebug($"OperatorMatrix_{level}.txt");
					op_lv.OperatorMatrix.SaveToTextFileSparse($"OperatorMatrix_{level}.txt");

					op_lv.GetPrologonationOperator.SaveToTextFileSparseDebug($"ProlongationMatrix_{level}.txt");
					op_lv.GetPrologonationOperator.SaveToTextFileSparse($"ProlongationMatrix_{level}.txt");
				}

				var coarserTP = new TaskParallelMGOperator(op_lv.OperatorMatrix, op_lv.GetPrologonationOperator.CloneAs(), op_lv.Mapping, WorldSize - level, finerTP);

				finerTP.CoarserLevel = coarserTP;
                //coarserTP.FinerLevel = finerTP;
				finerTP = coarserTP;
				level++;
			}
			InitImpl(thisTP);
			csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
		}

		BlockMsrMatrix ChangeCommunicator(BlockMsrMatrix mtx, (long i0Global, int CellLen)[] localBlocks, MPI_Comm comm, int[] MPIRankMapping ) {
			if (mtx == null)
				return null;

			var newColPartition = GetPartitioning(localBlocks, comm);
			mtx.ChangeMPICommForColumns(MPIRankMapping, newColPartition); // changes the communicator to the new one (subComm)

			List<long> ColIndices = Enumerable.Range((int)mtx._ColPartitioning.i0, mtx._ColPartitioning.LocalLength).Select(i => (long)i).ToList();
			List<long> RowIndices = Enumerable.Range( (int)mtx._RowPartitioning.i0 , mtx._RowPartitioning.LocalLength).Select(i => (long)i).ToList();

			return mtx.GetSubMatrix(RowIndices, ColIndices, comm); // finally create a new matrix with the new Comm
		}


		//BlockMsrMatrix opCommRestrictionOperator = null;
		BlockMsrMatrix worldCommProlongationOperator => TpMapping.ProlongationMatrix;
		BlockMsrMatrix worldCommFromCoarseProlongationOperator => TpMapping.CoarserLevel.ProlongationMatrix;

		BlockMsrMatrix subCommFromCoarseProlongationOperator = null;
		BlockMsrMatrix subCommToCoarseRestrictionOperator = null;

		BlockMsrMatrix thisCommProlongationOperator = null;
		BlockMsrMatrix thisCommRestrictionOperator = null;

		//BlockMsrMatrix opLeftChangeOfBasisMatrix = null;
		//BlockMsrMatrix opRightChangeOfBasisMatrix = null;

		//BlockMsrMatrix subCommRestrictionOperator = null;
		//BlockMsrMatrix subCommProlongationOperator = null;
		BlockMsrMatrix thisCommOpMatrix = null;
		BlockMsrMatrix subCommSmootherOpMatrix = null;
		BlockMsrMatrix subCommCoarseOpMatrix = null;

		(BlockMsrMatrix Matrix, List<long> RowIndices, List<long> ColIndices, BlockMsrMatrix TransposeMtx) smootherPermutation = (null, null, null, null);
		(BlockMsrMatrix Matrix, List<long> RowIndices, List<long> ColIndices, BlockMsrMatrix TransposeMtx) coarsePermutation = (null, null, null, null);

		BlockMsrMatrix smootherPermutationMtx => smootherPermutation.Matrix;
		BlockMsrMatrix coarsePermutationMtx => coarsePermutation.Matrix;

		List<(long source, long target)> columnMappingWorldToThis => TpMapping.ThisNewCellMapping;
		//{
		//	get {
		//		var tpFine = TpMapping.FinerLevel as TaskParallelMGOperator;
		//		if (tpFine != null) { 
		//			return tpFine.CoarseNewCellMapping;
		//		} else { 
		//			var ret = Enumerable.Range(0, columnMappingWorldToSmoother.Count).Select(i => ((long)i, (long)i)).ToList();
		//			Debug.Assert( ret.Count == columnMappingWorldToSmoother.Count && ret.Zip(columnMappingWorldToSmoother, (a, b) => a == b).All(equal => equal));
		//			return ret;
		//		}
		//	}
		//}


		List<(long source, long target)> columnMappingWorldToSmoother => TpMapping.SmootherNewCellMapping;
        //  new List<(long source, long target)>(); // cell-based indexes for the mapping from the operator to the smoother matrix (source idx, target idx)
        List<(long source, long target)> columnMappingWorldToCoarse => TpMapping.CoarseNewCellMapping; 
            //= new List<(long source, long target)>();   // cell-based indexes for the mapping from the operator to the coarse matrix (source idx, target idx)

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


		// To Do:
		// - implement the restriction operator to be called from finer level
		// - in old version, it is done my mg operator, this time it should be done by the solver

		void InitImpl(TaskParallelMGOperator op) {
            Debugger.Launch();  
			using (var tr = new FuncTrace()) {
                if (object.ReferenceEquals(op, m_OpMapPair))
                    return; // already initialized

				else
                    this.Dispose();

				csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);


				if (FinerLevelCoarseComm == csMPI.Raw._COMM.NULL)
					return; // this init called by the smoother of the finer, which is not needed. // it can be still needed for smoother


				this.m_OpMapPair = op;
                var Mtx = op.OperatorMatrix;

				var MgMap = op.DgMapping;

                //if (!Mtx.RowPartitioning.EqualsPartition(MgMap))
                //    throw new ArgumentException("Row partitioning mismatch.");
                //if (!Mtx.ColPartition.EqualsPartition(MgMap))
                //    throw new ArgumentException("Column partitioning mismatch.");


				NoOfCoarseProcs = op.NoOfCoarseProcs;
				NoOfSmootherProcs = op.NoOfSmootherProcs;
				SplitCommunicator(NoOfCoarseProcs);
				//PermutateOpMatrixFromScratch();
				MigrateFromWorldToSubComms();
				ortho = new CoreOrthonormalizationProcedureTP(OpMatrix);
				CreatePermutationMatrices();

				//CommitCoarseningOperators(opCommRestrictionOperator, worldCommProlongationOperator);

				// set operator
				// ============

				// initiate coarser level
				// ======================
				if (TpLevel == 0)
					csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
				InitCoarse();
				if (TpLevel == 0) {
					csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
					double[] a = new double[] { worldCommRank };
					a = a.MPIMax(csMPI.Raw._COMM.WORLD);
				}



			}
        }

		void InitSmoothers(TaskParallelMGOperator op) {
			var SmootherOpMappingPairOnSubComm = new StandAloneOperatorMappingPairWithGridData(subCommSmootherOpMatrix, columnMappingWorldToSmoother, TpMapping.m_xadj, TpMapping.m_adj, TpMapping.m_NoOfSpecies);

			// init smoother
			// =============
			if (PreSmoother != null) {
				if (PreSmoother is SchwarzForTaskParallel schwarzTp) {
					schwarzTp.Init(SmootherOpMappingPairOnSubComm);
				} else if (PreSmoother is ISubsystemSolver ssPreSmother) {
					ssPreSmother.Init(SmootherOpMappingPairOnSubComm);
				} else {
					throw new NotSupportedException($"Unable to initialize pre-smoother if it is not a {typeof(ISubsystemSolver)} and until multigrid operator is designed to work on sub communicators, this won't work");
				}
			}

			// init post smoother
			// =============
			if (PostSmoother != null && !object.ReferenceEquals(PreSmoother, PostSmoother)) {
				if (PostSmoother is SchwarzForTaskParallel schwarzTp) {
					schwarzTp.Init(SmootherOpMappingPairOnSubComm);
				} else if (PostSmoother is ISubsystemSolver ssPostSmother) {
					ssPostSmother.Init(SmootherOpMappingPairOnSubComm);
				} else {
					throw new NotSupportedException($"Unable to initialize pre-smoother if it is not a {typeof(ISubsystemSolver)} and until multigrid operator is designed to work on sub communicators, this won't work");
				}
			}
		}

		double[] Restrict<T1>(T1 IN_fine)
				where T1 : IList<double> {
			if (TpMapping.CoarserLevel == null)
				throw new NotSupportedException("Already on finest level -- no finer level to restrict from.");

			if (IN_fine.Count != subCommToCoarseRestrictionOperator.ColPartition.LocalLength)
				throw new ArgumentException("Mismatch in length of fine grid vector (input).", "IN_fine");

			double[] OUT_coarse = new double[subCommToCoarseRestrictionOperator.RowPartitioning.LocalLength];
			this.subCommToCoarseRestrictionOperator.SpMV(1.0, IN_fine, 0.0, OUT_coarse);

			//if (this.LeftChangeOfBasis != null) {
			//	double[] LB = new double[OUT_coarse.Count];
			//	this.LeftChangeOfBasis.SpMV(1.0, OUT_coarse, 0.0, LB);
			//	OUT_coarse.SetV(LB);
			//}
			return OUT_coarse;
		}

		double[] Prolongate<T1>(T1 IN_coarse)
		where T1 : IList<double> {
			if (TpMapping.CoarserLevel == null)
				throw new NotSupportedException("Already on finest level -- no finer level to prolongate from.");

			if (IN_coarse.Count != subCommFromCoarseProlongationOperator.ColPartition.LocalLength)
				throw new ArgumentException("Mismatch in length of coarse grid vector (input).", "IN_coarse");

			double[] OUT_fine = new double[subCommFromCoarseProlongationOperator.RowPartitioning.LocalLength];
			this.subCommToCoarseRestrictionOperator.SpMV(1.0, IN_coarse, 0.0, OUT_fine);

			//if (this.LeftChangeOfBasis != null) {
			//	double[] LB = new double[OUT_coarse.Count];
			//	this.LeftChangeOfBasis.SpMV(1.0, OUT_coarse, 0.0, LB);
			//	OUT_coarse.SetV(LB);
			//}
			return OUT_fine;
		}

		ICoordinateMapping MgMap => m_OpMapPair.DgMapping;

		/// <summary>
		/// uses permutation matrix to permutate the operator matrix, only works for the matrix in this level of MG and if it is square
		/// </summary>
		void PermutateOpMatrixFromScratch() {
			using (var tr = new FuncTrace()) {
				subCommSmootherOpMatrix = ChangeDistributionAndCommunicators(smootherBlocks, WorldCommOpMatrix, smootherPermutation, columnMappingWorldToSmoother, true, $"s_lvl{TpLevel}");
				subCommCoarseOpMatrix = ChangeDistributionAndCommunicators(coarseBlocks, WorldCommOpMatrix, coarsePermutation, columnMappingWorldToCoarse, true, $"c_lvl{TpLevel}");
			}
		}

		int TpLevel => TpMapping.level;

		/// <summary>
		/// Changes the communicators fo the operator matrix (to smoother subComm) and the prolongation operator (to this levelComm - to the coarse subComm of the finer level)
		/// (the matrix entries should be already distributed via <see cref="TaskParallelMGOperator"/>
		/// It is not need to deal with the coarse mesh as it will be carried out by the coarser level 
		/// (except the coarsest level <see cref="InitiateCoarsestSolver"/>)
		/// </summary>
		void MigrateFromWorldToSubComms() {
			using (var tr = new FuncTrace()) {
				thisCommOpMatrix = ChangeCommunicator(WorldCommOpMatrix.CloneAs(), TpMapping.localBlocksForThisLevel, thisComm, WorldToThisCommMapping);

				subCommSmootherOpMatrix = ChangeCommunicator(TpMapping.m_OpMtx_smoother.CloneAs(), TpMapping.localBlocksForSmoother, subComm, WorldToSubCommMapping);
				thisCommProlongationOperator = ChangeCommunicator(worldCommProlongationOperator, TpMapping.localBlocksForThisLevel, thisComm, thisCommMap);
				thisCommRestrictionOperator = thisCommProlongationOperator?.Transpose();

				// In MG operator, the prolongation and restriction operators are stored at the coarse level and restriction performed by operator
				// However, this is not the case here. The MG operator is dedicated for world level communications and key informations.
				// This then must be done by this level solver here for the coarser solver.
				// So each level solver gets data already restricted and then return without prolongating.
				subCommFromCoarseProlongationOperator = ChangeCommunicator(worldCommFromCoarseProlongationOperator, TpMapping.CoarserLevel.localBlocksForThisLevel, subComm, WorldToSubCommMapping);
				subCommToCoarseRestrictionOperator = subCommFromCoarseProlongationOperator?.Transpose();

				if (verbose) {
					subCommSmootherOpMatrix.SaveToTextFileSparseDebug($"lvl_{TpLevel}_o_smoothOp.txt");
					subCommSmootherOpMatrix.SaveToTextFileSparse($"lvl_{TpLevel}_o_smoothOp.txt");
					thisCommProlongationOperator?.SaveToTextFileSparseDebug($"lvl_{TpLevel}_thisProl.txt");
					thisCommProlongationOperator?.SaveToTextFileSparse($"lvl_{TpLevel}_thisProl.txt");
					thisCommRestrictionOperator?.SaveToTextFileSparseDebug($"lvl_{TpLevel}_thisRest.txt");
					thisCommRestrictionOperator?.SaveToTextFileSparse($"lvl_{TpLevel}_thisRest.txt");
					subCommFromCoarseProlongationOperator?.SaveToTextFileSparseDebug($"lvl_{TpLevel}_subProl.txt");
					subCommFromCoarseProlongationOperator?.SaveToTextFileSparse($"lvl_{TpLevel}_subProl.txt");
					subCommToCoarseRestrictionOperator?.SaveToTextFileSparseDebug($"lvl_{TpLevel}_subRest.txt");
					subCommToCoarseRestrictionOperator?.SaveToTextFileSparse($"lvl_{TpLevel}_subRest.txt");
					//tests the old distribution done from world all to world this procs in the TaskParallelMGOperator, not this level to smoother or 
					//TestMatrices(TpMapping.old_OpMtx, smootherPermutation.Matrix, OpMatrix, "test_", verbose); 
				}
			}
		}

		double[] ChangeVectorDistr(double[] vec, IBlockPartitioning newPart) {
			double[] ret = new double[newPart.LocalLength];
			return ret;
		}

		double[] ChangeVectorDistrBack(double[] vec, IBlockPartitioning newPart) {
			double[] ret = new double[newPart.LocalLength];
			return ret;
		}

		double[] PermutateVector(double[] vec, (BlockMsrMatrix Matrix, List<long> RowIndices, List<long> ColIndices, BlockMsrMatrix TransposeMtx) permutation) {
			double[] ret = new double[permutation.Matrix.RowPartitioning.LocalLength];
			permutation.Matrix.SpMV(1.0, vec, 0.0, ret);
			return ret;
		}

		double[] PermutateVectorBack(double[] vec, (BlockMsrMatrix Matrix, List<long> RowIndices, List<long> ColIndices, BlockMsrMatrix TransposeMtx) permutation) {
			double[] ret = new double[permutation.TransposeMtx.RowPartitioning.LocalLength];
			permutation.TransposeMtx.SpMV(1.0, vec, 0.0, ret);
			return ret;
		}

		/// <summary>
		/// Changes the communicator of <paramref name="Mtx"/> w.r.t. <paramref name="permutation"/> and distributes them to communicator of <paramref name="subComm"/>
		/// matrices are the same but permutated so that they are stored only on their communicators
		/// </summary>
		/// <param name="Mtx"></param>
		/// <param name="permutation"></param>
		/// <returns></returns>
		BlockMsrMatrix ChangeDistributionAndCommunicators((long i0Global, int CellLen)[] localBlocks, BlockMsrMatrix Mtx, (BlockMsrMatrix Matrix, List<long> RowIndices, List<long> ColIndices, BlockMsrMatrix TransposeMtx) permutation, IList<(long Source, long Target)> columnMapping, bool test = true, string tag = "d_") {
			//smoother (should be called by all procs in opComm)
			var LocalizedMatrix = BlockMsrMatrix.Multiply(permutation.Matrix, Mtx); // transfer matrix to the target procs (changed the position of rows not columns)
			LocalizedMatrix.ChangeColumnIndices(columnMapping); //switches column indices with respect to the mapping
			LocalizedMatrix.ChangeColumnPartitioning(LocalizedMatrix._RowPartitioning); // calculates new external indices with respec to the new partitioning (still opComm)
			var newColMapping = GetPartitioning(localBlocks, subComm);
			LocalizedMatrix.ChangeMPICommForColumns(ThisCommToSubCommMapping, newColMapping); // changes the communicator to the new one (subComm)
			var permutatedMatrix = LocalizedMatrix.GetSubMatrix(permutation.RowIndices, permutation.RowIndices, subComm); // finally create a new matrix with the new Comm
#if DEBUG
			if (test)
				TestMatrices(Mtx, permutation.Matrix, permutatedMatrix, tag, verbose);
#endif
			return permutatedMatrix;
		}

		(long i0Cell, int lenCell)[] smootherBlocks;
		(long i0Cell, int lenCell)[] coarseBlocks;

		IBlockPartitioning thisPartitioningInThisComm => thisCommOpMatrix._RowPartitioning; //thisCommProlongationOperator == null ? TpMapping.ThisTargetPartitioning : thisCommProlongationOperator._RowPartitioning; // if this instance is for the finest level without prolongation operator, then we can use the this Partitioning at the world level. If not, get the this partitioning from the prolongation operator (comm operator has changed).

		/// <summary>
		/// Permutation matrices from thisComm to subComms (to disribute vectors)
		/// </summary>
		void CreatePermutationMatrices() {
			var worldToCoarseDict = columnMappingWorldToCoarse.ToDictionary(x => x.source, x => x.target);
			var worldToSmootherDict = columnMappingWorldToSmoother.ToDictionary(x => x.source, x => x.target);

			List<(long source, long target)> colMapThisToSmoother = columnMappingWorldToThis
				.Select(x => (x.target, worldToSmootherDict[x.source]))
				.ToList();

			List<(long source, long target)> colMapThisToCoarse = columnMappingWorldToThis
				.Select(x => (x.target, worldToCoarseDict[x.source]))
				.ToList();

			smootherBlocks = GetLocalDistribution(TpMapping.CellToDOFdata, colMapThisToSmoother, TpMapping.SmootherCellI0s, 0, NoOfSmootherProcs);
			Console.WriteLine("I was here");
			csMPI.Raw.Barrier(thisComm);
			smootherPermutation = GetPermutationMatrix(thisPartitioningInThisComm, smootherBlocks); //this is technically a permutation matrix but also distributes 

			coarseBlocks = GetLocalDistribution(TpMapping.CellToDOFdata, colMapThisToCoarse, TpMapping.CoarseCellI0s, NoOfSmootherProcs, NoOfCoarseProcs);
			coarsePermutation = GetPermutationMatrix(thisPartitioningInThisComm, coarseBlocks);
		}

  //      List<(long source, long target)> columnMappingOpToSmootherOld;
		//List<(long source, long target)> columnMappingOpToCoarseOld;

        /// <summary>
        /// tests if the permutation matrix is correct with respect to the original matrix and the subComm matrix
        /// </summary>
        /// <param name="opCommMtx">a matrix defined on opComm</param>
        /// <param name="PermutationMtx">permutation matrix to be tested</param>
        /// <param name="subCommMtx">the subComm version of <paramref name="opCommMtx"/> </param>
        /// <param name="tag">tag for i/o outputs</param>
        /// <param name="WriteToFiles">true or false</param>
        /// <exception cref="Exception"></exception>
        void TestMatrices(BlockMsrMatrix opCommMtx, BlockMsrMatrix PermutationMtx, BlockMsrMatrix subCommMtx, string tag = "t_", bool WriteToFiles = false) {
			//Create a test vector
			Random rnd = new Random(44); //seed is 44
			var TestVector = new double[opCommMtx.ColPartition.LocalLength];
			for (int i = 0; i < TestVector.Length; i++) {
				TestVector[i] = rnd.NextDouble();
			}

			// first test SpMV with opCommMatrix
			var opResult = new double[opCommMtx.RowPartitioning.LocalLength];
			opCommMtx.SpMV(1.0, TestVector, 0.0, opResult);
				
			// test vector to subComm
			var TestVectorSub = new double[PermutationMtx._RowPartitioning.LocalLength];
			PermutationMtx.SpMV(1.0, TestVector, 0.0, TestVectorSub);

			// result at the subComm
			var subOpResult = new double[subCommMtx._RowPartitioning.LocalLength];
			subCommMtx.SpMV(1.0, TestVectorSub, 0.0, subOpResult);

			// map back the subComm result to the opComm
			var TransposeRedist = PermutationMtx.Transpose();
			var subOpResultBackToOp = new double[TransposeRedist._RowPartitioning.LocalLength];
			TransposeRedist.SpMV(1.0, subOpResult, 0.0, subOpResultBackToOp);

			if (WriteToFiles) {
				//dist.SaveToTextFileDebug($"{tag}localDistForSmoother.txt");
				//GlobDOFidx.SaveToTextFileDebug($"{tag}GlobDOFidx.txt");

				opCommMtx.SaveToTextFileSparseDebug($"{tag}OpMatrix.txt"); //parallel
				opCommMtx.SaveToTextFileSparse($"{tag}OpMatrix.txt"); /// combine and write as a single file

				TestVector.SaveToTextFileDebug($"{tag}TestVector", ".txt");
				opResult.SaveToTextFileDebug($"{tag}opResult", ".txt");

				subCommMtx.SaveToTextFileSparseDebug($"{tag}localBlock.txt");
				subCommMtx.SaveToTextFileSparse($"{tag}localBlock.txt");
				
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

		(long i0Cell, int lenCell)[] GetLocalDistribution((long i0Cell, int lenCell)[] globalDOFs, List<(long Source, long Target)> cellColumnMapping, long[] targeti0s, int procOffset, int procSize) {
            int newRank = thisCommRank - procOffset;
            if (newRank < 0 || newRank >= procSize)
                return new (long i0Cell, int lenCell)[0];

            var ret = new List<(long i0Cell, int lenCell)>();
			var newBlocki0 = (int)targeti0s[newRank];
            var newBlockiE = (int)targeti0s[newRank + 1];

			// if cell is designated to be on this processor, add it to the list
			for (long iCell = 0; iCell < cellColumnMapping.Count; iCell++) {
				var iTargetCell = cellColumnMapping[(int)iCell].Target; //designated cell index
				var iSourceCell = cellColumnMapping[(int)iCell].Source;

				if (iTargetCell >= newBlocki0 && iTargetCell < newBlockiE) { // if designated cell index falls into the range of this processor
					var sourceBlock = globalDOFs[iSourceCell]; // get the current block index (source)
					ret.Add(sourceBlock);
				}
			}
			ret.Sort((x, y) => x.i0Cell.CompareTo(y.i0Cell)); // sort the list according to the cell index

            Debug.Assert(newBlockiE - newBlocki0 == ret.Count);
			return ret.ToArray();
		}

        /// <summary>
        /// Creates permutation matrix and calculates Row and col indixes of each DOF (not block as in cellColumnMapping)
        /// </summary>
        /// <param name="mapping">MG mapping for grid data at the current level</param>
        /// <param name="locDOFsData">DOFs data with global indices</param>
        /// <returns></returns>
		private (BlockMsrMatrix Matrix, List<long> RowIndices, List<long> ColIndices, BlockMsrMatrix TransposeMatrix) GetPermutationMatrix(IBlockPartitioning mapping, (long i0Global, int CellLen)[] locDOFsData) {
			using (new FuncTrace()) {

				BlockPartitioning TargetPartitioning = GetPartitioning(locDOFsData, thisComm);
				BlockMsrMatrix PermutationMatrix = new BlockMsrMatrix(TargetPartitioning, mapping);
                List<long> RowIndices = new List<long>(); 
                List<long> ColIndices= new List<long>();
				{
					int cnt = 0;
						if (locDOFsData != null && locDOFsData.Length > 0) {
							int NoCells = locDOFsData.Length;

							for (int j = 0; j < NoCells; j++) {
								//int Len = part.lenCell[j];
								int Len = locDOFsData[j].CellLen;
								long i0Row = TargetPartitioning.i0 + cnt;
								long i0Col = locDOFsData[j].i0Global;   //part.i0Cell[j];

								PermutationMatrix.AccBlock(i0Row, i0Col, 1.0, MultidimensionalArray.CreateEye(Len));

								for (int k = 0; k < Len; k++) {
									RowIndices.Add(i0Row + k);
									ColIndices.Add(i0Col + k);
									//columnMappingWorldToSmoother.Add((i0Col + k, i0Row + k));
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
                var TransposeMtx = PermutationMatrix.Transpose();

				return (PermutationMatrix, RowIndices, ColIndices, TransposeMtx);
			}
		}

		private int GetRank(MPI_Comm op_comm) {
			csMPI.Raw.Comm_Rank(op_comm, out int m_rank);
			return m_rank;
		}

		MPI_Comm thisComm => FinerLevelCoarseComm != csMPI.Raw._COMM.NULL ? FinerLevelCoarseComm : csMPI.Raw._COMM.WORLD;
		int thisCommRank => FinerLevelCoarseComm != csMPI.Raw._COMM.NULL ? FinerLevelCoarseCommRank : TpMapping.worldCommRank;
		int[] thisCommMap => FinerLevelCoarseComm != csMPI.Raw._COMM.NULL ? FinerLevelCoarseCommMap : Enumerable.Range(0, TpMapping.worldCommSize).ToArray();

		int thisCommsize {
			get {
				csMPI.Raw.Comm_Size(thisComm, out int size);
				return size;
			}
		}

		MPI_Comm subComm;

		public MPI_Comm FinerLevelCoarseComm = csMPI.Raw._COMM.NULL;
		public int FinerLevelCoarseCommRank;
		public int[] FinerLevelCoarseCommMap;

		int subCommRank; 
		int subCommSize; 

		MPI_Comm opComm => m_OpMapPair.OperatorMatrix.MPI_Comm;
		int opCommRank => m_OpMapPair.OperatorMatrix._RowPartitioning.MpiRank;
		int opCommSize => m_OpMapPair.OperatorMatrix._RowPartitioning.MpiSize;

		readonly int worldCommRank;

		/// <summary>
		/// Creates a new communicator for the coarse level solver.
		/// </summary>
		/// <param name="NoOfCoarseProcs"></param>
		void SplitCommunicator(int NumberOfCoarseProcs) {
			if (NoOfCoarseProcs < 1 || NoOfSmootherProcs < 1)
                throw new ArgumentOutOfRangeException($"Number of coarse processors and smoother processors must be greater than 0. Coarse: {NoOfCoarseProcs}, Smoother: {NoOfSmootherProcs}");

			if (thisCommRank < NoOfSmootherProcs) {
				myTask = TpTaskType.Smoother;
				csMPI.Raw.CommSplit(thisComm, 0, FinerLevelCoarseCommRank, out subComm);
			} else {
				myTask = TpTaskType.Coarse;
				csMPI.Raw.CommSplit(thisComm, 1, FinerLevelCoarseCommRank, out subComm);
			}

			csMPI.Raw.Comm_Rank(subComm, out subCommRank);
			csMPI.Raw.Comm_Size(subComm, out subCommSize);

            if(verbose)
			    Console.WriteLine($"The proc with worldRank-{worldCommRank} with opRank{opCommRank} is assigned to {subCommRank}of{subCommSize} new com{subComm.m1}");

            CreateMpiRankMappings();
		}


		/// <summary>
		/// ?we do not really need this as MPI_Split reorders new communicators w.r.t. old ranks
		/// </summary>
		void CreateMpiRankMappings() {
			// Each process sends its rank as its value.
			Debug.Assert(thisCommsize == TpMapping.NoOfThisProcs);

			int totalCount = thisCommsize; //  number of processes in the communicator

			// Prepare managed buffers.
			int[] sendBuffer = new int[] { subCommRank }; //workaround for fixed statement
			ThisCommToSubCommMapping = new int[totalCount];

			// Pin the buffers to obtain their addresses.
			GCHandle sendHandle = GCHandle.Alloc(sendBuffer, GCHandleType.Pinned);
			GCHandle recvHandle = GCHandle.Alloc(ThisCommToSubCommMapping, GCHandleType.Pinned);

			IntPtr sendPtr = sendHandle.AddrOfPinnedObject();
			IntPtr recvPtr = recvHandle.AddrOfPinnedObject();

			// Perform the Allgather call.
			csMPI.Raw.Allgather(sendPtr, 1, csMPI.Raw._DATATYPE.INT, recvPtr, 1, csMPI.Raw._DATATYPE.INT, thisComm);

			// Display gathered values on each process.
			Console.WriteLine($"Process {thisCommRank} gathered:");
			for (int i = 0; i < totalCount; i++) {
				Console.WriteLine($"Rank {i}: {ThisCommToSubCommMapping[i]}");
			}

			var CommToSubCommMappingD = new int[totalCount];
			var WorldToSubCommMappingD = new int[TpMapping.worldCommSize];
			WorldToSubCommMappingD.SetAll(-1);
			for (int i = 0; i < totalCount; i++) {
				CommToSubCommMappingD[i] = i < TpMapping.NoOfSmootherProcs ? i : i - TpMapping.NoOfSmootherProcs;
				WorldToSubCommMappingD[i + TpMapping.worldMPIOffset] = i < TpMapping.NoOfSmootherProcs ? i : i - TpMapping.NoOfSmootherProcs;
			}

			WorldToSubCommMapping = new int[TpMapping.worldCommSize];
			WorldToSubCommMapping.SetAll(-1);
			WorldToSubCommMapping = new int[TpMapping.worldCommSize];
			WorldToSubCommMapping.SetAll(-1); // keep empty ranks -1 to see if an error comes out, as the default value (0) might cause some unvisible bugs
			WorldToSubCommMapping.SetSubVector(ThisCommToSubCommMapping, TpMapping.worldMPIOffset, TpMapping.NoOfThisProcs);

			Debug.Assert(WorldToSubCommMappingD.SequenceEqual(WorldToSubCommMapping));

			WorldToThisCommMapping = Enumerable.Range(0, TpMapping.worldCommSize).Select(i => i < TpMapping.worldMPIOffset ? -1 : i - TpMapping.worldMPIOffset).ToArray();
			//WorldToSubCommMapping = ThisCommToSubCommMapping;// ??
		}

		bool verbose = true;

        int[] ThisCommToSubCommMapping;

		int[] WorldToSubCommMapping;
		int[] WorldToThisCommMapping;

		TpTaskType myTask = TpTaskType.All;

		BlockPartitioning GetPartitioning((long i0Global, int CellLen)[] DOFs, MPI_Comm comm) {
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
				var partitioning = new BlockPartitioning(LocalLength, i0Cell, LnCell, comm, i0isLocal: true);
				return partitioning;
			}
		}
		
        private void InitiateCoarsestSolver() {
			var CoarserTpMapping = (TpMapping.CoarserLevel as TaskParallelMGOperator);
			subCommCoarseOpMatrix = ChangeCommunicator(CoarserTpMapping.OperatorMatrix.CloneAs(), CoarserTpMapping.localBlocksForThisLevel, subComm, WorldToSubCommMapping);
			if (CoarserLevelSolver is PARDISOSolver PARSolver) {

				PARSolver.DefineMatrix(subCommCoarseOpMatrix);
			} else if (CoarserLevelSolver is ISubsystemSolver subsystemSolver) {
				var SmootherOpMappingPairOnSubComm = new StandAloneOperatorMappingPairWithGridData(subCommCoarseOpMatrix, TpMapping.CoarserLevel.CoarseNewCellMapping, TpMapping.m_xadj, TpMapping.m_adj, TpMapping.m_NoOfSpecies, false);
				subsystemSolver.Init(SmootherOpMappingPairOnSubComm);
			} else {
				throw new NotSupportedException();
			}
			//csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
			// for the final coarse level solver (probabaly direct solver), we need to actually change the communicator of its matrix
			//throw new NotImplementedException("Not yet implemented.");
		}

		private void InitCoarse() {
			using (var tr = new FuncTrace()) {
				if (CoarserLevelSolver == null)
					throw new NotSupportedException("Missing coarse level solver.");

				if (TpMapping.CoarserLevel == null)
					throw new NotSupportedException("Unexpected null CoarserLevel.");

				if (CoarserLevelSolver is TaskParallelOrthoMG ssCoarse) { 
					ssCoarse.FinerLevelCoarseComm = myTask != TpTaskType.Smoother ? subComm : csMPI.Raw._COMM.NULL;
					ssCoarse.FinerLevelCoarseCommRank = myTask != TpTaskType.Smoother ? subCommRank : -1;
					ssCoarse.FinerLevelCoarseCommMap = myTask != TpTaskType.Smoother ? FinerLevelCoarseCommMap : null;

					ssCoarse.Init(TpMapping.CoarserLevel);
				} else
					InitiateCoarsestSolver();
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
                Debug.Assert(Res.Length == OpMatrix.ColPartition.LocalLength);
                Debug.Assert(X.Length == OpMatrix.ColPartition.LocalLength);
                Debug.Assert(B.Length == OpMatrix.ColPartition.LocalLength);
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
                f.StdoutOnAllRanks();

				ThisLevelTime.Start();

				// Initialize vectors
				double[] B = InitializeVector(_B);
				double[] X = InitializeVector(_xl);
				int L = X.Length;
				double[] Res = new double[L];

				// Initialize residual
				double[] Sol0 = X.CloneAs();
				double[] Res0 = InitializeResidual(X, B, Res);
				Residual(Res0, Sol0, B);
				Array.Copy(Res0, Res, L);

				// Initialize task-specific data
				double[] XforSub, BforSub, ResforSub;
				InitializeTaskSpecificData(X, B, Res, out XforSub, out BforSub, out ResforSub);


				csMPI.Raw.Barrier(thisComm);

                if(ortho.Norm(Res0) <= 0) {
                    double normB = ortho.Norm(B);
                    double normX = ortho.Norm(X);
                    f.Error($" **** Residual is 0.0: |X| = {normB}; |partGlob| = {normB}; |Res0| = {ortho.Norm(Res0)}");
                }

                int iLevel = TpLevel;

                void WriteDebug(int iter, double res, string text) {
                    if (iLevel >= 0)
					    Console.WriteLine($"{string.Concat(Enumerable.Repeat("-", iLevel))} OrthoMG, current level={iLevel}, iteration={iter} {(text != null ? " - " + text : "")} and res norm: {res}");
					
                    return;
                }

                double iter0_resNorm = ortho.Norm(Res0);
                double resNorm = iter0_resNorm;
                //this.IterationCallback?.Invoke(0, Sol0, Res0, this.m_OpMapPair as MultigridOperator);
                                
                // clear history of coarse solvers
                ortho.Clear();



				if (myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {
					InitSmoothers(TpMapping);
				}


				PerformIterations(X, B, Res, Res0, XforSub, BforSub, ResforSub);



				int LforSub = XforSub.Length;

				int iIter;
				//           for (iIter = 1; bIterate; iIter++) {
				//               WriteDebug(iIter, resNorm, "initial start");
				//csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

				//var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
				//               if (!termState.bNotTerminate) {
				//                   Converged = termState.bSuccess;
				//                   break;
				//               } else {

				//               }


				//// pre-smoother
				//// ------------
				//if (myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {

				//	VerivyCurrentResidual(XforSub, BforSub, Res, iIter);
				//	double[] PreCorr = new double[L];

				//	if (thisCommRank == SmootherGroupLeader) {
				//		byte completionSignal = 0;
				//		MPI_Request req;
				//		unsafe {
				//			byte* pSignal = &completionSignal;
				//			IntPtr ptr = (IntPtr)pSignal;
				//			csMPI.Raw.Irecv(ptr, 1,
				//							csMPI.Raw._DATATYPE.BYTE,
				//							CoarseGroupLeader, 22,
				//							thisComm, out req);

				//		}
				//		int k = 1;
				//		while (!done) {
				//			PreSmoother.Solve(PreCorr, Res); // Vorglättung
				//			Thread.Sleep(1000);
				//			Console.WriteLine($"{k}-second waiting for completionSignal={completionSignal}");
				//			csMPI.Raw.Test(ref req, out done, out MPI_Status status);
				//			done.MPIOr(subComm);
				//			k++;
				//		}
				//		Console.WriteLine($"waited {k}-second for completionSignal={completionSignal}");

				//	} else {
				//		while (!done) {
				//			PreSmoother.Solve(PreCorr, Res); // Vorglättung
				//			Thread.Sleep(1000);
				//			done.MPIOr(subComm);
				//		}

				//	}

				//                   if (PreSmoother is SchwarzForCoarseMesh schwarz2) { 
				//		schwarz2.PropagateBackToOpComm(PreCorr);
				//	}
				//	// orthonormalization and residual minimization
				//	resNorm = ortho.AddSolAndMinimizeResidual(ref PreCorr, X, Sol0, Res0, Res, "presmoothL" + iLevel);

				//	WriteDebug(iIter, resNorm, " pre-smoother applied");

				//	//SpecAnalysisSample(iIter, X, "ortho1");
				//	var termState2 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
				//	if (!termState2.bNotTerminate) {
				//		Converged = termState2.bSuccess;
				//		break;
				//	}

				//}





				//               // coarse grid correction
				//               // ----------------------
				//               if (myTask == TpTaskType.All || myTask == TpTaskType.Coarse) {

				//	if (thisCommRank == CoarseGroupLeader) {
				//		Console.WriteLine("Waiting for 20 seconds...");
				//		Thread.Sleep(20000);  // 20 seconds delay
				//		Console.WriteLine("Done waiting.");

				//		//while (!done) {
				//		//	Thread.Sleep(1000);
				//		//	done.MPIOr(newComm);
				//		//}

				//		byte completionSignal = 1;
				//		unsafe {
				//			byte* pSignal = &completionSignal;
				//			IntPtr ptr = (IntPtr)pSignal;
				//			csMPI.Raw.Send(ptr, 1, csMPI.Raw._DATATYPE.BYTE, SmootherGroupLeader, 22, thisComm);
				//		}
				//		Console.WriteLine("Sent signal");
				//	} else {
				//		Thread.Sleep(19000);
				//		done = true;

				//	}
				//	Thread.Sleep(2000000);

				//	CrseLevelTime.Start();
				//                   // Test: Residual on this level / already computed by 'MinimizeResidual' above
				//                   VerivyCurrentResidual(X, B, Res, iIter);

				//                   double resNorm_b4Coarse = resNorm;
				//                   for (int i = 0; i < myConfig.m_omega; i++) {
				//                       if (this.CoarserLevelSolver != null) {

				//                           double[] vl = new double[L];
				//                           if (myConfig.CoarseOnLowerLevel) {

				//                               var _MgOperator = m_OpMapPair as MultigridOperator;

				//                               // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//                               // coarse grid solver defined on COARSER MESH LEVEL:
				//                               // this solver must perform restriction and prolongation
				//                               // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//                               using (new BlockTrace("Restriction", f)) {
				//                                   // restriction of residual
				//                                   _MgOperator.CoarserLevel.Restrict(Res, ResCoarse);
				//                               }
				//                               // Berechnung der Grobgitterkorrektur
				//                               double[] vlc = new double[Lc];
				//                               this.CoarserLevelSolver.Solve(vlc, ResCoarse);
				//                               using (new BlockTrace("Prolongation", f)) {
				//                                   // Prolongation der Grobgitterkorrektur
				//                                   _MgOperator.CoarserLevel.Prolongate(1.0, vl, 1.0, vlc);
				//                               }

				//                           } else {
				//                               // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//                               // coarse grid solver defined on the SAME MESH LEVEL:
				//                               // performs (probably) its own restriction/prolongation
				//                               // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

				//                               // Berechnung der Grobgitterkorrektur
				//                               this.CoarserLevelSolver.Solve(vl, Res);
				//                           }

				//                           // orthonormalization and residual minimization
				//                           if (vl.ContainsNanOrInf().MPIOr() && CoarseArithmeticExceptionCount < 20) {
				//                               CoarseArithmeticExceptionCount++;
				//                               Console.WriteLine("Coarse solver failed " + CoarseArithmeticExceptionCount);
				//                               if (CoarseArithmeticExceptionCount == 1) {
				//                                   BlockMsrMatrix coarseMtx = myConfig.CoarseOnLowerLevel ? (m_OpMapPair as MultigridOperator).CoarserLevel.OperatorMatrix : m_OpMapPair.OperatorMatrix;
				//                                   coarseMtx.SaveToTextFileSparse("FailedCoarseMatrix.txt");
				//                               }
				//                               vl = null;
				//                               this.CoarserLevelSolver.Dispose();
				//                               this.InitCoarse();
				//                           } else {
				//                               resNorm = ortho.AddSolAndMinimizeResidual(ref vl, X, Sol0, Res0, Res, "coarsecorL" + iLevel);
				//                           }

				//                       }
				//                   } // end of coarse-solver loop
				//                   WriteDebug(iIter, resNorm, "coarse-solver applied");

				//                   var termState3 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
				//                   if (!termState3.bNotTerminate) {
				//                       Converged = termState3.bSuccess;
				//                       break;
				//                   }
				//                   CrseLevelTime.Stop();                       
				//}

				//               if (myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {

				//                   // post-smoother
				//                   // -------------
				//                   if (PostSmoother != null || AdditionalPostSmoothers != null) {
				//                       bool termPost = false;



				//                       for (int g = 0; g < config.NoOfPostSmootherSweeps; g++) { // Test: Residual on this level / already computed by 'MinimizeResidual' above

				//                           ISolverSmootherTemplate _PostSmoother = allSmooters[iPostSmooter];// g % allSmooters.Length];
				//                           if (iPostSmooter > 0)
				//                               InitAdditionalPostSmooters();
				//                           VerivyCurrentResidual(X, B, Res, iIter); // 
				//                           double[] PostCorr = new double[L];

				//                           bool fail = false;
				//                           try {
				//                               _PostSmoother.Solve(PostCorr, Res); // compute correction (Nachglättung)
				//                           } catch (ArithmeticException ae) {
				//                               f.Error($"Smoother fail on Rank {m_OpMapPair.DgMapping.MpiRank}: " + ae.ToString());
				//                               fail = true;
				//                           }
				//                           fail = fail.MPIOr();

				//                           if ((PostCorr.ContainsNanOrInf() || fail).MPIOr()) {
				//                               PostSmootherArithmeticExceptionCount++;
				//                               Console.Error.WriteLine("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g + " for " + PostSmootherArithmeticExceptionCount + "-th time.");

				//                               if (Res.ContainsNanOrInf()) {
				//                                   Console.WriteLine("... so does RHS");
				//                                   throw new ArithmeticException("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g + " RHS already corrupted.");
				//                               } else
				//                                   Console.WriteLine("... although RHS is regular");

				//                               //var viz = new MGViz(m_OpMapPair as MultigridOperator);
				//                               //var __RHS = viz.ProlongateRhsToDg(Res, "RES");
				//                               //var __SOL = viz.ProlongateRhsToDg(PostCorr, "SOL");
				//                               //Tecplot.Tecplot.PlotFields(__SOL.Cat(__RHS), "iilufail", 0.0, 0);
				//                               //

				//                               if (PostSmootherArithmeticExceptionCount < 20) {
				//                                   PostCorr = null;
				//                                   _PostSmoother.Dispose();
				//                                   if (PostSmoother is ISubsystemSolver ssPostSmother) {
				//                                       ssPostSmother.Init(this.m_OpMapPair);
				//                                   } else {
				//                                       if (this.m_OpMapPair is MultigridOperator mgOp) {
				//                                           PostSmoother.Init(mgOp);
				//                                       } else {
				//                                           throw new NotSupportedException($"Unable to initialize post-smoother if it is not a {typeof(ISubsystemSolver)} and operator is not a {typeof(MultigridOperator)}");
				//                                       }
				//                                   }
				//                               } else {
				//                                   throw new ArithmeticException("Post-Smoother " + _PostSmoother.GetType().Name + " produces NAN/INF at sweep " + g);
				//                               }

				//                           } else {
				//                               resNorm = ortho.AddSolAndMinimizeResidual(ref PostCorr, X, Sol0, Res0, Res, "pstsmthL" + iLevel + "-sw" + g);
				//                               WriteDebug(iIter, resNorm, "post-smoother applied");



				//                               var termState4 = TerminationCriterion(iIter, iter0_resNorm, resNorm);
				//                               if (!termState4.bNotTerminate) {
				//                                   Converged = termState4.bSuccess;
				//                                   termPost = true;
				//                                   break;
				//                               }

				//                               if (ortho.CancellationTriggered) {
				//                                   //skipPreSmoothInFollowingIters = true; // most of the time, the pre-smoother does nothing different than the post-smoother; So, if we cancel post-smoothing there is no need to do pre-smoothing in the next loop.

				//                                   iPostSmooter++;
				//                                   if (iPostSmooter >= allSmooters.Length)
				//                                       iPostSmooter = 0;




				//                                   break;
				//                               }
				//                           }
				//                       }


				//                       if (termPost)
				//                           break;

				//                   } // end of post-smoother loop
				//               }


				//               // iteration callback
				//               // ------------------
				//               this.ThisLevelIterations++;
				//               IterationCallback?.Invoke(iIter, X, Res, this.m_OpMapPair as MultigridOperator);

				//           } // end of solver iterations
				iIter = 0;
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

		private double[] InitializeVector<T>(T input) where T : IList<double> {
			if (input is double[] array) {
				return array;
			}
			return input.ToArray();
		}

		private double[] InitializeResidual(double[] X, double[] B, double[] Res) {
			double[] Res0 = new double[X.Length];
			Residual(Res0, X, B);
			Array.Copy(Res0, Res, Res0.Length);

			if (ortho.Norm(Res0) <= 0) {
				double normB = ortho.Norm(B);
				double normX = ortho.Norm(X);
				throw new ArithmeticException($"Residual is 0.0: |X| = {normX}; |B| = {normB}; |Res0| = {ortho.Norm(Res0)}");
			}

			return Res0;
		}

		private void InitializeTaskSpecificData(double[] X, double[] B, double[] Res, out double[] XforSub, out double[] BforSub, out double[] ResforSub) {
			using (var trace = new FuncTrace()) {
				// Validate input
				if (X == null || B == null || Res == null) {
					throw new ArgumentNullException("Input vectors X, B, or Res cannot be null.");
				}

				// Permutation matrices live on thisComm, which means that they should be called here
				//var XforSubCoarse = PermutateVector(X, coarsePermutation);
				//var BforSubCoarse = PermutateVector(B, coarsePermutation);
				var ResforSubCoarse = PermutateVector(Res, coarsePermutation);

				var XforSubSmoother = PermutateVector(X, smootherPermutation);
				var BforSubSmoother = PermutateVector(B, smootherPermutation);
				//var ResforSubSmoother = PermutateVector(Res, smootherPermutation);

				// Initialize task-specific data based on the task type
				switch (myTask) {
					case TpTaskType.Coarse:
						// For coarse tasks, only the residual is permuted
						ResforSub = ResforSubCoarse;
						XforSub = null; // Not used for coarse tasks
						BforSub = null; // Not used for coarse tasks
						break;

					case TpTaskType.Smoother:
						// For smoother tasks, all vectors are permuted
						XforSub = XforSubSmoother;
						BforSub = BforSubSmoother;
						ResforSub = null;
						break;

					case TpTaskType.All:
						// For all tasks, use the original vectors
						XforSub = X;
						BforSub = B;
						ResforSub = Res;
						break;

					default:
						// Handle unexpected task types
						throw new InvalidOperationException($"Unexpected task type: {myTask}");
				}

				// Log the initialization
				trace.Info($"Task-specific data initialized for task type: {myTask}");
			}
		}


		private void PerformIterations(double[] X, double[] B, double[] Res, double[] Res0, double[] XforSub, double[] BforSub, double[] ResforSub) {
			int iIter;
			double iter0_resNorm = ortho.Norm(Res0);
			double resNorm = iter0_resNorm;
			csMPI.Raw.Barrier(thisComm);
			var X0 = X.CloneAs();
			for (iIter = 1; ; iIter++) {
				bool done = false;
				int SmootherGroupLeader = 0 + NoOfCoarseProcs - NoOfCoarseProcs;
				int CoarseGroupLeader = NoOfCoarseProcs + NoOfSmootherProcs - 1;
				done.MPIAnd(thisComm);

				// Check termination
				var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
				if (!termState.bNotTerminate) {
					Converged = termState.bSuccess;
					break;
				}

				double[] PreCorrSub = new double[smootherPermutationMtx._RowPartitioning.LocalLength];
				// Pre-smoother
				if (myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {
					PreCorrSub = ApplyPreSmoother(XforSub, BforSub);
				}
				var PreCorr = PermutateVectorBack(PreCorrSub, smootherPermutation);
				resNorm = ortho.AddSolAndMinimizeResidual(ref PreCorr, X, X0, Res0, Res, "Tp-presmooth" + TpLevel);
				// Log the result
				Console.WriteLine($"Pre-smoother is applied successfully. Residual norm: {resNorm}");

				double[] CoarseCorreectionOnSub = new double[coarsePermutationMtx._RowPartitioning.LocalLength];
				// Coarse grid correction
				if (myTask == TpTaskType.All || myTask == TpTaskType.Coarse) {
					CoarseCorreectionOnSub =  ApplyCoarseGridCorrection(ResforSub);
				}
				var CoarseCorrection = PermutateVectorBack(CoarseCorreectionOnSub, coarsePermutation);
				resNorm = ortho.AddSolAndMinimizeResidual(ref CoarseCorrection, X, X0, Res0, Res, "Tp-presmooth" + TpLevel);
				Console.WriteLine($"Coarse-correction is applied successfully. Residual norm: {resNorm}");


				// Post-smoother
				if (myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {
					//ApplyPostSmoother(X, Res, Res0, ref resNorm, iIter);
				}

				csMPI.Raw.Barrier(thisComm);


				// Iteration callback
				IterationCallback?.Invoke(iIter, X, Res, m_OpMapPair as MultigridOperator);
			}
		}

		private double[] ApplyPreSmoother(double[] X, double[] B) {
			using (var trace = new FuncTrace()) {
				// Verify that the residual is up-to-date

				// Apply the pre-smoother
				PreSmoother?.Solve(X, B); // Pre-smoother modifies PreCorr based on Res

				Console.WriteLine($"Pre-smoother applied on rank {thisCommRank}");

				return X;
			}
		}


		private double[] ApplyCoarseGridCorrection (double[] Res) {
			using (var trace = new FuncTrace()) {
				if (CoarserLevelSolver == null) {
					throw new InvalidOperationException("Coarse level solver is not initialized.");
				}

				// Start timing for coarse grid correction
				CrseLevelTime.Start();

				// Restrict the residual to the coarse grid
				double[] ResCoarse = Restrict(Res);

				// Solve the coarse grid problem
				double[] CoarseCorrection = new double[ResCoarse.Length];
				CoarserLevelSolver.Solve(CoarseCorrection, ResCoarse);

				// Prolongate the correction back to the fine grid
				double[] FineCorrection = Prolongate(CoarseCorrection);

				// Stop timing for coarse grid correction
				CrseLevelTime.Stop();

				return FineCorrection;
			}
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

        //int PostSmootherArithmeticExceptionCount = 0;
        //int CoarseArithmeticExceptionCount = 0;


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

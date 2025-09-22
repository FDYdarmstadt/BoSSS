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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Control;
using BoSSS.Solution.Gnuplot;
using BoSSS.Solution.Statistic;
using BoSSS.Solution.Tecplot;
using ilPSP;
using ilPSP.Kraypis;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.monkey;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using ilPSP.Utils;
using log4net.Core;
using MPI.Wrappers;
using Newtonsoft.Json.Linq;
using NUnit.Common;
using NUnit.Framework.Internal;
using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Net.NetworkInformation;
using System.Numerics;
using System.Reflection.Emit;
using System.Runtime.InteropServices;
using System.Runtime.Serialization;
using System.Security.Policy;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Linq;
using static System.Net.Mime.MediaTypeNames;
using static System.Reflection.Metadata.BlobBuilder;

namespace BoSSS.Solution.AdvancedSolvers {

	[StructLayout(LayoutKind.Sequential)]
	internal struct BlockInfo {
		public long iBlock;
		public long i0Cell;
		public int lenCell;
		internal BlockInfo(long iBlock, long i0Cell, int lenCell) {
			this.iBlock = iBlock;
			this.i0Cell = i0Cell;
			this.lenCell = lenCell;
		}
	}

	// To Do: write a CoreOrthonormalization class for the task parallel variant that works only on smoother processors.
	// There is no need for ThisLevel Matrices. They are only used for SpMV but the matrix operations are expensive.
	// Some room for improvement is there.

	/// <summary>
	/// Helper class designed to hold information getting from the WORLD.
	/// It distributes all the necessary information to the responsible processors on the world communicator.
	/// Later on, they can be easily copied/transferred to the sub communicators.
	/// </summary>
	public class TaskParallelMGOperator : IOperatorMappingPair {

		/// <summary>
		/// An ad-hoc information storing class for World-Level MG operator
		/// Constructor (responsible for world to sub communicators)
		/// should be defined on world communicator and then will be used on sub communicators with the solver
		/// Basically this defines the operator matrix and the prolongation matrix 
		/// and transfers them to the responsible processors without changing the communicator
		/// Technically each level is responsible for the prolongation matrix to the finer level and the operator matrix of the current level
		/// Every level works with three different communicators: level from finer level (coarse comm of the finer), 
		/// and the smoother and the coarse for this level (still on World but the processors are designated)
		/// In each level first processors are for the smoother and the rest is for the coarse,
		/// leading that 0-th rank is always the smoother with the finest and last rank is for the coarse on the coarsest level
		/// This can be called from the main program or from a finer level of another solver, so the finest level might have or not a prolongation op.
		/// </summary>
		/// <param name="OperatorMatrix"></param>
		/// <param name="ProlongationMatrix"></param>
		/// <param name="Mapping"></param>
		/// <param name="TotProcSize"></param>
		/// <param name="CoarseProcSize"></param>
		/// <param name="finerLevel"></param>
		/// <param name="LeftChangeOfBasis"></param>
		/// <param name="RightChangeOfBasis"></param>
		/// <param name="IsThereACoarserSolver"></param>
		public TaskParallelMGOperator(BlockMsrMatrix OperatorMatrix, BlockMsrMatrix ProlongationMatrix, MultigridMapping Mapping, int TotProcSize = 0, int CoarseProcSize = 0, IOperatorMappingPair finerLevel = null, BlockMsrMatrix LeftChangeOfBasis = null, BlockMsrMatrix RightChangeOfBasis = null, bool IsThereACoarserSolver = true) { 
            m_OpMtx = OperatorMatrix;
			m_ProlMtx = ProlongationMatrix;
			m_MultigridMapping = Mapping;
			NoOfThisProcs = TotProcSize != 0 ? TotProcSize :
										finerLevel is TaskParallelMGOperator tp ? // if zero check if finer level TP
										tp.NoOfCoarseProcs : finerLevel.OperatorMatrix.RowPartitioning.MpiSize;
			NoOfCoarseProcs = CoarseProcSize != 0 ? CoarseProcSize : Math.Max(1, Math.Min(3 * NoOfThisProcs / 4, NoOfThisProcs - 1));
			m_IsThereCoarserLevel = IsThereACoarserSolver;
			FinerLevel = finerLevel;
			m_LeftChangeOfBasis = LeftChangeOfBasis;
			m_RightChangeOfBasis = RightChangeOfBasis;	
			CheckMPICorrectness();

			(m_xadj,m_adj) = GetCurrentAggGridGraphForMetis(m_MultigridMapping);
            m_NoOfSpecies = GetNoOfSpeciesList(m_MultigridMapping);

			//distribution at cell/block level
			(ThisCellI0s, ThisNewCellMapping) = DistributeMapping(m_MultigridMapping, NoOfThisProcs); 
			(SmootherCellI0s, SmootherNewCellMapping) = DistributeMapping(m_MultigridMapping, NoOfSmootherProcs);
			(CoarseCellI0s, CoarseNewCellMapping) = DistributeMapping(m_MultigridMapping, NoOfCoarseProcs);
			CalculateWorldToSubDistribution();
		}

		public IOperatorMappingPair FinerLevel;
		bool m_IsThereCoarserLevel; //this must be know at the constructor time,
									//as it checks the number of processors for the coarser level (if false, no need for sub smoother level)
		public TaskParallelMGOperator CoarserLevel; // can be attaached later on
		public BlockMsrMatrix m_OpMtx;
		public BlockMsrMatrix m_OpMtx_smoother;
		public BlockMsrMatrix m_ProlMtx;
		public MultigridMapping m_MultigridMapping;
		public BlockMsrMatrix m_LeftChangeOfBasis;
		public BlockMsrMatrix m_RightChangeOfBasis;
		public int Level => FinerLevel is TaskParallelMGOperator fine ? fine.Level + 1 : InitLevel; //level for this type of solver not overall level
																							//if this is attached to some finer level
		public int InitLevel = 0;

		public int NoOfThisProcs;
		public int NoOfSmootherProcs => NoOfThisProcs - NoOfCoarseProcs;
		public int NoOfCoarseProcs;

		public IBlockPartitioning ThisTargetPartitioning;               // this level partitioning on the world communicator
		public IBlockPartitioning CoarseTargetPartitioning;             // partitioning for coarse solver on the world communicator
                                                                        // (this is equal to this level partitioning of the coarser level)
		public IBlockPartitioning SmootherTargetPartitioning;           // partitioning for smoother on the world communicator
		internal (long i0Cell, int lenCell)[] localBlocksForThisLevel;  //block data from the original matrix on the world communicator
		internal (long i0Cell, int lenCell)[] localBlocksForSmoother;   //block data from the original matrix on the world communicator
		internal (long i0Cell, int lenCell)[] localBlocksForCoarse;     //block data from the original matrix on the world communicator

		public Dictionary<long, long> ThisNewCellMapping;     //global data (all of them) and global indices
		public Dictionary<long, long> CoarseNewCellMapping;   //global data (all of them) and global indices
		public Dictionary<long, long> SmootherNewCellMapping; //global data (all of them) and global indices

		public long[] ThisCellI0s;     //stores the global indices of starting cell indices of the processors for this level (length = NoOfThisProcs + 1)
		public long[] CoarseCellI0s;   //stores the global indices of starting cell indices of the processors for this level (length = NoOfCoarseProcs + 1)
		public long[] SmootherCellI0s; //stores the global indices of starting cell indices of the processors for this level (length = NoOfSmootherProcs + 1)

		// the following three arrays are used for the metis distribution
		// and stored only on the processor with rank=worldOffset (will be the master on the subcommunicator)
		public int[] m_xadj;
		public int[] m_adj;
		public int[] m_NoOfSpecies; 

		// properties for the necessary information
		public IBlockPartitioning FinerLevelCoarsePartitiong => FinerLevel is TaskParallelMGOperator fine ? fine.CoarseTargetPartitioning : null;
        public IBlockPartitioning FinerLevelThisPartitiong => FinerLevel is TaskParallelMGOperator fine ? fine.ThisTargetPartitioning : null;

        public (long i0Cell, int lenCell)[] FinerLevelCoarseBlocks => FinerLevel is TaskParallelMGOperator fine ? fine.localBlocksForCoarse : null;

        public (long i0Cell, int lenCell)[] FinerLevelThisBlocks => FinerLevel is TaskParallelMGOperator fine ? fine.localBlocksForThisLevel : null;

        public MPI_Comm currentComm => m_MultigridMapping.MPI_Comm;

		/// <summary>
		/// the offset of this operator matrix in the world communicator, the mpi rank of the first processor in this level
		/// </summary>
		public int worldMPIOffset => worldCommSize - NoOfThisProcs; 
		public int worldCommRank => m_MultigridMapping.MpiRank;
		public int worldCommSize => m_MultigridMapping.MpiSize;
		public BlockMsrMatrix OperatorMatrix => m_OpMtx;
		public BlockMsrMatrix ProlongationMatrix => m_ProlMtx;
		public ICoordinateMapping DgMapping => m_MultigridMapping;


		bool verbose = false;

		void CheckMPICorrectness() {
			Console.WriteLine($"TaskParallel-MG level {Level}: {NoOfSmootherProcs} smoother + {NoOfCoarseProcs} coarse = {NoOfThisProcs} total procs");
			if (m_MultigridMapping.MPI_Comm != csMPI.Raw._COMM.WORLD)
				throw new Exception("The mg mapping must be defined on the WORLD communicator.");

			if (NoOfThisProcs < 1)
				throw new Exception("The number of processors must be equal or greater than 1.");

			if (NoOfThisProcs < 2 && m_IsThereCoarserLevel)
				throw new Exception("Not enough number of processors left for the coarser level.");

			if (NoOfSmootherProcs < 1 && m_IsThereCoarserLevel)
				throw new Exception("The number of processors for the smoother must be greater than 0. " +
					"Only allowed if this level is direct solver (i.e., this level is the coarsest level)");

			if (FinerLevel != null)
				if (FinerLevel is TaskParallelMGOperator TpFine) {
					if (TpFine.NoOfCoarseProcs != NoOfThisProcs)
						throw new Exception("Inconsistent no of processors in the nested initialization");

					if ((TpFine.NoOfSmootherProcs + TpFine.worldMPIOffset) != worldMPIOffset)
						throw new Exception("Inconsistent no of processors resulting incorrect world offset");
				} else {
					if (FinerLevel.OperatorMatrix.RowPartitioning.MpiSize != NoOfThisProcs)
						throw new Exception("Inconsistent no of processors in the nested initialization");
				}

		}

		void CalculateWorldToSubDistribution() {
			localBlocksForThisLevel = GetLocalDistribution(ThisNewCellMapping, ThisCellI0s, worldMPIOffset, NoOfThisProcs);
			localBlocksForCoarse = GetLocalDistribution(CoarseNewCellMapping, CoarseCellI0s, worldMPIOffset + NoOfSmootherProcs, NoOfCoarseProcs);

			ThisTargetPartitioning = GetPartitioning(localBlocksForThisLevel, currentComm);
			CoarseTargetPartitioning = GetPartitioning(localBlocksForCoarse, currentComm);

			if (verbose) {
				ThisNewCellMapping.SaveToTextFileDebug($"lvl{Level}_thisNewMapping", ".txt");
				SmootherNewCellMapping.SaveToTextFileDebug($"lvl{Level}_SmootherNewCellMapping", ".txt");
				CoarseNewCellMapping.SaveToTextFileDebug($"lvl{Level}_CoarseNewCellMapping", ".txt");
			}
				
			if (m_IsThereCoarserLevel) { //this level of solving is only for smoother so op should partitioned for smoother (reserved for optimization)
				localBlocksForSmoother = GetLocalDistribution(SmootherNewCellMapping, SmootherCellI0s, worldMPIOffset, NoOfSmootherProcs);
				SmootherTargetPartitioning = GetPartitioning(localBlocksForSmoother, currentComm);
				m_OpMtx_smoother = ChangeThisLevelPartitioning(m_OpMtx, SmootherTargetPartitioning, localBlocksForSmoother, "OpSmooth");
				thisLevelPermutationMtx = null;
				thisLevelPermutationMtxTranspose = null;
			} // [ToprakToDo]: this got nasty. Solve the problem here.

			m_OpMtx = ChangeThisLevelPartitioning(m_OpMtx, ThisTargetPartitioning, localBlocksForThisLevel, "Op"); //this level of 
			// if null, these will return null too.
			m_LeftChangeOfBasis = ChangeThisLevelPartitioning(m_LeftChangeOfBasis, ThisTargetPartitioning, localBlocksForThisLevel, "LeftCofB");
			m_RightChangeOfBasis = ChangeThisLevelPartitioning(m_RightChangeOfBasis, ThisTargetPartitioning, localBlocksForThisLevel, "RightCofB");

            m_ProlMtx = ChangeProlPartitioning(ThisTargetPartitioning, localBlocksForThisLevel, FinerLevelCoarsePartitiong, FinerLevelCoarseBlocks, "Pro"); //same processors different level of mg operator (different number of cells for row and column)

            thisLevelPermutationMtx = null;
			thisLevelPermutationMtxTranspose = null;
		}

		BlockMsrMatrix thisLevelPermutationMtx = null;
		BlockMsrMatrix thisLevelPermutationMtxTranspose = null;

		//reserved for optimization [ToprakToDo] (even permutation matrices are not needed. it must be jsut a copy/paste operation with mpi exchange)
		BlockMsrMatrix ChangeThisLevelPartitioning(BlockMsrMatrix Mtx, IBlockPartitioning targetPartitioning, (long i0Cell, int lenCell)[] blockData = null, string tag = "Op") {
			if (Mtx is null)
				return null;


			if (verbose) {
				Mtx.SaveToTextFileSparseDebug($"lvl{Level}_{tag}_oldOp.txt");
				Mtx.SaveToTextFileSparse($"lvl{Level}_{tag}_oldOp.txt");
			}

			BlockMsrMatrix ret; //notice that a null check here might disturb the mpi behavior as 

			if (Mtx.RowPartitioning == targetPartitioning && Mtx.ColPartition == targetPartitioning && blockData == null)
				//this can happen at the first level of tp
				ret = Mtx.CloneAs();
			else {
				 
				bool needsNewPermutation = thisLevelPermutationMtx is null ||  //if null
					!targetPartitioning.Equals(thisLevelPermutationMtx._RowPartitioning) || //if not designed for the same target(can happen for smoother)
					!Mtx._RowPartitioning.Equals(thisLevelPermutationMtx._ColPartitioning); //if not designed for the same source

				if (needsNewPermutation) {
					thisLevelPermutationMtx = CreateThisLevelPermutationMtx(Mtx._RowPartitioning, targetPartitioning, blockData);
					thisLevelPermutationMtxTranspose = thisLevelPermutationMtx.Transpose();
				}

				var temp = BlockMsrMatrix.Multiply(thisLevelPermutationMtx, Mtx);
				ret = BlockMsrMatrix.Multiply(temp, thisLevelPermutationMtxTranspose);
			}

			if (verbose) {
				ret.SaveToTextFileSparseDebug($"lvl{Level}_{tag}_newOp.txt");
				ret.SaveToTextFileSparse($"lvl{Level}_{tag}_newOp.txt");
			}

			return ret;
		}

		BlockMsrMatrix CreateThisLevelPermutationMtx(IBlockPartitioning sourcePartitioning, IBlockPartitioning targetPartitioning, (long i0Cell, int lenCell)[]  locDOFsData ) {
			BlockMsrMatrix PermutationMatrix = new BlockMsrMatrix(targetPartitioning, sourcePartitioning);
			{
				int cnt = 0;
				if (locDOFsData != null && locDOFsData.Length > 0) {
					int NoCells = locDOFsData.Length;

					for (int j = 0; j < NoCells; j++) {
						//int Len = part.lenCell[j];
						int Len = locDOFsData[j].lenCell;
						long i0Row = targetPartitioning.i0 + cnt;
						long i0Col = locDOFsData[j].i0Cell;

						PermutationMatrix.AccBlock(i0Row, i0Col, 1.0, MultidimensionalArray.CreateEye(Len));

						cnt += Len;
					}
				}
			}
			return PermutationMatrix;
		}

		BlockMsrMatrix ChangeProlPartitioning(IBlockPartitioning thisPartitioning, (long i0Cell, int lenCell)[] thisBlockData, IBlockPartitioning finerPartitioning, (long i0Cell, int lenCell)[] finerBlockData, string tag = "c") {
			if (m_ProlMtx == null)
				return null;

			if (verbose) {
				m_ProlMtx.SaveToTextFileSparse($"lvl{Level}_{tag}_oldPro.txt");
				m_ProlMtx.SaveToTextFileSparseDebug($"lvl{Level}_{tag}_oldPro.txt");
			}

			BlockMsrMatrix ret;
			if (finerPartitioning is null)
				ret = m_ProlMtx.ChangeColumnPartitioning2(thisPartitioning, thisBlockData);
			else
				ret = m_ProlMtx.ChangePartitioning(finerPartitioning, finerBlockData, thisPartitioning, thisBlockData);


			if (verbose) {
				ret.SaveToTextFileSparse($"lvl{Level}_{tag}_newPro.txt");
				ret.SaveToTextFileSparseDebug($"lvl{Level}_{tag}_newPro.txt");
			}

			return ret;
		}

		BlockPartitioning GetPartitioning((long i0Global, int CellLen)[] DOFs, MPI_Comm comm) {
			using (new FuncTrace()) {
				var LnCell = new List<int>();
				var i0Cell = new List<long>();

				int cnt = 0;

				if (DOFs != null && DOFs.Length > 0) {
					int NoCells = DOFs.Length;// part.jGlobal.Length;
					for (int j = 0; j < NoCells; j++) {
						int Len = DOFs[j].CellLen; // part.lenCell[j];
						if (Len > 0) {
							i0Cell.Add(cnt);
							LnCell.Add(Len);
							cnt += Len;
						} else { //keep block information to be consistent with global indices but assign them zero length and some suitable i0
							if (j > 0) {
								i0Cell.Add(i0Cell[j - 1]);
							} else {
								i0Cell.Add(0);
							}

							LnCell.Add(Len);
							cnt += Len;
						}
					}
				}

				int LocalLength = cnt;
				var partitioning = new BlockPartitioning(LocalLength, i0Cell, LnCell, comm, i0isLocal: true);
				return partitioning;
			}
		}

		(long i0Cell, int lenCell)[] GetLocalDistribution(Dictionary<long, long> cellColumnDict, long[] targeti0s, int procOffset, int procSize) {
			int B = m_MultigridMapping.LocalNoOfBlocks;

			List<BlockInfo> thisProcBlocks = new List<BlockInfo>();
			Dictionary<int, List<BlockInfo>> sendList = new Dictionary<int, List<BlockInfo>>();

			for (int b = 0; b < B; b++) { //loop for each block, and prepare necessary lists
				long iBlock = m_MultigridMapping.FirstBlock + b;
				long iTargetBlock = cellColumnDict[iBlock];
				int localRank = Array.BinarySearch(targeti0s, iTargetBlock) is var p && p >= 0 ? p : (~p) - 1;  // this is the rank of the target block in the target partitioning
				int TargetRank = procOffset + localRank; // added offset since we are still on world processor
				Debug.Assert(TargetRank >= procOffset && TargetRank < procOffset + procSize);

				BlockInfo block = new BlockInfo(iTargetBlock, m_MultigridMapping.GetBlockI0(iBlock), m_MultigridMapping.GetBlockLen(iBlock));

				if (TargetRank == worldCommRank) { 
					thisProcBlocks.Add(block);
				} else { //if the target block is another processor, add it to the send list
					if (!sendList.TryGetValue(TargetRank, out var list))
						sendList[TargetRank] = list = new List<BlockInfo>();
					list.Add(block);
				}
			}

			var sendArray = sendList.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.ToArray());
			IDictionary<int, BlockInfo[]> data = ArrayMessenger<BlockInfo>.ExchangeData(sendArray);

			foreach (var message in data) {
				Debug.Assert(message.Key != worldCommRank);
				thisProcBlocks.AddRange(message.Value);
			}

			thisProcBlocks.Sort((x, y) => x.iBlock.CompareTo(y.iBlock));
			return thisProcBlocks.Select(b => (b.i0Cell, b.lenCell)).ToArray();
		}

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

		public (long[] cellI0s, Dictionary<long, long> colMapping) DistributeMapping(MultigridMapping MGMapping, int NoOfParts) {
			using (var f = new FuncTrace()) {
				var columnMapping = new List<(long sourceIdx, long targetIdx)>();
                int totalLength = (int)MGMapping.AggGrid.CellPartitioning.TotalLength;

                if(NoOfParts == 1) { //nothing to distribute, just return the identity map
                    var identityMap = Enumerable.Range(0, totalLength).ToDictionary(i => (long)i, i => (long)i);
                    return (new long[] { 0, totalLength }, identityMap);
                }

                if(NoOfParts == worldCommSize) { //no need to re-distribute (current dist. should be fine), just return the identity map
                    var identityMap = Enumerable.Range(0, totalLength).ToDictionary(i => (long)i, i => (long)i);
                    return (MGMapping.AggGrid.CellPartitioning.GetOffsetArray(), identityMap);
                }

                int MPIrnk = MGMapping.MpiRank;
				int[] part;
				if (MPIrnk == worldMPIOffset) { //call metis on only the rank 0 (opComm)
                    f.Info($"{MPIrnk}-rank is calculating the re-distribution for the level-{Level} into {NoOfParts} parts");

                    int ncon = 1;
					int edgecut = 0;
					int[] options = new int[METIS.METIS_NOPTIONS];
					METIS.SETDEFAULTOPTIONS(options);

					options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = 1; // 
					options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 5; // This is the default refinement iterations
					options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 100; // Maximum imbalance of 10 percent (this is the default kway clustering)
                                                                               // 3 percent seems to be to strict for some test cases
					options[(int)METIS.OptionCodes.METIS_OPTION_NUMBERING] = 0;
                    options[(int)METIS.OptionCodes.METIS_OPTION_SEED] = 22; // fixed seed

                    int J = m_xadj.Length - 1;
					part = new int[J];
					Debug.Assert(m_xadj.Where(idx => idx > m_adj.Length).Count() == 0);
					Debug.Assert(m_adj.Where(j => j >= J).Count() == 0);

					METIS.PARTGRAPHKWAY(
							ref J, ref ncon,
							m_xadj,
							m_adj.ToArray(),
                            NoOfParts < 150 ? m_NoOfSpecies : null,
							null,
							null,
							ref NoOfParts,
							null,
							null,
							options,
							ref edgecut,
							part);
                    f.Info($"{MPIrnk}-rank is now broadcasting the distribution map for the level-{Level} into other ranks");
                } else {
                    f.Info($"{MPIrnk}-rank is waiting the distribution map for the level-{Level}");
                    part = null;
				}
                //broadcast to every processor. (this is necessary to create column mapping)
                var partGlob = part.MPIBroadcast(worldMPIOffset, MGMapping.MPI_Comm);
                f.Info($"Distribution is complete over the comm");

                if (partGlob.Length >= NoOfParts)
                    for(int p = 0; p < NoOfParts; p++)
					    if (!partGlob.Contains(p))
						    throw new InvalidOperationException("There are empty procs with TaskParallelOrthoMG. Either you are using too many processors such that there is not enough DOFs left for a processor or something odd is happening");

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

				return (cell_i0s, columnMapping.ToDictionary(pair => pair.sourceIdx, pair => pair.targetIdx));
			}
		}

        /// <summary>
        /// Prints the graph partitioning statistics for debugging purposes.
        /// </summary>
        /// <param name="xadj"></param>
        /// <param name="adj"></param>
        /// <param name="nodeWeights"></param>
        /// <param name="noOfParts"></param>
        /// <param name="f"></param>
        /// <param name="level"></param>
        public static void PrintGraphPartitioningStats(int[] xadj, int[] adj, int[] nodeWeights, int noOfParts, FuncTrace f, int level) {
            int numVertices = xadj.Length - 1;
            int numEdges = adj.Length;
            double avgDegree = numVertices > 0 ? (double)numEdges / numVertices : 0;

            f.Info("==== METIS Graph Partitioning Stats ====");
            f.Info($"Number of vertices (nodes): {numVertices}");
            f.Info($"Number of edges (entries in adj): {numEdges}");
            f.Info($"Average degree: {avgDegree:F2}");
            f.Info($"Number of partitions requested: {noOfParts}");

            xadj.SaveToTextFileDebug($"xadj_lvl{level}_{noOfParts}_parts.txt", ".txt");
            adj.SaveToTextFileDebug($"adj_lvl{level}_{noOfParts}_parts.txt", ".txt");
            nodeWeights.SaveToTextFileDebug($"nodeWeights_lvl{level}_{noOfParts}_parts.txt", ".txt");


            if(nodeWeights != null) {
                if(nodeWeights.Length != numVertices) {
                    f.Info("⚠️  Mismatch: nodeWeights.Length != number of vertices!");
                } else {
                    int minW = nodeWeights.Min();
                    int maxW = nodeWeights.Max();
                    f.Info($"Node weights: min = {minW}, max = {maxW}");
                }
            } else {
                f.Info("Node weights: none");
            }

            // Check if any xadj index is out of bounds
            bool invalidXadj = xadj.Any(idx => idx < 0 || idx > adj.Length);
            f.Info($"xadj indices valid: {!invalidXadj}");

            // Check if adj entries refer to valid node indices
            bool invalidAdj = adj.Any(j => j < 0 || j >= numVertices);
            f.Info($"adj indices valid: {!invalidAdj}");

            // Connectedness quick check (not full BFS/DFS)
            if(avgDegree < 1)
                f.Info("⚠️  Graph might be disconnected or too sparse.");

            if(noOfParts > numVertices)
                f.Info("⚠️  More partitions than nodes — this may slow down Metis or cause poor results.");

            f.Info("========================================");
        }


        public void ClearMemory() { 
			m_OpMtx = null;
			m_ProlMtx = null;
			m_LeftChangeOfBasis = null;
			m_RightChangeOfBasis = null;
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
		/// <param name="GridDataOnRank0"> if true, the grid data is stored on Rank0 explicitly</param>
		public StandAloneOperatorMappingPairWithGridData(BlockMsrMatrix OperatorMtx, Dictionary<long, long> cellIndexMapping, int[] xadj, int[] adj, int[] NoOfSpecies, bool GridDataOnRank0 = true) {
			m_OpMtx = OperatorMtx;

			Debug.Assert(cellIndexMapping.Count == OperatorMtx._RowPartitioning.TotalNoOfBlocks);
			m_Mapping = new PseudoCoordinateMapping(OperatorMtx._RowPartitioning);

			// if the master/leader proc, get the grid/csr data and convert to the new cell indices
			if (OperatorMtx._RowPartitioning.MpiRank == 0 && GridDataOnRank0) {
				m_xadj = new int[xadj.Length];
				m_adj = new int[adj.Length];
				m_NoOfSpecies = new int[NoOfSpecies.Length];
				RemapCSRAndSpieces(xadj, adj, NoOfSpecies, cellIndexMapping, out m_xadj, out m_adj, out m_NoOfSpecies);
			}

			if (GridDataOnRank0)
				CellToDOFdata = CellIndexToDOFs(m_Mapping);
		}

		/// <summary>
		/// Get block info from 
		/// </summary>
		/// <param name="map"></param>
		/// <returns></returns>
		(long i0Cell, int lenCell)[] CellIndexToDOFs(IBlockPartitioning map) {
			int J = map.LocalNoOfBlocks;
			var myList = new BlockInfo[J]; //for cell currently on this proc (may be distributed to another or not)
			
			List<BlockInfo> globList = new List<BlockInfo>();

			// Create matrix entry information from cell indices
			for (int jCellLoc = 0; jCellLoc < J; jCellLoc++) {
				long jCellGlob = map.FirstBlock + jCellLoc;
				myList[jCellLoc] = new BlockInfo(jCellGlob, map.GetBlockI0(jCellGlob), map.GetBlockLen(jCellGlob));
			}

			Dictionary<int, BlockInfo[]> sendList = new Dictionary<int, BlockInfo[]>();

			if (map.MpiRank != 0)
				sendList.Add(0, myList); //should be sent to 0
			else
				globList.AddRange(myList); //directly add the glob list on rank0

			var data = ArrayMessenger<BlockInfo>.ExchangeData(sendList, map.MPI_Comm);

			if (map.MpiRank == 0) {
				foreach (var message in data) {
					globList.AddRange(message.Value);
				}
				globList.Sort((x, y) => x.iBlock.CompareTo(y.iBlock));
				return globList.Select(b => (b.i0Cell, b.lenCell)).ToArray();
			} 

			return null;			
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
                                Dictionary<long, long> cellIndexMapping,
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

    public class PermutateAndDistribute {

        /// <summary>
        /// - 1st index: MPI rank of target processor
        /// - 2nd index:
        /// - content: index into local vector before permutation, which is to be sent to target processor
        /// </summary>
        int[][] IdxsToBeSent;


        /// <summary>
        /// - index: MPI rank of origin process
        /// - content: number of items received from origin process
        /// </summary>
        int[] PackagesSize;


        ArrayMessenger<double> Permutate; // Source To Target   
        ArrayMessenger<double> BackPermutate; // Target To Source

        // for the permutation from source to target partitioning
        long[] TargetIndices; //indices of cells to be on this processor, according to target partitioning
        long[] SourceIndices; //indices of cells to be on this processor, according to source partitioning

        public int LengthBeforePermutation;
        public int LengthAfterPermutation;

        /// <summary>
        /// Vector permutation class
        /// </summary>
        /// <param name="sourcePartitioning"></param>
        /// <param name="targetPartitioning"></param>
        /// <param name="blockData"></param>
        /// <exception cref="ArgumentNullException"></exception>
        /// <exception cref="ArgumentException"></exception>
        internal PermutateAndDistribute(IBlockPartitioning sourcePartitioning, IBlockPartitioning targetPartitioning, (long i0Cell, int lenCell)[] blockData = null) {
            if(sourcePartitioning == null || targetPartitioning == null)
                throw new ArgumentNullException("Source and target partitionings cannot be null.");

            if(sourcePartitioning.MPI_Comm != targetPartitioning.MPI_Comm)
                throw new ArgumentException("Source and target partitionings must have the same Communicator.");

            if(sourcePartitioning.TotalNoOfBlocks != targetPartitioning.TotalNoOfBlocks)
                throw new ArgumentException("Source and target partitionings must have the same number of blocks.");


            (SourceIndices, TargetIndices) = GetIndxes(targetPartitioning, blockData);
            externalSourceIdx = new List<int>[sourcePartitioning.MpiSize];
            externalTargetIdx = new List<int>[sourcePartitioning.MpiSize];

            CalculateSendAndReceiveLists(SourceIndices, sourcePartitioning, TargetIndices, targetPartitioning);    
        }

        // Cells which are both in source and target partitioning on this processor
        List<int> localSourceIdx = new List<int>();
        List<int> localTargetIdx = new List<int>();

        List<int>[] externalSourceIdx; //Receive list to other processors (local indices on their respective partition)
        List<int>[] externalTargetIdx; //Receive list from other processors (local indices on the target partitioning)

        void CalculateSendAndReceiveLists(long[] SourceIndices, IBlockPartitioning sourcePartitioning, long[] TargetIndices, IBlockPartitioning targetPartitioning) {

            Debug.Assert(SourceIndices.Length == TargetIndices.Length, "Source and target indices must have the same length.");

            for(int i = 0; i < SourceIndices.Length; i++) {
                long srcGlobIdx = SourceIndices[i];
                long trgGlobIdx = TargetIndices[i]; 

                if(sourcePartitioning.IsInLocalRange(srcGlobIdx)) {
                    localSourceIdx.Add(sourcePartitioning.Global2Local(srcGlobIdx));
                    localTargetIdx.Add(targetPartitioning.Global2Local(trgGlobIdx));
                } else {
                    int originRank = sourcePartitioning.FindProcess(srcGlobIdx); // data should be received from the `originRank`
                    // add index of source data to the list for the origin rank
                    (externalSourceIdx[originRank] ??= new List<int>()).Add(checked((int)(srcGlobIdx - sourcePartitioning.GetI0Offest(originRank))));
                    (externalTargetIdx[originRank] ??= new List<int>()).Add(targetPartitioning.Global2Local(trgGlobIdx));
                }
            }

            PackagesSize = externalSourceIdx.Select(l => l?.Count ?? -271897).ToArray();
            if(PackagesSize.Any(sz => sz == 0))
                throw new Exception("internal error");


            using(var f = new ArrayMessenger<int>(sourcePartitioning.MPI_Comm)) {
                for(int r = 0; r < sourcePartitioning.MpiSize; r++) {
                    if(externalSourceIdx[r] != null)
                        f.SetCommPath(r, externalSourceIdx[r].Count);
                }
                f.CommitCommPaths();
                for(int r = 0; r < sourcePartitioning.MpiSize; r++) {
                    if(externalSourceIdx[r] != null)
                        f.Transmit(r, externalSourceIdx[r].ToArray());
                }
                IdxsToBeSent = new int[sourcePartitioning.MpiSize][];
                while(f.GetNext(out int p, out int[] data)) {
                    IdxsToBeSent[p] = data;
                    Debug.Assert(IdxsToBeSent[p].Min() >= 0);
                    Debug.Assert(IdxsToBeSent[p].Max() < sourcePartitioning.LocalLength);
                }
            }

            {
                Permutate = new ArrayMessenger<double>(sourcePartitioning.MPI_Comm);
                for(int r = 0; r < sourcePartitioning.MpiSize; r++) {
                    if(IdxsToBeSent[r] != null)
                        Permutate.SetCommPath(r, IdxsToBeSent[r].Length);
                }
                Permutate.CommitCommPaths();
            }

            {
                BackPermutate = new ArrayMessenger<double>(sourcePartitioning.MPI_Comm);
                for(int r = 0; r < sourcePartitioning.MpiSize; r++) {
                    if(externalSourceIdx[r] != null)
                        BackPermutate.SetCommPath(r, PackagesSize[r]);
                }
                BackPermutate.CommitCommPaths();
            }

            LengthBeforePermutation = IdxsToBeSent.Sum(arr => arr?.Length ?? 0) + localSourceIdx.Count; // total number of items sent + local items
            LengthAfterPermutation = SourceIndices.Length;
        }

        public double[] PermutateVector<V>(V input) where V : IList<double> {
            using(new FuncTrace()) {
                Debug.Assert(Permutate.Size == IdxsToBeSent.Length);

                double[] output = new double[LengthAfterPermutation];

                // Start transmitting data
                // =======================
                { 
                    for(int iTargetProc = 0; iTargetProc < Permutate.Size; iTargetProc++) {
                        if(IdxsToBeSent[iTargetProc] != null)
                            Permutate.Transmit(iTargetProc, input.GetSubVector<int[], int[], double>(IdxsToBeSent[iTargetProc]));
                    }
                }

                // local data exchange
                // ===================
                { 
                    int L = localSourceIdx.Count; // Loop over local indices
                    for(int l = 0; l < L; l++) {
                        output[localTargetIdx[l]] = input[localSourceIdx[l]];
                    }
                }

                // insert received data
                // ====================
                { 
                    while(Permutate.GetNext(out int OriginProc, out double[] data)) {
                        var eTargetIdxs = externalTargetIdx[OriginProc];

                        Debug.Assert(externalSourceIdx[OriginProc].Count == eTargetIdxs.Count);
                        Debug.Assert(data.Length == PackagesSize[OriginProc]);

                        if(eTargetIdxs != null) {
                            int L = eTargetIdxs.Count;
                            for(int l = 0; l < L; l++) {
                                output[eTargetIdxs[l]] = data[l];
                            }
                        }
                    }
                }

                return output;
            }
        }


        public double[] PermutateVectorBack<V>(V input) where V : IList<double> {
            using(new FuncTrace()) {
                Debug.Assert(BackPermutate.Size == IdxsToBeSent.Length);
                int initialSize = LengthBeforePermutation; 
                double[] output = new double[initialSize];

                // Start transmitting data
                // =======================
                {
                    for(int iTargetProc = 0; iTargetProc < BackPermutate.Size; iTargetProc++) {
                        if(externalTargetIdx[iTargetProc] != null)
                            BackPermutate.Transmit(iTargetProc, input.GetSubVector<int[], int[], double>(externalTargetIdx[iTargetProc].ToArray()));
                    }
                }

                // local data exchange
                // ===================
                {
                    int L = localSourceIdx.Count; // Loop over local indices
                    for(int l = 0; l < L; l++) {
                        output[localSourceIdx[l]] = input[localTargetIdx[l]];
                    }
                }

                // insert received data
                // ====================
                {
                    while(BackPermutate.GetNext(out int OriginProc, out double[] data)) {
                        var eTargetIdxs = IdxsToBeSent[OriginProc];

                        if(eTargetIdxs != null) {
                            int L = eTargetIdxs.Length;
                            Debug.Assert(data.Length == L);

                            for(int l = 0; l < L; l++) {
                                output[eTargetIdxs[l]] = data[l];
                            }
                        }
                    }
                }

                return output;

            }
        }
        /// <summary>
        /// Creates permutation matrix and calculates Row and col indixes of each DOF (not block as in cellColumnMapping)
        /// </summary>
        /// <param name="targetPartitioning"></param>
        /// <param name="locDOFsData">DOFs data with global indices</param>
        /// <returns></returns>
        private (long[] SourceIndices, long[] TargetIndices) GetIndxes(IBlockPartitioning targetPartitioning, (long i0Global, int CellLen)[] locDOFsData) {
            using(new FuncTrace()) {

                List<long> TargetIndices = new List<long>();
                List<long> SourceIndices = new List<long>();
                {
                    int cnt = 0;
                    if(locDOFsData != null && locDOFsData.Length > 0) {
                        int NoCells = locDOFsData.Length;

                        for(int j = 0; j < NoCells; j++) {
                            //int Len = part.lenCell[j];
                            int Len = locDOFsData[j].CellLen;
                            long i0Row = targetPartitioning.i0 + cnt;
                            long i0Col = locDOFsData[j].i0Global; //


                            for(int k = 0; k < Len; k++) {
                                TargetIndices.Add(i0Row + k);
                                SourceIndices.Add(i0Col + k);
                            }

                            cnt += Len;
                        }
                    }
                }

#if DEBUG
        				{
        						if (TargetIndices != null) {

        							int L = TargetIndices.Count;
        							for (int l = 0; l < L; l++) {
        								if (l > 0) {
        									Debug.Assert(TargetIndices[l] > TargetIndices[l - 1], "Error, Row indexing is not strictly increasing for some reason.");
        									//Debug.Assert(ColIndices[l] > ColIndices[l - 1], "Error, Column indexing is not strictly increasing for some reason.");
        								}
        							}

        						}


        				}
#endif

                return (SourceIndices.ToArray(), TargetIndices.ToArray());
            }
        }

    }

    /// <summary>
    /// Task parallel (additive) variant of the <see cref="OrthonormalizationMultigrid"/>. 
    /// The coarse solver and smoothers are defined on different MPI communicators.
    /// Therefore, they are independent of each other.
    /// Useful when one wants to utilize large number of processors as smoothers are not idle while waiting for the coarse correction.
    /// They will be working on their own communicators and loop until coarse correction is done.
    /// This is first implementation, there is a large room for improvement here.
    /// </summary>
    public class TaskParallelOrthoMG : ISolverSmootherTemplate, ISolverWithCallback, IProgrammableTermination, ISubsystemSolver {

		enum TpTaskType {
			All,
			Smoother,
			Coarse
		}

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
                get { return "TaskParallelOrthoMG"; }
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
            public int MinimumNoOfPostSmootherSweeps = 1;
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
		/// When smoother or coarse operations are done before than each other, this option enables that they sweep until the other one is done
		/// When enabled, it possibly increases the convergence speed, but becomes run-time dependant.
		/// So, not deterministic anymore and different runs may give different iteration numbers and results.
		/// </summary>
		bool AdvancedParallelism = true;

        /// <summary>
        /// entry point for the TaskParallelOrthoMG (finest level)
        /// </summary>
        public void Init(MultigridOperator op) {
			if (op.OperatorMatrix.MPI_Comm != csMPI.Raw._COMM.WORLD)
				throw new Exception("Task parallel OrthoMG (finest level) should be initiated with an operator in world communicator");

            var ThisAndCoarserLevels = GetSubOperatorChain(op);
			var NoOfProcs = CalculateProcessorDistribution(ThisAndCoarserLevels);
			int WorldSize = op.Mapping.MpiSize;
			Debug.Assert(NoOfProcs[op.LevelIndex] == WorldSize);

			this.FinerLevelCoarseComm = csMPI.Raw._COMM.WORLD;
			this.FinerLevelCoarseCommRank = op.Mapping.MpiRank;

			// Initiate the finest TP MG operator
			var thisTP = new TaskParallelMGOperator(op.OperatorMatrix, op.GetPrologonationOperator, op.Mapping, 
				NoOfProcs[op.LevelIndex], NoOfProcs.ElementAtOrDefault(op.LevelIndex + 1))
				{ InitLevel = op.LevelIndex };

			TaskParallelMGOperator finerTP = thisTP;

			var op_lv_solver = this.CoarserLevelSolver;
			int maxLevel = NoOfProcs.Length;	
			
			for (int TpLevel = 1; TpLevel < ThisAndCoarserLevels.Length; TpLevel++) {
				var op_lv = ThisAndCoarserLevels[TpLevel];
				var coarserTP = new TaskParallelMGOperator(op_lv.OperatorMatrix, op_lv.GetPrologonationOperator, op_lv.Mapping, 
					NoOfProcs[op_lv.LevelIndex], NoOfProcs.ElementAtOrDefault(op_lv.LevelIndex+1), finerTP, 
					op_lv.LeftChangeOfBasis, op_lv.RightChangeOfBasis, 
					ThisAndCoarserLevels.ElementAtOrDefault(TpLevel + 1) != null);

				finerTP.CoarserLevel = coarserTP;
				finerTP = coarserTP;
			}
			InitImpl(thisTP);
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD); //To keep tracing and summary info separate from the Solve phase
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="root"></param>
        /// <returns></returns>
        MultigridOperator[] GetSubOperatorChain(MultigridOperator root) {
			MultigridOperator op_lv = root;
			ISolverSmootherTemplate op_lv_solver = this;
			List<MultigridOperator> SubChain = new List<MultigridOperator>();
			// find the coarsest level
			while (op_lv != null && op_lv_solver != null) {
				SubChain.Add(op_lv);
				op_lv = op_lv.CoarserLevel;
				op_lv_solver = (op_lv_solver is TaskParallelOrthoMG Tp) ? Tp.CoarserLevelSolver : null;
			}

			return SubChain.ToArray();
		}

		/// <summary>
		/// Calculate number of processors for task parallel distribution
		/// </summary>
		/// <param name="Chain">Chain of MG operators</param>
		/// <returns>An int array for each level with the number of processors (returns for whole MG but keeps zero until task parallel starts)</returns>
		int[] CalculateProcessorDistribution(MultigridOperator[] Chain) {
			using (var tr = new FuncTrace()) {
				tr.InfoToConsole = true;
				int WorldSize = Chain.First().Mapping.MpiSize;

				int StartOfTpLevel = Chain.First().LevelIndex;
				int EndOfTpLevel = Chain.Last().LevelIndex;

				long coarsestDOFs = Chain.Last().Mapping.TotalLength;
				long TotDOFs = Chain.Select(c => c.Mapping.TotalLength).Sum();

				// all the smoother operators should have more or less the same DOFs/proc
				int SmootherDOFsPerProc = (int)((TotDOFs - coarsestDOFs) / WorldSize); 
				int[] NoProcs = new int[EndOfTpLevel + 1]; // use the original level index and keep 0 for unused ones (better for exception catching)

				if (Chain.Length < 2) // at least one for this one and one for coarser/direct solver
					throw new ArgumentException("TaskParallelOrthoMG: MG chain must have at least 2 levels."); 

				// the coarsest level has always 1 processor
				NoProcs[EndOfTpLevel] = 1;

				// start with one finer level than the coarsest
				for (int l = EndOfTpLevel - 1; l >= StartOfTpLevel; l--) {
					var op_lv = Chain[l- StartOfTpLevel];

					int NoOfThisLevelSmootherProcs = Math.Max((int)(op_lv.Mapping.TotalLength / SmootherDOFsPerProc), 1);

					if (op_lv.Mapping.TotalNoOfBlocks < NoOfThisLevelSmootherProcs )
						tr.Warning($"Optimization error: task-parallel MG at level {l} has {NoOfThisLevelSmootherProcs} processors assigned but only {op_lv.Mapping.TotalNoOfBlocks} blocks available – you may be using too many processors.");


					int NoOfThisLevelProcs = NoOfThisLevelSmootherProcs + NoProcs[l + 1];

					if (NoOfThisLevelProcs > WorldSize || StartOfTpLevel == l)
						NoOfThisLevelProcs = WorldSize;

					NoProcs[l] = NoOfThisLevelProcs;
					op_lv = op_lv.FinerLevel;
				}


				tr.Info("\nTask Parallel Multigrid hierarchy (o = smoother, x = coarse/direct)");
				tr.Info($"Levels from {StartOfTpLevel} to {EndOfTpLevel}   |  level {EndOfTpLevel} (coarsest) uses 1 processor\n");
				tr.Info("Level        DOFs        # of procs   smoother  coarse    DOFs/proc");
				tr.Info("-----  -------------   ----------  --------  ------  -------------");

				for (int l = StartOfTpLevel; l <= EndOfTpLevel; l++) {
					var op_lv = Chain[l - StartOfTpLevel];
					int p = NoProcs[l];
					int coarse = (l < EndOfTpLevel) ? NoProcs[l + 1] : 1;
					int smooth = p - coarse;
					double ratio = (double)op_lv.Mapping.TotalLength / p;
					tr.Info(
						$"{l,3}   {op_lv.Mapping.TotalLength,12:N0}   {p,10}     {smooth,8}   {coarse,6}   {ratio,13:F1}"
					);

				}

#if DEBUG
				tr.Info("\nPer-core assignment");
				for (int j = StartOfTpLevel; j < EndOfTpLevel+1; j++) {
					// compute indent = sum of smoother procs on all finer levels
					int indent = 0;
					for (int k = StartOfTpLevel; k < j; k++) {
						int coarseK = (k < EndOfTpLevel) ? NoProcs[k + 1] : 1;
						indent += NoProcs[k] - coarseK;
					}
					// current level’s split
					int coarseJ = (j < EndOfTpLevel) ? NoProcs[j + 1] : 1;
					int smoothJ = NoProcs[j] - coarseJ;
					char leaf = (j < EndOfTpLevel) ? 'x' : 'c';   // ‘x’ = coarse, ‘c’ = coarsest
														   // build and print the bar
					string bar = new string(' ', indent)
							   + new string('o', smoothJ)
							   + new string(leaf, coarseJ);
					tr.Info($"L {j} | {bar}");
				}
#endif
				return NoProcs;
			}
		}

		BlockMsrMatrix ChangeCommunicator(BlockMsrMatrix mtx, (long i0Global, int CellLen)[] localBlocks, MPI_Comm comm, int[] MPIRankMapping ) {
			if (mtx == null)
				return null;

			// Row partitioning must be already in the desired form, otherwise, this is not only changing communicator but distributing the matrix

			// Get the local row block information
			var RowlocalBlocks = GetLocalBlocks(mtx._RowPartitioning); 
			var newRowPartition = GetPartitioning(RowlocalBlocks, comm);

			// Get the local column block information
			var newColPartition = GetPartitioning(localBlocks, comm);
			mtx.ChangeMPICommForColumns(MPIRankMapping, newColPartition); // changes the communicator to the new one (subComm)

			List<long> ColIndices = Enumerable.Range((int)mtx._ColPartitioning.i0, mtx._ColPartitioning.LocalLength).Select(i => (long)i).ToList();
			List<long> RowIndices = Enumerable.Range( (int)mtx._RowPartitioning.i0 , mtx._RowPartitioning.LocalLength).Select(i => (long)i).ToList();


			long[] Tlist1 = default(long[]);
			long[] Tlist2 = default(long[]);

			var ret = new BlockMsrMatrix(newRowPartition, newColPartition);
			mtx.AccSubMatrixTo(1.0, ret, RowIndices, Tlist1, ColIndices, Tlist2);

			return ret; // finally create a new matrix with the new Comm
		}

		(long i0Global, int CellLen)[] GetLocalBlocks(IBlockPartitioning part) { 
			int J = part.LocalNoOfBlocks;
			var ret = new (long i0Global, int CellLen)[J];
			long i0Block = part.FirstBlock;
			for (int b = 0; b < J; b++)
				ret[b] = (part.GetBlockI0(i0Block + b), part.GetBlockLen(i0Block+b));

			return ret;
		}

		//BlockMsrMatrix opCommRestrictionOperator = null;
		BlockMsrMatrix worldCommFromCoarseProlongationOperator => TpMapping.CoarserLevel.ProlongationMatrix;

		BlockMsrMatrix worldCommFromCoarseLeftChangeOfBasisMatrix => TpMapping.CoarserLevel.m_LeftChangeOfBasis;
		BlockMsrMatrix worldCommFromCoarseRightChangeOfBasisMatrix => TpMapping.CoarserLevel.m_RightChangeOfBasis;

		BlockMsrMatrix subCommFromCoarseProlongationOperator = null;
		BlockMsrMatrix subCommToCoarseRestrictionOperator = null;

		BlockMsrMatrix subCommLeftChangeOfBasisMatrix = null;
		BlockMsrMatrix subCommRightChangeOfBasisMatrix = null;

		BlockMsrMatrix thisCommOpMatrix = null;
		BlockMsrMatrix subCommSmootherOpMatrix = null;
		BlockMsrMatrix subCommCoarseOpMatrix = null;

        PermutateAndDistribute smootherPermutation = null;
        PermutateAndDistribute coarsePermutation = null;
        //PermutateAndDistribute coarseToSmootherPermutation = null;

        Dictionary<long, long> columnMappingWorldToThis => TpMapping.ThisNewCellMapping;
        Dictionary<long, long> columnMappingWorldToSmoother => TpMapping.SmootherNewCellMapping;
        //  new List<(long source, long target)>(); // cell-based indexes for the mapping from the operator to the smoother matrix (source idx, target idx)
        Dictionary<long, long> columnMappingWorldToCoarse => TpMapping.CoarseNewCellMapping; 
            //= new List<(long source, long target)>();   // cell-based indexes for the mapping from the operator to the coarse matrix (source idx, target idx)

		CoreOrthonormalizationProcedure ortho;

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

		/// <summary>
		/// Initiates the operator, transform all the necessary matrices to the sub communicators
		/// </summary>
		/// <param name="op"></param>
		void InitImpl(TaskParallelMGOperator op) {
			using (var tr = new FuncTrace()) {
                if (object.ReferenceEquals(op, m_OpMapPair))
                    return; // already initialized
				else
                    this.Dispose();

				if (FinerLevelCoarseComm == csMPI.Raw._COMM.NULL)
					return; // this init called by the smoother of the finer, which is not needed.

				this.m_OpMapPair = op;

				NoOfCoarseProcs = op.NoOfCoarseProcs;
				NoOfSmootherProcs = op.NoOfSmootherProcs;
				SplitCommunicator();

				// set operator
				// ============
				MigrateFromWorldToSubComms();
				ortho = new CoreOrthonormalizationProcedure(subCommSmootherOpMatrix);
				CreatePermutations();

#if DEBUG
				TestMatrices(OpMatrix, smootherPermutation, subCommSmootherOpMatrix, $"test_lvl{TpLevel}_s",verbose);
#endif

				// initiate smoother and coarser level (notice initializing them are also task parallel)
				// ======================
				InitSmoothers();
				InitCoarse();
			}
        }

		void InitSmoothers() {
			if (myTask != TpTaskType.All && myTask != TpTaskType.Smoother) return;
			if (Smoothers is null) return;

			var SmootherOpMappingPairOnSubComm = new StandAloneOperatorMappingPairWithGridData(subCommSmootherOpMatrix, columnMappingWorldToSmoother, TpMapping.m_xadj, TpMapping.m_adj, TpMapping.m_NoOfSpecies);

			// init smoothers
			// =============
			foreach (var smoother in Smoothers) { 
				if (smoother != null) {
					if (smoother is SchwarzForTaskParallel schwarzTp) {
						schwarzTp.Init(SmootherOpMappingPairOnSubComm);
					} else if (smoother is ISubsystemSolver ssPreSmother) {
						ssPreSmother.Init(SmootherOpMappingPairOnSubComm);
					} else {
						throw new NotSupportedException($"Unable to initialize pre-smoother if it is not a {typeof(ISubsystemSolver)} and until multigrid operator is designed to work on sub communicators, this won't work");
					}
				}
			}

		}

		double[] Restrict<T1>(T1 IN_fine)
				where T1 : IList<double> {
            using(var tr = new FuncTrace("RestrictFromLvl" + TpLevel)) {
                if(TpMapping.CoarserLevel == null)
                    throw new NotSupportedException("Already on finest level -- no finer level to restrict from.");

                if(IN_fine.Count != subCommToCoarseRestrictionOperator.ColPartition.LocalLength)
                    throw new ArgumentException("Mismatch in length of fine grid vector (input).", "IN_fine");

                double[] OUT_coarse = new double[subCommToCoarseRestrictionOperator.RowPartitioning.LocalLength];

                this.subCommToCoarseRestrictionOperator.SpMV(1.0, IN_fine, 0.0, OUT_coarse);

                if(this.subCommLeftChangeOfBasisMatrix != null) {
                    double[] LB = new double[OUT_coarse.Length];
                    this.subCommLeftChangeOfBasisMatrix.SpMV(1.0, OUT_coarse, 0.0, LB);
                    OUT_coarse.SetV(LB);
                }

                return OUT_coarse;
            }
		}

        double[] Prolongate<T1>(T1 IN_coarse)
        where T1 : IList<double> {
            using(var tr = new FuncTrace("ProlongateBackToLvl"+TpLevel)) {
                if(TpMapping.CoarserLevel == null)
                    throw new NotSupportedException("Already on finest level -- no finer level to prolongate from.");

                if(IN_coarse.Count != subCommFromCoarseProlongationOperator.ColPartition.LocalLength)
                    throw new ArgumentException("Mismatch in length of coarse grid vector (input).", "IN_coarse");

                double[] OUT_fine = new double[subCommFromCoarseProlongationOperator.RowPartitioning.LocalLength];

                if(this.subCommRightChangeOfBasisMatrix != null) {
                    double[] RX = new double[IN_coarse.Count];
                    this.subCommRightChangeOfBasisMatrix.SpMV(1.0, IN_coarse, 0.0, RX);

                    this.subCommFromCoarseProlongationOperator.SpMV(1.0, RX, 0.0, OUT_fine);
                } else {
                    this.subCommFromCoarseProlongationOperator.SpMV(1.0, IN_coarse, 0.0, OUT_fine);
                }

                return OUT_fine;
            }
        }

		ICoordinateMapping MgMap => m_OpMapPair.DgMapping;

		int TpLevel => TpMapping.Level;

		/// <summary>
		/// Changes the communicators fo the operator matrix (to smoother subComm) and the prolongation operator (to this levelComm - to the coarse subComm of the finer level)
		/// (the matrix entries should be already distributed via <see cref="TaskParallelMGOperator"/>
		/// It is not need to deal with the coarse mesh as it will be carried out by the coarser level 
		/// (except the coarsest level <see cref="InitiateCoarsestSolver"/>)
		/// </summary>
		void MigrateFromWorldToSubComms() {
			using (var tr = new FuncTrace()) {
				thisCommOpMatrix = ChangeCommunicator(WorldCommOpMatrix, TpMapping.localBlocksForThisLevel, thisComm, WorldToThisCommMapping);
				subCommSmootherOpMatrix = ChangeCommunicator(TpMapping.m_OpMtx_smoother, TpMapping.localBlocksForSmoother, subComm, WorldToSubCommMapping);

                // In MG operator, the prolongation and restriction operators are stored at the coarse level and restriction performed by operator
                // However, this is not the case here. The MG operator is dedicated for world level communications and key informations.
                // This then must be done by this level solver here for the coarser solver.
                // So each level solver gets data already restricted and then return without prolongating.
                subCommFromCoarseProlongationOperator = ChangeCommunicator(worldCommFromCoarseProlongationOperator, TpMapping.CoarserLevel.localBlocksForThisLevel, subComm, WorldToSubCommMapping); 
                subCommToCoarseRestrictionOperator = subCommFromCoarseProlongationOperator?.Transpose();

				if (worldCommFromCoarseLeftChangeOfBasisMatrix != null)
					subCommLeftChangeOfBasisMatrix = ChangeCommunicator(worldCommFromCoarseLeftChangeOfBasisMatrix,
						TpMapping.CoarserLevel.localBlocksForThisLevel, subComm, WorldToSubCommMapping);

				if (worldCommFromCoarseRightChangeOfBasisMatrix != null)
					subCommRightChangeOfBasisMatrix = ChangeCommunicator(worldCommFromCoarseRightChangeOfBasisMatrix,
					TpMapping.CoarserLevel.localBlocksForThisLevel, subComm, WorldToSubCommMapping);



				TpMapping.ClearMemory();
			}
		}

		(long i0Cell, int lenCell)[] smootherBlocks;
		(long i0Cell, int lenCell)[] coarseBlocks;

		IBlockPartitioning thisPartitioningInThisComm => thisCommOpMatrix._RowPartitioning; //thisCommProlongationOperator == null ? TpMapping.ThisTargetPartitioning : thisCommProlongationOperator._RowPartitioning; // if this instance is for the finest level without prolongation operator, then we can use the this Partitioning at the world level. If not, get the this partitioning from the prolongation operator (comm operator has changed).

		/// <summary>
		/// Permutation matrices from thisComm to subComms (to disribute vectors)
		/// </summary>
		void CreatePermutations() {
            var colMapThisToSmoother = columnMappingWorldToThis
                //.Where(kvp => columnMappingWorldToSmoother.TryGetValue(kvp.Key, out _)) // to avoid key not found exception
                .ToDictionary(
                    kvp => kvp.Value,                             // local index in "this"
                    kvp => columnMappingWorldToSmoother[kvp.Key]  // mapped index in smoother
                );

            var colMapThisToCoarse = columnMappingWorldToThis
                // .Where(kvp => columnMappingWorldToCoarse.TryGetValue(kvp.Key, out _))
                .ToDictionary(
                    kvp => kvp.Value,                           // local index in "this"
                    kvp => columnMappingWorldToCoarse[kvp.Key]  // mapped index in coarse
                );

            //var colMapCoarseToSmoother = columnMappingWorldToCoarse
            //    .ToDictionary(
            //        kvp => kvp.Value,                           // local index in coarse
            //        kvp => columnMappingWorldToSmoother[kvp.Key] // mapped index in smoother
            //    );
            //var SmootherglobalDOFs = CellIndexToDOFs2(smootherTargetPartitioning, TpMapping.localBlocksForSmoother);
            //coarseToSmootherPermutation = new PermutateAndDistribute(coarseTargetPartitioning, smootherTargetPartitioning, smootherBlocks);

            DebugPrint(colMapThisToSmoother, $"lvl_{TpLevel}_colMapThisToSmoother", ".txt");
            DebugPrint(colMapThisToCoarse, $"lvl_{TpLevel}_colMapThisToCoarse", ".txt");

			//var ThisglobalDOFsOld = CellIndexToDOFs(thisPartitioningInThisComm);

			var ThisglobalDOFs = CellIndexToDOFs2(thisPartitioningInThisComm, TpMapping.localBlocksForThisLevel);

			smootherBlocks = GetLocalDistribution(ThisglobalDOFs, colMapThisToSmoother, TpMapping.SmootherCellI0s, 0, NoOfSmootherProcs);
            BlockPartitioning smootherTargetPartitioning = GetPartitioning(smootherBlocks, thisComm);
            smootherPermutation = new PermutateAndDistribute(thisPartitioningInThisComm,smootherTargetPartitioning, smootherBlocks);

			coarseBlocks = GetLocalDistribution(ThisglobalDOFs, colMapThisToCoarse, TpMapping.CoarseCellI0s, NoOfSmootherProcs, NoOfCoarseProcs);
            BlockPartitioning coarseTargetPartitioning = GetPartitioning(coarseBlocks, thisComm);
            coarsePermutation = new PermutateAndDistribute(thisPartitioningInThisComm, coarseTargetPartitioning, coarseBlocks);

#if DEBUG
            TestPermutation(smootherPermutation, $"lvl_{TpLevel}_smootherPermutationTest_", verbose);
            TestPermutation(coarsePermutation, $"lvl_{TpLevel}_coarsePermutationTest_", verbose);
#endif
        }

        (long i0Cell, int lenCell)[] CellIndexToDOFs(IBlockPartitioning map) {
			int J = map.LocalNoOfBlocks;
			var myList = new (long i0Cell, int lenCell)[J]; //for cell currently on this proc (may be distributed to another or not)

			// Create matrix entry information from cell indices
			for (int jCellLoc = 0; jCellLoc < J; jCellLoc++) {
				long jCellGlob = map.FirstBlock + jCellLoc;
				myList[jCellLoc] = (map.GetBlockI0(jCellGlob), map.GetBlockLen(jCellGlob));
			}
			(long i0Cell, int lenCell)[][] globList = myList.MPIAllGatherO(map.MPI_Comm);
			var flatGlobArray = globList.SelectMany(x => x).ToArray();

			return flatGlobArray;
		}

		(long i0Cell, int lenCell)[] CellIndexToDOFs2(IBlockPartitioning map, (long i0Cell, int lenCell)[] worldLevelOriginalBlocks) {
			int B = worldLevelOriginalBlocks.Length;
			var myList = new (long i0Cell, int lenCell)[worldLevelOriginalBlocks.Length]; //for cell currently on this proc 

			long cnt = map.GetBlockI0(map.FirstBlock);
			for (int b = 0; b < B; b++) {
				var block = worldLevelOriginalBlocks[b];
				var updatedBlock = (cnt, block.lenCell);
				myList[b] = updatedBlock;
				cnt += block.lenCell;
			}

			(long i0Cell, int lenCell)[][] globList = myList.MPIAllGatherO(map.MPI_Comm); //[Toprak] ToDo: this should be improved
			var flatGlobArray = globList.SelectMany(x => x).ToArray();

			return flatGlobArray;
		}

        //      List<(long source, long target)> columnMappingOpToSmootherOld;
        //List<(long source, long target)> columnMappingOpToCoarseOld;

        /// <summary>
        /// tests if the permutation matrix is correct with respect to the original matrix and the subComm matrix
        /// </summary>
        /// <param name="opCommMtx">a matrix defined on opComm</param>
        /// <param name="Permutation">permutation to be tested</param>
        /// <param name="subCommMtx">the subComm version of <paramref name="opCommMtx"/> </param>
        /// <param name="tag">tag for i/o outputs</param>
        /// <param name="WriteToFiles">true or false</param>
        /// <exception cref="Exception"></exception>
        void TestMatrices(BlockMsrMatrix opCommMtx, PermutateAndDistribute Permutation, BlockMsrMatrix subCommMtx, string tag = "t_", bool WriteToFiles = false) {
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
			var TestVectorSub = Permutation.PermutateVector(TestVector);

			// result at the subComm
			var subOpResult = new double[subCommMtx._RowPartitioning.LocalLength];
			subCommMtx.SpMV(1.0, TestVectorSub, 0.0, subOpResult);

            // map back the subComm result to the opComm
            var subOpResultBackToOp = Permutation.PermutateVectorBack(subOpResult);

			if (WriteToFiles) {
				//dist.SaveToTextFileDebug($"{tag}localDistForSmoother.txt");
				//GlobDOFidx.SaveToTextFileDebug($"{tag}GlobDOFidx.txt");

				opCommMtx.SaveToTextFileSparseDebug($"{tag}OpMatrix.txt"); //parallel
				opCommMtx.SaveToTextFileSparse($"{tag}OpMatrix.txt"); // combine and write as a single file

				TestVector.SaveToTextFileDebug($"{tag}TestVector", ".txt");
				opResult.SaveToTextFileDebug($"{tag}opResult", ".txt");

				subCommMtx.SaveToTextFileSparseDebug($"{tag}localBlock.txt");
				subCommMtx.SaveToTextFileSparse($"{tag}localBlock.txt");
				
				subOpResult.SaveToTextFileDebug($"{tag}subOpResult", ".txt");
				subOpResultBackToOp.SaveToTextFileDebug($"{tag}backPermutation", ".txt");
			}

			for (int i = 0; i < subOpResultBackToOp.Length; i++) {
				if (Math.Abs(opResult[i] - subOpResultBackToOp[i]) > Math.Pow(10, -8))
					throw new Exception("Something odd with redistribution matrix, the solutions do not hold");
			}
		}

		/// <summary>
		/// tests if the permutation matrix is correct with respect to the original matrix and the subComm matrix
		/// </summary>
		/// <param name="Per">permutation to be tested</param>
		/// <param name="tag">tag for i/o outputs</param>
		/// <param name="WriteToFiles">true or false</param>
		/// <exception cref="Exception"></exception>
		void TestPermutation(PermutateAndDistribute Per, string tag = "t_", bool WriteToFiles = false) {
			//Create a test vector
			Random rnd = new Random(44); //seed is 44
			var TestVector = new double[Per.LengthBeforePermutation];
			for (int i = 0; i < TestVector.Length; i++) {
				TestVector[i] = rnd.NextDouble();
			}

			// test vector to subComm
			var TestVectorSub = Per.PermutateVector(TestVector);

            // map back the subComm result to the opComm
            var backPermutation = Per.PermutateVectorBack(TestVectorSub);


            if (WriteToFiles) {
				TestVector.SaveToTextFileDebug($"{tag}TestVector", ".txt");
                TestVectorSub.SaveToTextFileDebug($"{tag}TestVectorSub", ".txt");
				backPermutation.SaveToTextFileDebug($"{tag}backPermutation", ".txt");
            }

            for (int i = 0; i < backPermutation.Length; i++) {
				if (Math.Abs(TestVector[i] - backPermutation[i]) > Math.Pow(10, -12))
					throw new Exception("Something odd with redistribution matrix, the solutions do not hold");
			}
        }

		(long i0Cell, int lenCell)[] GetLocalDistribution((long i0Cell, int lenCell)[] globalDOFs, Dictionary<long, long> cellColumnMapping, long[] targeti0s, int procOffset, int procSize) {
            int newRank = thisCommRank - procOffset;
            if (newRank < 0 || newRank >= procSize)
                return new (long i0Cell, int lenCell)[0];

            var ret = new List<(long i0Cell, int lenCell)>();
			var newBlocki0 = (int)targeti0s[newRank];
            var newBlockiE = (int)targeti0s[newRank + 1];

            // if cell is designated to be on this processor, add it to the list
            foreach(var kvp in cellColumnMapping) {
                long iSourceCell = kvp.Key;
                long iTargetCell = kvp.Value;

                if (iTargetCell >= newBlocki0 && iTargetCell < newBlockiE) { // if designated cell index falls into the range of this processor
					var sourceBlock = globalDOFs[iSourceCell]; // get the current block index (source)
					ret.Add(sourceBlock);
				}
			}
			//ret.Sort((x, y) => x.i0Cell.CompareTo(y.i0Cell)); // sort the list according to the cell index

            Debug.Assert(newBlockiE - newBlocki0 == ret.Count);
			return ret.ToArray();
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
		int[] ThisCommToSubCommMapping;
		int[] WorldToSubCommMapping;
		int[] WorldToThisCommMapping;

		/// <summary>
		/// Creates a new communicator for the coarse level solver.
		/// </summary>
		void SplitCommunicator() {
			using (var f = new FuncTrace()) {
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

				if (verbose)
					f.Info($"The proc with worldRank-{worldCommRank} with opRank{opCommRank} is assigned to {subCommRank}of{subCommSize} new com{subComm.m1}");

				CreateMpiRankMappings();
			}
		}


		/// <summary>
		/// Creates mapings between MPIRanks on different communicators.
		/// Notice that MPI_Split reorders new communicators w.r.t. old ranks, so this is a basic operation
		/// </summary>
		void CreateMpiRankMappings() {
			Debug.Assert(thisCommsize == TpMapping.NoOfThisProcs);
			int totalCount = thisCommsize; //  number of processes in the communicator

			ThisCommToSubCommMapping = new int[totalCount];
			WorldToSubCommMapping = new int[TpMapping.worldCommSize];
			WorldToSubCommMapping.SetAll(-1);
			for (int i = 0; i < totalCount; i++) {
				ThisCommToSubCommMapping[i] = i < TpMapping.NoOfSmootherProcs ? i : i - TpMapping.NoOfSmootherProcs;
				WorldToSubCommMapping[i + TpMapping.worldMPIOffset] = i < TpMapping.NoOfSmootherProcs ? i : i - TpMapping.NoOfSmootherProcs;
			}

			WorldToThisCommMapping = Enumerable.Range(0, TpMapping.worldCommSize).
				Select(i => i < TpMapping.worldMPIOffset ? -1 : i - TpMapping.worldMPIOffset).ToArray();
		}

		bool verbose = false;

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

						if (Len > 0) {
							i0Cell.Add(cnt);
							LnCell.Add(Len);
							cnt += Len;
						} else { //keep block information to be consistent with global indices but assign them zero length and some suitable i0
							if (j > 0) {
								i0Cell.Add(i0Cell[j - 1]);
							} else {
								i0Cell.Add(0);
							}

							LnCell.Add(Len);
							cnt += Len;
						}
					}
				}

				int LocalLength = cnt;
				var partitioning = new BlockPartitioning(LocalLength, i0Cell, LnCell, comm, i0isLocal: true);
				return partitioning;
			}
		}

        /// <summary>
        /// For the final coarse level solver (probabaly direct solver), we need to actually change the communicator of its matrix at this level
        /// (because the coarsest solver is not a part of the multigrid, it is a stand-alone solver).
        /// </summary>
        /// <exception cref="NotSupportedException">Works only for certain coarse solvers, not with Mg operator</exception>
        private void InitiateCoarsestSolver() {
            using(var tr = new FuncTrace("InitiateCoarsestSolver")) {
                var CoarserTpMapping = (TpMapping.CoarserLevel as TaskParallelMGOperator);
                subCommCoarseOpMatrix = ChangeCommunicator(CoarserTpMapping.OperatorMatrix, CoarserTpMapping.localBlocksForThisLevel, subComm, WorldToSubCommMapping);

                if (verbose) {
                    subCommCoarseOpMatrix.SaveToTextFileSparseDebug($"lvl_{TpLevel}_coarseOpMatrix.txt");
                    subCommCoarseOpMatrix.SaveToTextFileSparse($"lvl_{TpLevel}_coarseOpMatrix.txt");
                }

                if (CoarserLevelSolver is PARDISOSolver PARSolver) {
                    PARSolver.DefineMatrix(subCommCoarseOpMatrix);
                } else if (CoarserLevelSolver is ISubsystemSolver subsystemSolver) {
                    var SmootherOpMappingPairOnSubComm = new StandAloneOperatorMappingPairWithGridData(subCommCoarseOpMatrix, TpMapping.CoarserLevel.CoarseNewCellMapping, TpMapping.m_xadj, TpMapping.m_adj, TpMapping.m_NoOfSpecies, false);
                    subsystemSolver.Init(SmootherOpMappingPairOnSubComm);
                } else {
                    throw new NotSupportedException();
                }

                tr.Info($"MinimumNoOfPostSmootherSweeps was equal to {this.config.MinimumNoOfPostSmootherSweeps} at TPLvl-{TpLevel} but changed to 1 ");
                this.config.MinimumNoOfPostSmootherSweeps = 1;
                CoarserTpMapping.ClearMemory();
            }
        }

		/// <summary>
		/// Initialize the coarse level solver.
		/// </summary>
		/// <exception cref="NotSupportedException"></exception>
		private void InitCoarse() {
			if (myTask != TpTaskType.All && myTask != TpTaskType.Coarse) return;

			using (var tr = new FuncTrace()) {
				if (CoarserLevelSolver == null)
					throw new NotSupportedException("Missing coarse level solver.");

				if (TpMapping.CoarserLevel == null)
					throw new NotSupportedException("Unexpected null CoarserLevel.");

				if (CoarserLevelSolver is TaskParallelOrthoMG ssCoarse) {
					ssCoarse.FinerLevelCoarseComm = myTask != TpTaskType.Smoother ? subComm : csMPI.Raw._COMM.NULL; //[ToprakToDo]: these are not need anymore
					ssCoarse.FinerLevelCoarseCommRank = myTask != TpTaskType.Smoother ? subCommRank : -1;
					ssCoarse.FinerLevelCoarseCommMap = myTask != TpTaskType.Smoother ? thisCommMap : null; //[ToprakToDo]: something might be problematic here. Check the difference between thisCommMap and FinerLevelCoarseCommMap

					ssCoarse.Init(TpMapping.CoarserLevel);
				} else {
					InitiateCoarsestSolver();
				}
			}
		}

		/// <summary>
		/// coarse-level correction; can be defined either
		/// - this option is not supported for this variant (on this level (then the coarse solver may perform its of prolongation/restriction)), or 
		/// - on coarser level, then prolongation/restriction is handled by this solver.
		/// </summary>
		public ISolverSmootherTemplate CoarserLevelSolver;

		/// <summary>
		/// high frequency solver parallel to coarse grid correction
		/// </summary>
		public ISolverSmootherTemplate[] Smoothers;

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
                Array.Copy(B, Res, L);
                OpMatrix.SpMV(-1.0, X, 1.0, Res);
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

		void WriteDebug(int iter, double res, string text) {
			//CurrentTrace.StdoutOnAllRanks();

            int iLevel = TpLevel;			
			if (iLevel >= 0)
				CurrentTrace.Info($"{string.Concat(Enumerable.Repeat("-", iLevel))} OrthoMG, current level={iLevel}, " +
					$"iteration={iter} {(text != null ? " - " + text : "")} and res norm: {res}");

			return;
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
			using (var f = new FuncTrace("TaskParallelOrthoMGSolveLvl" + TpLevel)) {
				CurrentTrace = f;
				f.InfoToConsole = true;

				ThisLevelTime.Start();

				// Initialize vectors
				double[] B = InitializeVector(_B);
				double[] X = InitializeVector(_xl);
				int L = X.Length;
				double[] Res = new double[L];

				// Initialize residual
				double[] Res0 = InitializeResidual(X, B, Res);

                //this.IterationCallback?.Invoke(0, Sol0, Res0, this.m_OpMapPair as MultigridOperator);
                                
                // clear history of coarse solvers
                ortho.Clear();

                PerformIterations(ref X, B, Res, Res0);

				// solution copy
				// =============
				if (!ReferenceEquals(_xl, X)) {
                    _xl.SetV(X);
                }
                ThisLevelTime.Stop();

                if(TpLevel == 0)
                    Console.WriteLine($"Time spent at TpLevel-0 is {ThisLevelTime.Elapsed}");
            } // end of functrace
        }

		FuncTrace CurrentTrace;

		private double[] InitializeVector<T>(T input) where T : IList<double> {
			if (input is double[] array) {
				return array;
			}
			return input.ToArray();
		}

		private double[] InitializeResidual(double[] X, double[] B, double[] Res) {
			int L = X.Length;
			double[] Res0 = new double[L];
			Residual(Res0, X, B); //res0 with x0
			Array.Copy(Res0, Res, L);

            if ((myTask == TpTaskType.Smoother || myTask == TpTaskType.All) && ortho.Norm(Res0) <= 0) {
				double normB = ortho.Norm(B);
				double normX = ortho.Norm(X);
				throw new ArithmeticException($"Residual is 0.0: |X| = {normX}; |B| = {normB}; |Res0| = {ortho.Norm(Res0)}; task={myTask}");
			}

			return Res0;
		}

        private void DebugPrint<T>(T input, string filename, string extension) { 
           if(verbose) {
                if (input is double[] array) {
                    array.SaveToTextFileDebug(filename, extension);
                } else if (input is IList<double> list) {
                    list.ToArray().SaveToTextFileDebug(filename, extension);
                } else {
                    throw new ArgumentException("Input must be of type double[] or IList<double>.");
                }
            }        
        }

        /// <summary>
        /// Distributes the input vectors according to the task type and initializes task-specific data. (must be called on thisComm level)
        /// </summary>
        /// <param name="X"></param>
        /// <param name="B"></param>
        /// <param name="Res"></param>
        /// <param name="XforCoarse"></param>
        /// <param name="XforSmoother"></param>
        /// <param name="ResforCoarse"></param>
        /// <param name="ResforSmoother"></param>
        /// <exception cref="ArgumentNullException"></exception>
        /// <exception cref="InvalidOperationException"></exception>
		private void InitializeTaskSpecificData(double[] X, double[] B, double[] Res, out double[] XforCoarse, out double[] XforSmoother, out double[] ResforCoarse, out double[] ResforSmoother) {
			using (var trace = new FuncTrace("InitializeTaskSpecificData")) {
				// Validate input
				if (X == null || B == null || Res == null) {
					throw new ArgumentNullException("Input vectors X, B, or Res cannot be null.");
				}

                DebugPrint(Res, "Res", ".txt");

                // Permutation operators live on thisComm, which means that they should be called on thisComm
                ResforSmoother = smootherPermutation.PermutateVector(Res);
                XforSmoother = smootherPermutation.PermutateVector(X); // new double[smootherPermutation.LengthAfterPermutation];

                ResforCoarse = coarsePermutation.PermutateVector(Res); 
                XforCoarse = new double[ResforCoarse.Length];  // Coarse correction will be applied to this vector

                Debug.Assert(myTask != TpTaskType.Smoother || ResforCoarse.Length == 0);
                Debug.Assert(myTask != TpTaskType.Coarse || ResforSmoother.Length == 0);

                // Log the initialization
                trace.Info($"Task-specific data initialized for task type: {myTask}");
			}
		}

        /// <summary>
        /// The original version without using optimizations (reference)
        /// </summary>
        /// <param name="X"></param>
        /// <param name="B"></param>
        /// <param name="Res"></param>
        /// <param name="Res0"></param>
		private void PerformIterationsReference(double[] X, double[] B, double[] Res, double[] Res0) {
			int iIter;
			double iter0_resNorm = ortho.Norm(Res0);
			double resNorm = iter0_resNorm;
			var X0 = X.CloneAs();
			for (iIter = 1; ; iIter++) {
				// Check termination
				var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
				if (!termState.bNotTerminate) {
					Converged = termState.bSuccess;
					break;
				}
                ThisLevelIterations++;

                // Initialize task-specific data
                double[] XforCoarse, XforSmoother, ResforCoarse, ResforSmoother;
                InitializeTaskSpecificData(X, B, Res, out XforCoarse, out XforSmoother, out ResforCoarse, out ResforSmoother);

				// Start listening to signal (see AdvancedParallelism)
				ListenSignal(iIter);
				IsEndSignalSent = false;

				// Start timing for coarse grid correction and smoother (as they are parallel and wait for each other in case of AdvancedParallelism)
				CrseLevelTime.Start();
				// from now, it is not serial as it looks on the screen.
				{ 
					// Coarse grid correction
					if (myTask == TpTaskType.All || myTask == TpTaskType.Coarse) {
						double[] XGuess = XforCoarse;
                        XforCoarse = ApplyCoarseGridCorrection(XGuess, ResforCoarse, iIter);
					}

					// smoother
					if (myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {
						double[] XGuess = XforSmoother;
                        XforSmoother = ApplySmoothersReference(XGuess, ResforSmoother, iIter);
					}
				}
				// Stop timing for coarse grid correction
				CrseLevelTime.Stop();

                using(new FuncTrace("PermutateAndMinimize")) {
                    var PreCorr = smootherPermutation.PermutateVectorBack(XforSmoother); //this can be further optimized as smooth part has the full matrix
                    resNorm = ortho.AddSolAndMinimizeResidual(ref PreCorr, X, X0, Res0, Res, "Tp-smoother" + TpLevel);
                    WriteDebug(iIter, resNorm, "Tp-smoother");

                    var CoarseCorrection = Prolongate(XforCoarse);
                    resNorm = ortho.AddSolAndMinimizeResidual(ref CoarseCorrection, X, X0, Res0, Res, "Tp-coarse" + TpLevel);
                    WriteDebug(iIter, resNorm, "Tp-coarse");
                }
				// Iteration callback
				//IterationCallback?.Invoke(iIter, X, Res, m_OpMapPair as MultigridOperator);
			}

            WriteDebug(iIter, resNorm, "final of this level");
		}

        /// <summary>
        /// The main method that performs the multigrid iterations.
        /// The behavior is controlled by the smoother ranks, which store this level matrix and add the solutions to the memory.
        /// This level solution is stored at the Smoother ranks, and the coarse grid correction is calculated at the coarse ranks.
        /// </summary>
        /// <param name="X"></param>
        /// <param name="B"></param>
        /// <param name="Res"></param>
        /// <param name="Res0"></param>
        private void PerformIterations(ref double[] X, double[] B, double[] Res, double[] Res0) {

            // Initialize task-specific data
            double[] XforCoarse, XforSmoother, ResforCoarse, ResforSmoother;
            InitializeTaskSpecificData(X, B, Res, out XforCoarse, out XforSmoother, out ResforCoarse, out ResforSmoother);

            double[] Res0forSmoother = new double[smootherPermutation.LengthAfterPermutation]; 
            Array.Copy(ResforSmoother, Res0forSmoother, smootherPermutation.LengthAfterPermutation);

            double resNorm, iter0_resNorm = -1;
            double[] X0forSmoother = null;
            if(myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {
                iter0_resNorm = ortho.Norm(Res0forSmoother);
                X0forSmoother = XforSmoother.CloneAs();
            }
            iter0_resNorm = iter0_resNorm.MPIBroadcast(0, thisComm); // broadcast the result to all ranks
            resNorm = iter0_resNorm;
            WriteDebug(0, resNorm, $"initial start with MinimumNoOfPostSmootherSweeps={config.MinimumNoOfPostSmootherSweeps}");

            int iIter;
            for(iIter = 1; ; iIter++) {
                using(var f = new FuncTrace("TPLvl" + TpLevel + "Iter" + iIter)) {
                    // Check termination
                    var termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if(!termState.bNotTerminate) {
                        Converged = termState.bSuccess;
                        break;
                    }
                    ThisLevelIterations++;

                    // Start listening to signal (see AdvancedParallelism)
                    ListenSignal(iIter);
                    IsEndSignalSent = false;

                    // Start timing for coarse grid correction and smoother (as they are parallel and wait for each other in case of AdvancedParallelism)
                    CrseLevelTime.Start();
                    // from now, it is not serial as it looks on the screen.
                    {
                        // Coarse grid correction
                        if(myTask == TpTaskType.All || myTask == TpTaskType.Coarse) {
                            XforCoarse = ApplyCoarseGridCorrection(null, ResforCoarse, iIter); //this is coarse correction, so we do not need to pass XGuess
                        }

                        // smoother
                        if(myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {
                            ApplySmoothers(XforSmoother, ResforSmoother, X0forSmoother, Res0forSmoother, iIter); // This will update XforSmoother and ResforSmoother
                        }
                    }
                    // Stop timing for coarse grid correction
                    CrseLevelTime.Stop();

                    using(new FuncTrace("PermutateAndMinimize")) {
                        // [Toprak] ToDo: this not ideal: there can be a better permutation between smoother and coarse, but we do not have it yet
                        var CoarseCorrectionToThisComm = coarsePermutation.PermutateVectorBack(XforCoarse);  
                        var CoarseCorrectionToSmoother = smootherPermutation.PermutateVector(CoarseCorrectionToThisComm); //from this to smoother

                        if(myTask == TpTaskType.All || myTask == TpTaskType.Smoother) {
                            resNorm = ortho.AddSolAndMinimizeResidual(ref CoarseCorrectionToSmoother, XforSmoother, X0forSmoother, Res0forSmoother, ResforSmoother, "Tp-coarse" + TpLevel);
                        }

                        resNorm = resNorm.MPIBroadcast(0, thisComm); // broadcast the result to all ranks
                        WriteDebug(iIter, resNorm, "Tp-coarse");
                    }

                    // Check termination again
                    termState = TerminationCriterion(iIter, iter0_resNorm, resNorm);
                    if(!termState.bNotTerminate) {
                        Converged = termState.bSuccess;
                        break;
                    }

                    using(var tr = new FuncTrace("InitArrays")) {
                        Res = smootherPermutation.PermutateVectorBack(ResforSmoother);
                        ResforCoarse = coarsePermutation.PermutateVector(Res);
                    }
                }
            }

            X = smootherPermutation.PermutateVectorBack(XforSmoother);
            WriteDebug(iIter, resNorm, "final of this level");
        }

        int ThisSmootherGroupLeader => 0;
		int ThisCoarseGroupLeader => NoOfCoarseProcs + NoOfSmootherProcs - 1;

		const int signalTag = 2;

		MPI_Request RecvRequest;
		bool IsEndSignalSent = false;

		static readonly byte[] _signalBuffer = { 1 };
		static readonly GCHandle _signalHandle = GCHandle.Alloc(_signalBuffer, GCHandleType.Pinned);
		static readonly IntPtr _signalPtr = _signalHandle.AddrOfPinnedObject();
		static readonly byte[] _recvBuffer = new byte[1];
		static readonly GCHandle _recvHandle = GCHandle.Alloc(_recvBuffer, GCHandleType.Pinned);
		static readonly IntPtr _recvPtr = _recvHandle.AddrOfPinnedObject();

		int GetTargetRankForSignal() {
			if (thisCommRank == ThisSmootherGroupLeader) {
				return ThisCoarseGroupLeader;
			} else if (thisCommRank == ThisCoarseGroupLeader) {
				return ThisSmootherGroupLeader;
			} else {
				return -1;
			}
		}

		void SendSignal(bool InitiateSignal, int Iter) {
			if (!AdvancedParallelism) return;
			if (!InitiateSignal && IsEndSignalSent) return;

			int targetRank = GetTargetRankForSignal();
			if (targetRank < 0) return;
			int tag = (TpLevel << 16) | Iter;
			csMPI.Raw.Isend(_signalPtr, 1,	csMPI.Raw._DATATYPE.BYTE, targetRank, (signalTag+ tag), thisComm, out MPI_Request req);
			//CurrentTrace.Info($"Sent signal from {thisCommRank} to {targetRank} on {thisComm} with size {thisCommsize} on TpLevel{TpLevel}");

			IsEndSignalSent = true;
			return;
		}

		void ListenSignal(int Iter) {
			if (!AdvancedParallelism) return;
			int targetRank = GetTargetRankForSignal();
			if (targetRank < 0) return;
			int tag = (TpLevel << 16) | Iter;
			csMPI.Raw.Irecv(_recvPtr, 1, csMPI.Raw._DATATYPE.BYTE, targetRank, (signalTag + tag), thisComm, out RecvRequest);

			byte completionSignal = _recvBuffer[0];

			//if (completionSignal == 1)
			//	CurrentTrace.Info($"Got signal on rank {thisCommRank} from {targetRank} on {thisComm} with size {thisCommsize} on TpLevel{TpLevel}");
		}

		bool CheckSignal() {
			if (!AdvancedParallelism) return true; //if not enabled, bypass this feature by returning true
			bool done = false;
			int targetRank = GetTargetRankForSignal();
			if (targetRank > -1)  csMPI.Raw.Test(ref RecvRequest, out done, out MPI_Status status);

			done = done.MPIOr(subComm);
			return done;
		}

        private double[] ApplySmoothers(double[] XforSmoother, double[] ResforSmoother, double[] X0forSmoother, double[] Res0forSmoother, int iIter) {
            using(var trace = new FuncTrace("SmootherTpLvl" + TpLevel)) {
                int i = 0;
                void RunAllSmoothers() {
                    foreach(var smoother in Smoothers) {
                        double[] PreCorrSmoother = new double[smootherPermutation.LengthAfterPermutation]; //XforSmoother;
                        smoother.Solve(PreCorrSmoother, ResforSmoother);

                        var resNorm = ortho.AddSolAndMinimizeResidual(ref PreCorrSmoother, XforSmoother, X0forSmoother, Res0forSmoother, ResforSmoother, "Tp-smoother" + TpLevel);
                        WriteDebug(iIter, resNorm, "Tp-smootherS" + i);
                        i++;
                    }
                }

                RunAllSmoothers(); // Apply smoothers
                for(int sweep = 1; sweep < this.config.MinimumNoOfPostSmootherSweeps; sweep++)
                    RunAllSmoothers();

                if(AdvancedParallelism) {
                    SendSignal(true, iIter);
                    bool done = CheckSignal();
                    int k = 0;
                    while(!done) { //until the coarse solver is done
                        RunAllSmoothers(); // Pre-smoother modifies PreCorr based on Res

                        k++;
                        done = CheckSignal();
                    }
                    CurrentTrace.Info($"{string.Concat(Enumerable.Repeat("-", TpLevel))} OrthoMG, current level={TpLevel}, " +
                    $"iteration={iIter} - All smoothers cycled extra {k}-times while waiting the coarse solver");
                }

                return XforSmoother;
            }
        }

        private double[] ApplySmoothersReference(double[] X, double[] B, int iIter) {
            using (var trace = new FuncTrace("SmootherTpLvl"+TpLevel)) {
                List<double[]> Xsolutions = new List<double[]>();
				void RunAllSmoothers() {
					foreach (var smoother in Smoothers)
						smoother.Solve(X, B);
				}
				
				RunAllSmoothers();                                                              // Apply smoothers (do the first anyway)
				for (int sweep = 1; sweep < this.config.MinimumNoOfPostSmootherSweeps; sweep++) // do rest if required
                    RunAllSmoothers();

				if (AdvancedParallelism) {
					SendSignal(true, iIter);
					bool done = CheckSignal();
					int k = 0;
					while (!done) { //until the coarse solver is done
						RunAllSmoothers(); // Pre-smoother modifies PreCorr based on Res

						k++;
						done = CheckSignal();
					}
					CurrentTrace.Info($"{string.Concat(Enumerable.Repeat("-", TpLevel))} OrthoMG, current level={TpLevel}, " +
					$"iteration={iIter} - All smoothers cycled extra {k}-times while waiting the coarse solver");
				}

				return X;
			}
		}

        private double[] ApplySmoothers2(double[] X, double[] B, int iIter) {
            using(var trace = new FuncTrace("SmootherTpLvl" + TpLevel)) {
                List<double[]> Xsolutions = new List<double[]>();
                void RunAllSmoothers() {
                    foreach(var smoother in Smoothers) {
                        smoother.Solve(X, B);

                    }
                }

                RunAllSmoothers(); // Apply smoothers

                for(int sweep = 1; sweep < this.config.MinimumNoOfPostSmootherSweeps; sweep++)
                    RunAllSmoothers();

                if(AdvancedParallelism) {
                    SendSignal(true, iIter);
                    bool done = CheckSignal();
                    int k = 0;
                    while(!done) { //until the coarse solver is done
                        RunAllSmoothers(); // Pre-smoother modifies PreCorr based on Res

                        k++;
                        done = CheckSignal();
                    }
                    CurrentTrace.Info($"{string.Concat(Enumerable.Repeat("-", TpLevel))} OrthoMG, current level={TpLevel}, " +
                    $"iteration={iIter} - All smoothers cycled extra {k}-times while waiting the coarse solver");
                }

                return X;
            }
        }

        /// <summary>
        /// Solve the coarse grid problem
        /// </summary>
        /// <param name="X">X is both input (initial guess) and output (approximation) (if left null in input, taken 0 array)</param>
        /// <param name="B">RHS</param>
        /// <param name="iIter"></param>
        /// <returns></returns>
        /// <exception cref="InvalidOperationException"></exception>
        private double[] ApplyCoarseGridCorrection (double[] X, double[] B, int iIter) {
			using (var trace = new FuncTrace("CoarseCorrrectionTpLvl"+TpLevel)) {
				if (CoarserLevelSolver == null) {
					throw new InvalidOperationException("Coarse level solver is not initialized.");
				}

                //[Toprak]: this operation is done on subComm level, but it is not necessary to do it on the sub.Instead, it should be done on this level and then permutated.
                //Accordingly, two variants can be used. (discussed in the paper)

                double[] ResCoarse = Restrict(B);
                double[] XCoarse = X is null ? new double[ResCoarse.Length] : Restrict(X); // This is a zero vector, no need for restriction Restrict(X);

                // Solve the coarse grid problem
                CoarserLevelSolver.Solve(XCoarse, ResCoarse);

				if (AdvancedParallelism) {
					SendSignal(true, iIter);

					if (CoarserLevelSolver is TaskParallelOrthoMG TpCoarse) {
						int k = 0;
						bool done = CheckSignal();
						while (!done) {
							CoarserLevelSolver.Solve(XCoarse, ResCoarse);
							k++;
							done = CheckSignal();
						}

                        trace.InfoToConsole = CurrentTrace.InfoToConsole;
                        trace.Info($"{string.Concat(Enumerable.Repeat("-", TpLevel))} OrthoMG, current level={TpLevel}, " +
						$"iteration={iIter} - Coarse cycled extra {k}-times while waiting the coarse solver");
					}
				}

                // Prolongate the correction back to the fine grid
                double[] FineCorrection = Prolongate(XCoarse);

				return FineCorrection;
			}
		}

		/// <summary>
		/// For performance optimization, the <see cref="OrthonormalizationMultigrid"/>
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
        Stopwatch CrseLevelTime = new Stopwatch(); //this includes both smoother + coarse solver time as they are now executed parallel and wait for each other

		/// <summary>
		/// ~
		/// </summary>
		public int IterationsInNested {
            get {
                int iter = 0;

				foreach (var smoother in Smoothers) {
					iter += (smoother?.IterationsInNested ?? 0) + (smoother?.ThisLevelIterations ?? 0);
				}
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

			if (Smoothers != null)
				foreach (var smoother in Smoothers) {
					if (smoother != null) 
						smoother.ResetStat();	
				}

			if (this.CoarserLevelSolver != null)
                this.CoarserLevelSolver.ResetStat();
        }

        /// <summary>
        /// %
        /// </summary>
        public object Clone() {
            throw new NotImplementedException("Clone of " + this.ToString() + " TODO");
        }

        bool m_verbose = false;

        public long UsedMemory() {
            long Memory = 0;
            Memory += MemoryOfMultigrid();
            Memory += MemoryOfSmoother();
            return Memory;
        }

        public long MemoryOfSmoother() {
            long Memory = 0;
            if (this.CoarserLevelSolver is TaskParallelOrthoMG)
                Memory += (this.CoarserLevelSolver as TaskParallelOrthoMG).MemoryOfSmoother();

			if (this.Smoothers != null) {
				foreach (var smoother in Smoothers) {
					if (smoother != null) Memory += smoother.UsedMemory();
				}
			}

            return Memory;
        }

        public long MemoryOfMultigrid() {
            long Memory = 0;
            if (this.CoarserLevelSolver is TaskParallelOrthoMG)
                Memory += (this.CoarserLevelSolver as TaskParallelOrthoMG).MemoryOfMultigrid();

            Memory += ortho?.GetMemUsage() ?? 0;
            return Memory;
        }

        public void Dispose() {
            //if (m_MTracker != null) m_MTracker.Dispose();
            if (m_verbose && m_OpMapPair != null && m_OpMapPair is TaskParallelMGOperator _mgop) {
                int lv = _mgop.Level;
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

			if (Smoothers != null)
				foreach (var smoother in Smoothers)
					smoother?.Dispose();

            //this.PreSmoother = null; // don't delete - we need this again for the next init
            //this.PostSmoother = null;  // don't delete - we need this again for the next init

            
        }

    }
}

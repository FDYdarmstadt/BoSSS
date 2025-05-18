using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Voronoi;
using BoSSS.Solution.AdvancedSolvers.Testing;
using ilPSP;
using ilPSP.Kraypis;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.monkey.CL;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Configuration;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text;
using System.Xml.Schema;
using System.Xml.Serialization;
using static System.Reflection.Metadata.BlobBuilder;


namespace BoSSS.Solution.AdvancedSolvers {


	/// <summary>
	/// Additive Schwarz method modified w.r.t. task parallelization. Similar to <see cref="SchwarzForCoarseMesh"/> but does not require a multigrid operator.
	/// Instead it uses <see cref="StandAloneOperatorMappingPairWithGridData"/> to get the operator matrix and the grid data."/>
	/// </summary>
	public class SchwarzForTaskParallel : ISolverSmootherTemplate {

        /// <summary>
        /// 
        /// </summary>
        [Serializable]
        public class Config : ISolverFactory {

            public string Name => "Additive Schwarz Preconditioner for TaskParallel";

            public string Shortname => "AddSwzTP";

            public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
               throw new NotImplementedException("This solver is not designed to be used with MultigridOperator.");
			}

            public bool Equals(ISolverFactory _other) {
                var other = _other as Config;

                if (other == null)
                    return false;


                if (other.NoOfBlocks != this.NoOfBlocks)
                    return false;
                if (other.Overlap != this.Overlap)
                    return false;
                if (other.EnableOverlapScaling != this.EnableOverlapScaling)
                    return false;

                return true;
            }

            /// <summary>
            /// <see cref="Overlap"/>
            /// </summary>
            private int m_Overlap = 1;

            /// <summary>
            /// Overlap of the Schwarz blocks, in number-of-cells.
            /// - in the case of 0, Additive Schwarz degenerates to Block-Jacobi
            /// - an overlap of 1 is recommended for most cases
            /// - a higher overlap may accelerate convergence, but along the MPI boundaries, overlap will always be limited to 1
            /// - values higher than 2 are currently not supported (its easy to unlock, but there might be no point
            /// </summary>
            public int Overlap {
                get {
                    return m_Overlap;
                }
                set {
                    if (value < 0) {
                        throw new ArgumentException("overlap cannot be negative");
                    }
                    if (value > 5) {
                        throw new ArgumentException($"overlap of {value} is not supported - maximum is 5.");
                    }
                    m_Overlap = value;
                }
            }

            /// <summary>
            /// Number of blocks on all MPI processors
            /// if negative or zero, auto-determined.
            /// </summary>
            public int NoOfBlocks = -1;

            /// <summary>
            /// If <see cref="Overlap"/> > 0, the solution, in cells which are covered by multiple blocks, 
            /// is scaled in the overlapping region by one over the multiplicity.
            /// This option might be useful in some applications but may also fail in others:
            /// - seems to **fail** e.g. if this is used as a preconditioner for PCG (<see cref="SoftPCG"/>)
            /// - improves number of iterations if used e.g. as a smoother for <see cref="OrthonormalizationMultigrid"/>
            /// </summary>
            public bool EnableOverlapScaling = true;
        }

        Config m_config = new Config();

        /// <summary>
        /// Solver configuration
        /// </summary>
        public Config config {
            get {
                return m_config;
            }
        }

        int NoIter = 0;

        /// <summary>
        /// <see cref="ISolverSmootherTemplate.ThisLevelIterations"/>
        /// </summary>
        public int ThisLevelIterations {
            get {
                return this.NoIter;
            }
        }

        /// <summary>
        /// ~
        /// </summary>
        public int IterationsInNested {
            get {
                return 0;
            }
        }

        /// <summary>
        /// <see cref="ISolverSmootherTemplate.Converged"/>
        /// </summary>
        public bool Converged {
            get {
                return m_Converged;
            }
        }

        bool m_Converged = false;

        /// <summary>
        /// <see cref="ISolverSmootherTemplate.ResetStat"/> 
        /// </summary>
        public void ResetStat() {
            this.m_Converged = false;
            this.NoIter = 0;
        }

        public object Clone() {
            throw new NotImplementedException();
        }

        public void Dispose() {
            if (m_BlockSolvers != null) {
                m_BlockSolvers.ForEach(solver => solver?.Dispose());
                m_BlockSolvers = null;
            }
            m_comm?.Dispose();
            m_comm = null;
            m_op = null;
            m_BlockSizes = null;
            m_OverlapScaling = null;
        }

		/// <summary>
		/// Uses METIS to assign a Schwarz block index to each cell on this processor.
		/// Note that the Schwarz block index is **not** the MPI rank, since the number of Schwarz blocks here is typically different than the MPI size
		/// </summary>
		/// <returns>
		/// Mapping from global cell index to Schwarz block
		/// - index: global cell index `j
		/// - content: for cell `j`, the global (among all MPI processors) Schwarz block index to which this cell belongs
		/// Note: only the processor with rank 0 will compute the mappinng, the exchange is done later on.
		/// </returns>

		long[] ComputeSchwarzBlockIndexMETISGlobal(StandAloneOperatorMappingPairWithGridData op) {
			using (var tr = new FuncTrace("TaskParallelSchwarzMetisDist")) {
				(int[] xadj, int[] adjncy) = (op.m_xadj, op.m_adj);
				int[] NoOfSpecies = op.m_NoOfSpecies;

				int[] part;
				if (thisCommRank == 0) {
					int NoOfParts = this.config.NoOfBlocks;
					int J = xadj.Length - 1;

                    if (NoOfParts == 1)
                        return new long[J];

					int ncon = 1;
					int edgecut = 0;
					int[] options = new int[METIS.METIS_NOPTIONS];
					METIS.SETDEFAULTOPTIONS(options);

					options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = 1; // 
					options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 10; // This is the default refinement iterations
					options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 30; // Maximum imbalance of 3 percent (this is the default kway clustering)
					options[(int)METIS.OptionCodes.METIS_OPTION_NUMBERING] = 0;

					part = new int[J];
					Debug.Assert(xadj.Where(idx => idx > adjncy.Length).Count() == 0);
					Debug.Assert(adjncy.Where(j => j >= J).Count() == 0);

                    int[] Weights = NoOfSpecies.Select(i => i * 100 + 1).ToArray(); // avoid zero weights

                    int k = 0;
                    while (k < 3) {
						bool ok = true;

						METIS.PARTGRAPHKWAY(ref J, ref ncon,
							                xadj, adjncy.ToArray(),
							                Weights, null,
							                null, ref NoOfParts,
							                null, null,
							                options, ref edgecut,
							                part);

						for (int p = 0; p < NoOfParts; p++) {
							ok &= part.Contains(p);
                        }
                        if (ok)
                            break;
                        else
                            k++;

                        tr.StdoutOnAllRanks();
                        tr.Warning($"METIS failed to assign all {NoOfParts} parts to the cells. Trying again with different weights.");
					}
				} else {
					part = null;
				}

				return part?.Select(i => (long) i).ToArray();
			}
		}

		/// <summary>
		/// Using CSR format, get the neighbors of a cell. (notice that CSR is stored at the rank 0)
		/// </summary>
		/// <param name="j"></param>
		/// <param name="xadj"></param>
		/// <param name="adj"></param>
		/// <returns></returns>
		long[] getNeighborsCSR(long j, int[] xadj, int[] adj) {
			int start = xadj[j];
			int end = xadj[j + 1];
			long[] neighbors = new long[end - start];
			for (int k = start; k < end; k++) {
				neighbors[k - start] = (long)adj[k];
			}
			return neighbors;
		}

		/// <summary>
		/// Enlarge the Schwarz block by adding the neighbors of the cells in the block.
		/// </summary>
		/// <param name="op"></param>
		/// <param name="block"></param>
		/// <exception cref="ArgumentException"></exception>
		void EnlargeSchwarzBlock(StandAloneOperatorMappingPairWithGridData op, SchwarzBlock block) {
            var Mapping = op.DgMapping;
			using (var tr = new FuncTrace()) {

				var blockCells = block.jGlobal.ToList();

				if (config.Overlap < 0)
					throw new ArgumentException();

				if (config.Overlap > 0) {
					if (config.Overlap > 1 && Mapping.MpiSize > 1) {
						//throw new NotSupportedException("In MPI parallel runs, the maximum supported overlap for the Schwarz preconditioner is 1.");
						tr.Warning("In MPI parallel runs, the overlap for the Schwarz preconditioner is reduced to 1 at MPI boundaries.");
					}

					var bi = block.jLocal;

					// determine overlap regions
					for (int k = 0; k < config.Overlap; k++) { // overlap sweeps
						int Jblock = blockCells.Count;
						for (int j = 0; j < Jblock; j++) { // loop over parts of block
							var jCell = blockCells[j];
							var neighbors = getNeighborsCSR(jCell, op.m_xadj, op.m_adj);
							foreach (var neighbor in neighbors)
								if (!blockCells.Contains(neighbor))
									blockCells.Add(neighbor);
						}
					}

				} else {
					tr.Info("Running Schwarz TP without overlap (at level not known )");
				}

				Debug.Assert(blockCells.Distinct().Count() == blockCells.Count, "Duplicate cells in Schwarz block");

				blockCells.Sort();
				block.jGlobal = blockCells.ToArray();
			}
		}

        SchwarzBlock[][] ExchangeBlocks(SchwarzBlock[] localPartsOfBlocks) {
            using (new FuncTrace()) {
                Debug.Assert(localPartsOfBlocks.Length == m_config.NoOfBlocks);
                int NoOfBlocks = localPartsOfBlocks.Length;
                int MyRank = thisCommRank;
                int MPIsize = thisCommSize;

                List<SchwarzBlock>[] ret = NoOfBlocks.ForLoop(iBlk => new List<SchwarzBlock>());
                foreach (var b in localPartsOfBlocks) {
                    if (b.jGlobal.Length > 0) {
                        b.iOriginRank = MyRank;
                        ret[b.iBlock].Add(b);
                    }
                }

                using (ArrayMessenger<long> messenger = new ArrayMessenger<long>(op_comm)) {

                    // prepare data and setup communication
                    // ====================================

                    List<long>[] commBuffers = MPIsize.ForLoop(i => new List<long>());
                    {
                        foreach (var block in localPartsOfBlocks) {
                            Debug.Assert(block.iOwnerProc >= 0);
                            Debug.Assert(block.iOwnerProc < MPIsize);

                            if (block.iOwnerProc != MyRank && block.jGlobal.Length > 0) {
                                block.Serialize(commBuffers[block.iOwnerProc]);
                            }
                        }
                        Debug.Assert(commBuffers[MyRank].Count == 0);
                        for (int iRank = 0; iRank < MPIsize; iRank++) {
                            if (commBuffers[iRank].Count > 0)
                                messenger.SetCommPath(iRank, commBuffers[iRank].Count);
                        }

                        messenger.CommitCommPaths();
                    }

                    // send data
                    // =========
                    {
                        for (int iRank = 0; iRank < MPIsize; iRank++) {
                            if (commBuffers[iRank].Count > 0)
                                messenger.Transmit(iRank, commBuffers[iRank].ToArray());

                        }
                    }

                    // receive
                    // =======
                    {
                        while (messenger.GetNext(out int rcv_rank, out long[] buffer)) {
                            var rcvBlocks = SchwarzBlock.Deserialze(buffer);
                            foreach (var b in rcvBlocks) {
                                b.iOriginRank = rcv_rank;
                                ret[b.iBlock].Add(b);
                            }
                        }
                    }
                }



                // return
                // ======
                foreach (var list in ret) {
                    list.Sort(new FuncComparer<SchwarzBlock>((SchwarzBlock A, SchwarzBlock B) => A.iOriginRank - B.iOriginRank));
                    Debug.Assert(list.Any(block => block.jGlobal.Length <= 0) == false);
                }
                return ret.Select(list => list.ToArray()).ToArray();
            }
        }

        (long i0Global, int CellLen)[][] SanitizeSchwarzBlocks(SchwarzBlock[][] SchwarzBlocksGlobal) {
            using (new FuncTrace()) {
                int NoOfBlocks = SchwarzBlocksGlobal.Length;
                Debug.Assert(NoOfBlocks == config.NoOfBlocks);

                var ret = new (long i0Global, int CellLen)[NoOfBlocks][];

                for (int iBlk = 0; iBlk < NoOfBlocks; iBlk++) {
                    var SchwarzBlock_iBlk = SchwarzBlocksGlobal[iBlk];
                    if (SchwarzBlock_iBlk != null && SchwarzBlock_iBlk.Length > 0) {
                        int OwnerProcess = SchwarzBlock_iBlk.First().iOwnerProc;

                        if (SchwarzBlock_iBlk.Where(part => part.iOwnerProc != OwnerProcess).Count() > 0) // check that all `SchwarzBlock` in a row have the same owner process
                            throw new ApplicationException("mismatch in owner process ranks");

#if DEBUG
                        for (int i = 0; i < SchwarzBlock_iBlk.Length; i++) {
                            if (i > 0)
                                Debug.Assert(SchwarzBlock_iBlk[i - 1].iOriginRank < SchwarzBlock_iBlk[i].iOriginRank);
                            Debug.Assert(SchwarzBlock_iBlk[i].iOriginRank >= 0);
                            Debug.Assert(SchwarzBlock_iBlk[i].iOriginRank < thisCommSize);
                        }
#endif

                        if (OwnerProcess == thisCommRank) {

                            // concat all blocks
                            // =================

                            long[] jGlobal = new long[0];
                            var helperIdx = new List<(int iPart, int idx)>();
                            int iPart = 0;
                            foreach (var part in SchwarzBlock_iBlk) {
                                jGlobal = jGlobal.Cat(part.jGlobal);
                                for (int i = 0; i < part.jGlobal.Length; i++)
                                    helperIdx.Add((iPart, i));
                                //int NoCells = part.jGlobal.Length;
                                //for (int j = 0; j < NoCells; j++) {
                                //    int Len = part.lenCell[j];
                                //
                                //    i0Cell.Add(cnt);
                                //    LnCell.Add(Len);
                                //    cnt += Len;
                                //}
                                iPart++;
                            }

                            // sort and remove duplicates
                            // ==========================

                            // Note: duplicates might have been added due to Schwarz block overlap into ghost cells.

                            var _hleperIdx = helperIdx.ToArray();
                            Array.Sort(jGlobal, _hleperIdx);

                            var finalBlock = new List<(long i0Global, int CellLen)>();
                            for (int i = 0; i < jGlobal.Length; i++) {
                                var idxPair = _hleperIdx[i];
                                if (finalBlock.Count > 0) {
                                    if (jGlobal[i] == jGlobal[i - 1])
                                        continue; // duplicate cell found
                                    Debug.Assert(SchwarzBlock_iBlk[idxPair.iPart].i0Cell[idxPair.idx] >= finalBlock[finalBlock.Count - 1].i0Global + finalBlock[finalBlock.Count - 1].CellLen, "Cell start indices not in strictly ascending order");
                                }
                                if (SchwarzBlock_iBlk[idxPair.iPart].lenCell[idxPair.idx] <= 0)
                                    continue; // empty cell

                                finalBlock.Add((SchwarzBlock_iBlk[idxPair.iPart].i0Cell[idxPair.idx], SchwarzBlock_iBlk[idxPair.iPart].lenCell[idxPair.idx]));

                            }

                            ret[iBlk] = finalBlock.ToArray();
                        }
                    }
                }

                return ret;
            }

        }

        BlockPartitioning GetPartitioning((long i0Global, int CellLen)[][] SchwarzBlocksGlobal) {
            using (new FuncTrace()) {
                var LnCell = new List<int>();
                var i0Cell = new List<long>();

                int cnt = 0;

                int NoOfBlocks = SchwarzBlocksGlobal.Length;
                Debug.Assert(NoOfBlocks == config.NoOfBlocks);

                for (int iBlk = 0; iBlk < NoOfBlocks; iBlk++) {
                    var SchwarzBlock_iBlk = SchwarzBlocksGlobal[iBlk];
                    if (SchwarzBlock_iBlk != null && SchwarzBlock_iBlk.Length > 0) {

                        {
                            {// foreach (var part in SchwarzBlock_iBlk) {
                                int NoCells = SchwarzBlock_iBlk.Length;// part.jGlobal.Length;
                                for (int j = 0; j < NoCells; j++) {
                                    int Len = SchwarzBlock_iBlk[j].CellLen; // part.lenCell[j];

                                    i0Cell.Add(cnt);
                                    LnCell.Add(Len);
                                    cnt += Len;
                                }
                            }
                        }
                    }
                }

                int LocalLength = cnt;
                var partitioning = new BlockPartitioning(LocalLength, i0Cell, LnCell, op_comm, i0isLocal: true);
                return partitioning;
            }
        }

        /// <summary>
        /// Re-distribution of operator matrix to the Schwarz blocks on the respective processors.
        /// </summary>
        /// <param name="op"></param>
        /// <param name="SchwarzBlocksGlobal"></param>
        /// <returns>
        /// - `Redist`: redistribution matrix; The local Schwarz blocks can than be extracted from the product `Redist*OperatorMatrix`;
        /// - `RowIndices`, `ColIndices`: index lists into the product `Redist*OperatorMatrix`, in order to extract Schwarz blocks,
        ///    (e.g., by <see cref="BlockMsrMatrix.AccSubMatrixTo{V1, V2, V3, V4}(double, IMutableMatrixEx, V1, V2, V3, V4)"/>).
        ///   - 1st index: Schwarz block index `i`; the respective list at `i` is only non-null, if block `i` is beeing owned by the current process
        ///   - 2nd index: Enumeration of matrix rows
        /// </returns>
        private (BlockMsrMatrix Redist, List<long>[] RowIndices, List<long>[] ColIndices) GetRedistributionMatrix(StandAloneOperatorMappingPairWithGridData op, (long i0Global, int CellLen)[][] SchwarzBlocksGlobal) {
            using (new FuncTrace()) {
                int NoOfBlocks = SchwarzBlocksGlobal.Length;
                Debug.Assert(NoOfBlocks == m_config.NoOfBlocks);

                BlockPartitioning TargetPartitioning = GetPartitioning(SchwarzBlocksGlobal);
                BlockMsrMatrix RedistMatrix = new BlockMsrMatrix(TargetPartitioning, OpMtx._RowPartitioning);
                var RowIndices = new List<long>[NoOfBlocks];
                var ColIndices = new List<long>[NoOfBlocks];
                {
                    int cnt = 0;
                    for (int iBlk = 0; iBlk < NoOfBlocks; iBlk++) {
                        (long i0Global, int CellLen)[] SchwarzBlock_iBlk = SchwarzBlocksGlobal[iBlk];
                        if (SchwarzBlock_iBlk != null && SchwarzBlock_iBlk.Length > 0) {

                            int NoCells = SchwarzBlock_iBlk.Length;
                            if (NoCells > 0 && RowIndices[iBlk] == null) {
                                RowIndices[iBlk] = new List<long>();
                                ColIndices[iBlk] = new List<long>();
                            }

                            for (int j = 0; j < NoCells; j++) {
                                //int Len = part.lenCell[j];
                                int Len = SchwarzBlock_iBlk[j].CellLen;
                                long i0Row = TargetPartitioning.i0 + cnt;
                                long i0Col = SchwarzBlock_iBlk[j].i0Global;   //part.i0Cell[j];

                                RedistMatrix.AccBlock(i0Row, i0Col, 1.0, MultidimensionalArray.CreateEye(Len));

                                for (int k = 0; k < Len; k++) {
                                    RowIndices[iBlk].Add(i0Row + k);
                                    ColIndices[iBlk].Add(i0Col + k);
                                }

                                cnt += Len;
                            }
                        }


                    }
                }

#if DEBUG
                {
                    for (int iBlk = 0; iBlk < NoOfBlocks; iBlk++) {
                        if (RowIndices[iBlk] != null) {

                            int L = RowIndices[iBlk].Count;
                            for (int l = 0; l < L; l++) {
                                if (l > 0) {
                                    Debug.Assert(RowIndices[iBlk][l] > RowIndices[iBlk][l - 1], "Error, Row indexing is not strictly increasing for some reason.");
                                    Debug.Assert(ColIndices[iBlk][l] > ColIndices[iBlk][l - 1], "Error, Column indexing is not strictly increasing for some reason.");
                                }
                            }

                        }

                    }
                }
#endif


                return (RedistMatrix, RowIndices, ColIndices);
            }
        }

        PARDISOSolver[] GetBlockSolvers(BlockMsrMatrix Redistributed, List<long>[] RowIndices, List<long>[] ColIndices, int[] metisDOFsPerBlockGlobal) {
            Debug.Assert(RowIndices.Length == ColIndices.Length);
            int NoOfBlocks = RowIndices.Length;
            Debug.Assert(NoOfBlocks == m_config.NoOfBlocks);

            #if DEBUG
            {
                // verify that each block is owned by exactly one process.
                int[] local_OwnedByProc = RowIndices.Select(ary => ary != null ? 1 : 0).ToArray();
                int[] OwnedByProc = local_OwnedByProc.MPISum(Redistributed.MPI_Comm);
                for (int i = 0; i < OwnedByProc.Length; i++) {
                    if (!(OwnedByProc[i] == 1 || metisDOFsPerBlockGlobal[i] == 0)) {
                        throw new ApplicationException($"Block {i} is owned by {OwnedByProc[i]} process(es); (expecting that each block is owned by exactly one processor, if not initially empty).");
                    }
                    Debug.Assert((RowIndices[i] != null) == (ColIndices[i] != null));
                }
            }

            #endif


            PARDISOSolver[] ret = new PARDISOSolver[NoOfBlocks];

            for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {
                if (RowIndices[iBlock] != null) {
                    if (RowIndices[iBlock].Count <= 0)
                        throw new ArgumentException("empty blocks are not allowed");

                    var blockPart = Redistributed._RowPartitioning.GetSubBlocking(RowIndices[iBlock], csMPI.Raw._COMM.SELF, -1);
#if DEBUG
                    Debug.Assert(blockPart.MpiSize == 1);
                    Debug.Assert(blockPart.MpiRank == 0);
                    Debug.Assert(blockPart.LocalLength == RowIndices[iBlock].Count);
                    Debug.Assert(blockPart.i0 == 0);
                    Debug.Assert(blockPart.TotalLength == RowIndices[iBlock].Count);
                    for (int j = 0; j < blockPart.LocalNoOfBlocks; j++) {
                        if (j > 0)
                            Debug.Assert(blockPart.GetBlockI0(j) == blockPart.GetBlockI0(j - 1) + blockPart.GetBlockLen(j - 1));
                        int bt = blockPart.GetBlockType(j);
                        Debug.Assert(blockPart.GetSubblkLen(bt).Sum() == blockPart.GetBlockLen(j));
                    }

                    Debug.Assert(RowIndices[iBlock].Count == RowIndices[iBlock].Count);
#endif
                    var ColIndices_iBlock = ColIndices[iBlock];
                    int L = blockPart.LocalLength;
                    int firstLocal = 0;
                    int lastLocal = 0;
                    for (int i = 0; i < L; i++) {
                        if (i > 0) {
                            Debug.Assert(ColIndices_iBlock[i - 1] < ColIndices_iBlock[i], "column indices must be strictly increasing");
                            Debug.Assert(RowIndices[iBlock][i - 1] < RowIndices[iBlock][i], "row indices must be strictly increasing");
                        }
                        if (ColIndices_iBlock[firstLocal] < Redistributed._ColPartitioning.i0)
                            firstLocal++;
                        if (ColIndices_iBlock[lastLocal] < Redistributed._ColPartitioning.iE)
                            lastLocal++;
                    }

                    if (lastLocal > 0)
                        Debug.Assert(ColIndices_iBlock[lastLocal - 1] < Redistributed._ColPartitioning.iE);
                    if (lastLocal < L)
                        Debug.Assert(ColIndices_iBlock[lastLocal] >= Redistributed._ColPartitioning.iE);
                    if (firstLocal < L)
                        Debug.Assert(ColIndices_iBlock[firstLocal] >= Redistributed._ColPartitioning.i0);
                    if (firstLocal > 0)
                        Debug.Assert(ColIndices_iBlock[firstLocal - 1] < Redistributed._ColPartitioning.i0);





                    var Mtx_iBlk = new BlockMsrMatrix(blockPart);
					Redistributed.AccSubMatrixTo(1.0, Mtx_iBlk, RowIndices[iBlock], default(long[]),
                        ColIndices_iBlock.GetSubVector(firstLocal, lastLocal - firstLocal), (lastLocal - firstLocal).ForLoop(i => (long)i + firstLocal),
                        ColIndices_iBlock.GetSubVector(0, firstLocal).Cat(ColIndices_iBlock.GetSubVector(lastLocal, L - lastLocal)), firstLocal.ForLoop(i => (long)i).Cat((L - lastLocal).ForLoop(i => (long)i + lastLocal))
                        );

                    if (verbose)
                        Mtx_iBlk.SaveToTextFileSparseDebug($"MtxLocBlock{Mapping.TotalNoOfBlocks}_{iBlock}.txt");

					var slv_iBlk = new PARDISOSolver() {
                        CacheFactorization = true,
                        UseDoublePrecision = false,
                        Parallelism = Parallelism.SEQ // hugely important!
                    };
                    slv_iBlk.DefineMatrix(Mtx_iBlk);

                    ret[iBlock] = slv_iBlk;

                }
            }


            return ret;
        }

        /// <summary>
        /// temporary data structure for assembling Schwarz blocks
        /// </summary>
        class SchwarzBlock {
            public int iOriginRank = -21313; // rank on which the block was assembled.

            public int iBlock = -21380934;

            public int iOwnerProc = -123456; // rank, on which the block should be solved (negative init value causes exception if forgotten)

            public List<int> jLocal = new List<int>();

            public long[] jGlobal;

            public long[] i0Cell;

            public int[] lenCell;


            public int Serialize(List<long> stream) {
                int inCount = stream.Count;
                Debug.Assert(jGlobal.Length == lenCell.Length);
                Debug.Assert(lenCell.Length == i0Cell.Length);

                int Len = jGlobal.Length;

                stream.Add(Len);
                stream.Add(iOriginRank);
                stream.Add(iBlock);
                stream.Add(iOwnerProc);

                for (int k = 0; k < Len; k++) {
                    //stream.Add(jLocal[k]); // local cell indices are not required on other processors
                    stream.Add(jGlobal[k]);
                    stream.Add(i0Cell[k]);
                    stream.Add(lenCell[k]);
                }

                return stream.Count - inCount;
            }

            static public IEnumerable<SchwarzBlock> Deserialze(long[] buffer) {
                int cnt = 0;
                List<SchwarzBlock> ret = new List<SchwarzBlock>();
                while (cnt < buffer.Length) {
                    SchwarzBlock b = new SchwarzBlock();
                    ret.Add(b);

                    int Len = checked((int)buffer[cnt]); cnt++;
                    b.iOriginRank = checked((int)buffer[cnt]); cnt++;
                    b.iBlock = checked((int)buffer[cnt]); cnt++;
                    b.iOwnerProc = checked((int)buffer[cnt]); cnt++;

                    b.jGlobal = new long[Len];
                    b.i0Cell = new long[Len];
                    b.lenCell = new int[Len];
                    for (int k = 0; k < Len; k++) {
                        b.jGlobal[k] = buffer[cnt]; cnt++;
                        b.i0Cell[k] = buffer[cnt]; cnt++;
                        b.lenCell[k] = checked((int)buffer[cnt]); cnt++;
                    }

                }
                return ret;
            }
        }

        int GetBlockOwnerRank(int iBlock) {
            Debug.Assert(thisCommSize == m_op.OperatorMatrix._RowPartitioning.MpiSize);
            return (iBlock % thisCommSize);
        }

		StandAloneOperatorMappingPairWithGridData m_op;
		BlockMsrMatrix OpMtx => m_op.OperatorMatrix;
        ICoordinateMapping Mapping => m_op.DgMapping;
        
		MPI_Comm op_comm => Mapping.MPI_Comm;

		public int thisCommRank => Mapping.MpiRank;
		public int thisCommSize => Mapping.MpiSize;

		public void Init(MultigridOperator op) {
			throw new NotSupportedException("Use Init(StandAloneOperatorMappingPairWithGridData op) instead.");
		}

		public void Init(StandAloneOperatorMappingPairWithGridData op) {
#if DEBUG
            InitWithTest(op, true);
#else
            InitWithTest(op, false);
#endif
        }

        bool verbose = false;

		/// <summary>
		/// this step can be performed on world communicator how ever, it should return results on smoother operator communicator
		/// </summary>
		/// <param name="op"></param>
		/// <param name="doTest"></param>
		/// <exception cref="ApplicationException"></exception>
		public void InitWithTest(StandAloneOperatorMappingPairWithGridData op, bool doTest) {
  			using (var f = new FuncTrace()) {
                m_op = op;
                f.StdoutOnAllRanks();

				if (this.config.NoOfBlocks < thisCommSize) 
					f.Warning("!! Warning !! Task parallel Schwarz does not have a block per processor. Either you are using too many cores or there is something wrong.");

				var locBlocks = CalculateBlocks(op); //by utilizing mostly rank0, calculate blocks and get them on respective procs
                var locDOFs = SanitizeSchwarzBlocks(locBlocks);


                int[] DOFsPerBlock;
#if DEBUG
				var DOFsPerBlockLocal = new int[config.NoOfBlocks];
				for (int iBlk = 0; iBlk < locDOFs.Length; iBlk++)
					DOFsPerBlockLocal[iBlk] = locDOFs[iBlk]?.Length ?? 0;

				DOFsPerBlock = DOFsPerBlockLocal.MPISum(op_comm);
#endif

                var RedistAndIndices = GetRedistributionMatrix(op, locDOFs);

				// obtain the part of the matrix which should be solved on this processor via multiplication with the redistribution matrix.
				var LocalBlocks = BlockMsrMatrix.Multiply(RedistAndIndices.Redist, OpMtx);

				// create the local block solvers
				m_BlockSolvers = GetBlockSolvers(LocalBlocks, RedistAndIndices.RowIndices, RedistAndIndices.ColIndices, DOFsPerBlock);
                m_BlockSizes = RedistAndIndices.RowIndices.Select(IDXs => IDXs?.Count ?? -12345).ToArray();

                m_comm = new CommunicationStuff(this, OpMtx._ColPartitioning, RedistAndIndices.ColIndices);

                if (doTest)
                    TestCommunication(RedistAndIndices.Redist);

                if (config.EnableOverlapScaling && config.Overlap > 0) 
                    CalculateScaling();

				if (verbose) {
					RedistAndIndices.Redist.SaveToTextFileSparse("RedistSz" + op.DgMapping.TotalNoOfBlocks + ".txt");
					for (int i = 0; i < RedistAndIndices.ColIndices.Length; i++) {
						RedistAndIndices.RowIndices[i].SaveToTextFileDebug($"RowIndicesSz{op.DgMapping.TotalNoOfBlocks}_{i}", ".txt");
						RedistAndIndices.ColIndices[i].SaveToTextFileDebug($"ColIndicesdistSz{op.DgMapping.TotalNoOfBlocks}_{i}", ".txt");
					}
					LocalBlocks.SaveToTextFileSparse("LocalBlockTogetherSz" + op.DgMapping.TotalNoOfBlocks + ".txt");
				}
			}
        }

		/// <summary>
		/// Calculates the scaling for the overlap region.
		/// </summary>
		void CalculateScaling() {
			int L = Mapping.LocalLength;
			double[] temp = new double[L];
			temp.SetAll(1.0);

			double[][] fwd = AllocBlockMem();
			m_comm.CommRHStoBlocks(temp, fwd);
			temp.Clear();
			m_comm.AccBlockSol(temp, fwd);

			for (int l = 0; l < L; l++)
				temp[l] = 1.0 / temp[l];

			//temp.SaveToTextFile("OverlapScaling.txt");

			m_comm.CommRHStoBlocks(temp, fwd);
			m_OverlapScaling = fwd;
		}

		/// <summary>
		/// Calculates the Schwarz blocks for the given operator 
        /// ad-hoc distribution of dofs to blocks from 0 to others. Can be improved for better scaling between procs.
		/// </summary>
		/// <param name="op"></param>
		/// <returns></returns>
		SchwarzBlock[][] CalculateBlocks(StandAloneOperatorMappingPairWithGridData op) {
			int MPIrnk = thisCommRank;
			var cellToBlock = ComputeSchwarzBlockIndexMETISGlobal(op); //each cell to target block index

			SchwarzBlock[] Blocks = new SchwarzBlock[this.config.NoOfBlocks]; //there is only block info at rank0, the rests are null
			for (int iBlck = 0; iBlck < this.config.NoOfBlocks; iBlck++) {
				Blocks[iBlck] = new SchwarzBlock() { iBlock = iBlck };

				List<long> indices = new List<long>();
                if (cellToBlock != null)
				    for (int i = 0; i < cellToBlock.Length; i++) {
					    if (cellToBlock[i] == iBlck)
						    indices.Add((long)i);
				    }

				Blocks[iBlck].jGlobal = indices.ToArray();
				Blocks[iBlck].i0Cell = new long[] { };
				Blocks[iBlck].lenCell = new int[] { };
			}

			// determine the Owner processor for the respective Schwarz block:
			for (int i = 0; i < Blocks.Length; i++) {
				Blocks[i].iOwnerProc = GetBlockOwnerRank(i);
			}

			if (MPIrnk == 0) { //notice that this is under master rank as it is the only one which has csr graph. so it should also decide on the overlap
				for (int iBlock = 0; iBlock < Blocks.Length; iBlock++) {
					EnlargeSchwarzBlock(op, Blocks[iBlock]); 
					CalculateBlockCellToDOFdata(op, Blocks[iBlock]);
				}
			}

			var locBlocks = ExchangeBlocks(Blocks);
			return locBlocks;
		}

		/// <summary>
		/// Filters/calculates the matrix data for the given block.
		/// ad-hoc distribution of dofs to blocks. Can be improved for better scaling between procs.
        /// but since CSR is gathered already on proc 0. this should not kill the performance.
		/// </summary>
		/// <param name="op"></param>
		/// <param name="Block"></param>
		void CalculateBlockCellToDOFdata(StandAloneOperatorMappingPairWithGridData op, SchwarzBlock Block) {
            Debug.Assert(thisCommRank == 0); //designed only for the master rank as the necessary information is gathered there

			var blockCells = Block.jGlobal;
            Block.i0Cell = new long[blockCells.Length];
            Block.lenCell = new int[blockCells.Length];

            (long i0Cell, int lenCell)[] globalDOFs = op.CellToDOFdata;
			for (int j = 0; j < blockCells.Length; j++) { // loop over parts of block
				var jCell = blockCells[j];
				var blockDOFs = globalDOFs[jCell];
                (Block.i0Cell[j], Block.lenCell[j]) = blockDOFs;
			}
		}

		double[][] m_OverlapScaling = null;

		/// <summary>
		/// Test the communication of the redistribution matrix.
		/// </summary>
		/// <param name="Redist"></param>
		/// <exception cref="ApplicationException"></exception>
		private void TestCommunication(BlockMsrMatrix Redist) {
            using (new FuncTrace()) {
                int NoOfBlocks = this.m_BlockSolvers.Length;

                // create test data
                int LL = m_op.OperatorMatrix.RowPartitioning.LocalLength;
                double[] Origin = new double[LL];
                Origin.FillRandom(seed: thisCommRank);


                // communicate data forward and backward using the Redistribution matrix:
                double[] Dist = new double[Redist.RowPartitioning.LocalLength];
                Redist.SpMV(1.0, Origin, 0.0, Dist);

                double[] Back = new double[LL];
                Redist.Transpose().SpMV(1.0, Dist, 0.0, Back);

                // now, do the same with the Communication vector and compare; from global to blocks...
                var distComm = AllocBlockMem();
                m_comm.CommRHStoBlocks(Origin, distComm);

                double[] distCommCat = new double[0];
                for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {
                    distCommCat = distCommCat.Cat(distComm[iBlock] ?? new double[0]);
                }
                double fwdErr = distCommCat.MPI_L2DistPow2(Dist,op_comm);
                if (fwdErr != 0)
                    throw new ApplicationException("Forward communication error; err = " + fwdErr);


                // ... and back from blocks to global.
                double[] backComm = new double[LL];
                m_comm.AccBlockSol(backComm, distComm);

                double bckErr = Back.MPI_L2DistPow2(backComm, m_op.OperatorMatrix.MPI_Comm);
                if (bckErr != 0)
                    throw new ApplicationException("Backward communication error; err = " + bckErr);
            }
        }

		/// <summary>
		/// Allocates the memory for the blocks.
		/// </summary>
		/// <returns></returns>
		double[][] AllocBlockMem() {
            return m_BlockSizes.Select(sz => sz > 0 ? new double[sz] : null).ToArray();
        }

        int[] m_BlockSizes;

        PARDISOSolver[] m_BlockSolvers;

        CommunicationStuff m_comm;

        class CommunicationStuff : IDisposable {

            readonly SchwarzForTaskParallel m_owner;

            public CommunicationStuff(SchwarzForTaskParallel owner, IBlockPartitioning Part, List<long>[] BlockGlobalIdx) {
                using (new FuncTrace()) {
                    m_owner = owner;
                    int NoOfBlocks = m_owner.m_BlockSolvers.Length;

                    BlockIdxs = new int[NoOfBlocks][];
                    GlobalIdxs = new int[NoOfBlocks][];
                    var _externBlockIdxs = new List<int>[Part.MpiSize][];
                    var _externPackageIdxs = new List<int>[Part.MpiSize][];
                    externBlockIdxs = new int[Part.MpiSize][][];
                    externPackageIdxs = new int[Part.MpiSize][][];
                    var _Packets = new List<int>[Part.MpiSize];
                    PackagesSize = new int[Part.MpiSize];


                    for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) { // loop over Schwarz blocks...
                        if (m_owner.m_BlockSolvers[iBlock] != null) {
                            var locRHSinsertIdxSrc = new List<int>();
                            var locRHSinsertIdxTrg = new List<int>();
                            var glIdxS = BlockGlobalIdx[iBlock];
                            int LneBlock = glIdxS.Count;

                            for (int i = 0; i < LneBlock; i++) { // loop over rows/columns of Schwarz block...
                                long glIdx = glIdxS[i]; // global index which corresponds to block index 'i'

                                if (Part.IsInLocalRange(glIdx)) {
                                    locRHSinsertIdxSrc.Add(Part.Global2Local(glIdx));
                                    locRHSinsertIdxTrg.Add(i);
                                } else {
                                    int originRank = Part.FindProcess(glIdx); // data is received from the `originRank`
                                                                              // 
                                    if (_externBlockIdxs[originRank] == null) {
                                        _externBlockIdxs[originRank] = new List<int>[NoOfBlocks];
                                        _externPackageIdxs[originRank] = new List<int>[NoOfBlocks];
                                    }
                                    Debug.Assert((_externPackageIdxs[originRank] == null) == (_externBlockIdxs[originRank] == null));
                                    if (_externBlockIdxs[originRank][iBlock] == null) {
                                        _externBlockIdxs[originRank][iBlock] = new List<int>();
                                        _externPackageIdxs[originRank][iBlock] = new List<int>();
                                    }
                                    Debug.Assert((_externPackageIdxs[originRank][iBlock] == null) == (_externBlockIdxs[originRank][iBlock] == null));
                                    if (_Packets[originRank] == null) {
                                        _Packets[originRank] = new List<int>();
                                    }

                                    _externBlockIdxs[originRank][iBlock].Add(i);
                                    _externPackageIdxs[originRank][iBlock].Add(_Packets[originRank].Count);

                                    _Packets[originRank].Add(checked((int)(glIdx - Part.GetI0Offest(originRank))));
                                }

                            }


                            GlobalIdxs[iBlock] = locRHSinsertIdxSrc.ToArray();
                            BlockIdxs[iBlock] = locRHSinsertIdxTrg.ToArray();

                            if (GlobalIdxs[iBlock].Length > 0) {
                                Debug.Assert(GlobalIdxs[iBlock].Min() >= 0);
                                Debug.Assert(GlobalIdxs[iBlock].Max() < Part.LocalLength);

                                Debug.Assert(BlockIdxs[iBlock].Min() >= 0);
                                Debug.Assert(BlockIdxs[iBlock].Max() < LneBlock);
                            }
                        }
                    }
                    for (int rnk = 0; rnk < Part.MpiSize; rnk++) {
                        if (_externBlockIdxs[rnk] != null) {
                            externBlockIdxs[rnk] = new int[NoOfBlocks][];
                            externPackageIdxs[rnk] = new int[NoOfBlocks][];
                            for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {

                                if (_externBlockIdxs[rnk][iBlock] != null) {
                                    int LneBlock = BlockGlobalIdx[iBlock].Count;
                                    externBlockIdxs[rnk][iBlock] = _externBlockIdxs[rnk][iBlock].ToArray();
                                    externPackageIdxs[rnk][iBlock] = _externPackageIdxs[rnk][iBlock].ToArray();

                                    Debug.Assert(externBlockIdxs[rnk][iBlock].Length > 0);
                                    Debug.Assert(externPackageIdxs[rnk][iBlock].Length > 0);

                                    Debug.Assert(externBlockIdxs[rnk][iBlock].Min() >= 0);
                                    Debug.Assert(externBlockIdxs[rnk][iBlock].Max() < LneBlock);
                                    Debug.Assert(externPackageIdxs[rnk][iBlock].Min() >= 0);
                                    Debug.Assert(externPackageIdxs[rnk][iBlock].Max() < _Packets[rnk].Count);
                                }
                            }
                        }
                    }

                    PackagesSize = _Packets.Select(l => l?.Count ?? -271897).ToArray();
                    if (PackagesSize.Any(sz => sz == 0))
                        throw new Exception("internal error");


                    using (var f = new ArrayMessenger<int>(Part.MPI_Comm)) {
                        for (int r = 0; r < Part.MpiSize; r++) {
                            if (_Packets[r] != null)
                                f.SetCommPath(r, _Packets[r].Count);
                        }
                        f.CommitCommPaths();
                        for (int r = 0; r < Part.MpiSize; r++) {
                            if (_Packets[r] != null)
                                f.Transmit(r, _Packets[r].ToArray());
                        }
                        Global2BlocksPackages = new int[Part.MpiSize][];
                        while (f.GetNext(out int p, out int[] data)) {
                            Global2BlocksPackages[p] = data;
                            Debug.Assert(Global2BlocksPackages[p].Min() >= 0);
                            Debug.Assert(Global2BlocksPackages[p].Max() < Part.LocalLength);
                        }
                    }

                    {
                        Global2Blocks = new ArrayMessenger<double>(Part.MPI_Comm);
                        for (int r = 0; r < Part.MpiSize; r++) {
                            if (Global2BlocksPackages[r] != null)
                                Global2Blocks.SetCommPath(r, Global2BlocksPackages[r].Length);
                        }
                        Global2Blocks.CommitCommPaths();
                    }

                    {
                        Blocks2Global = new ArrayMessenger<double>(Part.MPI_Comm);
                        for (int r = 0; r < Part.MpiSize; r++) {
                            if (PackagesSize[r] > 0)
                                Blocks2Global.SetCommPath(r, PackagesSize[r]);
                        }
                        Blocks2Global.CommitCommPaths();
                    }
                }
            }




            /// <summary>
            /// Data copy indices into Schwarz blocks, for each Schwarz block, for data exchange on this MPI processor
            /// - 1st index: Schwarz block index
            /// - 2nd index: enumeration
            /// </summary>
            int[][] BlockIdxs;

            /// <summary>
            /// Data copy indices into global vectors, for each Schwarz block, for data exchange on this MPI processor
            /// - 1st index: Schwarz block index
            /// - 2nd index: enumeration
            /// </summary>
            int[][] GlobalIdxs;


            /// <summary>
            /// Data copy indices into Schwarz blocks, for each MPI processor, for each Schwarz block, for data exchange between processors
            /// - 1st index: MPI rang of data origin process
            /// - 2nd index: Schwarz block index
            /// - 3rd index: enumeration
            /// </summary>
            int[][][] externBlockIdxs;

            /// <summary>
            /// Data copy indices into data package for origin process, for each MPI processor, for each Schwarz block, for data between processors
            /// - 1st index: MPI rank of origin process
            /// - 2nd index: Schwarz block index
            /// - 3rd index: enumeration
            /// </summary>
            int[][][] externPackageIdxs;

            /// <summary>
            /// - 1st index: MPI rank of target processor
            /// - 2nd index:
            /// - content: index into global vector
            /// </summary>
            int[][] Global2BlocksPackages;


            /// <summary>
            /// - index: MPI rank of origin process
            /// - content: number of items received from origin process
            /// </summary>
            int[] PackagesSize;


            ArrayMessenger<double> Global2Blocks;

            ArrayMessenger<double> Blocks2Global;


            /// <summary>
            /// Distributes the RHS from the global system (<paramref name="globalRHS"/>) to the RHS-vectors for the Schwarz blocks on the respective MPI processors (<paramref name="localRHSs"/>)
            /// </summary>
            /// <param name="globalRHS">
            /// input
            /// </param>
            /// <param name="localRHSs">
            /// output; RHS for the respective Schwarz block
            /// </param>
            public void CommRHStoBlocks<V>(V globalRHS, double[][] localRHSs) where V : IList<double> {
                using (new FuncTrace()) {
                    int NoOfBlocks = m_owner.m_BlockSolvers.Length;
                    Debug.Assert(NoOfBlocks == BlockIdxs.Length);
                    Debug.Assert(NoOfBlocks == GlobalIdxs.Length);

                    Debug.Assert(Global2Blocks.Size == Global2BlocksPackages.Length);

                    // Start transmitting data
                    // =======================

                    for (int iTargetProc = 0; iTargetProc < Global2Blocks.Size; iTargetProc++) {
                        if (Global2BlocksPackages[iTargetProc] != null)
                            Global2Blocks.Transmit(iTargetProc, globalRHS.GetSubVector<int[], int[], double>(Global2BlocksPackages[iTargetProc]));
                    }

                    // local data exchange
                    // ===================

                    for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {
                        Debug.Assert((BlockIdxs[iBlock] == null) == (m_owner.m_BlockSolvers[iBlock] == null));
                        Debug.Assert((GlobalIdxs[iBlock] == null) == (m_owner.m_BlockSolvers[iBlock] == null));

                        if (m_owner.m_BlockSolvers[iBlock] != null) {
                            double[] localRHS = localRHSs[iBlock];
                            var locRHSinsSrc = BlockIdxs[iBlock];
                            var locRHSinsTrg = GlobalIdxs[iBlock];
                            Debug.Assert(locRHSinsSrc.Length == locRHSinsTrg.Length);

                            int L = locRHSinsSrc.Length;
                            for (int l = 0; l < L; l++) {
                                localRHS[locRHSinsSrc[l]] = globalRHS[locRHSinsTrg[l]];
                            }
                        }
                    }

                    // insert received data
                    // ====================

                    while (Global2Blocks.GetNext(out int OriginProc, out double[] data)) {
                        var eBlockIdxs = externBlockIdxs[OriginProc];
                        var ePacketIdxs = externPackageIdxs[OriginProc];

                        Debug.Assert(data.Length == PackagesSize[OriginProc]);

                        for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {
                            Debug.Assert((eBlockIdxs[iBlock] == null) == (ePacketIdxs[iBlock] == null));

                            if (eBlockIdxs[iBlock] != null) {
                                Debug.Assert(m_owner.m_BlockSolvers[iBlock] != null, "received data for a block which is not solved on this processor");
                                double[] localRHS = localRHSs[iBlock];
                                var locRHSinsSrc = eBlockIdxs[iBlock];
                                var locRHSinsTrg = ePacketIdxs[iBlock];
                                Debug.Assert(locRHSinsSrc.Length == locRHSinsTrg.Length);

                                int L = locRHSinsSrc.Length;
                                for (int l = 0; l < L; l++) {
                                    localRHS[locRHSinsSrc[l]] = data[locRHSinsTrg[l]];
                                }
                            }
                        }
                    }
                }
            }


            /// <summary>
            /// Accumulates the Solution from all Schwarz blocks (<paramref name="localSolS"/>) onto the global solution vector <paramref name="globalSol"/>
            /// </summary>
            /// <typeparam name="U"></typeparam>
            /// <param name="globalSol">
            /// output/accumulator
            /// </param>
            /// <param name="localSolS">
            /// input; 
            /// </param>
            /// <remarks>
            /// Note: this is the quasi-inverse operation to <see cref="CommRHStoBlocks{V}(V, double[][])"/>
            /// </remarks>
            public void AccBlockSol<U>(U globalSol, double[][] localSolS) where U : IList<double> {
                using (new FuncTrace()) {
                    int NoOfBlocks = m_owner.m_BlockSolvers.Length;
                    Debug.Assert(NoOfBlocks == BlockIdxs.Length);
                    Debug.Assert(NoOfBlocks == GlobalIdxs.Length);


                    // note: Since this is the inverse, the role of origin and target rank are vertauscht

                    // Start transmitting data
                    // =======================

                    for (int iOriginProc = 0; iOriginProc < Global2Blocks.Size; iOriginProc++) { // loop over MPI ranks
                        var eBlockInsertIdxs = externBlockIdxs[iOriginProc];
                        var ePacketInsertIdxs = externPackageIdxs[iOriginProc];
                        Debug.Assert((eBlockInsertIdxs == null) == (ePacketInsertIdxs == null));

                        if (eBlockInsertIdxs != null) {
                            double[] sendBackBuffer = new double[PackagesSize[iOriginProc]];

                            for (int iBlock = 0; iBlock < eBlockInsertIdxs.Length; iBlock++) {
                                Debug.Assert((eBlockInsertIdxs[iBlock] == null) == (ePacketInsertIdxs[iBlock] == null));
                                if (eBlockInsertIdxs[iBlock] != null) {
                                    Debug.Assert(eBlockInsertIdxs[iBlock].Length == ePacketInsertIdxs[iBlock].Length);

                                    // sendBackBuffer[ePacketInsertIdxs[iBlock]] +=  localSolS[iBlock][eBlockInsertIdxs[iBlock]]
                                    GenericBlas.AccV(sendBackBuffer, 1.0, localSolS[iBlock], ePacketInsertIdxs[iBlock], eBlockInsertIdxs[iBlock]);
                                }
                            }

                            Blocks2Global.Transmit(iOriginProc, sendBackBuffer);
                        }
                    }


                    // local data exchange
                    // ===================

                    for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {
                        Debug.Assert((BlockIdxs[iBlock] == null) == (m_owner.m_BlockSolvers[iBlock] == null));
                        Debug.Assert((GlobalIdxs[iBlock] == null) == (m_owner.m_BlockSolvers[iBlock] == null));

                        if (m_owner.m_BlockSolvers[iBlock] != null) {
                            double[] localSol = localSolS[iBlock];
                            var locRHSinsSrc = BlockIdxs[iBlock];
                            var locRHSinsTrg = GlobalIdxs[iBlock];
                            Debug.Assert(locRHSinsSrc.Length == locRHSinsTrg.Length);


                            GenericBlas.AccV(globalSol, 1.0, localSol, locRHSinsTrg, locRHSinsSrc);
                            //int L = locRHSinsSrc.Length;
                            //for (int l = 0; l < L; l++) {
                            //    localSol[locRHSinsTrg[l]] = globalRHS[locRHSinsSrc[l]];
                            //}
                        }
                    }

                    // insert received data
                    // ====================

                    while (Blocks2Global.GetNext(out int TargetProc, out double[] packet)) {
                        Debug.Assert(packet.Length == Global2BlocksPackages[TargetProc].Length);

                        // globalSol[Global2BlocksPackages[TargetProc]] += packet;
                        GenericBlas.AccV(globalSol, 1.0, packet, Global2BlocksPackages[TargetProc], default(int[]));
                    }
                }
            }

            public void Dispose() {
                Global2Blocks.Dispose();
                Blocks2Global.Dispose();
            }
        }

        double[][] RHSblocks;
        double[][] Xblocks;

        void CommitInitialXandB<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {

            RHSblocks = AllocBlockMem();
            Xblocks = AllocBlockMem();

            if (X.MPI_L2NormPow2(m_op.OperatorMatrix.MPI_Comm) == 0.0) {
                m_comm.CommRHStoBlocks(B, RHSblocks);
            } else {
                double[] RES = B.ToArray();
                m_op.OperatorMatrix.SpMV(-1.0, X, 1.0, RES);
                m_comm.CommRHStoBlocks(RES, RHSblocks);
            }

        }

		/// <summary>
		/// A special case to test the convergence of the Schwarz method.
		/// </summary>
		/// <typeparam name="U"></typeparam>
		/// <typeparam name="V"></typeparam>
		/// <param name="X"></param>
		/// <param name="B"></param>
		public void TestConvergence<U, V>(U X, V B)
			where U : IList<double>
			where V : IList<double> {
            double[] XFull = X.ToArray().CloneAs();
			double[] BFull = B.ToArray().CloneAs();

			double[] XBlock= X.ToArray().CloneAs();
			double[] BBlock = B.ToArray().CloneAs();


			var FullSolver = new PARDISOSolver();
			FullSolver.DefineMatrix(m_op.OperatorMatrix);
            FullSolver.Solve(XFull, BFull);
            double[] diff = new double[XFull.Length];

			for (int k = 0; k < 100; k++) { 
                this.Solve(XBlock, BBlock);
                diff.Clear();
                diff.AccV(1.0, XFull);
				diff.AccV(-1.0, XBlock);
				Console.WriteLine(k + "X Norm: " + diff.MPI_L2NormPow2(m_op.OperatorMatrix.MPI_Comm));
			}

			Debug.Assert(X.MPI_L2NormPow2(m_op.OperatorMatrix.MPI_Comm) < Math.Pow(10,-2), "Schwarz method does not converge well");
		}

		public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {
            using (new FuncTrace()) {

                CommitInitialXandB(X, B);

                for (int iBlock = 0; iBlock < m_BlockSolvers.Length; iBlock++) {
                    if (m_BlockSolvers[iBlock] != null) {
                        m_BlockSolvers[iBlock].Solve(Xblocks[iBlock], RHSblocks[iBlock]);
      //                  Xblocks[iBlock].SaveToTextFileDebug($"XblocksSz{Mapping.TotalNoOfBlocks}_{iBlock}.txt");
						//RHSblocks[iBlock].SaveToTextFileDebug($"RHSblocksSz{Mapping.TotalNoOfBlocks}_{iBlock}.txt");

						if (m_OverlapScaling != null) {
                            var scale = m_OverlapScaling[iBlock];
							//scale.SaveToTextFileDebug($"scaleSz{Mapping.TotalNoOfBlocks}_{iBlock}.txt");

							var Xb = Xblocks[iBlock];
                            int L = scale.Length;
                            for (int l = 0; l < L; l++)
                                Xb[l] *= scale[l];
                        }
                    }
                }

                m_comm.AccBlockSol(X, Xblocks);

				this.NoIter++;
            }
        }

        public long UsedMemory() {
            long beits = 0;
            if (m_BlockSolvers != null) {
                beits += m_BlockSolvers.Sum(pardiso => pardiso?.UsedMemory() ?? 0);
            }

            return beits;
        }
    }
}

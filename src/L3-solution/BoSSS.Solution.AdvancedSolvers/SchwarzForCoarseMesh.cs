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
using System.Security.Cryptography;
using System.Text;
using System.Xml.Schema;


namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// Additive Schwarz method with optional, multiplicative coarse-grid correction.
    /// 
    /// In this class, we assume to have a relatively low number of DOFs per MPI rank, therefore we can have Schwarz blocks which span more than one process.
    /// So, it is intended to be used at the coarser end of the multigrid structure.
    /// 
    /// In contrast to the <see cref="Schwarz"/> implementation of Schwarz, this implementation:
    /// - computes the blocking globally using <see cref="METIS"/>, and use less blocks than MPI processors
    /// - is hard-coded to use PARDISO in the blocks, i.e., can't use p-multigrid.
    /// </summary>
    public class SchwarzForCoarseMesh : ISolverSmootherTemplate {

        /// <summary>
        /// 
        /// </summary>
        [Serializable]
        public class Config : ISolverFactory {

            public string Name => "Additive Schwarz Preconditioner for Coarse Meshes";

            public string Shortname => "AddSwzCrs";

            public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
                var R = new SchwarzForCoarseMesh();
                R.m_config = this;
                R.Init(level);
                return R;
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
                    if (value > 2) {
                        throw new ArgumentException($"overlap of {value} is not supported - maximum is 2.");
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
            if(m_BlockSolvers != null) {
                m_BlockSolvers.ForEach(solver => solver?.Dispose());
                m_BlockSolvers = null;
            }
            m_comm?.Dispose();
            m_comm = null;
            m_op = null;
            m_BlockSizes = null;
            m_OverlapScaling = null;
        }

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
                for (int j = 0; j<Jloc; j++) {
                    xadj[j] = cnt;
                    int[] Neighbors_j = Neighbors[j];
                    for (int iN = 0; iN<Neighbors_j.Length; iN++) {
                        int jNeigh = Neighbors_j[iN];
                        if (jNeigh<Jloc) {
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
                    xadjOffsets[i] = xadjOffsets[i-1] + locLengths[i-1];
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

        int[] GetNoOfSpeciesList(MultigridOperator op) {

            int J = op.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int[] NoOfSpecies = new int[J];

            XdgAggregationBasis xb = (XdgAggregationBasis)(op.Mapping.AggBasis.FirstOrDefault(b => b is XdgAggregationBasis));

            if (xb != null) {
                for (int jCell = 0; jCell < J; jCell++) {
                    NoOfSpecies[jCell] = xb.GetNoOfSpecies(jCell);
                }
            } else {
                NoOfSpecies.SetAll(1);
            }

            // MPI gather on rank 0
            int MPIsz = op.Mapping.MpiSize;
            int[] rcvCount = MPIsz.ForLoop(r => op.GridData.CellPartitioning.GetLocalLength(r));
            return NoOfSpecies.MPIGatherv(rcvCount, 0, op.OperatorMatrix.MPI_Comm);
        }

        /// <summary>
        /// Uses METIS to assign a Schwarz block index to each cell on this processor.
        /// Note that the Schwarz block index is **not** the MPI rank, since the number of Schwarz blocks here is typically different than theMPI size
        /// </summary>
        /// <returns>
        /// - index: global (among all MPI processors) Schwarz block index
        /// - content: semi-initialized Schwarz block
        /// Note: at first, all processors have all blocks, although the block might be empty on the respective processor.
        /// </returns>
        int[] ComputeSchwarzBlockIndexMETIS(MultigridOperator op) {
            using (new FuncTrace()) {
                (int[] xadj, int[] adjncy) = GetGridGraphForMetis(op);
                int[] NoOfSpecies = GetNoOfSpeciesList(op);

                int MPIrnk = op.Mapping.MpiRank;
                int[] part;
                if (MPIrnk == 0) {
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

                    int NoOfParts = this.config.NoOfBlocks;

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
                } else {
                    part = null;
                }

                // scatter back to processors
                int MPIsz = op.Mapping.MpiSize;
                int[] sendCount = MPIsz.ForLoop(r => op.GridData.CellPartitioning.GetLocalLength(r));
                var partLoc = part.MPIScatterv(sendCount, 0, op.OperatorMatrix.MPI_Comm);
                return partLoc;
            }
        }

        SchwarzBlock[] GetCellsFromPart(int[] part, int NoOfBlocks) {
            var ret = NoOfBlocks.ForLoop(iBlck => new SchwarzBlock() { iBlock = iBlck });

            int L = part.Length;
            for (int j = 0; j < L; j++) {
                ret[part[j]].jLocal.Add(j);
            }

            return ret;
        }

        SchwarzBlock EnlargeSchwarzBlock(SchwarzBlock block, MultigridOperator Mop) {
            using (var tr = new FuncTrace()) {
                int JComp = Mop.Mapping.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                int JGhost = Mop.Mapping.AggGrid.iLogicalCells.NoOfExternalCells;
                int[][] CellNeighbours = Mop.Mapping.AggGrid.iLogicalCells.CellNeighbours;
                
                BitArray marker = new BitArray(JComp + JGhost);

                if (config.Overlap < 0)
                    throw new ArgumentException();
                if (config.Overlap > 0) {
                    if (config.Overlap > 1 && Mop.OperatorMatrix.RowPartitioning.MpiSize > 1) {
                        //throw new NotSupportedException("In MPI parallel runs, the maximum supported overlap for the Schwarz preconditioner is 1.");
                        tr.Warning("In MPI parallel runs, the overlap for the Schwarz preconditioner is reduced to 1 at MPI boundaries.");
                    }

                    var bi = block.jLocal;

                    foreach (int jcomp in bi)
                        marker[jcomp] = true;

                    // determine overlap regions
                    for (int k = 0; k < config.Overlap; k++) { // overlap sweeps
                        int Jblock = bi.Count;
                        for (int j = 0; j < Jblock; j++) { // loop over parts of block
                            int jCell = bi[j];
                            if (jCell < JComp) {


                                int[] Neighs = CellNeighbours[jCell];
                                foreach (int jNeigh in Neighs) {
                                    if (marker[jNeigh] == false) {
                                        // neighbor cell is not already a member of the block
                                        // => add it.
                                        bi.Add(jNeigh);
                                        marker[jNeigh] = true;
                                    }
                                }
                            } else {

                            }
                        }
                    }

                } else {
                    tr.Info("Running Schwarz without overlap (level " + Mop.LevelIndex + ")");
                }

                
                return block;
            }
        }

        /// <summary>
        /// assigns:
        /// - <see cref="SchwarzBlock.jGlobal"/>
        /// - <see cref="SchwarzBlock.i0Cell"/>
        /// - <see cref="SchwarzBlock.lenCell"/>
        /// </summary>
        SchwarzBlock GetBlockPartProperties(SchwarzBlock block, MultigridMapping map) {
            IGridData g = map.AggGrid;
            int NoOfCells = block.jLocal.Count;
            int J = g.iLogicalCells.NoOfLocalUpdatedCells;
            long j0 = g.CellPartitioning.i0;
            long[] gIdxExt = g.iParallel.GlobalIndicesExternalCells;

            long[] jGlobal = new long[NoOfCells];
            block.jGlobal = jGlobal;
            long[] i0Cell = new long[NoOfCells];
            block.i0Cell = i0Cell;
            int[] lenCell = new int[NoOfCells];
            block.lenCell = lenCell;

//            var Temp = new (long i0Cell, int lenCell)[NoOfCells];
            for (int i = 0; i < NoOfCells; i++) {
                int jCellLoc = block.jLocal[i];

                long jCellGlob;
                if (jCellLoc < J)
                    jCellGlob = jCellLoc + j0;
                else
                    jCellGlob = gIdxExt[jCellLoc - J];

                jGlobal[i] = jCellGlob;
                //Temp[i] = (map.GlobalUniqueIndex(0, jCellLoc, 0), map.GetLength(jCellLoc));
                i0Cell[i] = map.GetCellI0(jCellLoc);
                lenCell[i] = map.GetLength(jCellLoc);
            }

            //            // sort according to global cell index
            //            Array.Sort(jGlobal, Temp);
            //            for (int i = 0; i < NoOfCells; i++) {
            //                if(i > 0) {
            //                    Debug.Assert(Temp[i-1].i0Cell < Temp[i].i0Cell);
            //                    Debug.Assert(jGlobal[i-1] < jGlobal[i]);
            //                }
            //#if DEBUG
            //                if (g.CellPartitioning.IsInLocalRange(jGlobal[i])) {
            //                    int jLocal = g.CellPartitioning.TransformIndexToLocal(jGlobal[i]);
            //                    Debug.Assert(Temp[i].i0Cell == map.GlobalUniqueIndex(0, jLocal, 0), "sorting messed things up");
            //                }
            //#endif

            //                i0Cell[i] = Temp[i].i0Cell;
            //                lenCell[i] = Temp[i].lenCell;
            //            }

            block.jLocal = null; // not needed anymore

            return block;
        }


        SchwarzBlock[][] ExchangeBlocks(SchwarzBlock[] localPartsOfBlocks, MultigridOperator op) {
            using (new FuncTrace()) {
                Debug.Assert(localPartsOfBlocks.Length == m_config.NoOfBlocks);
                int NoOfBlocks = localPartsOfBlocks.Length;
                int MyRank = op.Mapping.MpiRank;
                int MPIsize = op.Mapping.MpiSize;

                List<SchwarzBlock>[] ret = NoOfBlocks.ForLoop(iBlk => new List<SchwarzBlock>());
                foreach (var b in localPartsOfBlocks) {
                    if (b.jGlobal.Length > 0) {
                        b.iOriginRank = MyRank;
                        ret[b.iBlock].Add(b);
                    }
                }

                using (ArrayMessenger<long> messenger = new ArrayMessenger<long>(op.OperatorMatrix.MPI_Comm)) {

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

        (long i0Global, int CellLen)[][] SanitizeSchwarzBlocks(SchwarzBlock[][] SchwarzBlocksGlobal, MultigridOperator op) {
            using (new FuncTrace()) {
                int NoOfBlocks = SchwarzBlocksGlobal.Length;
                Debug.Assert(NoOfBlocks == config.NoOfBlocks);

                //var LnCell = new List<int>();
                //var i0Cell = new List<long>();

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
                            Debug.Assert(SchwarzBlock_iBlk[i].iOriginRank < op.Mapping.MpiSize);
                        }
#endif

                        if (OwnerProcess == op.Mapping.MpiRank) {

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
                                    if (jGlobal[i] == jGlobal[i-1])
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

        BlockPartitioning GetPartitioning((long i0Global, int CellLen)[][] SchwarzBlocksGlobal, MultigridOperator op) {
            using (new FuncTrace()) {
                var LnCell = new List<int>();
                var i0Cell = new List<long>();

                int cnt = 0;

                int NoOfBlocks = SchwarzBlocksGlobal.Length;
                Debug.Assert(NoOfBlocks == config.NoOfBlocks);

                for (int iBlk = 0; iBlk < NoOfBlocks; iBlk++) {
                    var SchwarzBlock_iBlk = SchwarzBlocksGlobal[iBlk];
                    if (SchwarzBlock_iBlk != null && SchwarzBlock_iBlk.Length > 0) {
                        //                        int OwnerProcess = SchwarzBlock_iBlk.First().iOwnerProc;

                        //                        if (SchwarzBlock_iBlk.Where(part => part.iOwnerProc != OwnerProcess).Count() > 0) // check that all `SchwarzBlock` in a row have the same owner process
                        //                            throw new ApplicationException("mismatch in owner process ranks");

                        //#if DEBUG
                        //                        for (int i = 0; i < SchwarzBlock_iBlk.Length; i++) {
                        //                            if (i > 0)
                        //                                Debug.Assert(SchwarzBlock_iBlk[i - 1].iOriginRank < SchwarzBlock_iBlk[i].iOriginRank);
                        //                            Debug.Assert(SchwarzBlock_iBlk[i].iOriginRank >= 0);
                        //                            Debug.Assert(SchwarzBlock_iBlk[i].iOriginRank < op.Mapping.MpiSize);
                        //                        }
                        //#endif

                        //                        if (OwnerProcess == op.Mapping.MpiRank) {
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
                var partitioning = new BlockPartitioning(LocalLength, i0Cell, LnCell, op.OperatorMatrix.MPI_Comm, i0isLocal: true);
                return partitioning;
            }
        }


        private (BlockMsrMatrix Redist, List<long>[] RowIndices, List<long>[] ColIndices) GetRedistributionMatrix(MultigridOperator op, (long i0Global, int CellLen)[][] SchwarzBlocksGlobal) {
            using (new FuncTrace()) {
                int NoOfBlocks = SchwarzBlocksGlobal.Length;
                Debug.Assert(NoOfBlocks == m_config.NoOfBlocks);

                BlockPartitioning TargetPartitioning = GetPartitioning(SchwarzBlocksGlobal, op);
                BlockMsrMatrix RedistMatrix = new BlockMsrMatrix(TargetPartitioning, op.OperatorMatrix._RowPartitioning);
                var RowIndices = new List<long>[NoOfBlocks];
                var ColIndices = new List<long>[NoOfBlocks];
                {
                    int cnt = 0;
                    for (int iBlk = 0; iBlk < NoOfBlocks; iBlk++) {
                        var SchwarzBlock_iBlk = SchwarzBlocksGlobal[iBlk];
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

        public PARDISOSolver[] GetBlockSolvers(BlockMsrMatrix Redist, List<long>[] RowIndices, List<long>[] ColIndices) {
            Debug.Assert(RowIndices.Length == ColIndices.Length);
            int NoOfBlocks = RowIndices.Length;
            Debug.Assert(NoOfBlocks == m_config.NoOfBlocks);

#if DEBUG
            int[] OwnedByProc = RowIndices.Select(ary => ary != null ? 1 : 0).ToArray().MPISum(Redist.MPI_Comm);
            for(int i = 0; i < OwnedByProc.Length; i++) {
                Debug.Assert(OwnedByProc[i] == 1);
                Debug.Assert((RowIndices[i] != null) == (ColIndices[i] != null));
            }

#endif


            PARDISOSolver[] ret = new PARDISOSolver[NoOfBlocks];

            for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {
                if (RowIndices[iBlock] != null) {
                    if (RowIndices[iBlock].Count <= 0)
                        throw new ArgumentException("empty blocks are not allowed");

                    var blockPart = Redist._RowPartitioning.GetSubBlocking(RowIndices[iBlock], csMPI.Raw._COMM.SELF, -1);
#if DEBUG
                    Debug.Assert(blockPart.MpiSize == 1);
                    Debug.Assert(blockPart.MpiRank == 0);
                    Debug.Assert(blockPart.LocalLength == RowIndices[iBlock].Count);
                    Debug.Assert(blockPart.i0 == 0);
                    Debug.Assert(blockPart.TotalLength == RowIndices[iBlock].Count);
                    for(int j = 0; j < blockPart.LocalNoOfBlocks; j++) {
                        if(j > 0)
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
                            Debug.Assert(ColIndices_iBlock[i-1] < ColIndices_iBlock[i], "column indices must be strictly increasing");
                            Debug.Assert(RowIndices[iBlock][i-1] < RowIndices[iBlock][i], "row indices must be strictly increasing");
                        }
                        if (ColIndices_iBlock[firstLocal] < Redist._ColPartitioning.i0)
                            firstLocal++;
                        if(ColIndices_iBlock[lastLocal] < Redist._ColPartitioning.iE)
                            lastLocal++;
                    }

                    if(lastLocal > 0)
                        Debug.Assert(ColIndices_iBlock[lastLocal - 1] < Redist._ColPartitioning.iE);
                    if(lastLocal < L)
                        Debug.Assert(ColIndices_iBlock[lastLocal] >= Redist._ColPartitioning.iE);
                    if(firstLocal < L)
                        Debug.Assert(ColIndices_iBlock[firstLocal] >= Redist._ColPartitioning.i0);
                    if (firstLocal > 0)
                        Debug.Assert(ColIndices_iBlock[firstLocal - 1] < Redist._ColPartitioning.i0);
                    




                    var Mtx_iBlk = new BlockMsrMatrix(blockPart);
                    Redist.AccSubMatrixTo(1.0, Mtx_iBlk, RowIndices[iBlock], default(long[]), 
                        ColIndices_iBlock.GetSubVector(firstLocal, lastLocal - firstLocal), (lastLocal - firstLocal).ForLoop(i => (long)i + firstLocal),
                        ColIndices_iBlock.GetSubVector(0, firstLocal).Cat(ColIndices_iBlock.GetSubVector(lastLocal, L - lastLocal)), firstLocal.ForLoop(i => (long)i).Cat((L - lastLocal).ForLoop(i => (long)i + lastLocal))
                        );

                 
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

                for(int k = 0; k < Len; k++) {
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
                while(cnt < buffer.Length) {
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
            int MPIsz = m_op.Mapping.MpiSize;
            return (iBlock % (MPIsz));
            //return (iBlock % (MPIsz - 1)) + 1; // we avoid rank 0, because rank 0 is doing the coarse solve
        }

        MultigridOperator m_op;

        public void Init(MultigridOperator op) {
#if DEBUG
            InitWithTest(op, true);
#else
            InitWithTest(op, false);
#endif
        }

        public void InitWithTest(MultigridOperator op, bool doTest) {
            using (new FuncTrace()) {
                m_op = op;

                
                //
                int[] SwzBlkLocal = ComputeSchwarzBlockIndexMETIS(op); // index: global Schwarz block index
                                                                       // at first, all processors have all blocks, although the block might be empty on the respective processor.
                Debug.Assert(SwzBlkLocal.Min() >= 0);
                Debug.Assert(SwzBlkLocal.Max() < config.NoOfBlocks);

                int NoOfBlocks = config.NoOfBlocks;
                SchwarzBlock[] Blocks = GetCellsFromPart(SwzBlkLocal, NoOfBlocks);

                // overlap
                for (int iBlock = 0; iBlock < Blocks.Length; iBlock++)
                    Blocks[iBlock] = EnlargeSchwarzBlock(Blocks[iBlock], op);

                // determine the Owner processor for the respective Schwarz block:
                for (int i = 0; i < Blocks.Length; i++) {
                    Blocks[i].iOwnerProc = GetBlockOwnerRank(i);
                    Blocks[i] = GetBlockPartProperties(Blocks[i], op.Mapping);
                }

                // perform the data exchange
                var SchwarzBlocksGlobal = SanitizeSchwarzBlocks(ExchangeBlocks(Blocks, op), op);
                
                //  BlockPartitioning(int LocalLength, IEnumerable<long> BlockI0, IEnumerable<int> BlockLen, MPI_Comm MpiComm, bool i0isLocal = false) 
                var RedistAndIndices = GetRedistributionMatrix(op, SchwarzBlocksGlobal);

                // obtain the part of the matrix which should be solved on this processor via multiplication with the redistribution matrix.
                var LocalBlocks = BlockMsrMatrix.Multiply(RedistAndIndices.Redist, op.OperatorMatrix);

                // create the local block solvers
                m_BlockSolvers = GetBlockSolvers(LocalBlocks, RedistAndIndices.RowIndices, RedistAndIndices.ColIndices);
                m_BlockSizes = RedistAndIndices.RowIndices.Select(IDXs => IDXs?.Count ?? -12345).ToArray();

                m_comm = new CommunicationStuff(this, op.OperatorMatrix._ColPartitioning, RedistAndIndices.ColIndices);

                if (doTest)
                    TestCommunication(RedistAndIndices.Redist);

                if(config.EnableOverlapScaling && config.Overlap > 0) {
                    int L = op.Mapping.LocalLength;
                    double[] temp = new double[L];
                    temp.SetAll(1.0);

                    double[][] fwd = AllocBlockMem();
                    m_comm.CommRHStoBlocks(temp, fwd);
                    temp.Clear();
                    m_comm.AccBlockSol(temp, fwd);

                    for(int l = 0; l < L; l++)
                        temp[l] = 1.0/temp[l];

                    //temp.SaveToTextFile("OverlapScaling.txt");

                    m_comm.CommRHStoBlocks(temp, fwd);
                    m_OverlapScaling = fwd;
                }
            }
        }

        double[][] m_OverlapScaling = null;

        private void TestCommunication(BlockMsrMatrix Redist) {
            using (new FuncTrace()) {
                int NoOfBlocks = this.m_BlockSolvers.Length;

                // create test data
                int LL = m_op.OperatorMatrix.RowPartitioning.LocalLength;
                double[] Origin = new double[LL];
                Origin.FillRandom(seed: m_op.DgMapping.MpiRank);

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
                double fwdErr = distCommCat.MPI_L2DistPow2(Dist);
                Console.WriteLine("Forward communication difference = " + fwdErr);
                if (fwdErr != 0)
                    throw new ApplicationException("Forward communication error; err = " + fwdErr);


                // ... and back from blocks to global.
                double[] backComm = new double[LL];
                m_comm.AccBlockSol(backComm, distComm);

                double bckErr = Back.MPI_L2DistPow2(backComm, m_op.OperatorMatrix.MPI_Comm);
                Console.WriteLine("Backward communication difference = " + bckErr);
                if (bckErr != 0)
                    throw new ApplicationException("Backward communication error; err = " + bckErr);
            }
        }


        double[][] AllocBlockMem() {
            return m_BlockSizes.Select(sz => sz > 0 ? new double[sz] : null).ToArray();
        }


        int[] m_BlockSizes; 

        PARDISOSolver[] m_BlockSolvers;

        CommunicationStuff m_comm; 

        class CommunicationStuff : IDisposable {

            readonly SchwarzForCoarseMesh m_owner;

            public CommunicationStuff(SchwarzForCoarseMesh owner, IBlockPartitioning Part, List<long>[] BlockGlobalIdx) {
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


                    using(var f = new ArrayMessenger<int>(Part.MPI_Comm)) {
                        for(int r = 0; r < Part.MpiSize; r++) {
                            if (_Packets[r] != null)
                                f.SetCommPath(r, _Packets[r].Count);
                        }
                        f.CommitCommPaths();
                        for (int r = 0; r < Part.MpiSize; r++) {
                            if (_Packets[r] != null)
                                f.Transmit(r, _Packets[r].ToArray());
                        }
                        Global2BlocksPackages = new int[Part.MpiSize][];
                        while(f.GetNext(out int p, out int[] data)) {
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



        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {
            using (new FuncTrace()) {
                double[][] RHSblocks = AllocBlockMem();
                double[][] Xblocks = AllocBlockMem();

                if(X.MPI_L2NormPow2(m_op.OperatorMatrix.MPI_Comm) == 0.0) {
                    m_comm.CommRHStoBlocks(B, RHSblocks);
                } else {
                    double[] RES = B.ToArray();
                    m_op.OperatorMatrix.SpMV(-1.0, X, 1.0, RES);
                    m_comm.CommRHStoBlocks(RES, RHSblocks);
                }

                for(int iBlock = 0; iBlock < m_BlockSolvers.Length; iBlock++) {
                    if (m_BlockSolvers[iBlock] != null) {
                        m_BlockSolvers[iBlock].Solve(Xblocks[iBlock], RHSblocks[iBlock]);
                        if(m_OverlapScaling != null) {
                            var scale = m_OverlapScaling[iBlock];
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
            if(m_BlockSolvers == null) {
                beits += m_BlockSolvers.Sum(pardiso => pardiso?.UsedMemory() ?? 0);
            }

            return beits;
        }
    }
}

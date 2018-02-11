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

// extensive testing
#define MATLAB_CHECK

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using ilPSP.Kraypis;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP.LinSolvers.PARDISO;
using System.Collections;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid;
using MPI.Wrappers;
using System.Runtime.InteropServices;
using ilPSP.Tracing;

namespace BoSSS.Solution.Multigrid {


    public class Schwarz : ISolverSmootherTemplate, ISolverWithCallback {

        /// <summary>
        /// Abstract base class, template for different strategies 
        /// </summary>
        public abstract class BlockingStrategy {

            /// <summary>
            /// Returns lists which form the blocks of the additive-Schwarz domain decomposition.
            /// </summary>
            /// <param name="op"></param>
            /// <returns>
            /// - outer enumeration: correcponts to domain-decomposition blocks
            /// - inner index: indices within the sub-blocks
            /// - content: local cell indices which form the respective additive-Schwarz block (<see cref="MultigridOperator."/>
            /// </returns>
            abstract internal IEnumerable<List<int>> GetBlocking(MultigridOperator op);

            /// <summary>
            /// Number of blocs returned by <see cref="GetBlocking(MultigridOperator)"/>
            /// </summary>
            internal abstract int GetNoOfBlocks(MultigridOperator op);
        }

        /// <summary>
        /// Additive-Schwarz blocks which are formed from coarser multigrid-levels.
        /// </summary>
        public class MultigridBlocks : BlockingStrategy {

            /// <summary>
            /// Number of multigrid levels to step down.
            /// </summary>
            public int Depth {
                get {
                    return m_Depht;
                }
                set {
                    if (value < 0) {
                        throw new ArgumentException();
                    }
                    m_Depht = value;
                }
            }

            int m_Depht = 1;

            /// <summary>
            /// Returns the multigrid blocking.
            /// </summary>
            internal override IEnumerable<List<int>> GetBlocking(MultigridOperator op) {

                AggregationGrid thisLevel = op.Mapping.AggGrid;

                List<AggregationGrid> blockLevelS = new List<AggregationGrid>();
                blockLevelS.Add(thisLevel);
                MultigridOperator blokOp = op;
                for (int i = 0; i < this.Depth; i++) {
                    if (blokOp.CoarserLevel == null)
                        throw new NotSupportedException("Not enough multigrid levels set to support a depth of " + m_Depht + ".");
                    blokOp = blokOp.CoarserLevel;
                    blockLevelS.Add(blokOp.Mapping.AggGrid);
                }
                AggregationGrid blckLevel = blockLevelS.Last(); // the cells of this level form the additive-Schwarz blocks
                int NoBlocks = blckLevel.iLogicalCells.NoOfLocalUpdatedCells; // each cell of 'blckLevel' forms a block
                List<int>[] Blocks = NoBlocks.ForLoop(l => new List<int>());

#if DEBUG
                bool[] checkOnce = new bool[thisLevel.iLogicalCells.NoOfLocalUpdatedCells];
#endif
                for (int iBlk = 0; iBlk < NoBlocks; iBlk++) {
                    if (blockLevelS.Count == 0) {
                        Blocks[iBlk].Add(iBlk); // the cell itself is the multigrid block (either Depth is 0, or no more MG level available).
                    } else {
                        int[] CoarseCell = blckLevel.jCellCoarse2jCellFine[iBlk];
                        CollectBlock(Blocks[iBlk], blockLevelS, 0, CoarseCell);
                    }
#if DEBUG
                    foreach (int j in Blocks[iBlk]) {
                        Debug.Assert(j >= 0);
                        Debug.Assert(j < checkOnce.Length);
                        Debug.Assert(checkOnce[j] == false);
                        checkOnce[j] = true;
                    }

#endif
                }

#if DEBUG
                for (int j = 0; j < checkOnce.Length; j++) {
                    Debug.Assert(checkOnce[j] == true);
                }
#endif

                return Blocks;
            }

            /// <summary>
            /// %
            /// </summary>
            internal override int GetNoOfBlocks(MultigridOperator op) {
                return GetBlocking(op).Count();
            }



            void CollectBlock(List<int> output, List<AggregationGrid> blockLevelS, int RecDepth, int[] CoarseCell) {

                if (RecDepth == blockLevelS.Count - 2) {
#if DEBUG
                    foreach (int jFine in CoarseCell)
                        Debug.Assert(output.Contains(jFine) == false);
#endif
                    output.AddRange(CoarseCell);
                } else {
                    AggregationGrid blockLevel = blockLevelS[blockLevelS.Count - 2 - RecDepth];
                    int[][] C2F = blockLevel.jCellCoarse2jCellFine;
                    foreach (int jFine in CoarseCell) {
                        CollectBlock(output, blockLevelS, RecDepth + 1, C2F[jFine]);
                    }
                }
            }


        }


        /// <summary>
        /// creates a fixed number of blocks by using METIS
        /// </summary>
        public class METISBlockingStrategy : BlockingStrategy {

            /// <summary>
            /// Number of parts/additive Schwarz blocks on current MPI process.
            /// </summary>
            public int NoOfPartsPerProcess = 4;

            internal override IEnumerable<List<int>> GetBlocking(MultigridOperator op) {
                var MgMap = op.Mapping;

                //if(!M.RowPartitioning.Equals(MgMap.Partitioning))
                //    throw new ArgumentException("Row partitioning mismatch.");
                //if(!M.ColPartition.Equals(MgMap.Partitioning))
                //    throw new ArgumentException("Column partitioning mismatch.");

                var ag = MgMap.AggGrid;

                int JComp = ag.iLogicalCells.NoOfLocalUpdatedCells; // number of local cells
                int[] xadj = new int[JComp + 1];
                List<int> adjncy = new List<int>();
                for (int j = 0; j < JComp; j++) {
                    int[] neigh_j = ag.iLogicalCells.CellNeighbours[j];
                    foreach (int jNeigh in neigh_j) {
                        //adjncy.AddRange(neigh_j);
                        if (jNeigh < JComp)
                            adjncy.Add(jNeigh);
                        else
                            Console.WriteLine("Skipping external cell");
                    }
                    xadj[j + 1] = xadj[j] + neigh_j.Length;
                }

                int MPIrank, MPIsize;
                MPI.Wrappers.csMPI.Raw.Comm_Rank(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out MPIrank);
                MPI.Wrappers.csMPI.Raw.Comm_Size(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out MPIsize);
                //if (MPIrank == 1)
                //    NoOfParts = 1;
                //Debugger.Launch();

                int[] part = new int[JComp];
                {
                    if (NoOfPartsPerProcess > 1) {
                        int ncon = 1;
                        int edgecut = 0;
                        int[] options = null; //new int[] { 0, 0, 0, 0, 0 };
                        METIS.PartGraphKway(
                            ref JComp, ref ncon,
                            xadj,
                            adjncy.ToArray(),
                            null,
                            null,
                            null,
                            ref NoOfPartsPerProcess,
                            null,
                            null,
                            options,
                            ref edgecut,
                            part);
                    } else {
                        part.SetAll(0);
                    }
                }

                {
                    var _Blocks = NoOfPartsPerProcess.ForLoop(i => new List<int>((int)Math.Ceiling(1.1 * JComp / NoOfPartsPerProcess)));
                    for (int j = 0; j < JComp; j++) {
                        _Blocks[part[j]].Add(j);

                    }

                    return _Blocks;
                }
            }


            /// <summary>
            /// %
            /// </summary>
            internal override int GetNoOfBlocks(MultigridOperator op) {
                return NoOfPartsPerProcess;
            }
        }


        public BlockingStrategy m_BlockingStrategy;
        MultigridOperator m_MgOp;

#if DEBUG
        bool m_MatlabParalellizationCheck = false;
#endif

        /// <summary>
        /// Debugging and checking of algorithm parallelization using the <see cref="ilPSP.Connectors.Matlab.BatchmodeConnector"/>.
        /// - only supported in DEBUG configuration
        /// - the checks are serial, no scaling can be expected
        /// - very expensive, only for debugging of small systems
        /// </summary>
        public bool MatlabParalellizationCheck {
            get {
#if DEBUG
                return m_MatlabParalellizationCheck;
#else
                return false;
#endif
            }
            set {
#if DEBUG
                m_MatlabParalellizationCheck = value;
#else
                if(value == true)
                    throw new NotSupportedException("Only supported in DEBUG mode.");
#endif
            }
        }


        public void Init(MultigridOperator op) {
            using (new FuncTrace()) {
                if (m_MgOp != null) {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // someone is trying to re-use this solver: see if the settings permit that
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    if (op.LevelIndex != m_MgOp.LevelIndex)
                        throw new ArgumentException("Re-use on different level not possible.");
                    if (!this.MtxFull._RowPartitioning.EqualsPartition(op.OperatorMatrix._RowPartitioning))
                        throw new ArgumentException("Matrix has changed, unable to re-use");
                    if (!this.MtxFull._ColPartitioning.EqualsPartition(op.OperatorMatrix._ColPartitioning))
                        throw new ArgumentException("Matrix has changed, unable to re-use");
#if DEBUG
                    if (!object.ReferenceEquals(this.MtxFull, op.OperatorMatrix)) {
                        BlockMsrMatrix Check = this.MtxFull.CloneAs();
                        Check.Acc(-1.0, op.OperatorMatrix);
                        if (Check.InfNorm() != 0.0) {
                            throw new ArgumentException("Matrix has changed, unable to re-use");
                        }
                    }
#endif
                    if(this.m_BlockingStrategy.GetNoOfBlocks(op) != this.blockSolvers.Count()) {
                        throw new ArgumentException("Blocking, unable to re-use");
                    }


                    return;
                }

                var Mop = op.OperatorMatrix;
                var MgMap = op.Mapping;
                this.m_MgOp = op;
                int myMpiRank = MgMap.MpiRank;
                int myMpisize = MgMap.MpiSize;

                if (!Mop.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!Mop.ColPartition.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Column partitioning mismatch.");

                var ag = MgMap.AggGrid;

                int JComp = ag.iLogicalCells.NoOfLocalUpdatedCells;
                int JGhost = ag.iLogicalCells.NoOfExternalCells;

#if DEBUG
                ilPSP.Connectors.Matlab.BatchmodeConnector matlab;
                if (m_MatlabParalellizationCheck)
                    matlab = new ilPSP.Connectors.Matlab.BatchmodeConnector();
                else
                    matlab = null;
#endif

                //Mop.Clear();
                //for(int i = Mop.RowPartitioning.i0; i < Mop.RowPartitioning.iE; i++) {
                //    Mop[i, i] = i + 1;
                //}


                // get cell blocks
                // ===============

                var _Blocks = this.m_BlockingStrategy.GetBlocking(op);
                int NoOfSchwzBlocks = _Blocks.Count();

                // test cell blocks
                // ================
#if DEBUG
                {
                    // ensure that each cell is used exactly once, among all blocks
                    bool[] test = new bool[ag.iLogicalCells.NoOfLocalUpdatedCells];
                    foreach (var bi in _Blocks) {
                        foreach (int j in bi) {
                            Debug.Assert(test[j] == false);
                            test[j] = true;
                        };
                    }
                    for (int i = 0; i < test.Length; i++)
                        Debug.Assert(test[i] == true);
                }
#endif

                // extend blocks according to desired overlap
                // ==========================================
                {
                    BitArray marker = new BitArray(JComp + JGhost);

                    if (Overlap < 0)
                        throw new ArgumentException();
                    if (Overlap > 0) {
                        if (Overlap > 1 && Mop.RowPartitioning.MpiSize > 1) {
                            throw new NotSupportedException("In MPI parallel runs, the maximum supported overlap for the Schwarz preconditioner is 1.");
                        }

                        foreach (List<int> bi in _Blocks) { // loop over blocks...
                            marker.SetAll(false); // marks all cells which are members of the block
                            foreach (int jcomp in bi)
                                marker[jcomp] = true;

                            // determine overlap regions
                            for (int k = 0; k < Overlap; k++) {
                                int Jblock = bi.Count;
                                for (int j = 0; j < Jblock; j++) {
                                    int jCell = bi[j];
                                    int[] Neighs = ag.iLogicalCells.CellNeighbours[jCell];
                                    foreach (int jNeigh in Neighs) {
                                        if (marker[jNeigh] == false) {
                                            // neighbor cell is not already a member of the block
                                            // => add it.
                                            bi.Add(jNeigh);
                                            marker[jNeigh] = true;
                                        }
                                    }

                                }
                            }

                            bi.Sort();
                        }
                    }

                    BlockCells = _Blocks.Select(list => list.ToArray()).ToArray();
                }


                // convert cell blocks to DOF blocks
                // =================================

                List<int>[] BlkIdx_gI_lR; //  for each Schwarz block, (global) indices in the local range 
                List<int>[] BlkIdx_gI_eR; //  for each Schwarz block, (global) indices of external rows and columns
                List<int>[] TempRowIdx_gI; // for each Schwarz block, (global) indices into the temporary matrix
                List<int>[] BlkIdx_lI_eR; //  for each Schwarz block, (local)  indices of external rows and columns
                List<int>[] LocalBlocks_i0, LocalBlocks_N; // blocking of the Schwarz-Blocks.

                // for matrix 'ExternalRowsTemp': which rows of 'Mop' are required locally
                List<int> ExternalRowsIndices, ExternalRows_BlockI0, ExternalRows_BlockN;
                {
                    int Jup = MgMap.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                    int Jgh = MgMap.AggGrid.iLogicalCells.NoOfExternalCells;

                    int LocalizedBlockCounter = 0;

                    BlkIdx_gI_lR = NoOfSchwzBlocks.ForLoop(iPart => new List<int>(BlockCells[iPart].Length * MgMap.MaximalLength));
                    BlkIdx_gI_eR = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                    LocalBlocks_i0 = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                    LocalBlocks_N = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());

                    TempRowIdx_gI = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                    BlkIdx_lI_eR = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());


                    ExternalRowsIndices = new List<int>();
                    ExternalRows_BlockI0 = new List<int>();
                    ExternalRows_BlockN = new List<int>();

                    for (int iPart = 0; iPart < NoOfSchwzBlocks; iPart++) { // loop over parts...
                        int[] bc = BlockCells[iPart];
                        var biI = BlkIdx_gI_lR[iPart];
                        var biE = BlkIdx_gI_eR[iPart];
                        var l1 = TempRowIdx_gI[iPart];
                        var l2 = BlkIdx_lI_eR[iPart];
                        var LBBi0 = LocalBlocks_i0[iPart];
                        var LBBN = LocalBlocks_N[iPart];


                        int Jblock = bc.Length;
                        int anotherCounter = 0;

                        for (int jblk = 0; jblk < Jblock; jblk++) { // loop over cells in blocks...
                            int j = bc[jblk];
                            int N = MgMap.GetLength(j);

                            if (j < Jup) {
                                // locally updated cell
                                int i0 = MgMap.GlobalUniqueIndex(0, j, 0);

                                for (int n = 0; n < N; n++) {
                                    biI.Add(i0 + n);
                                }
                            } else {
                                // external cell
                                int i0E = MgMap.GlobalUniqueIndex(0, j, 0); // 
                                int i0L = MgMap.LocalUniqueIndex(0, j, 0); // 
                                ExternalRows_BlockI0.Add(LocalizedBlockCounter);
                                ExternalRows_BlockN.Add(N);
                                //LEBi0.Add(LocalizedBlockCounter);
                                //LEBn.Add(N);
                                for (int n = 0; n < N; n++) {
                                    biE.Add(i0E + n);
                                    ExternalRowsIndices.Add(i0E + n);
                                    l1.Add(LocalizedBlockCounter + n);
                                    l2.Add(i0L + n);
                                    Debug.Assert(Mop._RowPartitioning.FindProcess(i0E + n) != myMpiRank);
                                }

                                LocalizedBlockCounter += N;
                            }

                            LBBi0.Add(anotherCounter);
                            LBBN.Add(N);

                            anotherCounter += N;
                        }
                    }

                    //this.BlockIndices = _LocallyStoredBlockIndices.Select(bi => bi.ToArray()).ToArray();
                }


                // get rows for blocks that use external cells
                // ===========================================

#if DEBUG
                {
                    if (Overlap == 0) {
                        Debug.Assert(ExternalRowsIndices.Count == 0);
                        Debug.Assert(ExternalRows_BlockI0.Count == 0);
                        Debug.Assert(ExternalRows_BlockN.Count == 0);
                    }

                    foreach (var bi in BlkIdx_gI_lR) {
                        foreach (int idx in bi) {
                            Debug.Assert(idx >= m_MgOp.Mapping.i0);
                            Debug.Assert(idx < m_MgOp.Mapping.iE);
                        }
                    }

                    foreach (var ei in BlkIdx_gI_eR) {
                        foreach (int idx in ei) {
                            Debug.Assert(idx < m_MgOp.Mapping.i0 || idx >= m_MgOp.Mapping.iE);
                        }
                    }


                    int LL = m_MgOp.Mapping.LocalLength;
                    int jMax = m_MgOp.Mapping.AggGrid.iLogicalCells.NoOfCells - 1;
                    int LE = m_MgOp.Mapping.LocalUniqueIndex(0, jMax, 0) + m_MgOp.Mapping.GetLength(jMax);


                    foreach (var ci in BlkIdx_lI_eR) {
                        foreach (int idx in ci) {
                            Debug.Assert(idx >= LL);
                            Debug.Assert(idx < LE);
                        }
                    }

                    if (m_MatlabParalellizationCheck) {
                        int globalBlockCounter = 0;
                        for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                            int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                            if (rankCounter == myMpiRank)
                                Debug.Assert(rank_NoBlks == NoOfSchwzBlocks);

                            for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {
                                double[] vec;
                                if (rankCounter == myMpiRank) {
                                    vec = ArrayTools.Cat(BlkIdx_gI_lR[iBlock], BlkIdx_gI_eR[iBlock]).Select(ii => ((double)(ii + 1))).ToArray();
                                } else {
                                    vec = new double[0];
                                }

                                matlab.PutVector(vec, string.Format("BlockIdx{0}", globalBlockCounter));

                                globalBlockCounter++;
                                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                            }

                        }
                    }
                }
#endif


                BlockMsrMatrix ExternalRowsTemp;
                if (myMpisize > 1 && Overlap > 0) {
                    //int NoOfLocalRows = _ExternalBlockIndices.Sum(L => L.Count);

                    BlockPartitioning PermRow = new BlockPartitioning(ExternalRowsIndices.Count, ExternalRows_BlockI0, ExternalRows_BlockN, Mop.MPI_Comm, i0isLocal: true);

                    // Remark: we use a permutation matrix for MPI-exchange of rows

                    BlockMsrMatrix Perm = new BlockMsrMatrix(PermRow, Mop._RowPartitioning);
                    for (int iRow = 0; iRow < ExternalRowsIndices.Count; iRow++) {
                        Debug.Assert(Mop._RowPartitioning.IsInLocalRange(ExternalRowsIndices[iRow]) == false);
                        Perm[iRow + PermRow.i0, ExternalRowsIndices[iRow]] = 1;
                    }

                    ExternalRowsTemp = BlockMsrMatrix.Multiply(Perm, Mop);

#if DEBUG
                    if (m_MatlabParalellizationCheck) {
                        matlab.PutSparseMatrix(Perm, "Perm");
                        matlab.PutSparseMatrix(ExternalRowsTemp, "ExternalRowsTemp");
                    }
#endif
                } else {
                    ExternalRowsTemp = null;
                }

                ExternalRowsIndices = null;
                ExternalRows_BlockI0 = null;
                ExternalRows_BlockN = null;




                // create solvers
                // ==============


                {
                    blockSolvers = new ISparseSolver[NoOfSchwzBlocks];

#if DEBUG
                    List<BlockMsrMatrix> Blocks = new List<BlockMsrMatrix>();
#endif
                    for (int iPart = 0; iPart < NoOfSchwzBlocks; iPart++) {
                        var bi = BlkIdx_gI_lR[iPart];

                        int Bsz;
                        if (MgMap.MinimalLength == MgMap.MaximalLength)
                            Bsz = MgMap.MaximalLength;
                        else
                            Bsz = 1;

                        var l1 = TempRowIdx_gI[iPart];

                        //if (M.RowPartitioning.MpiSize > 1) {
                        //    int i0Proc = M.RowPartitioning.i0;
                        //    bi = bi.CloneAs();
                        //    for (int i = 0; i < bi.Length; i++) {
                        //        bi[i] += i0Proc;
                        //    }
                        //}

                        BlockPartitioning localBlocking = new BlockPartitioning(bi.Count + l1.Count, LocalBlocks_i0[iPart], LocalBlocks_N[iPart], csMPI.Raw._COMM.SELF);

                        if (l1.Count > 0) {
                            // convert the indices into 'ExternalRowsTemp' to global indices
                            int l1L = l1.Count;
                            int offset = ExternalRowsTemp._RowPartitioning.i0;
                            for (int i = 0; i < l1L; i++)
                                l1[i] += offset;
                        }

                        BlockMsrMatrix Block = new BlockMsrMatrix(localBlocking, localBlocking);// bi.Length, bi.Length, Bsz, Bsz);
                        Mop.WriteSubMatrixTo(Block, bi, default(int[]), bi, default(int[]));
                        if (l1.Count > 0) {
                            int offset = bi.Count;
                            int[] targRows = l1.Count.ForLoop(i => i + offset);

                            var biE = BlkIdx_gI_eR[iPart];
                            int[] extTargCols = biE.Count.ForLoop(i => i + offset);

                            Mop.AccSubMatrixTo(1.0, Block, bi, default(int[]), new int[0], default(int[]), biE, extTargCols);
                            ExternalRowsTemp.AccSubMatrixTo(1.0, Block, l1, targRows, bi, default(int[]), biE, extTargCols);
                        }
#if DEBUG
                        if (m_MatlabParalellizationCheck != null) {
                            Blocks.Add(Block);
                        }
#endif
                        blockSolvers[iPart] = new PARDISOSolver() {
                            CacheFactorization = true
                        };
                        //blockSolvers[iPart] = new FullDirectSolver();
                        //blockSolvers[iPart] = new ilPSP.LinSolvers.MUMPS.MUMPSSolver(MPI: false);
                        blockSolvers[iPart].DefineMatrix(Block);
                    }

#if DEBUG
                    if (m_MatlabParalellizationCheck) {
                        int globalBlockCounter = 0;
                        for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                            int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                            for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {
                                BlockMsrMatrix Block;
                                if (rankCounter == myMpiRank) {
                                    Block = Blocks[iBlock];
                                } else {
                                    Block = null;
                                }

                                matlab.PutSparseMatrix(Block, string.Format("Block{0}", globalBlockCounter));

                                globalBlockCounter++;
                                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

                            }

                        }
                    }
#endif
                }

                // Record required indices
                // =======================
                {
                    this.BlockIndices_Local = new int[NoOfSchwzBlocks][];
                    this.BlockIndices_External = new int[NoOfSchwzBlocks][];
                    int LocalI0 = MgMap.i0;
                    int LocalLength = MgMap.LocalLength;

                    for (int iBlock = 0; iBlock < NoOfSchwzBlocks; iBlock++) {
                        var _bi = BlkIdx_gI_lR[iBlock];
                        int L = _bi.Count;
                        int[] bil = new int[L];
                        this.BlockIndices_Local[iBlock] = bil;

                        for (int l = 0; l < L; l++) {
                            bil[l] = _bi[l] - LocalI0;
                            Debug.Assert(bil[l] >= 0);
                            Debug.Assert(bil[l] < MgMap.LocalLength);
                        }

                        var _biE = BlkIdx_lI_eR[iBlock];
                        if (_biE.Count > 0) {
                            this.BlockIndices_External[iBlock] = _biE.ToArray();
                        }
                    }
                }


                this.MtxFull = Mop;

                if (CoarseSolver != null) {
                    CoarseSolver.Init(op.CoarserLevel);
                }

                // Debug & Test-Code 
                // =================
#if DEBUG
                if (m_MatlabParalellizationCheck) {
                    Console.WriteLine("Matlab dir: " + matlab.WorkingDirectory);

                    matlab.PutSparseMatrix(Mop, "Full");
                    int GlobalNoOfBlocks = NoOfSchwzBlocks.MPISum();



                    for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                        matlab.Cmd("BlockErr({0} + 1, 1) = norm( Block{0} - Full( BlockIdx{0}, BlockIdx{0} ), inf );", iGlbBlock);
                    }

                    Random rnd = new Random(myMpiRank);
                    double[] testRHS = new double[MgMap.LocalLength];
                    for (int i = 0; i < testRHS.Length; i++) {
                        testRHS[i] = rnd.NextDouble();
                    }
                    matlab.PutVector(testRHS, "testRHS");

                    MPIexchange<double[]> ResExchange = new MPIexchange<double[]>(MgMap, testRHS);
                    ResExchange.TransceiveStartImReturn();
                    ResExchange.TransceiveFinish(0.0);

                    int offset = MgMap.LocalLength;

                    int g = 0;
                    for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                        int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                        for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {
                            double[] SubVec;
                            if (rankCounter == myMpiRank) {
                                int LL = this.BlockIndices_Local[iBlock].Length;
                                int LE;
                                if (this.BlockIndices_External[iBlock] != null) {
                                    LE = this.BlockIndices_External[iBlock].Length;
                                } else {
                                    LE = 0;
                                }
                                int L = LL + LE;

                                SubVec = new double[L];
                                for (int i = 0; i < LL; i++) {
                                    SubVec[i] = testRHS[this.BlockIndices_Local[iBlock][i]];
                                }
                                if (LE > 0) {
                                    for (int i = 0; i < LE; i++) {
                                        SubVec[i + LL] = ResExchange.Vector_Ext[this.BlockIndices_External[iBlock][i] - offset];
                                    }
                                }
                            } else {
                                SubVec = new double[0];
                            }

                            matlab.PutVector(SubVec, "SubVec" + g);

                            g++;
                        }
                    }

                    for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                        matlab.Cmd("RhsErr({0} + 1, 1) = norm( SubVec{0} - testRHS( BlockIdx{0} ), inf );", iGlbBlock);
                    }

                    double[] testX = new double[testRHS.Length];
                    MPIexchangeInverse<double[]> XExchange = new MPIexchangeInverse<double[]>(MgMap, testX);

                    g = 0;
                    for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                        int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                        for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {

                            if (rankCounter == myMpiRank) {
                                int LL = this.BlockIndices_Local[iBlock].Length;
                                int LE;
                                if (this.BlockIndices_External[iBlock] != null) {
                                    LE = this.BlockIndices_External[iBlock].Length;
                                } else {
                                    LE = 0;
                                }
                                int L = LL + LE;


                                for (int i = 0; i < LL; i++) {
                                    testX[this.BlockIndices_Local[iBlock][i]] += (g + 1);
                                }
                                if (LE > 0) {
                                    for (int i = 0; i < LE; i++) {
                                        XExchange.Vector_Ext[this.BlockIndices_External[iBlock][i] - offset] += (g + 1);
                                    }
                                }
                            } else {
                                //nop
                            }

                            g++;
                        }
                    }
                    XExchange.TransceiveStartImReturn();
                    XExchange.TransceiveFinish(1.0);

                    matlab.Cmd("testXref = zeros({0},1);", MgMap.TotalLength);
                    for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                        matlab.Cmd("testXref(BlockIdx{0},1) = testXref(BlockIdx{0},1) + ({0} + 1);", iGlbBlock);
                    }

                    matlab.PutVector(testX, "testX");
                    matlab.Cmd("testXErr = norm(testX - testXref, inf);");

                    MultidimensionalArray BlockErr = MultidimensionalArray.Create(GlobalNoOfBlocks, 1);
                    MultidimensionalArray RhsErr = MultidimensionalArray.Create(GlobalNoOfBlocks, 1);
                    MultidimensionalArray testXErr = MultidimensionalArray.Create(1, 1);

                    matlab.GetMatrix(BlockErr, "BlockErr");
                    matlab.GetMatrix(RhsErr, "RhsErr");
                    matlab.GetMatrix(testXErr, "testXErr");

                    matlab.Execute();

                    for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                        Console.WriteLine("Block #{0} Error (external? ) " + BlockErr[iGlbBlock, 0], iGlbBlock);
                        Console.WriteLine("RHS #{0} Error " + RhsErr[iGlbBlock, 0], iGlbBlock);
                        Debug.Assert(BlockErr[iGlbBlock, 0] == 0);
                        Debug.Assert(RhsErr[iGlbBlock, 0] == 0);
                    }

                    Console.WriteLine("X Error " + testXErr[0, 0]);
                    Debug.Assert(testXErr[0, 0] == 0.0);

                    matlab.Dispose();
                }
#endif
            }
        }

        /*
        class FullDirectSolver : ISparseSolver {
            public void DefineMatrix(IMutableMatrixEx M) {
                fullMtx = M.ToFullMatrixOnProc0();
            }

            MultidimensionalArray fullMtx;

            public void Dispose() {
            }

            public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
                where Tunknowns : IList<double>
                where Trhs : IList<double> {
                double[] _x = x.ToArray();
                double[] _rhs = rhs.ToArray();
                fullMtx.Solve(_x, _rhs);
                x.SetV(_x);
                return new SolverResult() { Converged = true, NoOfIterations = 1, };
            }
        }
        */

        BlockMsrMatrix MtxFull;


        ISparseSolver[] blockSolvers;

        /// <summary>
        /// List of local cell indices for each Schwarz block.
        /// - 1st index: block
        /// - 2nd index: enumeration; list of composite cells which make up the respective block
        /// - content: local cell index
        /// </summary>
        int[][] BlockCells;


        /// <summary>
        /// List of local row and column indices for each Schwarz block.
        /// - 1st index: local block index
        /// - 2nd index: enumeration
        /// - content: local row and column indices (i.e. in between 0 and <see cref="MultigridMapping.LocalLength"/>).
        /// </summary>
        int[][] BlockIndices_Local;

        /// <summary>
        /// List of external row and column indices for each Schwarz block.
        /// - 1st index: local block index
        /// - 2nd index: enumeration
        /// - content: external index (can be used e.g. as index into <see cref="MPIexchange{T}.Vector_Ext"/>, <see cref="MPIexchangeInverse{T}{T}.Vector_Ext"/>, but taking into account a local offset).
        /// </summary>
        int[][] BlockIndices_External;




        private int m_Overlap = 1;


        /// <summary>
        /// Overlap of the Schwarz blocks, in number-of-cells.
        /// </summary>
        public int Overlap {
            get {
                return m_Overlap;
            }
            set {
                if (value < 0) {
                    throw new ArgumentException();
                }
                if (value > 2) {
                    throw new ArgumentException();
                }
                m_Overlap = value;
            }
        }

        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }

        public int m_MaxIterations = 1;


        public ISolverSmootherTemplate CoarseSolver;

        public bool CoarseSolverIsMultiplicative = true;

        int NoIter = 0;


        //static void SingleFilter(double[] V) {
        //    for(int i = 0; i < V.Length; i++) {
        //        float vi = (float)(V[i]);
        //        V[i] = vi;
        //    }
        //}


        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (var tr = new FuncTrace()) {
                int NoParts = this.BlockIndices_Local.Length;

                // --------
                // Reminder: we solve for a correction, the basic idea is:
                //
                // X0 is current solution, state of 'X' at input time
                //    Residual = B - M*X0;
                // now solve (approximately)
                //    M*Xc = Residual,
                // if we would solve exactly, i.e.
                //    Xc = M^{-1}*Residual
                // then X=(X0+Xc) is an exact solution of 
                //    M*X = B
                // --------

                double[] Res = new double[B.Count];
                MPIexchange<double[]> ResExchange;
                MPIexchangeInverse<U> XExchange;
                if (Overlap > 0) {
                    ResExchange = new MPIexchange<double[]>(this.m_MgOp.Mapping, Res);
                    XExchange = new MPIexchangeInverse<U>(this.m_MgOp.Mapping, X);
                } else {
                    ResExchange = null;
                    XExchange = null;
#if DEBUG
                    foreach (var ciE in BlockIndices_External) {
                        Debug.Assert(ciE == null || ciE.Length <= 0);
                    }
#endif
                }


                int LocLength = m_MgOp.Mapping.LocalLength;

                for (int iIter = 0; iIter < m_MaxIterations; iIter++) {
                    this.NoIter++;

                    Res.SetV(B);
                    this.MtxFull.SpMV(-1.0, X, 1.0, Res);

                    if (IterationCallback != null)
                        IterationCallback(iIter, X.ToArray(), Res.CloneAs(), this.m_MgOp);

                    using (new BlockTrace("coarse_solve_level" + this.m_MgOp.LevelIndex, tr)) {
                        if (CoarseSolver != null) {
                            var XC = X.ToArray().CloneAs();
                            double[] bc = new double[m_MgOp.CoarserLevel.Mapping.LocalLength];// = Res.CloneAs();
                            m_MgOp.CoarserLevel.Restrict(Res.CloneAs(), bc);
                            double[] xc = new double[bc.Length];
                            CoarseSolver.Solve(xc, bc);
                            //SingleFilter(xc);
                            m_MgOp.CoarserLevel.Prolongate(1, XC, 1, xc);
                            X.AccV(1.0, XC);

                            if (CoarseSolverIsMultiplicative) {
                                Res.SetV(B);
                                this.MtxFull.SpMV(-1.0, X, 1.0, Res);
                            }
                        }
                    }

                    if (Overlap > 0) {
                        ResExchange.TransceiveStartImReturn();
                        ResExchange.TransceiveFinish(0.0);
                    }

                    using (new BlockTrace("block_solve_level" + this.m_MgOp.LevelIndex, tr)) {
                        for (int iPart = 0; iPart < NoParts; iPart++) {
                            int[] ci = BlockIndices_Local[iPart];
                            int[] ciE = BlockIndices_External[iPart];
                            int L = ci.Length;
                            if (ciE != null)
                                L += ciE.Length;

                            double[] bi = new double[L];
                            double[] xi = new double[L];

                            // extract block part of residual
                            bi.AccV(1.0, Res, default(int[]), ci);
                            if (ciE != null && ciE.Length > 0)
                                bi.AccV(1.0, ResExchange.Vector_Ext, default(int[]), ciE, acc_index_shift: ci.Length, b_index_shift: (-LocLength));

                            blockSolvers[iPart].Solve(xi, bi);
                            //SingleFilter(xi);

                            // accumulate block solution 'xi' to global solution 'X'
                            X.AccV(1.0, xi, ci, default(int[]));
                            if (ciE != null && ciE.Length > 0)
                                XExchange.Vector_Ext.AccV(1.0, xi, ciE, default(int[]), acc_index_shift: (-LocLength), b_index_shift: ci.Length);
                        }
                    }

                    if (Overlap > 0) {
                        // block solutions stored on *external* indices will be accumulated on other processors.
                        XExchange.TransceiveStartImReturn();
                        XExchange.TransceiveFinish(1.0);

                        if (iIter < m_MaxIterations - 1)
                            XExchange.Vector_Ext.ClearEntries();
                    }
                }
            }
        }


        public int IterationsInNested {
            get {
                if (this.CoarseSolver != null)
                    return this.CoarseSolver.IterationsInNested + this.CoarseSolver.ThisLevelIterations;
                else
                    return 0;
            }
        }

        public int ThisLevelIterations {
            get {
                return this.m_MaxIterations;
            }
        }

        public bool Converged {
            get {
                return m_Converged;
            }
        }

        bool m_Converged = false;

        public void ResetStat() {
            this.m_Converged = false;
            this.NoIter = 0;
            if (this.CoarseSolver != null)
                this.CoarseSolver.ResetStat();
        }
    }


}

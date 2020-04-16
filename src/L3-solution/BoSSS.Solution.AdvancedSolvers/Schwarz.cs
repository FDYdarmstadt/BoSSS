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

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// Additive Schwarz method with optional, multiplicative coarse-grid correction.
    /// </summary>
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
            /// - outer enumeration: corresponds to domain-decomposition blocks
            /// - inner index: indices within the sub-blocks
            /// - content: local cell indices which form the respective additive-Schwarz block (<see cref="MultigridOperator"/>
            /// </returns>
            abstract internal IEnumerable<List<int>> GetBlocking(MultigridOperator op);

            /// <summary>
            /// Number of blocs returned by <see cref="GetBlocking(MultigridOperator)"/>
            /// </summary>
            internal abstract int GetNoOfBlocks(MultigridOperator op);
        }

        /// <summary>
        /// Additive-Schwarz blocks which are formed from cells of coarser multigrid-levels.
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
                AggregationGridData thisLevel = op.Mapping.AggGrid;

                List<AggregationGridData> blockLevelS = new List<AggregationGridData>();
                blockLevelS.Add(thisLevel);
                MultigridOperator blokOp = op;
                for (int i = 0; i < this.Depth; i++) {
                    if (blokOp.CoarserLevel == null)
                        throw new NotSupportedException("Not enough multigrid levels set to support a depth of " + m_Depht + ".");
                    blokOp = blokOp.CoarserLevel;
                    blockLevelS.Add(blokOp.Mapping.AggGrid);
                }
                AggregationGridData blckLevel = blockLevelS.Last(); // the cells of this level form the additive-Schwarz blocks
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



            void CollectBlock(List<int> output, List<AggregationGridData> blockLevelS, int RecDepth, int[] CoarseCell) {

                if (RecDepth == blockLevelS.Count - 2) {
#if DEBUG
                    foreach (int jFine in CoarseCell)
                        Debug.Assert(output.Contains(jFine) == false);
#endif
                    output.AddRange(CoarseCell);
                } else {
                    AggregationGridData blockLevel = blockLevelS[blockLevelS.Count - 2 - RecDepth];
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
                    Debug.Assert(xadj[j] == adjncy.Count);

                    int[] neigh_j = ag.iLogicalCells.CellNeighbours[j];
                    int nCnt = 0;
                    foreach (int jNeigh in neigh_j) {
                        //adjncy.AddRange(neigh_j);
                        if (jNeigh < JComp) {
                            adjncy.Add(jNeigh);
                            nCnt++;
                        } else {
                            //Console.WriteLine("Skipping external cell");
                        }
                    }
                    xadj[j + 1] = xadj[j] + nCnt;
                }
                Debug.Assert(xadj[JComp] == adjncy.Count);

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
                        int[] options = new int[METIS.METIS_NOPTIONS]; 
                        METIS.SETDEFAULTOPTIONS(options);

                        options[(int)METIS.OptionCodes.METIS_OPTION_NCUTS] = 1; // 
                        options[(int)METIS.OptionCodes.METIS_OPTION_NITER] = 10; // This is the default refinement iterations
                        options[(int)METIS.OptionCodes.METIS_OPTION_UFACTOR] = 30; // Maximum imbalance of 3 percent (this is the default kway clustering)
                        options[(int)METIS.OptionCodes.METIS_OPTION_NUMBERING] = 0;

                        Debug.Assert(xadj.Where(idx => idx > adjncy.Count).Count() == 0);
                        Debug.Assert(adjncy.Where(j => j >= JComp).Count() == 0);

                        METIS.PARTGRAPHKWAY(
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
                    List<List<int>> _Blocks = NoOfPartsPerProcess.ForLoop(i => new List<int>((int)Math.Ceiling(1.1 * JComp / NoOfPartsPerProcess))).ToList();
                    for (int j = 0; j < JComp; j++) {
                        _Blocks[part[j]].Add(j);
                    }
                    
                    for(int iB = 0; iB < _Blocks.Count; iB++) {
                        if(_Blocks[iB].Count <= 0) {
                            _Blocks.RemoveAt(iB);
                            iB--;
                        }
                    }

                    if(_Blocks.Count < NoOfPartsPerProcess)
                        Console.WriteLine("METIS WARNING: requested " + NoOfPartsPerProcess + " blocks, but got " + _Blocks.Count);

                    return _Blocks.ToArray();
                }
            }


            /// <summary>
            /// %
            /// </summary>
            internal override int GetNoOfBlocks(MultigridOperator op) {
                return NoOfPartsPerProcess;
            }
        }


        public class SimpleBlocking : BlockingStrategy {

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
                        var temp = (int)Math.Ceiling(part.Length / (double)NoOfPartsPerProcess);
                        int[] localBlocks = new int[NoOfPartsPerProcess];
                        int sum = 0;
                        for (int i = 0; i < (NoOfPartsPerProcess - 1); i++) {
                            localBlocks[i] = temp;
                            sum += temp;
                        }
                        localBlocks[NoOfPartsPerProcess - 1] = (part.Length - sum);

                        sum = 0;
                        for (int blkIdx = 0; blkIdx < NoOfPartsPerProcess; blkIdx++) {
                            for (int i = 0; i < localBlocks[blkIdx]; i++) {
                                part.SetValue(blkIdx, i + sum);
                            }
                            sum += localBlocks[blkIdx];
                        }

                    } else {
                        part.SetAll(0);
                    }
                }

                var _Blocks = NoOfPartsPerProcess.ForLoop(i => new List<int>((int)Math.Ceiling(1.1 * JComp / NoOfPartsPerProcess)));
                for (int j = 0; j < JComp; j++) {
                    _Blocks[part[j]].Add(j);

                }

                return _Blocks;
            }

            internal override int GetNoOfBlocks(MultigridOperator op) {
                return NoOfPartsPerProcess;
            }
        }

        /// <summary>
        /// Strategy for finding the Schwarz blocks.
        /// </summary>
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
                if (value == true)
                    throw new NotSupportedException("Only supported in DEBUG mode.");
#endif
            }
        }

        /// <summary>
        /// turn P-multigrid for block solvers on/off
        /// </summary>
        public bool UsePMGinBlocks {
            get;
            set;
        }

        /// <summary>
        /// DG degree at low order sub-blocks; If p-multi-grid is used (<see cref="UsePMGinBlocks"/>), 
        /// this degree is the boundary which divides into low order and high order blocks.
        /// </summary>
        public int pLow = 1;


        /// <summary>
        /// ~
        /// </summary>
        public void Init(MultigridOperator op) {
            using (new FuncTrace()) {
                if (UsePMGinBlocks) {
                    Console.WriteLine("Schwarz: pmg is used as blocksolve");
                }

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
                    if (this.m_BlockingStrategy.GetNoOfBlocks(op) != this.blockSolvers.Count()) {
                        throw new ArgumentException("Blocking, unable to re-use");
                    }


                    return;
                }

                var Mop = op.OperatorMatrix;
                var MgMap = op.Mapping;
                
                this.m_MgOp = op;
                int myMpiRank = MgMap.MpiRank;
                int myMpisize = MgMap.MpiSize;
                int D = m_MgOp.GridData.SpatialDimension;

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

                foreach(var b in _Blocks) {
                    if (b.Count <= 0)
                        throw new ArithmeticException("Empty Schwarz-Block found");
                }
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
                    /*
                    // hack for visualization of blocks, pt. 1/2
                    double[] blockIndex = new double[ag.iLogicalCells.NoOfLocalUpdatedCells];
                    for (int idx_bi = 0; idx_bi < _Blocks.Count(); idx_bi++) {
                        var bi = _Blocks.ElementAt(idx_bi);
                        foreach (int j in bi) {
                            blockIndex[j] = idx_bi;
                        };
                    }
                    //*/

                    BitArray marker = new BitArray(JComp + JGhost);

                    if (Overlap < 0)
                        throw new ArgumentException();
                    if (Overlap > 0) {
                        if (Overlap > 1 && Mop.RowPartitioning.MpiSize > 1) {
                            //throw new NotSupportedException("In MPI parallel runs, the maximum supported overlap for the Schwarz preconditioner is 1.");
                            Console.WriteLine("In MPI parallel runs, the overlap for the Schwarz preconditioner is reduced to 1 at MPI boundaries.");
                        }

                        foreach (List<int> bi in _Blocks) { // loop over blocks...
                            marker.SetAll(false); // marks all cells which are members of the block
                            foreach (int jcomp in bi)
                                marker[jcomp] = true;

                            // determine overlap regions
                            for (int k = 0; k < Overlap; k++) { // overlap sweeps
                                int Jblock = bi.Count;
                                for (int j = 0; j < Jblock; j++) { // loop over parts of block
                                    int jCell = bi[j];
                                    if (jCell < JComp) {


                                        int[] Neighs = ag.iLogicalCells.CellNeighbours[jCell];
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

                            bi.Sort();
                        }
                    }

                    BlockCells = _Blocks.Select(list => list.ToArray()).ToArray();


                    /*
                    // hack for visualization of blocks, pt. 2/2
                    double[][] ExtendedBlocks = new double[BlockCells.Length][];
                    string[] Names = new string[BlockCells.Length];
                    for(int iBlock = 0; iBlock < ExtendedBlocks.Length; iBlock++) {
                        ExtendedBlocks[iBlock] = new double[blockIndex.Length];
                        foreach (int j in BlockCells[iBlock]) {
                            ExtendedBlocks[iBlock][j] = 1.0;
                            
                        }
                        Names[iBlock] = "block" + iBlock;
                    }

                    if (op.LevelIndex == 0) {

                        var gdata = op.BaseGridProblemMapping.GridDat;
                        

                        double[][] datas = ArrayTools.Cat<double[]>(blockIndex, ExtendedBlocks);
                        string[] names = ArrayTools.Cat<string>("indices", Names);
                        DGField[] vizData = new DGField[datas.Length];

                        for(int i = 0; i < vizData.Length; i++) {
                            vizData[i] = new SinglePhaseField(new Basis(gdata, 0), names[i]);

                            for(int j = 0; j < datas[i].Length; j++)
                                vizData[i].SetMeanValue(j, datas[i][j]);
                        }

                        Tecplot.Tecplot.PlotFields(vizData, "blox", 0.0, 0);
                    }
                    //*/
                }


                // convert cell blocks to DOF blocks
                // =================================

                /*
                {
                    Console.WriteLine("test code");
                    Mop = Mop.CloneAs();
                    Mop.Clear();
                    int Jup = MgMap.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;

                    AggregationGridBasis[] BS = MgMap.AggBasis;
                    for(int j = 0; j < Jup; j++) {
                        int NoOfSpc = BS[0].GetNoOfSpecies(j);
                        int p = MgMap.DgDegree[0];
                        int Np = BS[0].GetLength(j, p) / NoOfSpc;


                        for (int iSpc = 0; iSpc < NoOfSpc; iSpc++) {
                            for (int jSpc = 0; jSpc < NoOfSpc; jSpc++) {
                                for (int n = 0; n < Np; n++) {
                                    for (int m = 0; m < Np; m++) {
                                        double marker = (j + 1) + (iSpc + 1) * 0.1 + (jSpc + 1) * 0.01 + n * 0.001 + m * 0.0001;

                                        int irow = MgMap.LocalUniqueIndex(0, j, n + iSpc * Np);
                                        int icol = MgMap.LocalUniqueIndex(0, j, m + jSpc * Np);
                                        Mop[irow, icol] = marker;                                        
                                    }
                                }
                            }
                        }
                    }
                }
                */
                // Index lists and sub-blocking for the Schwarz blocks:
                List<int>[] BlkIdx_gI_lR; //  (global) indices in the local range 
                List<int>[] BlkIdx_gI_eR; //  (global) indices of external rows and columns
                List<int>[] TempRowIdx_gI; // (global) row indices into the temporary matrix
                List<int>[] BlkIdx_lI_eR; //  (local)  indices of external rows and columns
                List<int>[] LocalBlocks_i0, LocalBlocks_N; // sub-blocking of the Schwarz-Blocks.

                // Index lists and sub-blocking specifically if P-multigrid within Schwarz blocks is used: 
                // Rem.: sub-block (within Schwarz-blocks) == DG cells
                List<int>[] Blk_LoModes; //  indices of the low-order modes
                List<int>[] Blk_HiModes; //  indices of the high-order modes
                List<int>[] Blk_i0LoModes, Blk_NLoModes; // Blocking for low-order modes
                List<int>[] Blk_i0HiModes, Blk_NHiModes; // Blocking for high-order modes
                List<int>[] Blk_NspcHiModes; // number of species in sub-blocks/cells


                // for matrix 'ExternalRowsTemp': which rows of 'Mop' are required locally
                List<int> ExternalRowsIndices, ExternalRows_BlockI0, ExternalRows_BlockN;
                {
                    int Jup = MgMap.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                    int Jgh = MgMap.AggGrid.iLogicalCells.NoOfExternalCells;

                    AggregationGridBasis[] BS = MgMap.AggBasis;
                    int NoOfVars = BS.Length;
                    int[] Degrees = MgMap.DgDegree;


                    int LocalizedBlockCounter = 0;

                    BlkIdx_gI_lR = NoOfSchwzBlocks.ForLoop(iPart => new List<int>(BlockCells[iPart].Length * MgMap.MaximalLength));
                    BlkIdx_gI_eR = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                    LocalBlocks_i0 = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                    LocalBlocks_N = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());

                    TempRowIdx_gI = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                    BlkIdx_lI_eR = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());

                    if (this.UsePMGinBlocks) {
                        Blk_LoModes = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                        Blk_HiModes = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                        Blk_i0LoModes = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                        Blk_NLoModes = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                        Blk_i0HiModes = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                        Blk_NHiModes = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                        Blk_NspcHiModes = NoOfSchwzBlocks.ForLoop(iPart => new List<int>());
                    } else {
                        Blk_LoModes = null;
                        Blk_HiModes = null;
                        Blk_i0LoModes = null;
                        Blk_NLoModes = null;
                        Blk_i0HiModes = null;
                        Blk_NHiModes = null;
                        Blk_NspcHiModes = null;
                    }

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
                        int loBlockCounter = 0;
                        int hiBlockCounter = 0;

                        for (int jblk = 0; jblk < Jblock; jblk++) { // loop over cells in blocks...
                            int j = bc[jblk];
                            int Nj = MgMap.GetLength(j);


                            if (j < Jup) {
                                // locally updated cell
                                int i0 = MgMap.GlobalUniqueIndex(0, j, 0);

                                for (int n = 0; n < Nj; n++) {
                                    biI.Add(i0 + n);
                                }
                            } else {
                                // external cell
                                int i0E = MgMap.GlobalUniqueIndex(0, j, 0); // 
                                int i0L = MgMap.LocalUniqueIndex(0, j, 0); // 
                                                                           //LEBi0.Add(LocalizedBlockCounter);
                                                                           //LEBn.Add(N);

                                ExternalRows_BlockI0.Add(ExternalRowsIndices.Count);
                                ExternalRows_BlockN.Add(Nj);

                                for (int n = 0; n < Nj; n++) {
                                    biE.Add(i0E + n);
                                    ExternalRowsIndices.Add(i0E + n);
                                    l1.Add(LocalizedBlockCounter + n);
                                    l2.Add(i0L + n);
                                    Debug.Assert(Mop._RowPartitioning.FindProcess(i0E + n) != myMpiRank);
                                }
                                LocalizedBlockCounter += Nj;

                            }

                            Debug.Assert(ExternalRows_BlockI0.Count == ExternalRows_BlockN.Count);
                            Debug.Assert(ExternalRows_BlockI0.Count <= 1 ||
                                ExternalRows_BlockI0[ExternalRows_BlockI0.Count - 1] == ExternalRows_BlockI0[ExternalRows_BlockI0.Count - 2] + ExternalRows_BlockN[ExternalRows_BlockI0.Count - 2]);



                            LBBi0.Add(anotherCounter);
                            LBBN.Add(Nj);


                            if (UsePMGinBlocks) {
                                int cellOffset = anotherCounter;
                                var PMGloBlock = Blk_LoModes[iPart];
                                var PMGhiBlock = Blk_HiModes[iPart];
                                var PMGi0LoModes = Blk_i0LoModes[iPart];
                                var PMGNLoModes = Blk_NLoModes[iPart];
                                var PMGi0HiModes = Blk_i0HiModes[iPart];
                                var PMGNHiModes = Blk_NHiModes[iPart];
                                var PMGNspcHiModes = Blk_NspcHiModes[iPart];

                                PMGi0LoModes.Add(loBlockCounter);
                                int NpLoTot = 0;
                                for (int iVar = 0; iVar < NoOfVars; iVar++) {
                                    int p = Degrees[iVar];
                                    int NoOfSpc = BS[iVar].GetNoOfSpecies(j);
                                    int Np = BS[iVar].GetLength(j, p);

                                    int _pLow = this.pLow;
                                    if (iVar == D)
                                        _pLow -= 1; // Quick hack for Stokes systems

                                    int NpLo = BS[iVar].GetLength(j, _pLow);
                                    NpLoTot += NpLo;

                                    int NpBase = Np / NoOfSpc; // DOFs in cell *per species*
                                    int NpBaseLo = NpLo / NoOfSpc; // DOFs in cell *per species* for low-order modes

                                    PMGi0HiModes.Add(hiBlockCounter + NpBaseLo);
                                    PMGNHiModes.Add(Np - NpLo);
                                    PMGNspcHiModes.Add(NoOfSpc);

                                    for (int iSpc = 0; iSpc < NoOfSpc; iSpc++) {
                                        int n = 0;
                                        for (; n < NpBaseLo; n++) { // loop over low-order modes
                                            PMGloBlock.Add(cellOffset + n);
                                        }
                                        for (; n < NpBase; n++) { // loop over high-order modes
                                            PMGhiBlock.Add(cellOffset + n);
                                        }

                                        cellOffset += NpBase;
                                    }

                                    hiBlockCounter += Np;
                                }
                                loBlockCounter += NpLoTot;
                                PMGNLoModes.Add(NpLoTot);

                                Debug.Assert(cellOffset - anotherCounter == Nj);
                            }
                            anotherCounter += Nj;
                        }

                        if (UsePMGinBlocks) {
                            Debug.Assert(Blk_i0LoModes[iPart].Count == Jblock);
                            Debug.Assert(Blk_NLoModes[iPart].Count == Jblock);
                            Debug.Assert(Blk_i0HiModes[iPart].Count == Jblock * NoOfVars);
                            Debug.Assert(Blk_NHiModes[iPart].Count == Jblock * NoOfVars);
                            Debug.Assert(Blk_NspcHiModes[iPart].Count == Jblock * NoOfVars);
                        }
                    }
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
                    int jMax = m_MgOp.Mapping.AggGrid.iLogicalCells.Count - 1;
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
                    //  set   Perm[l,m] = I;
                    //  then  ExternalRowsTemp[l,-]  =  Mop[m,-]
                    //
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

                BlockMsrMatrix[] LoBlocks = new BlockMsrMatrix[NoOfSchwzBlocks];
                {
                    blockSolvers = new ISparseSolver[NoOfSchwzBlocks];
                    if (UsePMGinBlocks) {
                        BlockMatrices = new BlockMsrMatrix[NoOfSchwzBlocks];
                        //BlockMatrices = new ilPSP.LinSolvers.monkey.CPU.RefMatrix[NoOfSchwzBlocks];
                        PmgBlock_HiModeSolvers = new MultidimensionalArray[NoOfSchwzBlocks][];
                    }
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

                        int BlockSize = bi.Count + l1.Count;
                        if (BlockSize <= 0)
                            throw new ArithmeticException("Empty Schwarz-Block found");

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
                            for (int i = 0; i < l1L; i++) {
                                Debug.Assert(l1[i] >= 0);
                                Debug.Assert(l1[i] < ExternalRowsTemp._RowPartitioning.LocalLength);
                                l1[i] += offset;
                                Debug.Assert(ExternalRowsTemp._RowPartitioning.IsInLocalRange(l1[i]));
                            }
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
                        if (m_MatlabParalellizationCheck) {
                            Blocks.Add(Block);
                        }
#endif
                        blockSolvers[iPart] = new PARDISOSolver() {
                            CacheFactorization = true,
                            UseDoublePrecision = false
                        };
                        

                        if (UsePMGinBlocks) {
                            // define low-order solver
                            var lo_block = new BlockPartitioning(Blk_LoModes[iPart].Count, Blk_i0LoModes[iPart], Blk_NLoModes[iPart], csMPI.Raw._COMM.SELF);
                            var LoBlock = new BlockMsrMatrix(lo_block, lo_block);
                            Block.AccSubMatrixTo(1.0, LoBlock, Blk_LoModes[iPart], default(int[]), Blk_LoModes[iPart], default(int[]));
                            blockSolvers[iPart].DefineMatrix(LoBlock);
                            LoBlocks[iPart] = LoBlock;

                            // record matrix
                            BlockMatrices[iPart] = Block;
                            //BlockMatrices[iPart] = new ilPSP.LinSolvers.monkey.CPU.RefMatrix(Block.ToMsrMatrix());

                            // define high-oder sub-block solvers
                            var Nhi = Blk_NHiModes[iPart];
                            var i0Hi = Blk_i0HiModes[iPart];
                            var NspcHi = Blk_NspcHiModes[iPart];
                            var HiModes = Blk_HiModes[iPart];

                            Debug.Assert(Nhi.Count == i0Hi.Count);
                            int NoOfSubblocks = Nhi.Count;

                            PmgBlock_HiModeSolvers[iPart] = new MultidimensionalArray[NoOfSubblocks];
                            var HiModeSolvers = PmgBlock_HiModeSolvers[iPart];

                            int ptrHiModes = 0;
                            for (int iSubBlock = 0; iSubBlock < NoOfSubblocks; iSubBlock++) {
                                int Np_hi = Nhi[iSubBlock];
                                int i0 = i0Hi[iSubBlock];
                                int NoOfSpc = NspcHi[iSubBlock];
                                int Np_hiBase = Np_hi / NoOfSpc;
                                Debug.Assert(i0 == HiModes[ptrHiModes]);

                                HiModeSolvers[iSubBlock] = MultidimensionalArray.Create(Np_hi, Np_hi);

                                for (int iSpcRow = 0; iSpcRow < NoOfSpc; iSpcRow++) {
                                    int i0_spc = HiModes[ptrHiModes + Np_hiBase * iSpcRow];

                                    for (int iSpcCol = 0; iSpcCol < NoOfSpc; iSpcCol++) {
                                    //{ int iSpcCol = iSpcRow;
                                        int j0_spc = HiModes[ptrHiModes + Np_hiBase * iSpcCol];
#if DEBUG
                                        for(int n = 0; n < Np_hiBase; n++) {
                                            Debug.Assert(HiModes[ptrHiModes + Np_hiBase * iSpcRow + n] - i0_spc == n);
                                            Debug.Assert(HiModes[ptrHiModes + Np_hiBase * iSpcCol + n] - j0_spc == n);
                                        }
#endif
                                        var HiModeSolver_spc = HiModeSolvers[iSubBlock].ExtractSubArrayShallow(new[] { iSpcRow * Np_hiBase, iSpcCol * Np_hiBase },
                                            new[] { (iSpcRow + 1) * Np_hiBase - 1, (iSpcCol + 1) * Np_hiBase - 1 });

                                        Block.ReadBlock(i0_spc, j0_spc, HiModeSolver_spc);

                                    }
                                }
                                ptrHiModes += Np_hi;

                                HiModeSolvers[iSubBlock].Invert();
                            }

                        } else {
                            blockSolvers[iPart].DefineMatrix(Block);
                        }
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


                    if (UsePMGinBlocks) {
                        this.PmgBlock_LoModes = NoOfSchwzBlocks.ForLoop(iblk => Blk_LoModes[iblk].ToArray());
                        this.PmgBlock_HiModes = NoOfSchwzBlocks.ForLoop(iblk => Blk_HiModes[iblk].ToArray());
                    } else {
                        this.PmgBlock_LoModes = null;
                        this.PmgBlock_HiModes = null;
                    }
                }


                if (CoarseSolver != null) {
                    CoarseSolver.Init(op.CoarserLevel);
                }

                // solution scaling in overlapping regions
                // =======================================

                if (Overlap > 0 && EnableOverlapScaling) {
                    int LocalLength = MgMap.LocalLength;
                    this.SolutionScaling = new double[LocalLength];
                    var SolScale = this.SolutionScaling;

                    var XExchange = new MPIexchangeInverse<double[]>(this.m_MgOp.Mapping, SolScale);

                    for (int iPart = 0; iPart < NoOfSchwzBlocks; iPart++) {
                        int[] ci = BlockIndices_Local[iPart];
                        int[] ciE = BlockIndices_External[iPart];
                        int L = ci.Length;
                        int Le = ciE != null ? ciE.Length : 0;


                        // accumulate block solution 'xi' to global solution 'X'
                        //X.AccV(1.0, xi, ci, default(int[]));
                        for (int l = 0; l < L; l++) {
                            SolScale[ci[l]] += 1.0;
                        }

                        if (ciE != null && ciE.Length > 0) {
                            //XExchange.Vector_Ext.AccV(1.0, xi, ciE, default(int[]), acc_index_shift: (-LocLength), b_index_shift: ci.Length);
                            for (int l = 0; l < Le; l++) {
                                XExchange.Vector_Ext[ciE[l] - LocalLength] += 1.0;
                            }
                        }
                    }

                    XExchange.TransceiveStartImReturn();
                    XExchange.TransceiveFinish(1.0);

                    for (int l = 0; l < LocalLength; l++) {
                        SolScale[l] = 1.0 / SolScale[l];
                    }

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

                //Mop.SaveToTextFileSparse("Mop.txt");
                //for(int iBlock = 0; iBlock < NoOfSchwzBlocks; iBlock++) {
                //    this.BlockMatrices[iBlock].SaveToTextFileSparse("Block" + iBlock + ".txt");
                //    LoBlocks[iBlock].SaveToTextFileSparse("LoBlock" + iBlock + ".txt");
                //    for(int i = 0; i < PmgBlock_HiModeSolvers[iBlock].Length; i++) {
                //        PmgBlock_HiModeSolvers[iBlock][i].SaveToTextFile("HiBlock" + iBlock + "-" + i + ".txt");
                //    }
                //}



            }
        }

        /// <summary>
        /// scaling of blocks in the overlapping regions (<see cref="EnableOverlapScaling"/>).
        /// </summary>
        double[] SolutionScaling;

        /// <summary>
        /// If <see cref="Overlap"/> > 0, the solution, in cells which are covered by multiple blocks, 
        /// is scaled in the overlapping region by one over the multiplicity.
        /// This option might be useful in some applications but may also fail in others:
        /// - seems to **fail** e.g. if this is used as a preconditioner for PCG (<see cref="SoftPCG"/>)
        /// - improves number of iterations if used e.g. as a smoother for <see cref="OrthonormalizationMultigrid"/>
        /// </summary>
        public bool EnableOverlapScaling = false;

        /// <summary>
        /// the full operator matrix
        /// </summary>
        BlockMsrMatrix MtxFull {
            get {
                return m_MgOp.OperatorMatrix;
            }
        }

        /// <summary>
        /// Linear solver for each block
        /// </summary>
        ISparseSolver[] blockSolvers;

        /// <summary>
        /// In the case of P-multigrid for each level (<see cref="UsePMGinBlocks"/>), the matrices for the block
        /// - index: Schwarz block
        /// - content: matrix 
        /// </summary>
        BlockMsrMatrix[] BlockMatrices;
        //ilPSP.LinSolvers.monkey.CPU.RefMatrix[] BlockMatrices;

        /// <summary>
        /// List of low-order modes local indices in block (only used for <see cref="UsePMGinBlocks"/>).
        /// - 1st index: Schwarz block
        /// - 2nd index: enumeration; list of block-local indices which belong to low-order modes
        /// </summary>
        int[][] PmgBlock_LoModes;

        /// <summary>
        /// List of high-order modes local indices in block (only used for <see cref="UsePMGinBlocks"/>).
        /// - 1st index: block
        /// - 2nd index: enumeration; list of block-local indices which belong to high-order modes
        /// </summary>
        int[][] PmgBlock_HiModes;

        // <summary>
        // - 1st index: Schwarz block
        // - 2nd index: sub-block within the respective Schwarz block (there is one sub-block for each cell and each variable)
        // - content: the offset of the respective sub-block within the Schwarz-block
        // </summary>
        //int[][] PmgBlock_HiModeSolverI0;

        /// <summary>
        /// Cell-local solvers for the high-order modes 
        /// - 1st index: Schwarz block
        /// - 2nd index: sub-block within the respective Schwarz block (there is one sub-block for each cell and each variable)
        /// - content:some matrix inverse
        /// </summary>
        MultidimensionalArray[][] PmgBlock_HiModeSolvers;


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


        /*
        class PMGsolver : ISparseSolver {

            MsrMatrix Mtx;

            public void DefineMatrix(IMutableMatrixEx M) {
                Mtx = M.ToMsrMatrix();
            }

            public void Dispose() {
                Mtx = null;
            }

            public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
                where Tunknowns : IList<double>
                where Trhs : IList<double> {
                var StartTime = DateTime.Now;

                Mtx.SolveMATLAB(x, rhs);

                return new SolverResult() {
                    Converged = true,
                    NoOfIterations = 1,
                    RunTime = DateTime.Now - StartTime
                };
            }
        }
        */



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

        /// <summary>
        /// ~
        /// </summary>
        public Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }


        /// <summary>
        /// The fixed number of iteration on this level
        /// </summary>
        public int FixedNoOfIterations = 1;
              

        /// <summary>
        /// Coarse-grid correction
        /// </summary>
        public ISolverSmootherTemplate CoarseSolver;


        public bool CoarseSolverIsMultiplicative = true;

        int NoIter = 0;


        //static void SingleFilter(double[] V) {
        //    for(int i = 0; i < V.Length; i++) {
        //        float vi = (float)(V[i]);
        //        V[i] = vi;
        //    }
        //}


        /// <summary>
        /// ~
        /// </summary>
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

                for (int iIter = 0; iIter < FixedNoOfIterations; iIter++) {
                    this.NoIter++;

                    Res.SetV(B);
                    //Console.WriteLine("norm on swz entry: " + X.L2Norm());
                    this.MtxFull.SpMV(-1.0, X, 1.0, Res);

                    IterationCallback?.Invoke(iIter, X.ToArray(), Res.CloneAs(), this.m_MgOp);

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
                            //chackScheiss_iPart.AccV(1.0, chackScheiss, default(int[]), ci);
                            //chackScheiss_iPart.SaveToTextFile("chackScheiss" + iPart + ".txt");
                            if (ciE != null && ciE.Length > 0)
                                bi.AccV(1.0, ResExchange.Vector_Ext, default(int[]), ciE, acc_index_shift: ci.Length, b_index_shift: (-LocLength));

                            if (UsePMGinBlocks) {
                                // +++++++++++++++++++++++++++++++++
                                // P-multigrid in each Schwarz block
                                // +++++++++++++++++++++++++++++++++

                                // solve the low-order system
                                int[] ciLo = PmgBlock_LoModes[iPart];
                                int Llo = ciLo.Length;

                                double[] biLo = new double[Llo];
                                biLo.AccV(1.0, bi, default(int[]), ciLo);
                                double[] xiLo = new double[Llo];
                                try {
                                    blockSolvers[iPart].Solve(xiLo, biLo);
                                } catch (ArithmeticException ae) {
                                    Console.Error.WriteLine(ae.Message);
                                    throw ae;
                                }
                                xi.AccV(1.0, xiLo, ciLo, default(int[]));

                                {
                                    // re-evaluate the residual
                                    this.BlockMatrices[iPart].SpMV(-1.0, xi, 1.0, bi);

                                    // solve the high-order system
                                    int[] ciHi = PmgBlock_HiModes[iPart];
                                    var HiModeSolvers = PmgBlock_HiModeSolvers[iPart];
                                    int NoCells = HiModeSolvers.Length;

                                    double[] xiHi = null;
                                    double[] biHi = null;

                                    int ptr_CiHi = 0;
                                    for (int j = 0; j < NoCells; j++) {
                                        var HiModeSolver = HiModeSolvers[j];
                                        int Np = HiModeSolver.NoOfRows;


                                        if (xiHi == null || xiHi.Length != Np)
                                            xiHi = new double[Np];
                                        if (biHi == null || biHi.Length != Np)
                                            biHi = new double[Np];

                                        //double[] chackScheiss_iPart_j = new double[Np];
                                        for (int n = 0; n < Np; n++) {
                                            biHi[n] = bi[ciHi[ptr_CiHi + n]];
                                            //chackScheiss_iPart_j[n] = chackScheiss_iPart[ciHi[ptr_CiHi + n]];
                                        }
                                        //chackScheiss_iPart_j.SaveToTextFile("chackScheiss" + iPart + "-" + j + ".txt");


                                        HiModeSolver.GEMV(1.0, biHi, 0.0, xiHi);

                                        for (int n = 0; n < Np; n++) {
                                            xi[ciHi[ptr_CiHi + n]] = xiHi[n];
                                        }
                                        ptr_CiHi += Np;
                                    }
                                    Debug.Assert(ptr_CiHi == ciHi.Length);
                                }

                            } else {
                                // ++++++++++++++++++++++++++++++
                                // use block solver for all modes
                                // ++++++++++++++++++++++++++++++

                                blockSolvers[iPart].Solve(xi, bi);
                                //SingleFilter(xi);
                            }

                            // accumulate block solution 'xi' to global solution 'X'
                            X.AccV(1.0, xi, ci, default(int[]));
                            if (ciE != null && ciE.Length > 0)
                                XExchange.Vector_Ext.AccV(1.0, xi, ciE, default(int[]), acc_index_shift: (-LocLength), b_index_shift: ci.Length);
                        }
                    }


                    if (Overlap > 0 && EnableOverlapScaling) {
                        // block solutions stored on *external* indices will be accumulated on other processors.
                        XExchange.TransceiveStartImReturn();
                        XExchange.TransceiveFinish(1.0);

                        if (iIter < FixedNoOfIterations - 1)
                            XExchange.Vector_Ext.ClearEntries();

                        var SolScale = this.SolutionScaling;
                        for (int l = 0; l < LocLength; l++) {
                            X[l] *= SolScale[l];
                        }

                        //SolScale.SaveToTextFile("SolScale.txt");
                    }

                    /*
                     * Note (fk, 04aug19:)
                     * Doing the PMG after the low-p Schawrz blocks does not work very well
                     * - it would be faster, since one global SpMV replaces multiple SpMVs in blocks
                     * - but somehow,  oscillations at block boundaries cause oscillations and slower convergence in general
                    if (vi == null)
                        vi = new MGViz(this.m_MgOp);
                    var Xb4 = X.ToArray();
                    
                    double ResiNormB4 = 0.0;
                    if (UsePMGinBlocks) {
                        int Jloc = this.m_MgOp.Mapping.LocalNoOfBlocks;
                        bool[] hackMarker = new bool[Jloc];

                        using (new BlockTrace("himode_smoother_level" + this.m_MgOp.LevelIndex, tr)) {
                            Res.SetV(B);
                            this.MtxFull.SpMV(-1.0, X, 1.0, Res);
                            ResiNormB4 = Res.L2Norm();

                            for (int iPart = 0; iPart < NoParts; iPart++) {

                                int[] BlockCells_iPart = this.BlockCells[iPart];
                                
                                int[] ci = BlockIndices_Local[iPart];
                                int[] ciE = BlockIndices_External[iPart];
                                int L = ci.Length;
                                if (ciE != null)
                                    L += ciE.Length;

                                double[] bi = new double[L];
                                double[] xi = new double[L];

                                // extract block part of residual
                                bi.AccV(1.0, Res, default(int[]), ci);
                                //if (ciE != null && ciE.Length > 0)
                                //    bi.AccV(1.0, ResExchange.Vector_Ext, default(int[]), ciE, acc_index_shift: ci.Length, b_index_shift: (-LocLength));

                                // solve the high-order system
                                int[] ciHi = PmgBlock_HiModes[iPart];
                                var HiModeSolvers = PmgBlock_HiModeSolvers[iPart];
                                int NoCells = HiModeSolvers.Length;
                                Debug.Assert(NoCells == BlockCells_iPart.Length);

                                double[] xiHi = null;
                                double[] biHi = null;

                                int ptr_CiHi = 0;
                                for (int j = 0; j < NoCells; j++) {
                                    int jCell = BlockCells_iPart[j];
                                    var HiModeSolver = HiModeSolvers[j];
                                    int Np = HiModeSolver.NoOfRows;

                                    if (jCell < Jloc && hackMarker[jCell] == false) {


                                        if (xiHi == null || xiHi.Length != Np)
                                            xiHi = new double[Np];
                                        if (biHi == null || biHi.Length != Np)
                                            biHi = new double[Np];

                                        for (int n = 0; n < Np; n++) {
                                            biHi[n] = bi[ciHi[ptr_CiHi + n]];
                                        }

                                        HiModeSolver.GEMV(1.0, biHi, 0.0, xiHi);

                                        for (int n = 0; n < Np; n++) {
                                            xi[ciHi[ptr_CiHi + n]] = xiHi[n];
                                        }
                                    }
                                    ptr_CiHi += Np;
                                }
                                Debug.Assert(ptr_CiHi == ciHi.Length);

                                // accumulate block solution 'xi' to global solution 'X'
                                X.AccV(-1.0, xi, ci, default(int[]));
                            }
                        }
                        

                    }//*/

                } // end loop Schwarz iterations

            } // end FuncTrace
        }


        /// <summary>
        /// ~
        /// </summary>
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
                return this.NoIter;
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

        /// <summary>
        /// Forget any factorization stored for blocks.
        /// </summary>
        public void DisposeBlocks() {
            foreach(var b in this.blockSolvers) {
                b.Dispose();
            }
        }

        /// <summary>
        /// ~
        /// </summary>
        public object Clone() {
            Schwarz Clone = new Schwarz();
            if (this.m_BlockingStrategy is METISBlockingStrategy) {
                Clone.m_BlockingStrategy = new METISBlockingStrategy() {
                    NoOfPartsPerProcess = ((METISBlockingStrategy)this.m_BlockingStrategy).NoOfPartsPerProcess
                };
            } else if (this.m_BlockingStrategy is SimpleBlocking) {
                Clone.m_BlockingStrategy = new SimpleBlocking() {
                    NoOfPartsPerProcess = ((SimpleBlocking)this.m_BlockingStrategy).NoOfPartsPerProcess
                };
            } else if (this.m_BlockingStrategy is MultigridBlocks) {
                Clone.m_BlockingStrategy = new MultigridBlocks() {
                    Depth = ((MultigridBlocks)this.m_BlockingStrategy).Depth
                };
            } else {
                throw new NotImplementedException();
            }

            if (Clone.m_BlockingStrategy.GetType() != this.m_BlockingStrategy.GetType())
                throw new NotImplementedException();


            Clone.FixedNoOfIterations = this.FixedNoOfIterations;
            Clone.m_Overlap = this.m_Overlap;
            Clone.IterationCallback = this.IterationCallback;
            if (this.CoarseSolver != null)
                Clone.CoarseSolver = this.CoarseSolver.CloneAs();
            return Clone;
        }

        /// <summary>
        /// ~
        /// </summary>
        public long UsedMem {
            get {
                long s = 0;

                if(BlockMatrices != null) {
                    s += BlockMatrices.Sum(mtx => mtx.UsedMemory);
                }

                if(PmgBlock_HiModeSolvers != null) {
                    s += PmgBlock_HiModeSolvers.Sum(mdaa => mdaa.Sum(mda => (long)mda.Length * sizeof(double)));
                }
                
                return s;
            }
        }
    }
}


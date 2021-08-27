﻿/* =======================================================================
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
using ilPSP.LinSolvers.MUMPS;

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
            /// Number of parts/additive Schwarz blocks on current MPI process (can be different on other processors)
            /// </summary>
            public int NoOfPartsOnCurrentProcess = 4;

            internal override IEnumerable<List<int>> GetBlocking(MultigridOperator op) {

                if(cache != null) {
                    return cache.Select(orgList => new List<int>(orgList)).ToArray();
                }
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
                    if (NoOfPartsOnCurrentProcess > 1) {
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
                                ref NoOfPartsOnCurrentProcess,
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
                    List<List<int>> _Blocks = NoOfPartsOnCurrentProcess.ForLoop(i => new List<int>((int)Math.Ceiling(1.1 * JComp / NoOfPartsOnCurrentProcess))).ToList();
                    for (int j = 0; j < JComp; j++) {
                        _Blocks[part[j]].Add(j);
                    }
                    
                    for(int iB = 0; iB < _Blocks.Count; iB++) {
                        if(_Blocks[iB].Count <= 0) {
                            _Blocks.RemoveAt(iB);
                            iB--;
                        }
                    }

                    if(_Blocks.Count < NoOfPartsOnCurrentProcess)
                        Console.WriteLine("METIS WARNING: requested " + NoOfPartsOnCurrentProcess + " blocks, but got " + _Blocks.Count);

                    cache = _Blocks.ToArray();
                    Console.WriteLine("MetisBlocking Testcode active !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                    return cache.Select(orgList => new List<int>(orgList)).ToArray();
                }
            }

            List<int>[] cache;


            /// <summary>
            /// %
            /// </summary>
            internal override int GetNoOfBlocks(MultigridOperator op) {
                return NoOfPartsOnCurrentProcess;
            }
        }


        /// <summary>
        /// 
        /// </summary>
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

        // Hack the hack, if pressure is equal order ...
        // Only viable, if p-two-grid used
        public bool EqualOrder = false;

        /// <summary>
        /// turn P-multigrid for block solvers on/off.
        /// Not recommended: This may cause bad convergence in the presence of pressure.
        /// </summary>
        public bool UsePMGinBlocks = false;
        
        private bool AnyHighOrderTerms {
            get {
                Debug.Assert(m_MgOp != null, "there is no matrix given yet!");
                return m_MgOp.Mapping.DgDegree.Any(p => p > pLow);
            }
        }

        /// <summary>
        /// The maximum order of the coarse system, which is solved by a direct solver.
        /// NOTE: there is a hack, which consideres <see cref="CoarseLowOrder"/>-1 for pressure.
        /// pressure is assumed to be the Dimension-1-th variable
        /// </summary>
        public int CoarseLowOrder {
            get { return pLow; }
            set { pLow = value; }
        }

        /// <summary>
        /// Determines, if cutcells are fully assigned (<see cref="CoarseLowOrder"/>=p) to the coarse solver; only applicable, if p-two-grid is used as block solver
        /// </summary>
        public bool CoarseSolveOfCutcells = true;

        /// <summary>
        /// DG degree at low order sub-blocks; If p-two-grid is used (<see cref="UsePMGinBlocks"/>), 
        /// this degree is the boundary which divides into low order and high order blocks.
        /// </summary>
        private int pLow = 1;

        /// <summary>
        /// ~
        /// </summary>
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
                    //if (!object.ReferenceEquals(this.MtxFull, op.OperatorMatrix)) {
                    //    BlockMsrMatrix Check = this.MtxFull.CloneAs();
                    //    Check.Acc(-1.0, op.OperatorMatrix);
                    //    if (Check.InfNorm() != 0.0) {
                    //        throw new ArgumentException("Matrix has changed, unable to re-use");
                    //    }
                    //}
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

                //if (UsePMGinBlocks && AnyHighOrderTerms)
                //    Console.WriteLine("Schwarz: pmg is used as blocksolve");
                //if (UsePMGinBlocks && !AnyHighOrderTerms)
                //    Console.WriteLine("Schwarz: Only low order terms present. Schwarz blocks will be solved direct instead of PMG");
                if (!Mop.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Row partitioning mismatch.");
                if (!Mop.ColPartition.EqualsPartition(MgMap.Partitioning))
                    throw new ArgumentException("Column partitioning mismatch.");

                var ag = MgMap.AggGrid;

                int JComp = ag.iLogicalCells.NoOfLocalUpdatedCells;
                int JGhost = ag.iLogicalCells.NoOfExternalCells;

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

                int[][] BlockCells = null;

                // extend blocks according to desired overlap
                // ==========================================
                {

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
                    } else {
                        //Console.WriteLine("Running Schwarz without overlap (level " + this.m_MgOp.LevelIndex + ")");
                    }

                    BlockCells = _Blocks.Select(list => list.ToArray()).ToArray();

                }

                // Get all the External rows at once, for performance sake!
                BlockMsrMatrix ExtRows = null;
                if(m_Overlap > 0)
                    ExtRows = BlockMask.GetAllExternalRows(m_MgOp.Mapping, m_MgOp.OperatorMatrix);
#if TEST
                ExtRows.SaveToTextFileSparseDebug("ExtRows");
#endif

                blockSolvers = new ISparseSolver[NoOfSchwzBlocks];
                BMfullBlocks = new BlockMask[NoOfSchwzBlocks];
                BlockMatrices = new BlockMsrMatrix[NoOfSchwzBlocks];
                Levelpmgsolvers = new BlockLevelPmg[NoOfSchwzBlocks];

                // Factory to provide the level pmg solver for blocks
                // cannot be merged with LevelPmg right now,
                // because sub-mapping is not implemented right now
                var PTGFactory = new BlockLevelPmgFactory(op, ExtRows) {
                    pLow = CoarseLowOrder,
                    FullSolveOfCutcells = CoarseSolveOfCutcells
                };

                for (int iPart = 0; iPart < NoOfSchwzBlocks; iPart++) { // loop over parts...
                    Debug.Assert(BlockCells != null);
                    int[] bc = BlockCells[iPart];

                    BlockMask fullMask;
                    BlockMsrMatrix fullBlock;

                    if (UsePMGinBlocks && AnyHighOrderTerms) {
                        // Init of level pmg block solvers
                        Levelpmgsolvers[iPart] = PTGFactory.CreateAndInit(bc.ToList(), out fullBlock, out fullMask);
                    } else {
                        // direct solver for blocks

                        // generates the block mask
                        var fullSel = new SubBlockSelector(op.Mapping);
                        fullSel.CellSelector(bc.ToList(), false);
                        fullMask = new BlockMask(fullSel, ExtRows);
                        fullBlock = fullMask.GetSubBlockMatrix(op.OperatorMatrix);
                        Debug.Assert(fullBlock.RowPartitioning.MPI_Comm == csMPI.Raw._COMM.SELF);

                        //blockSolvers[iPart] = new PARDISOSolver() {
                        //    CacheFactorization = true,
                        //    UseDoublePrecision = true,
                        //    Parallelism = Parallelism.SEQ,
                        //};

                        blockSolvers[iPart] = new MUMPSSolver() {
                            Parallelism = Parallelism.SEQ,
                        };

                        // ILU nicht ratsam, viel mehr Iterationen nötig, als mit PARDISO
                        //blockSolvers[iPart] = new ilPSP.LinSolvers.HYPRE.Euclid() {
                        //    Level = 4,
                        //    Comm = csMPI.Raw._COMM.SELF
                        //};

                        blockSolvers[iPart].DefineMatrix(fullBlock);
                        BlockMatrices[iPart] = fullBlock; // just used to calculate memory consumption

                        
                    }
                    BMfullBlocks[iPart] = fullMask;
                }

                // Watchdog bomb!
                // ==============

                // multiple watchdogs to detect out-of-sync                
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();
                MPICollectiveWatchDog.Watch();


                // solution scaling in overlapping regions
                // =======================================

                if (Overlap > 0 && EnableOverlapScaling) {
                    int LocalLength = MgMap.LocalLength;
                    this.SolutionScaling = new double[LocalLength];
                    var SolScale = this.SolutionScaling;

                    var XExchange = new MPIexchangeInverse<double[]>(this.m_MgOp.Mapping, SolScale);

                    //int extLength = MgMap.AggGrid.iLogicalCells.NoOfExternalCells + MgMap.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                    //double[] sammel = new double[extLength];

                    for (int iPart = 0; iPart < NoOfSchwzBlocks; iPart++) {

                        int rows = BMfullBlocks[iPart].GetNoOfMaskedRows;
                        double[] druffdamit = rows.ForLoop<double>(i=>1.0);

                        BMfullBlocks[iPart].AccSubVec(druffdamit, XExchange.Vector_Ext, SolScale);

                    }

                    XExchange.TransceiveStartImReturn();
                    XExchange.TransceiveFinish(1.0);

                    for (int l = 0; l < LocalLength; l++) {
                        SolScale[l] = 1.0 / SolScale[l];
                    }

                }
            }
        }

        private void ModifyLowSelector(SubBlockSelector sbs, MultigridOperator op) {
            AssignXdgBlocksModification(sbs, op, true);
        }

        private void ModifyHighSelector(SubBlockSelector sbs, MultigridOperator op) {
            AssignXdgBlocksModification(sbs, op, false);
        }

        private void AssignXdgBlocksModification(SubBlockSelector sbs, MultigridOperator op, bool IsLowSelector) {
            var Filter = sbs.ModeFilter;
            Func<int, int, int, int, bool> Modification = delegate (int iCell, int iVar, int iSpec, int pDeg) {
                int NoOfSpec = op.Mapping.AggBasis[0].GetNoOfSpecies(iCell);
                if (NoOfSpec >= 2)
                    return IsLowSelector;
                else
                    return Filter(iCell, iVar, iSpec, pDeg);
            };
            sbs.ModeSelector(Modification);
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
        /// Level pmg solver, which are used instead of direct solvers,
        /// if <see cref="UsePMGinBlocks"/> is set.
        /// </summary>
        BlockLevelPmg[] Levelpmgsolvers;

        /// <summary>
        /// In the case of P-multigrid for each level (<see cref="UsePMGinBlocks"/>), the matrices for the block
        /// - index: Schwarz block
        /// - content: matrix 
        /// </summary>
        BlockMsrMatrix[] BlockMatrices;


        /// <summary>
        /// masks for the Schwarz blocks
        /// - index: Schwarz block
        /// </summary>
        BlockMask[] BMfullBlocks;

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
                    throw new ArgumentException("overlap cannot be negative");
                }
                if (value > 2) {
                    throw new ArgumentException($"overlap of {value} is not supported - maximum is 2.");
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


        /// <summary>
        /// If the <see cref="CoarseSolver"/> is specified, setting this to true
        /// causes the residual to be re-evaluated before the coarse correction; in this fashion, the Schwarz method can be written multiplicative.
        /// </summary>
        public bool CoarseSolverIsMultiplicative = true;

        int NoIter = 0;

        //private double[] Xdummy, Resdummy;

        /// <summary>
        /// ~
        /// </summary>
        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (var tr = new FuncTrace()) {
                int NoParts = this.BMfullBlocks.Length;

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
                // Reminder: even though no overlap considered, these Vectors have to be there ...
                ResExchange = new MPIexchange<double[]>(this.m_MgOp.Mapping, Res);
                XExchange = new MPIexchangeInverse<U>(this.m_MgOp.Mapping, X);

                int LocLength = m_MgOp.Mapping.LocalLength;

                for (int iIter = 0; iIter < FixedNoOfIterations; iIter++) {
                    this.NoIter++;

                    Res.SetV(B);
                    //Console.WriteLine("norm on swz entry: " + X.L2Norm());
                    this.MtxFull.SpMV(-1.0, X, 1.0, Res);

                    IterationCallback?.Invoke(iIter, X.ToArray(), Res.CloneAs(), this.m_MgOp);

                    //using (new BlockTrace("coarse_solve_level" + this.m_MgOp.LevelIndex, tr)) {
                    //    if (CoarseSolver != null) {
                    //        var XC = X.ToArray().CloneAs();
                    //        double[] bc = new double[m_MgOp.CoarserLevel.Mapping.LocalLength];// = Res.CloneAs();
                    //        m_MgOp.CoarserLevel.Restrict(Res.CloneAs(), bc);
                    //        double[] xc = new double[bc.Length];
                    //        CoarseSolver.Solve(xc, bc);
                    //        //SingleFilter(xc);
                    //        m_MgOp.CoarserLevel.Prolongate(1, XC, 1, xc);
                    //        X.AccV(1.0, XC);

                    //        if (CoarseSolverIsMultiplicative) {
                    //            Res.SetV(B);
                    //            this.MtxFull.SpMV(-1.0, X, 1.0, Res);
                    //        }
                    //    }
                    //}

                    if (Overlap > 0) {
                        ResExchange.TransceiveStartImReturn();
                        ResExchange.TransceiveFinish(0.0);
                    }


                    using (new BlockTrace("block_solve_level", tr)) {

                        for (int iPart = 0; iPart < NoParts; iPart++) {

                            var bi = BMfullBlocks[iPart].GetSubVec(ResExchange.Vector_Ext, Res);
                            double[] xi = new double[bi.Length];
                            if (UsePMGinBlocks && AnyHighOrderTerms) {
                                // +++++++++++++++++++++++++++++++++
                                // P-multigrid in each Schwarz block
                                // +++++++++++++++++++++++++++++++++

                                Levelpmgsolvers[iPart].Solve(xi, bi);

                            } else {
                                // ++++++++++++++++++++++++++++++
                                // use block solver for all modes
                                // ++++++++++++++++++++++++++++++

                                blockSolvers[iPart].Solve(xi, bi);
                                //SingleFilter(xi);
                            }

                            //xi.SaveToTextFileDebug(String.Format("x{0}",iPart));

                            // accumulate block solution 'xi' to global solution 'X'
                            BMfullBlocks[iPart].AccSubVec(xi, XExchange.Vector_Ext, X);
                        }
                    }

                    using (new BlockTrace("overlap_scaling", tr)) {
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
                        }
                    }
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

        /// <summary>
        /// <see cref="ISolverSmootherTemplate.ThisLevelIterations"/>
        /// </summary>
        public int ThisLevelIterations {
            get {
                return this.NoIter;
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
            if (this.CoarseSolver != null)
                this.CoarseSolver.ResetStat();
        }

        /// <summary>
        /// Forget any factorization stored for blocks.
        /// </summary>
        private void DisposeBlockSolver() {
            if (this.blockSolvers == null || this.blockSolvers.Count() <= 0)
                return;
            foreach (var b in this.blockSolvers) {
                if(b!=null) b.Dispose();
            }
        }

        private void DisposePMGSolvers() {
            if (UsePMGinBlocks && Levelpmgsolvers != null) {
                foreach (var solver in this.Levelpmgsolvers) {
                    solver.Dispose();
                }
                this.Levelpmgsolvers = null;
            }
        }

        /// <summary>
        /// ~
        /// </summary>
        public object Clone() {
            Schwarz Clone = new Schwarz();
            if (this.m_BlockingStrategy is METISBlockingStrategy) {
                Clone.m_BlockingStrategy = new METISBlockingStrategy() {
                    NoOfPartsOnCurrentProcess = ((METISBlockingStrategy)this.m_BlockingStrategy).NoOfPartsOnCurrentProcess
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
                if (UsePMGinBlocks && AnyHighOrderTerms) {
                    s += Levelpmgsolvers.Sum(solver => solver.UsedMem);
                } else {
                    if (BlockMatrices != null) {
                        s += BlockMatrices.Sum(mtx => mtx.UsedMemory);
                    }
                }
                return s;
            }
        }

        /// <summary>
        /// This can be used for testing MPI parallel execution
        /// </summary>
        /// <param name="BMs">Block Mask for schwarz blocks</param>
        /// <param name="Mop"></param>
        public static void MatlabDebugging(BlockMask[] BMs, MultigridOperator Mop) {

            int NoOfSchwzBlocks = BMs.Length;
            int myMpisize = Mop.Mapping.MpiSize;
            int myMpiRank = Mop.Mapping.MpiRank;
            var MopMap = Mop.Mapping;

            ilPSP.Connectors.Matlab.BatchmodeConnector matlab;
            matlab = new ilPSP.Connectors.Matlab.BatchmodeConnector();


            List<BlockMsrMatrix> Blocks = new List<BlockMsrMatrix>();
            var BlkIdx_gI_lR = NoOfSchwzBlocks.ForLoop(b => new List<long>());
            var BlkIdx_gI_eR = NoOfSchwzBlocks.ForLoop(b => new List<long>());
            int[][] BlockIndices_Local = new int[NoOfSchwzBlocks][];
            int[][] BlockIndices_External = new int[NoOfSchwzBlocks][];

            long LocalI0 = MopMap.i0;
            for (int iPart = 0; iPart < NoOfSchwzBlocks; iPart++) {
                BlkIdx_gI_lR[iPart] = BMs[iPart].GlobalIList_Internal;
                BlkIdx_gI_eR[iPart] = BMs[iPart].GlobalIList_External;
                var locallist = new List<int>();
                var extlist = new List<int>();
                foreach (int lIdx in BlkIdx_gI_lR[iPart]) {
                    locallist.Add((int)(lIdx - LocalI0));
                }
                foreach (int eIdx in BlkIdx_gI_eR[iPart]) {
                    extlist.Add((int)(eIdx - LocalI0));
                }
                BlockIndices_Local[iPart] = locallist.ToArray();
                BlockIndices_External[iPart] = extlist.ToArray();
                Blocks.Add(BMs[iPart].GetSubBlockMatrix(Mop.OperatorMatrix));
            }


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

            //matlab.PutSparseMatrix(Perm, "Perm");
            //matlab.PutSparseMatrix(ExternalRowsTemp, "ExternalRowsTemp");

            globalBlockCounter = 0;
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


            Console.WriteLine("Matlab dir: " + matlab.WorkingDirectory);

            matlab.PutSparseMatrix(Mop.OperatorMatrix, "Full");
            int GlobalNoOfBlocks = NoOfSchwzBlocks.MPISum();



            for (int iGlbBlock = 0; iGlbBlock < GlobalNoOfBlocks; iGlbBlock++) {
                matlab.Cmd("BlockErr({0} + 1, 1) = norm( Block{0} - Full( BlockIdx{0}, BlockIdx{0} ), inf );", iGlbBlock);
            }

            Random rnd = new Random(myMpiRank);
            double[] testRHS = new double[MopMap.LocalLength];
            for (int i = 0; i < testRHS.Length; i++) {
                testRHS[i] = rnd.NextDouble();
            }
            matlab.PutVector(testRHS, "testRHS");

            MPIexchange<double[]> ResExchange = new MPIexchange<double[]>(Mop.Mapping, testRHS);
            ResExchange.TransceiveStartImReturn();
            ResExchange.TransceiveFinish(0.0);

            int offset = MopMap.LocalLength;

            int g = 0;
            for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {
                    double[] SubVec;
                    if (rankCounter == myMpiRank) {
                        int LL = BlockIndices_Local[iBlock].Length;
                        int LE;
                        if (BlockIndices_External[iBlock] != null) {
                            LE = BlockIndices_External[iBlock].Length;
                        } else {
                            LE = 0;
                        }
                        int L = LL + LE;

                        SubVec = new double[L];
                        for (int i = 0; i < LL; i++) {
                            SubVec[i] = testRHS[BlockIndices_Local[iBlock][i]];
                        }
                        if (LE > 0) {
                            for (int i = 0; i < LE; i++) {
                                SubVec[i + LL] = ResExchange.Vector_Ext[BlockIndices_External[iBlock][i] - offset];
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
            MPIexchangeInverse<double[]> XXExchange = new MPIexchangeInverse<double[]>(MopMap, testX);

            g = 0;
            for (int rankCounter = 0; rankCounter < myMpisize; rankCounter++) {
                int rank_NoBlks = NoOfSchwzBlocks.MPIBroadcast(rankCounter);
                for (int iBlock = 0; iBlock < rank_NoBlks; iBlock++) {

                    if (rankCounter == myMpiRank) {
                        int LL = BlockIndices_Local[iBlock].Length;
                        int LE;
                        if (BlockIndices_External[iBlock] != null) {
                            LE = BlockIndices_External[iBlock].Length;
                        } else {
                            LE = 0;
                        }
                        int L = LL + LE;


                        for (int i = 0; i < LL; i++) {
                            testX[BlockIndices_Local[iBlock][i]] += (g + 1);
                        }
                        if (LE > 0) {
                            for (int i = 0; i < LE; i++) {
                                XXExchange.Vector_Ext[BlockIndices_External[iBlock][i] - offset] += (g + 1);
                            }
                        }
                    } else {
                        //nop
                    }

                    g++;
                }
            }
            XXExchange.TransceiveStartImReturn();
            XXExchange.TransceiveFinish(1.0);

            matlab.Cmd("testXref = zeros({0},1);", MopMap.TotalLength);
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

        public void Dispose() {
            this.DisposeBlockSolver();
            this.DisposePMGSolvers();
            this.SolutionScaling = null;
            this.BlockMatrices = null;
            this.BMfullBlocks = null;
        }

        public double UsedMemory() {
            double LScaling = this.SolutionScaling.Length * sizeof(double);
            double MemoryOfBlocks = 0.0;
            foreach (var block in BlockMatrices)
                MemoryOfBlocks+=block.UsedMemory;
            return (LScaling + MemoryOfBlocks) / (1024 * 1024);
        }
    }
}


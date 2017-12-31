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
        }

        /// <summary>
        /// Additive-Schwarz blocks which are fromed from coarser multigrid-levels.
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
                    if(value < 0) {
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
                    if(blokOp.CoarserLevel == null)
                        throw new NotSupportedException("Not enough multigrid levels set to support a depth of "+ m_Depht + ".");
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
                    foreach(int j in Blocks[iBlk]) {
                        Debug.Assert(j >= 0);
                        Debug.Assert(j < checkOnce.Length);
                        Debug.Assert(checkOnce[j] == false);
                        checkOnce[j] = true;
                    }

#endif
                }

#if DEBUG
                for(int j = 0; j < checkOnce.Length; j++) {
                    Debug.Assert(checkOnce[j] == true);
                }
#endif

                return Blocks;
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
                    foreach(int jFine in CoarseCell) {
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
            public int NoOfParts = 4;

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
                    if (NoOfParts > 1) {
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
                            ref NoOfParts,
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
                    var _Blocks = NoOfParts.ForLoop(i => new List<int>((int)Math.Ceiling(1.1 * JComp / NoOfParts)));
                    for (int j = 0; j < JComp; j++) {
                        _Blocks[part[j]].Add(j);

                    }

                    return _Blocks;
                }
            }
        }


        public BlockingStrategy m_BlockingStrategy;
        MultigridOperator m_MgOp;

        public void Init(MultigridOperator op) {
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

            // get cell blocks
            // ===============

            var _Blocks = this.m_BlockingStrategy.GetBlocking(op);
            int NoParts = _Blocks.Count();

            // test cell blocks
            // ================
#if DEBUG
            {
                // ensure that each cell is used exactly once, among all blocks
                bool[] test = new bool[ag.iLogicalCells.NoOfLocalUpdatedCells];
                foreach (var bi in _Blocks) {
                    foreach(int j in bi) { 
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
                    if(Overlap > 1 && Mop.RowPartitioning.MpiSize > 1) {
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

            List<int>[] _LocallyStoredBlockIndices; // for each Schwarz block, (global) indices in the local range 
            List<int>[] ExternalBlockIndices; //       for each Schwarz block, (global) indices of external rows and columns
            List<int>[] TempExternalRowIndices; //     for each Schwarz block, (global) indices into the temporary matrix
            List<int>[] L2; //                         for each Schwarz block, (local) indices of external rows and columns
            List<int>[] LocalBlocks_i0, LocalBlocks_N; // blocking of the Schwarz-Blocks.

            // for matrix 'ExternalRowsTemp': which rows of 'Mop' are required locally
            List<int> ExternalRowsIndices, ExternalRows_BlockI0, ExternalRows_BlockN;
            {
                int Jup = MgMap.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                int Jgh = MgMap.AggGrid.iLogicalCells.NoOfExternalCells;

                int LocalizedBlockCounter = 0;
                

                _LocallyStoredBlockIndices = NoParts.ForLoop(iPart => new List<int>(BlockCells[iPart].Length * MgMap.MaximalLength));
                ExternalBlockIndices = NoParts.ForLoop(iPart => new List<int>());
                LocalBlocks_i0 = NoParts.ForLoop(iPart => new List<int>());
                LocalBlocks_N = NoParts.ForLoop(iPart => new List<int>());

                TempExternalRowIndices = NoParts.ForLoop(iPart => new List<int>());
                L2 = NoParts.ForLoop(iPart => new List<int>());

                
                ExternalRowsIndices = new List<int>();
                ExternalRows_BlockI0 = new List<int>();
                ExternalRows_BlockN = new List<int>();

                for (int iPart = 0; iPart < NoParts; iPart++) { // loop over parts...
                    int[] bc = BlockCells[iPart];
                    var biI = _LocallyStoredBlockIndices[iPart];
                    var biE = ExternalBlockIndices[iPart];
                    var l1 = TempExternalRowIndices[iPart];
                    var l2 = L2[iPart];
                    var LBBi0 = LocalBlocks_i0[iPart];
                    var LBBN = LocalBlocks_N[iPart];
                    //var LEBidx = LocalizedExternalBlockIndices[iPart];
                    //var LEBi0 = LocalizedExternalBlockI0[iPart];
                    //var LEBn = LocalizedExternalBlockN[iPart];

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
                                //LEBidx.Add(LocalizedBlockCounter); LocalizedBlockCounter++;
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

            BlockMsrMatrix ExternalRowsTemp;
            if(myMpisize > 1) {
                //int NoOfLocalRows = _ExternalBlockIndices.Sum(L => L.Count);

                BlockPartitioning PermRow = new BlockPartitioning(ExternalRowsIndices.Count, ExternalRows_BlockI0, ExternalRows_BlockN, Mop.MPI_Comm, i0isLocal: true);

                // Remark: we use a permutation matrix for MPI-exchange of rows

                BlockMsrMatrix Perm = new BlockMsrMatrix(PermRow, Mop._RowPartitioning);
                for(int iRow = 0; iRow < ExternalRowsIndices.Count; iRow++) {
                    Perm[iRow + PermRow.i0, ExternalRowsIndices[iRow]] = 1;
                }
                
                ExternalRowsTemp = BlockMsrMatrix.Multiply(Perm, Mop);
            } else {
                ExternalRowsTemp = null;
            }


            
            
            //Console.WriteLine("Schwarz, blocking strategy {0}", this.m_BlockingStrategy.GetType().FullName);
            //double avg = 0;
            //for(int iblock  = 0; iblock < BlockIndices.Length; iblock++) {
            //    Console.WriteLine("Schwarz, blocking strategy {0}, block {0}, DOFs per block {1}", iblock, BlockIndices[iblock].Length);
            //    avg += BlockIndices[iblock].Length;
            //}
            //avg = Math.Round(avg / BlockIndices.Length);
            //Console.WriteLine("Average block size: " + avg);
            
            // create solvers
            // ==============

            {
                blockSolvers = new ISparseSolver[NoParts];

                for (int iPart = 0; iPart < NoParts; iPart++) {
                    var bi = _LocallyStoredBlockIndices[iPart];

                    int Bsz;
                    if (MgMap.MinimalLength == MgMap.MaximalLength)
                        Bsz = MgMap.MaximalLength;
                    else
                        Bsz = 1;

                    var l1 = TempExternalRowIndices[iPart];

                    //if (M.RowPartitioning.MpiSize > 1) {
                    //    int i0Proc = M.RowPartitioning.i0;
                    //    bi = bi.CloneAs();
                    //    for (int i = 0; i < bi.Length; i++) {
                    //        bi[i] += i0Proc;
                    //    }
                    //}

                    BlockPartitioning localBlocking = new BlockPartitioning(bi.Count + l1.Count, LocalBlocks_i0[iPart], LocalBlocks_N[iPart], csMPI.Raw._COMM.SELF);

                    BlockMsrMatrix Block = new BlockMsrMatrix(localBlocking, localBlocking);// bi.Length, bi.Length, Bsz, Bsz);
                    Mop.WriteSubMatrixTo(Block, bi, default(int[]), bi, default(int[]));
                    if (l1.Count > 0) {
                        int offset = bi.Count;
                        int[] targRows = l1.Count.ForLoop(i => i + offset);
                        ExternalRowsTemp.WriteSubMatrixTo(Block, l1, targRows, bi, default(int[]));
                    }

                    blockSolvers[iPart] = new PARDISOSolver() {
                        CacheFactorization = true
                    };
                    //blockSolvers[iPart].DefineMatrix(Block);
                    //blockSolvers[iPart] = new FullDirectSolver();
                    //blockSolvers[iPart].DefineMatrix(Block);

                    blockSolvers[iPart] = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
                    blockSolvers[iPart].DefineMatrix(Block);
                }
            }

            // Record required indices
            // =======================
            {
                this.BlockIndices_Local = new int[NoParts][];
                int LocalI0 = MgMap.i0;
                int LocalLength = MgMap.LocalLength;

                for(int iBlock = 0; iBlock < NoParts; iBlock++) {
                    var _bi = _LocallyStoredBlockIndices[iBlock];
                    int L = _bi.Count;
                    int[] bil = new int[L];
                    this.BlockIndices_Local[iBlock] = bil;

                    for(int l = 0; l < L; l++) {
                        bil[l] = _bi[l] - LocalI0;
                        Debug.Assert(bil[l] >= 0);
                        Debug.Assert(bil[l] < MgMap.LocalLength);
                    }

                    var _biE = L2[iBlock];
                    if(_biE.Count > 0) {
                        L = _biE.Count;
                        int[] bie = new int[L];
                        this.BlockIndices_External[iBlock] = bie;

                        for(int l = 0; l < L; l++) {
                            bie[l] = _biE[l];
#if DEBUG
                            Debug.Assert(bie[l] >= MgMap.LocalLength);
                            int jMax = MgMap.AggGrid.iLogicalCells.NoOfCells - 1;
                            Debug.Assert(bie[l] < MgMap.LocalUniqueIndex(0, jMax, 0) + MgMap.GetLength(jMax));
#endif
                            bie[l] -= LocalLength;
                        }
                    }
                }
            }


            this.MtxFull = new ilPSP.LinSolvers.monkey.CPU.RefMatrix(Mop.ToMsrMatrix());

            if (CoarseSolver != null) {
                CoarseSolver.Init(op.CoarserLevel);
            }
        }

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

        ilPSP.LinSolvers.monkey.MatrixBase MtxFull;


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
        /// - content: external index, minus the local offset (i.e. can be used directly as index into <see cref="MPIexchange{T}.Vector_Ext"/>, <see cref="MPIexchangeInverse{T}{T}.Vector_Ext"/>).
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
                if(value < 0) {
                    throw new ArgumentException();
                }
                if(value > 2) {
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


        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            int NoParts = this.BlockIndices_Local.Length;

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

            double[] Res = new double[B.Count];
            var ResExchange = new MPIexchange<double[]>(this.m_MgOp.Mapping, Res);
            var XExchange = new MPIexchangeInverse<U>(this.m_MgOp.Mapping, X);

            for (int iIter = 0; iIter < m_MaxIterations; iIter++) {
                this.NoIter++;
                //double[] Res = B.ToArray();
                Res.SetV(B);
                this.MtxFull.SpMV(-1.0, X, 1.0, Res);

                if (IterationCallback != null)
                    IterationCallback(iIter, X.ToArray(), Res.CloneAs(), this.m_MgOp);

                if (CoarseSolver != null) {
                    var XC = X.ToArray().CloneAs();
                    double[] bc = new double[m_MgOp.CoarserLevel.Mapping.TotalLength];// = Res.CloneAs();
                    m_MgOp.CoarserLevel.Restrict(Res.CloneAs(), bc);
                    double[] xc = new double[bc.Length];                  
                    CoarseSolver.Solve(xc, bc);
                    m_MgOp.CoarserLevel.Prolongate(1,XC,1, xc);
                    X.AccV(1.0, XC);

                    if (CoarseSolverIsMultiplicative) {
                        Res.SetV(B);
                        this.MtxFull.SpMV(-1.0, X, 1.0, Res);
                    }
                }

                ResExchange.TransceiveStartImReturn();
                ResExchange.TransceiveFinish(0.0);

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
                    if(ciE != null)
                        bi.AccV(1.0, ResExchange.Vector_Ext, default(int[]), ciE, acc_index_shift: ci.Length);
                    
                    blockSolvers[iPart].Solve(xi, bi);

                    // accumulate block solution 'xi' to global solution 'X'
                    X.AccV(1.0, xi, ci, default(int[]));
                    if (ciE != null)
                        XExchange.Vector_Ext.AccV(1.0, xi, ciE, default(int[]), b_index_shift: ci.Length);
                }

                XExchange.TransceiveStartImReturn();
                XExchange.TransceiveFinish(1.0);
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

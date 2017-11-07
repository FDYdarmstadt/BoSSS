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
                    if (blokOp.CoarserLevel == null)
                        break;
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
            var M = op.OperatorMatrix;
            var MgMap = op.Mapping;
            this.m_MgOp = op;

            if (!M.RowPartitioning.EqualsPartition(MgMap.Partitioning))
                throw new ArgumentException("Row partitioning mismatch.");
            if (!M.ColPartition.EqualsPartition(MgMap.Partitioning))
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

                if (overlap < 0)
                    throw new ArgumentException();
                if (overlap > 0) {
                    foreach (List<int> bi in _Blocks) { // loop over blocks...
                        marker.SetAll(false); // marks all cells which are members of the block
                        foreach (int jcomp in bi)
                            marker[jcomp] = true;

                        // determine overlap regions
                        for (int k = 0; k < overlap; k++) {
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

            {
                int Jup = MgMap.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
                int Jgh = MgMap.AggGrid.iLogicalCells.NoOfExternalCells;


                var _BlockIndices = NoParts.ForLoop(iPart => new List<int>(BlockCells[iPart].Length * MgMap.MaximalLength));
                for (int iPart = 0; iPart < NoParts; iPart++) { // loop over parts...
                    int[] bc = BlockCells[iPart];
                    var bi = _BlockIndices[iPart];
                    int Jblock = bc.Length;


                    for (int jblk = 0; jblk < Jblock; jblk++) {
                        int j = bc[jblk];

                        if (j < Jup) {

                            int N = MgMap.GetLength(j);
                            int i0 = MgMap.LocalUniqueIndex(0, j, 0);

                            for (int n = 0; n < N; n++) {
                                bi.Add(i0 + n);
                            }
                        } else {
                            throw new NotImplementedException("todo: MPI parallelization;");
                        }
                    }

                }

                this.BlockIndices = _BlockIndices.Select(bi => bi.ToArray()).ToArray();
            }


            // test DOF blocks
            // ===============
            {
                int L = M.RowPartitioning.LocalLength;
                int i0 = M.RowPartitioning.i0;

                // ensure that each cell is used exactly once, among all blocks

                BitArray test = new BitArray(L);
                foreach (var bi in BlockIndices)
                    bi.ForEach(delegate (int i) {
                        if (i < L) {
                            Debug.Assert(test[i] == false || this.overlap > 0);
                            test[i] = true;
                        }
                    });
                for (int i = 0; i < L; i++)
                    Debug.Assert(test[i] == true);

            }
            

            // create solvers
            // ==============

            {
                blockSolvers = new ISparseSolver[NoParts];

                for (int iPart = 0; iPart < NoParts; iPart++) {
                    int[] bi = BlockIndices[iPart];

                    int Bsz;
                    if (MgMap.MinimalLength == MgMap.MaximalLength)
                        Bsz = MgMap.MaximalLength;
                    else
                        Bsz = 1;

                    if (M.RowPartitioning.MpiSize > 1) {
                        int i0Proc = M.RowPartitioning.i0;
                        bi = bi.CloneAs();
                        for (int i = 0; i < bi.Length; i++) {
                            bi[i] += i0Proc;
                        }
                    }


                    MsrMatrix Block = new MsrMatrix(bi.Length, bi.Length, Bsz, Bsz);
                    M.WriteSubMatrixTo(Block, bi, default(int[]), bi, default(int[]));

                    //blockSolvers[iPart] = new PARDISOSolver();
                    //blockSolvers[iPart].CacheFactorization = true;
                    //blockSolvers[iPart].DefineMatrix(Block);
                    //blockSolvers[iPart] = new FullDirectSolver();
                    //blockSolvers[iPart].DefineMatrix(Block);

                    blockSolvers[iPart] = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
                    blockSolvers[iPart].DefineMatrix(Block);
                }
            }

            this.MtxFull = new ilPSP.LinSolvers.monkey.CPU.RefMatrix(M.ToMsrMatrix());

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
        /// 1st index: block
        /// 2nd index: enumeration; list of composite cells which make up the respective block
        /// </summary>
        int[][] BlockCells;


        int[][] BlockIndices;

        public int overlap = 2;


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
            int NoParts = this.BlockIndices.Length;

            for (int iIter = 0; iIter < m_MaxIterations; iIter++) {
                this.NoIter++;
                double[] Res = B.ToArray();
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

                for (int iPart = 0; iPart < NoParts; iPart++) {
                    int[] ci = BlockIndices[iPart];
                    int L = ci.Length;

                    double[] bi = new double[L];
                    double[] xi = new double[L];
                    bi.AccV(1.0, Res, default(int[]), ci);

                    blockSolvers[iPart].Solve(xi, bi);


                    X.AccV(1.0, xi, ci, default(int[]));
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

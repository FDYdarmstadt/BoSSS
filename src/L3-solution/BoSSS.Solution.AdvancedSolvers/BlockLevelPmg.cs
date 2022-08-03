using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.AdvancedSolvers {
    /*
    /// <summary>
    /// The main purpose of this factory is to wrap up input parameters for the block solvers and to provide buffer arrays.
    /// All block solver share these buffers, which reduces the pressure on the garbage collector.
    /// </summary>
    class BlockLevelPmgFactory {

        /// <summary>
        /// Factory for block solvers, which based on level pmg alias p-two-grid.
        /// Sets parameters, which all block solver share.
        /// Allocates buffer arrays, for index translation between block - sub-block and block - local Matrix.
        /// </summary>
        /// <param name="mgo">the basis mapping and matrix is derived from this</param>
        /// <param name="ExtMatrix">if blocks exceed local range, which is stored on this processor, it is reasonable to exchange external rows first and hand them over</param>
        public BlockLevelPmgFactory(MultigridOperator mgo, BlockMsrMatrix ExtMatrix) {
            m_op = mgo;
            m_ExtMatrix = ExtMatrix;

            // set dummy arrays, to reduce pressure on garbage collector
            bool ExtCellsIncluded = ExtMatrix != null && ExtMatrix.NoOfRows > 0;
            AllocateBuffer(ExtCellsIncluded, mgo);
            Debug.Assert(m_XBuffer != null && m_XBuffer.Length > 0);
            Debug.Assert(m_ResBuffer != null && m_ResBuffer.Length > 0);
        }

        public bool FullSolveOfCutcells = true;
        public bool EqualOrder = false;
        public int pLow = 1;
        private MultigridOperator m_op;
        private BlockMsrMatrix m_ExtMatrix;
        private double[] m_XBuffer;
        private double[] m_ResBuffer;

        /// <summary>
        /// Instantiates and initializes the <see cref="LevelPmg"/> block solvers.
        /// Use this one, if you have to calculate block mask and matrix anyway, to avoid redundancy.
        /// </summary>
        /// <param name="BlockCellIdc">input cell indices</param>
        /// <param name="fullBlock">the matrix block corresponding to input</param>
        /// <param name="fullMask">the corresponding </param>
        /// <returns></returns>
        public BlockLevelPmg CreateAndInit(List<int> BlockCellIdc, out BlockMsrMatrix fullBlock, out BlockMask fullMask) {
            var solver = new BlockLevelPmg() {
                m_FullSolveOfCutcells = FullSolveOfCutcells,
                m_pLow = pLow
            }; 
            Debug.Assert(m_XBuffer.Length == m_ExtMatrix.NoOfRows + m_op.Mapping.LocalLength);
            Debug.Assert(m_ResBuffer.Length == m_ExtMatrix.NoOfRows + m_op.Mapping.LocalLength);
            solver.Xdummy = m_XBuffer;
            solver.Resdummy = m_ResBuffer;

            var fullSel = new SubBlockSelector(m_op.Mapping);
            fullSel.CellSelector(BlockCellIdc.ToList(), false);
            var ExtRows = BlockMask.GetAllExternalRows(m_op.Mapping, m_op.OperatorMatrix);
            fullMask = new BlockMask(fullSel, ExtRows);
            fullBlock = fullMask.GetSubBlockMatrix_MpiSelf(m_op.OperatorMatrix);

            solver.Init(m_op,BlockCellIdc, m_ExtMatrix, fullBlock, fullMask);
            return solver;
        }

        /// <summary>
        /// these buffers are used for index translation: block - sub-block and block - local Matrix.
        /// </summary>
        /// <param name="EnhancedBuffer"></param>
        /// <param name="mgo"></param>
        private void AllocateBuffer(bool EnhancedBuffer, MultigridOperator mgo) {
            if (m_XBuffer == null || m_ResBuffer == null) {
                int ExtLen = EnhancedBuffer ? BlockMask.GetLocalAndExternalDOF(mgo.Mapping) : mgo.Mapping.LocalLength;
                m_XBuffer = new double[ExtLen];
                m_ResBuffer = new double[ExtLen];
            } else {
                m_XBuffer.ClearEntries();
                m_ResBuffer.ClearEntries();
            }
        }
    }

    */

    /*
    class BlockLevelPmg : IDisposable {

        public bool m_FullSolveOfCutcells = true;
        //public bool m_EqualOrder = false;
        public int m_pLow = 1;
        MultigridOperator m_op = null;

        // dummy arrays used for index translation
        // working around the lack of sub-mapping
        public double[] Resdummy = null;
        public double[] Xdummy = null;

        /// <summary>
        /// masks for the Schwarz blocks, high order modes, only initialized if PMG is used, <see cref="Schwarz.Config.UsePMGinBlocks"/>
        /// - index: Schwarz block
        /// </summary>
        BlockMask BMhiBlocks;

        /// <summary>
        /// masks for the Schwarz blocks, low order modes, only initialized if PMG is used, <see cref="Schwarz.Config.UsePMGinBlocks"/>
        /// - index: Schwarz block
        /// </summary>
        BlockMask BMloBlock;

        /// <summary>
        /// masks for the Schwarz blocks
        /// - index: Schwarz block
        /// </summary>
        BlockMask m_BMfullBlock;

        /// <summary>
        /// Cell-local solvers for the high-order modes 
        /// - 1st index: Schwarz block
        /// - 2nd index: sub-block within the respective Schwarz block (there is one sub-block for each cell)
        /// - content: matrix inverse of the diagonal block for high-order DG modes
        /// </summary>
        MultidimensionalArray[] HiModeBlocks;

        /// <summary>
        /// Linear solver for each block
        /// </summary>
        ISparseSolver loSolver;

        /// <summary>
        /// contains the Block MSR matrix of the full block.
        /// This should be a matrix small than the local size of the operator matrix
        /// </summary>
        BlockMsrMatrix m_fullBlock;

        /// <summary>
        /// only used for memory calculation
        /// </summary>
        BlockMsrMatrix loModeBlock;

        ///// <summary>
        ///// experimental. LU Pivoting. Can be dismissed ...
        ///// </summary>
        //int[][] HighOrderBlocks_LUpivots;

        private bool AnyHighOrderTerms {
            get {
                Debug.Assert(m_op != null, "there is no matrix given yet!");
                return m_op.Mapping.DgDegree.Any(p => p > m_pLow);
            }
        }

        
        /// <summary>
        /// For performance sake, just build the objects, which are passed to this function, once.
        /// </summary>
        /// <param name="op"></param>
        /// <param name="BlockCellIdc"></param>
        /// <param name="ExtMatrix"></param>
        /// <param name="fullBlock">used for residual re-evaluation</param>
        /// <param name="_BMfullBlock"></param>
        public void Init(MultigridOperator op, List<int> BlockCellIdc, BlockMsrMatrix ExtMatrix, BlockMsrMatrix fullBlock, BlockMask _BMfullBlock) {
            //generate selector instructions
            m_op = op;
            m_BMfullBlock = _BMfullBlock;
            m_fullBlock = fullBlock;
            int D = m_op.GridData.SpatialDimension;

            var lowSel = new SubBlockSelector(op.Mapping);
            lowSel.CellSelector(BlockCellIdc, false);

            int[] lowDegs = op.GetBestFitLowOrder(m_pLow);
            bool LowSelector(int iCell, int iVar, int iSpec, int pDeg) {
                return pDeg <= lowDegs[iVar];
            }


            //lowSel.SetModeSelector((int iCell, int iVar, int iSpec, int pDeg) => pDeg <= (iVar != D && !m_EqualOrder ? m_pLow : m_pLow - 1));
            lowSel.SetModeSelector(LowSelector);
            if (m_FullSolveOfCutcells)
                ModifyLowSelector(lowSel, op);

            var HiSel = new SubBlockSelector(op.Mapping);
            HiSel.CellSelector(BlockCellIdc, false);
            HiSel.SetModeSelector((int iCell, int iVar, int iSpec, int pDeg) => !LowSelector(iCell, iVar, iSpec, pDeg));
            if (m_FullSolveOfCutcells)
                ModifyHighSelector(HiSel, op);

            //generate Blockmasking
            var lowMask = new BlockMask(lowSel, ExtMatrix);
            Debug.Assert(lowMask != null);
            BMloBlock = lowMask;
            var HiMask = new BlockMask(HiSel, ExtMatrix);
            Debug.Assert(HiMask.NoOfMaskedCells == lowMask.NoOfMaskedCells);
            Debug.Assert(HiMask.NoOfMaskedCells == fullBlock._RowPartitioning.LocalNoOfBlocks);

            //Console.WriteLine("Testcode in Schwarz.");
            //Debug.Assert(HiMask.GetNoOfMaskedCells == BlockCellIdc.Count() || m_FullSolveOfCutcells); // Probably not fulfilled for IBM, there maybe emtpy cells, due to no species
            //Debug.Assert(lowMask.GetNoOfMaskedCells == BlockCellIdc.Count()); // same here
            BMhiBlocks = HiMask;

            //get subblocks from masking
            MultidimensionalArray[] hiBlocks = HiMask.GetDiagonalBlocks(op.OperatorMatrix, false, false); //gets diagonal-blocks only        
            var loBlock = lowMask.GetSubBlockMatrix_MpiSelf(op.OperatorMatrix);

            //get inverse of high-order blocks
            if (hiBlocks != null) {
                foreach (var block in hiBlocks) {
                    try {
                        block.InvertInPlace();
                    } catch (Exception e) {
                        Console.WriteLine(e);
                    }
                }
            }

            // LU Pivoting showed no advantage vs Lapack invert ... 
            //HighOrderBlocks_LUpivots = new int[hiBlocks.Length][];
            //if (hiBlocks != null) {
            //    for (int jBlock=0;jBlock<hiBlocks.Length;jBlock++) {
            //        int len = hiBlocks[jBlock].NoOfRows;
            //        HighOrderBlocks_LUpivots[jBlock] = new int[len];
            //        hiBlocks[jBlock].FactorizeLU(HighOrderBlocks_LUpivots[jBlock]);
            //    }
            //}
            
            HiModeBlocks = hiBlocks;


#if TEST
                        int cnt = 0;
                        foreach (var block in PmgBlock_HiModeSolvers[iPart]) {
                            block.SaveToTextFile(String.Format("{0}_block{1}_part{2}", op.Mapping.MpiRank, cnt, iPart));
                            cnt++;
                        }

                        loBlock.SaveToTextFileSparseDebug(String.Format("{0}_block_part{1}", op.Mapping.MpiRank, iPart));
#endif
            loModeBlock = loBlock;
            loSolver = new PARDISOSolver() {
                CacheFactorization = true,
                UseDoublePrecision = true,
                Parallelism = Parallelism.SEQ
            };
            loSolver.DefineMatrix(loBlock);

        }
        

        private void ModifyLowSelector(SubBlockSelector sbs, MultigridOperator op) {
            AssignXdgBlocksModification(sbs, op, true);
        }

        private void ModifyHighSelector(SubBlockSelector sbs, MultigridOperator op) {
            AssignXdgBlocksModification(sbs, op, false);
        }

        /// <summary>
        /// Executes the assignment of cut-cell blocks to the low solver.
        /// </summary>
        /// <param name="sbs"></param>
        /// <param name="op"></param>
        /// <param name="IsLowSelector"></param>
        private void AssignXdgBlocksModification(SubBlockSelector sbs, MultigridOperator op, bool IsLowSelector) {
            var Filter = sbs.ModeFilter;
            Func<int, int, int, int, bool> Modification = delegate (int iCell, int iVar, int iSpec, int pDeg) {
                int NoOfSpec = op.Mapping.AggBasis[0].GetNoOfSpecies(iCell);
                if (NoOfSpec >= 2)
                    return IsLowSelector;
                else
                    return Filter(iCell, iVar, iSpec, pDeg);
            };
            sbs.SetModeSelector(Modification);
        }



        public void Solve<U, V>(U X, V B)
           where U : IList<double>
           where V : IList<double> // 
       {
            Resdummy.ClearEntries();
            Xdummy.ClearEntries();

            var bi = new double[B.Count()];
            bi.SetV(B);
            m_BMfullBlock.AccSubVec(bi, Resdummy);

            // solve low-order system
            SolveLoSystem();
            // re-evaluate residual
            ReEvaluateRes(bi);
            // hi order smoother
            SolveHiSystem();

            //// this gives better results, but is more costly
            //// experiments have shown,
            //// it is probably not worth the additional effort
            //ReEvaluateRes(iPart, bi);
            //SolveLoSystem(iPart);

            var xi = m_BMfullBlock.GetSubVec(Xdummy);
            Debug.Assert(xi.Length == X.Count());
            X.SetV(xi);
            B.SetV(bi);
        }

        private void ReEvaluateRes(double[] res) {
            using (new FuncTrace()) {
                var ri = res.CloneAs();
                Resdummy.ClearEntries();
                var xi = m_BMfullBlock.GetSubVec(Xdummy);
                this.m_fullBlock.SpMV(-1.0, xi, 1.0, ri);
                m_BMfullBlock.AccSubVec(ri, Resdummy);
            }
        }

        private void SolveHiSystem() {
            using (new FuncTrace()) {
                int NoCells = HiModeBlocks.Length;

                double[] xiHi = null;
                double[] biHi = null;

                for (int j = 0; j < NoCells; j++) {
                    var HiModeSolver = HiModeBlocks[j];
                    int Np = HiModeSolver.NoOfRows;
                    if (xiHi == null || xiHi.Length != Np)
                        xiHi = new double[Np];
                    if (biHi == null || biHi.Length != Np)
                        biHi = new double[Np];
                    biHi = BMhiBlocks.GetSubVecOfCell(Resdummy, j);
                    //HiModeSolver.BacksubsLU(HighOrderBlocks_LUpivots[j], xiHi, biHi);
                    HiModeSolver.GEMV(1.0, biHi, 0.0, xiHi);
                    BMhiBlocks.AccSubVecOfCell(xiHi, j, Xdummy);
                }
            }
        }

        private void SolveLoSystem() {
            using (new FuncTrace()) {
                var rLo = BMloBlock.GetSubVec(Resdummy);
                double[] xiLo = new double[rLo.Length];
                try {
                    loSolver.Solve(xiLo, rLo);
                } catch (ArithmeticException ae) {
                    Console.Error.WriteLine(ae.Message);
                    throw ae;
                }
                BMloBlock.AccSubVec(xiLo, Xdummy);
            }
        }

        /// <summary>
        /// Gives memory consumption of used matrices
        /// </summary>
        public long UsedMem {
            get {
                long s = 0;

                if (loModeBlock != null) {
                    s += loModeBlock.UsedMemory;
                }

                if (HiModeBlocks != null) {
                    s += HiModeBlocks.Sum(mda => (long)mda.Length * sizeof(double));
                }

                return s;
            }
        }

        public void Dispose() {
            throw new NotImplementedException();
        }
    }

    */
    
}

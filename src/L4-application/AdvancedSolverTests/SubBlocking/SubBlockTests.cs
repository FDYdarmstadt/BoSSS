using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Foundation;

namespace AdvancedSolverTests
{
    internal class TestMasking : BlockMask
    {
        public TestMasking(SubBlockSelector sbs, BlockMsrMatrix ExtRows = null) : base(sbs, ExtRows) {
            Global_IList_LocalCells = base.GlobalIList_Internal;
            Global_IList_ExternalCells = base.GlobalIList_External;
        }

        public List<int> Global_IList_LocalCells;
        public List<int> Global_IList_ExternalCells;

        public int[] GetGlobalIdx(int iCell) {
            return base.GetGlobalIdxOfCell(iCell);
        }
        public int[] GetLocalIdx(int iCell) {
            return base.GetLocalIdxOfCell(iCell);
        }
    }
    static public class SubBlockTests
    {

        public static void TestInit(params int[] Testparameters) {
            int L=Testparameters.Length;
            unsafe {
                int[] Params = new int[L*2], ParamsGlob = new int[L*2];
                fixed (int* pParams = Params, pParamsGlob = ParamsGlob) {
                    for(int i = 0; i < L; i++) {
                        pParams[i] = (int)Testparameters[i];
                        pParams[i + L] = -pParams[i];
                    }
                    

                    csMPI.Raw.Allreduce((IntPtr)pParams, (IntPtr)pParamsGlob, L*2, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                }

                int[] ParamsMin = ParamsGlob.GetSubVector(0, L);
                int[] ParamsMax = ParamsGlob.GetSubVector(L, L);
                for (int i = 0; i < L; i++) {
                    if (Params[i] != ParamsMin[i])
                        throw new ApplicationException();
                    if (Params[i] != -ParamsMax[i])
                        throw new ApplicationException();
                }
            }
        }

        public static void GetExternalRowsTest() {

        }


        public static void 

        public static void IndexingTest(
            [Values(XDGusage.none, XDGusage.mixed1, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder) {

            TestInit((int)UseXdg, DGOrder);

            Console.WriteLine("TestIndexing({0},{1})", UseXdg, DGOrder);

            using (var solver = new SubBlockTestSolver() { m_UseXdg = UseXdg, m_DGorder = DGOrder }) {
                solver.Init(null);
                solver.RunSolverMode();

                //Arrange --- Get global index by mapping
                var map = solver.GetMapping;
                int[] fields = map.NoOfVariables.ForLoop(i=>i);
                int[] GlobalIdxMap_loc = map.GetSubvectorIndices(fields);
                int[] GlobalIdxMap_ext = map.GetSubvectorIndices_Ext(fields); // maybe not sorted, so do it ...
                Array.Sort(GlobalIdxMap_ext);

                //Arrange --- Prepare stuff for mask
                var selector = new SubBlockSelector(map);
                var dummy = new BlockMsrMatrix(map); // we are not interested in any operational stuff
                var stw = new Stopwatch();
                stw.Reset();

                //Act --- do the masking to get index lists
                stw.Start();
                var mask = new TestMasking(selector, dummy);
                stw.Stop();
                int[] GlobalIdxMask_loc = mask.Global_IList_LocalCells.ToArray();
                int[] GlobalIdxMask_ext = mask.Global_IList_ExternalCells.ToArray();

                //Act --- Get global

                //Assert --- Idx lists are of same length
                Assert.IsTrue(GlobalIdxMap_loc.Length == GlobalIdxMask_loc.Length);
                Assert.IsTrue(GlobalIdxMap_ext.Length == GlobalIdxMask_ext.Length);

                //Assert --- local idx list contains local idx
                Assert.IsTrue(map.IsInLocalRange(GlobalIdxMask_loc[0]));
                Assert.IsTrue(map.IsInLocalRange(GlobalIdxMask_loc.Last()));

                //Assert --- external idx list contains external idx
                Assert.IsTrue(!map.IsInLocalRange(GlobalIdxMask_ext[0]));
                Assert.IsTrue(!map.IsInLocalRange(GlobalIdxMask_ext.Last()));

                //Act --- Compare Idx
                bool TestLoc = true;
                for (int iLoc = 0; iLoc < GlobalIdxMask_loc.Length; iLoc++) {
                    GlobalIdxMap_loc[iLoc]==
                }
                


                mask.Global_IList_LocalCells
                mask.Global_IList_ExternalCells
            }
        }

        /// <summary>
        /// Test for extracting masked diagonal cell blocks and respective vector matrix multiplication
        /// </summary>
        [Test]
        public static void ExtractDiagonalBlocks(
            [Values(XDGusage.none, XDGusage.mixed1, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal,MatrixShape.diagonal_var,MatrixShape.diagnoal_spec,MatrixShape.diagnoal_var_spec)] MatrixShape MShape
            ) {

            TestInit((int)UseXdg, DGOrder,(int)MShape);

            Console.WriteLine("ExtractDiagonalBlocks({0},{1},{2})", UseXdg, DGOrder, MShape);

                using (var solver = new SubBlockTestSolver() { m_UseXdg = UseXdg, m_DGorder = DGOrder, m_Mshape = MShape }) {

                    //Arrange --- test matrix and subblock masking
                    solver.Init(null);
                    solver.RunSolverMode();

                    Stopwatch stw = new Stopwatch();
                    stw.Reset();
                    

                    BlockMsrMatrix M = solver.OperatorMatrix;
                    
                    SubBlockSelector SBS = new SubBlockSelector(solver.GetMapping);
                    
                //masking operates on row and col indices,
                //but we want to ignore block coupling ...
                bool[] coupling = new bool[3];
                switch (MShape) {
                    case MatrixShape.diagonal:
                        coupling = new bool[]{ true, false , false};
                        break;
                    case MatrixShape.diagonal_var:
                        coupling = new bool[] { true, true, false };
                        break;
                    case MatrixShape.diagnoal_spec:
                        coupling = new bool[] { true, false, true };
                        break;
                    case MatrixShape.diagnoal_var_spec:
                        coupling = new bool[] { true, true, true };
                        break;
                    default:
                        throw new NotSupportedException(String.Format("{0} is not supported by this test",MShape));
                }

                TestMasking mask = new TestMasking(SBS, null);

                //Act --- diagonal subblock extraction
                stw.Start();
                var blocks = mask.GetSubBlocks(M, coupling[0], coupling[1], coupling[2]);
                stw.Stop();

                //Assert --- number of extracted blocks equal to subblocks total 
                Assert.IsTrue(blocks.Length == solver.GetMapping.LocalNoOfBlocks);

                //Arrange --- build temporal matrix to store subblocks
                mask.Global_IList_LocalCells;

                for (int i = 0; i < blocks.Length; i++) {
                    //Arrange --- Get diagonal subblock of ith cell
                    stw.Start();
                    int[] MblockLIdx = mask.GetGlobalIdx(i);
                    stw.Stop();
                    int MblockL = MblockLIdx.Length;
                    //var Mblock = MultidimensionalArray.Create(MblockL, MblockL);
                    //M.ReadBlock(MblockLIdx[0], MblockLIdx[0], Mblock);


                    //Act --- diff of entries of ith diag. block in matrix and resp. subblock
                    Mblock.Acc(-1.0, blocks[i]);

                    //Act --- multiply subvectors with subblocks 
                    mask.GetVectorCellwise(i);

                    //Assert --- Are subblock and ith diagonal block entries equal?
                    Assert.IsTrue(Mblock.InfNorm() == 0.0, String.Format("infNorm of block {0} neq 0!", i));
                }
                
                



            }

        }

        public static void ExtractSubMatrix(
            [Values(XDGusage.none, XDGusage.mixed1, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal, MatrixShape.diagonal_var, MatrixShape.diagnoal_spec, MatrixShape.diagnoal_var_spec)] MatrixShape MShape
            ) {

            TestInit((int)UseXdg, DGOrder, (int)MShape);

            Console.WriteLine("ExtractDiagonalBlocks({0},{1},{2})", UseXdg, DGOrder, MShape);

            using (var solver = new SubBlockTestSolver() { m_UseXdg = UseXdg, m_DGorder = DGOrder, m_Mshape = MShape }) {

                //Arrange --- test matrix and subblock masking
                solver.Init(null);
                solver.RunSolverMode();

                Stopwatch stw = new Stopwatch();
                stw.Reset();


                BlockMsrMatrix M = solver.OperatorMatrix;

                SubBlockSelector SBS = new SubBlockSelector(solver.GetMapping);


            }
        }
}

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

namespace AdvancedSolverTests
{
    internal class TestMasking : BlockMask
    {
        public TestMasking(SubBlockSelector sbs, BlockMsrMatrix ExtRows = null) : base(sbs, ExtRows) {
            Global_IList_LocalCells = this.GlobalIList_Internal;
            Global_IList_ExternalCells = this.GlobalIList_External;
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

        /// <summary>
        /// Test for submatrix extraction.
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

                    solver.Init(null);
                    solver.RunSolverMode();

                    Stopwatch stw = new Stopwatch();
                    stw.Reset();
                    stw.Start();

                    BlockMsrMatrix M = solver.OperatorMatrix;
                    
                    SubBlockSelector SBS = new SubBlockSelector(solver.GetMapping);
                    TestMasking mask = new TestMasking(SBS,null); // testing mpi-local extraction

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
                    var blocks = mask.GetSubBlocks(M, coupling[0], coupling[1], coupling[2]); // get diagonal cell blocks only
                    Assert.IsTrue(blocks.Length == solver.GetMapping.LocalNoOfBlocks);

                    for (int i = 0; i < blocks.Length; i++) {
                        int MblockL = mask.GetGlobalIdx(i).Length;
                        var Mblock = MultidimensionalArray.Create(MblockL, MblockL);
                        M.ReadBlock(mask.GetGlobalIdx(i)[0], mask.GetGlobalIdx(i)[0], Mblock);
                        Mblock.Acc(-1.0, blocks[i]);
                        Assert.IsTrue(Mblock.InfNorm() == 0.0, String.Format("infNorm of block {0} neq 0!", i));
                    }

                    stw.Stop();

                if (MShape == MatrixShape.diagnoal_spec) {

                }
                }

        }
    }
}

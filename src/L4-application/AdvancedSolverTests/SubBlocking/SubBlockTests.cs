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

        [Test]
        public static void GetExternalRowsTest(
            ) {
            //Matlabaufruf --> gesamte Matrix nach Matlab schreiben
            //Teilmatritzen gemäß Globalid extrahieren
            //Mit ExternalRows vergleichen
            //Die große Frage: funktioniert der batchmode connector parallel? Beim rausschreiben beachten
        }

        [Test]
        public static void LocalIndexTest(
            [Values(XDGusage.none, XDGusage.mixed1, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder) {

            TestInit((int)UseXdg, DGOrder);

            Console.WriteLine("TestIndexing({0},{1})", UseXdg, DGOrder);

            using (var solver = new SubBlockTestSolver() { m_UseXdg = UseXdg, m_DGorder = DGOrder }) {
                solver.Init(null);
                solver.RunSolverMode();

                //Arrange --- Get global index by mapping
                var map = solver.GetMapping;
                int[] fields = map.NoOfVariables.ForLoop(i => i);
                int[] GlobalIdxMap_loc = map.GetSubvectorIndices(fields);
               
                //Arrange --- Prepare stuff for mask
                var selector = new SubBlockSelector(map);
                var stw = new Stopwatch();
                stw.Reset();

                //Act --- do the masking to get index lists
                stw.Start();
                var mask = new TestMasking(selector);
                stw.Stop();
                int[] GlobalIdxMask_loc = mask.Global_IList_LocalCells.ToArray();

                //Assert --- Idx lists are of same length
                Assert.IsTrue(GlobalIdxMap_loc.Length == GlobalIdxMask_loc.Length);

                //Assert --- Compare map and mask indices
                for (int iLoc = 0; iLoc < GlobalIdxMask_loc.Length; iLoc++) {
                    Assert.True(GlobalIdxMap_loc[iLoc] == GlobalIdxMask_loc[iLoc]);
                }
            }
        }

        [Test]
        public static void ExternalIndexTest(
            [Values(XDGusage.none, XDGusage.mixed1, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder) {

            TestInit((int) UseXdg, DGOrder);

            Console.WriteLine("TestIndexing({0},{1})", UseXdg, DGOrder);

            using (var solver = new SubBlockTestSolver() { m_UseXdg = UseXdg, m_DGorder = DGOrder }) {
                solver.Init(null);
                solver.RunSolverMode();

                //Arrange --- Get global index by mapping
                var map = solver.GetMapping;
                int[] fields = map.NoOfVariables.ForLoop(i => i);
                int[] GlobalIdxMap_ext = map.GetSubvectorIndices_Ext(fields); // maybe not sorted, so do it ...
                Array.Sort(GlobalIdxMap_ext);

                //Arrange --- Prepare stuff for mask
                var selector = new SubBlockSelector(map);
                var dummy = new BlockMsrMatrix(map); // we do not need the external rows yet, because we are not testing any operations
                var stw = new Stopwatch();
                stw.Reset();

                //Act --- do the masking to get index lists
                stw.Start();
                var mask = new TestMasking(selector, dummy);
                stw.Stop();
                int[] GlobalIdxMask_ext = mask.Global_IList_ExternalCells.ToArray();

                //Assert --- Idx lists are of same length
                Assert.IsTrue(GlobalIdxMap_ext.Length == GlobalIdxMask_ext.Length);

                //Assert --- Compare map and mask indices
                for (int iLoc = 0; iLoc< GlobalIdxMask_ext.Length; iLoc++) {
                    Assert.IsTrue(GlobalIdxMask_ext[iLoc] == GlobalIdxMap_ext[iLoc]);
                }
            }
        }

        /// <summary>
        /// Test for extracting masked diagonal cell blocks and respective vector matrix multiplication
        /// </summary>
        [Test]
        public static void ExtractDiagonalCellBlocks(
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



                //BlockMsrMatrix all1;
                //all1.SetAll(1);
                //Generate broken diagonal matrix, die zur Maske passt: M
                //M+all1=M_prep
                //Wende Extraction auf M_prep an, Man sollte nun M bekommen
                //Test: M_prep-extract(M_prep)=all1
                //Test-crit: Result.SumEntries=DOF^2 oder Result.Max()==Result.Min()==1
                //oder (besser)
                //Test: M-extract(M_prep)=zeros
                //Test-crit: Result.InfNorm()==0

                //Der Test kann für ExtractSubMatrix mit ignore coupling wiederholt werden
                //eventuell: Testmatrix finden mit brauchbaren Nebendiagonalen für einen Fall

                //Was wird getestet: funktioniert ignorecoupling richtig?

            }

        }

        public static void AccSubVectorCellwise() {
            //matrix Erzeugung wie in ExtractDiagonalCellBlocks...
            //Auf der HierarchieEbene, auf der Kopplung ausgesetzt wird kann Auswahl vorgenommen werden
            //bei var: 0 / 1, bei DG: <=1 / >1, bei spec: A / B, bei Cells: odd / even
            //accumulierte Teilergebnisse sind dann == fullM*fullX 
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

                //Hier interessiert uns kein ignorecoupling
                //Die Frage ist: werden blöcke abseits der Diagonalen richtig behandelt?
                //Bauen 2er Testmatritzen durch solver:
                //M1: enthält nur Elemente mit Var1 / SpecA / DG<=1 (Einstellen über Manipulation der Flüsse/Sourceterme)
                //M2: keine Einschränkung aber gleiche Flüsse und Sourceterme wie M1
                //Maskierung gemäß M1: nur Elemente mit Var1 / SpecA / DG<=1
                //Test: extract(M2)-M1=zeros
                //Test-crit: Result.InfNorm()==0

                //Eventuell auch möglich für Zellen
            }
        }

        public static void AccSubVector() {
            //matrix Erzeugung wie in ExtractSubMatrix ...
            //Auf der HierarchieEbene, auf der Kopplung ausgesetzt wird kann Auswahl vorgenommen werden
            //bei var: 0 / 1, bei DG: <=1 / >1, bei spec: A / B, bei Cells: odd / even
            //accumulierte Teilergebnisse sind dann == fullM*fullX 
        }

        public static void ExtractSubMatrix_Ext() {
            //Kombination aus GetExternalRows und SubMatrix extraction
            //Es werden nur external Row Einträge berücksichtigt
            //Berechnung eventuell in Matlab durchführen
        }

        public static void AccSubVec_Ext() {
            //Kombination aus GetExternalRows und AccSubVec_Ext
            //Es werden nur external Row Einträge berücksichtigt
            //Berechnung eventuell in Matlab durchführen
        }
    }

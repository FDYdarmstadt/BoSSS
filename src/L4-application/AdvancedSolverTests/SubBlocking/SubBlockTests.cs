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
using System.Threading;
using System.Globalization;

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

    [TestFixture]
    static public class SubBlockTests
    {

        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);
            Thread.CurrentThread.CurrentCulture = CultureInfo.CurrentCulture;
        }

        /// <summary>
        /// MPI shutdown.
        /// </summary>
        [TestFixtureTearDown]
        public static void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        public static void TestInit(params int[] Testparameters) {
            int L = Testparameters.Length;
            unsafe {
                int[] Params = new int[L * 2], ParamsGlob = new int[L * 2];
                fixed (int* pParams = Params, pParamsGlob = ParamsGlob) {
                    for (int i = 0; i < L; i++) {
                        pParams[i] = (int)Testparameters[i];
                        pParams[i + L] = -pParams[i];
                    }


                    csMPI.Raw.Allreduce((IntPtr)pParams, (IntPtr)pParamsGlob, L * 2, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
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

        public static void WriteOutTestMatrices() {
            MultigridOperator mgo;
            mgo=CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal);
            mgo.OperatorMatrix.SaveToTextFileSparseDebug("M");
            mgo = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_var);
            mgo.OperatorMatrix.SaveToTextFileSparseDebug("M_var");
            mgo = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_spec);
            mgo.OperatorMatrix.SaveToTextFileSparseDebug("M_spec");
            mgo = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_var_spec);
            mgo.OperatorMatrix.SaveToTextFileSparseDebug("M_var_spec");
        }

        public static MultigridOperator CreateTestMGOperator(XDGusage UseXdg = XDGusage.none, int DGOrder = 2, MatrixShape MShape = MatrixShape.standard, int Resolution = 2) {
            return CreateTestMGOperator(out double[] Vec, UseXdg, DGOrder, MShape, Resolution);
        }

        public static MultigridOperator CreateTestMGOperator(out double[] Vec, XDGusage UseXdg = XDGusage.none, int DGOrder = 2, MatrixShape MShape = MatrixShape.standard, int Resolution = 2) {
            MultigridOperator retMGOp;
            using (var solver = new SubBlockTestSolver2Var() { m_UseXdg = UseXdg, m_DGorder = DGOrder, m_Mshape = MShape, m_Res = Resolution }) {
                solver.Init(null);
                solver.RunSolverMode();
                retMGOp = solver.MGOp;
                Vec=solver.someVec;
            }
            return retMGOp;
        }

        public static void SetAll(this BlockMsrMatrix A, double val) {
            var rowmap = A._RowPartitioning;
            var colmap = A._ColPartitioning;
            int RowBlocks = rowmap.LocalNoOfBlocks;
            int ColBlocks = colmap.LocalNoOfBlocks;
            
            Partitioning rowpart = new Partitioning(RowBlocks);

            for (int iBlock = rowpart.i0; iBlock < rowpart.iE; iBlock++) {
                for (int jBlock = rowpart.i0; jBlock < rowpart.iE; jBlock++) {
                    int i0 = rowmap.GetBlockI0(iBlock);
                    int j0 = colmap.GetBlockI0(jBlock);
                    int iL = rowmap.GetBlockLen(iBlock);
                    int jL = colmap.GetBlockLen(jBlock);
                    var subM = MultidimensionalArray.Create(iL, jL);
                    subM.SetAll(val);
                    A.AccBlock(i0, j0, 1.0, subM);
                }
            }
            double min, max;
            int minc, minr, maxc, maxr;
            A.GetMinimumAndMaximum_MPILocal(out min, out minr, out minc, out max, out maxr, out maxc);
            Debug.Assert(min == max);
            Debug.Assert(min == val);
        }

        public static BlockMsrMatrix CreateShapeOfOnes(BlockMsrMatrix A) {
            var rowmap = A._RowPartitioning;
            var colmap = A._ColPartitioning;
            int RowBlocks = rowmap.LocalNoOfBlocks;
            int ColBlocks = colmap.LocalNoOfBlocks;

            BlockMsrMatrix B = new BlockMsrMatrix(rowmap, colmap);

            Partitioning rowpart = new Partitioning(RowBlocks);
            for (int iBlock = rowpart.i0; iBlock < rowpart.iE; iBlock++) {
                for (int jBlock = rowpart.i0; jBlock < rowpart.iE; jBlock++) {
                    int i0 = rowmap.GetBlockI0(iBlock);
                    int j0 = colmap.GetBlockI0(jBlock);
                    int iL = rowmap.GetBlockLen(iBlock);
                    int jL = colmap.GetBlockLen(jBlock);
                    var subM = MultidimensionalArray.Create(iL, jL);
                    A.ReadBlock(i0, j0, subM);
                    subM.ApplyAll(i => i != 0.0 ? 1 : 0);
                    B.AccBlock(i0, j0, 1.0, subM);
                }
            }
            double min, max;
            int minc, minr, maxc, maxr;
            B.GetMinimumAndMaximum_MPILocal(out min, out minr, out minc, out max, out maxr, out maxc);
            Debug.Assert(min == 0);
            Debug.Assert(max == 1);
            return B;
        }

        static public bool IsLocallyEqual(this IPartitioning t, IPartitioning o) {
            if (t == null && o == null)
                return true;
            if (o == null)
                return false;
            if (t == null)
                return false;
            if (object.ReferenceEquals(t, o))
                return true;
            if (o.LocalLength != t.LocalLength)
                return false;
            if (o.iE-o.i0 != t.iE)
                return false;
            return true;
        }

    [Test]
        public static void GetExternalRowsTest(
            ) {
            throw new NotImplementedException();
            //Matlabaufruf --> gesamte Matrix nach Matlab schreiben
            //Teilmatritzen gemäß Globalid extrahieren
            //Mit ExternalRows vergleichen
            //Die große Frage: funktioniert der batchmode connector parallel? Beim rausschreiben beachten
        }

        [Test]
        public static void LocalIndexTest(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder) {

            TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("LocalLIndexTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- Get global index by mapping
            MultigridOperator MGOp = CreateTestMGOperator(UseXdg, DGOrder);
            var map = MGOp.Mapping;
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

        [Test]
        public static void ExternalIndexTest(
            [Values(XDGusage.none, XDGusage.mixed1, XDGusage.mixed2, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder) {

            TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("ExternalIndexTest({0},{1})", UseXdg, DGOrder);

            //Arrange --- Get global index by mapping
            MultigridOperator MGOp = CreateTestMGOperator(UseXdg, DGOrder);
            var map = MGOp.Mapping;
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
            for (int iLoc = 0; iLoc < GlobalIdxMask_ext.Length; iLoc++) {
                Assert.IsTrue(GlobalIdxMask_ext[iLoc] == GlobalIdxMap_ext[iLoc]);
            }
        }

        [Test]
        public static void MapConsistencyTest(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder
            ) {

            //if no selection is chosen mapping should be the same as of origin
            //assume: indexing is right 'cause of other tests
            TestInit((int)UseXdg, DGOrder);
            Console.WriteLine("MapConsistencyTest({0},{1})", UseXdg, DGOrder);

            //Arrange
            MultigridOperator MGOp = CreateTestMGOperator(UseXdg, DGOrder);
            var sbs = new SubBlockSelector(MGOp.Mapping);
            var mask = new TestMasking(sbs);
            var stw = new Stopwatch();
            stw.Restart();

            //Act --- Create Mapping from mask
            stw.Start();
            var submatrix = mask.GetSubBlockMatrix(MGOp.OperatorMatrix);
            stw.Stop();
            var rowpart = submatrix._RowPartitioning;
            var colpart = submatrix._ColPartitioning;

            //Assert --- Equal Partition of mask and origin
            Assert.AreEqual(rowpart, colpart);
            Assert.IsTrue(rowpart.IsLocallyEqual(MGOp.Mapping));
        }

        private static bool[] SetCoupling(MatrixShape shape) {
            bool[] coupling = new bool[3];
            switch (shape) {
                case MatrixShape.diagonal:
                    coupling = new bool[] { true, true, true };
                    break;
                case MatrixShape.diagonal_var:
                    coupling = new bool[] { true, false, true };
                    break;
                case MatrixShape.diagonal_spec:
                    coupling = new bool[] { true, true, false };
                    break;
                case MatrixShape.diagonal_var_spec:
                    coupling = new bool[] { true, false, false };
                    break;
                default:
                    throw new NotSupportedException(String.Format("{0} is not supported by this test", shape));
            }
            return coupling;
        }

        /// <summary>
        /// Test for extracting masked diagonal cell blocks and respective vector matrix multiplication
        /// </summary>
        [Test]
        public static void SubBlockExtractionWithCoupling(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal, MatrixShape.diagonal_var, MatrixShape.diagonal_spec, MatrixShape.diagonal_var_spec)] MatrixShape MShape
            ) {

            TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("ExtractDiagonalBlocks({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- get multigridoperator
            MultigridOperator MGOp = CreateTestMGOperator(UseXdg, DGOrder, MShape);
            BlockMsrMatrix M = MGOp.OperatorMatrix;
            MultigridMapping map = MGOp.Mapping;

            //Arrange --- setup masking
            SubBlockSelector SBS = new SubBlockSelector(map);
            TestMasking mask = new TestMasking(SBS, null);
            bool[] coupling = SetCoupling(MShape);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Arrange --- setup auxiliary matrix
            //this will show us if more is extracted, than it should ...
            var Mprep = new BlockMsrMatrix(map);
            Mprep.Acc(1.0, M);

            //Act --- diagonal subblock extraction
            stw.Start();
            var blocks = mask.GetSubBlocks(Mprep, coupling[0], coupling[1], coupling[2]);
            stw.Stop();

            //Assert --- all diagonal blocks are extracted
            Assert.IsTrue(blocks.Length == map.LocalNoOfBlocks);


            for (int i = 0; i < map.LocalNoOfBlocks; i++) {

                //Arrange --- get ith diagonal block of M: M_i
                int iBlock = i + map.AggGrid.CellPartitioning.i0;
                int L = map.GetBlockLen(iBlock);
                int i0 = map.GetBlockI0(iBlock);
                var Mblock = MultidimensionalArray.Create(L, L);
                M.ReadBlock(i0, i0, Mblock);

                //Act --- M_i-Mones_i
                Mblock.Acc(-1.0, blocks[i]);

                //Assert --- are extracted blocks and 
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

        [Test]
        public static void SubMatrixExtractionWithCoupling(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal, MatrixShape.diagonal_var, MatrixShape.diagonal_spec, MatrixShape.diagonal_var_spec)] MatrixShape MShape
            ) {

            TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("ExtractSubMatrixAndIgnoreCoupling({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- get multigridoperator
            MultigridOperator MGOp = CreateTestMGOperator(UseXdg, DGOrder, MShape);
            BlockMsrMatrix M = MGOp.OperatorMatrix;
            MultigridMapping map = MGOp.Mapping;

            //Arrange --- setup masking
            SubBlockSelector SBS = new SubBlockSelector(map);
            TestMasking mask = new TestMasking(SBS, null);
            bool[] coupling = SetCoupling(MShape);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- establish submatrix
            stw.Start();
            var Mext = mask.GetSubBlockMatrix(M, coupling[0], coupling[1], coupling[2]);
            stw.Stop();
            Mext.Acc(-1.0, M);

            //Assert --- Mext conains only diagonal blocks of M
            Assert.IsTrue(Mext.InfNorm() == 0);
        }

        [Test]
        public static void FastSubMatrixExtraction(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal_var_spec)] MatrixShape MShape
            ) {

            TestInit((int)UseXdg, DGOrder, (int)MShape);
            Console.WriteLine("ExtractSubMatrixAndIgnoreCoupling({0},{1},{2})", UseXdg, DGOrder, MShape);

            //Arrange --- get multigridoperator
            MultigridOperator MGOp = CreateTestMGOperator(UseXdg, DGOrder, MShape);
            BlockMsrMatrix M = MGOp.OperatorMatrix;
            MultigridMapping map = MGOp.Mapping;

            //Arrange --- setup masking
            SubBlockSelector SBS = new SubBlockSelector(map);
            TestMasking mask = new TestMasking(SBS, null);
            bool[] coupling = SetCoupling(MShape);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- establish submatrix
            stw.Start();
            var Mext = mask.GetSubBlockMatrix(M);
            stw.Stop();
            Mext.Acc(-1.0, M);

            //Assert --- Mext conains only diagonal blocks of M
            Assert.IsTrue(Mext.InfNorm() == 0);
        }

        [Test]
        public static void CellBlockVectorOperation(
            [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
            [Values(2)] int DGOrder,
            [Values(MatrixShape.diagonal, MatrixShape.diagonal_var, MatrixShape.diagonal_spec, MatrixShape.diagonal_var_spec)] MatrixShape MShape
            ) {
            //matrix Erzeugung wie in ExtractDiagonalCellBlocks...
            //Auf der HierarchieEbene, auf der Kopplung ausgesetzt wird kann Auswahl vorgenommen werden
            //bei var: 0 / 1, bei DG: <=1 / >1, bei spec: A / B, bei Cells: odd / even
            //accumulierte Teilergebnisse sind dann == fullM*fullX
            double[] Vec;
            var mop = CreateTestMGOperator(out Vec, UseXdg, DGOrder, MShape);
            var M = mop.OperatorMatrix;
            var map = mop.Mapping;

            //Arrange --- setup masking
            SubBlockSelector SBS = new SubBlockSelector(map);
            
            TestMasking mask = new TestMasking(SBS, null);
            bool[] coupling = SetCoupling(MShape);
            var blocks = mask.GetSubBlocks(M, coupling[0], coupling[1], coupling[2]);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- diagonal subblock extraction
            

            //Assert --- all diagonal blocks are extracted
            Assert.IsTrue(blocks.Length == map.LocalNoOfBlocks);


            double[] Res = new double[map.LocalLength];
            double[] MRes = new double[map.LocalLength];

            for (int i = 0; i < map.LocalNoOfBlocks; i++) {
                stw.Start();
                double[] Vec_i = mask.GetVectorCellwise(Vec,i);
                double[] Res_i = new double[Vec_i.Length];
                blocks[i].MatVecMul(1.0,Vec_i,1.0,Res_i);
                mask.AccVecCellwiseToFull(Res_i,i, Res);
                stw.Stop();
            }
            M.SpMV(1.0, Vec, 1.0, MRes);
            Debug.Assert(MRes.L2Norm() != 0);
            Res.AccV(-1.0,MRes);

            //Assert --- are extracted blocks and 
            Assert.IsTrue(Res.L2Norm() == 0.0, String.Format("L2Norm neq 0!"));
        }

        private static void SetDefaultSplitSelection(this SubBlockSelector sbs, MatrixShape shape, bool upper) {
            switch (shape) {
                case MatrixShape.diagonal:
                    sbs.SetDefaultCellSelection(upper);
                    break;
                case MatrixShape.diagonal_var:
                    
                    break;
                case MatrixShape.diagonal_spec:
                    
                    break;
                case MatrixShape.diagonal_var_spec:
                    
                    break;
                default:
                    throw new NotSupportedException(String.Format("{0} is not supported by this test", shape));
            }
        }

        private static void SetDefaultCellSelection(this SubBlockSelector sbs, bool upper) {
            List<int> odds = new List<int>();
            List<int> even = new List<int>();
            for (int i = 0; i < sbs.GetMapping.LocalNoOfBlocks; i++) {
                if (i % 2 != 0)
                    odds.Add(i);
                else
                    even.Add(i);
            }
            sbs.CellSelector(upper?odds:even);
        }

        [Test]
        public static void SplitVectorOperation(
        [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
        [Values(2)] int DGOrder,
        [Values(MatrixShape.diagonal)] MatrixShape MShape
        ) {
            //matrix Erzeugung wie in ExtractDiagonalCellBlocks...
            //Auf der HierarchieEbene, auf der Kopplung ausgesetzt wird kann Auswahl vorgenommen werden
            //bei var: 0 / 1, bei DG: <=1 / >1, bei spec: A / B, bei Cells: odd / even
            //accumulierte Teilergebnisse sind dann == fullM*fullX
            double[] Vec;
            var mop = CreateTestMGOperator(out Vec, UseXdg, DGOrder, MShape);
            var M = mop.OperatorMatrix;
            var map = mop.Mapping;

            //Arrange --- setup masking
            SubBlockSelector sbsA = new SubBlockSelector(map);
            sbsA.SetDefaultSplitSelection(MShape, true);
            TestMasking maskA = new TestMasking(sbsA, null);
            
            SubBlockSelector sbsB = new SubBlockSelector(map);
            sbsB.SetDefaultSplitSelection(MShape, false);
            TestMasking maskB = new TestMasking(sbsB, null);


            var subA = maskA.GetSubBlockMatrix(M);
            var subB = maskB.GetSubBlockMatrix(M);

            //Arrange --- some time measurement
            Stopwatch stw = new Stopwatch();
            stw.Reset();

            //Act --- 
            stw.Start();
            var VecA=maskA.GetSubBlockVec(Vec);
            var VecB = maskB.GetSubBlockVec(Vec);
            stw.Stop();
            
            double[] Res = new double[map.LocalLength];
            double[] ResA = new double[subA.RowPartitioning.LocalLength];
            double[] ResB = new double[subB.RowPartitioning.LocalLength];
            double[] MRes = new double[map.LocalLength];

            subA.SpMV(1.0, VecA, 1.0, ResA);
            subB.SpMV(1.0, VecB, 1.0, ResB);

            maskA.AccVecToFull(ResA, Res);
            maskB.AccVecToFull(ResB, Res);

            M.SpMV(1.0, Vec, 1.0, MRes);
            Debug.Assert(MRes.L2Norm() != 0);
            Res.AccV(-1.0, MRes);

            //Assert --- are extracted blocks and 
            Assert.IsTrue(Res.L2Norm() == 0.0, String.Format("L2Norm neq 0!"));
        }

        public static void MatrixSplitting(
        [Values(XDGusage.none, XDGusage.all)] XDGusage UseXdg,
        [Values(2)] int DGOrder,
        [Values(MatrixShape.diagonal_var_spec)] MatrixShape MShape
        ) {


            //Hier interessiert uns kein ignorecoupling
            //Die Frage ist: werden blöcke abseits der Diagonalen richtig behandelt?
            //Bauen 2er Testmatritzen durch solver:
            //M1: enthält nur Elemente mit Var1 / SpecA / DG<=1 (Einstellen über Manipulation der Flüsse/Sourceterme)
            //M2: keine Einschränkung aber gleiche Flüsse und Sourceterme wie M1
            //Maskierung gemäß M1: nur Elemente mit Var1 / SpecA / DG<=1
            //Test: extract(M2)-M1=zeros
            //Test-crit: Result.InfNorm()==0
            throw new NotImplementedException();
            //Eventuell auch möglich für Zellen

        }
        public static void AccSubVector() {
            throw new NotImplementedException();
            //matrix Erzeugung wie in ExtractSubMatrix ...
            //Auf der HierarchieEbene, auf der Kopplung ausgesetzt wird kann Auswahl vorgenommen werden
            //bei var: 0 / 1, bei DG: <=1 / >1, bei spec: A / B, bei Cells: odd / even
            //accumulierte Teilergebnisse sind dann == fullM*fullX 
        }

        public static void ExtractSubMatrix_Ext() {
            throw new NotImplementedException();
            //Kombination aus GetExternalRows und SubMatrix extraction
            //Es werden nur external Row Einträge berücksichtigt
            //Berechnung eventuell in Matlab durchführen
        }

        public static void AccSubVec_Ext() {
            throw new NotImplementedException();
            //Kombination aus GetExternalRows und AccSubVec_Ext
            //Es werden nur external Row Einträge berücksichtigt
            //Berechnung eventuell in Matlab durchführen
        }
    }
}

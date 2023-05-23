using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
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

namespace AdvancedSolverTests.SubBlocking
{
    public enum SelectionType
    {
        variables,
        species,
        all_combined,
        degrees,
    }

    internal class TestMask : BlockMask {
        public TestMask(SubBlockSelector SBS, BlockMsrMatrix ExtRows) : base(SBS, ExtRows){
            Global_IList_LocalCells = base.GlobalIndices_Internal.ToArray();
            Global_IList_ExternalCells = base.GlobalIndices_External.ToArray();
        }

        public long[] Global_IList_ExternalCells;
        public long[] Global_IList_LocalCells;
    }

    public class MgoSolverPair : IDisposable {

        public MgoSolverPair(SubBlockTestSolver2Var solver) {
            MGOp = solver.MGOp;
            Vec = solver.someVec;
            this.Solver = solver;
            MGSeq = solver.MgSeq;
            ilPSP.Tracing.Tracer.NamespacesToLog = new string[] { "" };
            AssertWatch = new Stopwatch();
            AssertWatch.Start();
        }

        public MultigridOperator MGOp {
            get;
            private set;
        }

        public double[] Vec {
            get;
            private set;
        }

        public AggregationGridData[] MGSeq {
            get;
            private set;
        }

        public SubBlockTestSolver2Var Solver {
            get;
            private set;
        }

        public void Dispose() {
            Solver.Dispose();
            CheckAssertWatch();
        }

        private Stopwatch AssertWatch;

        private void CheckAssertWatch() {
            AssertWatch.Stop();
            double time = AssertWatch.Elapsed.TotalSeconds;
            double timelimit = 240; //sec
            // FK to Jens: what are you doing???
            //Assert.IsTrue(time < timelimit, "time limit of " + timelimit + " seconds exceeded. There is something rotten, plz check ...");
            if (time > timelimit)
                Console.WriteLine("Warning: time limit of " + timelimit + " seconds exceeded. There is something rotten, plz check ...");
        }
    }


    public static class Utils
    {
        
        public static MgoSolverPair CreateTestMGOperator(XDGusage UseXdg = XDGusage.none, int DGOrder = 2, MatrixShape MShape = MatrixShape.full, int Resolution = 4) {
            //using (var solver = new SubBlockTestSolver2Var() { m_UseXdg = UseXdg, m_DGorder = DGOrder, m_Mshape = MShape, m_Res = Resolution }) {
            var solver = new SubBlockTestSolver2Var() { m_UseXdg = UseXdg, m_DGorder = DGOrder, m_Mshape = MShape, m_Res = Resolution };
            ilPSP.Tracing.Tracer.NamespacesToLog = new string[] { "" };
            solver.Init(null);
            ilPSP.Tracing.Tracer.NamespacesToLog = new string[] { "" };
            solver.RunSolverMode();
            ilPSP.Tracing.Tracer.NamespacesToLog = new string[] { "" };
            //}
            // Note to Jens:
            // the "using"-block calls the Dispose()-Method, so the return value may depend on disposed object
            return new MgoSolverPair(solver);
        }




        public static void TestInit(params int[] Testparameters) {
            int L = Testparameters.Length;
            int rank, size;
            var com = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(com, out rank);
            csMPI.Raw.Comm_Size(com, out size);

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
            using(var T1 = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal)) {
                T1.MGOp.OperatorMatrix.SaveToTextFileSparseDebug("M");
            }
            using(var T2 = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_var)) {
                T2.MGOp.OperatorMatrix.SaveToTextFileSparseDebug("M_var");
            }
            using(var T3 = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_spec)) {
                T3.MGOp.OperatorMatrix.SaveToTextFileSparseDebug("M_spec");
            }
            using(var T4 = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_var_spec)) {
                T4.MGOp.OperatorMatrix.SaveToTextFileSparseDebug("M_var_spec");
            }
        }

        public static void SetAll(this BlockMsrMatrix A, double val) {
            var rowmap = A._RowPartitioning;
            var colmap = A._ColPartitioning;
            int RowBlocks = rowmap.LocalNoOfBlocks;
            int ColBlocks = colmap.LocalNoOfBlocks;

            Partitioning rowpart = new Partitioning(RowBlocks);

            for (long iBlock = rowpart.i0; iBlock < rowpart.iE; iBlock++) {
                for (long jBlock = rowpart.i0; jBlock < rowpart.iE; jBlock++) {
                    long i0 = rowmap.GetBlockI0(iBlock);
                    long j0 = colmap.GetBlockI0(jBlock);
                    int iL = rowmap.GetBlockLen(iBlock);
                    int jL = colmap.GetBlockLen(jBlock);
                    var subM = MultidimensionalArray.Create(iL, jL);
                    subM.SetAll(val);
                    A.AccBlock(i0, j0, 1.0, subM);
                }
            }
            double min, max;
            long minc, minr, maxc, maxr;
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
            for (long iBlock = rowpart.i0; iBlock < rowpart.iE; iBlock++) {
                for (long jBlock = rowpart.i0; jBlock < rowpart.iE; jBlock++) {
                    long i0 = rowmap.GetBlockI0(iBlock);
                    long j0 = colmap.GetBlockI0(jBlock);
                    int iL = rowmap.GetBlockLen(iBlock);
                    int jL = colmap.GetBlockLen(jBlock);
                    var subM = MultidimensionalArray.Create(iL, jL);
                    A.ReadBlock(i0, j0, subM);
                    subM.ApplyAll(i => i != 0.0 ? 1 : 0);
                    B.AccBlock(i0, j0, 1.0, subM);
                }
            }
            double min, max;
            long minc, minr, maxc, maxr;
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
            if (o.iE - o.i0 != t.iE)
                return false;
            return true;
        }

        public static bool[] SetCoupling(MatrixShape shape) {
            bool[] coupling = new bool[2];
            switch (shape) {
                case MatrixShape.diagonal:
                case MatrixShape.full:
                    coupling = new bool[] { true, true };
                    break;
                case MatrixShape.diagonal_var:
                case MatrixShape.full_var:
                    coupling = new bool[] { false, true };
                    break;
                case MatrixShape.diagonal_spec:
                case MatrixShape.full_spec:
                    coupling = new bool[] { true, false };
                    break;
                case MatrixShape.diagonal_var_spec:
                case MatrixShape.full_var_spec:
                    coupling = new bool[] { false, false };
                    break;
                default:
                    throw new NotSupportedException(String.Format("{0} is not supported by this test", shape));
            }
            return coupling;
        }


        public static double[] GetRandomVector(int Length) {
            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);
            double[] vec = new double[Length];
            var rndgen = new Random(rank);
            for (int i = 0; i < Length; i++) {
                vec[i] = rndgen.NextDouble();
                //vec[i] = (rank + 1);
                Debug.Assert(vec[i]!=0);
            }
            return vec;
        }

        public static void SetDefaultSplitSelection(this SubBlockSelector sbs, MatrixShape shape, bool upper, bool islocal=true) {
            switch (shape) {
                case MatrixShape.diagonal:
                case MatrixShape.full:
                    sbs.DefaultCellSplit(upper, islocal);
                    break;
                case MatrixShape.diagonal_var:
                case MatrixShape.full_var:
                    sbs.DefaultSpeciesSplit(upper);
                    if(!islocal) sbs.AllExternalCellsSelection();
                    break;
                case MatrixShape.diagonal_spec:
                case MatrixShape.full_spec:
                    sbs.DefaultVarSplit(upper);
                    if (!islocal) sbs.AllExternalCellsSelection();
                    break;
                case MatrixShape.diagonal_var_spec:
                case MatrixShape.full_var_spec:
                    sbs.DefaultCellSplit(upper, islocal);
                    break;
                default:
                    throw new NotSupportedException(String.Format("{0} is not supported by this test", shape));
            }
        }

        private static void DefaultSpeciesSplit(this SubBlockSelector sbs, bool upper) {
            if (sbs.Mapping.IsXDGvariable(0)) {
                SpeciesId[] SIdc = sbs.Mapping.UsedSpecies;
                SpeciesId[] OtherSpec = (SIdc.Length - 1).ForLoop(i => SIdc[i + 1]);
                if (upper)
                    sbs.SetSpeciesSelector(SIdc[0]);
                else
                    sbs.SetSpeciesSelector(OtherSpec);
            } else {
                throw new NotSupportedException();
            }
        }

        private static void DefaultVarSplit(this SubBlockSelector sbs, bool upper) {
            sbs.SetVariableSelector(upper ? 0 : 1);
            int NoOfVar = sbs.Mapping.NoOfVariables;
            int[] OtherVars = (NoOfVar - 1).ForLoop(i => i + 1);
            if (upper)
                sbs.SetVariableSelector(0);
            else
                sbs.SetVariableSelector(OtherVars);
        }

        private static void DefaultCellSplit(this SubBlockSelector sbs, bool upper, bool islocal = true) {
            List<int> odds = new List<int>();
            List<int> even = new List<int>();
            int i0, iE;
            if(islocal) {
                i0 = 0;
                iE = sbs.Mapping.LocalNoOfBlocks;
            } else {
                i0 = sbs.Mapping.LocalNoOfBlocks;
                iE = sbs.Mapping.NoOfExternalCells + sbs.Mapping.LocalNoOfBlocks;
            }

            for(int i = i0; i < iE; i++) {
                if(i % 2 != 0)
                    odds.Add(i);
                else
                    even.Add(i);
            }
            sbs.CellSelector(upper ? odds : even, false);
        }


        public static void GetDefaultSelection(this SubBlockSelector sbs, SelectionType SType, int iCell) {
            SpeciesId A = sbs.Mapping.UsedSpecies[0];
            SpeciesId B = sbs.Mapping.UsedSpecies[1];

            sbs.CellSelector(iCell);
            //do not change this, selection corresponds to hardcoded masking
            //see GetSubIndices
            switch (SType) {
                case SelectionType.degrees:
                    sbs.SetModeSelector(p => p == 1);
                    break;
                case SelectionType.species:
                    sbs.SetSpeciesSelector(A);
                    break;
                case SelectionType.variables:
                    sbs.SetVariableSelector(1);
                    break;
                case SelectionType.all_combined:
                    sbs.SetModeSelector(p => p == 1);
                    sbs.SetSpeciesSelector(A);
                    sbs.SetVariableSelector(1);
                    break;
            }
        }

        public static BlockMsrMatrix GetCellCompMatrix(SelectionType SType, MultigridOperator mop, int iB) {

            int rank = mop.Mapping.MpiRank;
            long iBlock = mop.Mapping.AggGrid.CellPartitioning.i0 + iB;
            long i0 = mop.Mapping.GetBlockI0(iBlock);

            //int jBlock = i0 + jB;
            int R = mop.Mapping.GetBlockLen(iBlock);
            //int C = mop.Mapping.GetBlockLen(jBlock);

            bool ZwoSpecR = Math.Max(mop.Mapping.GetSubblkLen(0)[0], mop.Mapping.GetSubblkLen(1)[0]) == R;
            //bool ZwoSpecC = (mop.Mapping.AggBasis[0].GetMaximalLength(DGdegree) + mop.Mapping.AggBasis[0].GetMaximalLength(DGdegree - 1)) == C;

            SpeciesId A = ((XdgAggregationBasis)mop.Mapping.AggBasis[0]).UsedSpecies[0];
            int Specpos = ((XdgAggregationBasis)mop.Mapping.AggBasis[0]).GetSpeciesIndex(iB, A);

            long[] SubIdcR = GetSubIndices(SType, ZwoSpecR, Specpos).Select((int i) => (long)i).ToArray();
            //int[] SubIdcC = GetSubIndices(SType, ZwoSpecC);

            for (int i = 0; i < SubIdcR.Length; i++) {
                Debug.Assert(SubIdcR[i] < R);
                SubIdcR[i] += i0;
            }
            //for (int i = 0; i < SubIdcC.Length; i++) {
            //    Debug.Assert(SubIdcC[i] < C);
            //    SubIdcC[i] += i0;
            //}
            //return mop.OperatorMatrix.GetSubMatrix(SubIdcR, SubIdcC);
            var part = new BlockPartitioning(SubIdcR.Length, new long[] { 0 }, new int[] { SubIdcR.Length }, csMPI.Raw._COMM.SELF);

            BlockMsrMatrix sub = new BlockMsrMatrix(part);

            mop.OperatorMatrix.WriteSubMatrixTo(sub, SubIdcR, default(long[]), SubIdcR, default(long[]));
            return sub;

            //return mop.OperatorMatrix.GetSubMatrix(SubIdcR, SubIdcR);
        }

        public static int[] GetSubIndices(SelectionType SType, bool ZwoSpec, int lookatpos) {
            int[] SubIdc;
            switch (SType) {
                case SelectionType.degrees:
                    //dg==1
                    SubIdc = ZwoSpec ? new int[] { 1, 2, 7, 8, 13, 14, 16, 17 } : new int[] { 1, 2, 7, 8 };
                    break;
                case SelectionType.variables:
                    //var==u2
                    SubIdc = ZwoSpec ? new int[] { 12, 13, 14, 15, 16, 17 } : new int[] { 6, 7, 8 };
                    break;
                case SelectionType.species:
                    switch (lookatpos) {
                        case 0:
                            SubIdc = ZwoSpec ? new int[] { 0, 1, 2, 3, 4, 5, 12, 13, 14 } : new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
                            break;
                        case 1:
                            SubIdc = ZwoSpec ? new int[] { 6, 7, 8, 9, 10, 11, 15, 16, 17 } : new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
                            break;
                        default:
                            throw new NotSupportedException("This selection type is not supported");
                    }
                    break;
                case SelectionType.all_combined:
                    switch (lookatpos) {
                        case 0:
                            SubIdc = ZwoSpec ? new int[] { 13, 14 } : new int[] { 7, 8 };
                            break;
                        case 1:
                            SubIdc = ZwoSpec ? new int[] { 16, 17 } : new int[] { 7, 8 };
                            break;
                        default:
                            throw new NotSupportedException("This selection type is not supported");
                    }
                    break;
                default:
                    throw new NotSupportedException("This selection type is not supported");
            }
            return SubIdc;
        }

        public static int[] GetAllExternalCells(ICoordinateMapping map) {
            int NoOfExternalCells = map.NoOfExternalCells;
            int offset = map.NoOfLocalUpdatedCells;
            int[] extcells = NoOfExternalCells.ForLoop(i => i + offset);
            return extcells;
        }

        public static int[] AllExternalCellsSelection(this SubBlockSelector sbs) {
            var map = sbs.Mapping;
            var extcells = GetAllExternalCells(map);
            sbs.CellSelector(extcells, false);
            return extcells;
        }

        public static long[] GetIndcOfExtCell(MultigridMapping map, int jCell) {
            int Jup = map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            long i0 = map.GlobalUniqueIndex(0, jCell, 0);
            int fld = map.NoOfVariables;
            int N = 0;
            for (int iF = 0; iF < fld; iF++)
                N+=map.AggBasis[iF].GetLength(jCell, map.DgDegree[iF]);
            long[] ret = N.ForLoop(i => i + i0);
            return ret;
        }

        public static long[] GetAllExtCellIdc(MultigridMapping map) {
            var extC = GetAllExternalCells(map);
            List<long> extIdcL = new List<long>();
            foreach (int eC in extC) {
                Debug.Assert(eC < map.AggGrid.iLogicalCells.NoOfExternalCells + map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells);
                Debug.Assert(eC >= map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells);
                long[] Idc = GetIndcOfExtCell(map, eC);
                extIdcL.AddRange(Idc);
            }
            long[] extIdc = extIdcL.ToArray();
            Array.Sort(extIdc);
#if DEBUG
            int[] fields = map.NoOfVariables.ForLoop(i => i);
            long[] GlobalIdxMap_ext = map.GetSubvectorIndices_Ext(fields);
            Array.Sort(GlobalIdxMap_ext);
            Debug.Assert(extIdc.Length == GlobalIdxMap_ext.Length);
            for(int iCell=0;iCell<extIdc.Length; iCell++) {
                Debug.Assert(extIdc[iCell] == GlobalIdxMap_ext[iCell]);
            }
#endif
            return extIdc;
        }

        public static Dictionary<int,long[]> GetDictOfAllExtCellIdc(MultigridMapping map) {
            int[] cells = GetAllExternalCells(map);
            Dictionary<int, long[]> extDict = new Dictionary<int, long[]>();

            foreach(int eC in cells) {
                long[] idc = GetIndcOfExtCell(map, eC);
                extDict.Add(eC,idc);
            }
            return extDict;
        }

        public static long[] GimmeAllBlocksWithSpec(MultigridMapping map, int DOF) {
            long Bi0 = map.AggGrid.CellPartitioning.i0;
            List<long> hits = new List<long>();
            for(int i=0; i < map.LocalNoOfBlocks; i++) {
                long iBlock = Bi0 + i;
                int type = map.GetBlockType(iBlock);
                int BlockDOF = map.GetSubblkLen(type)[0];
                if (BlockDOF == DOF)
                    hits.Add(iBlock);
            }

            int[] recvcnt = hits.Count().MPIAllGather();
            Debug.Assert(recvcnt.Length==map.MpiSize);
            long[] VecOfAllhits = hits.ToArray().MPIAllGatherv(recvcnt);
            
            return VecOfAllhits;
        }

        public static MsrMatrix ConvertToMsr(this MultidimensionalArray M) {
            Partitioning rowpart = new Partitioning(M.Lengths[0], MPI.Wrappers.csMPI.Raw._COMM.SELF);
            Partitioning colpart = new Partitioning(M.Lengths[1], MPI.Wrappers.csMPI.Raw._COMM.SELF);
            var bla = new MsrMatrix(rowpart, colpart);
            bla.AccBlock(0, 0, 1.0, M);
            return bla;
        }

        public static BlockMsrMatrix ConvertToQuadraticBMsr(this BlockMsrMatrix M, long[] Colidx, bool isinternal) {

            Debug.Assert(M._RowPartitioning.LocalLength == Colidx.Length);

            int NoOfBlocks = M._RowPartitioning.LocalNoOfBlocks;
            long[] Offsets = new long[NoOfBlocks];
            int[] Lengths = new int[NoOfBlocks];
            long IdxOffset = M._RowPartitioning.i0;
            long ColIdxOffset = M._ColPartitioning.i0;
            for(int i = 0; i < NoOfBlocks; i++) {
                long iBlock = i + M._RowPartitioning.FirstBlock;
                Offsets[i] = M._RowPartitioning.GetBlockI0(iBlock) - IdxOffset;
                Lengths[i] = M._RowPartitioning.GetBlockLen(iBlock);
            }
            BlockPartitioning part = new BlockPartitioning(M._RowPartitioning.LocalLength, Offsets, Lengths, csMPI.Raw._COMM.SELF, true);
            BlockMsrMatrix ret = new BlockMsrMatrix(part);

            long[] RowISrc = M._RowPartitioning.LocalLength.ForLoop(i => i + IdxOffset);
            //int[] ColISrc = M._ColPartitioning.LocalLength.ForLoop(i => Colidx[i]);
            //if (ColISrc.Length < RowISrc.Length)
            //    ExtISrc = (RowISrc.Length - ColISrc.Length).ForLoop(i => Colidx[i+ ColISrc.Length]);
            long[] ExtISrc = M._RowPartitioning.LocalLength.ForLoop(i => (long)Colidx[i]);
            long[] ExtITrg = M._RowPartitioning.LocalLength.ForLoop(i => (long)i);
            if(isinternal) {
                M.AccSubMatrixTo(1.0, ret, RowISrc, default(long[]), ExtISrc, default(long[]));
            } else {
                M.AccSubMatrixTo(1.0, ret, RowISrc, default(long[]), new long[0], default(long[]), ExtISrc, ExtITrg);
            }

            return ret;
        }

        public static int[] GetCellsOfOverlappingTestBlock(MultigridMapping map) {
            List<int> testcells = new List<int>();

            //find block with 4 neighbours
            for (int iCell = 0; iCell < map.LocalNoOfBlocks; iCell++) {
                int[] NC = map.AggGrid.iLogicalCells.CellNeighbours[iCell];
                if (NC.Length == 4) {
                    foreach (int c in NC) {
                        if (c >= map.LocalNoOfBlocks) {
                            testcells.Add(c);
                        } else if (!map.IsLocalBlock(map.FirstBlock + c)) {
                            testcells.Add(c);
                        }
                    }
                    if (testcells.Count > 0) {
                        testcells.Add(iCell);
                        break;
                    }

                }
            }
            Debug.Assert(testcells.Count > 0);
            return testcells.ToArray();
        }

        public static long[] GetIdcOfSubBlock(MultigridMapping map, int[] cells) {
            List<long> idc = new List<long>();
            foreach (int c in cells)
                idc.AddRange(GetIndcOfExtCell(map,c));
            return idc.ToArray();
        }
        public static int GetIdxOfFirstBlockWith(MultigridMapping map, bool ZwoSpec) {
            int maxLen = Math.Max(map.GetSubblkLen(0)[0], map.GetSubblkLen(1)[0]);
            int minLen = Math.Min(map.GetSubblkLen(0)[0], map.GetSubblkLen(1)[0]);
            int crit = ZwoSpec ? maxLen : minLen;
            for (int iCell = 0; iCell < map.LocalNoOfBlocks; iCell++) {
                long iBlock = iCell + map.AggGrid.CellPartitioning.i0;
                if(map.GetBlockLen(iBlock) == crit)
                    return iCell;
            }
            return -1;
        } 
    }
}

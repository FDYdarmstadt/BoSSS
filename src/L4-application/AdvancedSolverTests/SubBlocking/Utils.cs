using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MPI.Wrappers;
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
            Global_IList_LocalCells=base.GlobalIList_Internal.ToArray();
            Global_IList_ExternalCells=base.GlobalIList_External.ToArray();
        }

        public int[] Global_IList_ExternalCells;
        public int[] Global_IList_LocalCells;
    }

    internal static class Utils
    {
        public static MultigridOperator CreateTestMGOperator(XDGusage UseXdg = XDGusage.none, int DGOrder = 2, MatrixShape MShape = MatrixShape.full, int Resolution = 4) {
            return CreateTestMGOperator(out double[] Vec, UseXdg, DGOrder, MShape, Resolution);
        }

        public static MultigridOperator CreateTestMGOperator(out double[] Vec, XDGusage UseXdg = XDGusage.none, int DGOrder = 2, MatrixShape MShape = MatrixShape.full, int Resolution = 4) {
            MultigridOperator retMGOp;
            using (var solver = new SubBlockTestSolver2Var() { m_UseXdg = UseXdg, m_DGorder = DGOrder, m_Mshape = MShape, m_Res = Resolution }) {
                solver.Init(null);
                solver.RunSolverMode();
                retMGOp = solver.MGOp;
                Vec = solver.someVec;
            }
            return retMGOp;
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
            mgo = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal);
            mgo.OperatorMatrix.SaveToTextFileSparseDebug("M");
            mgo = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_var);
            mgo.OperatorMatrix.SaveToTextFileSparseDebug("M_var");
            mgo = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_spec);
            mgo.OperatorMatrix.SaveToTextFileSparseDebug("M_spec");
            mgo = CreateTestMGOperator(UseXdg: XDGusage.all, MShape: MatrixShape.diagonal_var_spec);
            mgo.OperatorMatrix.SaveToTextFileSparseDebug("M_var_spec");
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
            if (o.iE - o.i0 != t.iE)
                return false;
            return true;
        }

        public static bool[] SetCoupling(MatrixShape shape) {
            bool[] coupling = new bool[3];
            switch (shape) {
                case MatrixShape.diagonal:
                case MatrixShape.full:
                    coupling = new bool[] { true, true, true };
                    break;
                case MatrixShape.diagonal_var:
                case MatrixShape.full_var:
                    coupling = new bool[] { true, false, true };
                    break;
                case MatrixShape.diagonal_spec:
                case MatrixShape.full_spec:
                    coupling = new bool[] { true, true, false };
                    break;
                case MatrixShape.diagonal_var_spec:
                case MatrixShape.full_var_spec:
                    coupling = new bool[] { true, false, false };
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
            var rndgen = new Random();
            for (int i = 0; i < Length; i++) {
                //vec[i] = rndgen.NextDouble() * (rank+1);
                vec[i] = (rank + 1);
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
            if (sbs.GetMapping.AggBasis[0] is XdgAggregationBasis) {
                SpeciesId[] SIdc = ((XdgAggregationBasis)sbs.GetMapping.AggBasis[0]).UsedSpecies;
                SpeciesId[] OtherSpec = (SIdc.Length - 1).ForLoop(i => SIdc[i + 1]);
                if (upper)
                    sbs.SpeciesSelector(SIdc[0]);
                else
                    sbs.SpeciesSelector(OtherSpec);
            }
        }

        private static void DefaultVarSplit(this SubBlockSelector sbs, bool upper) {
            sbs.VariableSelector(upper ? 0 : 1);
            int NoOfVar = sbs.GetMapping.NoOfVariables;
            int[] OtherVars = (NoOfVar - 1).ForLoop(i => i + 1);
            if (upper)
                sbs.VariableSelector(0);
            else
                sbs.VariableSelector(OtherVars);
        }

        private static void DefaultCellSplit(this SubBlockSelector sbs, bool upper, bool islocal=true) {
            List<int> odds = new List<int>();
            List<int> even = new List<int>();
            int i0, iE;
            if (islocal) {
                i0 = 0;
                iE = sbs.GetMapping.LocalNoOfBlocks;
            } else {
                i0 = sbs.GetMapping.LocalNoOfBlocks;
                iE = sbs.GetMapping.AggGrid.iLogicalCells.NoOfExternalCells+ sbs.GetMapping.LocalNoOfBlocks;
            }

            for (int i = i0; i < iE; i++) {
                if (i % 2 != 0)
                    odds.Add(i);
                else
                    even.Add(i);
            }
            sbs.CellSelector(upper ? odds : even, false);
        }


        public static void GetDefaultSelection(this SubBlockSelector sbs, SelectionType SType, int iCell) {
            SpeciesId B = ((XdgAggregationBasis)sbs.GetMapping.AggBasis[0]).UsedSpecies[0];
            SpeciesId A = ((XdgAggregationBasis)sbs.GetMapping.AggBasis[0]).UsedSpecies[1];
            sbs.CellSelector(iCell);
            //do not change this, selection corresponds to hardcoded masking
            //see GetSubIndices
            switch (SType) {
                case SelectionType.degrees:
                    sbs.ModeSelector(p => p == 1);
                    break;
                case SelectionType.species:
                    sbs.SpeciesSelector(B);
                    break;
                case SelectionType.variables:
                    sbs.VariableSelector(1);
                    break;
                case SelectionType.all_combined:
                    sbs.ModeSelector(p => p == 1);
                    sbs.SpeciesSelector(B);
                    sbs.VariableSelector(1);
                    break;
            }
        }

        public static BlockMsrMatrix GetCellCompMatrix(SelectionType SType, MultigridOperator mop, int iB) {
            int DGdegree = 2;

            int rank = mop.Mapping.MpiRank;
            int iBlock = mop.Mapping.AggGrid.CellPartitioning.i0 + iB;
            int i0 = mop.Mapping.GetBlockI0(iBlock);

            //int jBlock = i0 + jB;
            int R = mop.Mapping.GetBlockLen(iBlock);
            //int C = mop.Mapping.GetBlockLen(jBlock);

            bool ZwoSpecR = (mop.Mapping.AggBasis[0].GetMaximalLength(DGdegree) + mop.Mapping.AggBasis[0].GetMaximalLength(DGdegree - 1)) == R;
            //bool ZwoSpecC = (mop.Mapping.AggBasis[0].GetMaximalLength(DGdegree) + mop.Mapping.AggBasis[0].GetMaximalLength(DGdegree - 1)) == C;

            int[] SubIdcR = GetSubIndices(SType, ZwoSpecR);
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
            return mop.OperatorMatrix.GetSubMatrix(SubIdcR, SubIdcR);
        }

        public static int[] GetSubIndices(SelectionType SType, bool ZwoSpec) {
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
                    //spec==B
                    SubIdc = ZwoSpec ? new int[] { 6, 7, 8, 9, 10, 11, 15, 16, 17 } : new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
                    break;
                case SelectionType.all_combined:
                    SubIdc = ZwoSpec ? new int[] { 16, 17 } : new int[] { 7, 8 };
                    break;
                default:
                    throw new NotSupportedException("This selection type is not supported");
            }
            return SubIdc;
        }

        public static int[] GetAllExternalCells(MultigridMapping map) {
            int NoOfExternalCells = map.AggGrid.iLogicalCells.NoOfExternalCells;
            int offset = map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int[] extcells = NoOfExternalCells.ForLoop(i => i + offset);
            return extcells;
        }

        public static int[] AllExternalCellsSelection(this SubBlockSelector sbs) {
            var map = sbs.GetMapping;
            var extcells=GetAllExternalCells(map);
            sbs.CellSelector(extcells,false);
            return extcells;
        }

        public static int[] GetIndcOfExtCell(MultigridMapping map, int jCell) {
            int Jup = map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int i0 = map.GlobalUniqueIndex(0, jCell, 0);
            int fld = map.NoOfVariables;
            int N = 0;
            for (int iF = 0; iF < fld; iF++)
                N+=map.AggBasis[iF].GetLength(jCell, map.DgDegree[iF]);
            int[] ret = N.ForLoop(i => i + i0);
            return ret;
        }

        public static int[] GetAllExtCellIdc(MultigridMapping map) {
            var extC = GetAllExternalCells(map);
            List<int> extIdcL = new List<int>();
            foreach (int eC in extC) {
                Debug.Assert(eC < map.AggGrid.iLogicalCells.NoOfExternalCells + map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells);
                Debug.Assert(eC >= map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells);
                int[] Idc = GetIndcOfExtCell(map, eC);
                extIdcL.AddRange(Idc);
            }
            int[] extIdc = extIdcL.ToArray();
            Array.Sort(extIdc);
#if Debug
            int[] fields = map.NoOfVariables.ForLoop(i => i);
            int[] GlobalIdxMap_ext = map.GetSubvectorIndices_Ext(fields);
            Array.Sort(GlobalIdxMap_ext);
            Debug.Assert(extIdc.Length == GlobalIdxMap_ext.Length);
            for(int iCell=0;iCell<extIdc.Length; iCell++) {
                Debug.Assert(extIdc[iCell] == GlobalIdxMap_ext[iCell]);
            }
#endif
            return extIdc;
        }

        public static Dictionary<int,int[]> GetDictOfAllExtCellIdc(MultigridMapping map) {
            int[] cells = GetAllExternalCells(map);
            Dictionary<int, int[]> extDict = new Dictionary<int, int[]>();

            foreach(int eC in cells) {
                int[] idc = GetIndcOfExtCell(map, eC);
                extDict.Add(eC,idc);
            }
            return extDict;
        }

        public static int[] GimmeAllBlocksWithSpec(MultigridMapping map, int DOF) {
            int Bi0 = map.AggGrid.CellPartitioning.i0;
            List<int> hits = new List<int>();
            for(int i=0; i < map.LocalNoOfBlocks; i++) {
                int iBlock = Bi0 + i;
                int type = map.GetBlockType(iBlock);
                int BlockDOF = map.GetSubblkLen(type)[0];
                if (BlockDOF == DOF)
                    hits.Add(iBlock);
            }

            int[] recvcnt = hits.Count().MPIAllGather();
            Debug.Assert(recvcnt.Length==map.MpiSize);
            int[] VecOfAllhits = hits.ToArray().MPIAllGatherv(recvcnt);
            
            return VecOfAllhits;
        }

        public static MsrMatrix ConvertToMsr(this MultidimensionalArray M) {
            Partitioning rowpart = new Partitioning(M.Lengths[0], MPI.Wrappers.csMPI.Raw._COMM.SELF);
            Partitioning colpart = new Partitioning(M.Lengths[1], MPI.Wrappers.csMPI.Raw._COMM.SELF);
            var bla = new MsrMatrix(rowpart, colpart);
            bla.AccBlock(0, 0, 1.0, M);
            return bla;
        }

        public static BlockMsrMatrix ConvertToQuadraticBMsr(this BlockMsrMatrix M, int[] Colidx) {

            Debug.Assert(M._RowPartitioning.LocalLength == Colidx.Length);
            
            int NoOfBlocks = M._RowPartitioning.LocalNoOfBlocks;
            int[] Offsets = new int[NoOfBlocks];
            int[] Lengths = new int[NoOfBlocks];
            int IdxOffset = M._RowPartitioning.i0;
            int ColIdxOffset = M._ColPartitioning.i0;
            for (int i=0;i < NoOfBlocks; i++) {
                int iBlock = i + M._RowPartitioning.FirstBlock;
                Offsets[i]=M._RowPartitioning.GetBlockI0(iBlock) - IdxOffset;
                Lengths[i]=M._RowPartitioning.GetBlockLen(iBlock);
            }
            BlockPartitioning part = new BlockPartitioning(M._RowPartitioning.LocalLength,Offsets,Lengths,csMPI.Raw._COMM.SELF,true);
            BlockMsrMatrix ret = new BlockMsrMatrix(part);

            int[] RowISrc = M._RowPartitioning.LocalLength.ForLoop(i=> i + IdxOffset);
            //int[] ColISrc = M._ColPartitioning.LocalLength.ForLoop(i => Colidx[i]);
            //if (ColISrc.Length < RowISrc.Length)
            //    ExtISrc = (RowISrc.Length - ColISrc.Length).ForLoop(i => Colidx[i+ ColISrc.Length]);
            int[] ExtISrc = M._RowPartitioning.LocalLength.ForLoop(i => Colidx[i]);
            int[] ExtITrg = M._RowPartitioning.LocalLength.ForLoop(i => i);

            M.AccSubMatrixTo(1.0,ret, RowISrc, default(int[]), new int[0], default(int[]), ExtISrc, ExtITrg);
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

        public static int[] GetIdcOfSubBlock(MultigridMapping map, int[] cells) {
            List<int> idc = new List<int>();
            foreach (int c in cells)
                idc.AddRange(GetIndcOfExtCell(map,c));
            return idc.ToArray();
        }
           
    }
}

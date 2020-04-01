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
        public static MultigridOperator CreateTestMGOperator(XDGusage UseXdg = XDGusage.none, int DGOrder = 2, MatrixShape MShape = MatrixShape.full, int Resolution = 2) {
            return CreateTestMGOperator(out double[] Vec, UseXdg, DGOrder, MShape, Resolution);
        }

        public static MultigridOperator CreateTestMGOperator(out double[] Vec, XDGusage UseXdg = XDGusage.none, int DGOrder = 2, MatrixShape MShape = MatrixShape.full, int Resolution = 2) {
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
            double[] vec = new double[Length];
            var rndgen = new Random();
            for (int i = 0; i < Length; i++) {
                vec[i] = rndgen.NextDouble();
            }
            return vec;
        }

        public static void SetDefaultSplitSelection(this SubBlockSelector sbs, MatrixShape shape, bool upper) {
            switch (shape) {
                case MatrixShape.diagonal:
                case MatrixShape.full:
                    sbs.DefaultCellSplit(upper);
                    break;
                case MatrixShape.diagonal_var:
                case MatrixShape.full_var:
                    sbs.DefaultSpeciesSplit(upper);
                    break;
                case MatrixShape.diagonal_spec:
                case MatrixShape.full_spec:
                    sbs.DefaultVarSplit(upper);
                    break;
                case MatrixShape.diagonal_var_spec:
                case MatrixShape.full_var_spec:
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

        private static void DefaultCellSplit(this SubBlockSelector sbs, bool upper) {
            List<int> odds = new List<int>();
            List<int> even = new List<int>();
            for (int i = 0; i < sbs.GetMapping.LocalNoOfBlocks; i++) {
                if (i % 2 != 0)
                    odds.Add(i);
                else
                    even.Add(i);
            }
            sbs.CellSelector(upper ? odds : even);
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

        public static void AllExternalCellsSelection(this SubBlockSelector sbs) {
            var map = sbs.GetMapping;
            int NoOfExternalCells = map.AggGrid.iLogicalCells.NoOfExternalCells;
            int offset = map.AggGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int[] extcells = NoOfExternalCells.ForLoop(i => i + offset);
            sbs.CellSelector(extcells,false);
        }
    }
}

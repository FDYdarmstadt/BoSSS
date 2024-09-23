using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace StokesHelical_Ak {
    public partial class HelicalMain : Application<HelicalControl> {

        /*
        BlockMsrMatrix DoBlockPrecond(BlockMsrMatrix Mtx) {

            AggregationGridBasis[][] MgBasis;
            {
                var MgSequence = this.MultigridSequence;
                if (MgSequence.Length <= 0)
                    throw new NotSupportedException();

                MgBasis = new AggregationGridBasis[MgSequence.Length][];
                Basis[] OrgBasisS = this.UnknownsMap.BasisS.ToArray();
                int NoOfVars = OrgBasisS.Length;

                Basis maxDGbasis = null;
                foreach (var b in OrgBasisS) {
                    if (maxDGbasis == null || maxDGbasis.Degree < b.Degree)
                        maxDGbasis = b;
                }

                //AggregationGridBasis[] mgDgB = MgSequence.Select(aggGrd => new AggregationGridBasis(maxDGbasis, aggGrd)).ToArray();

                AggregationGridBasis[][] mgDgB = AggregationGridBasis.CreateSequence(MgSequence, new Basis[] { maxDGbasis });

                for (int iLevel = 0; iLevel < MgSequence.Length; iLevel++) {
                    MgBasis[iLevel] = new AggregationGridBasis[NoOfVars];
                    for (int r = 0; r < MgBasis[iLevel].Length; r++) {
                        MgBasis[iLevel][r] = mgDgB[iLevel][0];
                    }
                }
            }



            MultigridOperator mgOp = new MultigridOperator(MgBasis, this.CurrentSolution.Mapping,
                Mtx, null, this.MultigridOperatorConfig, null);

            return mgOp.OperatorMatrix;
        }

        /*
        protected void MatrixAnalysis(BlockMsrMatrix RawMatrix) {
            using (new FuncTrace()) {

                BlockMsrMatrix PcdMatrix = DoBlockPrecond(RawMatrix);

                if (Globals.writeMtxBeforeR0fix == true) {
                    PcdMatrix.SaveToTextFileSparse("PCDmatrixBoSSSbeforeR0Fix.txt");
                } else {
                    PcdMatrix.SaveToTextFileSparse("PCDmatrixBoSSSafterR0Fix.txt");
                }

                MultidimensionalArray output = MultidimensionalArray.Create(10, 1);

                CellMask myInnerCells = this.GridData.GetBoundaryCells().Complement();


                SubGrid InnerCells = new SubGrid(myInnerCells);
                SubGrid OneCell = new SubGrid(
                    new CellMask(this.GridData,
                        new Chunk[] {
                            Chunk.GetSingleElementChunk(myInnerCells.First().i0)
                        }));

                ///
                // Full_0Vars ist Matrix der SecondOrder Viscosity Terms bei der ur mit r-Impuls koppelt, alle Zellen mit Rand
                // Full_12Vars ist Matrix der SecondOrder Viscosity Terms bei der uxi, ueta mit z-Impuls und eta-Impuls koppelt, alle Zellen mit Rand
                // Full_012Vars ist Matrix der SecondOrder Viscosity Terms (ganze Matrix) alle Zellen mit Rand
                // Inner: s.oben nur innere Zellen
                // OneCell: s.oben nur eine Zelle, nicht am Rand
                double[] Full_0Vars = this.UnknownsMap.GetSubvectorIndices(true, new int[] { 0 }).Select(i => i + 1.0).ToArray();
                double[] Full_12Vars = this.UnknownsMap.GetSubvectorIndices(true, new int[] { 1, 2 }).Select(i => i + 1.0).ToArray();
                double[] Full_012Vars = this.UnknownsMap.GetSubvectorIndices(true, new int[] { 0, 1, 2 }).Select(i => i + 1.0).ToArray();
                double[] Full_AllVars = this.UnknownsMap.GetSubvectorIndices(true, new int[] { 0, 1, 2, 3 }).Select(i => i + 1.0).ToArray();

                double[] Inner_0Vars = this.UnknownsMap.GetSubvectorIndices(InnerCells, true, new int[] { 0 }).Select(i => i + 1.0).ToArray();
                double[] Inner_12Vars = this.UnknownsMap.GetSubvectorIndices(InnerCells, true, new int[] { 1, 2 }).Select(i => i + 1.0).ToArray();
                double[] Inner_012Vars = this.UnknownsMap.GetSubvectorIndices(InnerCells, true, new int[] { 0, 1, 2 }).Select(i => i + 1.0).ToArray();

                double[] OneCell_0Vars = this.UnknownsMap.GetSubvectorIndices(OneCell, true, new int[] { 0 }).Select(i => i + 1.0).ToArray();
                double[] OneCell_12Vars = this.UnknownsMap.GetSubvectorIndices(OneCell, true, new int[] { 1, 2 }).Select(i => i + 1.0).ToArray();
                double[] OneCell_012Vars = this.UnknownsMap.GetSubvectorIndices(OneCell, true, new int[] { 0, 1, 2 }).Select(i => i + 1.0).ToArray();

                Console.WriteLine("Starting Evaluation of Condition in Matlab");
                //BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave_cygwin;
                string[] names = new string[] { "Full_0Vars", "Full_12Vars", "Full_012Vars", "Inner_0Vars", "Inner_12Vars", "Inner_012Vars", "OneCell_0Vars", "OneCell_12Vars", "OneCell_012Vars", "Full_AllVars" };
                using (var bmc = new BatchmodeConnector()) {
                    Console.WriteLine($"Using {bmc.m_Flav.ToString()}");
                    bmc.PutSparseMatrix(PcdMatrix, "Mtx");
                    bmc.PutVector(Full_0Vars, "Full_0Vars");
                    bmc.PutVector(Full_12Vars, "Full_12Vars");
                    bmc.PutVector(Full_012Vars, "Full_012Vars");
                    bmc.PutVector(Inner_0Vars, "Inner_0Vars");
                    bmc.PutVector(Inner_12Vars, "Inner_12Vars");
                    bmc.PutVector(Inner_012Vars, "Inner_012Vars");
                    bmc.PutVector(OneCell_0Vars, "OneCell_0Vars");
                    bmc.PutVector(OneCell_12Vars, "OneCell_12Vars");
                    bmc.PutVector(OneCell_012Vars, "OneCell_012Vars");
                    bmc.PutVector(Full_AllVars, "Full_AllVars");

                    //bmc.Cmd("retMatrix = [condNoFullMatrix,condNoDiffMatrix]");

                    bmc.Cmd("output = zeros(10,1)");
                    int k = 1;
                    foreach (var s in names) {
                        bmc.Cmd("output({1}) = condest(Mtx({0},{0}));", s, k);
                        k++;
                    }
                    bmc.GetMatrix(output, "output");

                    bmc.Execute(true);

                }
                Console.WriteLine("Finished Evaluation of Condition in Matlab");

                cond_Full_0Vars = output[0, 0];
                cond_Full_12Vars = output[1, 0];
                cond_Full_012Vars = output[2, 0];
                cond_Inner_0Vars = output[3, 0];
                cond_Inner_12Vars = output[4, 0];
                cond_Inner_012Vars = output[5, 0];
                cond_OneCell_0Vars = output[6, 0];
                cond_OneCell_12Vars = output[7, 0];
                cond_OneCell_012Vars = output[8, 0];
                cond_Full_AllVars = output[9, 0];

                //base.QueryHandler.ValueQuery("condNum", cond_Full_012Vars, true);

                Guid sessionID = this.CurrentSessionInfo.ID;
                TextWriter condfile = base.DatabaseDriver.FsDriver.GetNewLog("condition_numbers", sessionID);
                int jj = 0;
                foreach (var s in names) {
                    Console.WriteLine("Condition no {0}: {1:0.####e-00}", s, output[jj, 0]);
                    condfile.WriteLine("Condition no {0}: {1:0.####e-00}", s, output[jj, 0]);
                    condfile.Flush();
                    jj++;
                }

                //ExtractMatrixBlock(OneCell, new int[] { 0, 1, 2 }, PcdMatrix).ToFullMatrixOnProc0().SaveToTextFile(@"d:\Users\rieckmann\Helical\ErrorsAndMatrices\PcdMtxBlk.txt");

                //ExtractMatrixBlock(OneCell, new int[] { 0, 1, 2 }, RawMatrix).ToFullMatrixOnProc0().SaveToTextFile(@"d:\Users\rieckmann\Helical\ErrorsAndMatrices\RawMtxBlk.txt");


                #region Save operator matrix
                RawMatrix.SaveToTextFileSparse("matrixBoSSS.txt");
                PcdMatrix.SaveToTextFileSparse("PCDmatrixBoSSS.txt");
                #endregion
                Console.WriteLine("Operator Matrix written: matrixBoSSS.txt");



            }


        }
        */
        public double cond_Full_0Vars;
        public double cond_Full_12Vars;
        public double cond_Full_012Vars;
        public double cond_Inner_0Vars;
        public double cond_Inner_12Vars;
        public double cond_Inner_012Vars;
        public double cond_OneCell_0Vars;
        public double cond_OneCell_12Vars;
        public double cond_OneCell_012Vars;
        public double cond_Full_AllVars;

        /*
        private static void WriteInConsole(string whichMatrix, double condNumber) {
            Console.WriteLine(whichMatrix + condNumber);
        }

        */

        /*
        private MsrMatrix ExtractMatrixBlock(SubGrid InnerCells, int[] VarIdx, BlockMsrMatrix FullMatrix) {

            long[] innerCellIndices = this.UnknownsMap.GetSubvectorIndices(InnerCells, true, VarIdx);

            MsrMatrix InnerOpMatrix = new MsrMatrix(innerCellIndices.Length, innerCellIndices.Length, 1, 1);
            FullMatrix.WriteSubMatrixTo(InnerOpMatrix, innerCellIndices, default(long[]), innerCellIndices, default(long[]));

            return InnerOpMatrix;
        }


        private MsrMatrix ExtractMatrixBlock(int[] VarIdx, BlockMsrMatrix FullMatrix) {
            MsrMatrix DiffMatrix;

            long[] USubMatrixIdx_Row = this.UnknownsMap.GetSubvectorIndices(true, VarIdx);  // alle Zeilen in denen Geschwindigkeitsvariablen auftreten
            long[] USubMatrixIdx_Col = this.UnknownsMap.GetSubvectorIndices(true, VarIdx);  // alle Spalten in denen Geschwindigkeitsvariablen auftreten
            int L = USubMatrixIdx_Row.Length;

            DiffMatrix = new MsrMatrix(L, L, 1, 1);
            FullMatrix.WriteSubMatrixTo(DiffMatrix, USubMatrixIdx_Row, default(long[]), USubMatrixIdx_Col, default(long[]));

            return DiffMatrix;
        }
        */
    }
}

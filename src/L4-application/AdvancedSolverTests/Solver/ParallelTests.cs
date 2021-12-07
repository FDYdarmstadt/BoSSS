using AdvancedSolverTests.SolverChooser;
using BoSSS.Solution.AdvancedSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Solution;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using MPI.Wrappers;
using ilPSP.LinSolvers;

namespace AdvancedSolverTests.Solver {
    class ParallelTests  {


        public class BlockingStrat : Schwarz.BlockingStrategy {

            public BlockingStrat() {
            }

            private int IdentifyBlock(int iCell, MultigridOperator op) {
                //int NoOfNodes = cell.TransformationParams.NoOfRows;
                //int D = op.GridData.SpatialDimension;
                //var center = new double[D];
                //for (int k = 0; k < NoOfNodes; k++) {
                //    for (int d = 0; d < D; d++) {
                //        center[d] += cell.TransformationParams[k, d];
                //    }
                //}
                var center = op.Mapping.AggGrid.ParentGrid.iLogicalCells.GetCenter(iCell);
                if (center.x < 0 && center.y < 0)
                    return 0;
                if (center.x > 0 && center.y < 0)
                    return 1;
                if (center.x > 0 && center.y > 0)
                    return 2;
                if (center.x < 0 && center.y > 0)
                    return 3;
                return -1;
            }

            public override IEnumerable<List<int>> GetBlocking(MultigridOperator op) {
                if (op.Mapping.MpiSize > 1) {
                    //int MPIrank = op.Mapping.MpiRank;
                    //return new List<int>[] { SblockList[MPIrank] };
                    NoOfBlocks = 1;
                    return new List<int>[] { op.Mapping.LocalNoOfBlocks.ForLoop(i => i).ToList() };
                } else {
                    List<int>[] SblockList = NoOfBlocks.ForLoop(s => new List<int>());
                    long N_cells = op.Mapping.TotalNoOfBlocks;
                    for (int iCell = 0; iCell < N_cells; iCell++) {
                        int iBlock = IdentifyBlock(iCell, op);
                        SblockList[iBlock].Add(iCell);
                    }
                    return SblockList;
                }
            }

            private int NoOfBlocks = 4;

            public override int GetNoOfBlocks(MultigridOperator op) {
                return NoOfBlocks;
            }
        }

        public static bool RunTest() {

            System.Threading.Thread.Sleep(5000);

            var solver = new Schwarz() {
                //m_BlockingStrategy = ,
                m_BlockingStrategy = new BlockingStrat(),
                CoarseSolver = null,
                Overlap = 1,
                EnableOverlapScaling = false
            };

            var MgoPair = Utils.CreateTestMGOperator(XDGusage.all,3,MatrixShape.full_spec,8);
            solver.Init(MgoPair.MGOp);

            var op = MgoPair.MGOp;
            int NoOfBlocks = (solver as Schwarz).m_BlockingStrategy.GetNoOfBlocks(op);
            var SblockList = (solver as Schwarz).m_BlockingStrategy.GetBlocking(op).ToArray();
            var grid = op.Mapping.AggGrid;
            var basis = new Basis(op.BaseGridProblemMapping.GridDat, 0);
            var MPIranks = new SinglePhaseField(basis, "MPIrank");
            long L = grid.CellPartitioning.LocalLength;
            //for (int iCell = 0; iCell < L; iCell++) {
            //MPIranks.SetMeanValue(iCell, op.Mapping.MpiRank);
            //}
            for (int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {
                foreach (int iCell in SblockList[iBlock])
                    MPIranks.SetMeanValue(iCell, iBlock+op.Mapping.MpiRank);
            }
            var list = new List<SinglePhaseField>();
            list.Add(MPIranks);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(list, "BLargh.plt", 0.0, 0);
            int size = op.Mapping.MpiSize;
            long totLength = op.Mapping.TotalLength;
            double[] X, B;
            if (size > 1) {
                int locLength = op.Mapping.LocalLength;
                
                //X = new double[locLength];
                //B = new double[locLength];

                //B.LoadFromTextFile("B_dumbed");

                ////if (op.Mapping.MpiRank == 0) {
                ////X = totLength.ForLoop(s => new Random(0).NextDouble());
                ////B = totLength.ForLoop(s => new Random(1).NextDouble());
                ////}
                ////int[] loclenghs = size.ForLoop(iRank => op.Mapping.GetLocalLength(iRank));
                ////X = X.MPIScatterv(loclenghs);
                ////B = B.MPIScatterv(loclenghs);

                X = locLength.ForLoop(s => 1.0);
            } else {
                ////B = totLength.ForLoop(s => new Random().NextDouble());
                //X = new double[totLength];
                //////X = totLength.ForLoop(s => (double)s);
                //B = totLength.ForLoop(s => (double)s);
                //B.SaveToTextFile("B_dumbed");
                X = totLength.ForLoop(s => 1.0);
            }
            B = new double[X.Length];

            op.OperatorMatrix.Clear();
            op.OperatorMatrix.AccEyeSp();
            op.OperatorMatrix.SpMV(1.0, X, 0.0, B);
            X.Clear();
            Console.WriteLine("NormL2:" + B.MPI_L2NormPow2());
            solver.Solve(X, B);
            Console.WriteLine("NormL2:"+X.MPI_L2NormPow2());
            return true;
        }
    }
}

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using MPI.Wrappers;
using System;
using System.Diagnostics;

namespace IntersectingLevelSetTest
{

    /// <summary>
    /// some helper for manual debugging
    /// </summary>
    static class DebugTools
    {
        public static void MatrixTests(MsrMatrix OpMatrix, IGridData gridData, LevelSetTracker LsTrk)
        {
            int J = gridData.iLogicalCells.NoOfLocalUpdatedCells;
            int E = gridData.iLogicalEdges.Count;
            int[,] e2c = gridData.iLogicalEdges.CellIndices;

            MsrMatrix ConMatrix = new MsrMatrix(new Partitioning(J));
            MsrMatrix ConMatrix2 = new MsrMatrix(new Partitioning(J));
            var map = new UnsetteledCoordinateMapping(new Basis(gridData, 0));
            var FConMatrix = new XSpatialOperatorMk2.SpeciesFrameMatrix<MsrMatrix>(ConMatrix2, LsTrk.Regions, LsTrk.GetSpeciesId("B"), map, map);

            int jCell0 = gridData.CellPartitioning.i0;



            //for (int e = 0; e < E; e++) {
            foreach (int e in LsTrk.Regions.GetSpeciesSubGrid("B").InnerEdgesMask.ItemEnum)
            {
                int j0 = e2c[e, 0];
                int j1 = e2c[e, 1];

                int j0G;
                if (j0 >= J)
                {
                    j0G = (int)(gridData.iParallel.GlobalIndicesExternalCells[j0 - J]);
                }
                else
                {
                    j0G = j0 + jCell0;
                }

                int j1G;
                if (j1 >= 0)
                {
                    if (j1 >= J)
                    {
                        j1G = (int)(gridData.iParallel.GlobalIndicesExternalCells[j1 - J]);
                    }
                    else
                    {
                        j1G = j1 + jCell0;
                    }
                }
                else
                {
                    j1G = -1;
                }

                ConMatrix[j0G, j0G] += 1.0;
                FConMatrix[j0G, j0G] += 1.0;

                if (j1G >= 0)
                {
                    ConMatrix[j0G, j1G] += 1.0;
                    FConMatrix[j0G, j1G] += 1.0;

                    if (gridData.CellPartitioning.IsInLocalRange(j1G))
                    {
                        ConMatrix[j1G, j1G] += 1.0;
                        FConMatrix[j1G, j1G] += 1.0;
                        ConMatrix[j1G, j0G] += 1.0;
                        FConMatrix[j1G, j0G] += 1.0;
                    }
                }
            }


            for (int jCell = 0; jCell < J; jCell++)
            {
                Debug.Assert(jCell + jCell0 == gridData.iLogicalCells.GetGlobalID(jCell));
            }


            if (gridData.MpiSize == 1)
            {
                OpMatrix.SaveToFile("matrix.bin");
                ConMatrix.SaveToFile("conMtx.bin");
                ConMatrix2.SaveToFile("conMtx2.bin");
                OpMatrix.SaveToTextFile("matrix.txt");

            }
            else
            {
                var compOpMatrix = MsrMatrix.LoadFromFile("matrix.bin", csMPI.Raw._COMM.WORLD, OpMatrix.RowPartitioning, OpMatrix.ColPartition);
                var compConMatrix = MsrMatrix.LoadFromFile("conMtx.bin", csMPI.Raw._COMM.WORLD, ConMatrix.RowPartitioning, ConMatrix.ColPartition);
                var compConMatrix2 = MsrMatrix.LoadFromFile("conMtx2.bin", csMPI.Raw._COMM.WORLD, ConMatrix2.RowPartitioning, ConMatrix2.ColPartition);


                compOpMatrix.Acc(-1.0, OpMatrix);
                compConMatrix.Acc(-1.0, ConMatrix);
                compConMatrix2.Acc(-1.0, ConMatrix2);

                double OpErrNrm = compOpMatrix.InfNorm();
                double ConErrNrm = compConMatrix.InfNorm();
                double Con2ErrNrm = compConMatrix.InfNorm();


                Console.WriteLine("Matrix Comparison: {0}, {1}, {2}", ConErrNrm, Con2ErrNrm, OpErrNrm);

            }
        }
    }
}

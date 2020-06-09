using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace BoSSS.Foundation {
    public class TestingUtils {

        /// <summary>
        /// Enables the comparison of calculations with different number of MPI cores, 
        /// typically the comparison of a single-core vs. a parallel run.
        /// </summary>
        /// <param name="G"></param>
        /// <param name="FileName"></param>
        /// <param name="ReferenceMPISize">
        /// The number of MPI course for the reference computation, typically 1 (single core).
        /// - If the current MPI size is equal to this value, data is saved in file <paramref name="FileName"/>.
        /// - If the current MPI size is different, data is loaded from file <paramref name="FileName"/> and error norms are returned
        /// </param>
        /// <param name="infoFunc">
        /// Functions providing cell-wise data 
        /// </param>
        /// <returns></returns>
        static public (double[] AbsL2Errors, double[] RelL2Errors) Compare(IGridData G, string FileName, int ReferenceMPISize, params ExecutionMask.ItemInfo[] infoFunc) {
            if(ReferenceMPISize < 1)
                throw new ArgumentOutOfRangeException();

            var GatheredData = GatherAndSort(G, infoFunc);
            int D = G.SpatialDimension;

            string CoordName(int d) {
                return  "CellCenter" + (new[] { "X", "Y", "Z" })[d];
            }

            if(G.MpiSize == ReferenceMPISize) {

                if(G.MpiRank == 0) {
                    Dictionary<string, System.Collections.IEnumerable> table = new Dictionary<string, System.Collections.IEnumerable>();

                    table.Add("GlobalID", GatheredData.GlobalID);
                    for(int d = 0; d < D; d++) {
                        string name = CoordName(d);
                        table.Add(name, GatheredData.DataColumns[d]);
                    }

                    for(int iCol = D; iCol < GatheredData.DataColumns.Length; iCol++) {
                        table.Add("Col" + (iCol - G.SpatialDimension), GatheredData.DataColumns[iCol]);
                    }

                    table.SaveToCSVFile(FileName);
                }

                return (new double[infoFunc.Length],new double[infoFunc.Length]); 
            } else {

                double[] AbsErrors = null;
                double[] RelErrors = null;

                if(G.MpiRank == 0) {
                    Dictionary<string, System.Collections.IEnumerable> table = new Dictionary<string, System.Collections.IEnumerable>();

                    Func<string, object>[] parsers = new Func<string, object>[GatheredData.DataColumns.Length + 1];
                    parsers[0] = (s => long.Parse(s));
                    for(int i = 1; i < parsers.Length; i++)
                        parsers[i] = (s => double.Parse(s));

                    table.ReadFromCSVFile(filename:FileName, parsers:parsers);


                    AbsErrors = new double[infoFunc.Length];
                    RelErrors = new double[infoFunc.Length];

                    if(!GatheredData.GlobalID.ListEquals(table["GlobalID"] as IEnumerable<long>))
                        throw new IOException("Mismatch in GlobalID list - unable to compare.");

                    for(int d = 0; d < D; d++) {
                        IEnumerable<double> Coord_d_comp = table[CoordName(d)] as IEnumerable<double>;
                        IEnumerable<double> Coord_d = GatheredData.DataColumns[d];

                        double AbsErr = Coord_d.L2Distance(Coord_d_comp);
                        double MaxL2 = Math.Max(Coord_d.L2Norm(), Coord_d_comp.L2Norm());
                        double RelErr = AbsErr / MaxL2;
                        if(RelErr > BLAS.MachineEps.Sqrt())
                            throw new IOException("Mismatch for cell center coordinates - unable to compare.");

                    }

                    for(int iCol = D; iCol < infoFunc.Length; iCol++) {
                        IEnumerable<double> Coord_d_comp = table["Col" + iCol] as IEnumerable<double>;
                        IEnumerable<double> Coord_d = GatheredData.DataColumns[iCol];

                        double AbsErr = Coord_d.L2Distance(Coord_d_comp);
                        double MaxL2 = Math.Max(Coord_d.L2Norm(), Coord_d_comp.L2Norm());
                        double RelErr = AbsErr / MaxL2;

                        AbsErrors[iCol] = AbsErr;
                        RelErrors[iCol] = RelErr;
                    }

                }


                AbsErrors = AbsErrors.MPIBroadcast(0, G.CellPartitioning.MPI_Comm);
                RelErrors = RelErrors.MPIBroadcast(0, G.CellPartitioning.MPI_Comm);
                return (AbsErrors, RelErrors);
            }
        }



        /// <summary>
        /// Gathers cell-wise data on rank 0 in a sorted fashion (independent of grid permutation), which can be 
        /// compared among different runs, with different MPI sizes.
        /// </summary>
        static public (long[] GlobalID, double[][] DataColumns) GatherAndSort(IGridData G, params ExecutionMask.ItemInfo[] infoFunc) {

            // Build columns on each MPI process
            // =================================

            long[] GiD;
            double[][] Columns;
            {

                int D = G.SpatialDimension;
                int NoOfCol = infoFunc.Length + D;
                int J = G.iLogicalCells.NoOfLocalUpdatedCells;
                int[][] Agg2Part = G.iLogicalCells.AggregateCellToParts;

                GiD = G.CurrentGlobalIdPermutation.Values.CloneAs();
                Columns = NoOfCol.ForLoop(i => new double[J]);

                for(int j = 0; j < J; j++) {
                    int jGeom;
                    if(Agg2Part == null || Agg2Part[j] == null) {
                        jGeom = j;
                    } else {
                        jGeom = Agg2Part[j][0];
                    }

                    NodeSet localCenter = G.iGeomCells.GetRefElement(jGeom).Center;
                    var globalCenter = G.GlobalNodes.GetValue_Cell(localCenter, jGeom, 1);
                    double[] X = new double[D];


                    for(int d = 0; d < D; d++) {
                        Columns[d][j] = globalCenter[0, 0, d];
                        X[d] = globalCenter[0, 0, d];
                    }

                    for(int i = 0; i < NoOfCol - D; i++) {
                        Columns[i + D][j] = infoFunc[i](X, j, jGeom);
                    }
                }
            }

            // gather all data on MPI rank 0
            // =============================

            MPI_Comm comm = G.CellPartitioning.MPI_Comm;
            var AllGiD = GiD.MPIGatherO(0, comm);
            var AllCoumns = Columns.MPIGatherO(0, comm);
            int rank = G.MpiRank;
            int size = G.MpiSize;

            // write to stream on rank 0
            // =========================

            if(rank == 0) {

                // concat data
                // -----------

                long[] CatGiD = AllGiD[0];
                double[][] CatColumns = AllCoumns[0];

                int NoOfCols = Columns.Length;

                for(int r = 1; r < size; r++) {
                    CatGiD = CatGiD.Cat(AllGiD[r]);
                    for(int ic = 0; ic < NoOfCols; ic++) {
                        CatColumns[ic] = CatColumns[ic].Cat(AllCoumns[r][ic]);
                    }
                }
            

                // sort 
                // ----

                int[] SortIndices = CatGiD.Length.ForLoop(i => i);
                double[][] SortColumns= CatGiD.Length.ForLoop(i => new double[CatGiD.Length]);

                Array.Sort(CatGiD, SortIndices);

                for(int iRow = 0; iRow < CatGiD.Length; iRow++) {
                    Debug.Assert(iRow == CatGiD[iRow]);

                    int idx = SortIndices[iRow];

                    

                    for(int ic = 0; ic < NoOfCols; ic++) {
                        SortColumns[ic][iRow] = CatColumns[ic][idx];
                    }

                    
                }

                // return
                // ------

                return (CatGiD, SortColumns);

            } else {
                return default(ValueTuple<long[], double[][]>);
            }
        }
    }
}

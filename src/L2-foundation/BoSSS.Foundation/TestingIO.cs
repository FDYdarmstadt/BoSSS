using BoSSS.Foundation;
using BoSSS.Foundation.Comm;
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

    /// <summary>
    /// Utility class for testing: Enables the comparison of calculations with different number of MPI cores, 
    /// typically the comparison of a single-core vs. a parallel run, without using a BoSSS database
    /// </summary>
    public class TestingIO {
        /// <summary>
        /// Data file.
        /// </summary>
        public string CSVfile {
            get;
            private set;
        }

        /// <summary>
        /// The number of MPI course for the reference computation, typically 1 (single core).
        /// - (Reference mode): If the current MPI size is equal to this value, data from <see cref="CurrentData"/> is saved in file <paramref name="CSVfile"/> .
        /// - (Comparison mode): If the current MPI size is different, data is loaded from file <paramref name="CSVfile"/> into <see cref="ReferenceData"/>
        /// </summary>
        public int ReferenceMPISize {
            get;
            private set;
        }

        public IGridData GridDat {
            get;
            private set;
        }

        Dictionary<string, double[]> m_CurrentData = new Dictionary<string, double[]>();

        /// <summary>
        /// Data produced in current process 
        /// </summary>
        public IDictionary<string, double[]> CurrentData {
            get {
                return m_CurrentData;
            }
        }

        /// <summary>
        /// Reference data loaded from file 
        /// </summary>
        public IDictionary<string, double[]> ReferenceData {
            get;
            private set;
        }


        /// <summary>
        /// 
        /// </summary>
        public TestingIO(IGridData __g, string __CSVfile, int __ReferenceMPISize = 1) {
            if(__ReferenceMPISize < 1)
                throw new ArgumentOutOfRangeException();

            CSVfile = __CSVfile;
            GridDat = __g;
            this.ReferenceMPISize = __ReferenceMPISize;
            Setup();
        }

        const string ColName_Gid = "GlobalID";

        static string ColName_Coörd(int d) {
            return "CellCenter" + (new[] { "X", "Y", "Z" })[d];
        }

        void Setup() {
            var G = this.GridDat;

            int D = G.SpatialDimension;
            int J = G.iLogicalCells.NoOfLocalUpdatedCells;
            int[][] Agg2Part = G.iLogicalCells.AggregateCellToParts;

            double[] GiD = G.CurrentGlobalIdPermutation.Values.Select((long gid) => (double)gid).ToArray();
            Debug.Assert(GiD.Length == J);
            double[][] Coords = D.ForLoop(i => new double[J]);

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
                    Coords[d][j] = globalCenter[0, 0, d];
                    X[d] = globalCenter[0, 0, d];
                }
            }

            // ---

            m_CurrentData.Add(ColName_Gid, GiD);
            for(int d = 0; d < D; d++)
                m_CurrentData.Add(ColName_Coörd(d), Coords[d]);
        }

        /// <summary>
        /// Sorted column names
        /// </summary>
        public string[] ColumnNames {
            get {
                int D = GridDat.SpatialDimension;
                var colnames = this.m_CurrentData.Keys.ToList();
                colnames.Remove(ColName_Gid);
                for(int d = 0; d < D; d++)
                    colnames.Remove(ColName_Coörd(d));

                colnames.Sort();

                for(int d = D-1; d >=0 ; d--)
                    colnames.Insert(0, ColName_Coörd(d));
                colnames.Insert(0, ColName_Gid);

                return colnames.ToArray();
            }
        }

        /// <summary>
        /// Sorted column names
        /// </summary>
        public string[] ColumnNamesWithoutReserved {
            get {
                int D = GridDat.SpatialDimension;
                var colnames = this.m_CurrentData.Keys.ToList();
                colnames.Remove(ColName_Gid);
                for(int d = 0; d < D; d++)
                    colnames.Remove(ColName_Coörd(d));

                colnames.Sort();


                return colnames.ToArray();
            }
        }


        double[][] Dict2AoA<T>(IDictionary<string, T> table) where T : IEnumerable<double> {
            var coln = ColumnNames;
            if(!coln.SetEquals(table.Keys)) {
                throw new ApplicationException();
            }

            double[][] DataColumns = coln.Select(nmn => table[nmn].ToArray()).ToArray();

            return DataColumns;
        }

        Dictionary<string, double[]> AoA2Dict(double[][] AoA) {
            var coln = ColumnNames;
            if(coln.Length != AoA.Length) {
                throw new ApplicationException();
            }

            Dictionary<string, double[]> table = new Dictionary<string, double[]>();
            for(int iCol = 0; iCol < coln.Length; iCol++)
                table.Add(coln[iCol], AoA[iCol]);

            return table;
        }

        /// <summary>
        /// 
        /// </summary>
        public void DoIOnow() {
            if(GridDat.CellPartitioning.MpiSize == ReferenceMPISize) {
                SaveData();
                this.ReferenceData = this.CurrentData;
            } else {
                LoadData();
            }
        }



        /// <summary>
        /// saves data in <see cref="CSVfile"/>.
        /// </summary>
        void SaveData() {

            
            double[][] DataColumns = Dict2AoA(m_CurrentData);

            var GatheredData = GatherAndSort(this.GridDat, DataColumns);

            if(GridDat.CellPartitioning.MpiRank == 0) {
                var table = AoA2Dict(GatheredData);

                table.SaveToCSVFile(CSVfile);

            }
        }
        /// <summary>
        /// Loads data from <see cref="CSVfile"/> (on rank 0) and scatters it among MPI processors
        /// </summary>
        void LoadData() {

      
            double[][] DataColumns;

            if(GridDat.CellPartitioning.MpiRank == 0) {
                Dictionary<string, IEnumerable<double>> table = new Dictionary<string, IEnumerable<double>>();
                table.ReadFromCSVFile(filename: CSVfile);

                double[][] Rank0Columns = Dict2AoA(table);
                DataColumns = ScatterSortedData(this.GridDat, Rank0Columns);
            } else {
                DataColumns = ScatterSortedData(this.GridDat, null);
            }

            ReferenceData = AoA2Dict(DataColumns);

            if(AbsError(ColName_Gid) > 0.0)
                throw new ArgumentException("Mismatch in Global ID between reference data and current grid.");
            for(int d = 0; d < GridDat.SpatialDimension; d++) {
                if(RelError(ColName_Coörd(d)) > BLAS.MachineEps.Sqrt())
                    throw new ArgumentException("Mismatch in cell center coordinates.");
            }
        }

        /// <summary>
        /// Absolute L2 for DG fields
        /// </summary>
        public double AbsError(ConventionalDGField f) {
            if(GridDat.MpiSize == ReferenceMPISize)
                return 0.0;

            var err = f.CloneAs();
            OverwriteDGField(err);
            err.Acc(-1.0, f);

            return err.L2Norm();
        }

        /// <summary>
        /// Absolute L2 error between columns
        /// </summary>
        public double AbsError(string Colname) {
            if(GridDat.MpiSize == ReferenceMPISize)
                return 0.0;

            double[] CurGid = m_CurrentData[Colname];
            double[] RefGid = ReferenceData[Colname];

            return CurGid.L2Distance(RefGid).Pow2().MPISum(GridDat.CellPartitioning.MPI_Comm).Sqrt();
        }

        /// <summary>
        /// Difference between data in file and data in <paramref name="vec"/>
        /// </summary>
        public double[] LocError(string Colname) {
            
            double[] RefData = ReferenceData[Colname];
            double[] LoclErr = m_CurrentData[Colname].CloneAs();

            LoclErr.AccV(-1.0, RefData);

            return LoclErr;
        }



        /// <summary>
        /// Relative L2 error between columns
        /// </summary>
        public double RelError(string Colname) {
            if(GridDat.MpiSize == ReferenceMPISize)
                return 0.0;

            double[] CurGid = m_CurrentData[Colname];
            double[] RefGid = ReferenceData[Colname];

            double[] vals = new[] {
                CurGid.L2Distance(RefGid).Pow2(),
                CurGid.L2NormPow2(),
                RefGid.L2NormPow2()
                };

            vals = vals.MPISum(GridDat.CellPartitioning.MPI_Comm);

            return vals[0].Sqrt() / Math.Max(Math.Max(vals[1], vals[2]), double.Epsilon*1.0e20).Sqrt();

        }

        /// <summary>
        /// Relative Errors for all known columns
        /// </summary>
        public IDictionary<string, double> AllRelErr() {
            var R = new Dictionary<string, double>();
            foreach(string col in this.ColumnNames) {
                R.Add(col, RelError(col));
            }
            return R;
        }

        /// <summary>
        /// Relative Errors for all known columns
        /// </summary>
        public IDictionary<string, double> AllAbsErr() {
            var R = new Dictionary<string, double>();
            foreach(string col in this.ColumnNames) {
                R.Add(col, AbsError(col));
            }
            return R;
        }
        
        /*
        /// <summary>
        /// Comparison of reference data 
        /// </summary>
        public (double[] AbsL2Errors, double[] RelL2Errors) Compare() {
            if(GridDat.MpiSize == ReferenceMPISize) {
                int N = m_CurrentData.Count;
                return (new double[N], new double[N]);

            } else { 


            
                int D = this.GridDat.SpatialDimension;


                double[] AbsErrors = null;
                double[] RelErrors = null;

                double[][] ReferenceData;

                if(G.MpiRank == 0) {
                    




                    AbsErrors = new double[infoFunc.Length];
                    RelErrors = new double[infoFunc.Length];

                    if(!GatheredData.GlobalID.ListEquals(table["GlobalID"] as IEnumerable<long>))
                        throw new IOException("Mismatch in GlobalID list - unable to compare.");

                    for(int d = 0; d < D; d++) {
                        IEnumerable<double> Coord_d_ref = table[ColName_Coörd(d)] as IEnumerable<double>;
                        IEnumerable<double> Coord_d = GatheredData.DataColumns[d];

                        double AbsErr = Coord_d.L2Distance(Coord_d_ref);
                        double MaxL2 = Math.Max(Coord_d.L2Norm(), Coord_d_ref.L2Norm());
                        double RelErr = AbsErr / MaxL2;
                        if(RelErr > BLAS.MachineEps.Sqrt())
                            throw new IOException("Mismatch for cell center coordinates - unable to compare.");

                    }

                    double[][] Rank0Data = new double[infoFunc.Length][];
                    if(ComparisonFile != null) {
                        Dictionary<string, System.Collections.IEnumerable> CompTable = new Dictionary<string, System.Collections.IEnumerable>();
                        CompTable.Add("GlobalID", GatheredData.GlobalID);
                        for(int d = 0; d < D; d++) {
                            CompTable.Add(ColName_Coörd(d), table[ColName_Coörd(d)]);
                        }

                        for(int iCol = D; iCol < infoFunc.Length + D; iCol++) {
                            string colname = "Col" + (iCol - D);
                            IEnumerable<double> Ref = table[colname] as IEnumerable<double>;
                            IEnumerable<double> Cur = GatheredData.DataColumns[iCol];
                            double[] Err = Cur.ToArray().CloneAs();
                            Err.AccV(-1.0, Ref.ToArray());

                            CompTable.Add(colname + "(REF)", Ref);
                            CompTable.Add(colname + "(CUR)", Cur);
                            CompTable.Add(colname + "(ERR)", Err);

                            Rank0Data[iCol] = Ref.ToArray();
                        }

                        CompTable.SaveToCSVFile(ComparisonFile);
                    }

                    Debug.Assert(infoFunc.Length + D == GatheredData.DataColumns.Length);

                    for(int iCol = D; iCol < infoFunc.Length + D; iCol++) {
                        IEnumerable<double> Coord_d_comp = table["Col" + (iCol - D)] as IEnumerable<double>;
                        IEnumerable<double> Coord_d = GatheredData.DataColumns[iCol];

                        double AbsErr = Coord_d.L2Distance(Coord_d_comp);
                        double MaxL2 = Math.Max(Coord_d.L2Norm(), Coord_d_comp.L2Norm());
                        double RelErr = AbsErr / MaxL2;

                        AbsErrors[iCol - D] = AbsErr;
                        RelErrors[iCol - D] = RelErr;
                    }

                    ReferenceData = ScatterSortedData(G, Rank0Data);

                } else {
                    ReferenceData = ScatterSortedData(G, null);
                }


                AbsErrors = AbsErrors.MPIBroadcast(0, G.CellPartitioning.MPI_Comm);
                RelErrors = RelErrors.MPIBroadcast(0, G.CellPartitioning.MPI_Comm);
                return (AbsErrors, RelErrors, ReferenceData);
            }
        }
        */


        /// <summary>
        /// 
        /// </summary>
        public void AddVector(string ColName, IEnumerable<double> data) {
            int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            int L = data.Count();
            int K = L / J;

            
            if(J*K != L || K < 1)
                throw new ArgumentException("wrong length of input vector");
            if(L == 1) {
                m_CurrentData.Add(ColName, data.ToArray());
            } else {
                var M = MultidimensionalArray.CreateWrapper(data.ToArray(), J, K);
                for(int k = 0; k < K; k++) {
                    m_CurrentData.Add(ColName + "-k" + k, M.GetColumn(k));
                }

            }
        }


        /// <summary>
        /// Collects cell-wise data from <see cref="ExecutionMask.ItemInfo"/>
        /// </summary>
        public void AddColumn(string ColName, ExecutionMask.ItemInfo infoFunc) {

            long[] GiD;
            double[] Column;
            var G = this.GridDat;
            int D = G.SpatialDimension;

            int J = G.iLogicalCells.NoOfLocalUpdatedCells;
            int[][] Agg2Part = G.iLogicalCells.AggregateCellToParts;

            GiD = G.CurrentGlobalIdPermutation.Values.CloneAs();
            Column = new double[J];

            double[][] Coords = D.ForLoop(d => this.m_CurrentData[ColName_Coörd(d)]);


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

                    X[d] = Coords[d][j];
                }


                Column[j] = infoFunc(X, j, jGeom);

            }


            m_CurrentData.Add(ColName, Column);
        }

        /// <summary>
        /// Adds a DG field
        /// </summary>
        public void AddDGField(ConventionalDGField f) {
            if(!object.ReferenceEquals(GridDat, f.GridDat))
                throw new ArgumentException("DG field is defined on different mesh");

            int N = f.Basis.Length;
            int J = GridDat.CellPartitioning.LocalLength;

            for(int n = 0; n < N; n++) {
                var col = f.Coordinates.GetColumn(n).GetSubVector(0, J);
                m_CurrentData.Add(ColName_DGfield(f, n), col);
            }
        }

        /// <summary>
        /// Adds a DG field
        /// </summary>
        public ConventionalDGField LocalError(ConventionalDGField f) {
            if(!object.ReferenceEquals(GridDat, f.GridDat))
                throw new ArgumentException("DG field is defined on different mesh");

            var Error = f.CloneAs();
            Error.Clear();
            OverwriteDGField(Error);
            Error.Scale(-1);
            Error.Acc(1.0, f);

            Error.Identification = "Error-" + f.Identification;
            return Error;
        }




        static string ColName_DGfield(DGField f, int n) {
            return f.Identification + "_mode" + n;
        }

        /// <summary>
        /// Overwrites the memory of a DG field with the reference data 
        /// </summary>
        /// <param name="f"></param>
        public void OverwriteDGField(ConventionalDGField f) {
            if(!object.ReferenceEquals(GridDat, f.GridDat))
                throw new ArgumentException("DG field is defined on different mesh");

            int N = f.Basis.Length;
            int N1 = ReferenceData.Keys.Where(colName => colName.StartsWith(f.Identification)).Count();
            if(N != N1)
                throw new ArgumentException("DG degree seems different");

            int J = this.GridDat.CellPartitioning.LocalLength;

            for(int n = 0; n < N; n++) {
                var col = ReferenceData[ColName_DGfield(f, n)];
                //for(int j = 0; j < J; j++)
                //    f.Coordinates[j, n] = col[j];
                f.Coordinates.SetColumn(n, col);
            }
        }

        /// <summary>
        /// Gathers cell-wise data on rank 0 in a sorted fashion (independent of grid permutation), which can be 
        /// compared among different runs, with different MPI sizes.
        /// </summary>
        static public double[][] GatherAndSort(IGridData G, double[][] Columns) {

            // Build columns on each MPI process
            // =================================

            long[] GiD = G.CurrentGlobalIdPermutation.Values.CloneAs();
         

            // gather all data on MPI rank 0
            // =============================

            MPI_Comm comm = G.CellPartitioning.MPI_Comm;
            var AllGiD = GiD.MPIGatherO(0, comm);
            var AllCoumns = Columns.MPIGatherO(0, comm);
            int rank = G.MpiRank;
            int size = G.MpiSize;


            if(rank == 0) {
                // +++++++++++++++++++++
                // concat data on rank 0
                // +++++++++++++++++++++

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

                Debug.Assert(CatColumns.Length == NoOfCols);

                // sort 
                // ----

                int[] SortIndices = CatGiD.Length.ForLoop(i => i);
                double[][] SortColumns = NoOfCols.ForLoop(i => new double[CatGiD.Length]);

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

                return SortColumns;

            } else {
                return null;
            }
        }

        /// <summary>
        /// Scatters data from rank 0 (independent of grid permutation) to all processes, according to MPI partition.
        /// </summary>
        /// <param name="SortedDataColumns">Data, sorted according to GlobalId, at rank 0</param>
        /// <param name="G"></param>
        static public double[][] ScatterSortedData(IGridData G, double[][] SortedDataColumns) {
            //Debugger.Launch();
            int Jglob = G.CellPartitioning.TotalLength;
            int Jloc = G.CellPartitioning.LocalLength;
            int MPIrank = G.MpiRank;
            MPI_Comm comm = G.CellPartitioning.MPI_Comm;
            int NoOfColumns = MPIrank == 0 ? SortedDataColumns.Length : 0;
            NoOfColumns = NoOfColumns.MPIBroadcast(0);


            if(MPIrank == 0) {
                foreach(double[] vec in SortedDataColumns)
                    if(vec.Length != Jglob)
                        throw new ArgumentException();
            } else {
                SortedDataColumns = new double[NoOfColumns][];
                for(int iCol = 0; iCol < NoOfColumns; iCol++)
                    SortedDataColumns[iCol] = new double[0];
            }


            // Compute resorting permutation
            // =============================

            Permutation Resorting;
            {
                // id    is the GlobalID-permutation that we have for the loaded vector
                // sigma is the current GlobalID-permutation of the grid
                Permutation sigma = G.CurrentGlobalIdPermutation;
                Permutation id;// new Permutation(DataVec.Select(cd => cd.GlobalID).ToArray(), csMPI.Raw._COMM.WORLD);
                if(MPIrank == 0) {
                    id = new Permutation(Jglob.ForLoop((int i) => (long)i), comm);
                } else {
                    id = new Permutation(new long[0], comm);
                }

                // compute resorting permutation
                Permutation invSigma = sigma.Invert();
                Resorting = invSigma * id;
                id = null;
                invSigma = null;
            }

            // test
            {
                long[] Check;
                if(MPIrank == 0)
                    Check = Jglob.ForLoop((int i) => (long)i);
                else
                    Check = new long[0];

                long[] Resort = new long[Jloc];
                Resorting.ApplyToVector(Check, Resort, G.CellPartitioning);

                for(int j = 0; j < Jloc; j++) {
                    if(Resort[j] != G.CurrentGlobalIdPermutation.Values[j])
                        throw new ApplicationException("error in algorithm");
                }
            }


            // apply resorting
            // ===============

            double[][] RetColunms = new double[NoOfColumns][];
            for(int iCol = 0; iCol < NoOfColumns; iCol++) {
                RetColunms[iCol] = new double[Jloc];

                Resorting.ApplyToVector(SortedDataColumns[iCol], RetColunms[iCol], G.CellPartitioning);

            }

            return RetColunms;
        }


    }
}

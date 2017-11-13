using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Collections;


namespace BoSSS.Application.AdaptiveMeshRefinementTest {

    /// <summary>
    /// App which performs basic tests on adaptive mesh refinement and dynamic load balancing.
    /// </summary>
    class AdaptiveMeshRefinementTestMain : BoSSS.Solution.Application {

        static void RefineScheisse() {
            double[] nodes = GenericBlas.Linspace(-5, 5, 7);

            BoundingBox BL = new BoundingBox(new[] { new[] { -10.0, -10.0 }, new[] { 0.0, 0.0 } });
            GridCommons grd = Grid2D.Cartesian2DGrid(nodes, nodes, CutOuts: new BoundingBox[] { BL });
            GridData gdat = new GridData(grd);

            SinglePhaseField u = new SinglePhaseField(new Basis(gdat, 0), "u");
            Tecplot.PlotFields(new DGField[] { u }, "Refinment-0", 0.0, 0);

            int MaxRefinement = 5;

            for(int iLevel = 0; iLevel <= MaxRefinement; iLevel++) {
                Console.WriteLine("Level " + iLevel);

                HashSet<int> CellsToRefine = new HashSet<int>();
                int J = gdat.Cells.NoOfLocalUpdatedCells;
                for(int j = 0; j < J; j++) {
                    double xc = gdat.Cells.CellCenter[j, 0];
                    double yc = gdat.Cells.CellCenter[j, 1];
                    int RefinementLevel = gdat.Cells.GetCell(j).RefinementLevel;

                    double r = Math.Sqrt(xc * xc + yc * yc);

                    if((RefinementLevel < MaxRefinement) && (r < gdat.Cells.h_minGlobal*2)) {
                        CellsToRefine.Add(j);
                        RefineNeighboursRecursive(gdat, CellsToRefine, j, RefinementLevel);

                    }
                }

                Console.Write("  refine: ");
                foreach(int j in CellsToRefine) {
                    Console.Write(j + " ");
                }
                Console.WriteLine();


                grd = gdat.Adapt(CellsToRefine.ToArray(), null, out GridCorrelation Old2NewCorr);
                gdat = new GridData(grd);
                u = new SinglePhaseField(new Basis(gdat, 0), "u");
                Tecplot.PlotFields(new DGField[] { u }, "Refinment-" + iLevel, 0.0, 0);

            }
        }

        static void RefineNeighboursRecursive(GridData gdat, HashSet<int> CellsToRefine, int j, int DesiredLevel) {

            foreach(var jNeigh in gdat.Cells.CellNeighbours[j]) {
                var cl = gdat.Cells.GetCell(j);
                if(cl.RefinementLevel < DesiredLevel) {
                    CellsToRefine.Add(jNeigh);
                    RefineNeighboursRecursive(gdat, CellsToRefine, jNeigh, cl.RefinementLevel);
                }
            }

        }



        static void Main(string[] args) {
            //bool MpiInit;
            //ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out MpiInit);

            //RefineScheisse();

            //return;

            BoSSS.Solution.Application._Main(
                args,
                true,
                null,
                () => new AdaptiveMeshRefinementTestMain());
        }

        protected override GridCommons CreateOrLoadGrid() {
            double[] nodes = GenericBlas.Linspace(-5, 5, 21);
            var grd = Grid2D.Cartesian2DGrid(nodes, nodes);
            base.m_GridPartitioningType = GridPartType.none;
            return grd;
        }

        SinglePhaseField u;
        VectorField<SinglePhaseField> Grad_u;
        SinglePhaseField MagGrad_u;

        protected override void CreateFields() {
            var Basis = new Basis(base.GridData, DEGREE);
            u = new SinglePhaseField(Basis, "u");
            Grad_u = new VectorField<SinglePhaseField>(base.GridData.SpatialDimension, Basis, "Grad_u", (bs, nmn) => new SinglePhaseField(bs, nmn));
            MagGrad_u = new SinglePhaseField(new Basis(base.GridData, 0), "Magnitude_Grad_u");
            base.m_RegisteredFields.Add(u);
            base.m_RegisteredFields.AddRange(Grad_u);
            base.m_RegisteredFields.Add(MagGrad_u);
        }


        /// <summary>
        /// DG polynomial degree
        /// </summary>
        internal int DEGREE = 3;

        /// <summary>
        /// Setting initial value.
        /// </summary>
        protected override void SetInitial() {
            UpdateBaseGrid(0.0);

            RefinedGrid = this.GridData;
            Refined_u = this.u;
            Refined_TestData = new SinglePhaseField(this.u.Basis, "TestData");
            Refined_Grad_u = this.Grad_u;
            Refined_MagGrad_u = this.MagGrad_u;

            Refined_TestData.ProjectField((x, y) => x * x + y * y);
        }

        private void UpdateBaseGrid(double time) {

            u.Clear();
            u.ProjectField(NonVectorizedScalarFunction.Vectorize(uEx, time));
            
            Grad_u.Clear();
            Grad_u.Gradient(1.0, u);

            MagGrad_u.Clear();
            MagGrad_u.ProjectFunction(1.0,
                (double[] X, double[] U, int jCell) => Math.Sqrt(U[0].Pow2() + U[1].Pow2()),
                new Foundation.Quadrature.CellQuadratureScheme(),
                Grad_u.ToArray());
        }

        static double uEx(double[] X, double time) {
            double x = X[0];
            double y = X[1];
            
            double r = Math.Sqrt(x * x + y * y);
            double alpha = Math.Atan2(y, x);

            double alpha0 = alpha - time;

            double x0 = Math.Cos(alpha0) * r;
            double y0 = Math.Sin(alpha0) * r;
            
            double val = Math.Exp(- x0.Pow2() - (y0 - 2.5).Pow2()); 

            return val;
        }


        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
            
            if (L == null) {
                
            } else {
                throw new NotImplementedException("todo");
            }
        }

        GridData RefinedGrid = null;
         
        SinglePhaseField Refined_u;
        VectorField<SinglePhaseField> Refined_Grad_u;
        SinglePhaseField Refined_MagGrad_u;
        SinglePhaseField Refined_TestData;


        int[][] FindCoarseningClusters(BitArray Ok2Coarsen) {
            int JE = RefinedGrid.Cells.NoOfCells;
            int J = RefinedGrid.Cells.NoOfLocalUpdatedCells;

            if (Ok2Coarsen.Length != JE)
                throw new ArgumentException();

            int[][] CellNeighbours = RefinedGrid.Cells.CellNeighbours;

            List<Cell> temp = new List<Cell>();
            List<int> tempCC = new List<int>();

            int RecursionDepht = -1;
            if (RefinedGrid.SpatialDimension == 1)
                RecursionDepht = 1;
            else if (RefinedGrid.SpatialDimension == 2)
                RecursionDepht = 2;
            else if (RefinedGrid.SpatialDimension == 3)
                RecursionDepht = 4;
            else
                throw new NotSupportedException();

            int[][] CoarseningCluster = new int[JE][];


            BitArray marker = new BitArray(JE);
            for(int j = 0; j < J; j++) {
                if (marker[j])
                    continue;
                if (!Ok2Coarsen[j])
                    continue;

                Cell Cell_j = RefinedGrid.Cells.GetCell(j);
                int Level = Cell_j.RefinementLevel;

                if (Level == 0) {
                    marker[j] = true;
                    continue;
                }

                temp.Clear();
                temp.Add(Cell_j);

                tempCC.Clear();
                tempCC.Add(j);

                //SearchGids = new long[Cell_j.CoarseningPeers.Length + 1];
                //Array.Copy(Cell_j.CoarseningPeers, SearchGids, SearchGids.Length - 1);
                //SearchGids[SearchGids.Length - 1] = Cell_j.GlobalID;
                int SearchID = Cell_j.CoarseningClusterID;
                int ClusterSize = Cell_j.CoarseningClusterSize;
                Debug.Assert(SearchID > 0);


                bool complete = false;
                FindCoarseiningClusterRecursive(j, RecursionDepht, CellNeighbours, marker, RefinedGrid.Cells.GetCell, Level, SearchID, ClusterSize, tempCC, ref complete);
                foreach(int jC in tempCC) {
                    marker[jC] = true;
                }

                if(complete) {
                    if (!tempCC.Any(jC => Ok2Coarsen[jC] == false)) {
                        int[] CC = tempCC.ToArray();
                        foreach (int jC in CC) {
                            Debug.Assert(CoarseningCluster[jC] == null);
                            CoarseningCluster[jC] = CC;
                        }
                    } else {
                        //Console.WriteLine("Not ok to coarsen.");
                    }
                }
            }


            return CoarseningCluster;
        }

        static void FindCoarseiningClusterRecursive(int j, int MaxRecursionDeph,
            int[][] CellNeighbours, BitArray marker, Func<int, Cell> GetCell, int Level, int SearchID, int ClusterSize,
            List<int> CoarseiningCluster, ref bool complete) {

            Debug.Assert(CoarseiningCluster.Contains(j) == true);

            int[] Neighs = CellNeighbours[j];
            foreach (var jNeigh in Neighs) {
                if (marker[jNeigh] == true)
                    continue;
                if (CoarseiningCluster.Contains(jNeigh))
                    continue;

                Cell Cell_neigh = GetCell(jNeigh);
                if (Cell_neigh.RefinementLevel != Level)
                    continue;

                if (Cell_neigh.CoarseningClusterID != SearchID)
                    continue;

                CoarseiningCluster.Add(jNeigh);
                
                if(CoarseiningCluster.Count == ClusterSize) {
                    complete = true;
#if DEBUG
                    foreach(int j1 in CoarseiningCluster) {
                        Cell Cell_j1 = GetCell(j1);
                        Debug.Assert(Cell_j1.RefinementLevel > 0);
                        //Debug.Assert(Cell_j1.CoarseningPeers != null);
                        //Debug.Assert(Cell_j1.CoarseningPeers.Length == CoarseiningCluster.Count - 1);
                        Debug.Assert(Cell_j1.RefinementLevel == Level);
                        Debug.Assert(Cell_j1.CoarseningClusterID == SearchID);
                        Debug.Assert(Cell_j1.CoarseningClusterSize == ClusterSize);
                        //Debug.Assert(Array.IndexOf(SearchGids, Cell_j1.GlobalID) >= 0);
                        //foreach(long gid_cp in Cell_j1.CoarseningPeers) {
                        //    Debug.Assert(Array.IndexOf(SearchGids, gid_cp) >= 0);
                        //}
                    }
#endif
                    return;
                }

                if (MaxRecursionDeph > 0)
                    FindCoarseiningClusterRecursive(jNeigh, MaxRecursionDeph - 1,
                        CellNeighbours, marker, GetCell, Level, SearchID, ClusterSize, CoarseiningCluster, ref complete);

                if(complete) {
                    return;
                }
            }
            
        }



        void UpdateRefinedGrid(double time, int iTimestep) {
            int oldJ = RefinedGrid.Cells.NoOfLocalUpdatedCells;

            double[][] TestData_DGCoordinates = new double[oldJ][];
            for(int j = 0; j < oldJ; j++) {
                TestData_DGCoordinates[j] = Refined_TestData.Coordinates.GetRow(j);
            }


            bool NoRefinement = true;
            int[] DesiredLevel = new int[oldJ];
            for(int j = 0; j < oldJ; j++) {
                double GradMag = Refined_MagGrad_u.GetMeanValue(j);

                int DesiredLevel_j = 0;
                if(GradMag > 0.6)
                    DesiredLevel_j = 1;
                if(GradMag > 0.81)
                    DesiredLevel_j = 2;
                
                if(DesiredLevel[j] < DesiredLevel_j) {
                    DesiredLevel[j] = DesiredLevel_j;
                    NoRefinement = false;
                    RefineNeighboursRecursive(RefinedGrid, DesiredLevel, j, DesiredLevel_j - 1);
                }
            }

            SinglePhaseField Ok2CoarsenViz = new SinglePhaseField(new Basis(RefinedGrid, 0), "Ok2Coarsen");
            SinglePhaseField CombinedMarkr = new SinglePhaseField(new Basis(RefinedGrid, 0), "CombinedMarker");
            SinglePhaseField CClustersMViz = new SinglePhaseField(new Basis(RefinedGrid, 0), "CoarseningClusters");
            BitArray Ok2Coarsen = new BitArray(oldJ);
            var rnd = new Random();
            for(int j = 0; j < oldJ; j++) {
                int ActualLevel_j = RefinedGrid.Cells.GetCell(j).RefinementLevel;
                int DesiredLevel_j = DesiredLevel[j];

                if(ActualLevel_j > DesiredLevel_j) {
                    Ok2CoarsenViz.SetMeanValue(j, 1.0);
                    Ok2Coarsen[j] = true;
                }
                //Ok2Coarsen[j] = true;
            }

            int[][] CClusters = FindCoarseningClusters(Ok2Coarsen);

            List<int[]> Coarsening = new List<int[]>();
            int NoOfCellsToCoarsen = 0;
            for(int j = 0; j < oldJ; j++) {
                               
                if(CClusters[j] != null) {
                    NoOfCellsToCoarsen++;
                    double jitter = rnd.NextDouble() * 0.5 + 1.0;
                    foreach(int jc in CClusters[j]) {
                        if (CClustersMViz.GetMeanValue(jc) == 0.0) {
                            CClustersMViz.SetMeanValue(jc, jitter);
                        }
                    }

                    Debug.Assert(CClusters[j].Contains(j));
                    if(j == CClusters[j].Min()) {
                        Coarsening.Add(CClusters[j]);
                    }
                }
            }
            //if (iTimestep == 10) {
            //    Console.WriteLine("Now really coarsening");
            //} else {
            //    Coarsening.Clear();
            //}

            CombinedMarkr.ProjectProduct(1.0, Ok2CoarsenViz, CClustersMViz);

            Console.WriteLine("       No of cells to coarsen: " + NoOfCellsToCoarsen);

            Tecplot.PlotFields(ArrayTools.Cat<DGField>(RefinedGrid.BoundaryMark(), Ok2CoarsenViz, CClustersMViz, CombinedMarkr ), "Coarsen", time, 0);

            GridCorrelation Old2NewCorr;
            if((!NoRefinement) || (Coarsening.Count > 0)) {


                List<int> CellsToRefineList = new List<int>();
                for(int j = 0; j < oldJ; j++) {
                    int ActualLevel_j = RefinedGrid.Cells.GetCell(j).RefinementLevel;
                    int DesiredLevel_j = DesiredLevel[j];

                    if(ActualLevel_j < DesiredLevel_j)
                        CellsToRefineList.Add(j);
                }

                Console.WriteLine("       Refining " + CellsToRefineList.Count + " of " + oldJ + " cells");

                // create new Grid & DG fields
                // ===========================

                GridCommons newGrid = RefinedGrid.Adapt(CellsToRefineList, Coarsening, out Old2NewCorr);
                RefinedGrid = new GridData(newGrid);

                var Basis = new Basis(RefinedGrid, DEGREE);
                Refined_u = new SinglePhaseField(Basis, "u");
                Refined_TestData = new SinglePhaseField(Basis, "TestData");
                Refined_Grad_u = new VectorField<SinglePhaseField>(base.GridData.SpatialDimension, Basis, "Grad_u", (bs, nmn) => new SinglePhaseField(bs, nmn));
                Refined_MagGrad_u = new SinglePhaseField(new Basis(RefinedGrid, 0), "Magnitude_Grad_u");
            } else {
                Old2NewCorr = null;
            }


            {
                Refined_u.Clear();
                Refined_u.ProjectField(NonVectorizedScalarFunction.Vectorize(uEx, time));

                Refined_Grad_u.Clear();
                Refined_Grad_u.Gradient(1.0, Refined_u);

                Refined_MagGrad_u.Clear();
                Refined_MagGrad_u.ProjectFunction(1.0,
                    (double[] X, double[] U, int jCell) => Math.Sqrt(U[0].Pow2() + U[1].Pow2()),
                    new Foundation.Quadrature.CellQuadratureScheme(),
                    Refined_Grad_u.ToArray());


                if(Old2NewCorr != null) {
                    Old2NewCorr.ComputeDataRedist(RefinedGrid);

                    int newJ = RefinedGrid.Cells.NoOfLocalUpdatedCells;

                    double[][][] ReDistDGCoords = new double[newJ][][];
                    Old2NewCorr.ApplyToVector(TestData_DGCoordinates, ReDistDGCoords, RefinedGrid.CellPartitioning);

                    for(int j = 0; j < newJ; j++) {
                        dfdsf
                    }

                }

            }
        }

        static void RefineNeighboursRecursive(GridData gdat, int[] DesiredLevel, int j, int DesiredLevelNeigh) {

            if(DesiredLevelNeigh <= 0)
                return;
            
            foreach(var jNeigh in gdat.Cells.CellNeighbours[j]) {
                var cl = gdat.Cells.GetCell(j);
                if(cl.RefinementLevel < DesiredLevelNeigh) {
                    DesiredLevel[jNeigh] = DesiredLevelNeigh;
                    RefineNeighboursRecursive(gdat, DesiredLevel, jNeigh, DesiredLevelNeigh - 1);
                }
            }
        }



        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                // Elo
                Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime + " ... ");

                // timestepping params
                base.NoOfTimesteps = 100;
                dt = Math.PI*2 / base.NoOfTimesteps;

                UpdateBaseGrid(phystime + dt);
                UpdateRefinedGrid(phystime + dt, TimestepNo);

                // return
                //base.TerminationKey = true;
                return dt;
            }
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "BaseGrid." + timestepNo;
            Tecplot.PlotFields(base.m_RegisteredFields.ToArray(), filename, physTime, superSampling);

            DGField[] RefinedFields = new[] { Refined_u, Refined_TestData, Refined_Grad_u[0], Refined_Grad_u[1], Refined_MagGrad_u };
            string filename2 = "RefinedGrid." + timestepNo;
            Tecplot.PlotFields(RefinedFields, filename2, physTime, superSampling);
        }
    }
}

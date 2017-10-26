using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution;
using BoSSS.Solution.Multigrid;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

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


                grd = gdat.Adapt(CellsToRefine.ToArray(), null);
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
            bool MpiInit;
            ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out MpiInit);

            RefineScheisse();

            return;

            BoSSS.Solution.Application._Main(
                args,
                true,
                null,
                () => new AdaptiveMeshRefinementTestMain());
        }

        protected override GridCommons CreateOrLoadGrid() {
            double[] nodes = GenericBlas.Linspace(-5, 5, 61);
            var grd = Grid2D.Cartesian2DGrid(nodes, nodes);
            base.m_GridPartitioningType = GridPartType.none;
            return grd;
        }

        SinglePhaseField u;


        protected override void CreateFields() {
            
            var Basis = new Basis(base.GridData, DEGREE);
            u = new SinglePhaseField(Basis, "u");
            base.m_RegisteredFields.Add(u);
        }

        /// <summary>
        /// turn dynamic balancing on/off
        /// </summary>
        internal bool DynamicBalance = false;

        /// <summary>
        /// DG polynomial degree
        /// </summary>
        internal int DEGREE = 3;

        /// <summary>
        /// Setting initial value.
        /// </summary>
        protected override void SetInitial() {
        }


        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
            
            if (L == null) {
                
            } else {
                throw new NotImplementedException("todo");
            }
        }




        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (new FuncTrace()) {
                // Elo
                Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime + " ... ");

                // timestepping params
                base.NoOfTimesteps = 20;
                dt = 0.1;


                // return
                //base.TerminationKey = true;
                return dt;
            }
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "AdaptiveMeshRefinementTestMain." + timestepNo;
            Tecplot.PlotFields(base.m_RegisteredFields.ToArray(), filename, physTime, superSampling);
        }
    }
}

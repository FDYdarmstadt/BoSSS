using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
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

        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;

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

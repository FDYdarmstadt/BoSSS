using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using ilPSP;
using System.Diagnostics;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.Statistic {
    
    /// <summary>
    /// test for <see cref="CellLocalization"/>
    /// </summary>
    [TestFixture]
    class CellLocalizationTest : BoSSS.Solution.Application {

        static void Main() {
            var t = new CellLocalizationTest();
            t.SetUp();
            t.CheckNorm();
            t.Finish();
        }

        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public void SetUp() {
            bool dommy;
            Enviroment.Bootstrap(
                new string[0],
                GetBoSSSInstallDir(),
                //Path.Combine(GetBoSSSInstallDir(), Path.Combine("bin", Path.Combine("native", "win"))),
                out dommy);            

            base.Init(null, new CommandLineOptions(), "");
            base.CreateContextAndWhatFollows();
            this.SetInitial();
            this.CreateEquationsAndSolvers();
        }

        /// <summary>
        /// MPI shutdown
        /// </summary>
        [TestFixtureTearDown]
        public void Finish() {
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// checks the evaluation of the jump-seminorm
        /// </summary>
        [Test]
        public void CheckNorm() {
            double val = f.JumpSemiNorm();

            // Length of one edge
            double hF = _Context.GridDat.h_max_Edge[0];
            for (int e = 0; e < _Context.GridDat.NoOfEdges; e++) 
                Debug.Assert(Math.Abs(_Context.GridDat.h_max_Edge[e] - 2.0 / ((double)Res)) <= 1.0e-9);

            // Referenz-Lösung (für 
            
            double valRef = Math.Sqrt(Res * JumpHeight * JumpHeight * hF / hF); 

            // check
            double err = Math.Abs(val - valRef);
            double thres = 1.0e-10;
            Assert.IsTrue(err < thres);
        }


        /// <summary>
        /// must be even !
        /// </summary>
        const int Res = 100;


        const double JumpHeight = 2.0;

        protected override void CreateOrLoadGrid() {
            Debug.Assert(Res % 2 == 0);

            double[] xNodes = Grid1D.Linspace(-1, 1, Res + 1);
            double[] yNodes = Grid1D.Linspace(-1, 1, Res + 1);
            GridCommons grd = new Cartesian2DGrid(m_Context.CommMaster, xNodes, yNodes, null, null, false, false);
            _Context.SetGrid(grd);
        }

        
        Field f;

        protected override void SetInitial() {
            f.ProjectField((x, y) => ((x > 0 && y > 0) ? JumpHeight : 0.0));
        }

        protected override void CreateFields() {
            f = new SinglePhaseField(new Basis(_Context, 2));
        }


        protected override void CreateEquationsAndSolvers() {
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            return 0.0;
        }

        protected override void PlotCurrentState(double physTime, int timestepNo) {
        }
    }
}

}

using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using ilPSP.Utils;

namespace AdvancedSolverTests.SubBlocking {
    class SubBlockTests : TestBench {
        public override void Run() {
            AritmethicTest();
            
        }

        [Test]
        public void AritmethicTest() {
            
            //int Res=8;
            //double[] xNodes = GenericBlas.Linspace(-1, +1, Res + 1);
            //double[] yNodes = GenericBlas.Linspace(-1, +1, Res + 1);
            //double[] zNodes = GenericBlas.Linspace(-1, +1, Res + 1);
            //int J = (xNodes.Length - 1) * (yNodes.Length - 1) * (zNodes.Length - 1);

            //string GridName = string.Format(WorkflowMgm.CurrentProject + "_J" + J);

            //grids[cnt] = null;
            //foreach (IGridInfo grd in tempDB.Grids) {
            //    bool check = grd.Name.Contains(string.Format("_J" + J));
            //    if (check) {
            //        grids[cnt] = grd;
            //    }
            //}

            //AggregationGridData[] mgs = new AggregationGridData[] { };
            //XdgAggregationBasis xdgbasis = XdgAggregationBasis.CreateSequence(mgs, dgBasis); ;
            //XDGField field = new XDGField(xdgbasis.XDGBasis);
            
            //UnsetteledCoordinateMapping ucm = new UnsetteledCoordinateMapping();
            //MultigridMapping mgm = new MultigridMapping(ucm, new XdgAggregationBasis[] { xdgbasis }, new int[] { 2 });
            //MsrMatrix OperatorM;
            
        }

    }
}

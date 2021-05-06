using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CutEdgeQuadrature {
    class Tester3D {

        GridCommons grid;

        int innerEdgeIndex;

        public Tester3D() {
            double[] xNodes = new double[] { -1, 0, 1 };
            double[] yNodes = new double[] { -1, 1 };
            double[] zNodes = new double[] { -1, 1 };

            grid = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);
            innerEdgeIndex = 0;
        }

        public double Test(IEdgeQuadratureTest3D testcase, int order) {
            LevelSet phi = new LevelSet(new BoSSS.Foundation.Basis(grid.GridData, order), "phi");
            phi.ProjectField(testcase.LevelSet);
            Tecplot plotter = new Tecplot(grid.GridData, 4);
            plotter.PlotFields("Phi", 0, phi);
            LevelSetTracker lsTrkr = new LevelSetTracker(grid.GridData, testcase.MomentFittingVariant, 1, new string[] { "A", "B"}, phi);
            SpeciesId[] species = new SpeciesId[] { lsTrkr.GetSpeciesId("A"), lsTrkr.GetSpeciesId("B") };
            XDGSpaceMetrics metrics = lsTrkr.GetXDGSpaceMetrics(species, order);
            
            double s = EvaluateEdgeArea(metrics, species[0]);
            double error = Math.Abs(testcase.EdgeArea - s);
            return error;
        }

        double EvaluateEdgeArea(XDGSpaceMetrics metrics, SpeciesId species) {
            MultidimensionalArray areas = metrics.CutCellMetrics.CutEdgeAreas[species];
            return areas[innerEdgeIndex];
        }
    }
}

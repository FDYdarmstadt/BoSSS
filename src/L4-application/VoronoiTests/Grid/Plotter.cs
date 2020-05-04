using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Tecplot;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace VoronoiTests.Grid
{
    static class Plotter
    {
        public static void Plot(IGrid grid)
        {
            EdgeMask boundaryEdges = EdgeMask.GetFullMask(grid.iGridData, MaskType.Logical);
            boundaryEdges.SaveToTextFile(
                "edges.txt",
                false,
                (double[] CoordGlobal, int LogicalItemIndex, int GeomItemIndex) => grid.iGridData.iGeomEdges.EdgeTags[GeomItemIndex]);
            Tecplot plt1 = new Tecplot(grid.iGridData, true, false, 0);
            Basis b = new Basis(grid.iGridData, 0);
            SinglePhaseField field = new SinglePhaseField(b, "u");
            plt1.PlotFields("grid", 0, field);
        }
    }
}

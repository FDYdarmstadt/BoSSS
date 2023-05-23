using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using BoSSS.Foundation;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid.RefElements;

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
            EdgeMask edges = EdgeMask.GetFullMask(grid.iGridData, MaskType.Geometrical);
            MultidimensionalArray centers = GetCenters(edges);
            for(int i = 0; i < edges.NoOfItemsLocally; ++i) {
                Console.WriteLine($"Edge center x:{centers[i,0,0]} " +
                    $"y:{centers[i, 0, 1]} z:{centers[i, 0, 2]} with area: {areas[i]}");
            }
            return areas[innerEdgeIndex];
        }

        MultidimensionalArray GetCenters(EdgeMask edges) {
            int D = grid.SpatialDimension;
            MultidimensionalArray localCenterEdge = MultidimensionalArray.Create(1, D - 1);
            MultidimensionalArray localCenterVolume = MultidimensionalArray.Create(1, D);

            MultidimensionalArray centers = MultidimensionalArray.Create(edges.NoOfItemsLocally,1, D);
            GridData gridDat = grid.GridData;
            foreach (Chunk chunk in edges) {
                for (int edge = 0; edge < chunk.Len; edge++) {
                    int iTrafo = gridDat.iGeomEdges.Edge2CellTrafoIndex[edge, 0];
                    int localEdge = gridDat.iGeomEdges.FaceIndices[edge, 0];
                    int cell = gridDat.iGeomEdges.CellIndices[edge, 0];
                    RefElement KrefCell = gridDat.iGeomCells.GetRefElement(cell);

                    gridDat.iGeomEdges.Edge2CellTrafos[iTrafo].Transform(localCenterEdge, localCenterVolume);

                    gridDat.TransformLocal2Global(new NodeSet(KrefCell, localCenterVolume, false), cell, 1, centers, edge);

                }
            }
            return centers;
        }
    }
}

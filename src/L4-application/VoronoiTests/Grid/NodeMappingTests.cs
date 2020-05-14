using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Voronoi;
using BoSSS.Foundation.Grid.Voronoi.Meshing;
using BoSSS.Foundation.Voronoi;
using BoSSS.Solution.Tecplot;
using ilPSP;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace VoronoiTests.Grid
{
    class NodeMappingTests : BoSSSTestBench
    {
        public override void Run()
        {
            MapMultipleCells();
        }

        [Test]
        public void MapOneCellPeriodic()
        {
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 1, "Dirichlet" }
            };
            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames,
                BoundingBox = GridShapes.Rectangle(2, 2),
            };
            NodeTrackingVoronoiMesher mesher = new NodeTrackingVoronoiMesher( 
                new VoronoiMesher<TrackableNode>.Settings() 
                { 
                    Boundary = gridBoundary,
                    NumberOfLloydIterations = 0,
                });
            MultidimensionalArray array = MultidimensionalArray.Create(8, 2);
            array.SetRow(0, new double[] { -0.9, 0.7 });
            array.SetRow(1, new double[] { -0.8, -0.4 });
            array.SetRow(2, new double[] { 0, 0.5 }); 
            array.SetRow(3, new double[] { 0.05, 0 }); 
            array.SetRow(4, new double[] { 0.1, -0.5 });
            array.SetRow(5, new double[] { 0.8, 0.8 });
            array.SetRow(6, new double[] { 1.2, 0.2 });
            array.SetRow(7, new double[] { 0.5, -0.6 });

            VoronoiNodes nodes = new VoronoiNodes(array);
            MappedVoronoiGrid grid = mesher.CreateGrid(nodes, 0);
            Assert.IsTrue((
                grid.Result.Nodes.Nodes[grid.InputNodesToResultNodes.GetMapping(0)].Position 
                - new Vector(-0.9, 0.7)).Abs() < 1e-12);
        }

        [Test]
        public void MapMultipleCellsPeriodic()
        {
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 1, "Dirichlet" }
            };
            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames,
                BoundingBox = GridShapes.Rectangle(2, 2)
            };
            NodeTrackingVoronoiMesher mesher = new NodeTrackingVoronoiMesher(
                new VoronoiMesher<TrackableNode>.Settings()
                {
                    Boundary = gridBoundary,
                    NumberOfLloydIterations = 0,
                });
            int numberOfNodes = 200;
            MultidimensionalArray array = MultidimensionalArray.Create(numberOfNodes, 2);
            Random random = new Random(0);
            for(int i = 0; i < numberOfNodes; ++i)
            {
                array[i, 0] = 2.0 * random.NextDouble() - 1.0;
                array[i, 1] = 2.0 * random.NextDouble() - 1.0;
            }
            VoronoiNodes nodes = new VoronoiNodes(array);
            MappedVoronoiGrid grid = mesher.CreateGrid(nodes, 0);
            bool mapsNodesOntoItself = IsPermutation(grid.InputNodesToResultNodes, nodes, grid.Result.Nodes);
            Assert.IsTrue(mapsNodesOntoItself);
            MappedVoronoiGrid grid2 = mesher.CreateGrid(grid.Result.Nodes, 0);
            mapsNodesOntoItself = IsPermutation(grid2.InputNodesToResultNodes, grid.Result.Nodes, grid2.Result.Nodes);
            Assert.IsTrue(mapsNodesOntoItself);
        }

        bool IsPermutation(Map map, VoronoiNodes source, VoronoiNodes target)
        {
            for(int i = 0; i < source.Nodes.Count; ++i)
            {
                VoronoiNode sourceNode = source.Nodes[i];
                VoronoiNode targetNode = target.Nodes[map.GetMapping(i)];
                if((sourceNode.Position - targetNode.Position).Abs() > 1e-10)
                {
                    return false;
                }
            }

            return true;
        }
        
        [Test]
        public void MapOneCell()
        {
            byte[] tags = { 1, 1, 1, 1 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 1, "Dirichlet" }
            };
            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames,
                BoundingBox = GridShapes.Rectangle(2, 2)
            };
            NodeTrackingVoronoiMesher mesher = new NodeTrackingVoronoiMesher(
                new VoronoiMesher<TrackableNode>.Settings()
                {
                    Boundary = gridBoundary,
                    NumberOfLloydIterations = 0,
                });
            MultidimensionalArray array = MultidimensionalArray.Create(8, 2);
            array.SetRow(0, new double[] { -0.9, 0.7 });
            array.SetRow(1, new double[] { -0.8, -0.4 });
            array.SetRow(2, new double[] { 0, 0.5 });
            array.SetRow(3, new double[] { 0.05, 0 });
            array.SetRow(4, new double[] { 0.1, -0.5 });
            array.SetRow(5, new double[] { 0.8, 0.8 });
            array.SetRow(6, new double[] { 0.9, 0.2 });
            array.SetRow(7, new double[] { 0.5, -0.6 });

            VoronoiNodes nodes = new VoronoiNodes(array);
            MappedVoronoiGrid grid = mesher.CreateGrid(nodes, 0);
            Assert.IsTrue((
                grid.Result.Nodes.Nodes[grid.InputNodesToResultNodes.GetMapping(0)].Position
                - new Vector(-0.9, 0.7)).Abs() < 1e-12);
        }

        [Test]
        public void MapMultipleCells()
        {
            byte[] tags = { 1, 1, 1, 1 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(1)
            {
                { 1, "Dirichlet" }
            };
            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames,
                BoundingBox = GridShapes.Rectangle(2, 2)
            };
            NodeTrackingVoronoiMesher mesher = new NodeTrackingVoronoiMesher(
                new VoronoiMesher<TrackableNode>.Settings()
                {
                    Boundary = gridBoundary,
                    NumberOfLloydIterations = 0,
                });
            int numberOfNodes = 200;
            MultidimensionalArray array = MultidimensionalArray.Create(numberOfNodes, 2);
            Random random = new Random(0);
            for (int i = 0; i < numberOfNodes; ++i)
            {
                array[i, 0] = 2.0 * random.NextDouble() - 1.0;
                array[i, 1] = 2.0 * random.NextDouble() - 1.0;
            }
            VoronoiNodes nodes = new VoronoiNodes(array);
            MappedVoronoiGrid grid = mesher.CreateGrid(nodes, 0);
            bool mapsNodesOntoItself = IsPermutation(grid.InputNodesToResultNodes, nodes, grid.Result.Nodes);
            Assert.IsTrue(mapsNodesOntoItself);
            MappedVoronoiGrid grid2 = mesher.CreateGrid(grid.Result.Nodes, grid.Result.FirstCornerNodeIndice);
            mapsNodesOntoItself = IsPermutation(grid2.InputNodesToResultNodes, grid.Result.Nodes, grid2.Result.Nodes);
            Assert.IsTrue(mapsNodesOntoItself);
        }
    }
}

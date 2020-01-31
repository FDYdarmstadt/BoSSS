using BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler;
using ilPSP;
using System.Collections.Generic;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Voronoi.Meshing.Cutter;
using BoSSS.Foundation.Grid.Voronoi.Meshing.MICHMesher;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class MeshGenerator<T>
         where T : ICloneable<T>, new()
    {
        readonly Cutter<T> cutter;

        readonly Vector[] boundingBox;

        readonly BoundaryLine[] boundaryLines;

        readonly PeriodicBoundaryHandler<T> boundaryHandler;

        public MeshGenerator(MeshingAlgorithm.Settings settings)
        {
            cutter = new Cutter<T>();
            boundingBox = settings.BoundingBox;
            boundaryLines = BoundaryLine.ToLines(settings.Boundary);
            boundaryHandler = new PeriodicBoundaryHandler<T>(settings.PeriodicMap);
        }

        public Domain<T> Generate(IList<T> nodes, MeshCell<T> firstCorner)
        {
            return Generate(nodes, firstCorner.ID);
        }

        public Domain<T> Generate(IList<T> nodes, int firstCornerNodeIndice)
        {
            Debug.Assert(nodes.Count > 0);

            Domain<T> mesh = CreateMeshFrom(nodes, firstCornerNodeIndice);
            if (boundaryHandler.ContainsPeriodicBoundaries)
            {
                nodes = boundaryHandler.CloneNodesAlongPeriodicBoundaries(mesh);
                mesh = CreateMeshFrom(nodes, mesh.Boundary.FirstCorner);
                boundaryHandler.RecomposePeriodicEdges(mesh);
            }
            return mesh;
        }

        Domain<T> CreateMeshFrom(IList<T> nodes, int firstCornerNodeIndice)
        {
            AddDistantBoundingNodes(nodes, boundingBox);

            IDMesh<T> mesh = MIConvexHullMeshGenerator.CreateMesh(nodes);
            Boundary<T> boundary = new Boundary<T>
            {
                BoundaryLines = boundaryLines,
                FirstCorner = mesh.Cells[firstCornerNodeIndice]
            };

            cutter.CutOut(mesh, boundary);
            return new Domain<T>
            {
                Mesh = mesh,
                Boundary = boundary
            };
        }

        Domain<T> CreateMeshFrom(IList<T> nodes, MeshCell<T> firstCorner)
        {
            return CreateMeshFrom(nodes, firstCorner.ID);
        }

        static void AddDistantBoundingNodes(IList<T> nodes, Vector[] boundingBox)
        {
            foreach (Vector corner in boundingBox)
            {
                T cornerNode = new T()
                {
                    Position = corner * 10
                };
                nodes.Add(cornerNode);
            }
        }
    }
}
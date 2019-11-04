using BoSSS.Foundation.Grid.Voronoi.Meshing.Recomposer;
using BoSSS.Platform.LinAlg;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class MeshGenerator<T>
         where T : ILocatable, new()
    {
        readonly Cutter<T> cutter;

        readonly Vector[] boundingBox;

        readonly BoundaryLine[] boundary;

        readonly int firstCellNode_indice;

        readonly PeriodicBoundaryHandler<T> boundaryHandler;

        public MeshGenerator(MeshingAlgorithm.Settings settings)
        {
            cutter = new Cutter<T>();
            boundingBox = settings.BoundingBox;
            boundary = BoundaryLine.ToLines(settings.Boundary);
            firstCellNode_indice = settings.FirstCellNode_indice;
            boundaryHandler = new PeriodicBoundaryHandler<T>(
                boundary, 
                settings.FirstCellNode_indice,
                settings.PeriodicBoundaryMap,
                settings.PeriodicTransformations);
        }

        public Mesh<T> Generate(IList<T> nodes)
        {
            Debug.Assert(nodes.Count > 0);
            IDMesh<T> mesh = CreateMeshFrom(nodes);
            if (boundaryHandler.ContainsPeriodicBoundaries)
            {
                nodes = boundaryHandler.CloneNodesAlongPeriodicBoundaries(mesh);
                mesh = CreateMeshFrom(nodes);
                boundaryHandler.RecomposePeriodicEdges(mesh);
            }
            return mesh;
        }

        IDMesh<T> CreateMeshFrom(IList<T> nodes)
        {
            AddDistantBoundingNodes(nodes, boundingBox);
            
            IDMesh<T> mesh = MIConvexHullMeshGenerator.CreateMesh(nodes);
            cutter.CutOut(mesh, boundary, firstCellNode_indice);
            
            return mesh;
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
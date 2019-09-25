using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class MeshGenerator<T>
         where T : IMesherNode, new()
    {
        readonly BoundaryCutter<T> cutter;

        readonly Vector[] boundingBox;

        readonly BoundaryLine[] boundary;

        readonly int firstCellNode_indice;

        readonly BoundaryHandler<T> boundaryHandler;

        public MeshGenerator(MeshingAlgorithm.Settings settings)
        {
            cutter = new BoundaryCutter<T>();
            boundingBox = settings.BoundingBox;
            boundary = BoundaryLine.ToLines(settings.Boundary);
            firstCellNode_indice = settings.FirstCellNode_indice;
            boundaryHandler = new BoundaryHandler<T>(boundary, settings.PeriodicBoundaryMap, settings.FirstCellNode_indice);
        }

        public Mesh<T> Generate(IList<T> nodes)
        {
            Debug.Assert(nodes.Count > 0);
            Mesh<T> mesh = CreateMeshFrom(nodes);
            if (boundaryHandler.ContainsPeriodicBoundaries)
            {
                nodes = boundaryHandler.CloneNodesAlongPeriodicBoundaries(mesh);
                mesh = CreateMeshFrom(nodes);
                boundaryHandler.RecomposePeriodicEdges(mesh);
            }
            return mesh;
        }

        Mesh<T> CreateMeshFrom(IList<T> nodes)
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
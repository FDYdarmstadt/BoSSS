using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public interface IMesherNode
    {
        Vector Position { get; set; }
    }

    static class MeshingAlgorithm
    {
        public class Settings
        {
            public Vector[] BoundingBox;
            public Vector[] Boundary;
            public int NumberOfLloydIterations = 10;
            public int FirstCellNode_indice = 0;
        }

        public static BoundaryMesh<T> ComputeMesh<T>(IList<T> nodeList, Settings settings)
            where T : IMesherNode, new()
        {
            // Create Voronoi mesh
            // =================================
            IEnumerator<Line> boundaryLines = Line.GetEnumerator(settings.Boundary);
            IntersectionMesh<T> mesh = null;

            for (int iLloyd = 0; iLloyd <= settings.NumberOfLloydIterations; ++iLloyd)
            {
                // Mesh generation
                //-------------------------------------
                AddDistantBoundingNodes(nodeList, settings.BoundingBox);
                mesh = IntersectionMeshGenerator.CreateMesh(nodeList, settings.FirstCellNode_indice);
                //Clip
                //-------------------------------------
                Intersecter.Intersect(mesh, boundaryLines);

                // Lloyds algorithm (Voronoi relaxation)
                // -------------------------------------
                IReadOnlyList<MeshCell<T>> cells = mesh.GetCells(); //GetInsideCell : Return Cells in order of MeshArray.
                if (iLloyd != settings.NumberOfLloydIterations)
                {
                    MoveNodesTowardsCellCenter(cells, ref settings.FirstCellNode_indice);
                }
                nodeList = mesh.GetNodes();
            }
            return mesh;
        }

        static void AddDistantBoundingNodes<T>(IList<T> nodes, Vector[] boundingBox)
            where T : IMesherNode, new()
        {
            Debug.Assert(nodes.Count > 0);
            foreach (Vector corner in boundingBox)
            {
                T cornerNode = new T();
                cornerNode.Position = corner * 10;
                nodes.Add(cornerNode);
            }
        }

        static void MoveNodesTowardsCellCenter<T>(IReadOnlyList<MeshCell<T>> Cells, ref int FirstCellNode_indice)
            where T : IMesherNode, new()
        {
            //Mark inside nodes
            //Use Original Nodes List and update. Use LinkedList?! Only iterate and cut some nodes, or insert
            //Let's give it a try!
            for (int i = 0; i < Cells.Count; ++i)
            {
                MeshCell<T> cell = Cells[i];
                double relaxValue = 0.1;
                Vector CenterOfGravity = new Vector(0, 0);
                foreach (Vertex vertex in cell.Vertices)
                {
                    CenterOfGravity += vertex.Position;
                }
                CenterOfGravity.Scale(1.0 / cell.Vertices.Length);
                CenterOfGravity = CenterOfGravity * relaxValue + new Vector(cell.Node.Position) * (1 - relaxValue);

                cell.Node.Position = CenterOfGravity;
                if (cell.ID == FirstCellNode_indice)
                {
                    FirstCellNode_indice = i;
                }
            }
        }
    }
}

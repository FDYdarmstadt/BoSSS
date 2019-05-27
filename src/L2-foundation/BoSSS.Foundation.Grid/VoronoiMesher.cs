using System;
using System.Collections.Generic;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi
{
    static class VoronoiMesher<T>
        where T : INode, new()
    {
        public static BoundaryMesh<T> Create(IList<T> nodeList, VoronoiMesherInfo info)
        {
            // Create Voronoi mesh
            // =================================
            IEnumerator<Line> boundaryLines = Line.GetEnumerator(info.Boundary);
            IntersectionMesh<T> mesh = null;

            for (int iLloyd = 0; iLloyd <= info.NumberOfLloydIterations; ++iLloyd)
            {
                // Mesh generation
                //-------------------------------------
                AddDistantBoundingNodes(nodeList, info.BoundingBox);
                mesh = IntersectionMeshGenerator.CreateMesh(nodeList, info.FirstCellNode_indice);
                //Clip
                //-------------------------------------
                Intersecter.Intersect(mesh, boundaryLines);

                // Lloyds algorithm (Voronoi relaxation)
                // -------------------------------------
                IReadOnlyList<Cell<T>> cells = mesh.GetCells(); //GetInsideCell : Return Cells in order of MeshArray.
                if (iLloyd != info.NumberOfLloydIterations)
                {
                    MoveNodesTowardsCellCenter(cells, ref info.FirstCellNode_indice); 
                }
                nodeList = mesh.GetNodes();
            }
            return mesh;
        }

        static void AddDistantBoundingNodes(IList<T> nodes, Vector[] boundingBox)
        {
            Debug.Assert(nodes.Count > 0);
            foreach(Vector corner in boundingBox)
            {
                T cornerNode = new T();
                cornerNode.Position = corner * 10;
                nodes.Add(cornerNode);
            }
        }

        static void MoveNodesTowardsCellCenter(IReadOnlyList<Cell<T>> Cells, ref int FirstCellNode_indice)
        {
            //Mark inside nodes
            //Use Original Nodes List and update. Use LinkedList?! Only iterate and cut some nodes, or insert
            //Let's give it a try!
            for (int i = 0; i < Cells.Count; ++i)
            {
                Cell<T> cell = Cells[i];
                double relaxValue = 0.1;
                Vector CenterOfGravity = new Vector(0, 0);
                foreach (Vertex vertex in cell.Vertices)
                {
                    CenterOfGravity += vertex.Position;
                }
                CenterOfGravity.Scale(1.0 / cell.Vertices.Length);
                CenterOfGravity = CenterOfGravity * relaxValue + new Vector(cell.Position) * (1 - relaxValue);

                cell.Position = CenterOfGravity;
                cell.Node.Position = CenterOfGravity;
                if (cell.ID == FirstCellNode_indice)
                {
                    FirstCellNode_indice = i;
                }
            }
        }
    }
}

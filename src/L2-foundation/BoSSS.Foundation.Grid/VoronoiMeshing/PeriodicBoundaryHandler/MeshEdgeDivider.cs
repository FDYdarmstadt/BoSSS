using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class MeshEdgeDivider<T>
    {
        readonly IDMesh<T> mesh;

        public MeshEdgeDivider(IDMesh<T> mesh)
        {
            this.mesh = mesh;
        }

        public void DivideEdge(Edge<T> edge, int edgeIndiceInCell)
        {
            Edge<T> clone = edge = MeshElementCloner.Clone(edge);
            mesh.AddVertex(clone.Start);
            mesh.AddVertex(clone.End);
            SwitchEdges(edge, edgeIndiceInCell, clone);
        }

        static void SwitchEdges(Edge<T> old, int indiceInCell, Edge<T> newEdge)
        {
            MeshCell<T> cell = old.Cell;
            if(cell.ID != newEdge.Cell.ID)
            {
                throw new Exception("Edges are not registered to the same cell.");
            }
            Edge<T> before = cell.Edges[(indiceInCell - 1 + cell.Edges.Length) % cell.Edges.Length];
            Edge<T> after = cell.Edges[(indiceInCell + 1 + cell.Edges.Length) % cell.Edges.Length];

            cell.Edges[indiceInCell] = newEdge;
            before.End = newEdge.Start;
            after.Start = newEdge.End;
            UpdateVerticesOf(cell);

            newEdge.Twin = old.Twin;
            newEdge.Twin.Twin = newEdge;

            SwitchVertex(before.Twin.Start, before.Twin.Cell, newEdge.Start);
            SwitchVertex(after.Twin.End, after.Twin.Cell, newEdge.End);
        }

        static void SwitchVertex(Vertex old, MeshCell<T> cell, Vertex newVertex)
        {
            (Edge<T> before, Edge<T> ofVertex) = GetEdgePair(old, cell);
            ofVertex.Start = newVertex;
            before.End = newVertex;
            UpdateVerticesOf(cell);
        }

        static (Edge<T> before, Edge<T> ofVertex) GetEdgePair(Vertex vertex, MeshCell<T> cell)
        {
            for (int i = 0; i < cell.Edges.Length; ++i)
            {
                Edge<T> current = cell.Edges[i];
                if (current.Start.ID == vertex.ID)
                {
                    Edge<T> before = cell.Edges[(i - 1 + cell.Edges.Length) % cell.Edges.Length];
                    return (before, current);
                }
            }
            throw new Exception("Vertex not found in Cell");
        }

        static void UpdateVerticesOf(MeshCell<T> cell)
        {
            for (int i = 0; i < cell.Edges.Length; ++i)
            {
                cell.Vertices[i] = cell.Edges[i].Start;
            }
        }
    }
}

using BoSSS.Platform.LinAlg;
using System.Collections.Generic;
using System.Collections;
using System.Diagnostics;
using BoSSS.Foundation.Voronoi;
using System;

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
            public HashSet<int> NotchedBoundaryLines = new HashSet<int>();
            public int NumberOfLloydIterations = 10;
            public int FirstCellNode_indice = 0;
        }

        static void AssertCorrectness(Settings settings)
        {
            Debug.Assert(settings.BoundingBox.Length == 4);
            Debug.Assert(settings.Boundary.Length > 2);
        }

        public static IDMesh<T> ComputeMesh<T>(IList<T> nodeList, Settings settings)
            where T : IMesherNode, new()
        {
            AssertCorrectness(settings);
            BoundaryLineEnumerator boundaryLines = BoundaryLine.GetEnumerator(
                settings.Boundary, 
                settings.NotchedBoundaryLines);
            IDMesh<T> mesh = null;
            var cutter = new BoundaryCutter<T>();

            // Create Voronoi mesh
            // =================================
            for (int iLloyd = 0; iLloyd <= settings.NumberOfLloydIterations; ++iLloyd)
            {
                // Voronoi Mesh generation
                //-------------------------------------
                AddDistantBoundingNodes(nodeList, settings.BoundingBox);
                mesh = MIConvexHullMeshGenerator.CreateMesh(nodeList);

                //Clip: cut out boundary polygon
                //-------------------------------------
                cutter.CutOut(mesh, boundaryLines, settings.FirstCellNode_indice);

                // Lloyds algorithm (Voronoi relaxation)
                // -------------------------------------
                //Get inside cells : Return cells in order of Mesh Array.
                if (iLloyd != settings.NumberOfLloydIterations)
                {
                    MoveNodesTowardsCellCenter(mesh.Cells, ref settings.FirstCellNode_indice);
                }
                nodeList = mesh.Nodes;
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
                Vector centerOfGravity = CenterOf(cell);
                double relaxValue = 0.1;
                centerOfGravity = centerOfGravity * relaxValue + new Vector(cell.Node.Position) * (1 - relaxValue);
                cell.Node.Position = centerOfGravity;

                if (cell.ID == FirstCellNode_indice)
                {
                    FirstCellNode_indice = i;
                }
            }
        }

        static Vector CenterOf<T>(MeshCell<T> cell)
        {
            Vector centerOfGravity = new Vector(0, 0);
            foreach (Vertex vertex in cell.Vertices)
            {
                centerOfGravity += vertex.Position;
            }
            centerOfGravity.Scale(1.0 / cell.Vertices.Length);
            return centerOfGravity;
        }
    }
}

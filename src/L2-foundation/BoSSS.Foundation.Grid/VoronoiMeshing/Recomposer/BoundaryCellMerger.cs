using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryCellMerger<T>
       where T : IMesherNode
    {
        readonly LinkedList<Edge<T>> newEdges;

        IEnumerator<Edge<T>> sourceEdgeEnumerator;

        IEnumerator<Edge<T>> targetEdgeEnumerator;

        int targetBoundaryEdgeNumber;

        public BoundaryCellMerger()
        {
            newEdges = new LinkedList<Edge<T>>();
        }

        public static void MergeAtBoundary(
            IList<MeshCell<T>> source,
            IList<MeshCell<T>> target,
            IList<(int sourceEdgeIndice, int targetEdgeIndice)> glueMap)
        {
            BoundaryCellMerger<T> merger = new BoundaryCellMerger<T>();
            if (source.Count != target.Count && source.Count != glueMap.Count)
            {
                throw new Exception("lenghts must align.");
            }
            for (int i = 0; i < source.Count; ++i)
            {
                merger.Merge(source[i], glueMap[i].sourceEdgeIndice, target[i], glueMap[i].targetEdgeIndice);
            }
        }

        public void Merge(MeshCell<T> source, int sourceEdgeIndice, MeshCell<T> target, int targetEdgeIndice)
        {
            //Wer ist der coolste hier?
            AssertCorrectness(source, sourceEdgeIndice, target, targetEdgeIndice);
            newEdges.Clear();

            //Todo: Add convolution stuff, if weld is longer than 1.
            SetupBeforeWeld(source.Edges, sourceEdgeIndice, target.Edges, targetEdgeIndice);
            ConstructNewEdges();
            ReshapeCell(target);
        }

        void AssertCorrectness(MeshCell<T> source, int sourceIndice, MeshCell<T> target, int targetIndice)
        {
            Debug.Assert(source.Edges.Length > 2 && source.Edges.Length > sourceIndice);
            Debug.Assert(target.Edges.Length > 2 && target.Edges.Length > targetIndice);
        }

        void SetupBeforeWeld(Edge<T>[] sourceEdges, int sourceEdgeIndice, Edge<T>[] targetEdges, int targetEdgeIndice)
        {
            targetBoundaryEdgeNumber = targetEdges[targetEdgeIndice].BoundaryEdgeNumber;

            CyclicArray<Edge<T>> sourceEdgesStartingAtWeld = new CyclicArray<Edge<T>>(sourceEdges, sourceEdgeIndice);
            sourceEdgeEnumerator = new CyclicArrayEnumerator<Edge<T>>(sourceEdgesStartingAtWeld);

            Edge<T>[] reversedTargetEdges = ArrayMethods.GetCopyInReverseOrder(targetEdges);
            int reversedInidce = targetEdges.Length - 1 - targetEdgeIndice;
            CyclicArray<Edge<T>> targetEdgesStartingAtWeld = new CyclicArray<Edge<T>>(reversedTargetEdges, reversedInidce);
            targetEdgeEnumerator = new CyclicArrayEnumerator<Edge<T>>(targetEdgesStartingAtWeld);
        }

        void ConstructNewEdges()
        {
            MergeBoundaryEdges();
            WeldEdgesAfterMerge();
            AddFollowingEdges();
            WeldEdgesBeforeMerge();
        }

        void MergeBoundaryEdges()
        {
            //Remove boundary edges
            sourceEdgeEnumerator.MoveNext();
            targetEdgeEnumerator.MoveNext();
        }

        void WeldEdgesAfterMerge()
        {
            //Weld edges After weld
            sourceEdgeEnumerator.MoveNext();
            targetEdgeEnumerator.MoveNext();
            AddSourceEdge(sourceEdgeEnumerator.Current);
            AddTargetEdge(targetEdgeEnumerator.Current);
            LinkedListNode<Edge<T>> sourceNodeAfterWeld = newEdges.First;
            LinkedListNode<Edge<T>> targetNodeAfterWeld = newEdges.Last;
            WeldEdges(sourceNodeAfterWeld, targetNodeAfterWeld);
        }

        void AddTargetEdge(Edge<T> targetEdge)
        {
            newEdges.AddLast(targetEdge);
        }

        void AddSourceEdge(Edge<T> sourceEdge)
        {
            sourceEdge.IsBoundary = true;
            sourceEdge.BoundaryEdgeNumber = targetBoundaryEdgeNumber;
            newEdges.AddFirst(sourceEdge);
        }

        void AddFollowingEdges()
        {
            //Add common edges and set source to target boundary
            while (sourceEdgeEnumerator.MoveNext())
            {
                Edge<T> sourceEdge = sourceEdgeEnumerator.Current;
                AddSourceEdge(sourceEdge);
            }
            while (targetEdgeEnumerator.MoveNext())
            {
                Edge<T> targetEdge = targetEdgeEnumerator.Current;
                AddTargetEdge(targetEdge);
            }
        }

        void WeldEdgesBeforeMerge()
        {
            //Last edge and fist edge are edges before weld
            var firstNode = newEdges.First;
            newEdges.RemoveFirst();
            newEdges.AddLast(firstNode);
            WeldEdges(newEdges.Last.Previous, newEdges.Last);
        }

        void WeldEdges(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            if (OnLine(source.Value.Start, source.Value.End, target.Value.Start))
            {
                target.Value.IsBoundary = false;
                target.Value.BoundaryEdgeNumber = source.Value.BoundaryEdgeNumber;
                target.Value = EdgeMerge(source.Value, target.Value);
                MergeNodes(source, target);
            }
        }

        static void MergeNodes(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            Debug.Assert(source.List.First.Value.Start.ID == target.List.First.Value.Start.ID, 
                "Nodes must be from same list.");
            source.List.Remove(source);
        }

        static Edge<T> EdgeMerge(Edge<T> source, Edge<T> target)
        {
            Debug.Assert((source.Start.Position - target.End.Position).Abs() < 1e-12,
                "Edges do not touch.");
            target.End = source.End;
            return target;
        }

        const double accuracy = 1e-12;

        static bool OnLine(Vertex a, Vertex b, Vertex c)
        {
            Vector vA = a.Position;
            Vector vB = b.Position;
            Vector vC = c.Position;
            double orientation = (vB.x - vA.x) * (vC.y - vA.y) - (vB.y - vA.y) * (vC.x - vA.x);

            if (Math.Abs(orientation) < accuracy)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        void ReshapeCell(MeshCell<T> cell)
        {
            cell.Edges = new Edge<T>[newEdges.Count];
            LinkedListNode<Edge<T>> current = newEdges.Last;
            int j = 0;
            while (current != null)
            {
                cell.Edges[j] = current.Value;
                ++j;
                current = current.Previous;
            }

            Vertex[] vertices = new Vertex[cell.Edges.Length];
            for (int i = 0; i < vertices.Length; ++i)
            {
                vertices[i] = cell.Edges[i].Start;
            }
            cell.Vertices = vertices;
        }
    }
}

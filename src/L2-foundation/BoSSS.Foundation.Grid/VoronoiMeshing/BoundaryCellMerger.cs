using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryCellMerger<T>
       where T : IMesherNode
    {
        LinkedList<Edge<T>> newEdges;

        int sourceBoundaryEdgeNumber;

        int targetBoundaryEdgeNumber;

        IEnumerator<Edge<T>> sourceEdgeEnumerator;

        IEnumerator<Edge<T>> targetEdgeEnumerator;

        public BoundaryCellMerger(int sourceBoundaryEdgeNumber, int targetBoundaryEdgeNumber)
        {
            newEdges = new LinkedList<Edge<T>>();
            this.sourceBoundaryEdgeNumber = sourceBoundaryEdgeNumber;
            this.targetBoundaryEdgeNumber = targetBoundaryEdgeNumber;
        }
        public static void MergeAtBoundary(
        IList<MeshCell<T>> source,
        int sourceBoundaryEdgeNumber,
        IList<MeshCell<T>> target,
        int targetBoundaryEdgeNumber)
        {
            BoundaryCellMerger<T> merger = new BoundaryCellMerger<T>(sourceBoundaryEdgeNumber, targetBoundaryEdgeNumber);
            if (source.Count != target.Count)
            {
                throw new Exception("lenghts must align.");
            }
            for (int i = 0; i < source.Count; ++i)
            {
                merger.MergeAtBoundary(source[i], target[i]);
            }
        }

        public void MergeAtBoundary(MeshCell<T> source, MeshCell<T> target)
        {
            //Wer ist der coolste hier?
            Initialize(source, target);

            LinkedListNode<Edge<T>> beginningOfWeld = CollectEdgesBeforeWeld();
            Weld(beginningOfWeld);
            CollectEdgesAfterWeld();
            ReshapeCell(source);
        }

        void Initialize(MeshCell<T> source, MeshCell<T> target)
        {
            AssertCorrectness(source, target);
            newEdges.Clear();
            sourceEdgeEnumerator = new ArrayEnumerator<Edge<T>>(source.Edges);
            targetEdgeEnumerator = new BackwardsEnumerator<Edge<T>>(target.Edges);
        }

        void AssertCorrectness(MeshCell<T> source, MeshCell<T> target)
        {
            Debug.Assert(source.Edges.Count(x => (x.BoundaryEdgeNumber == sourceBoundaryEdgeNumber)) == 1);
            Debug.Assert(target.Edges.Count(x => (x.BoundaryEdgeNumber == targetBoundaryEdgeNumber)) == 1);
        }

        LinkedListNode<Edge<T>> CollectEdgesBeforeWeld()
        {
            while (sourceEdgeEnumerator.MoveNext())
            {
                Edge<T> edge = sourceEdgeEnumerator.Current;
                if (edge.BoundaryEdgeNumber != sourceBoundaryEdgeNumber)
                {
                    ToTargetBoundary(edge);
                    newEdges.AddLast(edge);
                }
                else
                {
                    break;
                }
            }
            LinkedListNode<Edge<T>> beginningOfWeld = newEdges.Last;
            while (targetEdgeEnumerator.MoveNext())
            {
                Edge<T> edge = targetEdgeEnumerator.Current;
                if (edge.BoundaryEdgeNumber != targetBoundaryEdgeNumber)
                {
                    newEdges.AddAfter(beginningOfWeld, edge);
                }
                else
                {
                    break;
                }
            }
            return beginningOfWeld;
        }

        void Weld(LinkedListNode<Edge<T>> beginningOfWeld)
        {
            //Weld before cut
            LinkedListNode<Edge<T>> sourceNodeBeforeWeld = beginningOfWeld;
            LinkedListNode<Edge<T>> targetNodeBeforeWeld = beginningOfWeld.Next;
            WeldEdges(sourceNodeBeforeWeld, targetNodeBeforeWeld);

            //Remove
            sourceEdgeEnumerator.MoveNext();
            sourceEdgeEnumerator.MoveNext();
            
            //Find edges after weld
            if (sourceEdgeEnumerator.MoveNext())
            {
                Edge<T> sourceEdgeAfterWeld = sourceEdgeEnumerator.Current;
                ToTargetBoundary(sourceEdgeAfterWeld);
                newEdges.AddFirst(sourceEdgeAfterWeld);
            }
            if (targetEdgeEnumerator.MoveNext())
            {
                Edge<T> targetEdgeAfterWeld = targetEdgeEnumerator.Current;
                newEdges.AddLast(targetEdgeAfterWeld);
            }
            LinkedListNode<Edge<T>> sourceNodeAfterWeld = newEdges.First;
            LinkedListNode<Edge<T>> targetNodeAfterWeld = newEdges.Last;

            //Weld after cut
            newEdges.AddLast(sourceNodeAfterWeld);
            newEdges.RemoveFirst();
            WeldEdges(sourceNodeAfterWeld, targetNodeAfterWeld);
        }

        void WeldEdges(LinkedListNode<Edge<T>> node, LinkedListNode<Edge<T>> nextNode)
        {
            Debug.Assert((node.Value.End.Position - nextNode.Value.Start.Position).Abs() > 1e-12, "Edges do not touch.");
            if (OnLine(node.Value.Start, node.Value.End, nextNode.Value.Start))
            {
                Merge(node, nextNode);
            }
        }

        static void Merge(LinkedListNode<Edge<T>> node, LinkedListNode<Edge<T>> nextNode)
        {
            Debug.Assert(!node.List.Equals(nextNode.List), "Nodes must be from same list.");
            node.Value = EdgeMerge(node.Value, nextNode.Value);
            node.List.Remove(nextNode);
        }

        static Edge<T> EdgeMerge(Edge<T> source, Edge<T> target)
        {
            source.End = target.Start;
            return source;
        }

        const double accuracy = 1e-12;

        static bool OnLine(Vertex a, Vertex b, Vertex c)
        {
            Vector start = a.Position;
            Vector end = b.Position;
            Vector node = c.Position;
            double orientation =  ((end.x - start.x) * (node.y - start.y) - (end.y - start.y) * (node.x - start.x));

            if (Math.Abs(orientation) < accuracy)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        void CollectEdgesAfterWeld()
        {
            //-------------------------Close Edges
            LinkedListNode<Edge<T>> lastNode = newEdges.Last;
            while (sourceEdgeEnumerator.MoveNext())
            {
                Edge<T> edge = sourceEdgeEnumerator.Current;
                newEdges.AddAfter(lastNode, edge);
            }

            while (targetEdgeEnumerator.MoveNext())
            {
                Edge<T> edge = targetEdgeEnumerator.Current;
                newEdges.AddBefore(lastNode, edge);
            }
        }

        void ToTargetBoundary(Edge<T> edge)
        {
            edge.BoundaryEdgeNumber = targetBoundaryEdgeNumber;
        }

        void ReshapeCell(MeshCell<T> cell)
        {
            cell.Edges = newEdges.ToArray();

            Vertex[] vertices = new Vertex[cell.Edges.Length];
            for (int i = 0; i < vertices.Length; ++i)
            {
                vertices[i] = cell.Edges[i].Start;
            }
            cell.Vertices = vertices;
        }
    }
}

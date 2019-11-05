
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class BoundaryCellMerger<T>
       where T : ILocatable
    {
        readonly LinkedList<Edge<T>> newEdges;

        IEnumerator<Edge<T>> sourceEdgeEnumerator;

        IEnumerator<Edge<T>> targetEdgeEnumerator;

        int targetBoundaryEdgeNumber;

        MeshCell<T> targetCell;

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
            targetCell = target;

            //Todo: Add convolution stuff, if weld is longer than 1.
            Setup(source.Edges, sourceEdgeIndice, target.Edges, targetEdgeIndice);
            ConstructNewEdges();
            ReshapeCell(target);
        }

        void AssertCorrectness(MeshCell<T> source, int sourceIndice, MeshCell<T> target, int targetIndice)
        {
            Debug.Assert(source.Edges.Length > 2 && source.Edges.Length > sourceIndice);
            Debug.Assert(target.Edges.Length > 2 && target.Edges.Length > targetIndice);
        }

        void Setup(Edge<T>[] sourceEdges, int sourceEdgeIndice, Edge<T>[] targetEdges, int targetEdgeIndice)
        {
            targetBoundaryEdgeNumber = targetEdges[targetEdgeIndice].BoundaryEdgeNumber;

            CyclicArray<Edge<T>> sourceEdgesStartingAtWeld = new CyclicArray<Edge<T>>(sourceEdges, sourceEdgeIndice);
            sourceEdgeEnumerator = new CyclicArrayEnumerator<Edge<T>>(sourceEdgesStartingAtWeld);

            Edge<T>[] reversedTargetEdges = ArrayMethods.GetReverseOrderArray(targetEdges);
            int reversedInidce = targetEdges.Length - 1 - targetEdgeIndice;
            CyclicArray<Edge<T>> targetEdgesStartingAtWeld = new CyclicArray<Edge<T>>(reversedTargetEdges, reversedInidce);
            targetEdgeEnumerator = new CyclicArrayEnumerator<Edge<T>>(targetEdgesStartingAtWeld);
        }

        void ConstructNewEdges()
        {
            SkipEdgesWhereCellsAreGlued();
            WeldFirstEdges();
            AddCommonEdges();
            WeldLastEdges();
        }

        void SkipEdgesWhereCellsAreGlued()
        {
            //Remove boundary edges
            sourceEdgeEnumerator.MoveNext();
            targetEdgeEnumerator.MoveNext();
        }

        void WeldFirstEdges()
        {
            //Weld edges After weld
            sourceEdgeEnumerator.MoveNext();
            targetEdgeEnumerator.MoveNext();
            AddSourceEdge(sourceEdgeEnumerator.Current);
            
            AddTargetEdge(targetEdgeEnumerator.Current);
            LinkedListNode<Edge<T>> sourceNodeAfterWeld = newEdges.First;
            LinkedListNode<Edge<T>> targetNodeAfterWeld = newEdges.Last;
            Welder.TargetFirstWeldEdges(sourceNodeAfterWeld, targetNodeAfterWeld);
        }

        void AddTargetEdge(Edge<T> targetEdge)
        {
            newEdges.AddLast(targetEdge);
        }

        void AddSourceEdge(Edge<T> sourceEdge)
        {
            sourceEdge.Cell = targetCell;
            if (newEdges.Count != 0)
            {
                sourceEdge.Start = newEdges.First.Value.End;
            }
            Edge<T> neighborInUntransformed = sourceEdge.Twin.Twin.Twin;
            Edge<T> neighbor = sourceEdge.Twin;
            if (neighborInUntransformed.Start.ID != neighbor.Start.ID)
            {
                //pointer magic
                sourceEdge.Twin.Twin.Twin = sourceEdge.Twin;

                sourceEdge.Twin = neighborInUntransformed;
                //sourceEdge.Cell = sourceEdge.Twin.Twin.Cell;
                sourceEdge.Twin.Twin = sourceEdge;

                sourceEdge.Start.ID = sourceEdge.Twin.End.ID;
                sourceEdge.End.ID = sourceEdge.Twin.Start.ID;

                sourceEdge.IsBoundary = false;
                sourceEdge.BoundaryEdgeNumber = -1;
                sourceEdge.Twin.IsBoundary = false;
                sourceEdge.Twin.BoundaryEdgeNumber = -1;
            }
            else
            {
                sourceEdge.Twin.Twin = sourceEdge;
                sourceEdge.BoundaryEdgeNumber = targetBoundaryEdgeNumber;
                sourceEdge.IsBoundary = true;
            }
            newEdges.AddFirst(sourceEdge);
        }

        void AddCommonEdges()
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

        void WeldLastEdges()
        {
            //Attach last source Edge behind last target Edge 
            var firstNode = newEdges.First;
            newEdges.RemoveFirst();
            newEdges.AddLast(firstNode);
            //Weld last source Edge and last target Edge
            Welder.SourceFirstWeldEdges(newEdges.Last, newEdges.Last.Previous);
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

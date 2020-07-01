
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

        MeshCell<T> targetCell;

        public BoundaryCellMerger()
        {
            newEdges = new LinkedList<Edge<T>>();
        }

        public static void MergeAtBoundary(
            IList<(MeshCell<T> cell, bool)> source,
            IList<(MeshCell<T> cell, bool)> target,
            IList<(int sourceEdgeIndice, int targetEdgeIndice, bool glue)> glueMap)
        {
            BoundaryCellMerger<T> merger = new BoundaryCellMerger<T>();
            if (source.Count != target.Count && source.Count != glueMap.Count)
            {
                throw new Exception("lenghts must align.");
            }
            for (int i = 0; i < source.Count; ++i)
            {
                if (glueMap[i].glue)
                {
                    merger.Merge(source[i].cell, glueMap[i].sourceEdgeIndice, target[i].cell, glueMap[i].targetEdgeIndice);
                }
                else
                {
                    merger.Link(source[i].cell, glueMap[i].sourceEdgeIndice, target[i].cell, glueMap[i].targetEdgeIndice);
                }
            }
        }

        public void Link(MeshCell<T> source, int sourceEdgeIndice, MeshCell<T> target, int targetEdgeIndice)
        {
            source.Edges[sourceEdgeIndice].Twin = target.Edges[targetEdgeIndice];
            target.Edges[targetEdgeIndice].Twin = source.Edges[sourceEdgeIndice];
        }

        public void Merge(MeshCell<T> source, int sourceEdgeIndice, MeshCell<T> target, int targetEdgeIndice)
        {
            //Wer ist der coolste hier?
            AssertCorrectness(source, sourceEdgeIndice, target, targetEdgeIndice);
            newEdges.Clear();
            targetCell = target;

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
            Welder.TryTargetFirstWeldEdges(sourceNodeAfterWeld, targetNodeAfterWeld);
        }

        void AddTargetEdge(Edge<T> targetEdge)
        {
            newEdges.AddLast(targetEdge);
            targetEdge.Twin.Twin = targetEdge;
        }

        void AddSourceEdge(Edge<T> sourceEdge)
        {
            sourceEdge.Cell = targetCell;
            newEdges.AddFirst(sourceEdge);
        }

        void AddCommonEdges()
        {
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
            Welder.TrySourceFirstWeldEdges(newEdges.Last, newEdges.Last.Previous);
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

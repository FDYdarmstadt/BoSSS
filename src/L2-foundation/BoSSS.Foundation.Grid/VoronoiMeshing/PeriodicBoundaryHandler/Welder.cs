using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    static class Welder
    {
        public static void TargetFirstWeldEdges<T>(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            if (OnLine(source.Value.Start, source.Value.End, target.Value.Start))
            {
                target.Value = TargetFirstEdgeMerge(source.Value, target.Value);
                MergeNodes(source, target);
            }
            else
            {
                Console.WriteLine("Not merging");
            }
        }

        public static void SourceFirstWeldEdges<T>(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            if (OnLine(source.Value.Start, source.Value.End, target.Value.Start))
            {
                target.Value = SourceFirstEdgeMerge(source.Value, target.Value);
                MergeNodes(source, target);
            }
            else
            {
                Console.WriteLine("Not merging");
            }
        }

        static void MergeNodes<T>(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            Debug.Assert(source.List.First.Value.Start.ID == target.List.First.Value.Start.ID,
                "Nodes must be from same list.");
            source.List.Remove(source);
        }

        static Edge<T> SourceFirstEdgeMerge<T>(Edge<T> source, Edge<T> target)
        {
            Debug.Assert((source.End.Position - target.Start.Position).Abs() < 1e-10,
                "Edges do not touch.");
            target.Start = source.Start;
            if (target.IsBoundary)
            {
                if (source.IsBoundary)
                {
                    target.Twin = source.Twin;
                }
                else
                {
                    target.Twin.End = target.Start;
                }
            }
            return target;
        }

        static Edge<T> TargetFirstEdgeMerge<T>(Edge<T> source, Edge<T> target)
        {
            Debug.Assert((target.End.Position - source.Start.Position).Abs() < 1e-10,
                "Edges do not touch.");
            target.End = source.End;
            if (target.IsBoundary)
            {
                if (source.IsBoundary)
                {
                    target.Twin = source.Twin;
                }
                else
                {
                    target.Twin.Start = target.End;
                }
            }
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
    }
}

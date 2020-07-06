using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    static class Welder
    {
        public static void TryTargetFirstWeldEdges<T>(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            if (OnLine(target.Value.Start, source.Value.Start, source.Value.End))
            {
                target.Value = TargetFirstEdgeMerge(source.Value, target.Value);
                MergeNodes(source, target);
            }
            else
            {
                Console.WriteLine("Not merging");
                source.Value.Start = target.Value.End;
            }
        }

        public static void TrySourceFirstWeldEdges<T>(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            if (OnLine(source.Value.Start, source.Value.End, target.Value.End))
            {
                target.Value = SourceFirstEdgeMerge(source.Value, target.Value);
                MergeNodes(source, target);
            }
            else
            {
                Console.WriteLine("Not merging");
                source.Value.End = target.Value.Start;
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
            Debug.Assert(BoSSS.Platform.FloatingPointArithmetic.IsEqual(source.End.Position, target.Start.Position, 1e-8),
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
            Debug.Assert(BoSSS.Platform.FloatingPointArithmetic.IsEqual(target.End.Position, source.Start.Position, 1e-8),
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

        const double accuracy = 1e-5;

        static bool OnLine(Vertex a, Vertex b, Vertex c)
        {
            Vector vA = a.Position;
            Vector vB = b.Position;
            Vector vC = c.Position;
            double orientation = (vB.x - vA.x) * (vC.y - vA.y) - (vB.y - vA.y) * (vC.x - vA.x);

            double max = Math.Max(vA.AbsSquare(), Math.Max(vB.AbsSquare(), vC.AbsSquare()));
            if (Math.Abs(orientation) < accuracy * max)
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

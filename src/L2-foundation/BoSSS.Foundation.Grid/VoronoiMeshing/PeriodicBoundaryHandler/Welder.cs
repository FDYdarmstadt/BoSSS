using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class Welder
    {
        int targetBoundaryNumber;

        public Welder(int targetBoundaryNumber)
        {
            this.targetBoundaryNumber = targetBoundaryNumber;
        }

        public void TargetFirstWeldEdges<T>(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            if (OnLine(source.Value.Start, source.Value.End, target.Value.Start))
            {
                target.Value = TargetFirstEdgeMerge(source.Value, target.Value);
                MergeNodes(source, target);
            }
        }

        public void SourceFirstWeldEdges<T>(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            if (OnLine(source.Value.Start, source.Value.End, target.Value.Start))
            {
                target.Value = SourceFirstEdgeMerge(source.Value, target.Value);
                MergeNodes(source, target);
            }
        }

        void MergeNodes<T>(LinkedListNode<Edge<T>> source, LinkedListNode<Edge<T>> target)
        {
            Debug.Assert(source.List.First.Value.Start.ID == target.List.First.Value.Start.ID,
                "Nodes must be from same list.");
            source.List.Remove(source);
        }

        Edge<T> SourceFirstEdgeMerge<T>(Edge<T> source, Edge<T> target)
        {
            Debug.Assert((source.End.Position - target.Start.Position).Abs() < 1e-10,
                "Edges do not touch.");
            int ID = target.Start.ID; //hoho
            target.Start = source.Start;
            target.Start.ID = ID;
            if (target.IsBoundary && target.BoundaryEdgeNumber == targetBoundaryNumber)
            {
                source.Twin.Twin = target;
                target.Twin = source.Twin;
            }
            return target;
        }

        Edge<T> TargetFirstEdgeMerge<T>(Edge<T> source, Edge<T> target)
        {
            Debug.Assert((target.End.Position - source.Start.Position).Abs() < 1e-10,
                "Edges do not touch.");
            int ID = target.End.ID;
            target.End = source.End;
            target.End.ID = ID;

            if (target.IsBoundary && target.BoundaryEdgeNumber == targetBoundaryNumber)
            {
                source.Twin.Twin = target;
                target.Twin = source.Twin;
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

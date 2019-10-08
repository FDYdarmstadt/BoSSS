using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Recomposer
{
    class BoundaryEdgeHandler<T>
    {
        LinkedList<Edge<T>> candidates;

        int boundaryEdgeNumber; 

        public BoundaryEdgeHandler(int boundaryEdgeNumber)
        {
            candidates = new LinkedList<Edge<T>>();
            this.boundaryEdgeNumber = boundaryEdgeNumber;
        }

        public void AddCandidate(Edge<T> candidateForNewID)
        {
            candidates.AddFirst(candidateForNewID);
        }

        public void ConvertCandidatesToBoundary()
        {
            if (candidates.Count > 0)
            {
                FilterCandidates();
                if(candidates.Count > 0)
                {
                    CreateNewIDsForCandidates();
                }
            }
        }

        void FilterCandidates()
        {
            LinkedListNode<Edge<T>> current = candidates.First;
            while (current.Value.IsBoundary)
            {
                SetFirstCandidateToBoundary();
                current = candidates.First;
            }
            SetLastCandidateToBoundary();
        }

        void SetToBoundaryEdge(Edge<T> edge)
        {
            edge.IsBoundary = true;
            edge.BoundaryEdgeNumber = boundaryEdgeNumber;
        }

        void SetFirstCandidateToBoundary()
        {
            SetToBoundaryEdge(candidates.First.Value);
            candidates.RemoveFirst();
        }

        void SetLastCandidateToBoundary()
        {
            SetToBoundaryEdge(candidates.Last.Value);
            candidates.RemoveLast();
        }

        void CreateNewIDsForCandidates()
        {
            foreach (Edge<T> edge in candidates)
            {
                edge.End.ID = -100;
                edge.Start.ID = -50;
                SetToBoundaryEdge(edge);
            };
        }
    }
}

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class CornerCleaner
    {
        readonly ICollection<int> cornerCandidates;

        public CornerCleaner(int corners)
        {
            this.cornerCandidates = new List<int>(corners * 2);
        }

        public void RemoveAlreadyDealtWithCornerCellMergePairsFrom<T>(CellPairCollection<T>.EdgeCombo edge)
        {
            if (edge.Outer.Count > 0)
            {
                MeshCell<T> firstCandidate = edge.Outer[0].cell;
                if (edge.Inner.Count > 1)
                {
                    MeshCell<T> secondCandidate = edge.Outer[edge.Outer.Count - 1].cell;
                    if (cornerCandidates.Contains(secondCandidate.ID))
                    {
                        edge.Inner.RemoveAt(edge.Outer.Count - 1);
                        edge.Outer.RemoveAt(edge.Outer.Count - 1);
                    }
                    else
                    {
                        cornerCandidates.Add(secondCandidate.ID);
                    }
                }
                if (cornerCandidates.Contains(firstCandidate.ID))
                {
                    edge.Inner.RemoveAt(0);
                    edge.Outer.RemoveAt(0);
                }
                else
                {
                    cornerCandidates.Add(firstCandidate.ID);
                }
            }
        }
    }
}

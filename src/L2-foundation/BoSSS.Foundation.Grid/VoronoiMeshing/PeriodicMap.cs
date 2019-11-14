using BoSSS.Foundation.Grid.Classic;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    struct Corner : IEquatable<Corner>
    {
        public int FirstEdge;

        public int SecondEdge;

        public bool Equals(Corner other)
        {
            bool orderA = (other.FirstEdge == FirstEdge) && (other.SecondEdge == SecondEdge); 
            bool orderB = (other.SecondEdge == FirstEdge) && (other.SecondEdge == FirstEdge);
            return orderA || orderB;
        }

        public override bool Equals(object otherObject)
        {
            if(otherObject is Corner other)
            {
                bool orderA = (other.FirstEdge == FirstEdge) && (other.SecondEdge == SecondEdge);
                bool orderB = (other.SecondEdge == FirstEdge) && (other.SecondEdge == FirstEdge);
                return orderA || orderB;
            }
            else
            {
                return false;
            }
        }

        public override int GetHashCode()
        {
            return FirstEdge.GetHashCode();
        }
    }

    class PeriodicMap
    {
        public IDictionary<int, int> PeriodicBoundaryCorrelation;

        public IDictionary<int, Transformation> PeriodicBoundaryTransformations;

        public IDictionary<Corner, int> PeriodicCornerCorrelation;
    }
}

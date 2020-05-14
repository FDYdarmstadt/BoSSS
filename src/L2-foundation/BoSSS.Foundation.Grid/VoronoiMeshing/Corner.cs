using System;
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
            bool orderA = (other.FirstEdge == SecondEdge) && (other.SecondEdge == FirstEdge);
            bool orderB = (other.SecondEdge == SecondEdge) && (other.FirstEdge == FirstEdge);
            return orderA || orderB;
        }

        public override bool Equals(object otherObject)
        {
            if (otherObject is Corner other)
            {
                return Equals(other);
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
}

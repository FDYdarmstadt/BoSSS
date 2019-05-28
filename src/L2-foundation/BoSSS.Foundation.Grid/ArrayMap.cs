using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi
{
    public class ArrayMap
    {
        ArrayConnection[] map;

        public ArrayMap(ArrayConnection[] map)
        {
            this.map = map;
        }

        public ArrayConnection To(int j)
        {
            return map[j];
        }

        public ArrayMap Combine(ArrayMap map)
        {
            throw new NotImplementedException();
        }

        public ArrayMap Reverse(ArrayMap map, ArrayMap map2)
        {
            throw new NotImplementedException();
        }
    }

    public enum Connection { Created, Removed, Remained };

    public class ArrayConnection
    {
        public Connection Type;

        public int J;
    }

    public class MovingGrid
    {
        public ArrayMap Grid2PredecessorGrid;
        public VoronoiGrid Grid;
    }
}

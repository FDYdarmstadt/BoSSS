using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi
{
    public class ArrayMap
    {
        public ArrayConnection[] Map;

        public ArrayMap() { }

        public ArrayMap(ArrayConnection[] map)
        {
            this.Map = map;
        }

        public ArrayConnection To(int j)
        {
            return Map[j];
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
}

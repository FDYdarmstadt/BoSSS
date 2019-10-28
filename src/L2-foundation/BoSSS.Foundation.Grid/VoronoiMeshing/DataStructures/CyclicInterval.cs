using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class CyclicInterval
    {
        readonly int start;
        readonly int end;
        int index;

        public CyclicInterval(int start, int end, int index)
        {
            this.start = start;
            this.end = end;
            this.index = index;
            if (index < start || index > end)
            {
                throw new IndexOutOfRangeException();
            }
        }

        public int Current()
        {
            return index;
        }

        public void Next()
        {
            ++index;
            if (index > end)
            {
                index = start;
            }
        }

        public void Previous()
        {
            --index;
            if(index < start)
            {
                index = end;
            }
        }
    }
}

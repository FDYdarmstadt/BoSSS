using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryLineEnumerator : CountingEnumerator<BoundaryLine>
    {
        int Length;

        public BoundaryLineEnumerator(IEnumerator<BoundaryLine> enumerator, int Length) 
            : base(enumerator)
        {
            this.Length = Length;
        }

        public new int Counter {
            get {
                return base.Counter % Length;
            }
        }

        public CyclicInterval LineIndex {
            get {
                return new CyclicInterval(0, Length - 1, Counter);
            }
        }
    }

    
}

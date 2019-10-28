using BoSSS.Platform.LinAlg;
using System.Collections.Generic;

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

        public CyclicInterval LineIndex()
        {
            return new CyclicInterval(0, Length - 1, Counter);
        }

        public static BoundaryLineEnumerator GetEnumerator(
            Vector[] polygon)
        {
            BoundaryLine[] lines = BoundaryLine.ToLines(polygon);
            BoundaryLineEnumerator enumerateAndCount = GetEnumerator(lines);
            return enumerateAndCount;
        }

        public static BoundaryLineEnumerator GetEnumerator(BoundaryLine[] polygon)
        {
            ArrayEnumerator<BoundaryLine> lineEnum = new ArrayEnumerator<BoundaryLine>(polygon);
            BoundaryLineEnumerator enumerateAndCount = new BoundaryLineEnumerator(lineEnum, polygon.Length);
            return enumerateAndCount;
        }
    }
}

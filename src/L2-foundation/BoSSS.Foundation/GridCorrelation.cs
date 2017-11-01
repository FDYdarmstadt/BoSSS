using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Classic {

    /// <summary>
    /// Mapping of an old to a new grid under adaptive mesh refinement (see <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]})"/>.
    /// </summary>
    public class GridCorrelation {

        public long[] OldGlobalId;



        public AffineTrafo[] GeometricMapping;

    }
}

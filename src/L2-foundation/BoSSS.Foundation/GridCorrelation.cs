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

        /// <summary>
        /// Mapping from old cell index to old GlobalId.
        /// - index: local cell index into old grid
        /// - content: GlobalId of respective cell
        /// </summary>
        public long[] OldGlobalId;

        /// <summary>
        /// Mapping form old cell index to GlobalId of new cells.
        /// - 1st index: local cell index into old grid
        /// - 2nd index: enumeration
        /// - content: 
        /// </summary>
        public long[][] DestGlobalId;

        /// <summary>
        /// Affine transformation between old and new cells.
        /// - indexing: correlates with <see cref="DestGlobalId"/>
        /// - content: index into <see cref="GeometricMapping"/>; if negative, the mapping is the identity.
        /// </summary>
        public int[][] MappingIndex;

        /// <summary>
        /// Affine transformation from cell-local coordinate system to cell-local coordinate system.
        /// </summary>
        public AffineTrafo[] GeometricMapping;

    }
}

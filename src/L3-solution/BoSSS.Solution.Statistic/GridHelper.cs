using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Aggregation;

namespace BoSSS.Solution.Statistic
{
    static class GridHelper
    {
        /// <summary>
        /// Extracts a GridData object from IGridData, if there is one.
        /// </summary>
        public static GridData ExtractGridData(Foundation.Grid.IGridData gridData)
        {
            if (gridData is GridData classicGridData)
            {
                return classicGridData;
            }
            else if (gridData is AggregationGridData aggGridData)
            {
                GridData data = aggGridData.AncestorGrid;
                return aggGridData.AncestorGrid;
            }
            else
            {
                throw new NotImplementedException("ToDo: Grid Type not yet implemented");
            }
        }
    }
}

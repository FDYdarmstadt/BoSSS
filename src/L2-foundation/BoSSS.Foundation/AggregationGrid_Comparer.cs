using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Aggregation
{
    static class AggregationGridComparer
    {
        public static IEqualityComparer<IGrid> ReferenceComparer {
            get {
                return new GridComparer<AggregationGrid>(AreReferencesEqual);
            }
        }

        public static IEqualityComparer<IGrid> CellComparer {
            get {
                return new GridComparer<AggregationGrid>(AreCellsEqual);
            }
        }

        static bool AreCellsEqual(AggregationGrid A, AggregationGrid B)
        {
            //To do: Compare aggregation grid
            IEqualityComparer<IGrid> parentCellComparer = A.ParentGrid.GridSerializationHandler.CellComparer;
            bool parentCellsAreEqual = parentCellComparer.Equals(A.ParentGrid, B.ParentGrid);
            return parentCellsAreEqual;
        }

        static bool AreReferencesEqual(AggregationGrid A, AggregationGrid B)
        {
            //To do: Compare aggregation grid
            IEqualityComparer<IGrid> parentReferenceComparer = A.ParentGrid.GridSerializationHandler.ReferenceComparer;
            bool parentReferencesAreEqual = parentReferenceComparer.Equals(A.ParentGrid, B.ParentGrid);
            return parentReferencesAreEqual;
        }
    }
}

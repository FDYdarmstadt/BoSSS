using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json;

namespace BoSSS.Foundation.Grid.Aggregation
{
    public partial class AggregationGrid
    {
        /// <summary>
        /// Creates a human-readable string representation of the grid info object.
        /// </summary>
        /// <returns>A string summary of the grid.</returns>
        public override string ToString()
        {
            return "{ Guid = " + ID + "; Name = " + Name + "; Cell Count = "
                + NumberOfCells + "; Dim = " + SpatialDimension + " }";
        }

        [NonSerialized]
        [JsonIgnore]
        AggregationGridDatabaseMethods dataBaseMethods;

        public IGridSerializationHandler GridSerializationHandler {
            get {
                if (dataBaseMethods == null)
                    dataBaseMethods = new AggregationGridDatabaseMethods(this);
                return dataBaseMethods;
            }
        }
    }
}

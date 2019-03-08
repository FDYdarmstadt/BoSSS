using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;

namespace BoSSS.Foundation.Grid.Aggregation
{
    public partial class AggregationGrid : Grid.Classic.ISerializableGrid
    {
        public Type[] GetDataTypes()
        {
            throw new NotImplementedException();
        }

        public object[,] GetVectorData()
        {
            throw new NotImplementedException();
        }

        public bool HasEqualCells(IGrid grid)
        {
            throw new NotImplementedException();
        }

        public bool HasEqualReferences(IGrid grid)
        {
            throw new NotImplementedException();
        }

        public void SetVectorData(object[] data)
        {
            throw new NotImplementedException();
        }
    }
}

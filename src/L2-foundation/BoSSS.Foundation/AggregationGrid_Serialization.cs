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
        public object[][] GetVectorData()
        {
            throw new NotImplementedException();
        }

        public Guid[] GetVectorGuids()
        {
            throw new NotImplementedException();
        }

        public Type[] GetVectorTypes()
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

        public void SetVectorData(object[][] vectorDatas)
        {
            throw new NotImplementedException();
        }

        public void SetVectorGuids(Guid[] guids)
        {
            throw new NotImplementedException();
        }
    }
}

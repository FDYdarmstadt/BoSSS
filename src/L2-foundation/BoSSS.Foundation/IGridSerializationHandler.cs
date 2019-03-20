using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.IO;

namespace BoSSS.Foundation.Grid
{
    public interface IGridSerializationHandler : IVectorDataGrid, IComparableGrid
    {
        IDatabaseInfo Database {
            get;
            set;
        }
    }

    public interface IVectorDataGrid
    {
        Guid[] GetVectorGuids();
        void SetVectorGuids(Guid[] guids);
        object[][] GetVectorData();
        void SetVectorData(object[][] vectorDatas);
        Type[] GetVectorTypes();
        void Update();
    }

    public interface IComparableGrid
    {
        IEqualityComparer<IGrid> CellComparer { get; }
        IEqualityComparer<IGrid> ReferenceComparer { get; }
    }
}

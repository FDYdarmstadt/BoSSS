using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.IO;

namespace BoSSS.Foundation.Grid {
    public interface IGridSerializationHandler : IVectorDataGrid, IComparableGrid {
        IDatabaseInfo Database {
            get;
            set;
        }
    }

    public interface IVectorDataGrid {
        Guid[] GetVectorGuids();
        void SetVectorGuids(Guid[] guids);
        object[][] GetVectorData();
        void SetVectorData(object[][] vectorDatas);
        Type[] GetVectorTypes();
    }

    /// <summary>
    /// Grid comparison, divided into two parts:
    /// - basic properties: <see cref="BasePropertiesComparer"/>, which is cheap
    /// - individual cells: <see cref="CellComparer"/>, expensive to run
    /// </summary>
    public interface IComparableGrid {
        /// <summary>
        /// Comparison of all cells
        /// </summary>
        IEqualityComparer<IGrid> CellComparer { get; }

        /// <summary>
        /// Comparison of spatial dimension, cell types, number of cells, etc.
        /// </summary>
        IEqualityComparer<IGrid> BasePropertiesComparer { get; }
    }
}

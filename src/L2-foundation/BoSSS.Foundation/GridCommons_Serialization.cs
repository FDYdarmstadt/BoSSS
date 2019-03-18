using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json;
using BoSSS.Foundation.IO;

namespace BoSSS.Foundation.Grid.Classic
{
    class GridCommonsDatabaseMethods : IGridSerializationHandler
    {
        readonly GridCommons grid; 
        public GridCommonsDatabaseMethods(GridCommons grid)
        {
            this.grid = grid;
        }

        public IDatabaseInfo Database {
            get {
                return grid.Database;
            }
            set {
                grid.Database = value;
            }
        }

        [JsonIgnore]
        object[][] data;

        public object[][] GetVectorData()
        {
            if (data == null)
            {
                int numberOfObjects = grid.m_PredefinedGridPartitioning.Count + 2;
                data = new object[numberOfObjects][];
                data[0] = grid.Cells;
                for (int i = 0; i < grid.m_PredefinedGridPartitioning.Count; ++i)
                {
                    var s = grid.m_PredefinedGridPartitioning.ElementAt(i);
                    int[] cellToRankMap = s.Value.CellToRankMap;
                    data[i + 1] = cellToRankMap.Cast<object>().ToArray();
                }
                data[numberOfObjects - 1] = grid.BcCells;
            }
            return data;
        }

        public void SetVectorData(object[][] data)
        {
            grid.Cells = data[0].Cast<Cell>().ToArray();

            for (int i = 0; i < grid.m_PredefinedGridPartitioning.Count; ++i)
            {
                var s = grid.m_PredefinedGridPartitioning.ElementAt(i);
                grid.m_PredefinedGridPartitioning[s.Key] = new GridCommons.GridPartitioningVector { Guid = guids[i + 1], CellToRankMap = data[i + 1].Cast<int>().ToArray() };
            }
            if (data.Last() != null)
                grid.BcCells = data.Last().Cast<BCElement>().ToArray();
        }

        Type[] types;

        public Type[] GetVectorTypes()
        {
            if (types == null)
            {
                int numberOfObjects = grid.m_PredefinedGridPartitioning.Count + 2;
                types = new Type[numberOfObjects];
                types[0] = typeof(Cell);

                for (int i = 0; i < grid.m_PredefinedGridPartitioning.Count; ++i)
                {
                    types[i + 1] = typeof(int);
                }
                types[numberOfObjects - 1] = typeof(BCElement);
            }
            return types;
        }

        Guid[] guids;

        public Guid[] GetVectorGuids()
        {
            if (guids == null)
            {
                guids = InitializeVectorGuids();
            }
            return guids;
        }

        Guid[] InitializeVectorGuids()
        {
            int numberOfObjects = grid.m_PredefinedGridPartitioning.Count + 2;
            Guid[] guids = new Guid[numberOfObjects];
            guids[0] = grid.StorageGuid;
            for (int i = 0; i < grid.m_PredefinedGridPartitioning.Count; ++i)
            {
                var s = grid.m_PredefinedGridPartitioning.ElementAt(i);
                guids[i + 1] = s.Value.Guid;
            }
            guids[numberOfObjects - 1] = Guid.Empty;
            return guids;
        }

        public void SetVectorGuids(Guid[] guids)
        {
            this.guids = guids;
            grid.StorageGuid = guids[0];
            grid.BcCellsStorageGuid = guids.Last();
        }

        public void Initialize()
        {
            grid.InitNumberOfCells();
        }

        public IEqualityComparer<IGrid> ReferenceComparer {
            get {
                return GridCommonsComparer.ReferenceComparer;
            }
        }

        public IEqualityComparer<IGrid> CellComparer {
            get {
                return GridCommonsComparer.CellComparer;
            }
        }

    }
}

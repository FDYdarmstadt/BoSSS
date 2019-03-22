using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using Newtonsoft.Json;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.Grid.Aggregation
{
    class AggregationGridDatabaseMethods : IGridSerializationHandler
    {
        readonly AggregationGrid grid;
        readonly IGridSerializationHandler parentGridHandler;

        public AggregationGridDatabaseMethods(AggregationGrid grid)
        {
            parentGridHandler = grid.ParentGrid.GridSerializationHandler;
            this.grid = grid;
        }

        public IDatabaseInfo Database { get => grid.Database; set => grid.Database = value; }


        Guid[] guids; 

        public Guid[] GetVectorGuids()
        {
            if(guids == null)
            {
                guids = InitializeVectorGuids();
            }
            return guids;
        }

        Guid[] InitializeVectorGuids()
        {
            Guid[] parentGuids = parentGridHandler.GetVectorGuids();
            Guid[] aggregationGrids = new Guid[] { grid.AggCellStorageGuid };
            Guid[] guids = CombineArrays( aggregationGrids, parentGuids);
            return guids;
        }

        public void SetVectorGuids(Guid[] guids)
        {
            this.guids = guids;
            grid.AggCellStorageGuid = guids[0];
            parentGridHandler.SetVectorGuids(new Guid[] { guids[1], guids[2]});
        }

        public IEqualityComparer<IGrid> CellComparer => throw new NotImplementedException();

        public IEqualityComparer<IGrid> ReferenceComparer => throw new NotImplementedException();

        [JsonIgnore]
        object[][] data;

        public object[][] GetVectorData()
        {
            if (data == null)
            {
                object[][] aggCells = new object[][] { grid.AggregationCells.Cast<object>().ToArray() };
                object[][] parentData = parentGridHandler.GetVectorData();
                data = CombineArrays (aggCells, parentData);
            }
            return data;
        }

        public void SetVectorData(object[][] vectorDatas)
        {
            grid.AggregationCells = vectorDatas[0].Cast<AggregationGrid.AggCell>().ToArray();
            object[][] parentGridVectorData = Slice(vectorDatas, 1, vectorDatas.Length);
            parentGridHandler.SetVectorData(parentGridVectorData);
        }

        Type[] types;

        public Type[] GetVectorTypes()
        {
            if (types == null)
            {
                Type[] parentTypes = parentGridHandler.GetVectorTypes();
                types = CombineArrays( new[] { typeof(AggregationGrid.AggCell) }, parentTypes);
            }
            return types;
        }

        T[] CombineArrays<T>(T[] first, params T[] second)
        {
            T[] compositeArray = new T[first.Length + second.Length];
            FillIn(compositeArray, first);
            for(int i = 0 ; i < second.Length; ++i)
            {
                compositeArray[i + first.Length] = second[i];
            }
            return (compositeArray);
        }

        void FillIn<T>(T[] to, T[] from)
        {
            if(to. Length < from.Length)
            {
                throw new IndexOutOfRangeException();
            }
            for (int i = 0; i < from.Length; ++i)
            {
                to[i] = from[i];
            }
        }

        T[][] Slice<T>(T[][] arr, int from, int to)
        {
            int slicedLength = to - from;
            T[][] slicedArray = new T[slicedLength][];
            for (int i = 0; i < slicedLength; ++i)
            {
                slicedArray[i] = arr[i + from]; 
            }
            return slicedArray;
        }
    }
}

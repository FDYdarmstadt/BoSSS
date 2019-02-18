using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Newtonsoft.Json;
using BoSSS.Foundation.IO;
using ilPSP;

namespace BoSSS.Foundation.Grid.Classic
{
    interface ICompoundSerializable
    {
        object[,] GetData();
        void SetData(object[] data);
        Type[] GetDataTypes();
        IDatabaseInfo Database { get; set; }
    }

    public partial class GridCommons : ICompoundSerializable
    {
        [JsonIgnore]
        object[,] data;

        public object[,] GetData()
        {
            if (data == null)
            {
                int numberOfObjects = m_PredefinedGridPartitioning.Count + 2;
                data = new object[numberOfObjects,2];
                data[0, 1] = StorageGuid;
                data[0, 0] = Cells;
                for(int i = 0; i < m_PredefinedGridPartitioning.Count ; ++i)
                {
                    var s = m_PredefinedGridPartitioning.ElementAt(i);
                    int[] cellToRankMap = s.Value.CellToRankMap;
                    data[i + 1, 0] = cellToRankMap;
                }
                data[numberOfObjects - 1, 0] = BcCells;
                data[numberOfObjects - 1, 1] = BcCellsStorageGuid;
                
            }
            return data;
        }

        Type[] types;

        public Type[] GetDataTypes()
        {
            if(types == null)
            {
                int numberOfObjects = m_PredefinedGridPartitioning.Count + 2;
                types = new Type[numberOfObjects];
                types[0] = typeof(Cell[]);

                for (int i = 0; i < m_PredefinedGridPartitioning.Count; ++i)
                {
                    var s = m_PredefinedGridPartitioning.ElementAt(i);
                    types[i + 1] = typeof(KeyValuePair<string,GridPartitioningVector>);
                }
                types[numberOfObjects - 1] = typeof(BCElement[]);
            }
            return types;
        }

        public void SetData(object[] data)
        {
            Cells = (Cell[])data[0];
            for (int i = 0; i < m_PredefinedGridPartitioning.Count; ++i)
            {
                KeyValuePair<string, GridPartitioningVector> pair = (KeyValuePair<string, GridPartitioningVector>)data[i + 1];
                m_PredefinedGridPartitioning.Add(pair.Key, pair.Value);
            }
            BcCells = (BCElement[])data.Last();
        }
    }
}

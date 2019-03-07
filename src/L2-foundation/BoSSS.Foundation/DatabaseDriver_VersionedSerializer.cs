using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI.Wrappers;
using System.IO;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using BoSSS.Foundation.Grid;
using ilPSP;

namespace BoSSS.Foundation.IO
{
    class GridSerializer: IDisposable
    {
        VectorDataSerializer activeSerializer;
        readonly VectorDataSerializer[] serializer;

        public GridSerializer(params VectorDataSerializer[] dataSerializers) 
        {
            serializer = dataSerializers;
            activeSerializer = dataSerializers[0];
        }

        public T DeserializeGrid<T>(Guid gridGuid)
        {
            using (Stream s = GetGridStream(false, gridGuid))
            {
                T grid = activeSerializer.Deserialize<T>(s);
                return grid;
            }
        }

        public void SerializeGrid(IGrid grid)
        {
            using (Stream stream = GetGridStream(true, grid.ID))
            {
                activeSerializer.Serialize(stream, grid);
            }
        }

        public Guid SaveVector<T>(IList<T> vector)
        {
            return activeSerializer.SaveVector(vector);
        }

        public void SaveVector<T>(IList<T> vector, Guid id)
        {
            activeSerializer.SaveVector(vector, id);
        }

        public IList<T> LoadVector<T>(Guid id, ref Partitioning part)
        {
            return activeSerializer.LoadVector<T>(id, ref part);
        }

        public Stream GetGridStream(bool create,Guid gID)
        {
            return activeSerializer.GetGridStream(create, gID);
        }

        public void Dispose()
        {
            foreach (var singleSerializer in serializer)
            {
                singleSerializer.Dispose();
            }
        }
    }
}

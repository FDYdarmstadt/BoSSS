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
    class VersionedSerializer : IVectorDataSerializer
    {
        IVectorDataSerializer activeSerializer;
        readonly IVectorDataSerializer[] serializer;

        public VersionedSerializer(params IVectorDataSerializer[] dataSerializers) 
        {
            serializer = dataSerializers;
            activeSerializer = dataSerializers[0];
        }

        public T Deserialize<T>(Stream stream)
        {
            T grid = activeSerializer.Deserialize<T>(stream);
            return grid;
        }

        public void Serialize<T>(Stream stream, T obj)
        {
            activeSerializer.Serialize(stream, obj);
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

    }

}

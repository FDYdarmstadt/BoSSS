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
        IVectorDataSerializer preferedSerializer;
        IVectorDataSerializer standardSerializer;
        readonly Dictionary<string, IVectorDataSerializer> allSerializers;
        
        public VersionedSerializer(
            IVectorDataSerializer prefered, 
            IVectorDataSerializer standard, 
            params IVectorDataSerializer[] additionalSerializers) 
        {
            this.preferedSerializer = prefered;
            this.standardSerializer = standard;
            allSerializers = new Dictionary<string, IVectorDataSerializer>(additionalSerializers.Length + 2);
            AddToAllSerializers(prefered);
            AddToAllSerializers(standard);  
            foreach(IVectorDataSerializer additionalSerializer in additionalSerializers)
            {
                AddToAllSerializers(additionalSerializer);
            }
        }

        void AddToAllSerializers( IVectorDataSerializer serializer)
        {
            if (!allSerializers.ContainsKey(serializer.Name))
            {
                allSerializers.Add(serializer.Name, serializer);
            }
        }

        IVectorDataSerializer GetDeserializer(Stream stream)
        {
            IVectorDataSerializer dataSerializer;
            var reader = new StreamReader(stream); //No Using, stream will be disposed elsewhere
            {
                string firstLine = reader.ReadLine();
                bool found = allSerializers.TryGetValue(firstLine, out dataSerializer);
                if (!found)
                {
                    stream.Position = 0;
                    dataSerializer = standardSerializer;
                }
            }
            return dataSerializer;
        }

        IVectorDataSerializer GetSerializer(Stream stream)
        {
            var writer = new StreamWriter(stream);
            {
                writer.WriteLine(preferedSerializer.Name);
            }
            return preferedSerializer;
        }

        public T Deserialize<T>(Stream stream)
        {
            T grid = GetDeserializer(stream).Deserialize<T>(stream);
            return grid;
        }

        public void Serialize<T>(Stream stream, T obj)
        {
            GetSerializer(stream).Serialize(stream, obj);
        }

        public Guid SaveVector<T>(IList<T> vector)
        {
            return standardSerializer.SaveVector(vector);
        }

        public void SaveVector<T>(IList<T> vector, Guid id)
        {
            standardSerializer.SaveVector(vector, id);
        }

        public IList<T> LoadVector<T>(Guid id, ref Partitioning part)
        {
            return standardSerializer.LoadVector<T>(id, ref part);
        }

        public string Name => throw new Exception("Does not have a Name");

    }

}

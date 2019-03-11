using System;
using MPI.Wrappers;
using System.IO;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using BoSSS.Foundation.Grid;

namespace BoSSS.Foundation.IO
{

    abstract class Serializer : ISerializer
    {
        public abstract string Name { get; }

        protected abstract JsonSerializer JsonFormatter { get; }

        /// <summary>
        /// Indicates whether the content of the database should be serialized
        /// in a human-readable (debugging) format, or in a significantly
        /// smaller binary format
        /// </summary>
        protected static readonly bool DebugSerialization = false;

        public virtual object Deserialize(Stream stream, Type objectType)
        {
            using (var reader = GetJsonReader(stream))
            {
                object obj = JsonFormatter.Deserialize(reader, objectType);
                if (obj == null)
                {
                    throw new Exception("Deserializing failed.");
                }
                return obj;
            }
        }

        public void Serialize(Stream s, object obj, Type objectType)
        {
            using (var writer = GetJsonWriter(s))
            {
                JsonFormatter.Serialize(writer, obj, objectType);
                writer.Close();
                s.Close();
            }
        }

        protected JsonReader GetJsonReader(Stream s)
        {
            if (DebugSerialization)
            {
                return new JsonTextReader(new StreamReader(s));
            }
            else
            {
                return new BsonReader(s);
            }
        }

        protected JsonWriter GetJsonWriter(Stream s)
        {
            if (DebugSerialization)
            {
                return new JsonTextWriter(new StreamWriter(s));
            }
            else
            {
                return new BsonWriter(s);
            }
        }
    }

    class SerializerVersion0 : Serializer, ISerializer
    {
        public override string Name => "Version 0";

        protected override JsonSerializer JsonFormatter { get { return jsonFormatter; } }

        JsonSerializer jsonFormatter = new JsonSerializer()
        {
            NullValueHandling = NullValueHandling.Ignore,
            TypeNameHandling = TypeNameHandling.Auto,
            ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            ReferenceLoopHandling = ReferenceLoopHandling.Serialize,
            Binder = new MySerializationBinder()
        };

        class MySerializationBinder : Newtonsoft.Json.Serialization.DefaultSerializationBinder
        {

            public override Type BindToType(string assemblyName, string typeName)
            {
                if (assemblyName.Equals("BoSSS.Foundation") && typeName.Equals("BoSSS.Foundation.Grid.Cell[]"))
                {
                    typeName = "BoSSS.Foundation.Grid.Classic.Cell[]";
                }

                if (assemblyName.Equals("BoSSS.Foundation") && typeName.Equals("BoSSS.Foundation.Grid.BCElement[]"))
                {
                    typeName = "BoSSS.Foundation.Grid.Classic.BCElement[]";
                }
                Type T = base.BindToType(assemblyName, typeName);
                return T;
            }
        }

        public override object Deserialize(Stream stream, Type objectType)
        {
            if(typeof(IGrid).Equals( objectType))
            {
                objectType = typeof(Grid.Classic.GridCommons);
            }
            return JsonDeserialize(stream, objectType);
        }

        object JsonDeserialize(Stream stream, Type objectType)
        {
            using (var reader = GetJsonReader(stream))
            {
                object obj = jsonFormatter.Deserialize(reader, objectType);
                if (obj == null)
                {
                    throw new Exception("Deserializing failed.");
                }
                return obj;
            }
        }
    }

    class SerializerVersion1 : Serializer, ISerializer
    {
        public override string Name => "Version 1";

        protected override JsonSerializer JsonFormatter => jsonFormatter;

        JsonSerializer jsonFormatter = new JsonSerializer()
        {
            NullValueHandling = NullValueHandling.Ignore,
            TypeNameHandling = TypeNameHandling.Objects,
            ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            ReferenceLoopHandling = ReferenceLoopHandling.Serialize,
        };

    }

    public class MPIProcess
    {
        /// <summary>
        /// MPI rank of actual process within the MPI world communicator
        /// </summary>
        public int MyRank {
            get {
                int rank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                return rank;
            }
        }

        /// <summary>
        /// Number of MPI processes within the MPI world communicator
        /// </summary>
        public int Size {
            get {
                int size;
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
                return size;
            }
        }
    }
}

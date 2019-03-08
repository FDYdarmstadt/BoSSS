using System;
using MPI.Wrappers;
using System.IO;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;

namespace BoSSS.Foundation.IO
{
    abstract class Serializer : MPIProcess, ISerializer
    {
        public abstract string Name { get;}

        /// <summary>
        /// Indicates whether the content of the database should be serialized
        /// in a human-readable (debugging) format, or in a significantly
        /// smaller binary format
        /// </summary>
        protected static readonly bool DebugSerialization = false;

        protected abstract JsonSerializer JsonFormatter { get ;}

        public T Deserialize<T>(Stream stream)
        {
            using (var reader = GetJsonReader(stream))
            {
                T obj = JsonFormatter.Deserialize<T>(reader);
                if(obj == null)
                {
                    throw new Exception("Deserializing failed.");
                }
                return obj;
            }
        }
        
        public void Serialize<T>(Stream s, T obj )
        {
            using (var writer = GetJsonWriter(s))
            {
                JsonFormatter.Serialize(writer, obj);
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

using System;
using MPI.Wrappers;
using System.IO;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;

namespace BoSSS.Foundation.IO
{
    class Serializer : MPIProcess, IDisposable
    {
        public Serializer(IFileSystemDriver driver)
        {
            m_fsDriver = driver;
        }

        /// <summary>
        /// the file system driver
        /// </summary>
        public IFileSystemDriver FsDriver {
            get {
                return m_fsDriver;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        protected IFileSystemDriver m_fsDriver;

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

        /// <summary>
        /// Indicates whether the content of the database should be serialized
        /// in a human-readable (debugging) format, or in a significantly
        /// smaller binary format
        /// </summary>
        protected static readonly bool DebugSerialization = false;

        /// <summary>
        /// serialization formatter used for all bigger (data) objects
        /// </summary>
        JsonSerializer jsonFormatter = new JsonSerializer()
        {
            NullValueHandling = NullValueHandling.Ignore,
            TypeNameHandling = TypeNameHandling.Auto,
            ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            ReferenceLoopHandling = ReferenceLoopHandling.Serialize,
            Binder = new MySerializationBinder()
        };

        protected virtual JsonSerializer JsonFormatter {
            get { return jsonFormatter; }
        }

        public T Deserialize<T>(Stream stream)
        {
            using (var reader = GetJsonReader(stream))
            {
                return JsonFormatter.Deserialize<T>(reader);
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

        public Stream GetGridStream(bool create, Guid gID)
        {
            return FsDriver.GetGridStream(create, gID);
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

        /// <summary>
        /// finalize
        /// </summary>
        public void Dispose()
        {
            if (this.FsDriver is IDisposable)
            {
                ((IDisposable)this.FsDriver).Dispose();

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

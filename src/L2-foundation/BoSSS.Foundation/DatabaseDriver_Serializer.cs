﻿using System;
using System.IO;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using BoSSS.Foundation.Grid;
using System.Collections.Generic;
using System.Reflection;
using System.Linq;
using Newtonsoft.Json.Serialization;
using System.Diagnostics;
using System.IO.Compression;

namespace BoSSS.Foundation.IO
{
    interface ISerializer
    {
        string Name { get; }

        object Deserialize(Stream stream, Type objectType);

        void Serialize(Stream stream, object obj, Type objectType);
    }

    abstract class Serializer : ISerializer {
        public abstract string Name { get; }

        protected abstract JsonSerializer JsonFormatter { get; }

        /// <summary>
        /// Indicates whether the content of the database should be serialized
        /// in a human-readable (debugging) format, or in a significantly
        /// smaller binary format
        /// </summary>
        protected static readonly bool DebugSerialization = false;

        public virtual object Deserialize(Stream stream, Type objectType) {
            // fk, 19aug23: Create a MemoryStream and write the buffer into it
            // for compressed streams, the deserialization is than much faster than directly deserializing from the compressed stream
            using (MemoryStream memoryStream = new MemoryStream()) {
                stream.CopyTo(memoryStream);
                memoryStream.Seek(0, SeekOrigin.Begin);

                using (var reader = GetJsonReader(memoryStream)) {
                    object obj = JsonFormatter.Deserialize(reader, objectType);
                    if (obj == null) {
                        throw new IOException("Deserializing failed.");
                    }
                    return obj;
                }
            }
        }

        public void Serialize(Stream _s, object obj, Type objectType) {
            using (var s = new MemoryStream()) {
                using (var writer = GetJsonWriter(s)) {
                    JsonFormatter.Serialize(writer, obj, objectType);


                    s.Seek(0, SeekOrigin.Begin);
                    s.CopyTo(_s);
                    writer.Close();
                }
            }
        }

        protected JsonReader GetJsonReader(Stream s) {
            if (DebugSerialization) {
                return new JsonTextReader(new StreamReader(s));
            } else {
                return new BsonDataReader(s);
            }
        }

        protected JsonWriter GetJsonWriter(Stream s) {
            if (DebugSerialization) {
                return new JsonTextWriter(new StreamWriter(s));
            } else {
                return new BsonDataWriter(s);
            }
        }
    }

    class SerializerVersion0 : Serializer, ISerializer {
        public override string Name => "Version 0";

        protected override JsonSerializer JsonFormatter { get { return jsonFormatter; } }

        JsonSerializer jsonFormatter = new JsonSerializer() {
            NullValueHandling = NullValueHandling.Ignore,
            TypeNameHandling = TypeNameHandling.Auto,
            ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            ReferenceLoopHandling = ReferenceLoopHandling.Serialize,
            //Binder = new MySerializationBinder()
            SerializationBinder = new MySerializationBinder()
        };

        /// <summary>
        /// Serialization binder to ensure compatibility with older file versions, where certain classes
        /// were organized in a different namespace.
        /// </summary>
        class MySerializationBinder : Newtonsoft.Json.Serialization.DefaultSerializationBinder {

            public override Type BindToType(string assemblyName, string typeName) {
                if (assemblyName.Equals("BoSSS.Foundation") && typeName.Equals("BoSSS.Foundation.Grid.Cell[]")) {
                    typeName = "BoSSS.Foundation.Grid.Classic.Cell[]";
                }

                if (assemblyName.Equals("BoSSS.Foundation") && typeName.Equals("BoSSS.Foundation.Grid.BCElement[]")) {
                    typeName = "BoSSS.Foundation.Grid.Classic.BCElement[]";
                }
                Type T = base.BindToType(assemblyName, typeName);
                return T;
            }
        }

        public override object Deserialize(Stream stream, Type objectType) {
            if (typeof(IGrid).Equals(objectType)) {
                objectType = typeof(Grid.Classic.GridCommons);
            }
            return JsonDeserialize(stream, objectType);
        }

        object JsonDeserialize(Stream stream, Type objectType) {
            // fk, 19aug23: Create a MemoryStream and write the buffer into it
            // for compressed streams, the deserialization is than much faster than directly deserializing from the compressed stream
            using (MemoryStream memoryStream = new MemoryStream()) {
                stream.CopyTo(memoryStream);
                memoryStream.Seek(0, SeekOrigin.Begin);
           
                using (var reader = GetJsonReader(memoryStream)) {




                    object obj = jsonFormatter.Deserialize(reader, objectType);
                
                    if (obj == null) {
                        throw new IOException("Deserializing failed; object is null");
                    }
                    return obj;
                }
            }
        }
    }

    class SerializerVersion1 : Serializer, ISerializer {
        public override string Name => "Version 1";

        protected override JsonSerializer JsonFormatter => jsonFormatter;

        JsonSerializer jsonFormatter = new JsonSerializer() {
            NullValueHandling = NullValueHandling.Ignore,
            TypeNameHandling = TypeNameHandling.Objects,
            ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            ReferenceLoopHandling = ReferenceLoopHandling.Serialize,
            //Binder = new AdvancedLookupBinder()
        };
    }
}

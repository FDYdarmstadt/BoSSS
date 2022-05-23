using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Container for a list of <see cref="BatchProcessorClient"/> objects.
    /// </summary>
    [DataContract]
    public class BatchProcessorConfig {

        /// <summary>
        /// All items.
        /// </summary>
        [DataMember]
        public BatchProcessorClient[] AllQueues;

        /// <summary>
        /// Default queue, 
        /// </summary>
        [DataMember]
        public int DefaultQueueIndex = 0;


        static string ConfigFilePath {
            get {
                string settingsDir = Foundation.IO.Utils.GetBoSSSUserSettingsPath();
                string filePath = Path.Combine(settingsDir, "etc", "BatchProcessorConfig.json");
                return filePath;
            }
        }

        static string DefaultQueuesProjectOverridePath {
            get {
                string settingsDir = Foundation.IO.Utils.GetBoSSSUserSettingsPath();
                string filePath = Path.Combine(settingsDir, "etc", "DefaultQueuesProjectOverride.txt");
                return filePath;
            }
        }

        /// <summary>
        /// Loading of a configuration file from settings directory; if file does not exist, a default 
        /// configuration is returned and the file is created.
        /// </summary>
        static public BatchProcessorConfig LoadOrDefault() {

            string p = ConfigFilePath;
            if(File.Exists(p)) {
                string text = File.ReadAllText(p);
                var bpc = Deserialize(text);
                return bpc;
            } else {
                var r = new BatchProcessorConfig() {
                    AllQueues = new BatchProcessorClient[] {
                        new MiniBatchProcessorClient()
                    }
                };

                SaveConfiguration(r);

                return r;
            }
        }

        /// <summary>
        /// IO for the <see cref="DefaultQueuesProjectOverridePath"/>-file.
        /// </summary>
        static internal string GetDefaultBatchnameForProject(string projectName) {
            if(projectName == null)
                return null;

            try {
                if(!File.Exists(DefaultQueuesProjectOverridePath)) {
                    string[] lines = new string[]{
                    "## This file contains the default batch queue for specific BoSSS projects. ",
                    "## Precisely, this is an override for the global `DefaultQueueIndex` set in file `BatchProcessorConfig.json`.",
                    "## ",
                    "## The format for each line is",
                    "## ```",
                    "## Project-name : queue-name",
                    "## ```",
                    "## where `Project-name` is set in e Jupyter Notebook using `wmg.Init(\"LinslvPerf_ConstPoissonMpi1\")`",
                    "## and `queue-name` correlates with `BatchProcessorClient.Name`;" +
                    "## (Note: string comparisons are case-insensitive)"
                };

                    File.WriteAllLines(DefaultQueuesProjectOverridePath, lines);
                    return null;
                } else {
                    string[] lines = File.ReadAllLines(DefaultQueuesProjectOverridePath);
                    foreach(var l in lines) {
                        var _l = l.TrimStart();
                        if(_l.StartsWith("##"))
                            continue;

                        var l12 = _l.Split(":", StringSplitOptions.RemoveEmptyEntries);
                        if(l12.Length < 2)
                            continue;

                        if(l12[0].Trim().Equals(projectName, StringComparison.InvariantCultureIgnoreCase))
                            return l12[1].Trim();
                    }


                    return null;

                }
            } catch(Exception e) {
                Console.Error.WriteLine($"Error during 'DefaultQueuesProjectOverride.txt'-io: {e.GetType()}: {e.Message}");
                return null;
            }
        }


        /// <summary>
        /// Saves the given configuration in the default location
        /// </summary>
        /// <param name="bpc"></param>
        static public void SaveConfiguration(BatchProcessorConfig bpc) {
            string p = ConfigFilePath;
            string text = bpc.Serialize();
            File.WriteAllText(p, text);
        }



        /// <summary>
        /// Serializing this object into a string.
        /// </summary>
        public string Serialize() {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error,
                Formatting = Formatting.Indented
                //                ObjectCreationHandling = ObjectCreationHandling.
            };

            using(var tw = new StringWriter()) {
                //tw.WriteLine(this.GetType().AssemblyQualifiedName);
                using(JsonWriter writer = new JsonTextWriter(tw)) {  // Alternative: binary writer: BsonWriter
                    formatter.Serialize(writer, this);
                }

                string Ret = tw.ToString();
                return Ret;
            }

        }

        /// <summary>
        /// Used for control objects in work-flow management, de-serializing from a string.
        /// </summary>
        public static BatchProcessorConfig Deserialize(string Str) {
            var formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error
            };

            using(var tr = new StringReader(Str)) {
                //string typeName = tr.ReadLine();
                Type ControlObjectType = typeof(BatchProcessorConfig); //Type.GetType(typeName);
                using(JsonReader reader = new JsonTextReader(tr)) {
                    var obj = formatter.Deserialize(reader, ControlObjectType);

                    BatchProcessorConfig ctrl = (BatchProcessorConfig)obj;
                    return ctrl;
                }

            }
        }
        /*

        /// <summary>
        /// The default serialization binder used when resolving and loading classes from type names.
        /// </summary>
        public class DefaultSerializationBinder :
            SerializationBinder,
            Newtonsoft.Json.Serialization.ISerializationBinder {
            internal static readonly DefaultSerializationBinder Instance = new DefaultSerializationBinder();

          
            /// <summary>
            /// Initializes a new instance of the <see cref="DefaultSerializationBinder"/> class.
            /// </summary>
            public DefaultSerializationBinder() {
            }

            private Type GetTypeFromTypeNameKey(string? assemblyName, string typeName) {

                if(assemblyName != null) {
                    Assembly assembly;

                    // look, I don't like using obsolete methods as much as you do but this is the only way
                    // Assembly.Load won't check the GAC for a partial name
                    assembly = Assembly.LoadWithPartialName(assemblyName);

                    if(assembly == null) {
                        // will find assemblies loaded with Assembly.LoadFile outside of the main directory
                        Assembly[] loadedAssemblies = AppDomain.CurrentDomain.GetAssemblies();
                        foreach(Assembly a in loadedAssemblies) {
                            // check for both full name or partial name match
                            if(a.FullName == assemblyName || a.GetName().Name == assemblyName) {
                                assembly = a;
                                break;
                            }
                        }
                    }

                    if(assembly == null) {
                        throw new JsonSerializationException($"Could not load assembly '{assemblyName}'.");
                    }

                    Type? type = assembly.GetType(typeName);
                    if(type == null) {
                        // if generic type, try manually parsing the type arguments for the case of dynamically loaded assemblies
                        // example generic typeName format: System.Collections.Generic.Dictionary`2[[System.String, mscorlib, Version=2.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089],[System.String, mscorlib, Version=2.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089]]
                        if(typeName.IndexOf('`') >= 0) {
                            try {
                                type = GetGenericTypeFromTypeName(typeName, assembly);
                            } catch(Exception ex) {
                                throw new JsonSerializationException($"Could not find type '{typeName}' in assembly '{assembly.FullName}'.", ex);
                            }
                        }

                        if(type == null) {
                            throw new JsonSerializationException($"Could not find type '{typeName}' in assembly '{assembly.FullName}'.");
                        }
                    }

                    return type;
                } else {
                    return Type.GetType(typeName);
                }
            }

            private Type? GetGenericTypeFromTypeName(string typeName, Assembly assembly) {
                Type? type = null;
                int openBracketIndex = typeName.IndexOf('[');
                if(openBracketIndex >= 0) {
                    string genericTypeDefName = typeName.Substring(0, openBracketIndex);
                    Type genericTypeDef = assembly.GetType(genericTypeDefName);
                    if(genericTypeDef != null) {
                        List<Type> genericTypeArguments = new List<Type>();
                        int scope = 0;
                        int typeArgStartIndex = 0;
                        int endIndex = typeName.Length - 1;
                        for(int i = openBracketIndex + 1; i < endIndex; ++i) {
                            char current = typeName[i];
                            switch(current) {
                                case '[':
                                if(scope == 0) {
                                    typeArgStartIndex = i + 1;
                                }
                                ++scope;
                                break;
                                case ']':
                                --scope;
                                if(scope == 0) {
                                    string typeArgAssemblyQualifiedName = typeName.Substring(typeArgStartIndex, i - typeArgStartIndex);

                                    StructMultiKey<string?, string> typeNameKey = ReflectionUtils.SplitFullyQualifiedTypeName(typeArgAssemblyQualifiedName);
                                    genericTypeArguments.Add(GetTypeByName(typeNameKey));
                                }
                                break;
                            }
                        }

                        type = genericTypeDef.MakeGenericType(genericTypeArguments.ToArray());
                    }
                }

                return type;
            }

            private Type GetTypeByName(string? assemblyName, string typeName) {
                return _typeCache.Get(typeNameKey);
            }

            /// <summary>
            /// When overridden in a derived class, controls the binding of a serialized object to a type.
            /// </summary>
            /// <param name="assemblyName">Specifies the <see cref="Assembly"/> name of the serialized object.</param>
            /// <param name="typeName">Specifies the <see cref="System.Type"/> name of the serialized object.</param>
            /// <returns>
            /// The type of the object the formatter creates a new instance of.
            /// </returns>
            public override Type BindToType(string? assemblyName, string typeName) {
                return GetTypeByName(assemblyName, typeName);
            }

            
        }
        */
    }
}

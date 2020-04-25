using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP {
    

    /// <summary>
    /// Savind and loading for configuration files/objects
    /// </summary>
    static public class ConfigFileIOUtils {

         
        /// <summary>
        /// Loading of a configuration file from settings directory; if file does not exist, a default 
        /// configuration is returned and the file is created.
        /// </summary>
        static public T LoadOrDefault<T>(string ConfigFilePath) where T : new() {

            string p = ConfigFilePath;
            if(File.Exists(p)) {
                string text = File.ReadAllText(p);
                var bpc = Deserialize<T>(text);
                return bpc;
            } else {
                var r = new T();

                SaveConfiguration(r, ConfigFilePath);

                return r;
            }

        }

        /// <summary>
        /// Saves the given configuration object <paramref name="bpc"/> in the location <paramref name="ConfigFilePath"/>
        /// </summary>
        static public void SaveConfiguration<T>(T bpc, string ConfigFilePath) {
            string p = ConfigFilePath;
            string text = Serialize(bpc);
            File.WriteAllText(p, text);
        }

        ///// <summary>
        ///// Serializing this object into a string.
        ///// </summary>
        //public string Serialize() {
        //    return Serialize(this);
        //}


        /// <summary>
        /// Serializing object <paramref name="dis"/> into a string.
        /// </summary>
        public static string Serialize<T>(T dis) {
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
                    formatter.Serialize(writer, dis);
                }

                string Ret = tw.ToString();
                return Ret;
            }
            
        }

        /// <summary>
        /// De-serializing from a JSON string.
        /// </summary>
        public static T Deserialize<T>(string Str) {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error
            };

            
            using(var tr = new StringReader(Str)) {
                //string typeName = tr.ReadLine();
                Type ControlObjectType = typeof(T); //Type.GetType(typeName);
                using(JsonReader reader = new JsonTextReader(tr)) {
                    var obj = formatter.Deserialize(reader, ControlObjectType);

                    T ctrl = (T)obj;
                    return ctrl;
                }
              
            }
        }
    }
}

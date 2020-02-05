using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// 
    /// </summary>
    public class BatchProcessorConfig {

        public BatchProcessorClient.Config[] AllQueus;


        static string ConfigFilePath {
            get {
                string settingsDir = Foundation.IO.Utils.GetBoSSSUserSettingsPath();
                string filePath = Path.Combine(settingsDir, "etc", "BatchProcessorConfig.json");
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
                    AllQueus = new BatchProcessorClient.Config[] {
                        new MiniBatchProcessorClient.Config()
                    }
                };

                SaveConfiguration(r);

                return r;
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
        /// Used for control objects in work-flow management, 
        /// </summary>
        public static BatchProcessorConfig Deserialize(string Str) {
            JsonSerializer formatter = new JsonSerializer() {
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





    }
}

using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution {

    /// <summary>
    /// Configuration of <see cref="ilPSP.Connectors.Matlab.BatchmodeConnector"/>
    /// </summary>
    [DataContract]
    public class MatlabConnectorConfig {
        /// <summary>
        /// <see cref="ilPSP.Connectors.Matlab.BatchmodeConnector.MatlabExecuteable"/>
        /// </summary>
        [DataMember]
        public string MatlabExecuteable;

        /// <summary>
        /// <see cref="ilPSP.Connectors.Matlab.BatchmodeConnector.Flav"/>
        /// </summary>
        [DataMember]
        public string Flav = ilPSP.Connectors.Matlab.BatchmodeConnector.Flavor.Matlab.ToString();


        /// <summary>
        /// Loads configuration from default location in user directory
        /// </summary>
        public static MatlabConnectorConfig LoadDefault(string userDir) {
            string ConfigFile = Path.Combine(userDir, "etc", "MatlabConnectorConfig.json");
            if (!File.Exists(ConfigFile)) {
                var r = new MatlabConnectorConfig();
                r.SaveDefault(userDir);
                return r;
            }
            string str = File.ReadAllText(ConfigFile);

            return Deserialize(str);
        }

        /// <summary>
        /// Saves configuration in default location in user directory
        /// </summary>
        public void SaveDefault(string userDir) {
            string ConfigFile = Path.Combine(userDir, "etc", "MatlabConnectorConfig.json");
            var Str = this.Serialize();
            File.WriteAllText(ConfigFile, Str);
        }


        /// <summary>
        /// JSON deserialization
        /// </summary>
        public static MatlabConnectorConfig Deserialize(string Str) {


            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error
            };


            using (var tr = new StringReader(Str)) {
                //string typeName = tr.ReadLine();
                Type ControlObjectType = typeof(MatlabConnectorConfig); //Type.GetType(typeName);
                using (JsonReader reader = new JsonTextReader(tr)) {
                    var obj = formatter.Deserialize(reader, ControlObjectType);

                    MatlabConnectorConfig ctrl = (MatlabConnectorConfig)obj;
                    return ctrl;
                }

            }
        }

        /// <summary>
        /// JSON serialization
        /// </summary>
        public string Serialize() {
            JsonSerializer formatter = new JsonSerializer() {
                NullValueHandling = NullValueHandling.Ignore,
                TypeNameHandling = TypeNameHandling.Auto,
                ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
                ReferenceLoopHandling = ReferenceLoopHandling.Error,
                Formatting = Newtonsoft.Json.Formatting.Indented
            };

            using (var tw = new StringWriter()) {
                //tw.WriteLine(this.GetType().AssemblyQualifiedName);
                using (JsonWriter writer = new JsonTextWriter(tw)) {  // Alternative: binary writer: BsonWriter
                    formatter.Serialize(writer, this);
                }

                string Ret = tw.ToString();
                return Ret;
            }
        }

    }
}

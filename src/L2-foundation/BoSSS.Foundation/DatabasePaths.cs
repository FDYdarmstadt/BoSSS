using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ilPSP;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// "magical" file `DatabasePaths.json` in the current working directory, (resp, under `~/.BoSSS/etc/DatabasePaths.json`)
    /// </summary>
    public class DatabasePaths :ConfigFileBase {
        public DatabasePaths() {
            PathsOfDatabases = new Dictionary< Guid, ISet<string> > ();
        }


        /// <summary>
        /// - key: database IF, <see cref="IDatabaseInfo.ID"/>, <see cref="DatabaseInfoJson.ID"/>
        /// . value: paths to access this database, <see cref="IDatabaseInfo.AlternateDbPaths"/>, <see cref="DatabaseInfoJson.Paths"/>
        /// </summary>
        public IDictionary<Guid, ISet<string>> PathsOfDatabases;


        /// <summary>
        /// default filename
        /// </summary>
        public static readonly string filename = "DatabasePaths.json";

        protected override string GetDefaultFilePath() {
            return filename; 
        }

        /// <summary>
        /// tries to load from working and user directory
        /// </summary>
        static public IDictionary<Guid, ISet<string>> LoadAll() {
            var ret = new Dictionary<Guid, ISet<string>>();

            string[] paths = [filename, Path.Combine(Utils.GetBoSSSUserSettingsPath(), "etc", filename)];
            foreach(string path in paths) {
                if(File.Exists(path)) {
                    try {
                        var txt = File.ReadAllText(path);
                        var dbp = ConfigFileIOUtils.Deserialize<DatabasePaths>(txt);

                        foreach(var id in dbp.PathsOfDatabases.Keys) {

                            if(!ret.TryGetValue(id, out ISet<string> d)) {
                                d = new HashSet<string>();
                                ret.Add(id, d);

                            }

                            d.AddRange(dbp.PathsOfDatabases[id]);
                        }


                    } catch (Exception e) {
                        Console.Error.WriteLine($"Loading file `{path}`: {e}");
                    }
                }
            }

            return ret;
        }



    }
}

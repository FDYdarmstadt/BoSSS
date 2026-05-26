using ilPSP;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace BoSSS.Foundation.IO {


    /// <summary>
    /// serialization of the `DatabaseInfo.json` in the root directory of a BoSSS da´tabase.
    /// </summary>
    [Serializable]
    public class DatabaseInfoJson :ConfigFileBase {

        /// <summary>
        /// ctor
        /// </summary>
        public DatabaseInfoJson() {
            ID = Guid.Empty;
            Paths = new HashSet<string>();
        }

        /// <summary>
        /// ctor
        /// </summary>
        public DatabaseInfoJson(string dbPath) {
            CurrentPath = dbPath;
            ID = Guid.NewGuid();
            Paths = new HashSet<string>();
            Paths.Add(CurrentPath);

            try {
                if(base.Load()) {
                    // file already existed 

                    // add current path if necessary
                    if(!Paths.Contains(CurrentPath)) {
                        Paths.Add(CurrentPath);
                        base.Save();
                    }


                } else {
                    // file was freshly created

                }
            } catch(Exception e) {
                Console.WriteLine("IO Error database file: " + e);
            }

        }

        /// <summary>
        /// database ID
        /// </summary>
        public Guid ID;

        /// <summary>
        /// paths to reach the database, e.g.,
        /// different mount-paths on different computers for network file sastems
        /// </summary>
        public ISet<string> Paths;

        /// <summary>
        /// filename 
        /// </summary>
        public static readonly string Filename = "DatabaseInfo.json";


        [NonSerialized]
        string CurrentPath;

        /// <summary>
        /// fielname for IO
        /// </summary>
        protected override string GetDefaultFilePath() {
            return Path.Combine(CurrentPath, Filename);
        }
    }



}

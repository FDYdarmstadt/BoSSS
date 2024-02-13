/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Standard implementation of a database info object.
    /// </summary>
    public class DatabaseInfo : IDatabaseInfo {

        static List<DatabaseInfo> DatabaseInfos;

        static object padlock_DatabaseInfos = new object();

        /// <summary>
        /// Tries to open a database if <paramref name="path"/> is existent;
        /// Otherwise, creates a new database and opens it.
        /// </summary>
        public static DatabaseInfo CreateOrOpen(string path) {


            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int MpiRank);
            if(MpiRank == 0) {
                
                if(!File.Exists(path) && !Directory.Exists(path)) {
                    DirectoryInfo targetDirectory = new DirectoryInfo(path);
                    if(!targetDirectory.Exists) {
                        targetDirectory.Create();
                    } else {
                        if(targetDirectory.GetFiles().Length > 0)
                            throw new ArgumentException("Must be empty.");
                        if(targetDirectory.GetDirectories().Length > 0)
                            throw new ArgumentException("Must be empty.");
                    }

                    // Create structure
                    Directory.CreateDirectory(System.IO.Path.Combine(targetDirectory.FullName, "data"));
                    Directory.CreateDirectory(System.IO.Path.Combine(targetDirectory.FullName, "timesteps"));
                    Directory.CreateDirectory(System.IO.Path.Combine(targetDirectory.FullName, "grids"));
                    Directory.CreateDirectory(System.IO.Path.Combine(targetDirectory.FullName, "sessions"));
                }
                
            }
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
            Thread.Sleep(5000); // allow Network file systems to catch up
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

            return Open(path);
        }



        /// <summary>
        /// Opens an existing 
        /// </summary>
        public static DatabaseInfo Open(string _path) {
            return Open(_path, null);
        }

        /// <summary>
        /// Open the database
        /// </summary>
        public static DatabaseInfo Open((string DbPath, string MachineFilter)[] AlternateDbPaths) {
            return Open(null, AlternateDbPaths);
        }

        /// <summary>
        /// One of the stupid hacks that we have to do due to the shitty convoluted design of the database.
        /// </summary>
        public static void Close(IDatabaseInfo _dbi) {
            DatabaseInfo dbi = _dbi as DatabaseInfo;
            if(dbi == null) {
                //Console.Error.WriteLine($"Reminder: some strange {typeof(IDatabaseInfo).Name} implementation seems to be around here - unable to close.");
                return;
            }
            
            dbi.Controller.DBDriver.Dispose();
            lock(padlock_DatabaseInfos) {
                if(DatabaseInfos == null)
                    DatabaseInfos = new List<DatabaseInfo>();

                for(int i = 0; i < DatabaseInfos.Count; i++) {
                    if(object.ReferenceEquals(dbi, DatabaseInfos[i])) {
                        DatabaseInfos.RemoveAt(i);
                        i--;
                    }
                }
            }
        }



        /// <summary>
        /// Open the database
        /// </summary>
        public static DatabaseInfo Open(string _path, (string DbPath, string MachineFilter)[] AlternateDbPaths = null) {
            if(_path == null && (AlternateDbPaths == null || AlternateDbPaths.Length <= 0))
                throw new ArgumentException("hey, buddy, you must provide some path to open");
            
            List<ValueTuple<string, string>> allPaths = new List<(string, string)>();

            allPaths.Add((_path, null));
            if(AlternateDbPaths != null)
                allPaths.AddRange(AlternateDbPaths);

            string mName = System.Environment.MachineName.ToLowerInvariant();

            string dbPath = null;
            foreach(var t in allPaths) {
                string __path = t.Item1;
                string filter = t.Item2;

                if(!filter.IsNullOrEmpty() && !filter.IsEmptyOrWhite()) {
                    if(!mName.Contains(filter)) {
                        continue;
                    }
                }

                if(Directory.Exists(__path) || File.Exists(__path)) { // the latter is for ZIP-file databases
                    dbPath = __path;
                    break;
                }

            }

            if(dbPath == null) {
                Console.Error.WriteLine("Unable to open database: ");
                Console.Error.WriteLine($"primary path: {_path}");
                if (AlternateDbPaths == null || AlternateDbPaths.Length <= 0) {
                    Console.Error.WriteLine("No alternative paths specified.");
                } else {
                    Console.Error.WriteLine("Got " + AlternateDbPaths.Length + " all other paths: ");
                    for(int i = 0; i < AlternateDbPaths.Length; i++) {
                        Console.Error.WriteLine($"  #{i}: {AlternateDbPaths[i].DbPath}, filter is {AlternateDbPaths[i].MachineFilter}");
                    }

                }

                throw new IOException("Unable to open database - all given paths either don't exist or are ruled out by the machine filter on this machine with name " + mName + " .");
            }


            lock(padlock_DatabaseInfos) {
                if(DatabaseInfos == null)
                    DatabaseInfos = new List<DatabaseInfo>();



                for(int i = 0; i < DatabaseInfos.Count; i++) {
                    var dbi = DatabaseInfos[i];
                    IFileSystemDriver dbdr = dbi.Controller.DBDriver.FsDriver;
                    if(dbdr is StandardFsDriver fsdr ) {
                        if(fsdr.IsDisposed) {
                            dbi.Controller.DBDriver.Dispose();
                        }
                        DatabaseInfos.RemoveAt(i);
                        i--;
                    }

                }

                foreach(var db in DatabaseInfos) {
                    if(db.PathMatch(dbPath))
                        return db;
                }

                var newDb = new DatabaseInfo(dbPath);
                DatabaseInfos.Add(newDb);
                return newDb;
            }
        }



        /// <summary>
        /// Stores the path
        /// </summary>
        /// <param name="path">Path to the database</param>
        private DatabaseInfo(string path) {
            this.Path = path;
            if(path == null) {
                Controller = NullDatabaseController.Instance;
            } else {
                Controller = new DatabaseController(this);
            }
        }

        /// <summary>
        /// Full path to the base directory of the database.
        /// </summary>
        public string Path {
            get;
            private set;
        }

        

        /// <summary>
        /// Provides functionality to copy/move/delete info objects stored in
        /// the database
        /// </summary>
        public IDatabaseController Controller {
            get;
            private set;
        }

        /// <summary>
        /// Returns a string representation of this database.
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            // TO DO: THIS REQUIRES IO OPERATIONS => BAD!
            //return "{ Session Count = " + Controller.Sessions.Count()
            //    + "; Grid Count = " + Controller.Grids.Count()
            //    + "; Path = " + Path + " }";

            // Temporary workaround
            return "{ Session Count = " + ((DatabaseController)Controller).SessionCount
                + "; Grid Count = " + ((DatabaseController)Controller).GridCount
                + "; Path = " + Path + " }";
        }

        /// <summary>
        /// detects if some other path actually also points to this database
        /// </summary>
        public bool PathMatch(string otherPath) {
            if(this.Path == otherPath)
                return true;
            if(!Directory.Exists(otherPath))
                return false;

            string TokenName = Guid.NewGuid().ToString() + ".token";

            string file1 = System.IO.Path.Combine(this.Path, TokenName);
            File.WriteAllText(file1, "this is a test file which can be safely deleted.");

            string file2 = System.IO.Path.Combine(otherPath, TokenName);

            return File.Exists(file2);
        }

        /// <summary>
        /// The sessions of this database.
        /// </summary>
        public IList<ISessionInfo> Sessions {
            get {
                using (new FuncTrace()) {
                    //Stopwatch stw = new Stopwatch();

                    //Console.WriteLine("aquire sessions...");
                    //stw.Reset();
                    //stw.Start();
                    var allsessions = Controller.Sessions;
                    //stw.Stop();
                    //Console.WriteLine("done. " + stw.ElapsedMilliseconds);

                    //Console.WriteLine("sorting sessions...");
                    //stw.Reset();
                    //stw.Start();
                    var R = allsessions.OrderByDescending(s => s.WriteTime).ToList();
                    //stw.Stop();
                    //Console.WriteLine("done. " + stw.ElapsedMilliseconds);

                    return R;
                }
            }
        }

        /// <summary>
        /// The grids of this database.
        /// </summary>
        public IList<IGridInfo> Grids {
            get {
                return Controller.Grids.OrderByDescending(g => g.WriteTime).ToList();
            }
        }

        /// <summary>
        /// Sessions sorted according to projects, see <see cref="ISessionInfo.ProjectName"/>.
        /// </summary>
        public IDictionary<string, IEnumerable<ISessionInfo>> Projects {
            get {
                Dictionary<string, IEnumerable<ISessionInfo>> R = new Dictionary<string, IEnumerable<ISessionInfo>>();

                foreach(var s in this.Sessions) {
                    string PrjNmn = s.ProjectName;
                    if(PrjNmn == null || PrjNmn.Length <= 0)
                        PrjNmn = "__unknown_project__";

                    IEnumerable<ISessionInfo> ProjectSessions;
                    if(!R.TryGetValue(PrjNmn, out ProjectSessions)) {
                        ProjectSessions = new List<ISessionInfo>();
                        R.Add(PrjNmn, ProjectSessions);
                    }
                    Debug.Assert(ProjectSessions != null);

                    ((List<ISessionInfo>)ProjectSessions).Add(s);
                }


                return R;
            }
        }

         /// <summary>
        /// Reference equality
        /// </summary>
        public bool Equals(IDatabaseInfo other) {
            if(other == null)
                return false;
            if(object.ReferenceEquals(this, other))
                return true;

            string mName = System.Environment.MachineName.ToLowerInvariant();

            List<ValueTuple<string, string>> allPaths = new List<(string, string)>();
            allPaths.Add((other.Path, null));
            allPaths.AddRange(other.AlternateDbPaths);

            foreach(var t in allPaths) {
                string path = t.Item1;
                string filter = t.Item2;

                if(!filter.IsNullOrEmpty() && !filter.IsEmptyOrWhite()) {
                    if(!mName.Contains(filter)) {
                        continue;
                    }
                }

                if(this.PathMatch(path))
                    return true;
            }

            return false;
        }

        /// <summary>
        /// 
        /// </summary>
        public override bool Equals(object obj) {
            if (object.ReferenceEquals(obj, this))
                return true;

            return this.Equals(obj as DatabaseInfo);
        }

        /// <summary>
        /// 
        /// </summary>
        public override int GetHashCode() {
            return 1; // deactivate hashing
        }


        /// <summary>
        /// Alternative paths to access the database, if the main path is not present on a given machine.
        /// This allows to use the same control file or object on different machines, where the database is located in a different path.
        /// - 1st entry: path into the local file system
        /// - 2nd entry: optional machine name filter
        /// </summary>
        public (string DbPath, string MachineFilter)[] AlternateDbPaths {
            get {
                return ReadAlternateDbPaths(this.Path);
            }    
        }

        static (string DbPath, string MachineFilter)[] ReadAlternateDbPaths(string dbPath) {
            if (dbPath == null)
                throw new ArgumentNullException();
            string p = System.IO.Path.Combine(dbPath, "AlternatePaths.txt");

            if (!File.Exists(p))
                return new ValueTuple<string, string>[0];

            string[] lines = File.ReadAllLines(p);

            var ret = new List<ValueTuple<string, string>>();
            foreach (var line in lines) {
                if (line.StartsWith(";;"))
                    continue;
                string[] parts = line.Split(new[] { "," }, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length >= 2) {
                    ret.Add((parts[0], parts[1]));
                } else if (parts.Length >= 1) {
                    ret.Add((parts[0], null));
                }
            }

            return ret.ToArray();
        }

        /// <summary>
        /// Adds a new <see cref="IDatabaseInfo.AlternateDbPaths"/> to the database configuration
        /// </summary>
        /// <return>
        /// - true: the path was added
        /// - false: the alternate path is already registered for the database and was not added again.
        /// </return>
        public static bool AddAlternateDbPaths(string dbPath, string altPath, string machFilter) {
            if (altPath == null)
                throw new ArgumentNullException();

            var PathsSoFar = ReadAlternateDbPaths(dbPath);
            foreach(var tt in PathsSoFar) {
                if (tt.DbPath == altPath && tt.MachineFilter == machFilter)
                    return false;
            }

            string p = System.IO.Path.Combine(dbPath, "AlternatePaths.txt");

            string LeadText = "";
            if (File.Exists(p))
                LeadText = File.ReadAllText(p);

            using(var stw = new StreamWriter(p)) {
                if(LeadText.Length > 0) {
                    stw.WriteLine(LeadText);
                }

                if(machFilter == null) {
                    stw.WriteLine(altPath + ",");
                } else {
                    stw.WriteLine(altPath + "," + machFilter);
                }
            }

            return true;
        }


    }
}

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
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
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
        /// <param name="DbPaths">
        /// multiple directories where to search for the datebease;
        /// this is intended for use in HPC cluster environments, where the database directory might have different mount paths on different computers.
        /// </param>
        public static DatabaseInfo Open(string[] DbPaths) {
            return Open(null, DbPaths);
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
        /// <param name="AlternateDbPaths">
        /// alternative paths, where to search for a BoSSS database if it is not found at the primary path <paramref name="_path"/>;
        /// this is intended for use in HPC cluster environments, where the database directory might have different mount paths on different computers.
        /// </param>
        /// <param name="_path">
        /// default path to load from
        /// </param>
        public static DatabaseInfo Open(string _path, IEnumerable<string> AlternateDbPaths = null) {
            if(_path == null && (AlternateDbPaths == null || AlternateDbPaths.Count() <= 0))
                throw new ArgumentException("hey, buddy, you must provide some path to open");

            HashSet<string> allPaths = new HashSet<string>() { _path };
            allPaths.AddRange(AlternateDbPaths);


            string dbPath = null;
            foreach(var __path in allPaths) {
                if(Directory.Exists(__path) || File.Exists(__path)) { // the latter is for ZIP-file databases
                    dbPath = __path;
                    break;
                }
            }

            if(dbPath == null) {
                using(var sw = new StringWriter()) {
                    sw.WriteLine("Unable to open database: ");
                    sw.WriteLine($"primary path: {_path}");
                    sw.WriteLine(" Alernative paths: " + AlternateDbPaths.ToConcatString());

                    throw new IOException(sw.ToString());
                }
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

        private DatabaseInfoJson databaseInfoJson;



        /// <summary>
        /// Stores the path
        /// </summary>
        /// <param name="path">Path to the database</param>
        private DatabaseInfo(string path) {
            this.Path = path;
            if(path == null) {
                Controller = NullDatabaseController.Instance;
                databaseInfoJson = new DatabaseInfoJson();
            } else {
                Controller = new DatabaseController(this);
                databaseInfoJson = new DatabaseInfoJson(path);
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
            /*
            if(this.Path == otherPath)
                return true;
            if(!Directory.Exists(otherPath))
                return false;

            string TokenName = Guid.NewGuid().ToString() + ".token";

            string file1 = System.IO.Path.Combine(this.Path, TokenName);
            File.WriteAllText(file1, "this is a test file which can be safely deleted.");

            string file2 = System.IO.Path.Combine(otherPath, TokenName);

            return File.Exists(file2);
            */

            if(this.Path == otherPath)
                return true;
            if(otherPath.IsEmptyOrWhite())
                return false;
            if(!Directory.Exists(otherPath))
                return false;

            if(!DatabaseUtils.IsValidBoSSSDatabase(otherPath)) 
                return false;

            var o = new DatabaseInfoJson(otherPath);
            return o.ID == this.ID;

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
        /// true, id the <see cref="IDatabaseInfo.ID"/> of <paramref name="other"/> matches with this database id
        /// </summary>
        public bool Equals(IDatabaseInfo other) {
            if(other == null)
                return false;
            if(object.ReferenceEquals(this, other))
                return true;

            string mName = System.Environment.MachineName.ToLowerInvariant();

            List<string> allPaths = new List<string>();
            allPaths.Add(other.Path);
            allPaths.AddRange(other.AlternateDbPaths);

            foreach(var path in allPaths) {
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
            return this.ID.GetHashCode();
        }


        /// <summary>
        /// Alternative paths to access the database, if the main path is not present on a given machine.
        /// This allows to use the same control file or object on different machines, where the database is located in a different path.
        /// </summary>
        public string[] AlternateDbPaths {
            get {
                var r = new List<string>();
                r.AddRange(databaseInfoJson.Paths);
                r.Remove(this.Path);
                return r.ToArray();
            }
        }


        public Guid ID {
            get {
                return databaseInfoJson.ID;
            }
        }

        /*
        static (string DbPath, string MachineFilter)[] ReadAlternateDbPaths(string dbPath) {
            if (dbPath == null)
                throw new ArgumentNullException();
            string p = System.IO.Path.Combine(dbPath, "AlternatePaths.txt");

            return ReadAlternateDbPathsFile(p);
        }
        */

        /*
        /// <summary>
        /// Reads the line-oriented alternate database path format used by database-local
        /// <c>AlternatePaths.txt</c> and job-local <see cref="AlternateDatabasePathsFileName"/>.
        /// </summary>
        public static (string DbPath, string MachineFilter)[] ReadAlternateDbPathsFile(string path) {
            if(path == null)
                throw new ArgumentNullException();
            if(!File.Exists(path))
                return new ValueTuple<string, string>[0];

            return ParseAlternateDbPathLines(File.ReadAllLines(path));
        }

        /*
        /// <summary>
        /// Parses alternate database path lines. Blank lines are ignored, and lines starting
        /// with <c>;;</c> after leading whitespace are treated as comments.
        /// </summary>
        public static (string DbPath, string MachineFilter)[] ParseAlternateDbPathLines(IEnumerable<string> lines) {
            if(lines == null)
                throw new ArgumentNullException();

            var ret = new List<ValueTuple<string, string>>();
            foreach(var line in lines) {
                if(line == null)
                    continue;
                if(line.IsEmptyOrWhite())
                    continue;
                if(line.TrimStart().StartsWith(";;"))
                    continue;

                string[] parts = line.Split(new[] { "," }, 2, StringSplitOptions.None);
                string dbPath = parts[0].Trim();
                if(dbPath.IsEmptyOrWhite())
                    continue;
                string machineFilter = null;
                if(parts.Length >= 2) {
                    machineFilter = parts[1].Trim();
                    if(machineFilter.IsEmptyOrWhite())
                        machineFilter = null;
                }
                ret.Add((dbPath, machineFilter));
            }

            return ret.ToArray();
        }

        static (string DbPath, string MachineFilter)[] ReadDeployedAlternateDatabasePaths() {
            string path = System.IO.Path.Combine(Directory.GetCurrentDirectory(), AlternateDatabasePathsFileName);
            return ReadAlternateDbPathsFile(path);
        }
        */

        static void AddAlternateDbPathsUnique(List<ValueTuple<string, string>> paths, IEnumerable<ValueTuple<string, string>> additionalPaths) {
            foreach(var p in additionalPaths) {
                if(p.Item1.IsEmptyOrWhite())
                    continue;
                if(paths.Any(q => q.Item1 == p.Item1 && q.Item2 == p.Item2))
                    continue;
                paths.Add(p);
            }
        }

        /// <summary>
        /// Adds a new <see cref="IDatabaseInfo.AlternateDbPaths"/> to the database configuration
        /// </summary>
        /// <return>
        /// - true: the path was added
        /// - false: the alternate path is already registered for the database and was not added again.
        /// </return>
        public bool AddAlternateDbPath(string altPath) {
            if (altPath.IsEmptyOrWhite())
                throw new ArgumentNullException();

            /*
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
            */

            if(!databaseInfoJson.Paths.Contains(altPath)) {
                databaseInfoJson.Paths.Add(altPath);
                databaseInfoJson.Save();
                return true;
            } else {
                return false;
            }
        }

        
    }
}

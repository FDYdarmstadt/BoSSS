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

using BoSSS.Foundation.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP.Utils;
using System.Data;
using System.Reflection;
using System.Threading;
using ilPSP;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using ilPSP.Tracing;
using System.Collections.ObjectModel;

namespace BoSSS.Application.BoSSSpad {




    /// <summary>
    /// Workflow management.
    /// </summary>
    public partial class WorkflowMgm {

        /// <summary>
        /// Not intended for user interaction.
        /// </summary>
        internal WorkflowMgm() {
            SetEqualityBasedSessionJobControlCorrelation();
        }

        string m_CurrentProject;

        /// <summary>
        /// Name of the current project;
        /// </summary>
        public string CurrentProject {
            get {
                return m_CurrentProject;
            }
        }

        /// <summary>
        /// Clears/Invalidates all cached data.
        /// </summary>
        public void InvalidateCaches() {
            m_Sessions = null;
            m_Grids = null;
        }


        TimeSpan m_UpdatePeriod = new TimeSpan(0, 5, 0);

        /// <summary>
        /// Data like <see cref="Sessions"/>, <see cref="Grids"/> is cached for performance reasons; after this time span is elapsed, the data is re-read from disk.
        /// </summary>
        public TimeSpan UpdatePeriod {
            get {
                return m_UpdatePeriod;
            }
            set {
                m_UpdatePeriod = value;
            }
        }

        internal static Stopwatch getKeys = new Stopwatch();

        /// <summary>
        /// Correlation of session, job and control object is done by name
        /// </summary>
        public void SetNameBasedSessionJobControlCorrelation() {
            SessionInfoJobCorrelation = delegate (ISessionInfo sinf, Job job) {
                //using (var tr = new FuncTrace("NameBasedSessionInfoJobControlCorrelation")) {
                try {
                    var kaq = new Dictionary<string, object>(sinf.KeysAndQueries); // caching for IO acceleration
                    
                    // compare project name
                    if (!kaq.ContainsKey(BoSSS.Solution.Application.PROJECTNAME_KEY))
                        return false;

                    if (!Convert.ToString(kaq[BoSSS.Solution.Application.PROJECTNAME_KEY]).Equals(this.CurrentProject))
                        return false;

                    // compare session name
                    if (!kaq.ContainsKey(BoSSS.Solution.Application.SESSIONNAME_KEY))
                        return false;

                    if (!Convert.ToString(kaq[BoSSS.Solution.Application.SESSIONNAME_KEY]).Equals(job.Name))
                        return false;

                    // fall tests passed
                    return true;
                } catch (Exception) {
                    return false;
                }
                //}

            };
            SessionInfoAppControlCorrelation = delegate (ISessionInfo sinf, AppControl ctrl) {
                try {

                    var kaq = new Dictionary<string, object>(sinf.KeysAndQueries); // for IO acceleration

                    // compare project name
                    if (!kaq.ContainsKey(BoSSS.Solution.Application.PROJECTNAME_KEY))
                        return false;

                    if(!Convert.ToString(kaq[BoSSS.Solution.Application.PROJECTNAME_KEY]).Equals(ctrl.ProjectName))
                        return false;

                    // compare session name
                    if(!kaq.ContainsKey(BoSSS.Solution.Application.SESSIONNAME_KEY))
                        return false;

                    if(!Convert.ToString(kaq[BoSSS.Solution.Application.SESSIONNAME_KEY]).Equals(ctrl.SessionName))
                        return false;

                    // fall tests passed
                    return true;

                } catch(Exception) {
                    return false;
                }
            };

            JobAppControlCorrelation = delegate (Job j, AppControl c) {
                Debug.Assert(j != null);
                var jc = j.GetControl();
                if(jc == null) {
                    return false;
                }
                if(c == null)
                    return false;
                return (jc.SessionName == c.SessionName) && (jc.ProjectName == c.ProjectName);
            };
        }


        /// <summary>
        /// Correlation of session, job and control object is done <see cref="AppControl.Equals(object)"/>
        /// </summary>
        public void SetEqualityBasedSessionJobControlCorrelation() {
            SessionInfoJobCorrelation = delegate (ISessionInfo sinf, Job job) {
                var c_job = job.GetControl();
                
                try {
                    var c_sinf = sinf.GetControl();
                    if(c_sinf != null)
                        return c_sinf.Equals(c_job);
                    else
                        return false;
                } catch(Exception) {
                    return false;
                }
            };
            SessionInfoAppControlCorrelation = delegate (ISessionInfo sinf, AppControl ctrl) {
                try {
                    var c_sinf = sinf.GetControl();
                    if(c_sinf != null)
                        return c_sinf.Equals(ctrl);
                    else
                        return false;
                } catch(Exception) {
                    return false;
                }
            };

            JobAppControlCorrelation = delegate (Job j, AppControl c) {
                Debug.Assert(j != null);
                var jc = j.GetControl();
                if(jc == null) {
                    return false;
                }
                return jc.Equals(c);
            };
        }


        /// <summary>
        /// Defines, globally for the entire workflow management, how session in the project correlate to jobs.
        /// </summary>
        public Func<ISessionInfo, Job, bool> SessionInfoJobCorrelation;

        /// <summary>
        /// Defines, globally for the entire workflow management, how session in the project correlate to control objects.
        /// </summary>
        public Func<ISessionInfo, AppControl, bool> SessionInfoAppControlCorrelation;

        /// <summary>
        /// Defines, globally for the entire workflow management, how jobs in the project correlate to control objects.
        /// </summary>
        public Func<Job, AppControl, bool> JobAppControlCorrelation;

        internal bool RunWorkflowFromBackup {
            get {
                using (var tr = new FuncTrace("RunWorkflowFromBackup")) {
                    //string runfromBackup = System.Environment.GetEnvironmentVariable("BOSSS_RUNTESTFROMBACKUP");
                    //return runfromBackup.IsEmptyOrWhite() == false;
                    const string magicFile = "BOSSS_RUNTESTFROMBACKUP.txt";
                    bool exists = File.Exists(magicFile);
                    tr.Info($"File {magicFile} exists? {exists} (current directory: {System.Environment.CurrentDirectory})");
                    return exists;
                }
            }
        }


        /// <summary>
        /// Defines the name of the current project; also creates a default database
        /// </summary>
        public void Init(string ProjectName, BatchProcessorClient ExecutionQueue = null) {
            if ((m_CurrentProject == null) || (!m_CurrentProject.Equals(ProjectName)))
                InvalidateCaches();
            m_CurrentProject = ProjectName;
            Console.WriteLine("Project name is set to '{0}'.", ProjectName);

            if (ExecutionQueue == null) {
                ExecutionQueue = BoSSSshell.GetDefaultQueue();
                Console.WriteLine("Default Execution queue is chosen for the database.");
            }

            //if(InteractiveShell.ExecutionQueues.Any(Q => Q is MiniBatchProcessorClient))
            //    MiniBatchProcessor.Server.StartIfNotRunning();

            if (RunWorkflowFromBackup == false) {
                try {
                    DefaultDatabase = ExecutionQueue.CreateOrOpenCompatibleDatabase(ProjectName);
                } catch (Exception e) {
                    Console.Error.WriteLine($"{e.GetType().Name} caught during creation/opening of default database: {e.Message}.");
                }
            } else {
                Console.WriteLine("trying to run from backup database...");
                var pp = ExecutionQueue.AllowedDatabasesPaths[0];
                var baseDir = new DirectoryInfo(pp.LocalMountPath);

                if (!Path.IsPathRooted(pp.LocalMountPath))
                    throw new IOException($"Illegal entry for `AllowedDatabasesPaths` for {this.ToString()}: only absolute/rooted paths are allowed, but {pp.LocalMountPath} is not.");

                var bkupDbs = baseDir.GetDirectories("bkup*." + ProjectName);
                Console.WriteLine("   Bkup Database dirs: " + bkupDbs.ToConcatString("", ", ", ";"));

                if (bkupDbs.Length <= 0) {
                    Console.Error.WriteLine("No Backups found; unable to run worksheet from backup database.");
                    Console.WriteLine("Trying to create/open default database.");

                    try {
                        DefaultDatabase = ExecutionQueue.CreateOrOpenCompatibleDatabase(ProjectName);
                    } catch (Exception e) {
                        Console.Error.WriteLine($"{e.GetType().Name} caught during creation/opening of default database: {e.Message}.");
                    }

                } else {

                    var dbDir = bkupDbs.OrderBy(dir => dir.CreationTime).Last(); // select newest available backup
                    Console.WriteLine("Selecting newest available backup: " + dbDir);


                    IDatabaseInfo dbi = BoSSSshell.OpenDatabase(dbDir.FullName);

                    if (!pp.PathAtRemote.IsEmptyOrWhite()) {
                        string fullPathAtRemote = pp.PathAtRemote.TrimEnd('/', '\\');
                        string remoteDirSep = pp.PathAtRemote.Contains('/') ? "/" : "\\";
                        fullPathAtRemote = fullPathAtRemote + remoteDirSep + dbDir;
                        DatabaseInfo.AddAlternateDbPaths(dbDir.FullName, fullPathAtRemote, null);
                    }

                    DefaultDatabase = dbi;
                }
            }

        }

        /// <summary>
        /// Removes any pre-existing results regarding this project; use with extreme care!
        /// </summary>
        public void ResetProject(bool ResetJobs = true, bool deleteDeployments = false, bool deleteSessions = false, bool deleteGrids = false) {
            if (deleteSessions) {
                Console.WriteLine("Deleting Sessions in projects...");
                foreach(var s in this.Sessions) {
                    s.Delete(true);
                }
            } else {
                Console.WriteLine("Not deleting any sessions, because not specified (`deleteSessions:false`).");
            }

            if (deleteGrids) {
                Console.WriteLine("Deleting Grids in projects...");
                foreach (var g in this.Grids) {
                    g.Delete(true);
                }
            } else {
                Console.WriteLine("Not deleting any grids, because not specified (`deleteGrids:false`).");
            }

            if (deleteDeployments) {
                Console.WriteLine("Deleting Job deployments in projects...");
                foreach (var j in this.AllJobs) {
                    j.Value.DeleteOldDeploymentsAndSessions(deleteSessions);
                }

            } else {
                Console.WriteLine("Not deleting any grids, because not specified (`deleteDeployments:false`).");
            }

            if (ResetJobs) {
                Console.WriteLine("Forgetting all Jobs defined in this notebook so far...");
                this.m_AllJobs.Clear();
            } else {
                Console.WriteLine("Job objects defined so far remain valid (`ResetJobs:true`). ");
            }
        }



        IDatabaseInfo m_DefaultDatabase;

        /// <summary>
        /// primary database to store objects use in this project.
        /// </summary>
        public IDatabaseInfo DefaultDatabase {
            get {
                return m_DefaultDatabase;
            }
            set {
                if (CurrentProject.IsEmptyOrWhite()) {
                    throw new NotSupportedException("Workflow management not initialized yet - call Init(...)!");
                }

                m_DefaultDatabase = value;
            }
        }

        DateTime m_AllDatabases_CacheTime = DateTime.Now;
        List<IDatabaseInfo> m_AllDatabases;


        /// <summary>
        /// all databases which match the <see cref="CurrentProject"/> name
        /// </summary>
        public IReadOnlyList<IDatabaseInfo> AllDatabases {
            get {
                if (RunWorkflowFromBackup == false) {
                    if (m_AllDatabases == null || ((DateTime.Now - m_AllDatabases_CacheTime) > UpdatePeriod)) {

                        var allDBs = new List<IDatabaseInfo>();

                        foreach (var q in BoSSSshell.ExecutionQueues) {
                            int cnt = 0;
                            foreach (var dbPath in q.AllowedDatabasesPaths) {
                                IDatabaseInfo dbi = null;
                                string db_path = Path.Combine(dbPath.LocalMountPath, this.CurrentProject);
                                if (cnt == 0) {
                                    try {
                                        dbi = BoSSSshell.OpenOrCreateDatabase(db_path);
                                    } catch (Exception e) {
                                        Console.Error.WriteLine($"{e.GetType().Name} caught during creation/opening of database: {e.Message}.");
                                    }
                                } else {
                                    try {
                                        dbi = BoSSSshell.OpenDatabase(db_path);
                                    } catch (Exception) {
                                        dbi = null;
                                    }
                                }
                                if (dbi != null)
                                    allDBs.Add(dbi);
                                cnt++;
                            }
                        }

                        if (m_DefaultDatabase != null) {
                            int DefaultDbMatch = -1;
                            int cnt = 0;
                            foreach (var dbi in allDBs) {
                                if (Object.ReferenceEquals(m_DefaultDatabase, dbi) || dbi.PathMatch(m_DefaultDatabase.Path)) {
                                    DefaultDbMatch = cnt;
                                    break;
                                }
                                cnt++;
                            }

                            if (DefaultDbMatch >= 0) {
                                allDBs.RemoveAt(DefaultDbMatch);
                                allDBs.Insert(0, m_DefaultDatabase);
                            } else {
                                allDBs.Insert(0, m_DefaultDatabase);
                            }
                        }

                        m_AllDatabases_CacheTime = DateTime.Now;
                        m_AllDatabases = allDBs;
                    }

                    return m_AllDatabases.AsReadOnly();
                } else {
                    return (new List<IDatabaseInfo>(new[] { DefaultDatabase })).AsReadOnly();
                }
            }
        }





        
        ISessionInfo[] m_Sessions;
        List<SessionAtomic> m_Atomics = new List<SessionAtomic>();

        /// <summary>
        /// Clears the cache for <see cref="Sessions"/> and enforces to re-read the database.
        /// </summary>
        public void ResetSessionsCache() {
            m_Sessions = null;
        }

        /// <summary>
        /// Prevents changes (through IO) in <see cref="Sessions"/> as long the atomic is in use;   with the exception of <see cref="ResetSessionsCache"/>
        /// This also helps IO performance;
        /// </summary>
        public class SessionAtomic : IDisposable {
            
            internal SessionAtomic(WorkflowMgm owner) {
                m_owner = owner;
            }

            internal readonly WorkflowMgm m_owner;

            /// <summary>
            /// releases the atomic
            /// </summary>
            public void Dispose() {
                m_owner.m_Atomics.Remove(this);
                if (m_owner.m_Atomics.Count >= 0)
                    m_owner.ResetSessionsCache();
            }
        }

        /// <summary>
        /// <see cref="SessionAtomic"/>
        /// </summary>
        public SessionAtomic EnterSessionAtomic() {
            m_Atomics.Add(new SessionAtomic(this));
            return m_Atomics.Last();
        }


        /// <summary>
        /// A list of all sessions in the current project.
        /// </summary>
        public ISessionInfo[] Sessions {
            get {
                Debugger.Launch();
                using (var tr = new FuncTrace()) {
                    if (CurrentProject.IsEmptyOrWhite()) {
                        Console.WriteLine("Workflow management not initialized yet - call Init(...)!");
                        return new ISessionInfo[0];
                    }

                    if(m_Atomics.Count > 0 && m_Sessions != null) {
                        return m_Sessions;
                    } else { 

                    
                        List<ISessionInfo> ret = new List<ISessionInfo>();

                        if (BoSSSshell.databases != null) {
                            foreach (var db in BoSSSshell.databases) {
                                var SS = db.Sessions.Where(delegate (ISessionInfo si) {
                                    //#if DEBUG 
                                    //                                return si.ProjectName.Equals(this.CurrentProject);
                                    //#else
                                    Guid g = Guid.Empty;
                                    try {
                                        g = si.ID;

                                        if(si.ProjectName == null) {
                                            if(si is SessionProxy sip) {
                                                sip.TriggerReload = true;
                                            }
                                        }

                                        return si.ProjectName.Equals(this.CurrentProject);
                                    } catch (Exception e) {
                                        string sessionString;
                                        try {
                                            sessionString =  (si?.GetType().ToString() ?? "X") + " // " + (si?.ToString() ?? "NULL");
                                        } catch(Exception e2) {
                                            sessionString = e2.ToString();
                                        }

                                        Console.WriteLine("Warning: " + e.Message + " reading session " + g + ".");
                                        tr.Warning(" reading session " + g + ": " + e + " (" + e.StackTrace + "); Session = " + sessionString);
                                        return false;
                                    }
                                    //#endif
                                });
                                ret.AddRange(SS);
                            }
                        }

                        if (m_Atomics.Count > 0)
                            m_Sessions = ret.ToArray();
                        else
                            m_Sessions = null;
                        

                        return ret.ToArray();
                    }

                }
            }
        }





        DataTable m_Projects;
        /// <summary>
        /// A list of all available Projects
        /// </summary>
        public DataTable Projects {
            get {


                m_Projects = new DataTable("Projects");
                // Create Project DataTable
                {
                    DataColumn column;

                    // Create new DataColumn, set DataType,
                    // ColumnName and add to DataTable.
                    column = new DataColumn();
                    column.DataType = typeof(string);
                    column.ColumnName = "Name";
                    column.ReadOnly = true;
                    column.Unique = true;
                    // Add the Column to the DataColumnCollection.
                    m_Projects.Columns.Add(column);

                    // Create second column.
                    column = new DataColumn();
                    column.DataType = typeof(int);
                    column.ColumnName = "SessionCount";
                    column.ReadOnly = false;
                    column.Unique = false;
                    // Add the column to the table.
                    m_Projects.Columns.Add(column);

                    // Create second column.
                    column = new DataColumn();
                    column.DataType = typeof(List<Guid>);
                    column.ColumnName = "SessionIds";
                    column.ReadOnly = false;
                    column.Unique = false;
                    // Add the column to the table.
                    m_Projects.Columns.Add(column);
                }

                // Fill with values
                { 
                    if (BoSSSshell.databases != null) {
                        foreach (var db in BoSSSshell.databases) {
                            foreach(var si in db.Sessions) {
                                DataRow[] foundProject = m_Projects.Select("Name = '"+si.ProjectName+"'");
                                if (foundProject.Length != 0) {
                                    foundProject[0]["SessionCount"] = (int)foundProject[0]["SessionCount"] + 1;
                                    ((List<Guid>)foundProject[0]["SessionIds"]).Add(si.ID);
                                } else {
                                    DataRow newProject = m_Projects.NewRow();
                                    newProject["Name"] = si.ProjectName;
                                    newProject["SessionCount"] = 1;
                                    newProject["SessionIds"] = new List<Guid> { si.ID };
                                    m_Projects.Rows.Add(newProject);
                                }
                            }                            
                        }
                    }

                }

                return m_Projects;
            }
        }        

        /// <summary>
        /// A list of all tags in all sessions.
        /// </summary>
        public string[] Tags {
            get {
                if (CurrentProject.IsEmptyOrWhite()) {
                    Console.WriteLine("Workflow management not initialized yet - call Init(...)!");
                    return new string[0];
                }

                HashSet<string> r = new HashSet<string>();

                foreach (var s in this.Sessions) {
                    r.AddRange(s.Tags);
                }

                return r.ToArray();
            }
        }

        IGridInfo[] m_Grids;
        DateTime m_Grids_CacheTime;

        /// <summary>
        /// A list of all grids which are used in the current project.
        /// </summary>
        public IGridInfo[] Grids {
            get {
                if (CurrentProject.IsEmptyOrWhite()) {
                    Console.WriteLine("Workflow management not initialized yet - call Init(...)!");
                    return new IGridInfo[0];
                }

                if (m_Grids == null || ((DateTime.Now - m_Grids_CacheTime) > UpdatePeriod)) {
                    HashSet<IGridInfo> grids = new HashSet<IGridInfo>(
                    new ilPSP.FuncEqualityComparer<IGridInfo>((g1, g2) => g1.ID.Equals(g2.ID), g => g.ID.GetHashCode()));

                    foreach(var dbi in this.AllDatabases) {
                        foreach(var g in dbi.Grids) {
                            grids.Add(g);
                        }
                    }


                    foreach (var s in this.Sessions) {
                        //Console.Write("Session " + s.ID + " ... ");
                        grids.AddRange(s.GetGrids());
                        //Console.WriteLine(" done.");
                    }

                    m_Grids = grids.ToArray();
                    m_Grids_CacheTime = DateTime.Now;
                }

                return m_Grids;
            }
        }

        /// <summary>
        /// <see cref="IDatabaseInfoExtensions.SaveGrid{TG}(IDatabaseInfo, ref TG, bool)"/>
        /// </summary>
        /// <param name="filename"></param>
        /// <param name="EdgeTagFunc">
        /// <see cref="IGrid_Extensions.DefineEdgeTags(IGrid, Func{double[], string})"/>
        /// </param>
        /// <returns></returns>
        public GridCommons ImportGrid(string filename, Func<double[],string> EdgeTagFunc = null) {
            //using(var md5 = System.Security.Cryptography.MD5.Create()) {
            //    using(var stream = File.OpenRead(filename)) {
            //        var Hasch = md5.ComputeHash(stream);
            //    }
            //}
            if(this.DefaultDatabase == null) {
                throw new NotImplementedException("Default database for project not set yet.");
            }

            GridCommons r = Solution.GridImport.GridImporter.Import(filename);

            if(EdgeTagFunc != null)
                r.DefineEdgeTags(EdgeTagFunc);

            this.DefaultDatabase.SaveGrid(ref r, force:false);
            m_Grids = null; // trigger re-read;
            return r;
        }

        /// <summary>
        /// Saves the grid <paramref name="g"/> in the <see cref="DefaultDatabase"/>;
        /// see <see cref="IDatabaseInfoExtensions.SaveGrid{TG}(IDatabaseInfo, ref TG, bool)"/>
        /// </summary>
        public GridCommons SaveGrid(GridCommons g) {
            this.DefaultDatabase.SaveGrid(ref g, force: false);
            m_Grids = null; // trigger re-read;
            return g;
        }
        

        /// <summary>
        /// The keys and queries <see cref="ISessionInfo.KeysAndQueries"/> of all sessions in the 
        /// project (see <see cref="Sessions"/>) in one table.
        /// </summary>
        public DataTable SessionTable {
            get {
                var adiColi = AdditionalSessionTableColums.Select(kv => new Tuple<string, Func<ISessionInfo, object>>(kv.Key, kv.Value)).ToArray();
                return this.Sessions.GetSessionTable(adiColi);
            }
        }

        Dictionary<string, Func<ISessionInfo, object>> m_AdditionalSessionTableColums = new Dictionary<string, Func<ISessionInfo, object>>();

        /// <summary>
        /// Custom, user-defined columns for the session table (<see cref="SessionTable"/>).
        /// - keys: column name
        /// - values: functions which map the session info to a column value.
        /// </summary>
        public IDictionary<string, Func<ISessionInfo, object>> AdditionalSessionTableColums {
            get {
                return m_AdditionalSessionTableColums;
            }
        }

        



        internal Dictionary<string, Job> m_AllJobs = new Dictionary<string, Job>();

        /// <summary>
        /// Lists all compute jobs which are currently known by the work flow management system.
        /// - key: job name
        /// - item 
        /// </summary>
        public ReadOnlyDictionary<string, Job> AllJobs {
            get {
                return new ReadOnlyDictionary<string, Job>(m_AllJobs);
            }
        }


        /// <summary>
        /// Blocks until all jobs in <see cref="AllJobs"/> are either <see cref="JobStatus.FailedOrCanceled"/>
        /// or <see cref="JobStatus.FinishedSuccessful"/>.
        /// </summary>
        /// <param name="TimeOutSeconds">
        /// If positive, this method should terminate at latest after approximately this time period.
        /// </param>
        /// <param name="PollingIntervallSeconds">
        /// Seconds to wait before checking the jobs status again; should be in the order of seconds, not to overload the IO.
        /// </param>
        public void BlockUntilAllJobsTerminate(double TimeOutSeconds = -1, double PollingIntervallSeconds = 10) {
            using(var tr = new FuncTrace()) {
                DateTime start = DateTime.Now;
                while(true) {

                    //if(InteractiveShell.ExecutionQueues.Any(Q => Q is MiniBatchProcessorClient))
                    //    MiniBatchProcessor.Server.StartIfNotRunning(false); // hack for parallel execution of tests

                    Thread.Sleep((int)PollingIntervallSeconds);

                    if(TimeOutSeconds > 0) {
                        double RuntimeSoFar = (DateTime.Now - start).TotalSeconds;
                        if(RuntimeSoFar > TimeOutSeconds) {
                            Console.WriteLine("Timeout.");
                            return;
                        }
                    }

                    // dbg_launch();

                    bool terminate = true;
                    foreach(var J in this.AllJobs) {
                        var s = J.Value.Status;
                        tr.Info("Testing job: " + J);
                        if(s != JobStatus.FailedOrCanceled && s != JobStatus.FinishedSuccessful && s != JobStatus.PreActivation && s != JobStatus.Unknown) {
                            tr.Info("not terminating because of job: " + J);
                            terminate = false;
                            break;
                        }
                    }

                    if(terminate) {
                        Console.WriteLine("All jobs finished.");
                        m_Sessions = null;
                        return;
                    }

                    Thread.Sleep((int)(1000.0 * PollingIntervallSeconds));
                }
            }
        }
        List<Tuple<AppControl, int>> RegisteredControls = new List<Tuple<AppControl, int>>();


        /// <summary>
        /// Blocks until any running or queued job in <see cref="AllJobs"/> reaches 
        /// either <see cref="JobStatus.FailedOrCanceled"/>
        /// or <see cref="JobStatus.FinishedSuccessful"/>.
        /// </summary>
        /// <param name="TimeOutSeconds">
        /// If positive, this method should terminate at latest after approximately this time period.
        /// </param>
        /// <param name="PollingIntervallSeconds">
        /// Seconds to wait before checking the jobs status again; should be in the order of seconds, not to overload the IO.
        /// </param>
        /// <param name="JustFinished">
        /// On exit, the recently finished job
        /// </param>
        /// <returns>
        /// Remaining number of queued/running jobs plus 1
        /// </returns>
        public int BlockUntilAnyJobTerminate(out Job JustFinished,  double TimeOutSeconds = -1, double PollingIntervallSeconds = 10) {
            DateTime start = DateTime.Now;

            var QueueAndRun = this.AllJobs.Select(kv => kv.Value).Where(delegate (Job j) {
                var s = j.Status;
                if (s == JobStatus.FailedOrCanceled)
                    return false;
                if (s == JobStatus.FinishedSuccessful)
                    return false;
                if (s == JobStatus.PreActivation)
                    return false;
                return true;
            }).ToArray();
            if(QueueAndRun.Length == 0) {
                JustFinished = null;
                return 0;
            }

            while(true) {
                
                if(TimeOutSeconds > 0) {
                    double RuntimeSoFar = (DateTime.Now - start).TotalSeconds;
                    if(RuntimeSoFar > TimeOutSeconds) {
                        Console.WriteLine("Timeout.");
                        JustFinished = null;
                        return QueueAndRun.Length;
                    }
                }

                foreach(var J in QueueAndRun) {
                    var s = J.Status;
                    if(s == JobStatus.FailedOrCanceled || s == JobStatus.FinishedSuccessful) {
                        JustFinished = J;
                        return QueueAndRun.Length;
                    }
                }

                Thread.Sleep((int)(1000*PollingIntervallSeconds));
            }
            
        }


        /// <summary>
        /// Records the control object <paramref name="C"/> in an internal list, for its entire lifetime,
        /// and provides an index for it. 
        /// </summary>
        public int RegisterControl(AppControl C) {
            int max = 0;
            foreach (var t in RegisteredControls) {
                if (object.ReferenceEquals(t.Item1, C))
                    return t.Item2;
                max = Math.Max(t.Item2, max);
            }

            RegisteredControls.Add(Tuple.Create(C, max + 1));
            return max + 1;
        }


    }


    /*
    public static class MetaJobManager {

        static Dictionary<string, BatchProcessorClient> m_Computers;

        /// <summary>
        /// 
        /// </summary>
        static public Dictionary<string, BatchProcessorClient> Computers {
            get {
                if (m_Computers == null) {
                    m_Computers = new Dictionary<string, BatchProcessorClient>();
                }

                return m_Computers;
            }
        }

    }
    */
}




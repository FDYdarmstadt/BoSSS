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
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// A compute-job for the meta-job-manager.
    /// </summary>
    public class Job {

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="name">
        /// See <see cref="Name"/>.
        /// </param>
        /// <param name="solver"></param>
        public Job(string name, Type solver) {
            if(name.IsEmptyOrWhite()) {
                int i = 1;
                string newName;
                do {
                    newName = "EmptyJobName_" + i;
                    i++;
                } while(InteractiveShell.WorkflowMgm.AllJobs.ContainsKey(newName));
                
                Console.WriteLine($"Empty job name - picking new name '{newName}'");
                name = newName;
            }
            this.Solver = solver;
            this.Name = name;
            this.SessionReqForSuccess = true;
            
            if (InteractiveShell.WorkflowMgm.AllJobs.ContainsKey(name)) {
                throw new ArgumentException("Job with name '" + name + "' is already defined in the workflow management.");
            }
            InteractiveShell.WorkflowMgm.AllJobs.Add(name, this);
            if (string.IsNullOrWhiteSpace(InteractiveShell.WorkflowMgm.CurrentProject)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }
        }

        /// <summary>
        /// The name for the job; note that it serves as a unique, persistent identifier within the work flow management.
        /// Usually, set equal to the session name <see cref="BoSSS.Solution.Control.AppControl.SessionName"/>.
        /// </summary>
        public string Name {
            private set;
            get;
        }

        
        /// <summary>
        /// The memory (in MB) that is reserved for every core
        /// </summary>
        public string MemPerCPU {
            set;
            get;
        }

        /*
        string m_BatchProcessorIdentifierToken;

        /// <summary>
        /// (Optional) object used by batch processor (after calling <see cref="BatchProcessorClient.Submit(Job)"/>)
        /// in order to identify the job.
        /// </summary>
        public string BatchProcessorIdentifierToken {
            private set {
                m_BatchProcessorIdentifierToken = value;
            }
            get {
                var bpcToken = m_BatchProcessorIdentifierToken;

                if(m_BatchProcessorIdentifierToken.IsNullOrEmpty()) {
                    var directories = GetAllDeploymantDirectories();

                    if(directories == null || directories.Length <= 0) {
                        return null;
                    }
                    Array.Sort(directories, FuncComparerExtensions.ToComparer((DirectoryInfo a, DirectoryInfo b) => DateTime.Compare(a.CreationTime, b.CreationTime)));
                    DirectoryInfo _DD = directories.Last();
                    var DD = _DD.FullName;

                    try {
                        var l = File.ReadAllText(Path.Combine(DD, "IdentifierToken.txt"));
                        m_BatchProcessorIdentifierToken = l.Trim();

                    } catch(Exception) {
                        // job was probably deployed, but never submitted
                        // ignore this.
                    }
                }

                return m_BatchProcessorIdentifierToken;
            }
        }

        /// <summary>
        /// Some internal object that the job keeps for the batch processor
        /// </summary>
        object BatchProcessorObject;
        */

        /// <summary>
        /// Class which contains the main-method of the solver (or general application to launch).
        /// </summary>
        public Type Solver {
            get;
            private set;
        }

        /// <summary>
        /// Assembly which defines <see cref="Solver"/>.
        /// </summary>
        public Assembly EntryAssembly {
            get {
                return Solver.Assembly;
            }
        }

        /// <summary>
        /// Add dependent assemblies, may also contain stuff like mscorlib.dll.
        /// </summary>
        public IEnumerable<Assembly> AllDependentAssemblies {
            get {
                HashSet<Assembly> assiList = new HashSet<Assembly>();
                GetAllAssemblies(this.EntryAssembly, assiList, Path.GetDirectoryName(EntryAssembly.Location));
                return assiList.ToArray();
            }
        }

        /// <summary>
        /// Recursive collection of all dependencies of some assembly.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="assiList">
        /// Output, list where all dependent assemblies are collected.
        /// </param>
        /// <param name="SearchPath">
        /// Path to search for assemblies
        /// </param>
        private static void GetAllAssemblies(Assembly a, HashSet<Assembly> assiList, string SearchPath) {
            if (assiList.Contains(a))
                return;
            assiList.Add(a);

            string fileName = Path.GetFileName(a.Location);
            var allMatch = assiList.Where(_a => Path.GetFileName(_a.Location).Equals(fileName)).ToArray();
            if(allMatch.Length > 1) {
                throw new ApplicationException("internal error in assembly collection.");
            }


            foreach (AssemblyName b in a.GetReferencedAssemblies()) {
                Assembly na;
                try {
                    na = Assembly.Load(b);
                } catch (FileNotFoundException) {
                    string[] AssiFiles = ArrayTools.Cat(Directory.GetFiles(SearchPath, b.Name + ".dll"), Directory.GetFiles(SearchPath, b.Name + ".exe"));
                    if(AssiFiles.Length != 1) {
                        //throw new FileNotFoundException("Unable to locate assembly '" + b.Name + "'.");
                        Console.WriteLine("Skipping: " + b.Name);
                        continue;
                    }
                    na = Assembly.LoadFile(AssiFiles[0]);

                }

                GetAllAssemblies(na, assiList, SearchPath);
            }
        }

        List<Tuple<byte[], string>> m_AdditionalDeploymentFiles = new List<Tuple<byte[], string>>();

        /// <summary>
        /// Additional data files which will be deployed in the <see cref="DeploymentDirectory"/> together with the
        /// assemblies.
        ///  - 1st item: file content
        ///  - 2nd item: file name
        /// </summary>
        public IList<Tuple<byte[], string>> AdditionalDeploymentFiles {
            get {
                return m_AdditionalDeploymentFiles;
            }
        }

        /// <summary>
        /// Adds a text file to <see cref="AdditionalDeploymentFiles"/>.
        /// </summary>
        /// <param name="Content">Content of the text file.</param>
        /// <param name="FileName">Designated name of the file in the deployment directory.</param>
        public void AddTextFile(string Content, string FileName) {
            TestActivation();

            byte[] Bytes = Encoding.UTF8.GetBytes(Content);
            AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(Bytes, FileName));
        }

        /// <summary>
        /// Specifies the '--control' statement for application startup; overrides any startup arguments (<see cref="CommandLineArguments"/>)
        /// set so far. Session/Job and project name are set automatically.
        /// </summary>
        /// <param name="Args">Path to control file, resp. C#-statement (starting with 'cs:')</param>
        /// <remarks>
        /// Note: we pass the startup-arguments through environment variables, which is 
        /// (a bit) more robust (with respect to escape-characters, etc.) than 
        /// command line arguments.
        /// </remarks>
        public void SetControlStatement(string Args) {
            TestActivation();

            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(PrjName)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }

            // we pass the startup-arguments through environment variables, which is 
            // (a bit) more robust (with respect to escape-characters, etc.) than 
            // command line arguments.
            string[] args = new string[] {
                "--control", Args,
                "--prjnmn" , PrjName,
                "--sesnmn", this.Name
            };

            m_EnvironmentVars.Clear();
            for (int i = 0; i < args.Length; i++) {
                m_EnvironmentVars.Add("BOSSS_ARG_" + i, args[i]);
            }
        }

        /// <summary>
        /// Specifies the command line arguments for application startup; overrides any startup arguments (<see cref="CommandLineArguments"/>)
        /// set so far. 
        /// </summary>
        /// <param name="Args">startup string/command line arguments</param>
        /// <remarks>
        /// Note: we pass the startup-arguments through environment variables, which is 
        /// (a bit) more robust (with respect to escape-characters, etc.) than 
        /// command line arguments.
        /// </remarks>
        public void MySetCommandLineArguments(params string[] Args) {
            TestActivation();

            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(PrjName)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }


            m_EnvironmentVars.Clear();
            for (int i = 0; i < Args.Length; i++) {
                m_EnvironmentVars.Add("BOSSS_ARG_" + i, Args[i]);
            }
        }

   
        BoSSS.Solution.Control.AppControl m_ctrl;
        int m_ctrl_index;
        string ControlName = null;

        private void FiddleControlFile(BatchProcessorClient bpc) {
            if(m_ctrl == null)
                return;

            // check the database 
            // ==================
            IDatabaseInfo ctrl_db = m_ctrl.GetDatabase();
            if(bpc.AllowedDatabases != null && bpc.AllowedDatabases.Count > 0) {
                
                IDatabaseInfo newDb = null;
                if(ctrl_db == null) {
                    newDb = bpc.AllowedDatabases[0];
                } else {
                    bool ok = false;
                    foreach(var allow_dba in bpc.AllowedDatabases) {
                        if(allow_dba.Equals(ctrl_db)) {
                            ok = true;
                            break;
                        }
                    }

                    if(!ok)
                        newDb = bpc.AllowedDatabases[0];
                }

                if(newDb != null) {
                    Console.WriteLine("Resetting database for control object to " + newDb.ToString());

                    //newDb.AlternateDbPaths



                    m_ctrl.SetDatabase(newDb);
                    ctrl_db = newDb;
                }

                Console.WriteLine("Submitting job with the following database info: ");
                Console.WriteLine("Primary: " + m_ctrl.DbPath);
                if(ctrl_db.AlternateDbPaths != null && ctrl_db.AlternateDbPaths.Length > 0) {
                    int cnt = 0;
                    foreach (var t in ctrl_db.AlternateDbPaths) {
                        Console.WriteLine($" Alternative[{cnt}]: {t.DbPath}, MachineFilter: '{t.MachineFilter}'");
                        cnt++;
                    }
                } else {
                    Console.WriteLine("No alternative paths specified.");
                }
            } else {
                Console.WriteLine("");
            } 

            // check grid & restart info
            // =========================

            if(ctrl_db != null) {
                if(!m_ctrl.GridGuid.Equals(Guid.Empty)) {

                    var GridIn_ctrl_db = ctrl_db.Grids.FirstOrDefault(GrdInf => GrdInf.ID.Equals(m_ctrl.GridGuid));

                    if(GridIn_ctrl_db == null) {
                        Console.WriteLine($"Grid {m_ctrl.GridGuid} is not present in database - copy to target system...");

                        var grid2copy = InteractiveShell.AllGrids.FirstOrDefault(dbGrid => dbGrid.ID.Equals(m_ctrl.GridGuid));
                        if(grid2copy == null) {
                            // maybe replace exception with a warning, if job should be tried anyway
                            throw new IOException($"Unable to find grid '{m_ctrl.GridGuid}' in any database - job will most likely crash.");
                        } else {
                            grid2copy.Copy(ctrl_db);
                        }
                        
                        Console.WriteLine("done.");
                    }
                } else {
                    Console.Error.WriteLine($"Warning: no grid seems to be specified for the job to submit.");
                }

                if(m_ctrl.RestartInfo != null) {
                    Guid Rstsess_guid = m_ctrl.RestartInfo.Item1;


                    var Rstsess_ctrl_db = ctrl_db.Sessions.FirstOrDefault(sinf => sinf.ID.Equals(Rstsess_guid));

                    if(Rstsess_ctrl_db == null) {
                        Console.WriteLine($"Session {m_ctrl.GridGuid} to restart from is not present in database - copy to target system...");

                        var sess_2copy = InteractiveShell.AllSessions.FirstOrDefault(sinf => sinf.ID.Equals(Rstsess_guid));
                        if(sess_2copy == null) {
                            // maybe replace exception with a warning, if job should be tried anyway
                            throw new IOException($"Unable to find session '{sess_2copy}' in any database - job will most likely crash.");
                        } else {
                            sess_2copy.Copy(ctrl_db);
                        }
                        
                        Console.WriteLine("done.");
                    }
                }

            } else {
                Console.Error.WriteLine($"Warning: no database is set for the job to submit; nothing ma be saved.");
            }

            // finally, serialize the object
            // =============================
            {
                string text;
                m_ctrl_index = -1;
                if(m_ctrl.GeneratedFromCode) {
                    text = m_ctrl.ControlFileText;
                    ControlName = "control.cs";
                    m_ctrl_index = m_ctrl.ControlFileText_Index;
                } else {
                    text = m_ctrl.Serialize();
                    ControlName = "control.obj";
                }
                byte[] buffer = Encoding.UTF8.GetBytes(text);

                int remIdx = AdditionalDeploymentFiles.IndexWhere(tt => tt.Item2 == ControlName);
                if(remIdx >= 0)
                    AdditionalDeploymentFiles.RemoveAt(remIdx);

                AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(buffer, ControlName));
            }
        }



        /// <summary>
        /// Specifies the control object for application startup; overrides any startup arguments (<see cref="CommandLineArguments"/>)
        /// set so far.
        /// </summary>
        public void SetControlObject(BoSSS.Solution.Control.AppControl ctrl) {
            TestActivation();

            // serialize control object
            // ========================

            // grid function hack:
            if (ctrl.GridFunc != null) {
                Console.WriteLine("Control object contains grid function. Trying to Serialize the grid...");
                var dbi = ctrl.GetDatabase();
                if (dbi == null) {
                    throw new NotSupportedException("If a gird function is specified (instead of a grid id), a database must be specified to save the gird (when using the job manager).");
                }

                Foundation.Grid.IGrid g = ctrl.GridFunc();
                Guid id = dbi.SaveGrid(ref g);

                ctrl.GridFunc = null;
                ctrl.GridGuid = id;
                Console.WriteLine("Control object modified.");

            }


            ctrl.VerifyEx();
            m_ctrl = ctrl;
            m_ctrl.ProjectName = InteractiveShell.WorkflowMgm.CurrentProject;

            // note: serialization is done later, immediately before deployment,
            // since we may need to fix database issues (path on batch system, evtl. transfer of grid)

            m_ctrl_index = -1;
            if(m_ctrl.GeneratedFromCode) {
                ControlName = "control.cs";
                m_ctrl_index = m_ctrl.ControlFileText_Index;
            } else {
                ControlName = "control.obj";
            }


            // Project & Session Name
            // ======================
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(PrjName)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }

            string[] args = new string[] {
                "--control", "control.obj",
                "--prjnmn", PrjName,
                "--sesnmn", this.Name
            };
            if (m_ctrl_index >= 0) {
                ArrayTools.Cat(args, "--pstudy_case", m_ctrl_index.ToString());
            }

            m_EnvironmentVars.Clear();
            for (int i = 0; i < args.Length; i++) {
                m_EnvironmentVars.Add("BOSSS_ARG_" + i, args[i]);
            }

            // 
        }

        /// <summary>
        /// Returns the control object linked to this job
        /// </summary>
        public BoSSS.Solution.Control.AppControl GetControl() {
           
            return m_ctrl;
        }

        
        string[] m_CommandLineArguments = new string[0];

        /// <summary>
        /// Startup arguments for the process. 
        /// </summary>
        public string[] CommandLineArguments {
            get {
                return m_CommandLineArguments;
            }
        }

        Dictionary<string, string> m_EnvironmentVars = new Dictionary<string, string>();

        /// <summary>
        /// Additional environment variables for the process. 
        /// </summary>
        public IDictionary<string, string> EnvironmentVars {
            get {
                return m_EnvironmentVars;
            }
        }

        /// <summary>
        /// Data record for each deployment, i.e. each attempt to run a job.
        /// </summary>
        public class Deployment {


            Job m_owner;

            internal Deployment(DirectoryInfo _DeploymentDirectory, Job _owner, object __optInfo = null) {
                m_owner = _owner;
                if(_DeploymentDirectory != null) {
                    if(!_DeploymentDirectory.Exists)
                        throw new ArgumentException($"Directory {_DeploymentDirectory} does not exist -- must be either an existing directory or null.");
                }
                m_DeploymentDirectory = _DeploymentDirectory;
                optInfo = __optInfo;

            }

            internal Deployment(ISessionInfo _Session, Job _owner) {
                m_owner = _owner;
                if(_Session == null) {
                    throw new ArgumentNullException();
                }
                this.Session = _Session;
            }


            DirectoryInfo m_DeploymentDirectory;

            /// <summary>
            /// The deployment directory, accessible from the local machine (e.g. a mounted path if the job is deployed to another computer);
            /// This should be an absolute path.
            /// </summary>
            public DirectoryInfo DeploymentDirectory {
                get {
                    return m_DeploymentDirectory;
                }
            }

            /// <summary>
            /// Creation time of <see cref="DeploymentDirectory"/>
            /// </summary>
            public DateTime CreationDate {
                get {
                    return DeploymentDirectory.CreationTime;
                }
            }


            string m_BatchProcessorIdentifierToken = null;

            /// <summary>
            /// Job ID uses by the respective batch processor
            /// </summary>
            public string BatchProcessorIdentifierToken {
                get {
                    if(m_BatchProcessorIdentifierToken.IsNullOrEmpty()) {
                        string DD = this.DeploymentDirectory.FullName;
                        if(DD != null && Directory.Exists(DD))
                        try {
                            var l = File.ReadAllText(Path.Combine(DD, "IdentifierToken.txt"));
                            m_BatchProcessorIdentifierToken = l.Trim();

                        } catch(Exception) {
                            // job was probably deployed, but never submitted
                            // ignore this.
                        }
                    }

                    return m_BatchProcessorIdentifierToken;
                }
            }

            /// <summary>
            /// The deepest-level directory of <see cref="DeploymentDirectory"/>.
            /// </summary>
            public string RelativeDeploymentDirectory {
                get {
                    if(DeploymentDirectory != null) {
                        var parts = DeploymentDirectory.FullName.Split(new[] { Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar }, StringSplitOptions.RemoveEmptyEntries);
                        return parts.Last();
                    }

                    if(Session != null) {
                        var parts = Session.DeployPath.Split(new[] { Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar }, StringSplitOptions.RemoveEmptyEntries);
                        return parts.Last();
                    }

                    return null; // unable to determine
                }
            }

            internal object optInfo;

            JobStatus? StatusCache = null;
            int? ExitCodeCache = null;

            /// <summary>
            /// Status of the deployment
            /// </summary>
            public JobStatus Status {
                get {
                    if(StatusCache != null)
                        return StatusCache.Value;

                    if(m_owner.AssignedBatchProc == null)
                        return JobStatus.PreActivation;

                    var s = m_owner.AssignedBatchProc.EvaluateStatus(this.BatchProcessorIdentifierToken, this.optInfo, this.DeploymentDirectory.FullName);
                    string ExitCodeStr = s.ExitCode.HasValue ? s.ExitCode.Value.ToString() : "null";

                    if(s.Item1 == JobStatus.FinishedSuccessful) {
                        StatusCache = s.Item1;
                        if(s.ExitCode == null || s.ExitCode != 0)
                            throw new ApplicationException($"Error in implementation of {m_owner.AssignedBatchProc}: job marked as {JobStatus.FinishedSuccessful}, but exit code is {ExitCodeStr} -- expecting 0.");
                    }
                        
                    if(s.Item1 == JobStatus.FailedOrCanceled) {
                        StatusCache = s.Item1;
                        if(s.ExitCode == null || s.ExitCode != 0)
                            throw new ApplicationException($"Error in implementation of {m_owner.AssignedBatchProc}: job marked as {JobStatus.FailedOrCanceled}, but exit code is {ExitCodeStr} -- expecting any number except 0.");

                    }

                    return s.Item1;
                }
            }

            /// <summary>
            /// If deployment terminated, the exit coede of the app; otherwise, null
            /// </summary>
            public int? ExitCode {
                get {
                    _ = Status; // happy side-effect programming
                    return this.ExitCodeCache;
                }
            }

            /// <summary>
            /// If available the database session
            /// </summary>
            public ISessionInfo Session {
                get;
                internal set;
            }

            /// <summary>
            /// 
            /// </summary>
            public override string ToString() {
                return "Job "
                    + (this.BatchProcessorIdentifierToken != null ? BatchProcessorIdentifierToken : "null") + ", "
                    + this.Status
                    + " '" + this.RelativeDeploymentDirectory + "'"
                    + (m_owner.AssignedBatchProc != null ? (" @ " + m_owner.AssignedBatchProc) : "");
            }
        }

        List<Deployment> m_Deployments = new List<Deployment>();


        string JobDirectoryBaseName() {
            string Exe = Path.GetFileNameWithoutExtension(EntryAssembly.Location);
            string Proj = InteractiveShell.WorkflowMgm.CurrentProject;
            string Sess = Name;

            return Proj
                //+ "-" + Sess 
                + "-" + Exe;
        }


        /// <summary>
        /// Returns the directory where the assemblies for <paramref name="myJob"/> 
        /// are deployed if <paramref name="myJob"/> is assigned to this batch processor.
        /// </summary>
        string GetNewDeploymentDir() {
            if(AssignedBatchProc == null)
                throw new NotSupportedException("Job is not activated yet.");

            if (!Path.IsPathRooted(AssignedBatchProc.DeploymentBaseDirectory))
                throw new IOException($"Deployment base directory for {AssignedBatchProc.ToString()} must be rooted/absolute, but '{AssignedBatchProc.DeploymentBaseDirectory}' is not.");

            string ShortName = JobDirectoryBaseName();
            string DeployDir;
            int Counter = 0;
            do {
                string Suffix = Counter > 0 ? "-" + Counter : "";
                string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss");
                DeployDir = Path.Combine(AssignedBatchProc.DeploymentBaseDirectory, ShortName + DateNtime + Suffix);
                Counter++;
            } while (Directory.Exists(DeployDir) == true);

            return DeployDir;
        }

        /// <summary>
        /// All deployment directories which potentially could match the job on the current batch processor.
        /// </summary>
        DirectoryInfo[] GetAllUnkonwnExistingDeployDirectories() {
            using (var tr = new FuncTrace()) {
                if(this.GetControl() == null && this.EntryAssembly == null)
                    throw new NotSupportedException("Insufficient information to find job deployments.");

                if (!Path.IsPathRooted(AssignedBatchProc.DeploymentBaseDirectory))
                    throw new IOException($"Deployment base directory for {AssignedBatchProc.ToString()} must be rooted/absolute, but '{AssignedBatchProc.DeploymentBaseDirectory}' is not.");

                var jobControl = this.GetControl();
                if (jobControl == null)
                    return new DirectoryInfo[0];

                bool DirMatch(DirectoryInfo dir1, DirectoryInfo dir2) {
                    string _dir1 = dir1.FullName.TrimEnd(Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar);
                    string _dir2 = dir2.FullName.TrimEnd(Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar);
                    return _dir1.Equals(_dir2);
                }

                bool DirIsKnown(DirectoryInfo dir) {
                    foreach(var dep in m_Deployments) {
                        if(DirMatch(dep.DeploymentDirectory, dir))
                            return true;
                    }
                    return false;
                }

                // find all deployment directories relevant for project & job
                // ==========================================================
                string ShortName = this.JobDirectoryBaseName();

                DirectoryInfo[] AllDirs;
                using (new BlockTrace("DIRECTORY_QUERY", tr)) {
                    AllDirs = Directory.GetDirectories(this.AssignedBatchProc.DeploymentBaseDirectory, ShortName + "*").Select(str => new DirectoryInfo(str)).ToArray();
                }

                // filter appropriate ones 
                // =======================
                var filtDirs = new List<DirectoryInfo>();
                using(new BlockTrace("DIRECTORY_FILTERING", tr)) {
                    foreach(DirectoryInfo dir in AllDirs) {
                        if(!DirIsKnown(dir)) {


                            string ControlObj = Path.Combine(dir.FullName, "control.obj");
                            if(File.Exists(ControlObj)) {
                                var ctrl = BoSSS.Solution.Control.AppControl.Deserialize(File.ReadAllText(ControlObj));
                                if(InteractiveShell.WorkflowMgm.JobAppControlCorrelation(this, ctrl)) {
                                    filtDirs.Add(new DirectoryInfo(dir.FullName));
                                    continue;
                                }
                            }

                            string ControlScript = Path.Combine(dir.FullName, "control.cs");
                            if(File.Exists(ControlScript)) {
                                int control_index = 0;
                                int i = CommandLineArguments.IndexWhere(arg => arg == "--pstudy_case");
                                if(i >= 0) {
                                    control_index = int.Parse(CommandLineArguments[i + 1]);
                                }

                                var ctrl = BoSSS.Solution.Control.AppControl.FromFile(ControlScript, jobControl.GetType(), control_index);
                                if(InteractiveShell.WorkflowMgm.JobAppControlCorrelation(this, ctrl)) {
                                    filtDirs.Add(dir);
                                    continue;
                                }
                            }
                        }
                    }
                }

                // return
                // ======
                return filtDirs.ToArray();
            }
        }

        void UpdateDeployments() {
            
            // determine all new sessions
            // ==========================
            HashSet<Guid> KnownSessionGuids = new HashSet<Guid>();
            foreach(var dep in m_Deployments) {
                if(dep.Session != null) {
                    KnownSessionGuids.Add(dep.Session.ID);
                }
            }

            ISessionInfo[] AllNewSessions = InteractiveShell.WorkflowMgm.Sessions
                .Where(sinf => !KnownSessionGuids.Contains(sinf.ID)) // for performance reasons, filter sessions that we already know
                .Where(sinf => InteractiveShell.WorkflowMgm.SessionInfoJobCorrelation(sinf, this)).ToArray();

            // add all new deployment directories
            // ==================================
            foreach(DirectoryInfo dirInfo in this.GetAllUnkonwnExistingDeployDirectories()) {
                m_Deployments.Add(new Deployment(dirInfo, this));
            }

            // finally, correlate the new sessions to the deployments
            // ======================================================
            for(int i = 0; i < AllNewSessions.Length; i++) {
                var sinf = AllNewSessions[i];

                string RelDeployPath = sinf.DeployPath.Split(new[] { Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar }, StringSplitOptions.RemoveEmptyEntries).Last();
                

                foreach(var dep in m_Deployments) {
                    if(dep.RelativeDeploymentDirectory == RelDeployPath) {
                        // bingo
                        dep.Session = sinf;
                        AllNewSessions[i] = null;
                    }
                }
            }

            for(int i = 0; i < AllNewSessions.Length; i++) {
                var sinf = AllNewSessions[i];
                if(sinf != null) {
                    // delpoyment was deleted, but some session was found in the database.
                    m_Deployments.Add(new Deployment(sinf, this));
                }
            }
        }


        /// <summary>
        /// All deployments/attempts for this job which can be found in the data sources
        /// (sessions in BoSSS databases, deployment directories and batch processor queues)
        /// </summary>
        public IReadOnlyList<Deployment> AllDeployments {
            get {
                UpdateDeployments();
                return m_Deployments.AsReadOnly();
            }
        }


        /// <summary>
        /// Deletes any trace of previous calculations for this job.
        /// </summary>
        void DeleteOldDeploymentsAndSessions() {
            foreach(var dep in AllDeployments) {
                if(dep.Session != null) {
                    dep.Session.Delete(true);
                }

                if(dep.DeploymentDirectory != null && dep.DeploymentDirectory.Exists) {
                    Exception op(int iTry) {
                        dep.DeploymentDirectory.Delete(true);
                        return null;
                    }


                    MetaJobMgrIO.RetryIOop(op, "deletion of directory '" + dep.DeploymentDirectory.FullName + "'", false);
                }
            }

            m_Deployments.Clear();
        }



        

        int m_NumberOfMPIProcs = 1;

        /// <summary>
        /// %
        /// </summary>
        public int NumberOfMPIProcs {
            get {
                return m_NumberOfMPIProcs;
            }
            set {
                TestActivation();
                m_NumberOfMPIProcs = value;
            }
        }

        bool m_UseComputeNodesExclusive = false;

        /// <summary>
        /// If true, the batch system should try not to run any other jobs in parallel on the assigned compute nodes.
        /// </summary>
        public bool UseComputeNodesExclusive {
            get {
                return m_UseComputeNodesExclusive;
            }
            set {
                TestActivation();
                m_UseComputeNodesExclusive = value;
            }
        }

        
        /// <summary>
        /// Returns all session which can be correlated to this job,
        /// <see cref="WorkflowMgm.SessionInfoJobCorrelation"/>.
        /// </summary>
        public ISessionInfo[] AllSessions {
            get {
                List<ISessionInfo> ret = new List<ISessionInfo>();
                foreach(var dep in AllDeployments) {
                    if(dep.Session != null)
                        ret.Add(dep.Session);
                }
                return ret.ToArray();
            }
        }

      


        /// <summary>
        /// Returns the session which is the result of this job.
        /// </summary>
        public ISessionInfo LatestSession {
            get {
                ISessionInfo[] RR = this.AllSessions;
                return RR.Length > 0 ? RR.OrderBy(si => si.CreationTime).Last() : null;
            }
        }

        /// <summary>
        /// most recent (according to date) entry in <see cref="AllDeployments"/>
        /// </summary>
        public Deployment LatestDeployment {
            get {
                var deps = this.AllDeployments.ToArray();
                if(deps.Length > 0)
                    return deps.OrderBy(d => d.CreationDate).Last();
                else
                    return null;
            }
        }



        /// <summary>
        /// The content of the standard output stream as a string.
        /// </summary>
        public string Stdout {
            get {
                using(new FuncTrace()) {
                    if(AssignedBatchProc == null)
                        throw new NotSupportedException("Job is not activated.");

                    

                    string stdout = "";

                    Exception op(int itry) {

                        string StdoutFile = null;
                        if(this.Status == JobStatus.FinishedSuccessful) {
                            if(this.LatestSession != null) {
                                StdoutFile = this.LatestSession.FilesInSessionDir("stdout.0.txt").FirstOrDefault();
                            }
                        }
                        if(StdoutFile == null) {
                            var ld = LatestDeployment;
                            if(ld != null) {
                                StdoutFile = AssignedBatchProc.GetStdoutFile(ld.BatchProcessorIdentifierToken, ld.DeploymentDirectory.FullName);
                            }
                        }
    
                        if(StdoutFile != null && File.Exists(StdoutFile)) {
                            using(FileStream stream = File.Open(StdoutFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)) {
                                using(StreamReader reader = new StreamReader(stream)) {
                                    stdout = reader.ReadToEnd();
                                }
                            }
                        } else {
                            stdout = "";
                        }
                        return null;
                    }

                    MetaJobMgrIO.RetryIOop(op, "reading of stdout file", true);

                    return stdout;
                }
            }
        }

        /// <summary>
        /// The content of the standard error stream as a string.
        /// </summary>
        public string Stderr {
            get {
                using(new FuncTrace()) {
                    if(AssignedBatchProc == null)
                        throw new NotSupportedException("Job is not activated.");

                    string stderr = "";

                    Exception op(int itry) {
                        string StderrFile = null;
                        if(this.Status == JobStatus.FinishedSuccessful) {
                            if(this.LatestSession != null) {
                                StderrFile = this.LatestSession.FilesInSessionDir("stderr.0.txt").FirstOrDefault();
                            }
                        }
                        if(StderrFile == null) {
                            var ld = LatestDeployment;
                            if(ld != null) {
                                StderrFile = AssignedBatchProc.GetStderrFile(ld.BatchProcessorIdentifierToken, ld.DeploymentDirectory.FullName);
                            }

                        }

                        if(StderrFile != null && File.Exists(StderrFile)) {
                            using(FileStream stream = File.Open(StderrFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)) {
                                using(StreamReader reader = new StreamReader(stream)) {
                                    stderr = reader.ReadToEnd();
                                }
                            }
                        } else {
                            stderr = "";
                        }
                        return null;
                    }

                    MetaJobMgrIO.RetryIOop(op, "reading of stderr file", true);

                    return stderr;
                }
            }
        }

        /// <summary>
        /// Should open a separate console window showing the rolling output text of the job.
        /// </summary>
        public void ShowOutput() {
            if (AssignedBatchProc == null)
                throw new NotSupportedException("Job is not activated.");

            var ld = LatestDeployment;
            if(ld == null) {
                Console.WriteLine("No deployment available so far.");
                return;
            }

            string StderrFile = AssignedBatchProc.GetStderrFile(ld.BatchProcessorIdentifierToken, ld.DeploymentDirectory.FullName);
            string StdoutFile = AssignedBatchProc.GetStdoutFile(ld.BatchProcessorIdentifierToken, ld.DeploymentDirectory.FullName);

            if (StdoutFile != null && StderrFile != null) {
                ProcessStartInfo psi = new ProcessStartInfo();
                psi.FileName = typeof(btail.TailMain).Assembly.Location;
                btail.TailMain.SetArgs(psi, StdoutFile, StderrFile);

                Console.WriteLine("Starting console...");
                Console.WriteLine("(You may close the new window at any time, the job will continue.)");

                Process p = Process.Start(psi);
            } else {
                Console.WriteLine("No known standard output/error path yet - try again later.");
            }
        }

        /// <summary>
        /// After calling <see cref="BatchProcessorClient.Submit(Job)"/>, this job
        /// is assigned to the respective batch processor, which is recorded in this member.
        /// </summary>
        public BatchProcessorClient AssignedBatchProc {
            get;
            private set;
        }

        JobStatus? statusCache;

        /// <summary>
        /// Whats up with this job?
        /// </summary>
        public JobStatus Status {
            get {
                return GetStatus(WriteHints: false);
            }
        }

        /// <summary>
        /// How often has this job been submitted to the batch system?
        /// </summary>
        public int SubmitCount {
            get {
                return m_Deployments.Count;
            }
        }

        /// <summary>
        /// - if true (default), a successful session in a database is required for this job to be a success
        /// - false must be set for jobs which do not write to a database.
        /// </summary>
        public bool SessionReqForSuccess {
            get;
            set;
        }



        /*
        /// <summary>
        /// All directories which can be found to which this job has ever been deployed
        /// </summary>
        public DirectoryInfo[] GetAllDeploymantDirectories() {
            DirectoryInfo[] directories = AssignedBatchProc.GetAllExistingDeployDirectories(this);
            if(this.DeploymentDirectory != null) {
                DirectoryInfo diAdd = new DirectoryInfo(this.DeploymentDirectory);
                if(diAdd.Exists) {
                    if(directories == null)
                        directories = new DirectoryInfo[0];

                    diAdd.AddToArray(ref directories);
                }
            }

            if(directories == null)
                directories = new DirectoryInfo[0];

            return directories;
        }
        */


        /// <summary>
        /// override of parameter in <see cref="Activate(BatchProcessorClient, bool)"/> for testing purpose.
        /// </summary>
        public static bool UndocumentedSuperHack = false;

        /// <summary>
        /// Activates the job; can only be called once in this objects lifetime.
        /// </summary>
        /// <remarks>
        /// The job is submitted to the batch processor if no 
        /// By using a unique name, see <see cref="Name"/>, the job becomes persistent: e.g., assume 
        /// that this method is called within a worksheet. Later, the worksheet is restarted, or closed and reopened,
        /// thus this method is called a second time. By identifying the job in the batch processor with a unique name and
        /// by tagging sessions with the job name, it is ensured that the job is _not_ submitted a second time,
        /// resp. every time the worksheet is restarted.
        /// </remarks>
        /// <param name="bpc"></param>
        /// <param name="DeleteOldDeploymentsAndSessions">
        /// Override job persistence: traces of old jobs will be deleted.
        /// </param>
        public void Activate(BatchProcessorClient bpc, bool DeleteOldDeploymentsAndSessions = false) {
            using (var tr = new FuncTrace()) {
                // ============================================
                // ensure that this method is only called once.
                // ============================================
                if (this.AssignedBatchProc != null)
                    throw new NotSupportedException("Job can only be activated once.");
                AssignedBatchProc = bpc;

                if(DeleteOldDeploymentsAndSessions || UndocumentedSuperHack)
                    this.DeleteOldDeploymentsAndSessions();


                // ================
                // status
                // ================

                var stat = GetStatus(true);
                if(stat != JobStatus.Unknown)
                    return;

                /*
                // ================
                // check job status
                // ================
                string DeployDir;
                int SubmitCount;
                var status = GetStatus(out SubmitCount, out DeployDir);
                this.DeploymentDirectory = DeployDir;

                // ==============
                // check Database
                // ==============
                var RR = this.AllSessions;

                // =================
                // decide what to do
                // =================
                switch (status) {
                    case JobStatus.PreActivation:
                        DeployDir = null;
                        this.DeploymentDirectory = null;
                        Console.WriteLine("Job not submitted yet, or no result session is known - starting submission.");
                        break;

                    case JobStatus.FailedOrCanceled:
                        if (RR.Length <= 0) {
                            DeployDir = null;
                            this.DeploymentDirectory = null;
                            Console.WriteLine("Job is marked as failed by job manager, no database entry is found; performing new deployment and submission.");
                            break;
                        }
                        if (RetryCount <= SubmitCount) {
                            Console.WriteLine("Job is failed {0} times, maximum number of tries reached, no further action.", SubmitCount);
                            return;
                        } else {
                            DeployDir = null;
                            this.DeploymentDirectory = null;
                            Console.WriteLine("Job is failed, retrying (Submitted {0} times so far); performing new deployment and submission.", SubmitCount);
                            break;
                        }
                    case JobStatus.PendingInExecutionQueue:
                        Console.WriteLine("Job submitted, waiting for launch - no further action.");
                        return;

                    case JobStatus.FinishedSuccessful:
                        if (RR.Length <= 0) {
                            DeployDir = null;
                            this.DeploymentDirectory = null;
                            Console.WriteLine("Job is marked as success by job manager, but no session info in database is found; performing new deployment and submission.");
                            break;
                        }
                        ISessionInfo LatestSession = RR.OrderBy(sinf => sinf.CreationTime).Last();
                        Console.WriteLine("Job was successful (according to job manager), latest session related to job is:");
                        Console.WriteLine(LatestSession.ToString());
                        Console.WriteLine("No further action.");
                        return;

                    case JobStatus.InProgress:
                        Console.Write("Job has been started (according to job manager), ");
                        if (RR.Length > 0) {
                            Console.WriteLine("latest known session is:");
                            ISessionInfo LatestSession2 = RR.OrderBy(sinf => sinf.CreationTime).Last();
                            Console.WriteLine(LatestSession2.ToString());
                        } else {
                            Console.WriteLine("no session information available at this point.");
                        }
                        Console.WriteLine("No further action.");
                        return;

                    default:
                        throw new NotImplementedException();
                }
                */

                // ========================================================================
                // finally, it might be necessary to submit the job to the batch processor. 
                // ========================================================================

                // some database syncing might be necessary 
                FiddleControlFile(bpc);

                // deploy additional files
                string DeploymentDirectory = this.DeployExecuteables();
                
                // submit job
                using (new BlockTrace("JOB_SUBMISSION", tr)) {
                    var rr = bpc.Submit(this, DeploymentDirectory);
                    File.WriteAllText(Path.Combine(DeploymentDirectory, "IdentifierToken.txt"), rr.id);

                    Deployment dep = AllDeployments.SingleOrDefault(d => d.BatchProcessorIdentifierToken.Equals(rr.id));
                    if(dep == null)
                        m_Deployments.Add(new Deployment(new DirectoryInfo(DeploymentDirectory), this, rr.optJobObj));
                    else
                        dep.optInfo = rr.optJobObj;
                }
            }
        }


        
        /// <summary>
        /// Status evaluation, with optional additional information.
        /// </summary>
        /// <param name="WriteHints">
        /// If true, prints a reason for the returned status.
        /// </param>
        /// <returns></returns>
        public JobStatus GetStatus(bool WriteHints = true) {
            using (new FuncTrace()) {
                //SubmitCount = 0;
                //DD = null;
                if(AssignedBatchProc == null)
                    return JobStatus.PreActivation;

                if(!WriteHints) {
                    if(statusCache.HasValue)
                        return statusCache.Value;
                }

                // ================
                // identify success
                // ================

                Deployment[] DeploymentsSoFar = this.AllDeployments.ToArray();
                Deployment[] Success = DeploymentsSoFar.Where(dep => dep.Status == JobStatus.FinishedSuccessful).ToArray();
                if(SessionReqForSuccess) {
                    ISessionInfo[] SuccessSessions = this.AllSessions.Where(si => si.SuccessfulTermination()).OrderByDescending(sess => sess.CreationTime).ToArray();
                    if(SuccessSessions.Length <= 0) {
                        // look twice
                        InteractiveShell.WorkflowMgm.ResetSessionsCache();
                        DeploymentsSoFar = this.AllDeployments.ToArray();
                        Success = DeploymentsSoFar.Where(dep => dep.Status == JobStatus.FinishedSuccessful).ToArray();
                        SuccessSessions = this.AllSessions.Where(si => si.SuccessfulTermination()).OrderByDescending(sess => sess.CreationTime).ToArray();
                    }
                    if(SuccessSessions.Length > 0) {
                        if(WriteHints)
                            Console.WriteLine($"Info: Found successful session \"{SuccessSessions.First()}\" -- job is marked as successful, no further action.");
                        this.statusCache = JobStatus.FinishedSuccessful;
                        return JobStatus.FinishedSuccessful;
                    }

                    if(Success.Length > 0) {
                        if(WriteHints)
                            Console.WriteLine($"Note: found {Success.Length} successful deployment(s), but job is configured to require a successful result session ('this.SessionReqForSuccess' is true), and none is found. {this.AllSessions.Length} sessions correlated to this job fount in total.");
                    }

                    if(WriteHints) {
                        var c = this.GetControl();
                        if(c != null && c.savetodb == false) {
                            throw new NotSupportedException("Configuration error: this job can never be successful: a successful result session ('this.SessionReqForSuccess' is true) is required, but control object is configured not to save in the database.");
                        }
                    }

                } else {
                    if(Success.Length > 0) {
                        if(WriteHints)
                            Console.WriteLine($"Info: Found successful deployment and no result session is expected ('this.SessionReqForSuccess' is false):  job is marked as successful, no further action.");
                        this.statusCache = JobStatus.FinishedSuccessful;
                        return JobStatus.FinishedSuccessful;
                    }
                }

                // ========================
                // identify running/waiting
                // ========================

                var inprog = DeploymentsSoFar.Where(dep => (dep.Status == JobStatus.InProgress)).ToArray();
                if(inprog.Length > 0) {
                    if(WriteHints)
                        Console.WriteLine($"Info: Job {inprog[0].BatchProcessorIdentifierToken} is currently executed on {this.AssignedBatchProc} -- no further action.");
                    return JobStatus.InProgress;
                }

                var inq = DeploymentsSoFar.Where(dep => (dep.Status == JobStatus.PendingInExecutionQueue)).ToArray();
                if(inq.Length > 0) {
                    if(WriteHints)
                        Console.WriteLine($"Info: Job {inq[0].BatchProcessorIdentifierToken} is currently waiting on {this.AssignedBatchProc} -- no further action.");
                    return JobStatus.PendingInExecutionQueue;
                }

                // =============
                // identify fail
                // =============

                // if we pass this point, we want to submit to a batch processor,
                // but only if still allowed.
                if(this.SubmitCount >= this.RetryCount) {
                    if(WriteHints) {
                        Console.WriteLine($"Note: Job has reached its maximum number of attempts to run ({this.RetryCount}) -- job is marked as fail, no further action.");
                        Console.WriteLine($"Hint: you might either remove old deployments or increase the 'RetryCount'.");
                    }
                    this.statusCache = JobStatus.FailedOrCanceled;
                    return JobStatus.FailedOrCanceled;
                }

                // ============
                // what now?
                // ============

                return JobStatus.Unknown;
            }
        }


        int m_RetryCount = 1;

        /// <summary>
        /// Number of times the job is submitted at maximum.
        /// </summary>
        public int RetryCount {
            get {
                return m_RetryCount;
            }
            set {
                if(RetryCount < 0)
                    throw new ArgumentOutOfRangeException("RetryCount must be positive");
                TestActivation();

                m_RetryCount = value;
            }
        }

        private void TestActivation() {
            if(this.Status != JobStatus.PreActivation)
                throw new NotSupportedException("Unable to change properties after job is activated.");
        }

        /// <summary>
        /// Human-readable summary about this job.
        /// </summary>
        public override string ToString() {
            using (var stw = new StringWriter()) {
                stw.Write(this.Name);
                stw.Write(": ");
                stw.Write(this.Status);

                if (AssignedBatchProc != null) {
                    stw.Write(" (");
                    stw.Write(this.AssignedBatchProc.ToString());
                    stw.Write(")");
                }

                return stw.ToString();
            }
        }

        string m_ExecutionTime = "00:05:00";

        /// <summary>
        /// Estimated execution time limit. Important for slurm queuing
        /// </summary>
        public string ExecutionTime {
            get {
                return m_ExecutionTime;
            }
            set {
                TestActivation();
                m_ExecutionTime = value;
            }
        }

        
        /// <summary>
        /// Copies the executable files to the <see cref="DeploymentBaseDirectory"/>, 
        /// but does not submit the job.
        /// </summary>
        string DeployExecuteables() {

            bool IsNotSystemAssembly(Assembly Ass, string MainAssemblyDir) {
                PlatformID CurrentSys = System.Environment.OSVersion.Platform;
                switch(CurrentSys) {
                    case PlatformID.Unix: { return Path.GetFullPath(Ass.Location).StartsWith("/home/"); }
                    case PlatformID.Win32S:
                    case PlatformID.Win32Windows:
                    default: {
                        return Path.GetDirectoryName(Ass.Location).Equals(MainAssemblyDir)
                            || Path.GetFileName(Ass.Location).StartsWith("BoSSS")
                            || Path.GetFileName(Ass.Location).StartsWith("ilPSP")
                            || !Ass.GlobalAssemblyCache;
                    }
                }
            }

            void TestWR() {
                using(new FuncTrace()) {
                    Exception OP(int iTry) {
                        if(this.GetControl() != null) {
                            var directories = this.AllDeployments.Select(dep => dep.DeploymentDirectory).Where(dep => dep != null).ToArray();
                            if(directories == null || directories.Length <= 0) {
                                return new IOException("Job is assigned to batch processor, but no deployment directory can be found.");
                            }
                        }
                        return null;
                    }

                    MetaJobMgrIO.RetryIOop(OP, "testing of job deployment", false);
                }
            }


            using (var tr = new FuncTrace()) {
                Console.WriteLine("Deploying executables and additional files ...");

                // Collect files
                List<string> files = new List<string>();
                using (new BlockTrace("ASSEMBLY_COLLECTION", tr)) {
                    //string SystemPath = Path.GetDirectoryName(typeof(object).Assembly.Location);
                    string MainAssemblyDir = Path.GetDirectoryName(EntryAssembly.Location);
                    foreach (var a in AllDependentAssemblies) {
                        if (IsNotSystemAssembly(a, MainAssemblyDir)) {
                            files.Add(a.Location);
                        }
                    }
               
                    // test for really strange errors
                    for (int i = 0; i < files.Count; i++) {
                        for (int j = i + 1; j < files.Count; j++) {
                            if (Path.GetFileName(files[i]).Equals(Path.GetFileName(files[j])))
                                throw new ApplicationException("strange internal error");
                        }
                    }
                }

                // create deployment directory.
                string DeployDir = this.GetNewDeploymentDir();
                MetaJobMgrIO.CreateDirectoryWR(DeployDir);

                Console.WriteLine("Deployment directory: " + DeployDir);
                string OriginDir = null;

                // copy files
                using (new BlockTrace("EXECOPY", tr)) {
                    foreach (var fOrg in files) {
                        string fNmn = Path.GetFileName(fOrg);
                        if (fNmn.Equals("mscorlib.dll"))
                            throw new ApplicationException("internal error - something went wrong during filtering.");

                        string fTarget = Path.Combine(DeployDir, fNmn);

                        MetaJobMgrIO.CopyFileWR(fOrg, fTarget);
                        if (OriginDir == null || !OriginDir.Equals(Path.GetDirectoryName(fOrg))) {
                            OriginDir = Path.GetDirectoryName(fOrg);
                        }
                    }
                }
                Console.WriteLine("copied " + files.Count + " files.");

                // additional files
                if (this.AdditionalDeploymentFiles != null) {
                    foreach (var t in this.AdditionalDeploymentFiles) {
                        string fTarget = Path.Combine(DeployDir, t.Item2);
                        byte[] Content = t.Item1;
                        MetaJobMgrIO.WriteFileWR(fTarget, Content);
                        Console.WriteLine("   written file: " + t.Item2);
                    }
                }

                // deploy runtime
                using (new BlockTrace("DEPLOY_RUNTIME", tr)) {
                    if (AssignedBatchProc.DeployRuntime) {
                        string BosssInstall = BoSSS.Foundation.IO.Utils.GetBoSSSInstallDir();
                        string BosssBinNative = Path.Combine(BosssInstall, "bin", Path.Combine("native", "win"));
                        MetaJobMgrIO.CopyDirectoryRec(BosssBinNative, DeployDir, "amd64");
                        Console.WriteLine("   copied 'amd64' runtime.");
                    }
                }

                // finally
                Console.WriteLine("deployment finished.");

                // test
                TestWR();

                // return
                return DeployDir;
            }
        }
    }

}

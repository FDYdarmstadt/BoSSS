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

        

        /// <summary>
        /// (Optional) object used by some batch processor (after calling <see cref="BatchProcessorClient.Submit(Job)"/>)
        /// in order to identify the job.
        /// </summary>
        public string BatchProcessorIdentifierToken {
            private set;
            get;
        }

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

        /*
        /// <summary>
        /// Specifies command line arguments for application startup; overrides any startup arguments (<see cref="CommandLineArguments"/>)
        /// set so far.
        /// </summary>
        public void SetControlCode(string code) {
            AddTextFile(code, "control.cs");

            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(PrjName)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }

            string[] args = new string[] {
                "--control", "control.cs",
                "--prjnmn", PrjName,
                "--sesnmn", this.Name
            };

            for (int i = 0; i < args.Length; i++) {
                m_EnvironmentVars.Add("BOSSS_ARG_" + i, args[i]);
            }
        }
        */


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
                    m_ctrl.SetDatabase(newDb);
                    ctrl_db = newDb;
                }
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
            /*
            if(ControlName == null)
                return null;

            var data = AdditionalDeploymentFiles.Single(tt => tt.Item2.Equals(ControlName));
            string ctrlfileContent = Encoding.UTF8.GetString(data.Item1);
            BoSSS.Solution.Control.AppControl ctrl;
            
            if (ControlName.EndsWith(".obj")) {

                ctrl = BoSSS.Solution.Control.AppControl.Deserialize(ctrlfileContent);

            } else if (ControlName.EndsWith(".cs")) {
                                
                BoSSS.Solution.Control.AppControl.FromCode(ctrlfileContent,m_ctrl.GetType(), out ctrl, out BoSSS.Solution.Control.AppControl[] ctrl_paramstudy);

                if (ctrl != null) {
                    // passt eh
                } else if (ctrl_paramstudy != null) {
                    

                    ctrl =  ctrl_paramstudy.ElementAt(m_ctrl_index);
                } else {
                    //throw new Exception(string.Format(
                    //    "Invalid control instruction: unable to cast the last"
                    //        + " result of the control file/cs-script of type {0} to type {1}",
                    //    controlObj.GetType().FullName,
                    //    typeof(T).FullName));
                    throw new NotSupportedException("unknown state.");
                }
            } else {
                throw new IOException("Unable to restore control object from given data.");
            }

            if(!m_ctrl.Equals(ctrl))
                throw new IOException("Control object mismatch after serialize/deserialize.");
                */


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
        /// The deployment directory, accessible from the local machine (e.g. a mounted path if the job is deployed to another computer);
        /// This should be an absolute path.
        /// </summary>
        public string DeploymentDirectory {
            get;
            private set;            
        }


        /// <summary>
        /// The last directory of <see cref="DeploymentDirectory"/>.
        /// </summary>
        public string RelativeDeploymentDirectory {
            get {
                var parts = DeploymentDirectory.Split(new[] { Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar }, StringSplitOptions.RemoveEmptyEntries);
                return parts.Last();
            }
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
                
                var AllCandidates = InteractiveShell.WorkflowMgm.Sessions.Where(
                    sinf => InteractiveShell.WorkflowMgm.SessionInfoJobCorrelation(sinf, this));

                var cnt = AllCandidates.Count();

                if (cnt <= 0)
                    return new ISessionInfo[0];


                return AllCandidates.ToArray();
            }
        }

        /// <summary>
        /// Alias for <see cref="AllSessions"/>
        /// </summary>
        public ISessionInfo[] GetAllSessions() {
            return AllSessions;
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
        /// The content of the standard output stream as a string.
        /// </summary>
        public string Stdout {
            get {
                if (AssignedBatchProc == null)
                    throw new NotSupportedException("Job is not activated.");

                string StdoutFile = AssignedBatchProc.GetStdoutFile(this);
                if (StdoutFile != null && File.Exists(StdoutFile)) {
                    using (FileStream stream = File.Open(StdoutFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)) {
                        using (StreamReader reader = new StreamReader(stream)) {
                            return reader.ReadToEnd();
                        }
                    }
                } else {
                    return "";
                }
            }
        }

        /// <summary>
        /// The content of the standard error stream as a string.
        /// </summary>
        public string Stderr {
            get {
                if (AssignedBatchProc == null)
                    throw new NotSupportedException("Job is not activated.");

                string StderrFile = AssignedBatchProc.GetStderrFile(this);
                if (StderrFile != null && File.Exists(StderrFile)) {
                    using (FileStream stream = File.Open(StderrFile, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)) {
                        using (StreamReader reader = new StreamReader(stream)) {
                            return reader.ReadToEnd();
                        }
                    }
                } else {
                    return "";
                }
            }
        }

        /// <summary>
        /// Should open a separate console window showing the rolling output text of the job.
        /// </summary>
        public void ShowOutput() {
            if (AssignedBatchProc == null)
                throw new NotSupportedException("Job is not activated.");

            string StderrFile = AssignedBatchProc.GetStderrFile(this);
            string StdoutFile = AssignedBatchProc.GetStdoutFile(this);

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
                if (statusCache == null) {
                    int SubmitCount;
                    string DD;
                    var s = GetStatus(out SubmitCount, out DD);

                    if(s == JobStatus.FinishedSuccessful) {
                        statusCache = s; // probably, not a lot things are happening status-wise with this job anymore
                    }

                    return s;
                } else {
                    return statusCache.Value;
                }
            }
        }

        /// <summary>
        /// How often has this job been submitted to the batch system?
        /// </summary>
        public int SubmitCount {
            get {
                int sc;
                string DD;
                GetStatus(out sc, out DD);
                return sc;
            }
        }

        private JobStatus GetStatus(out int SubmitCount, out string DD) {
            SubmitCount = 0;
            DD = null;
            if (AssignedBatchProc == null)
                return JobStatus.PreActivation;

            // test if session exists
            ISessionInfo[] RR = this.AllSessions;
            ISessionInfo R = RR.Length > 0 ? RR.OrderBy(si => si.CreationTime).Last() : null;

            if (RR.Any(si => !si.Tags.Contains(SessionInfo.NOT_TERMINATED_TAG)))
                return JobStatus.FinishedSuccessful;

            // find the deployment directory
            var directories = AssignedBatchProc.GetAllExistingDeployDirectories(this);
            if(this.DeploymentDirectory != null) {
                DirectoryInfo diAdd = new DirectoryInfo(this.DeploymentDirectory);
                if(diAdd.Exists) {
                    if(directories == null)
                        directories = new DirectoryInfo[0];

                    diAdd.AddToArray(ref directories);
                }
            }

            if(directories == null || directories.Length <= 0) {
                return JobStatus.PreActivation;
                //throw new IOException("Job is assigned to batch processor, but no deployment directory can be found.");
            }
            Array.Sort(directories, FuncComparerExtensions.ToComparer((DirectoryInfo a, DirectoryInfo b) => DateTime.Compare(a.CreationTime, b.CreationTime)));
            DirectoryInfo _DD = directories.Last();
            DD = _DD.FullName;
            SubmitCount = DD.Length;

            // Obtain token
            string bpcToken = this.BatchProcessorIdentifierToken;
            if(bpcToken.IsNullOrEmpty()) {
                try {
                    var l = File.ReadAllText(Path.Combine(DD, "IdentifierToken.txt"));
                    bpcToken = l.Trim();
                    this.BatchProcessorIdentifierToken = bpcToken;
                } catch(Exception) {
                    // job was probably deployed, but never submitted
                    // ignore this.
                }
            }
            
            if(bpcToken.IsEmptyOrWhite())
                return JobStatus.PreActivation;

            // ask the batch processor
            AssignedBatchProc.EvaluateStatus(bpcToken, DD, out bool isRunning, out bool isTerminated, out int ExitCode);
            bool isPending = (isRunning == false) && (isTerminated == false);

            if (isPending)
                return JobStatus.PendingInExecutionQueue;
            if (isRunning)
                return JobStatus.InProgress;

            if(isTerminated) {
                if(R != null && R.Tags.Contains(SessionInfo.NOT_TERMINATED_TAG))
                    return JobStatus.Failed;

                if(ExitCode != 0)
                    return JobStatus.Failed;

                return JobStatus.FinishedSuccessful;
            }

            throw new IOException("Unable to determine job status.");
        }

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
        public void Activate(BatchProcessorClient bpc) {

            // ============================================
            // ensure that this method is only called once.
            // ============================================
            if (this.AssignedBatchProc != null)
                throw new NotSupportedException("Job can only be activated once.");
            AssignedBatchProc = bpc;

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

                case JobStatus.Failed:
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


            // ========================================================================
            // finally, it might be necessary to submit the job to the batch processor. 
            // ========================================================================

            // deploy executables:
            bool RequiresDeploy = false;
            if (DeployDir == null) {
                if (!Directory.Exists(this.DeploymentDirectory) ||
                        !File.Exists(Path.Combine(DeployDir, Path.GetFileName(this.EntryAssembly.Location)))) {
                    RequiresDeploy = true;
                }
            } else {
                RequiresDeploy = false;
            }

            // some database syncing might be necessary 
            FiddleControlFile(bpc);

            // deploy additional files
            if (RequiresDeploy) {
                this.DeploymentDirectory = bpc.GetNewDeploymentDir(this);
                bpc.DeployExecuteables(this, AdditionalDeploymentFiles);
            }

            // submit job
            this.BatchProcessorIdentifierToken = bpc.Submit(this);
            File.WriteAllText(Path.Combine(this.DeploymentDirectory, "IdentifierToken.txt"), this.BatchProcessorIdentifierToken);
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

    }

}

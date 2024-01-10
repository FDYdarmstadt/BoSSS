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
using Microsoft.DotNet.Interactive.Formatting;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization.Formatters.Binary;
using System.Security.Cryptography.X509Certificates;
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
                } while(BoSSSshell.WorkflowMgm.AllJobs.ContainsKey(newName));
                
                Console.WriteLine($"Empty job name - picking new name '{newName}'");
                name = newName;
            }
            this.Solver = solver;
            this.Name = name;
            this.SessionReqForSuccess = true;
            
            if (BoSSSshell.WorkflowMgm.AllJobs.ContainsKey(name)) {
                throw new ArgumentException("Job with name '" + name + "' is already defined in the workflow management.");
            }
            BoSSSshell.WorkflowMgm.m_AllJobs.Add(name, this);
            if (string.IsNullOrWhiteSpace(BoSSSshell.WorkflowMgm.CurrentProject)) {
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
        
        /*
         * Fk, Anmerkung:
         * Sowas wie die folgenden Properties (MemPerCPU) sollten keine Eigenschaften des Job sein, weil es Scheduler-spezifisch ist.
         * Der ganze Quatsch soll in die BatchProcessorConfig.json:
         * ```
         * "AdditionalBatchCommands": [
         *          "#SBATCH -C avx512",
         *          "#SBATCH --mem-per-cpu=8000"
         * ]
         * ```
         * 

        ///// <summary>
        ///// The memory (in MB) that is reserved for every core
        ///// </summary>
        //public string MemPerCPU {
        //    set;
        //    get;
        //}

        
        private int m_NumberOfNodes = -1;

        /// <summary>
        /// overrides Memory per CPU criterion. MemoryPerCPU = memorypernode * nonode / cpupernode. Memory per node is architecture dependent.
        /// </summary>
        public int NumberOfNodes {
            get { return m_NumberOfNodes;  }
            set { m_NumberOfNodes = value;  }
        }
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
        /// path to the executable
        /// </summary>
        public string EntryAssemblyName {
            get {
                if(EntryAssemblyRedirection == null)
                    return Path.GetFileName(EntryAssembly.Location);
                else
                    return EntryAssemblyRedirection;
            }
        }

        /// <summary>
        /// Override to the assembly path;
        /// Note: typically, only used for the test-runners, where hundreds of jobs share the same assembly.
        /// </summary>
        public string EntryAssemblyRedirection {
            get;
            set;
        }


        /// <summary>
        /// Add dependent assemblies, may also contain stuff like mscorlib.dll.
        /// </summary>
        public IEnumerable<Assembly> AllDependentAssemblies {
            get {
                //HashSet<Assembly> assiList = new HashSet<Assembly>();
                //GetAllAssemblies(this.EntryAssembly, assiList, Path.GetDirectoryName(EntryAssembly.Location));
                var assiList = this.EntryAssembly.GetAllDependentAssemblies();
                return assiList.ToArray();
            }
        }


        Assembly[] GetRelevantAssemblies() {
            using (var tr = new FuncTrace()) {
                var files = new List<Assembly>();

                var allAssis = AllDependentAssemblies;
                //Debugger.Launch();
                Assembly mscorlib = allAssis.Single(a => a.Location.EndsWith("System.Runtime.dll"));
                string corLibPath = Path.GetDirectoryName(mscorlib.Location);


                string MainAssemblyDir = Path.GetDirectoryName(EntryAssembly.Location);
                tr.Info("MainAssemblyDir: " + MainAssemblyDir);
                tr.Info("mscorlib dir: " + corLibPath);

                //Debugger.Launch();

                foreach (var a in AllDependentAssemblies) {

                    // new rule for .NET6: if the file NOT located in the same directory as mscorlib.dll, it should be deployed;
                    // (in Jupyter, sometimes assemblies from some cache are used, therefore we cannot use the assembly location as a criterion)

                    if (Path.GetDirectoryName(a.Location) == corLibPath) {
                        tr.Info("ignoring (in corelib-path): " + a.Location);
                        continue;
                    }

                    files.Add(a);
                }

                return files.ToArray();
            }
        }

        /// <summary>
        /// All dependent assemblies which are not part of the dotnet SDK/runtime
        /// </summary>
        public IEnumerable<Assembly> RelevantDependentAssemblies {
            get {
                return GetRelevantAssemblies();

            }
        }


        /*
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
        internal static void GetAllAssemblies(Assembly a, HashSet<Assembly> assiList, string SearchPath) {
            if (assiList.Contains(a))
                return;
            assiList.Add(a);
            //if(a.FullName.Contains("codeanalysis", StringComparison.InvariantCultureIgnoreCase))
            //    Console.Write("");

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
        */

        List<Tuple<byte[], string>> m_AdditionalDeploymentFiles = new List<Tuple<byte[], string>>();

        /// <summary>
        /// Additional data files which will be deployed in the <see cref="Deployment.DeploymentDirectory"/> together with the assemblies.
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

            string PrjName = BoSSSshell.WorkflowMgm.CurrentProject;
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

            string PrjName = BoSSSshell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(PrjName)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }


            m_EnvironmentVars.Clear();
            for (int i = 0; i < Args.Length; i++) {
                m_EnvironmentVars.Add("BOSSS_ARG_" + i, Args[i]);
            }
        }


        /// <summary>
        /// Control, if set;
        /// </summary>
        public BoSSS.Solution.Control.AppControl Control {
            get {
                return m_ctrl;
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
            if(ctrl_db == null) {
                foreach(var db in BoSSSshell.wmg.AllDatabases) {
                    m_ctrl.SetDatabase(db);
                    if(bpc.IsDatabaseAllowed(m_ctrl)) {
                        Console.WriteLine("Set Database: " + db);
                        break;
                    }
                }
            }

            if(!bpc.IsDatabaseAllowed(m_ctrl)) {
                throw new IOException($"Database {ctrl_db} is not allowed for {this.ToString()}; You might either use a different database for this computation OR modify the 'AllowedDatabasesPaths' in '~/.BoSSS/etc/BatchProcessorConfig.json'.");
            }



            // check grid & restart info
            // =========================

            // grid function hack:
            if(m_ctrl.GridFunc != null) {
                Console.WriteLine("Control object contains grid function. Trying to Serialize the grid...");
                var dbi = m_ctrl.GetDatabase();
                if(dbi == null) {
                    throw new NotSupportedException("If a gird function is specified (instead of a grid id), a database must be specified to save the gird (when using the job manager).");
                }

                Foundation.Grid.IGrid g = m_ctrl.GridFunc();
                Guid id = dbi.SaveGrid(ref g);

                Console.WriteLine($"Using grid {id} at {dbi.Path}.");

                m_ctrl.GridFunc = null;
                m_ctrl.GridGuid = id;
                Console.WriteLine("Control object modified.");

            } else {
                if(m_ctrl.GridGuid != Guid.Empty) {
                    var db = m_ctrl.GetDatabase();
                    if(db != null) {
                        if(!db.Grids.Any(gi => gi.ID.Equals(m_ctrl.GridGuid))) {
                            Console.WriteLine("Grid is not in database yet...");
                            var grd = m_ctrl.m_Grid;
                            if(grd != null) {

                                Foundation.Grid.IGrid newInfo = db.SaveGrid(grd, true);
                                m_ctrl.SetGrid(newInfo);
                                Console.WriteLine("Grid successfully saved: " + newInfo.ID);
                            } else {
                                Console.WriteLine("Unable to save grid - a crash is likely.");
                            }
                        }
                    }
                }
            }

            if(ctrl_db != null) {
                if(!m_ctrl.GridGuid.Equals(Guid.Empty)) {

                    var GridIn_ctrl_db = ctrl_db.Grids.FirstOrDefault(GrdInf => GrdInf.ID.Equals(m_ctrl.GridGuid));

                    if(GridIn_ctrl_db == null) {
                        Console.WriteLine($"Grid {m_ctrl.GridGuid} is not present in database - copy to target system...");

                        var grid2copy = BoSSSshell.AllGrids.FirstOrDefault(dbGrid => dbGrid.ID.Equals(m_ctrl.GridGuid));
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

                        var sess_2copy = BoSSSshell.AllSessions.FirstOrDefault(sinf => sinf.ID.Equals(Rstsess_guid));
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
                Console.Error.WriteLine($"Warning: no database is set for the job to submit; nothing may be saved.");
            }

            // check
            // =====

            m_ctrl.VerifyEx();

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
        /// Specifies the control object for application startup; overrides any startup arguments (<see cref="CommandLineArguments"/>) set so far.
        /// The control object might be changed during <see cref="FiddleControlFile(BatchProcessorClient)"/>
        /// </summary>
        public void SetControlObject(BoSSS.Solution.Control.AppControl ctrl) {
            TestActivation();

            // serialize control object
            // ========================

            // Verification does not work before we execute `FiddleControlFile`
            //ctrl.VerifyEx();
            m_ctrl = ctrl;
            m_ctrl.ProjectName = BoSSSshell.WorkflowMgm.CurrentProject;

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
            string PrjName = BoSSSshell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(PrjName)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }

            string[] args = new string[] {
                "--control", "control.obj",
                "--prjnmn", PrjName,
                "--sesnmn", this.Name,
                "--num_threads", this.NumberOfThreads.ToString()
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
                m_EnvironmentVars["OMP_NUM_THREADS"] = this.NumberOfThreads.ToString();
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
                    if (DeploymentDirectory != null && DeploymentDirectory.Exists)
                        return DeploymentDirectory.CreationTime;
                    else
                        return new DateTime(1981, 7, 14);
                    
                }
            }


            string m_BatchProcessorIdentifierToken = null;

            /// <summary>
            /// Job ID uses by the respective batch processor
            /// </summary>
            public string BatchProcessorIdentifierToken {
                get {
                    if(m_BatchProcessorIdentifierToken.IsNullOrEmpty()) {
                        string DD = this.DeploymentDirectory?.FullName;
                        if (DD != null && Directory.Exists(DD)) {
                            try
                            {
                                var l = File.ReadAllText(Path.Combine(DD, "IdentifierToken.txt"));
                                m_BatchProcessorIdentifierToken = l.Trim();

                            }
                            catch (Exception)
                            {
                                // job was probably deployed, but never submitted
                                // ignore this.
                            }
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


            bool ReadExitCache(out JobStatus status, out int? ExitCode) {
                using (var tr = new FuncTrace()) {
                    ExitCode = null;
                    status = JobStatus.Unknown;

                    if (this.DeploymentDirectory == null)
                        return false;

                    string path = Path.Combine(this.DeploymentDirectory.FullName, "JobStatus_ExitCode.txt");
                    if (!File.Exists(path)) {
                        tr.Info("no job status cache file");
                        return false;
                    }

                    try {
                        using (var str = new StreamReader(path)) {
                            string l1 = str.ReadLine();
                            status = Enum.Parse<JobStatus>(l1);
                            tr.Info("status = " + ExitCode);
                            string l2 = str.ReadLine();
                            if (!l2.IsEmptyOrWhite()) {
                                ExitCode = int.Parse(l2);
                                tr.Info("Exit code = " + ExitCode);
                            }
                        }

                    } catch (Exception e) {
                        tr.Error($"{e.GetType()}: {e.Message}");
                        return false;
                    }

                    return true;
                }
            }

            /// <summary>
            /// If job deployment is finished (<see cref="JobStatus.FinishedSuccessful"/> or <see cref="JobStatus.FailedOrCanceled"/>),
            /// the <paramref name="status"/> and <paramref name="ExitCode"/> are stored in a file within the deployment directory.
            /// 
            /// This has two advantages:
            /// - less load for the job manager
            /// - most job managers forget jobs after a couple of days, so we better rememeber.
            /// </summary>
            void RememberCache(JobStatus status, int? ExitCode) {

                //Console.Error.WriteLine("Remembercache: " + ((this.DeploymentDirectory?.FullName)??"nix") + " exitsts? " + ((this.DeploymentDirectory?.Exists)??false));

                if (DeploymentDirectory != null && DeploymentDirectory.Exists) {
                    string path = Path.Combine(this.DeploymentDirectory.FullName, "JobStatus_ExitCode.txt");


                    try {
                        using (var str = new StreamWriter(path)) {
                            str.WriteLine(status.ToString());
                            if (ExitCode != null)
                                str.WriteLine(ExitCode.Value);
                            str.Flush();
                        }

                    } catch (Exception) {

                    }
                } 
            }



            /// <summary>
            /// Status of the deployment
            /// </summary>
            public JobStatus Status {
                get {
                    using(var tr = new FuncTrace()) {
                        tr.Info("Trying to get status of deployment: " + ((DeploymentDirectory?.FullName) ?? "no-path-avail"));
                        if (StatusCache != null) {
                            tr.Info("From chache: " + StatusCache.Value); 
                            return StatusCache.Value;
                        }

                        if (m_owner.AssignedBatchProc == null) {
                            tr.Info("No batch queue asigned: " + JobStatus.PreActivation);
                            return JobStatus.PreActivation;
                        }

                        JobStatus bpc_status;
                        int? ExitCode;
                        bool alreadyKnow = true;
                        string ExitCodeStr = null;
                        try {

                            alreadyKnow = ReadExitCache(out bpc_status, out ExitCode);
                            tr.Info($"From cache file: cached = {alreadyKnow}: {bpc_status}, exit code = {ExitCode}");
                            if (!alreadyKnow) {
                                if (this.BatchProcessorIdentifierToken == null) {
                                    if(this.Session != null) {
                                        if (this.Session.SuccessfulTermination) {
                                            bpc_status = JobStatus.FinishedSuccessful;
                                            ExitCode = 0;
                                        } else {
                                            bpc_status = JobStatus.Unknown;
                                            ExitCode = null;
                                        }
                                    } else {
                                        bpc_status = JobStatus.Unknown;
                                        ExitCode = null;
                                    }
                                } else {
                                    (bpc_status, ExitCode) = m_owner.AssignedBatchProc.EvaluateStatus(this.BatchProcessorIdentifierToken, this.optInfo, this.DeploymentDirectory?.FullName);
                                }
                            }
                            ExitCodeStr = ExitCode.HasValue ? ExitCode.Value.ToString() : "null";
                            this.ExitCodeCache = ExitCode;
                        } catch (Exception e) {
                           
                            tr.Error($"{e.GetType().Name} during Job.Deployment status evaluation: {e.Message}");
                            tr.Info("Exception trace: " + (e.StackTrace ?? ""));
                            bpc_status = JobStatus.Unknown;
                            ExitCode = null;
                        }
                        tr.Info("batch processor status: " + bpc_status);


                        if (bpc_status == JobStatus.FinishedSuccessful) {
                            StatusCache = bpc_status;
                            if (ExitCode == null || ExitCode != 0)
                                throw new ApplicationException($"Error in implementation of {m_owner.AssignedBatchProc}: job marked as {JobStatus.FinishedSuccessful}, but exit code is {ExitCodeStr} -- expecting 0.");
                            if (!alreadyKnow)
                                RememberCache(bpc_status, ExitCode);
                        }

                        if (bpc_status == JobStatus.FailedOrCanceled) {
                            StatusCache = bpc_status;
                            if (ExitCode == null || ExitCode == 0) {
                                tr.Info($"status is {bpc_status}, but exit code is {ExitCode}; resetting exit code to {-655321}");
                                ExitCode = -655321;
                                this.ExitCodeCache = ExitCode;
                            }
                            //    throw new ApplicationException($"Error in implementation of {m_owner.AssignedBatchProc}: job marked as {JobStatus.FailedOrCanceled}, but exit code is {ExitCodeStr} -- expecting any number except 0.");
                            if (!alreadyKnow)
                                RememberCache(bpc_status, ExitCode);
                        }

                        tr.Info($"Deployment: {bpc_status}, exit code = {ExitCodeCache}");
                        return bpc_status;
                    }
                }
            }

            /// <summary>
            /// If deployment terminated, the exit code of the app; otherwise, null
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
                return "Job token: "
                    + (this.BatchProcessorIdentifierToken ?? "unknown") + ", "
                    + this.Status
                    + " '" + (this.RelativeDeploymentDirectory ?? "unkown_DeploymentDirectory") + "'"
                    + (m_owner.AssignedBatchProc != null ? (" @ " + m_owner.AssignedBatchProc) : "");
            }
        }

        List<Deployment> m_Deployments = new List<Deployment>();


        string JobDirectoryBaseName() {
            string Exe = Path.GetFileNameWithoutExtension(EntryAssembly.Location);
            string Proj = BoSSSshell.WorkflowMgm.CurrentProject;
            string Sess = Name;

            return Proj
                //+ "-" + Sess 
                + "-" + Exe;
        }


        /// <summary>
        /// Creates a new directory where the assemblies for some job can be copied to for a deployment on this batch processor.
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
                string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss.ffffff");
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
                    if (dir1 == null && dir2 == null)
                        return true;
                    if (dir1 == null) // dir2 must be not null
                        return false;
                    if (dir2 == null) // dir1 must be not null
                        return false;
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
                try {
                    tr.Info("got " + AllDirs.Count() + " possible deployment directories: "
                        + AllDirs.Take(8).Select(dir => dir.FullName).ToConcatString("", ", ", AllDirs.Count() > 8 ? "..." : ""));
                } catch (Exception ex) {
                    tr.Warning("Exception during formatting of directory list: " + ex);
                }

                // filter appropriate ones 
                // =======================
                var filtDirs = new List<DirectoryInfo>();
                using(new BlockTrace("DIRECTORY_FILTERING", tr)) {
                    var mybind = new KnownTypesBinder(this);


                    foreach(DirectoryInfo dir in AllDirs) {
                        if(!DirIsKnown(dir)) {

                            try {
                                string ControlObj = Path.Combine(dir.FullName, "control.obj");
                                if(File.Exists(ControlObj)) {
                                    var ctrl = BoSSS.Solution.Control.AppControl.Deserialize(File.ReadAllText(ControlObj), mybind);
                                    if(BoSSSshell.WorkflowMgm.JobAppControlCorrelation(this, ctrl)) {
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
                                    if(BoSSSshell.WorkflowMgm.JobAppControlCorrelation(this, ctrl)) {
                                        filtDirs.Add(dir);
                                        continue;
                                    }
                                }
                            } catch (Exception e) {
                                tr.Error($"Warning: unable process deployment directory {dir}: " + e.Message);
                                Console.Error.WriteLine($"Warning: unable process deployment directory {dir}: " + e.Message);
                            }
                        }
                    }
                }

                try {
                    tr.Info("remaining with " + filtDirs.Count() + " deployment directories after filtering: "
                        + filtDirs.Take(8).Select(dir => dir.FullName).ToConcatString("", ", ", filtDirs.Count() > 8 ? "..." : ""));
                } catch (Exception ex) {
                    tr.Warning("Exception during formatting of directory list: " + ex);
                }

                // return
                // ======
                return filtDirs.ToArray();
            }
        }

        class KnownTypesBinder : Newtonsoft.Json.Serialization.DefaultSerializationBinder {

            Job m_owner;

            internal KnownTypesBinder(Job __owner) {
                m_owner = __owner;

                foreach(var a in __owner.AllDependentAssemblies) {
                    var tt = new Dictionary<string, Type>();
                    knownTypes.Add(a.GetName().Name, tt);
                    foreach(var t in a.GetExportedTypes()) {
                        tt.Add(t.FullName, t);
                    }
                }
            }

            Dictionary<string, Dictionary<string, Type>> knownTypes = new Dictionary<string, Dictionary<string, Type>>();

            /*
            public IList<Type> KnownTypes { get; set; }

            public Type BindToType(string assemblyName, string typeName) {
                return KnownTypes.SingleOrDefault(t => t.Name == typeName);
            }

            public void BindToName(Type serializedType, out string assemblyName, out string typeName) {
                assemblyName = null;
                typeName = serializedType.Name;
            }
            */
            public override Type BindToType(string assemblyName, string typeName) {
                var dd = knownTypes[assemblyName];
                var tt = dd[typeName];

                return tt;
            }
        }

        void UpdateDeployments() {
            using (var tr = new FuncTrace()) {
                

                // determine all new sessions
                // ==========================
                HashSet<Guid> KnownSessionGuids = new HashSet<Guid>();
                foreach (var dep in m_Deployments) {
                    if (dep.Session != null) {
                        KnownSessionGuids.Add(dep.Session.ID);
                    }
                }

               
                ISessionInfo[] AllNewSessionsTotal;
                using (new BlockTrace("NewSessionFiltering", tr)) {
                    AllNewSessionsTotal = BoSSSshell.WorkflowMgm.Sessions
                        .Where(sinf => !KnownSessionGuids.Contains(sinf.ID)).ToArray(); // for performance reasons, filter sessions that we already know
                }

                ISessionInfo[] AllNewSessions;
                using (new BlockTrace("SessionInfoJobCorrelation", tr)) {
                    AllNewSessions = AllNewSessionsTotal
                        .Where(sinf => BoSSSshell.WorkflowMgm.SessionInfoJobCorrelation(sinf, this)).ToArray();
                }

                // add all new deployment directories
                // ==================================
                foreach (DirectoryInfo dirInfo in this.GetAllUnkonwnExistingDeployDirectories()) {
                    m_Deployments.Add(new Deployment(dirInfo, this));
                }

                // finally, correlate the new sessions to the deployments
                // ======================================================
                for (int i = 0; i < AllNewSessions.Length; i++) {
                    var sinf = AllNewSessions[i];

                    if (sinf.DeployPath != null) {
                        string RelDeployPath = sinf.DeployPath.Split(new[] { Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar }, StringSplitOptions.RemoveEmptyEntries).Last();

                        foreach (var dep in m_Deployments) {
                            if (dep.RelativeDeploymentDirectory == RelDeployPath) {
                                // bingo
                                dep.Session = sinf;
                                AllNewSessions[i] = null;
                            }
                        }
                    }
                }

                for (int i = 0; i < AllNewSessions.Length; i++) {
                    var sinf = AllNewSessions[i];
                    if (sinf != null) {
                        // deployment was deleted, but some session was found in the database.
                        m_Deployments.Add(new Deployment(sinf, this));
                    }
                }

                m_Deployments.Sort(new FuncComparer<Deployment>((A, B) => A.CreationDate.CompareTo(B.CreationDate)));
            }
        }

        
        internal class DeploymentsAtomic : IDisposable {
            public DeploymentsAtomic(Job owner) {
                m_owner = owner;
                m_owner.UpdateDeployments();
            }

            readonly Job m_owner;

            public void Dispose() {
                m_owner.deploymentsAtomic = null;
            }
        }

        DeploymentsAtomic deploymentsAtomic = null;

        /// <summary>
        /// All deployments/attempts for this job which can be found in the data sources
        /// (sessions in BoSSS databases, deployment directories and batch processor queues)
        /// </summary>
        public IReadOnlyList<Deployment> AllDeployments {
            get {
                if (deploymentsAtomic != null) {
                    return m_Deployments.AsReadOnly();
                } else {
                    //maybe some caching here?
                    UpdateDeployments();
                    return m_Deployments.AsReadOnly();

                }
            }
        }


        /// <summary>
        /// Deletes any trace of previous calculations for this job.
        /// </summary>
        public void DeleteOldDeploymentsAndSessions(bool deleteSessions = true) {
            foreach(var dep in AllDeployments) {
                if(dep.Session != null) {
                    try {
                        dep.Session.Delete(true);
                    } catch (Exception e) {
                        Console.Error.WriteLine($"{e.GetType()} during deployment / session deletion: {e.Message}");
                    }

                }

                if(dep.DeploymentDirectory != null && dep.DeploymentDirectory.Exists) {
                    Exception op(int iTry) {
                        dep.DeploymentDirectory.Delete(true);
                        return null;
                    }


                    MetaJobMgrIO.RetryIOop(op, "deletion of directory '" + dep.DeploymentDirectory.FullName + "'", false);
                    Console.WriteLine($"Deleted deployment {dep.DeploymentDirectory.Name}");
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
                if (value <= 0)
                    throw new ArgumentOutOfRangeException("number of MPI processes must be at least 1");
                TestActivation();
                m_NumberOfMPIProcs = value;
            }
        }

        int m_NumberOfThreads = 4;

        /// <summary>
        /// Number of threads for each MPI rank
        /// </summary>
        public int NumberOfThreads {
            get {
                return m_NumberOfThreads;
            }
            set {
                if(value <= 0) 
                    throw new ArgumentOutOfRangeException("number of threads must be at least 1");
                TestActivation();
                m_NumberOfThreads = value;
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
                using (new FuncTrace()) {
                    List<ISessionInfo> ret = new List<ISessionInfo>();
                    foreach (var dep in AllDeployments) {
                        if (dep.Session != null)
                            ret.Add(dep.Session);
                    }
                    return ret.ToArray();
                }
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

            string EscapeKack(string p) {
                if (p.Contains(' ')) {
                    // i hate escaping
                    if (System.IO.Path.DirectorySeparatorChar == '\\') {
                        // probably windows --> use quotes
                        p = "\"" + p + "\"";
                    } else {
                        // Linux etc.
                        p = p.Replace(" ", "\\ ");
                    }

                }
                return p;
            }

            string StderrFile = AssignedBatchProc.GetStderrFile(ld.BatchProcessorIdentifierToken, ld.DeploymentDirectory.FullName);
            string StdoutFile = AssignedBatchProc.GetStdoutFile(ld.BatchProcessorIdentifierToken, ld.DeploymentDirectory.FullName);

            if (StdoutFile != null && StderrFile != null) {
                ProcessStartInfo psi = new ProcessStartInfo();
                psi.UseShellExecute = true;
                psi.WindowStyle = ProcessWindowStyle.Normal;
                psi.FileName = "dotnet";
                psi.Arguments = EscapeKack(typeof(btail.TailMain).Assembly.Location);
                btail.TailMain.SetArgs(psi, StdoutFile, StderrFile);

                Console.WriteLine($"Starting external console ...");
                Console.WriteLine("(You may close the new window at any time, the job will continue.)");

                Process p = Process.Start(psi);
            } else {
                Console.WriteLine("No known standard output/error path yet - try again later.");
            }
        }

        /// <summary>
        /// After calling <see cref="BatchProcessorClient.Submit"/>, this job
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

        /*
        /// <summary>
        /// triggers <see cref="DeleteOldDeploymentsAndSessions"/>
        /// </summary>
        public static bool UndocumentedSuperHack = false;
        */

        /// <summary>
        /// Activates the Job in the default queue (<see cref="BoSSSshell.GetDefaultQueue"/>)
        /// </summary>
        public void Activate() {
            var bpc = BoSSSshell.GetDefaultQueue();
            this.Activate(bpc);
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
        /// <param name="bpc"></param>
        public void Activate(BatchProcessorClient bpc) {
            using (var tr = new FuncTrace()) {

                // ============================================
                // ensure that this method is only called once.
                // ============================================
                if (this.AssignedBatchProc != null)
                    throw new NotSupportedException("Job can only be activated once.");
                AssignedBatchProc = bpc;

                Reactivate();
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public void Reactivate() {
            using (var tr = new FuncTrace()) {
                if (this.AssignedBatchProc == null)
                    throw new NotSupportedException("Job must be activated before.");

                // ================
                // status
                // ================

                var stat = GetStatus(true);
                if (stat != JobStatus.Unknown) {
                    int sc = this.SubmitCount;
                    if (stat == JobStatus.FailedOrCanceled && sc < RetryCount) {
                        Console.WriteLine($"Job is {stat}, but retry count is set to {RetryCount} and only {sc} tries yet - trying once more...");
                        this.statusCache = null;
                    } else {
                        Console.WriteLine("No submission, because job status is: " + stat.ToString());
                        return;
                    }
                }
             



                // ========================================================================
                // finally, it might be necessary to submit the job to the batch processor. 
                // ========================================================================

                Console.WriteLine($"Deploying job {this.Name} ... ");

                // some database syncing might be necessary 
                FiddleControlFile(AssignedBatchProc);

                // deploy additional files
                string DeploymentDirectory = this.DeployExecuteables();
                
                // submit job
                using (new BlockTrace("JOB_SUBMISSION", tr)) {
                    var rr = AssignedBatchProc.Submit(this, DeploymentDirectory);
                    File.WriteAllText(Path.Combine(DeploymentDirectory, "IdentifierToken.txt"), rr.id);

                    //Deployment dep = AllDeployments.SingleOrDefault(d => d?.BatchProcessorIdentifierToken == rr.id);
                    Deployment dep = AllDeployments.LastOrDefault(d => (d?.BatchProcessorIdentifierToken == rr.id) && PathMatch(d?.DeploymentDirectory?.FullName, DeploymentDirectory));
                    
                    if (dep == null)
                        m_Deployments.Add(new Deployment(new DirectoryInfo(DeploymentDirectory), this, rr.optJobObj));
                    else
                        dep.optInfo = rr.optJobObj;
                }

                Console.WriteLine();
            }
        }

        static bool PathMatch(string this_Path, string otherPath) {
            if (this_Path == otherPath)
                return true;
            if (this_Path == null && otherPath != null)
                return false;
            if (this_Path != null && otherPath == null)
                return false;

            if (!Directory.Exists(otherPath))
                return false;

            string TokenName = Guid.NewGuid().ToString() + ".token";

            string file1 = System.IO.Path.Combine(this_Path, TokenName);
            File.WriteAllText(file1, "this is a test file which can be safely deleted.");

            string file2 = System.IO.Path.Combine(otherPath, TokenName);

            return File.Exists(file2);
        }


        /// <summary>
        /// Status evaluation, with optional additional information.
        /// </summary>
        /// <param name="WriteHints">
        /// If true, prints a reason for the returned status.
        /// </param>
        /// <returns></returns>
        public JobStatus GetStatus(bool WriteHints = true) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = WriteHints;

                if (AssignedBatchProc == null) {
                    tr.Info($"No batch processor assigned to job jet; therefore, status is {JobStatus.PreActivation}");
                    return JobStatus.PreActivation;
                }

                if (!WriteHints) {
                    if (statusCache.HasValue) {
                        tr.Info("returning cached value: " + statusCache.Value);
                        return statusCache.Value;
                    }
                }

                using (new DeploymentsAtomic(this)) {
                    using (BoSSSshell.WorkflowMgm.EnterSessionAtomic()) {



                        // ================
                        // identify success
                        // ================

                        // we have to evaluate the status NOW, and work with this status in order for this method to work correctly.
                        // otherwise, the loophole describes below might happen!
                        (Deployment Depl, JobStatus fixedStatus)[] DeploymentsSoFar = this.AllDeployments.Select(dep => (dep, dep.Status)).ToArray();
                        tr.Info("Deployments so far (" + DeploymentsSoFar.Length + "): " + DeploymentsSoFar.ToConcatString("", ", ", ";"));

                        Deployment[] Success = DeploymentsSoFar.Where(dep => dep.fixedStatus == JobStatus.FinishedSuccessful).Select(TT => TT.Depl).ToArray();
                        tr.Info("Success: " + Success.Length);
                        if (SessionReqForSuccess) {
                            ISessionInfo[] SuccessSessions = this.AllSessions.Where(si => si.SuccessfulTermination()).OrderByDescending(sess => sess.CreationTime).ToArray();
                            if (SuccessSessions.Length <= 0) {
                                // look twice
                                BoSSSshell.WorkflowMgm.ResetSessionsCache();
                                DeploymentsSoFar = this.AllDeployments.Select(dep => (dep, dep.Status)).ToArray();
                                Success = DeploymentsSoFar.Where(dep => dep.fixedStatus == JobStatus.FinishedSuccessful).Select(TT => TT.Depl).ToArray();
                                SuccessSessions = this.AllSessions.Where(si => si.SuccessfulTermination()).OrderByDescending(sess => sess.CreationTime).ToArray();
                            }
                            if (SuccessSessions.Length > 0) {
                                tr.Info($"Info: Found successful session \"{SuccessSessions.First()}\" -- job is marked as successful, no further action.");
                                this.statusCache = JobStatus.FinishedSuccessful;
                                return JobStatus.FinishedSuccessful;
                            }

                            if (Success.Length > 0) {
                                tr.Info($"Note: found {Success.Length} successful deployment(s), but job is configured to require a successful result session ('this.SessionReqForSuccess' is true), and none is found. {this.AllSessions.Length} sessions correlated to this job fount in total.");
                            }

                            if (WriteHints) {
                                var c = this.GetControl();
                                if (c != null && c.savetodb == false) {
                                    throw new NotSupportedException("Configuration error: this job can never be successful: a successful result session ('this.SessionReqForSuccess' is true) is required, but control object is configured not to save in the database.");
                                }
                            }

                        } else {
                            if (Success.Length > 0) {
                                tr.Info($"Info: Found successful deployment and no result session is expected ('this.SessionReqForSuccess' is false):  job is marked as successful, no further action.");
                                this.statusCache = JobStatus.FinishedSuccessful;
                                return JobStatus.FinishedSuccessful;
                            }
                        }

                        // ========================
                        // identify running/waiting
                        // ========================

                        // If we would not use the status evaluated at the top of this method, the 
                        // potential async loophole could happen: job is still 'PendingInExecutionQueue' here...
                        var inprog = DeploymentsSoFar.Where(dep => (dep.fixedStatus == JobStatus.InProgress)).ToArray();
                        //tr.Info("inprog length: " + inprog.Length);
                        if (inprog.Length > 0) {
                            tr.Info($"Info: Job {inprog[0].Depl.BatchProcessorIdentifierToken} is currently executed on {this.AssignedBatchProc} -- no further action.");
                            return JobStatus.InProgress;
                        }

                        // ... but moves to 'InProgress' somewhere in between here ...
                        var inq = DeploymentsSoFar.Where(dep => (dep.fixedStatus == JobStatus.PendingInExecutionQueue)).ToArray();
                        //tr.Info("inq length: " + inq.Length);
                        if (inq.Length > 0) {
                            tr.Info($"Info: Job {inq[0].Depl.BatchProcessorIdentifierToken} is currently waiting on {this.AssignedBatchProc} -- no further action.");
                            return JobStatus.PendingInExecutionQueue;
                        }

                        // ... so we reach this line and receive an error here

                        // =============
                        // identify fail
                        // =============

                        // if we pass this point, we want to submit to a batch processor,
                        // but only if still allowed.
                        tr.Info("job submit count: " + this.SubmitCount);
                        if (this.SubmitCount >= this.RetryCount) {
                            tr.Info($"Note: Job has reached its maximum number of attempts to run ({this.RetryCount}) -- job is marked as fail, no further action.");
                            tr.Info($"Hint: you might either remove old deployments or increase the 'RetryCount'.");

                            this.statusCache = JobStatus.FailedOrCanceled;
                            return JobStatus.FailedOrCanceled;
                        }

                        if (this.SubmitCount > 0 && DeploymentsSoFar.All(dep => dep.fixedStatus == JobStatus.FailedOrCanceled)) {
                            tr.Info($"Note: Job was deployed ({this.SubmitCount}) number of times, all failed; RetryCount is {this.RetryCount}, so try once more.");
                            tr.Info($"Hint: want to re-activate the job.");

                            this.statusCache = JobStatus.FailedOrCanceled;
                            return JobStatus.FailedOrCanceled;
                        }

                        // ============
                        // what now?
                        // ============

                        tr.Error("unable to determine job status - unknown");
                        return JobStatus.Unknown;
                    }
                }
            }
        }

        int m_RetryCount = 3;

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

        /*

        // Note: Fk, 09feb22: this is a slurm-specific issue; therefore it should not be in this class.
        // If it is job-specific, it must be implemented for all batch-processors!


        string m_ExecutionTime = "05:00:00";

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
        */

        
        /// <summary>
        /// Copies the executable files to the <see cref="BatchProcessorClient.DeploymentBaseDirectory"/>, but does not submit the job.
        /// </summary>
        string DeployExecuteables() {
            /*
            bool IsNotSystemAssembly(Assembly Ass, string MainAssemblyDir) {
                PlatformID CurrentSys = System.Environment.OSVersion.Platform;
                switch(CurrentSys) {
                    case PlatformID.Unix: { return Path.GetFullPath(Ass.Location).StartsWith("/home/")
                            || Path.GetFullPath(Ass.Location).StartsWith("/Jenkins/"); }
                    case PlatformID.Win32S:
                    case PlatformID.Win32Windows:
                    default: {
                        return Path.GetDirectoryName(Ass.Location).Equals(MainAssemblyDir)
                            || Path.GetFileName(Ass.Location).StartsWith("BoSSS")
                            || Path.GetFileName(Ass.Location).StartsWith("ilPSP");
                            //|| !Ass.GlobalAssemblyCache;
                    }
                }
            }
            */

           
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

            static void CopyFilesRecursively(string sourcePath, string targetPath) {
                //Now Create all of the directories
                foreach(string dirPath in Directory.GetDirectories(sourcePath, "*", SearchOption.AllDirectories)) {
                    Directory.CreateDirectory(dirPath.Replace(sourcePath, targetPath));
                }

                //Copy all the files & Replaces any files with the same name
                foreach(string newPath in Directory.GetFiles(sourcePath, "*.*", SearchOption.AllDirectories)) {
                    File.Copy(newPath, newPath.Replace(sourcePath, targetPath), true);
                }
            }


            using (var tr = new FuncTrace()) {
                Console.WriteLine("Deploying executables and additional files ...");

                // Collect files
                List<string> files = new List<string>();
                using (new BlockTrace("ASSEMBLY_COLLECTION", tr)) {
                    if(EntryAssemblyRedirection == null) {
                        //string SystemPath = Path.GetDirectoryName(typeof(object).Assembly.Location);
                        files.AddRange(GetManagedFileList());
                    } else {
                        Console.WriteLine("Skipping assembly file copy, since 'EntryAssemblyRedirection' = " + EntryAssemblyRedirection);
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

                    if(EntryAssemblyRedirection == null) {
                        // copy "runtimes" directory from .NET core/.NET5
                        string runtimes_Src = Path.Combine(Path.GetDirectoryName(EntryAssembly.Location), "runtimes");
                        string runtimes_Dst = Path.Combine(DeployDir, "runtimes");
                        if(Directory.Exists(runtimes_Src)) {
                            CopyFilesRecursively(runtimes_Src, runtimes_Dst);
                        }
                    }
                }
                tr.Info("copied " + files.Count + " files.");
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
                    if (AssignedBatchProc.DeployRuntime == true) {
                        string BosssInstall = BoSSS.Foundation.IO.Utils.GetBoSSSInstallDir();
                        var BosssBinNative = new DirectoryInfo(Path.Combine(BosssInstall, "bin", "native", AssignedBatchProc.RuntimeLocation));
                        MetaJobMgrIO.CopyDirectoryRec(BosssBinNative.Parent.FullName, DeployDir, BosssBinNative.Name);
                        Console.WriteLine("   copied '" + AssignedBatchProc.RuntimeLocation + "' runtime.");
                        tr.Info("   copied '" + AssignedBatchProc.RuntimeLocation + "' runtime.");
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

        private string[] GetManagedFileList() {
            using (var tr = new FuncTrace()) {
                List<string> files = new List<string>();

                var allAssis = RelevantDependentAssemblies;

                string entry_dir = Path.GetDirectoryName(EntryAssembly.Location);

                foreach (var a in allAssis) {

                    string a_location = a.Location;
                    if (a_location != entry_dir) {
                        //
                        // FK, 24jan23:
                        // Try to take the file from the same directory as the entry assembly:
                        //  * problem: when using Jupyter, sometimes different assemblies are loaded.
                        //    - e.g. ~\.dotnet\tools\.store\microsoft.dotnet-interactive\1.0.360602\microsoft.dotnet-interactive\1.0.360602\tools\net7.0\any\System.CodeDom.dll
                        //      instead of the `System.CodeDom.dll` present in the current directory.
                        //    - this causes (sometimes) problems when the deployed application should be started.
                        //  * solution: test, whether the assembly with the same name can be found "locally"; if yes, prefer this one.
                        //  * future: if all packages are done via nuget, this will maybe be not an issue anymore
                        //
                        string alt_location = Path.Combine(entry_dir, Path.GetFileName(a_location));
                        if (File.Exists(alt_location)) {

                            tr.Info($"Take {alt_location} instead of {a_location}");
                            a_location = alt_location;
                        }
                    }

                    files.Add(a_location);

                    var additionalFiles = new string[] {
                        Path.Combine(Path.GetDirectoryName(a_location),
                                     System.IO.Path.GetFileNameWithoutExtension(a_location) + ".deps.json"),
                        Path.Combine(Path.GetDirectoryName(a_location),
                                     System.IO.Path.GetFileNameWithoutExtension(a_location) + ".runtimeconfig.json"),
                        a_location + ".config" // probably obsolete?; only for the old .NET Framework
                    };


                    foreach (var a_acc in additionalFiles) {
                        tr.Info("additional file " + a_acc);
                        if (File.Exists(a_acc)) {
                            tr.Info(" added.");
                            files.Add(a_acc);
                        } else {
                            tr.Info(" missing.");
                        }
                    }
                }

                return files.ToArray();
            }
        }
    }

}

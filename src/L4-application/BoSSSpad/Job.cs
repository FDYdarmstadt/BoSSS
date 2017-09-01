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
using ilPSP.Utils;
using System;
using System.Collections.Generic;
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
        /// <param name="controlObj"></param>
        public Job(string name, Type solver) { 
            this.Solver = solver;
            this.Name = name;
            if(InteractiveShell.WorkflowMgm.AllJobs.ContainsKey(name)) {
                throw new ArgumentException("Job with name '" + name + "' is already defined in the workflow management.");
            }
            InteractiveShell.WorkflowMgm.AllJobs.Add(name, this);
            if (string.IsNullOrWhiteSpace(InteractiveShell.WorkflowMgm.CurrentProject)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }
        }

        /// <summary>
        /// The name for the job; note that it serves as a unique, persistent identifier within the work flow management.
        /// </summary>
        public string Name {
            private set;
            get;
        }

        /// <summary>
        /// (Optional) object used by some batch processor (after calling <see cref="BatchProcessorClient.Submit(Job)"/>)
        /// in order to identify the job.
        /// </summary>
        public object BatchProcessorIdentifierToken {
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
                List<Assembly> assiList = new List<Assembly>();
                GetAllAssemblies(this.EntryAssembly, assiList, Path.GetDirectoryName(EntryAssembly.Location));
                return assiList.AsReadOnly();
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
        private static void GetAllAssemblies(Assembly a, List<Assembly> assiList, string SearchPath) {
            if (assiList.Contains(a))
                return;

            assiList.Add(a);

            foreach (AssemblyName b in a.GetReferencedAssemblies()) {
                Assembly na;
                try {
                    na = Assembly.Load(b);
                } catch (FileNotFoundException) {
                    string[] AssiFiles = ArrayTools.Cat(Directory.GetFiles(SearchPath, b.Name + ".dll"), Directory.GetFiles(SearchPath, b.Name + ".exe"));
                    if (AssiFiles.Length != 1)
                        throw new FileNotFoundException("Unable to locate assembly '" + b.Name + "'.");
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
        public void SetCommandLineArguments(string Args) {
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

            for (int i = 0; i < args.Length; i++) {
                m_EnvironmentVars.Add( "BOSSS_ARG_" + i, args[i]);
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
        public void MySetCommandLineArguments(string Args) {
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(PrjName)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }

            // we pass the startup-arguments through environment variables, which is 
            // (a bit) more robust (with respect to escape-characters, etc.) than 
            // command line arguments.
            string[] args = new string[] {
                Args
            };

            for (int i = 0; i < args.Length; i++) {
                m_EnvironmentVars.Add("BOSSS_ARG_" + i, args[i]);
            }
        }

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

        /// <summary>
        /// Specifies the control object for application startup; overrides any startup arguments (<see cref="CommandLineArguments"/>)
        /// set so far.
        /// </summary>
        public void SetControlObject(BoSSS.Solution.Control.AppControl ctrl) {
            byte[] buffer;
            using (var ms = new MemoryStream()) {
                var bf = new BinaryFormatter();
                bf.Serialize(ms, ctrl);
                buffer = ms.GetBuffer();
            }

            AdditionalDeploymentFiles.Add(new Tuple<byte[], string>(buffer, "control.obj"));

            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(PrjName)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }

            string[] args = new string[] {
                "--control", "control.obj",
                "--prjnmn", PrjName,
                "--sesnmn", this.Name
            };

            for (int i = 0; i < args.Length; i++) {
                m_EnvironmentVars.Add("BOSSS_ARG_" + i, args[i]);
            }
        }

        /// <summary>
        /// Specifies command line arguments for application startup; overrides any startup arguments (<see cref="CommandLineArguments"/>)
        /// set so far.
        /// </summary>
        public void SetControlFile(string FileName) {
            File.ReadAllText(FileName);
            SetControlCode(File.ReadAllText(FileName));
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
        /// Where are the job files?
        /// </summary>
        public string DeploymentDirectory {
            get;
            private set;
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
                if (AssignedBatchProc != null)
                    throw new NotSupportedException("Job is activated - no further change of parameters is possible.");
                m_NumberOfMPIProcs = value;
            }
        }

        /// <summary>
        /// Returns the session which is the result of this job.
        /// </summary>
        public ISessionInfo[] AllSessions {
            get {
                //string ProjectName = InteractiveShell.WorkflowMgm.CurrentProject;
                
                var AllCandidates = InteractiveShell.WorkflowMgm.Sessions.Where(
                    sinf => sinf.KeysAndQueries.ContainsKey(BoSSS.Solution.Application.SESSIONNAME_KEY) 
                         && Convert.ToString(sinf.KeysAndQueries[BoSSS.Solution.Application.SESSIONNAME_KEY]).Equals(this.Name)
                    );

                var cnt = AllCandidates.Count();

                if (cnt <= 0)
                    return new ISessionInfo[0];

                
                return AllCandidates.ToArray();
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
        /// The content of the standard output stream as a string.
        /// </summary>
        public string Stdout {
            get {
                if (AssignedBatchProc == null)
                    throw new NotSupportedException("Job is not activated.");

                string StdoutFile = AssignedBatchProc.GetStdoutFile(this);
                if (StdoutFile != null) {
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
                if (StderrFile != null) {
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
        /// After calling <see cref="BatchProcessorClient.Submit(Job)"/>, this job
        /// is assigned to the respective batch processor, which is recorded in this member.
        /// </summary>
        public BatchProcessorClient AssignedBatchProc {
            get;
            private set;
        }

        /// <summary>
        /// Whats up with this job?
        /// </summary>
        public JobStatus Status {
            get {
                int SubmitCount;
                string DD;
                return GetStatus(out SubmitCount, out DD);
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

            bool isSubmitted, isRunning, wasSuccessful, isFailed;
            AssignedBatchProc.EvaluateStatus(this, out SubmitCount, out isRunning, out wasSuccessful, out isFailed, out DD);
            isSubmitted = SubmitCount > 0;
            bool isPending = isSubmitted && !(isRunning || wasSuccessful || isFailed);

            if (isPending)
                return JobStatus.PendingInExecutionQueue;

            if (isRunning)
                return JobStatus.InProgress;

            ISessionInfo[] RR = this.AllSessions;
            ISessionInfo R = RR.Length > 0 ? RR.OrderBy(si => si.CreationTime).Last() : null;

            //if (RR.Length == 0)
            //    // maybe finished, but no result is known.
            //    return JobStatus.PreActivation;

            if (wasSuccessful || RR.Any(si => !si.Tags.Contains(BoSSS.Solution.Application.NOT_TERMINATED_TAG)))
                return JobStatus.FinishedSuccessful;

            if (isSubmitted && !(isFailed || wasSuccessful) && (R == null))
                return JobStatus.PendingInExecutionQueue;

            if (isFailed || (R == null || R.Tags.Contains(BoSSS.Solution.Application.NOT_TERMINATED_TAG)))
                return JobStatus.Failed;


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
                if(RR.Length <= 0) {
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
            if (DeployDir != null) {
                if (!Directory.Exists(this.DeploymentDirectory) ||
                        !File.Exists(Path.Combine(DeployDir, Path.GetFileName(this.EntryAssembly.Location)))) {
                    RequiresDeploy = true;
                }
            } else {
                RequiresDeploy = true;
            }

            if (RequiresDeploy) {
                this.DeploymentDirectory = bpc.GetNewDeploymentDir(this);
                bpc.DeployExecuteables(this, AdditionalDeploymentFiles);
            }

            // submit job
            this.BatchProcessorIdentifierToken = bpc.Submit(this);
        }

        int m_RetryCount = 1;

        /// <summary>
        /// Number of times the job is submitted at maximum.
        /// </summary>
        public int RetryCount {
            get {
                return m_RetryCount;
            }
        }

        /// <summary>
        /// Human-readable summary about this job.
        /// </summary>
        public string WriteInformation() {
            using (var stw = new StringWriter()) {


                return stw.ToString();
            }
        }

    }

}

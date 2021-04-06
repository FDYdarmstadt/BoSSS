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
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace BoSSS.Application.BoSSSpad {
    /// <summary>
    /// Abstraction for all kind of batch systems (aka. job managers).
    /// </summary>
    [DataContract]
    abstract public class BatchProcessorClient {

        /// <summary>
        /// Base directory where the executables should be deployed,
        /// accessible from the local machine (e.g. a mounted path if the batch processor deploys on another computer system)
        /// </summary>
        [DataMember]
        public string DeploymentBaseDirectory {
            get;
            protected set;
        }

        /// <summary>
        /// True, if it is required to deploy the runtime on the destination system,
        /// in order to execute a job.
        /// </summary>
        [DataMember]
        public bool DeployRuntime {
            get;
            set;
        }

        /// <summary>
        /// If not null, specifies paths to databases which are accessible to the computer system 
        /// on which this batch processor submits its jobs.
        /// This triggers data synchronization on job submission, if e.g. grid or restart timestep
        /// are in some other database.
        /// </summary>
        [DataMember]
        public string[] AllowedDatabasesPaths;


        List<IDatabaseInfo> m_AllowedDatabases;

        /// <summary>
        /// See <see cref="AllowedDatabasesPaths"/>
        /// </summary>
        public IReadOnlyList<IDatabaseInfo> AllowedDatabases {
            get {
                if(m_AllowedDatabases == null) {
                    m_AllowedDatabases = new List<IDatabaseInfo>();
                    if(AllowedDatabasesPaths != null) {

                        // fill up with empty entries
                        m_AllowedDatabases.AddRange(AllowedDatabasesPaths.Select(path => default(IDatabaseInfo)));

                        // pass 1: search in already-known databases
                        for(int i = 0; i < AllowedDatabasesPaths.Length; i++) {
                            foreach(var db in InteractiveShell.databases) {
                                if(db.PathMatch(AllowedDatabasesPaths[i])) {
                                    // bingo
                                    m_AllowedDatabases[i] = db;
                                    break;
                                }
                            }
                        }

                        // pass 2: try to open other databases
                        for(int i = 0; i < AllowedDatabasesPaths.Length; i++) {
                            if(m_AllowedDatabases[i] == null) {
                                try {
                                    var db = InteractiveShell.OpenDatabase(AllowedDatabasesPaths[i]);
                                    m_AllowedDatabases[i] = db;
                                    //if(!InteractiveShell.databases.Any(idb => object.ReferenceEquals(idb, db))) {
                                    //    InteractiveShell.databases.Add(db);
                                    //    Console.WriteLine($"Note: adding database {db} specified for batch queue {this}");
                                    //}
                                } catch(Exception e) {
                                    Console.Error.WriteLine($"Unable to open 'allowed database' for {this.ToString()} at path {AllowedDatabasesPaths[i]}. Check configuration file 'BatchProcessorConfig.json'. ({e.GetType().Name} : {e.Message})");
                                    Console.Error.WriteLine($"{this.ToString()} will continue to work, but database synchronization on job submission might not work correctly.");
                                }
                            }
                        }

                        // remove empty entries
                        for(int i = 0; i < m_AllowedDatabases.Count; i++) {
                            if(m_AllowedDatabases[i] == null) {
                                m_AllowedDatabases.RemoveAt(i);
                                i--;
                            }
                        }
                    }
                }

                var ret = new List<IDatabaseInfo>();
                ret.AddRange(m_AllowedDatabases);

                // add any local database which might be acceptable
                /*
                if(InteractiveShell.databases != null) {
                    foreach(var db in InteractiveShell.databases) {

                        bool found = false;
                        foreach(var odb in m_AllowedDatabases) {
                            if(odb.Equals(db)) {
                                found = true;
                                break;
                            }
                        }
                        if(found)
                            continue;

                        if(!found) {
                            // 'db' is not in the list to return, but is it located on the local machine?

                            db.
                        }
                    }
                }
                */
                return ret.AsReadOnly();
            }
        }





        /// <summary>
        /// Submits the job to the batch system.
        /// </summary>
        /// <param name="myJob">Job to submit.</param>
        /// <param name="DeploymentDirectory">Where the executable is.</param>
        /// <returns>
        /// An identifier token (<see cref="Job.Deployment.BatchProcessorIdentifierToken"/>)
        /// as well as an optional (internal) object
        /// </returns>
        abstract public (string id, object optJobObj) Submit(Job myJob, string DeploymentDirectory);

        /// <summary>
        /// Try to get some information about a job from the job manager.
        /// </summary>
        /// <param name="idToken">Identification within batch processor</param>
        /// <param name="optInfo">
        /// Optional internal job object, returned form <see cref="Submit"/>
        /// </param>
        /// <param name="DeployDir"></param>
        public abstract (BoSSSpad.JobStatus, int? ExitCode) EvaluateStatus(string idToken, object optInfo, string DeployDir);

        /// <summary>
        /// Path to standard output file, if present - otherwise null.
        /// </summary>
        public abstract string GetStdoutFile(string idToken, string DeployDir);

        /// <summary>
        /// Path to standard error file, if present - otherwise null.
        /// </summary>
        public abstract string GetStderrFile(string idToken, string DeployDir);


    }

    /// <summary>
    /// File IO utilities for unreliable file-systems (perform 10x re-try if operation fails)
    /// </summary>
    public static class MetaJobMgrIO { 

        /// <summary>
        /// Generic IO operation with re-try (seems to be necessary when working with network file systems).
        /// </summary>
        internal static void RetryIOop(Func<int,Exception> op, string Message, bool SurpressException) {
            int MaxTry = 10;
            Random rnd = null;
            Exception latest = null;
            for(int iTry = 0; iTry < MaxTry; iTry++) {
                // on network file-systems, there seem to be some rare hiccups, sometimes:
                // hundreds of files copied successfully, suddenly an IOException: file already exists.
                // File indeed exists, but is empty -- makes no sense, since deploy directory is freshly created.

                try  {

                    latest = op(iTry);
                    if(latest != null)
                        break;
                    else {
                        if(iTry > 0)
                            Console.WriteLine("success.");
                        return;
                    }
                } catch(IOException e) {
                    Console.Error.WriteLine(e.GetType().Name + " " + Message +  " : " + e.Message);
                    
                    Console.WriteLine($"Retrying {iTry + 1} of {MaxTry} (waiting for some time before) ...");
                    if(rnd == null)
                        rnd = new Random();
                    latest = e;
                    int iwait = rnd.Next(77 * 1000);
                    System.Threading.Thread.Sleep(iwait); // sleep for at most 77 seconds...
                }
            }

            if(SurpressException == false && latest != null)
                throw latest;
        }


        /// <summary>
        /// File write with re-try (seems to be necessary when working with network file systems).
        /// </summary>
        internal static void WriteFileWR(string fTarget, byte[] Content) {

            Exception OP(int iTry) {
                if(iTry == 0 && File.Exists(fTarget)) {
                    return new IOException("File '" + fTarget + "' already exists - wont over write.");
                }
                File.WriteAllBytes(fTarget, Content);
                return null;
            }

            RetryIOop(OP, "writing file '" + fTarget + "'", false);

          
        }

        /// <summary>
        /// File copy with re-try (seems to be necessary when working with network file systems).
        /// </summary>
        internal static void CopyFileWR(string fOrg, string fTarget, bool SurpressException = false) {
            
            Exception op(int iTry) {
                File.Copy(fOrg, fTarget, iTry > 0);
                return null;
            }

            RetryIOop(op, " copy of file '" + fOrg + "' --> '" + fTarget + "'", SurpressException);

       
        }

        /// <summary>
        /// File copy with re-try (seems to be necessary when working with network file systems).
        /// </summary>
        internal static void CreateDirectoryWR(string dstSubDir, bool SurpressException = false) {
            
            Exception op(int i) {
                if(!Directory.Exists(dstSubDir))
                    Directory.CreateDirectory(dstSubDir);
                return null;
            }

            RetryIOop(op, "creation of directory '" + dstSubDir + "'", SurpressException);

        }

        /// <summary>
        /// (tries to) do a recursive copy of a directory
        /// </summary>
        /// <param name="srcDir"></param>
        /// <param name="dstDir"></param>
        /// <param name="filter">search pattern/filter</param>
        public static void CopyDirectoryRec(string srcDir, string dstDir, string filter) {
            string[] srcFiles = Directory.GetFiles(srcDir);

            foreach (string srcFile in srcFiles) {
                //TryCopy(srcFile, Path.Combine(dstDir, Path.GetFileName(srcFile)));
                CopyFileWR(srcFile, Path.Combine(dstDir, Path.GetFileName(srcFile)), SurpressException:true);
            }

            string[] subDirs;
            if (filter == null)
                subDirs = Directory.GetDirectories(srcDir);
            else
                subDirs = Directory.GetDirectories(srcDir, filter);
            foreach (string srcAbsDir in subDirs) {
                string srcRelDir = Path.GetFileName(srcAbsDir);
                string dstSubDir = Path.Combine(dstDir, srcRelDir);
                CreateDirectoryWR(dstSubDir);
                CopyDirectoryRec(srcAbsDir, dstSubDir, null);
            }
        }

    }
}


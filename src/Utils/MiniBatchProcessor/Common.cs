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
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MiniBatchProcessor {

    /// <summary>
    /// Functions which are used by both, server and client.
    /// </summary>
    public abstract class ClientAndServer {

        static readonly Random rnd = new();

        /// <summary>
        /// ctor.
        /// </summary>
        internal ClientAndServer(string BatchInstructionDir) {
            this.config = new Configuration(BatchInstructionDir);
            this.config.Load();
        }

        /// <summary>
        /// %.
        /// </summary>
        public Configuration config;

        /// <summary>
        /// If an IO operation fails, a randomized time, in milliseconds, within a certain range, before the operation should be tried again.
        /// </summary>
        public static int IOwaitTime => rnd.Next(500, 1500);

        /// <summary>
        /// Maximum number of times an IO Operation is re-tryed.
        /// </summary>
        public const int IO_OPS_MAX_RETRY_COUNT = 5;

        /// <summary>
        /// Directory for jobs that are waiting.
        /// </summary>
        public const string QUEUE_DIR = "queue";

        /// <summary>
        /// Directory for jobs which are currently in progress, and finished ones.
        /// </summary>
        public const string WORK_FINISHED_DIR = "work";

        private void ReadWorkDir(Dictionary<int, Tuple<JobData, JobStatus>> R) {
            string dir = Path.Combine(this.config.BatchInstructionDir, WORK_FINISHED_DIR);
            if ( !Directory.Exists(dir) ) {
                _ = Directory.CreateDirectory(dir);
            }

            foreach ( string fName in Directory.GetFiles(dir, "*") ) {
                bool IsInt = int.TryParse(Path.GetFileName(fName), out int id);

                if ( !IsInt ) {
                    continue;
                }

                if ( !R.TryGetValue(id, out Tuple<JobData, JobStatus> JS) ) {

                    var Jdata = JobData.FromFile(Path.Combine(dir, fName));

                    // must be a new file
                    JS = Tuple.Create(
                        Jdata,
                        JobStatus.Undefined);
                    R.Add(id, JS);
                } else {
                    // job is already known
                    if ( JS.Item2 == JobStatus.Finished ) {
                        continue;
                    }
                }

                {
                    string exit_token = id + "_exit.txt";
                    if ( File.Exists(Path.Combine(dir, exit_token)) ) {
                        JS = Tuple.Create(JS.Item1, JobStatus.Finished);
                    } else {
                        if ( File.Exists(this.GetStdoutFile(id)) && File.Exists(this.GetStderrFile(id)) ) {
                            JS = Tuple.Create(JS.Item1, JobStatus.Working);
                        }
                    }
                }

                if ( JS.Item1.ID != id ) {
                    throw new ApplicationException("Mismatch between job id in file and file name.");
                }

                R[id] = JS;
            }
        }

        /// <summary>
        /// Path to a jobs standard output file.
        /// </summary>
        public string GetStdoutFile(int JobId) {
            string f = Path.Combine(this.config.BatchInstructionDir, ClientAndServer.WORK_FINISHED_DIR, JobId.ToString() + "_out.txt");
            return f;
        }

        /// <summary>
        /// Path to a jobs standard error file.
        /// </summary>
        public string GetStderrFile(int JobId) {
            string f = Path.Combine(this.config.BatchInstructionDir, ClientAndServer.WORK_FINISHED_DIR, JobId.ToString() + "_err.txt");
            return f;
        }

        private void ReadQueueDir(Dictionary<int, Tuple<JobData, JobStatus>> R) {
            string dir = Path.Combine(this.config.BatchInstructionDir, QUEUE_DIR);
            if ( !Directory.Exists(dir) ) {
                _ = Directory.CreateDirectory(dir);
            }

            foreach ( string fName in Directory.GetFiles(dir, "*") ) {
                bool IsInt = int.TryParse(Path.GetFileName(fName), out int id);

                if ( !IsInt ) {
                    continue;
                }

                if ( R.ContainsKey(id) ) {
                    continue;
                }

                var J = JobData.FromFile(Path.Combine(dir, fName));

                if ( J == null ) {
                    continue;
                }

                if ( J.ID != id ) {
                    //throw new ApplicationException("Mismatch between job id in file and file name.");
                    Console.Error.WriteLine("Mismatch between job id in file and file name.");
                    continue;
                }

                R.Add(J.ID, new Tuple<JobData, JobStatus>(J, JobStatus.Queued));
            }
        }

        void UpdateLists() {
            using ( new FuncTrace() ) {
                lock ( this.padlock_AllJobs ) {
                    var AllJobsNew = new Dictionary<int, Tuple<JobData, JobStatus>>();
                    AllJobsNew.AddRange(this.m_AllJobs);

                    this.ReadQueueDir(AllJobsNew);
                    this.ReadWorkDir(AllJobsNew);

                    //ReadDir(QUEUE_DIR, AllJobsNew, JobStatus.Queued);
                    //ReadDir(WORK_FINISHED_DIR, AllJobsNew, JobStatus.Working);

                    //foreach (var kv in AllJobsNew.ToArray()) {
                    //    var Ja = AllJobsNew[i];

                    //    for (int j = i + 1; j < AllJobsNew.Count; j++) {
                    //        var Jb = AllJobsNew[j];

                    //        if(Ja.Item1.ID == Jb.Item1.ID) {
                    //            Console.WriteLine("Warning: found job ID {0} more than once - ignoring second occurrence (first in list {1}, second in list {2}).", Ja.Item1, Ja.Item2, Jb.Item2);
                    //            AllJobsNew.RemoveAt(j);
                    //            j--;
                    //        }
                    //    }
                    //}

                    this.m_AllJobs = AllJobsNew;
                }
            }
        }

        /// <summary>
        /// %
        /// </summary>
        /// <param name="JobId">
        /// The jobs id, see <see cref="JobData.ID"/>.
        /// </param>
        /// <returns></returns>
        public (JobStatus stat, int ExitCode) GetStatusFromID(int JobId) {
            using ( new FuncTrace() ) {
                this.UpdateLists();
                Tuple<JobData, JobStatus> jd = null;
                lock ( this.padlock_AllJobs ) {
                    if ( !this.m_AllJobs.ContainsKey(JobId) ) {
                        return (JobStatus.Undefined, -1);
                    }

                    jd = this.m_AllJobs[JobId];
                }

                {
                    if ( jd.Item1.ID == JobId ) {
                        int ExitCode = 0;

                        if ( jd.Item2 == JobStatus.Finished ) {
                            try {
                                string ExitCodeFile = Path.Combine(this.config.BatchInstructionDir, WORK_FINISHED_DIR, JobId.ToString() + "_exit.txt");
                                using var exit = new StreamReader(ExitCodeFile);
                                ExitCode = int.Parse(exit.ReadLine());
                            } catch ( Exception ) {
                                ExitCode = -1;
                            }
                        }

                        return (jd.Item2, ExitCode);
                    }
                }

                throw new ArgumentException("Unable to find job with id '" + JobId + "'.");
            }
        }

        IDictionary<int, Tuple<JobData, JobStatus>> m_AllJobs;

        /// <summary>
        /// Exactly what it says.
        /// </summary>
        public IList<JobData> AllJobs {
            get {
                this.UpdateLists();
                lock ( this.padlock_AllJobs ) {
                    return this.m_AllJobs.Values.Select(job => job.Item1).ToList().AsReadOnly();
                }
            }
        }

        readonly object padlock_AllJobs = new();

        internal void UpdateStatus(int id, JobStatus newStat) {
            lock ( this.padlock_AllJobs ) {
                Tuple<JobData, JobStatus> JS = this.m_AllJobs[id];
                var nJS = Tuple.Create(JS.Item1, newStat);
                this.m_AllJobs[id] = nJS;
            }
        }

        /// <summary>
        /// Jobs which are in the waiting queue.
        /// </summary>
        public IList<JobData> Queue {
            get {
                this.UpdateLists();
                lock ( this.padlock_AllJobs ) {
                    return this.m_AllJobs.Values.Where(job => job.Item2 == JobStatus.Queued)
                        .Select(job => job.Item1)
                        .ToList().AsReadOnly();
                }
            }
        }

        /// <summary>
        /// Jobs which are currently being executed.
        /// </summary>
        public IList<JobData> Working {
            get {
                this.UpdateLists();
                lock ( this.padlock_AllJobs ) {
                    return this.m_AllJobs.Values.Where(job => job.Item2 == JobStatus.Working)
                        .Select(job => job.Item1)
                        .ToList().AsReadOnly();
                }
            }
        }

        /// <summary>
        /// Jobs which are Finished.
        /// </summary>
        public IList<JobData> Finished {
            get {
                this.UpdateLists();
                lock ( this.padlock_AllJobs ) {
                    return this.m_AllJobs.Values.Where(job => job.Item2 == JobStatus.Finished)
                        .Select(job => job.Item1)
                        .ToList().AsReadOnly();
                }
            }
        }
    }
}

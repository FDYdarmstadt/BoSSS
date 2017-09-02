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

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MiniBatchProcessor {

    /// <summary>
    /// Functions which are used by both, server and client.
    /// </summary>
    public static class ClientAndServer {

        static Random rnd;


        /// <summary>
        /// ctor.
        /// </summary>
        static ClientAndServer() {
            config = new Configuration();
            rnd = new Random();

        }
        
        /// <summary>
        /// %.
        /// </summary>
        public static Configuration config;

        /// <summary>
        /// If an IO operation fails, a randomized time, in milliseconds, within a certain range, before the operation should be tried again.
        /// </summary>
        static public int IOwaitTime {
            get {
                return 500 + (int)(Math.Round(rnd.NextDouble() * 500));
            }
        }
       

        /// <summary>
        /// Maximum number of times an IO Operation is re-tryed.
        /// </summary>
        public const int IO_OPS_MAX_RETRY_COUNT = 5;

        /// <summary>
        /// Directory for jobs that are waiting.
        /// </summary>
        public const string QUEUE_DIR = "queue";

        /// <summary>
        /// Directory for jobs which are currently in progress.
        /// </summary>
        public const string WORK_DIR = "work";

        /// <summary>
        /// Directory for finished jobs.
        /// </summary>
        public const string FINISHED_DIR = "finished";

        private static void ReadDir(string RelDir, List<Tuple<JobData, JobStatus>> R, JobStatus s) {
            string dir = Path.Combine(config.BatchInstructionDir, RelDir);

            foreach (var fName in Directory.GetFiles(dir, "*")) {
                int id;
                bool IsInt = Int32.TryParse(Path.GetFileName(fName), out id);

                if (!IsInt)
                    continue;

                var J = JobData.FromFile(Path.Combine(dir, fName));

                if (J == null)
                    continue;

                
                R.Add(new Tuple<JobData, JobStatus>(J, s));
            }
        }

        static void UpdateLists() {
            var AllJobsNew = new List<Tuple<JobData, JobStatus>>();
            ReadDir(QUEUE_DIR, AllJobsNew, JobStatus.Queued);
            ReadDir(WORK_DIR, AllJobsNew, JobStatus.Working);
            ReadDir(FINISHED_DIR, AllJobsNew, JobStatus.Finished);

            for(int i = 0; i < AllJobsNew.Count; i++) {
                var Ja = AllJobsNew[i];
                for (int j = i + 1; j < AllJobsNew.Count; j++) {
                    var Jb = AllJobsNew[j];

                    if(Ja.Item1.ID == Jb.Item1.ID) {
                        Console.WriteLine("Warning: found job ID {0} more than once - ignoring second occurrence (first in list {1}, second in list {2}).", Ja.Item1, Ja.Item2, Jb.Item2);
                        AllJobsNew.RemoveAt(j);
                        j--;
                    }
                }
            }

            m_AllJobs = AllJobsNew;
        }

        /// <summary>
        /// Returns a currently unused job id.
        /// </summary>
        public static int GetNewId() {
            UpdateLists();
            if (m_AllJobs.Count <= 0) {
                return 1;
            } else {
                int newId = m_AllJobs.Max(J => J.Item1.ID) + 1;
                return newId;
            }
        }


        /// <summary>
        /// %
        /// </summary>
        /// <param name="JobId">
        /// The jobs id, see <see cref="JobData.ID"/>.
        /// </param>
        /// <returns></returns>
        public static JobStatus GetStatusFromID(int JobId, out int ExitCode) {
            UpdateLists();
            foreach (var jd in m_AllJobs) {
                if(jd.Item1.ID == JobId) {
                    ExitCode = 0;

                    if(jd.Item2 == JobStatus.Finished) {
                        try {
                            string ExitCodeFile = Path.Combine(ClientAndServer.config.BatchInstructionDir, FINISHED_DIR, JobId.ToString() + "_exit.txt");
                            using (var exit = new StreamReader(ExitCodeFile)) {
                                ExitCode = int.Parse(exit.ReadLine());
                            }
                        } catch (Exception) {
                            ExitCode = -1;
                        }
                    }

                    return jd.Item2;
                }
            }
            throw new ArgumentException("Unable to find job with id '" + JobId + "'.");
        }

        


        static List<Tuple<JobData, JobStatus>> m_AllJobs;

        /// <summary>
        /// Exactly what it says.
        /// </summary>
        public static IEnumerable<JobData> AllJobs {
            get {
                UpdateLists();
                return m_AllJobs.Select(job => job.Item1).ToList().AsReadOnly();
            }
        }


        /// <summary>
        /// Jobs which are in the waiting queue.
        /// </summary>
        public static IEnumerable<JobData> Queue {
            get {
                UpdateLists();
                return m_AllJobs.Where(job => job.Item2 == JobStatus.Queued)
                    .Select(job => job.Item1)
                    .ToList().AsReadOnly();
            }
        }


        /// <summary>
        /// Jobs which are currently being executed.
        /// </summary>
        public static IEnumerable<JobData> Working {
            get {
                UpdateLists();
                return m_AllJobs.Where(job => job.Item2 == JobStatus.Working)
                    .Select(job => job.Item1)
                    .ToList().AsReadOnly();
            }
        }

        
        /// <summary>
        /// Jobs which are Finnish.
        /// </summary>
        public static IEnumerable<JobData> Finished {
            get {
                UpdateLists();
                return m_AllJobs.Where(job => job.Item2 == JobStatus.Finished)
                    .Select(job => job.Item1)
                    .ToList().AsReadOnly();
            }
        }
    }
}

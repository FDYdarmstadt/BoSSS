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

using BoSSS.Platform;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace MiniBatchProcessor {

    /// <summary>
    /// Exe dir, etc.
    /// </summary>
    public class JobData {

        internal int m_ID;

        /// <summary>
        /// Job identifier; unique among all currently known jobs.
        /// </summary>
        public int ID {
            get {
                return m_ID;
            }
        }

        /// <summary>
        /// Non-Unique job name.
        /// </summary>
        public string Name = "NoName";

        /// <summary>
        /// First line in the job file, 
        /// directory where the job should be executed
        /// </summary>
        public string ExeDir;

        /// <summary>
        /// Second line in the job file,
        /// number of MPI processes.
        /// </summary>
        public int NoOfProcs;

        /// <summary>
        /// Third line in the job file, program to execute
        /// </summary>
        public string exefile;

        /// <summary>
        /// Additional startup arguments.
        /// </summary>
        public string[] Arguments = new string[0];

        /// <summary>
        /// Additional environment variables.
        /// </summary>
        public Tuple<string, string>[] EnvVars = new Tuple<string, string>[0];

        /// <summary>
        /// Date at which the job file was written.
        /// </summary>
        public DateTime SubmitTime {
            get;
            private set;
        }


        /// <summary>
        /// Loads the job from some text file
        /// </summary>
        internal static JobData FromFile(string f) {
            int ReTryCount = 0;
            while (true) {
                try {
                    int id = int.Parse(Path.GetFileNameWithoutExtension(f));

                    using (var fStr = new StreamReader(new FileStream(f, FileMode.Open, FileAccess.Read, FileShare.None))) {
                        JobData J = new JobData();

                        J.m_ID = id;

                        J.Name = fStr.ReadLine();
                        J.ExeDir = fStr.ReadLine();
                        if (string.IsNullOrWhiteSpace(J.ExeDir)) {
                            J.ExeDir = null;
                        }
                        J.NoOfProcs = Convert.ToInt32(fStr.ReadLine());
                        J.exefile = fStr.ReadLine();

                        int NoOfArgs = int.Parse(fStr.ReadLine());
                        J.Arguments = new string[NoOfArgs];
                        for (int j = 0; j < NoOfArgs; j++) {
                            J.Arguments[j] = fStr.ReadLine();
                        }

                        int NoOfEnvVars = int.Parse(fStr.ReadLine());
                        J.EnvVars = new Tuple<string, string>[NoOfEnvVars];
                        for (int j = 0; j < NoOfEnvVars; j++) {
                            J.EnvVars[j] = new Tuple<string, string>(fStr.ReadLine(), fStr.ReadLine());
                        }

                        J.SubmitTime = File.GetCreationTime(f);

                        return J;
                    }
                } catch (Exception E) {
                    if (ReTryCount < ClientAndServer.IO_OPS_MAX_RETRY_COUNT) {
                        ReTryCount++;
                        Thread.Sleep(ClientAndServer.IOwaitTime);
                    } else {
                        Console.Error.WriteLine("Exception {0} reading job file '{1}': {2}", E.GetType().Name, f, E.Message);
                        return null;
                    }
                }
            }
        }

        


        /// <summary>
        /// Saves job data in queue-directory.
        /// </summary>
        internal void Save() {
            
            string RelDir = ClientAndServer.QUEUE_DIR;
            string f = Path.Combine(ClientAndServer.config.BatchInstructionDir, RelDir, this.ID.ToString());

            int ReTryCount = 0;
            while (true) {

                try {
                    using (var fStr = new StreamWriter(new FileStream(f, FileMode.CreateNew, FileAccess.Write, FileShare.None))) {

                        fStr.WriteLine(Name);
                        fStr.WriteLine(ExeDir != null ? ExeDir : "");
                        fStr.WriteLine(NoOfProcs);
                        fStr.WriteLine(exefile);

                        fStr.WriteLine(this.Arguments.Length);
                        for (int i = 0; i < this.Arguments.Length; i++) {
                            fStr.WriteLine(Arguments[i]);
                        }

                        fStr.WriteLine(this.EnvVars.Length);
                        for (int i = 0; i < this.EnvVars.Length; i++) {
                            fStr.WriteLine(this.EnvVars[i].Item1);
                            fStr.WriteLine(this.EnvVars[i].Item2);
                        }

                        return;
                    }
                } catch (Exception e) {
                    try {
                        if (File.Exists(f))
                            // try to avoid invalid files
                            File.Delete(f);
                    } catch (Exception) {

                    }
                    
                    if (ReTryCount < ClientAndServer.IO_OPS_MAX_RETRY_COUNT) {
                        ReTryCount++;
                        Thread.Sleep(ClientAndServer.IOwaitTime);
                    } else {
                        throw e;
                    }
                }
            }
        }
    }
}

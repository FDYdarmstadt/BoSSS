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
    /// Client interface: submitting a job, etc.
    /// </summary>
    public class Client : ClientAndServer {

        /// <summary>
        /// %
        /// </summary>
        /// <param name="BatchInstructionDir">
        /// if null, the default subdirectory within the user home (~/.BoSSS/batch)
        /// </param>
        public Client(string BatchInstructionDir = null) : base(BatchInstructionDir) {
            
        }

        /// <summary>
        /// Puts a job into the waiting queue.
        /// </summary>
        public int SubmitJob(JobData JD) {
            JD.m_ID = base.GetNewId();
            JD.Save(config.BatchInstructionDir);
            return JD.ID;
        }

        /// <summary>
        /// Path to a jobs standard output file.
        /// </summary>
        public string GetStdoutFile(int JobId) {
            string f = Path.Combine(config.BatchInstructionDir, ClientAndServer.WORK_FINISHED_DIR, JobId.ToString() + "_out.txt");
            return f;
        }

        /// <summary>
        /// Path to a jobs standard error file.
        /// </summary>
        public string GetStderrFile(int JobId) {
            string f = Path.Combine(config.BatchInstructionDir, ClientAndServer.WORK_FINISHED_DIR, JobId.ToString() + "_err.txt");
            return f;
        }
    }
}
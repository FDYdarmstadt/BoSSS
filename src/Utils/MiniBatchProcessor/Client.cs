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
using System.Threading;
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
            //JD.m_ID = base.GetNewId();
            //JD.Save(config.BatchInstructionDir);
            //return JD.ID;

            Exception lastExc = null; ;
            int Retry = 0;
            while (Retry < ClientAndServer.IO_OPS_MAX_RETRY_COUNT) {
                Retry++;
                try {
                    int newID;
                    if(base.AllJobs != null && base.AllJobs.Count > 0)
                        newID = base.AllJobs.Select(j => j.ID).Max() + 1;
                    else
                        newID = 1;

                    JD.m_ID = newID;
                    string f = Path.Combine(config.BatchInstructionDir, ClientAndServer.QUEUE_DIR, newID.ToString());

                    if(File.Exists(f)) {
                        // to slow -- race condition
                        Retry--;
                        Thread.Sleep(ClientAndServer.IOwaitTime);
                        continue;
                    }

                    using(var fStr = new StreamWriter(new FileStream(f, FileMode.CreateNew, FileAccess.Write, FileShare.None))) {

                        JD.Write(fStr);

                    }

                    return newID;
                } catch (Exception e) {
                    lastExc = e;
                }
            }

            throw lastExc;
        }

    }
}
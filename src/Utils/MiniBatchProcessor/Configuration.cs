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
    /// Configuration (directories, file paths, etc.) of the job manager.
    /// </summary>
    public class Configuration {

        /// <summary>
        /// ctor.
        /// </summary>
        public Configuration() {
            string BoSSSuserDir = BoSSS.Foundation.IO.Utils.GetBoSSSUserSettingsPath();
            BatchInstructionDir = Path.Combine(BoSSSuserDir, "batch");

            foreach (var dir in new string[] {
                Path.Combine(BatchInstructionDir, ClientAndServer.FINISHED_DIR),
                Path.Combine(BatchInstructionDir, ClientAndServer.QUEUE_DIR),
                Path.Combine(BatchInstructionDir, ClientAndServer.WORK_DIR)
            }) {
                if (!Directory.Exists(dir))
                    Directory.CreateDirectory(dir);
            }
        }

        /// <summary>
        /// Directory for job queues (see <see cref="ClientAndServer.WORK_DIR"/>, <see cref="ClientAndServer.QUEUE_DIR"/>,
        /// <see cref="ClientAndServer.FINISHED_DIR"/>).
        /// </summary>
        public string BatchInstructionDir;

        /// <summary>
        /// Maximum number of processors on this machine
        /// </summary>
        public int MaxProcessors {
            get {
                return 1;
                //return Math.Max(1, Environment.ProcessorCount / 2);  // most of our machines have hyper-threading
            }
        }

        /// <summary>
        /// If true, small jobs are allowed to overtake - this may delay or even prevent the start of large jobs.
        /// </summary>
        public bool BackFilling = false;
    }
}

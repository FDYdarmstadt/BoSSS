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
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MiniBatchProcessor {

    /// <summary>
    /// Exactly what it sounds to be.
    /// </summary>
    public enum JobStatus {

        /// <summary>
        /// Job status is not defined yet, e.g. if job-file has not been written to disk.
        /// </summary>
        Undefined = 0,

        /// <summary>
        /// Job is waiting to be executed, file is in directory <see cref="ClientAndServer.QUEUE_DIR"/>.
        /// </summary>
        Queued = 1,

        /// <summary>
        /// Job is currently being executed, file is in directory <see cref="ClientAndServer.WORK_DIR"/>.
        /// </summary>
        Working = 2,

        /// <summary>
        /// Job is finished, file is in directory <see cref="ClientAndServer.FINISHED_DIR"/>.
        /// </summary>
        Finished = 3
    }
}

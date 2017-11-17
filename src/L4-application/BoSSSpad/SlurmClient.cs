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

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// A <see cref="BatchProcessorClient"/> implementation for slurm systems on unix based hpc platforms
    /// </summary>
    class SlurmClient : BatchProcessorClient {
        public override void EvaluateStatus(Job myJob, out int SubmitCount, out bool isRunning, out bool wasSuccessful, out bool isFailed, out string DeployDir) {
            throw new NotImplementedException();
        }

        public override string GetStderrFile(Job myJob) {
            throw new NotImplementedException();
        }

        public override string GetStdoutFile(Job myJob) {
            throw new NotImplementedException();
        }

        public override object Submit(Job myJob) {
            throw new NotImplementedException();
        }
    }
}

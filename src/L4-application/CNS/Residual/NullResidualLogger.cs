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
using BoSSS.Foundation;
using System.Collections.Generic;

namespace CNS.Residual {

    /// <summary>
    /// No residuals are calculated. This class just serves as a realization of
    /// the GoF null object pattern.
    /// </summary>
    public class NullResidualLogger : ResidualLogger {

        /// <summary>
        /// <see cref="ResidualLogger"/>
        /// </summary>
        public NullResidualLogger(BoSSS.Solution.ResidualLogger baseLogger, SessionInfo currentSession, CNSFieldSet workingSet)
            : base(baseLogger, currentSession, workingSet, 0) {
        }

        /// <summary>
        /// Empty
        /// </summary>
        /// <param name="timeStepNumber">Not used</param>
        /// <param name="dt">Not used</param>
        /// <param name="phystime">Not used</param>
        /// <returns>Null</returns>
        public override IDictionary<string, double> LogTimeStep(int timeStepNumber, double dt, double phystime) {
            return new Dictionary<string, double>();
        }
    }
}

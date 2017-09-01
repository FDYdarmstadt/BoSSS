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
using BoSSS.Foundation;
using BoSSS.Foundation.IO;

namespace CNS.Residual {

    /// <summary>
    /// Base class for residual loggers
    /// </summary>
    public abstract class ResidualLogger : IResidualLogger {

        /// <summary>
        /// Access to the actual logger
        /// </summary>
        protected readonly BoSSS.Solution.ResidualLogger baseLogger;

        /// <summary>
        /// the current session; determines where log files are stored.
        /// </summary>
        private SessionInfo m_currentsession;

        /// <summary>
        /// the current session - this determines where log files will be stored.
        /// </summary>
        public SessionInfo CurrentSession {
            get {
                return m_currentsession;
            }
        }

        /// <summary>
        /// Reference to the working set that is updated every time-step and
        /// thus always reflects the current flow state.
        /// </summary>
        protected VectorField<DGField> CurrentState;

        /// <summary>
        /// Deep copy of the working created after the last call of
        /// <see cref="LogTimeStep"/> (Before the first call of this method, it
        /// is initialized with the initial value)
        /// </summary>
        protected VectorField<DGField> PreviousState;

        /// <summary>
        /// Interval at which residuals are calculated
        /// </summary>
        protected readonly int residualInterval;

        /// <summary>
        /// Initializes <see cref="PreviousState"/>, <see cref="CurrentState"/>
        /// and the log files.
        /// </summary>
        /// <param name="baseLogger">Access to the actual log file</param>
        /// <param name="currentsession">The current session</param>
        /// <param name="workingSet">
        /// The initial working set ("0-th" iteration)
        /// </param>
        /// <param name="residualInterval">
        /// Interval at which residuals are calculated
        /// </param>
        public ResidualLogger(BoSSS.Solution.ResidualLogger baseLogger, SessionInfo currentsession, CNSFieldSet workingSet, int residualInterval) {
            this.baseLogger = baseLogger;
            this.residualInterval = residualInterval;
            this.m_currentsession = currentsession;

            CurrentState = new VectorField<DGField>(workingSet.ConservativeVariables);
            PreviousState = CurrentState.CloneAs();
        }

        #region IResidualLogger Members

        /// <summary>
        /// Calculates and writes the residual information for one time-step to
        /// the log files
        /// </summary>
        /// <param name="timeStepNumber">
        /// The iteration number associated with <see cref="CurrentState"/>
        /// </param>
        /// <param name="dt">
        /// The time increment from <see cref="PreviousState"/> to
        /// <see cref="CurrentState"/>
        /// </param>
        /// <param name="phystime"></param>
        public abstract IDictionary<string, double> LogTimeStep(int timeStepNumber, double dt, double phystime);

        #endregion
    }
}

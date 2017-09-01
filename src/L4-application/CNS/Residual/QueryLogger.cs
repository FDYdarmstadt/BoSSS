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
using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Solution;
using BoSSS.Solution.Queries;

namespace CNS.Residual {

    /// <summary>
    /// A residual logger that logs results based on the output of the queries
    /// defined in the control file.
    /// </summary>
    public class QueryLogger : IResidualLogger {

        private BoSSS.Solution.ResidualLogger baseLogger;

        /// <summary>
        /// the current session - this determines where log files will be stored.
        /// </summary>
        public SessionInfo CurrentSession {
            get;
            private set;
        }

        /// <summary>
        /// Queries to be executed in every logged time-step
        /// </summary>
        private IDictionary<string, Query> queries;

        /// <summary>
        /// Interval at which residuals are calculated
        /// </summary>
        private int interval;

        /// <summary>
        /// Results from the last query evaluation
        /// </summary>
        private Dictionary<string, double> lastResult;

        private IApplication<CNSControl> application;

        /// <summary>
        /// Creates a new logger for the configured queries.
        /// </summary>
        /// <param name="baseLogger">
        /// The real logger
        /// </param>
        /// <param name="application"></param>
        public QueryLogger(BoSSS.Solution.ResidualLogger baseLogger, IApplication<CNSControl> application) {
            this.baseLogger = baseLogger;
            this.CurrentSession = application.CurrentSessionInfo;
            this.queries = application.QueryHandler.QueryMap;
            this.interval = application.Control.ResidualInterval;
            this.application = application;
        }

        #region IResidualLogger Members

        /// <summary>
        /// Calculates and writes the residual information based on the query
        /// results for one time-step to the log files
        /// </summary>
        /// <param name="timeStepNumber">The iteration number</param>
        /// <param name="dt">The time increment</param>
        /// <param name="phystime">Physical time of previous time-step</param>
        /// <remarks>
        /// The residual information based on the query results
        /// </remarks>
        public IDictionary<string, double> LogTimeStep(int timeStepNumber, double dt, double phystime) {
            Dictionary<string, double> fullResult = new Dictionary<string, double>();
            if (interval < 1) {
                return fullResult;
            }

            if (this.ShouldLog(timeStepNumber, interval)) {
                Dictionary<string, double> valueResult = new Dictionary<string, double>();
                Dictionary<string, double> valueChangeResult = new Dictionary<string, double>();

                foreach (var idQueryPair in queries) {
                    // Note the + dt since time-stepper has already performed at
                    // this point, but physical time has not been updated yet.
                    string id = "query_" + idQueryPair.Key;
                    double queryResult = idQueryPair.Value(application, phystime + dt);

                    valueResult.Add(id, queryResult);
                    fullResult.Add(id, queryResult);

                    double oldQueryResult = 0.0;
                    if (lastResult != null) {
                        oldQueryResult = lastResult[id];
                    }
                    double changeResult = Math.Abs(queryResult - oldQueryResult);

                    valueChangeResult.Add(id + "_change", changeResult);
                    fullResult.Add(id + "_change", changeResult);

                    baseLogger.CustomValue(queryResult, id);
                    baseLogger.CustomValue(changeResult, id + "_change");
                }

                lastResult = valueResult;
            }

            return fullResult;
        }

        #endregion
    }
}

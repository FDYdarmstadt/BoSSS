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

using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.IO;

namespace CNS.Residual {

    /// <summary>
    /// Common interface for residual loggers.
    /// </summary>
    public interface IResidualLogger {

        /// <summary>
        /// the current session - this determines where log files will be stored.
        /// </summary>
        SessionInfo CurrentSession {
            get;
        }

        /// <summary>
        /// Calculates and writes the residual information for one time-step to
        /// the log files
        /// </summary>
        /// <param name="timeStepNumber">
        /// The iteration number.
        /// </param>
        /// <param name="dt">
        /// The time increment.
        /// </param>
        /// <param name="phystime">
        /// Physical time
        /// </param>
        /// <returns>
        /// The residual information indexed by a unique residual name
        /// </returns>
        IDictionary<string, double> LogTimeStep(int timeStepNumber, double dt, double phystime);
    }

    /// <summary>
    /// Extension methods for <see cref="IResidualLogger"/>
    /// </summary>
    public static class IResidualLoggerExtensions {

        /// <summary>
        /// <see cref="ShouldLog(IResidualLogger, int, int)"/>
        /// </summary>
        /// <returns></returns>
        public static bool ShouldLog(this IEnumerable<IResidualLogger> loggers, int timeStepNumber, int residualInterval) {
            return (loggers != null && loggers.FirstOrDefault().ShouldLog(timeStepNumber, residualInterval));
        }

        /// <summary>
        /// Determines whether the residuals for the given
        /// <paramref name="timeStepNumber"/> should be logged according to the
        /// given <paramref name="residualInterval"/>
        /// </summary>
        /// <param name="logger">
        /// The logger to be tested
        /// </param>
        /// <param name="timeStepNumber">
        /// The number of the current time-step
        /// </param>
        /// <param name="residualInterval">
        /// The interval between residuals calculations
        /// </param>
        /// <returns></returns>
        public static bool ShouldLog(this IResidualLogger logger, int timeStepNumber, int residualInterval) {
            return (logger != null
                && residualInterval > 0
                && (timeStepNumber == 1 || timeStepNumber % residualInterval == 0));
        }

        /// <summary>
        /// Calls <see cref="IResidualLogger.LogTimeStep"/> for each defined
        /// residual logger.
        /// </summary>
        /// <param name="loggers">
        /// The loggers which should log a time-step
        /// </param>
        /// <param name="timeStepNumber">
        /// <see cref="IResidualLogger.LogTimeStep"/>
        /// </param>
        /// <param name="dt">
        /// <see cref="IResidualLogger.LogTimeStep"/>
        /// </param>
        /// <param name="phystime">
        /// <see cref="IResidualLogger.LogTimeStep"/>
        /// </param>
        public static IDictionary<string, double> LogTimeStep(this IEnumerable<IResidualLogger> loggers, int timeStepNumber, double dt, double phystime) {
            Dictionary<string, double> result = new Dictionary<string, double>();
            foreach (var logger in loggers) {
                IDictionary<string, double> localResult = logger.LogTimeStep(timeStepNumber, dt, phystime);
                foreach (var pair in localResult) {
                    result.Add(pair.Key, pair.Value);
                }
            }
            return result;
        }
    }
}

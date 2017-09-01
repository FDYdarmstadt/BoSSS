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
using BoSSS.Foundation.IO;

namespace CNS.Residual {

    /// <summary>
    /// Utility class for the calculation of "pseudo-residuals". That is, the
    /// change of the solution between the last states is calculated which may
    /// serve as a residual indicator in the case of steady-state problems.
    /// Compared to <see cref="RigorousResidualLogger{T}"/>, the evaluation should
    /// be much cheaper (both, in terms of memory and computing time).
    /// </summary>
    public class ChangeRateResidualLogger : ResidualLogger {

        /// <summary>
        /// <see cref="ResidualLogger"/>
        /// </summary>
        public ChangeRateResidualLogger(BoSSS.Solution.ResidualLogger baseLogger, SessionInfo currentSession, CNSFieldSet workingSet, int residualInterval)
            : base(baseLogger, currentSession, workingSet, residualInterval) {
        }

        /// <summary>
        /// Calculates the residuals solely based on the difference between
        /// <see cref="ResidualLogger.PreviousState"/> and
        /// <see cref="ResidualLogger.CurrentState"/> (which is adequate for
        /// steady-state problems)
        /// </summary>
        /// <param name="timeStepNumber"><see cref="ResidualLogger"/></param>
        /// <param name="dt">Not used</param>
        /// <param name="phystime">Not used</param>
        /// <returns>
        /// The difference between <see cref="ResidualLogger.PreviousState"/>
        /// and <see cref="ResidualLogger.CurrentState"/> int the L2 norm.
        /// </returns>
        public override IDictionary<string, double> LogTimeStep(int timeStepNumber, double dt, double phystime) {
            Dictionary<string, double> result = new Dictionary<string, double>();
            if (residualInterval < 1) {
                return result;
            }

            if (this.ShouldLog(timeStepNumber, residualInterval)) {
                for (int i = 0; i < CurrentState.Dim; i++) {
                    PreviousState[i].Acc(-1.0, CurrentState[i]);
                    double absoluteResidualsL2 = PreviousState[i].L2Norm();
                    double relativeResidualsL2 = absoluteResidualsL2 / CurrentState[i].L2Norm();

                    result.Add("changeRate_abs_" + CurrentState[i].Identification, absoluteResidualsL2);

                    baseLogger.CustomValue(absoluteResidualsL2, "changeRate_abs_" + CurrentState[i].Identification);
                    baseLogger.CustomValue(relativeResidualsL2, "changeRate_rel_" + CurrentState[i].Identification);
                }
            }

            if ((timeStepNumber + 1) % residualInterval == 0) {
                // Residuals will be calculated in next time-step
                PreviousState = CurrentState.CloneAs();
            }

            return result;
        }
    }
}

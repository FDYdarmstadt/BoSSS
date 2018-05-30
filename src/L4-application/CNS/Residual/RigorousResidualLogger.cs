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
using BoSSS.Foundation;
using BoSSS.Foundation.IO;

namespace CNS.Residual {

    /// <summary>
    /// Utility class for the calculation of residuals. In contrast to
    /// <see cref="ChangeRateResidualLogger"/>, it evaluates the
    /// <see cref="SpatialOperator"/> defining the equation system
    /// to be solved and approximates the time-derivative by a backward-
    /// difference. This makes it applicable to unsteady problems. 
    /// </summary>
    public class RigorousResidualLogger<T> : ResidualLogger
            where T : CNSControl, new() {

        /// <summary>
        /// The spatial differential operator defining (the spatial part of)
        /// the system of equations that are solved.
        /// </summary>
        private IEvaluatorNonLin evaluator;

        /// <summary>
        /// <see cref="ResidualLogger"/>
        /// </summary>
        /// <param name="program"></param>
        /// <param name="residualInterval">
        /// <see cref="ResidualLogger"/>
        /// </param>
        /// <param name="differentialOperator">
        /// The spatial differential operator defining (the spatial part of)
        /// the system of equations that are solved.
        /// </param>
        public RigorousResidualLogger(Program<T> program, int residualInterval, SpatialOperator differentialOperator)
            : base(program.ResLogger, program.CurrentSessionInfo, program.WorkingSet, residualInterval) {
            DGField[] primalFields = program.WorkingSet.ConservativeVariables;
            CoordinateMapping domainMapping = new CoordinateMapping(primalFields);
            UnsetteledCoordinateMapping codomainMapping = new UnsetteledCoordinateMapping(
                primalFields.Select((field) => field.Basis).ToArray());
            evaluator = differentialOperator.GetEvaluatorEx(
                domainMapping,
                program.ParameterMapping,
                codomainMapping);
        }

        /// <summary>
        /// For a system
        /// \f$ 
        /// \frac{\partial U}{dt} + F(U) = 0,
        /// \f$ 
        /// evaluates
        /// \f$ 
        /// || \frac{U_1 - U_0}{\Delta t} + F(U_1) ||
        /// \f$ 
        /// in the associated relative residuals in L2 and the Linf norm.
        /// </summary>
        /// <param name="timeStepNumber"></param>
        /// <param name="dt"></param>
        /// <param name="phystime">not used</param>
        public override IDictionary<string, double> LogTimeStep(int timeStepNumber, double dt, double phystime) {
            Dictionary<string, double> result = new Dictionary<string, double>();
            if (residualInterval < 1) {
                return result;
            }

            if (this.ShouldLog(timeStepNumber, residualInterval)) {
                double beta;
                if (dt == 0.0) {
                    // r = F(u1)
                    beta = 0.0;
                    PreviousState.Clear();
                } else {
                    // r = (u1-u0)/dt + F(u1)
                    beta = -1.0 / dt;
                    for (int i = 0; i < CurrentState.Dim; i++) {
                        PreviousState[i].Acc(-1.0, CurrentState[i]);
                    }
                }

                CoordinateVector resultVector = new CoordinateVector(
                    new CoordinateMapping(PreviousState.ToArray()));
                evaluator.Evaluate(1.0, beta, resultVector);

                for (int i = 0; i < CurrentState.Dim; i++) {
                    double abs = PreviousState[i].L2Norm();
                    double rel = abs / CurrentState[i].L2Norm();

                    result.Add("rigorous_L2_abs_" + CurrentState[i].Identification, abs);

                    baseLogger.CustomValue(abs, "residual_abs_" + CurrentState[i].Identification);
                    baseLogger.CustomValue(rel, "residual_rel_" + CurrentState[i].Identification);
                }
            }

            // If residuals will be calculated in next time-step
            if ((timeStepNumber + 1) % residualInterval == 0
                && dt != 0.0) {
                PreviousState = CurrentState.CloneAs();
            }

            return result;
        }
    }
}

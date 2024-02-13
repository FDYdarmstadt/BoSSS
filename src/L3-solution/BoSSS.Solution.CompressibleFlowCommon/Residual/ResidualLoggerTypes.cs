﻿/* =======================================================================
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

using BoSSS.Foundation;
using System;
using System.Collections.Generic;

namespace BoSSS.Solution.CompressibleFlowCommon.Residual {

    /// <summary>
    /// Available types of residual loggers.
    /// </summary>
    [Flags]
    public enum ResidualLoggerTypes {

        /// <summary>
        /// No residuals are logged, see <see cref="NullResidualLogger"/>.
        /// </summary>
        None = 0,

        /// <summary>
        /// Residuals are logged using <see cref="ChangeRateResidualLogger"/>
        /// which implies that the problem to be solved is a steady-state
        /// problem.
        /// </summary>
        ChangeRate = 1,

        /// <summary>
        /// Residuals are logged using <see cref="RigorousResidualLogger{T}"/>.
        /// </summary>
        Rigorous = 2,

        /// <summary>
        /// The results of the selected queries are taken as residual/error
        /// indicators.
        /// </summary>
        Query = 4
    }

    /// <summary>
    /// Extension methods for <see cref="ResidualLoggerTypes"/>.
    /// </summary>
    public static class ResidualLoggerExtensions {

        /// <summary>
        /// Factory for residual loggers. The instantiated objects depend on the
        /// defined <paramref name="loggerType"/>.
        /// </summary>
        /// <param name="loggerType">
        /// The type of logger to be instantiated
        /// </param>
        /// <param name="program">
        /// The program requesting the residual logger.
        /// </param>
        /// <param name="config">Configuration options</param>
        /// <param name="differentialOperator">
        /// The differential operator that defines the system of equations to
        /// be solved. May be null if <paramref name="loggerType"/> does
        /// <b>not</b> contain <see cref="ResidualLoggerTypes.Rigorous"/>.
        /// </param>
        /// <returns>
        /// A list of residual loggers, see <see cref="ResidualLogger"/>.
        /// </returns>
        public static IEnumerable<IResidualLogger> Instantiate<T>(
            this ResidualLoggerTypes loggerType,
            Application<T> program,
            CompressibleControl config,
            DifferentialOperator differentialOperator,
            DGField[] consVars,
            CoordinateMapping paramMap)
            where T : CompressibleControl, new() {

            if (loggerType.HasFlag(ResidualLoggerTypes.ChangeRate)
                && loggerType.HasFlag(ResidualLoggerTypes.Rigorous)) {
                throw new Exception(
                    "Residual types \"changeRate\" and \"rigorous\" are mutually exclusive");
            }

            if (loggerType == ResidualLoggerTypes.None) {
                yield return new NullResidualLogger(
                    program.ResLogger, program.CurrentSessionInfo, consVars);
            } else {
                if (loggerType.HasFlag(ResidualLoggerTypes.ChangeRate)) {
                    loggerType ^= ResidualLoggerTypes.ChangeRate;
                    yield return new ChangeRateResidualLogger(
                        program.ResLogger, program.CurrentSessionInfo, consVars, config.ResidualInterval);
                }

                if (loggerType.HasFlag(ResidualLoggerTypes.Rigorous)) {
                    loggerType ^= ResidualLoggerTypes.Rigorous;
                    yield return new RigorousResidualLogger<T>(
                        program,
                        consVars,
                        paramMap,
                        config.ResidualInterval,
                        differentialOperator);
                }

                if (loggerType.HasFlag(ResidualLoggerTypes.Query)) {
                    loggerType ^= ResidualLoggerTypes.Query;
                    yield return new QueryLogger(
                        program.ResLogger, program);
                }

                if (loggerType != ResidualLoggerTypes.None) {
                    throw new NotImplementedException(
                        "Residual logging for residual type " + loggerType + " not implemented");
                }
            }
        }
    }
}

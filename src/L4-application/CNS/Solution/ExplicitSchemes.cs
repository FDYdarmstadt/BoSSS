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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using BoSSS.Solution;
using BoSSS.Solution.Timestepping;
using CNS.EquationSystem;
using CNS.Exception;
using CNS.IBM;
using CNS.ShockCapturing;
using System;
using System.Collections;
using System.Linq;

namespace CNS.Solution {

    /// <summary>
    /// Supported classes of time-stepping schemes
    /// </summary>
    public enum ExplicitSchemes {

        /// <summary>
        /// Explicit time stepping deactivated.
        /// </summary>
        None = 0,

        /// <summary>
        /// Family of Runge-Kutta schemes, <see cref="RungeKutta"/>
        /// </summary>
        RungeKutta,

        /// <summary>
        /// Family of Adams-Bashforth multi-step schemes,
        /// <see cref="AdamsBashforth"/>
        /// </summary>
        AdamsBashforth,

        /// <summary>
        /// Family of Adams-Bashforth multi-step schemes with Local Time Stepping
        /// <see cref="AdamsBashforth"/>
        /// </summary>
        LTS,

        /// <summary>
        /// Family of adaptive stabilized Chebyshev Runge-Kutta schemes
        /// <see cref="Rock4"/>
        /// </summary>
        Rock4,

        /// <summary>
        /// Special fourth order Runge-Kutta scheme with 5 stages; see
        /// <see cref="RungeKutta.RungeKuttaScheme.SSP54"/>.
        /// </summary>
        SSP54,

        /// <summary>
        /// Special fourth order Runge-Kutta scheme with 8 stages; see
        /// <see cref="RungeKutta.RungeKuttaScheme.RKC84"/>.
        /// </summary>
        RKC84
    }

    /// <summary>
    /// Extension methods for <see cref="ExplicitSchemes"/>
    /// </summary>
    public static class ExplicitSchemesExtensions {

        /// <summary>
        /// Creates the appropriate time-stepper for a given
        /// <paramref name="timeStepperType"/>
        /// </summary>
        /// <param name="timeStepperType">
        /// The type of the time-stepper (e.g., Runge-Kutta) to be created
        /// </param>
        /// <param name="control">
        /// Configuration options
        /// </param>
        /// <param name="equationSystem">
        /// The equation system to be solved
        /// </param>
        /// <param name="fieldsMap">
        /// Fields affected by the time-stepper
        /// </param>
        /// <param name="parameterMap">
        /// Fields serving as parameter for the time-stepper.
        /// </param>
        /// <param name="speciesMap">
        /// Mapping of different species inside the domain
        /// </param>
        /// <param name="program"></param>
        /// <returns>
        /// Depending on <paramref name="timeStepperType"/>, an appropriate
        /// implementation of <see cref="ITimeStepper"/>.
        /// </returns>
        /// <remarks>
        /// Currently limiting is not supported for Adams-Bashforth methods
        /// </remarks>
        public static ITimeStepper Instantiate<T>(
            this ExplicitSchemes timeStepperType,
            CNSControl control,
            OperatorFactory equationSystem,
            CoordinateMapping fieldsMap,
            CoordinateMapping parameterMap,
            ISpeciesMap speciesMap,
            IProgram<T> program)
            where T : CNSControl, new() {

            ITimeStepper timeStepper;
            switch (timeStepperType) {
                case ExplicitSchemes.RungeKutta:
                    //RungeKutta timeStepper;
                    if (control.DomainType == DomainTypes.Standard) {
                        // HACK
                        //RungeKutta timeStepperBla = new RungeKutta(
                        timeStepper = new RungeKutta(
                        RungeKutta.GetDefaultScheme(control.ExplicitOrder),
                        equationSystem.GetJoinedOperator().ToSpatialOperator(),
                        fieldsMap,
                        parameterMap,
                        equationSystem.GetJoinedOperator().CFLConstraints);

                        //timeStepperBla.OnBeforeComputeChangeRate += (t1, t2) => program.WorkingSet.UpdateDerivedVariables(program, program.SpeciesMap.SubGrid.VolumeMask);
                        //timeStepper = timeStepperBla;
                    } else {
                        IBMControl ibmControl = (IBMControl)control;
                        timeStepper = ibmControl.TimesteppingStrategy.CreateRungeKuttaTimeStepper(
                            ibmControl, equationSystem, fieldsMap, parameterMap, speciesMap,
                            equationSystem.GetJoinedOperator().CFLConstraints);
                    }
                    break;

                case ExplicitSchemes.AdamsBashforth:
                    if (control.DomainType == DomainTypes.StaticImmersedBoundary
                        || control.DomainType == DomainTypes.MovingImmersedBoundary) {
                        timeStepper = new IBMAdamsBashforth(
                            equationSystem.GetConvectiveOperator().Union(equationSystem.GetDiffusiveOperator()).ToSpatialOperator(),
                            equationSystem.GetSourceTermOperator().ToSpatialOperator(),
                            fieldsMap,
                            parameterMap,
                            speciesMap,
                            (IBMControl)control,
                            equationSystem.GetJoinedOperator().CFLConstraints);
                    } else {
                        // HACK
                        //AdamsBashforth timeStepperBla = new AdamsBashforth(
                        timeStepper = new AdamsBashforth(
                        equationSystem.GetJoinedOperator().ToSpatialOperator(),
                        fieldsMap,
                        parameterMap,
                        control.ExplicitOrder,
                        equationSystem.GetJoinedOperator().CFLConstraints);

                        //timeStepperBla.OnBeforeComputeChangeRate += (t1, t2) => program.WorkingSet.UpdateDerivedVariables(program, program.SpeciesMap.SubGrid.VolumeMask);
                        //timeStepper = timeStepperBla;
                    }
                    break;

                case ExplicitSchemes.LTS:
                    if (control.DomainType == DomainTypes.StaticImmersedBoundary
                        || control.DomainType == DomainTypes.MovingImmersedBoundary) {
                        timeStepper = new IBMAdamsBashforthLTS(
                            equationSystem.GetConvectiveOperator().Union(equationSystem.GetDiffusiveOperator()).ToSpatialOperator(),
                            equationSystem.GetSourceTermOperator().ToSpatialOperator(),
                            fieldsMap,
                            parameterMap,
                            speciesMap,
                            (IBMControl)control,
                            equationSystem.GetJoinedOperator().CFLConstraints);
                    } else {
                        timeStepper = new AdamsBashforthLTS(
                            equationSystem.GetJoinedOperator().ToSpatialOperator(),
                            fieldsMap,
                            parameterMap,
                            control.ExplicitOrder,
                        control.NumberOfSubGrids,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    }
                    break;

                case ExplicitSchemes.Rock4:
                    timeStepper = new ROCK4(equationSystem.GetJoinedOperator().ToSpatialOperator(), new CoordinateVector(fieldsMap), null);
                    break;

                case ExplicitSchemes.SSP54:
                    timeStepper = new RungeKutta(
                        RungeKuttaScheme.SSP54,
                        equationSystem.GetJoinedOperator().ToSpatialOperator(),
                        fieldsMap,
                        parameterMap,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    break;

                case ExplicitSchemes.RKC84:
                    timeStepper = new RungeKutta(
                        RungeKuttaScheme.RKC84,
                        equationSystem.GetJoinedOperator().ToSpatialOperator(),
                        fieldsMap,
                        parameterMap,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    break;

                case ExplicitSchemes.None:
                    throw new System.ArgumentException("Cannot instantiate empty scheme");

                default:
                    throw new Exception.InternalErrorException(String.Format(
                        "Unknown explicit time stepper type \"{0}\"", timeStepperType));
            }

            if (control.Limiter != null) {
                ExplicitEuler explicitEulerBasedTimestepper = timeStepper as ExplicitEuler;
                if (explicitEulerBasedTimestepper == null) {
                    throw new ConfigurationException(String.Format(
                        "Limiting currently not implemented for time-steppers of type '{0}~",
                        timeStepperType));
                } else {
                    explicitEulerBasedTimestepper.OnAfterFieldUpdate +=
                        (t, f) => control.Limiter.LimitFieldValues(program);
                }
            }

            return timeStepper;
        }
    }
}

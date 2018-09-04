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
using BoSSS.Solution;
using BoSSS.Solution.Timestepping;
using CNS.EquationSystem;
using CNS.IBM;
using System;

namespace CNS {

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
        /// <see cref="AdamsBashforthLTS"/>
        /// </summary>
        LTS,

        /// <summary>
        /// Family of adaptive stabilized Chebyshev Runge-Kutta schemes
        /// <see cref="Rock4"/>
        /// </summary>
        Rock4,

        /// <summary>
        /// Special fourth order Runge-Kutta scheme with 5 stages; see
        /// <see cref="RungeKuttaScheme.SSP54"/>.
        /// </summary>
        SSP54,

        /// <summary>
        /// Special fourth order Runge-Kutta scheme with 8 stages; see
        /// <see cref="RungeKuttaScheme.RKC84"/>.
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
        /// <param name="fieldSet">
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
            CNSFieldSet fieldSet,
            CoordinateMapping parameterMap,
            ISpeciesMap speciesMap,
            IProgram<T> program)
            where T : CNSControl, new() {

            CoordinateMapping variableMap = new CoordinateMapping(fieldSet.ConservativeVariables);

            IBMControl ibmControl = control as IBMControl;
            IBMOperatorFactory ibmFactory = equationSystem as IBMOperatorFactory;
            if (control.DomainType != DomainTypes.Standard
                && (ibmFactory == null || ibmControl == null)) {
                throw new Exception();
            }

            ITimeStepper timeStepper;
            switch (timeStepperType) {
                case ExplicitSchemes.RungeKutta when control.DomainType == DomainTypes.Standard:
                    timeStepper = new RungeKutta(
                        RungeKutta.GetDefaultScheme(control.ExplicitOrder),
                        equationSystem.GetJoinedOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    break;

                case ExplicitSchemes.RungeKutta:
                    timeStepper = ibmControl.TimesteppingStrategy.CreateRungeKuttaTimeStepper(
                        ibmControl,
                        equationSystem,
                        fieldSet,
                        parameterMap,
                        speciesMap,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    break;

                case ExplicitSchemes.AdamsBashforth when control.DomainType == DomainTypes.Standard:
                    timeStepper = new AdamsBashforth(
                        equationSystem.GetJoinedOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        control.ExplicitOrder,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    break;

                case ExplicitSchemes.AdamsBashforth:
                    timeStepper = new IBMAdamsBashforth(
                        ibmFactory.GetJoinedOperator().ToSpatialOperator(fieldSet),
                        ibmFactory.GetImmersedBoundaryOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        speciesMap,
                        ibmControl,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    break;

                case ExplicitSchemes.LTS when control.DomainType == DomainTypes.Standard:
                    timeStepper = new AdamsBashforthLTS(
                        equationSystem.GetJoinedOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        control.ExplicitOrder,
                        control.NumberOfSubGrids,
                        equationSystem.GetJoinedOperator().CFLConstraints,
                        reclusteringInterval: control.ReclusteringInterval,
                        fluxCorrection: control.FluxCorrection,
                        saveToDBCallback: program.SaveToDatabase,
                        maxNumOfSubSteps: control.maxNumOfSubSteps,
                        forceReclustering: control.forceReclustering,
                        logging: control.WriteLTSLog,
                        consoleOutput: control.WriteLTSConsoleOutput);
                    break;

                case ExplicitSchemes.LTS:
                    timeStepper = new IBMAdamsBashforthLTS(
                        equationSystem.GetJoinedOperator().ToSpatialOperator(fieldSet),
                        ibmFactory.GetImmersedBoundaryOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        speciesMap,
                        ibmControl,
                        equationSystem.GetJoinedOperator().CFLConstraints,
                        reclusteringInterval: control.ReclusteringInterval,
                        fluxCorrection: control.FluxCorrection,
                        maxNumOfSubSteps: control.maxNumOfSubSteps,
                        forceReclustering: control.forceReclustering,
                        logging: control.WriteLTSLog,
                        consoleOutput: control.WriteLTSConsoleOutput);
                    break;

                case ExplicitSchemes.Rock4 when control.DomainType == DomainTypes.Standard:
                    timeStepper = new ROCK4(equationSystem.GetJoinedOperator().ToSpatialOperator(fieldSet), new CoordinateVector(variableMap), null);
                    break;

                case ExplicitSchemes.SSP54 when control.DomainType == DomainTypes.Standard:
                    timeStepper = new RungeKutta(
                        RungeKuttaScheme.SSP54,
                        equationSystem.GetJoinedOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    break;

                case ExplicitSchemes.RKC84 when control.DomainType == DomainTypes.Standard:
                    timeStepper = new RungeKutta(
                        RungeKuttaScheme.RKC84,
                        equationSystem.GetJoinedOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        equationSystem.GetJoinedOperator().CFLConstraints);
                    break;

                case ExplicitSchemes.None:
                    throw new System.ArgumentException("Cannot instantiate empty scheme");

                default:
                    throw new NotImplementedException(String.Format(
                        "Explicit time stepper type '{0}' not implemented for domain type '{1}'",
                        timeStepperType,
                        control.DomainType));
            }

            // Make sure shock sensor is updated before every flux evaluation
            if (control.ShockSensor != null) {
                ExplicitEuler explicitEulerBasedTimestepper = timeStepper as ExplicitEuler;
                if (explicitEulerBasedTimestepper == null) {
                    throw new Exception(String.Format(
                        "Shock-capturing currently not implemented for time-steppers of type '{0}'",
                        timeStepperType));
                }

                explicitEulerBasedTimestepper.OnBeforeComputeChangeRate += delegate (double absTime, double relTime) {
                    // Note: Only shock sensor is updated, _NOT_ the corresponding variable
                    program.Control.ShockSensor.UpdateSensorValues(
                        program.WorkingSet,
                        program.SpeciesMap,
                        explicitEulerBasedTimestepper.SubGrid.VolumeMask);
                    // Note: When being called, artificial viscosity is updated in the _ENTIRE_ (fluid) domain
                    var avField = program.WorkingSet.DerivedFields[Variables.ArtificialViscosity];
                    Variables.ArtificialViscosity.UpdateFunction(avField, program.SpeciesMap.SubGrid.VolumeMask, program);

                    // Test
                    //double sensorNorm = program.WorkingSet.DerivedFields[Variables.ShockSensor].L2Norm();
                    //double avNorm = program.WorkingSet.DerivedFields[Variables.ArtificialViscosity].L2Norm();
                    //Console.WriteLine("\r\nThis is OnBeforeComputeChangeRate");
                    //Console.WriteLine("SensorNeu: {0}", sensorNorm);
                    //Console.WriteLine("AVNeu: {0}", avNorm);
                };
            }

            // Make sure limiter is applied after each modification of conservative variables
            if (control.Limiter != null) {
                ExplicitEuler explicitEulerBasedTimestepper = timeStepper as ExplicitEuler;
                if (explicitEulerBasedTimestepper == null) {
                    throw new Exception(String.Format(
                        "Limiting currently not implemented for time-steppers of type '{0}~",
                        timeStepperType));
                }

                explicitEulerBasedTimestepper.OnAfterFieldUpdate +=
                    (t, f) => control.Limiter.LimitFieldValues(program);
            }

            return timeStepper;
        }
    }
}

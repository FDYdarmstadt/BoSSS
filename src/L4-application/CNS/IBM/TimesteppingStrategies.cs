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
using BoSSS.Solution.Timestepping;
using CNS.EquationSystem;
using System;
using BoSSS.Solution;
using System.Collections.Generic;

namespace CNS.IBM {

    public enum TimesteppingStrategies {

        LieSplitting,

        StrangSplitting,

        MovingFrameFlux
    }

    public static class TimesteppingStrategiesExtensions {

        public static IBMRungeKutta CreateRungeKuttaTimeStepper(
            this TimesteppingStrategies strategy,
            IBMControl control,
            OperatorFactory equationSystem,
            CNSFieldSet fieldSet,
            CoordinateMapping parameterMap,
            ISpeciesMap speciesMap,
            IList<TimeStepConstraint> timeStepConstraints) {

            ImmersedSpeciesMap ibmSpeciesMap = speciesMap as ImmersedSpeciesMap;
            if (ibmSpeciesMap == null) {
                throw new ArgumentException(
                    "Only supported for species maps of type 'ImmersedSpeciesMap'",
                    "speciesMap");
            }

            CoordinateMapping variableMap = new CoordinateMapping(fieldSet.ConservativeVariables);
            switch (strategy) {
                case TimesteppingStrategies.LieSplitting:
                case TimesteppingStrategies.StrangSplitting:
                    return new IBMSplitRungeKutta(
                        equationSystem.GetConvectiveOperator().Union(equationSystem.GetDiffusiveOperator()).ToSpatialOperator(fieldSet),
                        equationSystem.GetSourceTermOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        ibmSpeciesMap,
                        timeStepConstraints);

                case TimesteppingStrategies.MovingFrameFlux:
                    return new IBMMovingFrameRungeKutta(
                        equationSystem.GetConvectiveOperator().Union(equationSystem.GetDiffusiveOperator()).ToSpatialOperator(fieldSet),
                        equationSystem.GetSourceTermOperator().ToSpatialOperator(fieldSet),
                        variableMap,
                        parameterMap,
                        ibmSpeciesMap,
                        timeStepConstraints);

                default:
                    throw new System.NotImplementedException();
            }
        }
    }
}
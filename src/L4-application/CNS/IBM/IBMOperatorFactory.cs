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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using CNS.Boundary;
using CNS.EquationSystem;

namespace CNS.IBM {

    /// <summary>
    /// Specialization of <see cref="OperatorFactory"/> for immersed boundary
    /// flows.
    /// </summary>
    public class IBMOperatorFactory : OperatorFactory {

        /// <summary>
        /// Constructs a new operator factory which additionally implements
        /// the boundary conditions at immersed boundaries.
        /// </summary>
        /// <param name="control"></param>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        /// <param name="speciesMap"></param>
        /// <param name="boundaryMap"></param>
        public IBMOperatorFactory(
            IBMControl control,
            GridData gridData,
            CNSFieldSet workingSet,
            ISpeciesMap speciesMap,
            IBoundaryConditionMap boundaryMap)
            : base(control, gridData, workingSet, speciesMap, boundaryMap) {

            this.sourceTermBuilders.Add(new BoundaryConditionSourceFluxBuilder(
                control, boundaryMap, speciesMap, convectiveFluxBuilder, diffusiveFluxBuilder));
        }

        /// <summary>
        /// Overrides the behavior of
        /// <see cref="OperatorFactory.GetConvectiveOperator"/> by using
        /// <see cref="IBMConvectiveCFLConstraint"/> instead of
        /// <see cref="Convection.ConvectiveCFLConstraint"/>
        /// </summary>
        /// <returns></returns>
        public override Operator GetConvectiveOperator() {
            Operator op = new Operator(control);
            if (convectiveFluxBuilder != null) {
                convectiveFluxBuilder.BuildFluxes(op);
            }

            if (!op.IsEmpty) {
                op.CFLConstraints.Add(new IBMConvectiveCFLConstraint(
                    control,
                    gridData,
                    workingSet,
                    speciesMap));
            }

            return op;
        }

        /// <summary>
        /// Overrides the behavior of
        /// <see cref="OperatorFactory.GetDiffusiveOperator"/> by using
        /// <see cref="IBMDiffusiveCFLConstraint"/> instead of
        /// <see cref="Diffusion.DiffusiveCFLConstraint"/>
        /// </summary>
        /// <returns></returns>
        public override Operator GetDiffusiveOperator() {
            Operator op = new Operator(control);
            if (diffusiveFluxBuilder != null) {
                diffusiveFluxBuilder.BuildFluxes(op);
            }

            if (!op.IsEmpty) {
                op.CFLConstraints.Add(new IBMDiffusiveCFLConstraint(
                    control,
                    gridData,
                    workingSet,
                    speciesMap));
            }

            return op;
        }
    }
}

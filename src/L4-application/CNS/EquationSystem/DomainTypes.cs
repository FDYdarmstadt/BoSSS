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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.Convection;
using CNS.IBM;
using System;

namespace CNS.EquationSystem {

    /// <summary>
    /// Supported types of domains
    /// </summary>
    public enum DomainTypes {

        /// <summary>
        /// Standard domain bounded by boundary-fitted cells
        /// </summary>
        Standard,

        /// <summary>
        /// Domain immersed into a fixed background mesh, where the boundary of
        /// the domain is defined by the zero level set of a <b>fixed</b> level
        /// set function. 
        /// </summary>
        StaticImmersedBoundary,

        /// <summary>
        /// Domain immersed into a fixed background mesh, where the boundary of
        /// the domain is defined by the zero level set of an
        /// <b>in-stationary</b> level set function. 
        /// </summary>
        MovingImmersedBoundary
    }

    /// <summary>
    /// Extension methods for <see cref="DomainTypes"/>.
    /// </summary>
    public static class DomainTypesExtensions {

        /// <summary>
        /// Creates the appropriate instance of <see cref="CNSFieldSet"/>
        /// depending on the given <paramref name="domainType"/>
        /// </summary>
        /// <param name="domainType"></param>
        /// <param name="gridData"></param>
        /// <param name="control"></param>
        /// <returns></returns>
        public static CNSFieldSet CreateWorkingSet(this DomainTypes domainType, IGridData gridData, CNSControl control) {
            switch (domainType) {
                case DomainTypes.Standard:
                    return new CNSFieldSet(gridData, control);

                case DomainTypes.StaticImmersedBoundary:
                case DomainTypes.MovingImmersedBoundary:
                    // Make level set gradient exists in list of derived fields
                    for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                        if (!control.FieldOptions.ContainsKey(IBMVariables.LevelSetGradient[d])) {
                            control.AddVariable(
                                IBMVariables.LevelSetGradient[d],
                                control.VariableToDegreeMap[IBMVariables.LevelSet] - 1);
                        }
                    }
                    return new IBMFieldSet(gridData, (IBMControl)control);

                default:
                    throw new Exception();
            }
        }

        /// <summary>
        /// Returns the correct species map for the current domain type
        /// </summary>
        /// <param name="domainType">
        /// The considered domain type
        /// </param>
        /// <param name="workingSet"></param>
        /// <param name="control"></param>
        /// <param name="gridData"></param>
        /// <returns>
        /// A species map that is suitable for the current application.
        /// </returns>
        public static ISpeciesMap CreateSpeciesMap(this DomainTypes domainType, CNSFieldSet workingSet, CNSControl control, IGridData gridData) {
            Material material = new Material(control.EquationOfState, control.ViscosityLaw, control.MachNumber, control.ReynoldsNumber, control.PrandtlNumber, control.FroudeNumber, control.ViscosityRatio);

            switch (domainType) {
                case DomainTypes.Standard:
                    return new SingleSpeciesMap(gridData, material);

                case DomainTypes.StaticImmersedBoundary:
                case DomainTypes.MovingImmersedBoundary:
                    IBMFieldSet ibmWorkingSet = workingSet as IBMFieldSet;
                    return new ImmersedSpeciesMap(
                        (IBMControl)control, ibmWorkingSet.LevelSet, material);

                default:
                    throw new Exception();
            }
        }

        /// <summary>
        /// Creates an <see cref="OperatorFactory"/> with appropriate terms
        /// (convective/diffusive) for the selected <paramref name="formulation"/>.
        /// </summary>
        /// <param name="formulation">
        /// The chosen equation system
        /// </param>
        /// <param name="control"></param>
        /// <param name="gridData"></param>
        /// <param name="speciesMap"></param>
        /// <param name="workingSet"></param>
        /// <param name="boundaryMap">
        /// Boundary information
        /// </param>
        /// <returns>
        /// An instance of <see cref="OperatorFactory"/> that has been
        /// configured with the fluxes defined by
        /// <see cref="CNSControl.ConvectiveFluxType"/> and/or
        /// <see cref="CNSControl.DiffusiveFluxType"/>.
        /// </returns>
        public static OperatorFactory GetOperatorFactory(
            this DomainTypes formulation, CNSControl control, IGridData gridData, BoundaryConditionMap boundaryMap, CNSFieldSet workingSet, ISpeciesMap speciesMap) {
            switch (formulation) {
                case DomainTypes.Standard:
                    return new OperatorFactory(
                        control, gridData, workingSet, speciesMap, boundaryMap);

                case DomainTypes.StaticImmersedBoundary:
                case DomainTypes.MovingImmersedBoundary:
                    FluxBuilder convectiveBuilder = control.ConvectiveFluxType.GetBuilder(
                        control, boundaryMap, speciesMap);
                    return new IBMOperatorFactory(
                        (IBMControl)control,
                        gridData,
                        workingSet,
                        speciesMap,
                        boundaryMap);

                default:
                    throw new Exception(
                        "Unknown formulation \"" + control.DomainType + "\"");
            }
        }
    }
}

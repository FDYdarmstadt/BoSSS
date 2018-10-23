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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.Convection;
using CNS.Diffusion;
using CNS.IBM;
using CNS.ShockCapturing;
using CNS.Source;
using System;
using System.Collections.Generic;

namespace CNS.EquationSystem {

    /// <summary>
    /// Constructs the <see cref="SpatialOperator"/>s for derivatives of the
    /// Navier-Stokes equation in primitive variable formulation. The operators
    /// consist of convective terms (depending first derivatives only),
    /// diffusive terms (depending on first and second derivatives) and generic
    /// source terms that are defined by the flux associated with the supplied
    /// <see cref="FluxBuilder"/>s. For the available options, see
    /// <see cref="DomainTypes"/>
    /// </summary>
    public class OperatorFactory {

        /// <summary>
        /// Control file options
        /// </summary>
        protected readonly CNSControl control;

        /// <summary>
        /// Information about the grid
        /// </summary>
        protected readonly IGridData gridData;

        /// <summary>
        /// The DG fields
        /// </summary>
        protected readonly CNSFieldSet workingSet;

        /// <summary>
        /// Information about different species.
        /// </summary>
        protected readonly ISpeciesMap speciesMap;

        /// <summary>
        /// The builder for the convective fluxes (if required, i.e. for the
        /// Euler and the Navier-Stokes equations)
        /// </summary>
        protected readonly FluxBuilder convectiveFluxBuilder;

        /// <summary>
        /// The builder for the diffusive fluxes (if required, i.e. for the
        /// Stokes and Navier-Stokes equations)
        /// </summary>
        protected readonly FluxBuilder diffusiveFluxBuilder;

        /// <summary>
        /// The builders for the optional source terms
        /// </summary>
        protected readonly IList<FluxBuilder> sourceTermBuilders = new List<FluxBuilder>();

        /// <summary>
        /// Prepares the construction of a new equation system
        /// </summary>
        /// <param name="control"></param>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        /// <param name="speciesMap"></param>
        /// <param name="boundaryMap"></param>
        /// <remarks>
        /// Source terms are currently considered independent of the considered
        /// equation system and are thus constructed automatically from the
        /// control file (see <see cref="CustomSourceBuilder"/>).
        /// </remarks>
        public OperatorFactory(
            CNSControl control,
            IGridData gridData,
            CNSFieldSet workingSet,
            ISpeciesMap speciesMap,
            BoundaryConditionMap boundaryMap) {

            bool hasConvection = control.ActiveOperators.HasFlag(Operators.Convection);
            bool hasDiffusion = control.ActiveOperators.HasFlag(Operators.Diffusion);
            if (!hasConvection && !hasDiffusion) {
                throw new Exception(
                    "Either convective or diffusive terms must be active");
            }

            if (hasConvection) {
                this.convectiveFluxBuilder = control.ConvectiveFluxType.GetBuilder(
                    control, boundaryMap, speciesMap);
            }

            if (hasDiffusion) {
                this.diffusiveFluxBuilder = control.DiffusiveFluxType.GetBuilder(
                    control, boundaryMap, speciesMap, gridData);
            }

            if (control.ActiveOperators.HasFlag(Operators.Gravity)) {
                this.sourceTermBuilders.Add(
                    new GravityFluxBuilder(control, boundaryMap, speciesMap));
            }

            if (control.ActiveOperators.HasFlag(Operators.CustomSource)) {
                this.sourceTermBuilders.Add(
                    new CustomSourceBuilder(control, boundaryMap, speciesMap));
            }

            if (control.ActiveOperators.HasFlag(Operators.SpongeLayer)) {
                this.sourceTermBuilders.Add(
                    new SpongeLayerFluxBuilder(control, boundaryMap, speciesMap));
            }

            if (control.ActiveOperators.HasFlag(Operators.ArtificialViscosity)) {
                this.sourceTermBuilders.Add(
                    new LaplacianArtificialViscosityFluxBuilder(control, boundaryMap, speciesMap));
            }

            this.control = control;
            this.gridData = gridData;
            this.workingSet = workingSet;
            this.speciesMap = speciesMap;
        }

        /// <summary>
        /// Constructs the convective operator of the equation system.
        /// </summary>
        /// <returns>
        /// A new <see cref="SpatialOperator"/> with all fluxes
        /// defined by <see cref="convectiveFluxBuilder"/>.
        /// </returns>
        public virtual Operator GetConvectiveOperator() {
            Operator op = new Operator(control);
            if (convectiveFluxBuilder != null) {
                convectiveFluxBuilder.BuildFluxes(op);
            }

            if (!op.IsEmpty) {
                op.CFLConstraints.Add(new ConvectiveCFLConstraint(
                    control,
                    gridData,
                    workingSet,
                    speciesMap));
            }

            return op;
        }

        /// <summary>
        /// Constructs the diffusive operator of the equation system.
        /// </summary>
        /// <returns>
        /// A new <see cref="SpatialOperator"/> with all fluxes
        /// defined by <see cref="diffusiveFluxBuilder"/>.
        /// </returns>
        public virtual Operator GetDiffusiveOperator() {
            Operator op = new Operator(control);
            if (diffusiveFluxBuilder != null) {
                diffusiveFluxBuilder.BuildFluxes(op);
            }

            if (!op.IsEmpty) {
                op.CFLConstraints.Add(new DiffusiveCFLConstraint(
                    control,
                    gridData,
                    workingSet,
                    speciesMap));
            }

            return op;
        }

        /// <summary>
        /// Creates a joined operator consisting of all the convective and the
        /// diffusive fluxes defined by <see cref="convectiveFluxBuilder"/> and
        /// <see cref="diffusiveFluxBuilder"/>. Most useful for fully explicit
        /// schemes where a splitting between linear and nonlinear terms is
        /// futile.
        /// </summary>
        /// <returns>
        /// An operator for the full equations defined by this object
        /// </returns>
        public virtual Operator GetJoinedOperator() {
            return GetConvectiveOperator().
                Union(GetDiffusiveOperator()).
                Union(GetSourceTermOperator());
        }

        /// <summary>
        /// Constructs the operator containing all source term of the
        /// underlying equation system. By default, source terms are
        /// constructed from the configuration, sees
        /// <see cref="CustomSourceBuilder"/>.
        /// </summary>
        /// <returns>
        /// A spatial operator containing source terms.
        /// </returns>
        public virtual Operator GetSourceTermOperator() {
            Operator op = new Operator(control);
            foreach (FluxBuilder builder in sourceTermBuilders) {
                builder.BuildFluxes(op);
            }

            // TODO add IBM AV CFL constraint
            if (control.ActiveOperators.HasFlag(Operators.ArtificialViscosity)) {
                op.CFLConstraints.Add(new ArtificialViscosityCFLConstraint(
                    control,
                    gridData,
                    workingSet,
                    speciesMap));
            }

            return op;
        }
    }
}
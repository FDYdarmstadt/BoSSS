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
using ilPSP;
using BoSSS.Solution;
using System;

namespace CNS.EquationSystem {

    /// <summary>
    /// A generic operator that helps with the construction of a
    /// <see cref="SpatialOperator"/>. The main purpose of this class is the
    /// simplification of the process of error-prone merging different
    /// operators (like the convective and the diffusive operator)
    /// </summary>
    public class Operator {

        /// <summary>
        /// Configuration options
        /// </summary>
        private CNSControl config;

        /// <summary>
        /// The parameter ordering that is constructed from the parameter
        /// orderings of the different components.
        /// </summary>
        private IList<string> GetParameterOrdering(CNSFieldSet fieldSet) {
            //return DensityComponents.
            //    SelectMany(f => f.ParameterOrdering ?? new string[0]).
            //    Union(MomentumComponents.SelectMany(f => f.SelectMany(g => g.ParameterOrdering ?? new string[0]))).
            //    Union(EnergyComponents.SelectMany(f => f.ParameterOrdering ?? new string[0])).
            //    ToList();
            return fieldSet.ParameterFields.Select(f => f.Identification).ToList();
        }

        /// <summary>
        /// Constructs and empty operator.
        /// </summary>
        public Operator(CNSControl config) {
            this.config = config;

            DensityComponents = new List<IEquationComponent>();
            MomentumComponents = new IList<IEquationComponent>[
                CNSEnvironment.NumberOfDimensions];
            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                MomentumComponents[d] = new List<IEquationComponent>();
            }
            EnergyComponents = new List<IEquationComponent>();
            CFLConstraints = new List<TimeStepConstraint>();
        }

        /// <summary>
        /// A list of components related to the continuity equation.
        /// </summary>
        public IList<IEquationComponent> DensityComponents {
            get;
            private set;
        }

        /// <summary>
        /// A list of lists of components related to the [1-3] momentum
        /// equations.
        /// </summary>
        public IList<IEquationComponent>[] MomentumComponents {
            get;
            private set;
        }

        /// <summary>
        /// A list of components related to the energy equation.
        /// </summary>
        public IList<IEquationComponent> EnergyComponents {
            get;
            private set;
        }

        /// <summary>
        /// A list of CFL constraint provides (e.g., stemming from different
        /// operators that have been merged) that limit the admissible
        /// step-size for explicit time-integrators.
        /// </summary>
        public IList<TimeStepConstraint> CFLConstraints {
            get;
            private set;
        }

        /// <summary>
        /// Determines whether any components have been registered
        /// </summary>
        public bool IsEmpty {
            get {
                return !DensityComponents.Any() &&
                    !MomentumComponents.Any(m => m.Any()) &&
                    !EnergyComponents.Any();
            }
        }

        /// <summary>
        /// Merges <paramref name="other"/> into this operator but adding the
        /// respective compatible components (see
        /// <see cref="DensityComponents"/>, <see cref="MomentumComponents"/>
        /// and <see cref="EnergyComponents"/>) as well as the defined CFL
        /// constraints (if any) to this object
        /// </summary>
        /// <param name="other"></param>
        /// <returns>
        /// A new operator containing the components and CFL constraints of
        /// this object and the merged object <paramref name="other"/>.
        /// </returns>
        public Operator Union(Operator other) {
            Operator union = new Operator(other.config);
            union.DensityComponents =
                this.DensityComponents.Concat(other.DensityComponents).ToList();
            // It might not be readable, but lambdas are soooo tempting :)
            union.MomentumComponents = this.MomentumComponents.Select(
                (l, i) => l.Concat(other.MomentumComponents[i]).ToList()).ToArray();
            union.EnergyComponents =
                this.EnergyComponents.Concat(other.EnergyComponents).ToList();

            union.CFLConstraints =
                this.CFLConstraints.Concat(other.CFLConstraints).ToList();

            return union;
        }

        /// <summary>
        /// Creates a <see cref="SpatialOperator"/> with the equation
        /// components that have been assigned to this operator.
        /// </summary>
        /// <returns></returns>
        public SpatialOperator ToSpatialOperator(CNSFieldSet fieldSet) {
            SpatialOperator spatialOp = new SpatialOperator(
                CNSEnvironment.PrimalArgumentOrdering,
                GetParameterOrdering(fieldSet),
                CNSEnvironment.PrimalArgumentOrdering,
                QuadOrderFunc.NonLinearWithoutParameters(2));
            MapComponents(spatialOp);
            spatialOp.Commit();
            return spatialOp;
        }
        /// <summary>
        /// Maps the <see cref="IEquationComponent"/>s in
        /// <paramref name="op"/> to the relevant equation components in
        /// the <paramref name="op"/>.
        /// </summary>
        /// <param name="op">
        /// The operator onto which the components of
        /// <paramref name="op"/> should be mapped.
        /// </param>
        private void MapComponents(SpatialOperator op) {
            DensityComponents.ForEach(component =>
                op.EquationComponents[Variables.Density].Add(component));

            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                MomentumComponents[d].ForEach(component =>
                    op.EquationComponents[Variables.Momentum[d]].Add(component));
            }

            EnergyComponents.ForEach(component =>
                op.EquationComponents[Variables.Energy].Add(component));
        }
    }
}

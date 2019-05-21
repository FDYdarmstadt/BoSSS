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

using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;

namespace CNS.Source {

    /// <summary>
    /// Implementation of a generic source term given by a user-defined formula
    /// that may involve nonlinear terms
    /// </summary>
    public class AdHocSourceTerm : NonlinearSource {

        /// <summary>
        /// See constructor
        /// </summary>
        private ISpeciesMap speciesMap;

        /// <summary>
        /// <see cref="ArgumentOrdering"/>
        /// </summary>
        private string[] argumentOrdering;

        /// <summary>
        /// <see cref="AdHocSourceTerm.AdHocSourceTerm"/>
        /// </summary>
        private Func<double[], double, StateVector, double> formula;

        /// <summary>
        /// Constructs a source term represented by the given
        /// <paramref name="formula"/>.
        /// </summary>
        /// <param name="speciesMap"></param>
        /// <param name="formula">
        /// A formula representing the actual source term.
        /// </param>
        public AdHocSourceTerm(ISpeciesMap speciesMap, Func<double[], double, StateVector, double> formula) {
            this.speciesMap = speciesMap;
            this.formula = formula;
            this.argumentOrdering = CompressibleEnvironment.PrimalArgumentOrdering;
        }

        /// <summary>
        /// Evaluates the formula supplied to the constructor.
        /// </summary>
        /// <param name="time">
        /// <see cref="NonlinearSource.Source(double, int, double[], double[])"/>
        /// </param>
        /// <param name="j">
        /// <see cref="NonlinearSource.Source(double, int, double[], double[])"/>
        /// </param>
        /// <param name="x">
        /// <see cref="NonlinearSource.Source(double, int, double[], double[])"/>
        /// </param>
        /// <param name="U">
        /// <see cref="NonlinearSource.Source(double, int, double[], double[])"/>
        /// </param>
        /// <returns>
        /// The value returned by the supplied formula
        /// </returns>
        protected override double Source(double time, int j, double[] x, double[] U) {
            StateVector state = new StateVector(U, speciesMap.GetMaterial(double.NaN));
            return formula(x, time, state);
        }

        /// <summary>
        /// <see cref="NonlinearSource.ArgumentOrdering"/>
        /// </summary>
        public override IList<string> ArgumentOrdering {
            get {
                return argumentOrdering;
            }
        }
    }
}

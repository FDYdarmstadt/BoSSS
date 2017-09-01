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


using System;
using BoSSS.Foundation;

namespace CNS.Solution {

    /// <summary>
    /// Available schemes.
    /// </summary>
    public enum ImplicitSchemes {

        /// <summary>
        /// No implicit treatment
        /// </summary>
        None = 0,

        /// <summary>
        /// No time integration, only the steady state equations are
        /// considered
        /// </summary>
        Steady,

        /// <summary>
        /// Implicit Euler discretization of the temporal term. Temporal error
        /// should be of first order.
        /// </summary>
        ImplicitEuler,

        /// <summary>
        /// Crank-Nicolson discretization of the temporal term. Temporal error
        /// should be of second order, but the scheme is not unconditionally
        /// stable
        /// </summary>
        CrankNicolson,

        /// <summary>
        /// Temporal term is discretized with back differentiation formula of
        /// second order. Temporal error should be of second order.
        /// </summary>
        BDF2
    }

    /// <summary>
    /// Extension methods <see cref="ImplicitSchemes"/>.
    /// </summary>
    public static class ImplicitSchemesExtensions {

        /// <summary>
        /// Instantiates the appropriate instance of an implicit scheme
        /// according to <paramref name="scheme"/>.
        /// </summary>
        /// <param name="scheme">
        /// The selected implicit scheme
        /// </param>
        /// <param name="config">
        /// Configuration options
        /// </param>
        /// <param name="op">
        /// The operator that should be used in the solution process
        /// </param>
        /// <param name="mapping">
        /// Affected DG fields.
        /// </param>
        /// <param name="parameterMapping">
        /// DG fields serving as parameters for the spatial operator.
        /// </param>
        /// <returns>
        /// An instance of the selected implicit scheme
        /// </returns>
        public static IIterativeImplicitScheme Instantiate(this ImplicitSchemes scheme, CNSControl config, SpatialOperator op, CoordinateMapping mapping, CoordinateMapping parameterMapping) {
            switch (scheme) {
                case ImplicitSchemes.None:
                    throw new System.ArgumentException(
                        "Cannot instantiate empty scheme");

                case ImplicitSchemes.Steady:
                    return config.NonlinearSolver.Instantiate(
                        config,
                        op,
                        config.ImplicitSolverFactory(),
                        mapping,
                        parameterMapping);

                default:
                    throw new System.NotImplementedException(String.Format(
                        "Scheme \"{0}\" has not been implemented yet",
                        scheme));
            }
        }
    }
}

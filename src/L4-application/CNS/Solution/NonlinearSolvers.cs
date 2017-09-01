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
using ilPSP.LinSolvers;

namespace CNS.Solution {

    /// <summary>
    /// Available options for nonlinear solvers, i.e. implementations of 
    /// <see cref="IIterativeImplicitScheme"/>.
    /// </summary>
    public enum NonlinearSolvers {

        /// <summary>
        /// No nonlinear solver.
        /// </summary>
        None = 0,

        /// <summary>
        /// <see cref="FixedPointIteration"/>
        /// </summary>
        FixedPointIteration
    }

    /// <summary>
    /// Extension methods for <see cref="NonlinearSolvers"/>
    /// </summary>
    public static class NonlinearSolversExtensions {

        /// <summary>
        /// Instantiates the nonlinear solver indicated by
        /// <paramref name="solver"/>.
        /// </summary>
        /// <param name="solver">
        /// The type of the solver to be instantiated-
        /// </param>
        /// <param name="config">
        /// Configuration options.
        /// </param>
        /// <param name="originalOperator">
        /// The (nonlinear) spatial differential operator to be solved.
        /// </param>
        /// <param name="linearSolver">
        /// Linear solver required in the iterations of
        /// the nonlinear solver.
        /// </param>
        /// <param name="mapping">
        /// DG fields affected by the solution procedure.
        /// </param>
        /// <param name="parameterMapping">
        /// DG fields serving as parameters for the spatial operator.
        /// </param>
        /// <returns>
        /// An instance of an implementation of <see cref="IIterativeImplicitScheme"/>
        /// </returns>
        public static IIterativeImplicitScheme Instantiate(
            this NonlinearSolvers solver, CNSControl config, SpatialOperator originalOperator, ISparseSolver linearSolver, CoordinateMapping mapping, CoordinateMapping parameterMapping) {

            switch (solver) {
                case NonlinearSolvers.FixedPointIteration:
                    double damping = config.FixedPointDamping;
                    if (damping <= 0.0) {
                        throw new Exception.ConfigurationException(
                            "Damping factor must not be less than or equal to zero");
                    }

                    return new FixedPointIteration(
                        originalOperator, linearSolver, mapping, parameterMapping, damping);

                case NonlinearSolvers.None:
                    return null;

                default:
                    throw new System.NotImplementedException(
                        "Unknown nonlinear scheme \"" + solver + "\"");
            }
        }
    }
}

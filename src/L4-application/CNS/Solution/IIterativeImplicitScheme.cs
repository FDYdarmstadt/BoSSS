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
using BoSSS.Solution;

namespace CNS.Solution {

    /// <summary>
    /// Common interface for all iterative solvers for nonlinear systems of
    /// equations.
    /// </summary>
    public interface IIterativeImplicitScheme : ITimeStepper, IDisposable {

        /// <summary>
        /// The current iteration
        /// </summary>
        int CurrentIteration {
            get;
        }

        /// <summary>
        /// The spatial operator used by the nonlinear solution algorithm
        /// </summary>
        SpatialOperator Operator {
            get;
        }

        /// <summary>
        /// Performs a single iteration of the nonlinear scheme.
        /// </summary>
        void PerformIteration();
    }
}

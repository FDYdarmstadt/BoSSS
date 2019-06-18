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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ilPSP.LinSolvers;

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// The basic interface for a linear solver. 
    /// </summary>
    public interface ISolverSmootherTemplate {


        /// <summary>
        /// Initializes the linear solver.
        /// </summary>
        /// <param name="op">
        /// Provides the matrix for the solver.
        /// </param>
        void Init(MultigridOperator op);
 
        /// <summary>
        /// Call to solver/smoother.
        /// </summary>
        /// <param name="X">
        /// On entry: initial guess for the solution; on exit: aproximate solution after applying the solver/smoother <paramref name="iter"/> times.
        /// </param>
        /// <param name="B">
        /// On entry: the right-hand-side of the system.
        /// </param>
        void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double>;


        /// <summary>
        /// Total number of iterations (since the last call to <see cref="ResetStat"/>) in nested solvers.
        /// </summary>
        int IterationsInNested {
            get;
        }

        /// <summary>
        /// Total number of iterations (since the last call to <see cref="ResetStat"/>) on the multi-grid level associated with this solver.
        /// </summary>
        int ThisLevelIterations {
            get;
        }

        /// <summary>
        /// True, if the last call to <see cref="Solve{U, V}(U, V)"/> was successful.
        /// </summary>
        bool Converged {
            get;
        }
        
        /// <summary>
        /// Reset solver statistics, i.e. set <see cref="IterationsInNested"/>, <see cref="ThisLevelIterations"/>, etc. to 0.
        /// </summary>
        void ResetStat();

        ISolverSmootherTemplate Clone();
    }

    /// <summary>
    /// Interface for solvers which provide a call-back during the solver-iterations.
    /// </summary>
    public interface ISolverWithCallback : ISolverSmootherTemplate {

        /// <summary>
        ///  - 1st argument: iteration index<br/>
        ///  - 2nd argument: current solution<br/>
        ///  - 3rd argument: current residual<br/>
        ///  - 4th argument: the currently used operator<br/>
        /// </summary>      
        Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }
    }
}

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
using BoSSS.Solution.Control;
using ilPSP.LinSolvers;

namespace BoSSS.Solution.AdvancedSolvers {

    /// <summary>
    /// The basic interface for a linear solver. 
    /// </summary>
    public interface ISolverSmootherTemplate : ICloneable, IDisposable {

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
        /// On entry: initial guess for the solution; on exit: approximate solution after applying the solver/smoother.
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

        /// <summary>
        /// Estimate of used memory in bytes
        /// </summary>
        long UsedMemory();
    }

    /// <summary>
    /// Minimal information which a solver/smoother might require; This is
    /// - obviously, the matrix <see cref="OperatorMatrix"/>
    /// - minimal information about the underlying GD discretization, see <see cref="DgMapping"/>
    /// </summary>
    /// <remarks>
    /// Feb22: This was designed in an attempt to make the solvers defined by <see cref="ISolverSmootherTemplate"/>,
    /// which always work on the **entire** mesh and system usable on **sub-systems**.
    /// At this point, the functionality in <see cref="SubBlockSelector"/> and <see cref="BlockMask"/> allows
    /// to extract sub-matrices of the full system matrix. However, it is not possible to apply a <see cref="ISolverSmootherTemplate"/>
    /// onto this matrix, since those solvers require a full <see cref="MultigridOperator"/> 
    /// to be initialized (<see cref="ISolverSmootherTemplate.Init(MultigridOperator)"/>.
    /// </remarks>
    public interface IOperatorMappingPair {
        
        /// <summary>
        /// Returns the Operator matrix on the current multigrid level.
        /// </summary>
        BlockMsrMatrix OperatorMatrix {
            get;
        }

        /// <summary>
        /// Minimal information about the DG structure behind matrix <see cref="OperatorMatrix"/>
        /// </summary>
        ICoordinateMapping DgMapping {
            get;
        }

    }


    /// <summary>
    /// Defines a solver which is able to work on sub-system, i.e. some part of the matrix on the full mesh.
    /// </summary>
    public interface ISubsystemSolver : ISolverSmootherTemplate {
        
        /// <summary>
        /// Initializes the linear solver.
        /// </summary>
        /// <param name="op">
        /// Provides the matrix for the solver.
        /// </param>
        void Init(IOperatorMappingPair op);
        
    }
  

    /// <summary>
    /// 
    /// </summary>
    public interface ISolverFactory : IEquatable<ISolverFactory> {

        /// <summary>
        /// Name/Description of the solver configuration
        /// </summary>
        string Name { get; }

        /// <summary>
        /// short version of <see cref="Name"/>, e.g. to be used in plots or tables
        /// </summary>
        string Shortname { get; }

        /// <summary>
        /// Creates an Instance of the respective solver
        /// </summary>
        ISolverSmootherTemplate CreateInstance(MultigridOperator level);
    }



    /// <summary>
    /// For certain solvers, a programmable termination criterion seems handy.
    /// On the other hand, preconditioners most of the time run on a fixed number of iterations.
    /// </summary>
    public interface IProgrammableTermination : ISolverSmootherTemplate {
        
        /// <summary>
        /// User-Programmable termination criterion: 
        /// - 1st argument: iteration index
        /// - 2nd argument: l2-Norm of residual of initial solution 
        /// - 3rd argument: l2-Norm of residual of solution in current iteration
        /// - return value, 1st item: true to continue, false to terminate
        /// - return value, 2nd item: true for successful convergence (e.g. convergence criterion reached), false for failure (e.g. maximum number of iterations reached)
        /// </summary>
        /// <remarks>
        /// For an example implementation, see
        /// <see cref="IterativeSolverConfig.DefaultTermination"/>.
        /// </remarks>
        Func<int, double, double, (bool bNotTerminate, bool bSuccess)> TerminationCriterion {
            get;
            set;
        }
    }


    /// <summary>
    /// Interface for solvers which provide a call-back during the solver-iterations.
    /// </summary>
    public interface ISolverWithCallback : ISolverSmootherTemplate {

        /// <summary>
        ///  - 1st argument: iteration index
        ///  - 2nd argument: current solution
        ///  - 3rd argument: current residual
        ///  - 4th argument: the currently used operator
        /// </summary>      
        Action<int, double[], double[], MultigridOperator> IterationCallback {
            get;
            set;
        }
    }
}

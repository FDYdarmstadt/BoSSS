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
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution;

namespace BoSSS.Solution.Control {
    public class NonLinearSolverConfig {
        public enum Code {

            /// <summary>
            /// NewtonKrylov GMRES (<see cref="BoSSS.Solution.AdcancedSolvers.NonLinearSolver"/>) with linear solver (<see cref="LinearSolverConfig.Code"/>) used as preconditioner for matrix-free GMRES 
            /// </summary>
            NewtonGMRES = 0,

            /// <summary>
            /// Picard fixpoint solver (<see cref="BoSSS.Solution.AdvancedSolvers.FixpointIterator"/>) with linear solver (<see cref="LinearSolverConfig.Code"/>) for the linearized equation system
            /// </summary>
            Picard = 1,

            /// <summary>
            /// Newtons method (<see cref="BoSSS.Solution.AdvancedSolvers.Newton"/>) with linear solver (<see cref="LinearSolverConfig.Code"/>) used to approximate the inverse of the jacobian with the inverse operator matrix. 
            /// </summary>
            Newton = 2,

            PicardGMRES = 3,

        }

        public class _Precond : LinearSolverConfig {}

        public  _Precond Precond;

        ///// <summary>
        ///// If iterative saddle-point solvers like GMRES or Orthonormalization are used, the maximum number of basis vectors
        ///// that are used to construct the accelerated solution.
        ///// </summary>
        //[DataMember]
        //public int MaxKrylovDim = 30;

        /// <summary>
        /// If iterative solvers are used, the maximum number of iterations.
        /// </summary>
        [DataMember]
        public int MaxSolverIterations = 2000;

        /// <summary>
        /// If iterative solvers are used, the maximum number of iterations.
        /// </summary>
        [DataMember]
        public int MinSolverIterations = 2;

        /// <summary>
        /// Convergence criterion for nonlinear solver.
        /// </summary>
        [DataMember]
        public double Solver_ConvergenceCriterion = 1.0e-8;

        /// <summary>
        /// Sets the algorithm to use for nonlinear solving, e.g. Newton or Picard.
        /// </summary>
        [DataMember]
        public NonLinearSolverConfig.Code SolverCode = NonLinearSolverConfig.Code.Picard;
    }

}

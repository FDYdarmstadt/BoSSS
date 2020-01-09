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

    public enum NonLinearSolverCode {

        /// <summary>
        /// The bald guy from the Enterprise.
        /// Picard fixpoint solver (<see cref="BoSSS.Solution.AdvancedSolvers.FixpointIterator"/>) with linear solver (<see cref="LinearSolverCode"/>) for the linearized equation system
        /// </summary>
        Picard = 1,

        /// <summary>
        /// Newtons method (<see cref="BoSSS.Solution.AdvancedSolvers.Newton"/>) with linear solver (<see cref="LinearSolverCode"/>) used to approximate the inverse of the jacobian with the inverse operator matrix. 
        /// </summary>
        Newton = 2,

        /// <summary>
        /// Mixed sequence of nonlinear solvers. Useful for example if one wishes to start a simulation with Picard and then change to Newton.
        /// </summary>
        NLSolverSequence = 4, 

        /// <summary>
        /// Overgiven solver to Solver Chooser object. Weisse bescheid ...
        /// </summary>
        selfmade = 999,
    }

    /// <summary>
    /// User-Options for nonlinear solver configuration; 
    /// </summary>
    public class NonLinearSolverConfig : ICloneable, IEquatable<NonLinearSolverConfig>{

        /// <summary>
        /// This will print out more information about iterations.
        /// </summary>
        public bool verbose = false;

        
        /// <summary>
        /// If iterative solvers are used, the maximum number of iterations.
        /// </summary>
        [DataMember]
        public int MaxSolverIterations = 2000;

        /// <summary>
        /// If iterative solvers are used, the minimum number of iterations.
        /// </summary>
        [DataMember]
        public int MinSolverIterations = 2;

        /// <summary>
        /// Convergence criterion for nonlinear solver.
        /// </summary>
        [DataMember]
        public double ConvergenceCriterion = 1.0e-8;

        /// <summary>
        /// under relaxation, if a fixpoint iterator (e.g.Picard) is used.
        /// </summary>
        [DataMember]
        public double UnderRelax = 1;

        /// <summary>
        /// Sets the algorithm to use for nonlinear solving, e.g. Newton or Picard.
        /// </summary>
        [DataMember]
        public NonLinearSolverCode SolverCode = NonLinearSolverCode.Picard;

        /// <summary>
        /// Number of iterations, where jacobi is not updated. Also known as constant newton method. Default 1, means regular newton.
        /// </summary>
        [DataMember]
        public int constantNewtonIterations = 1;

        /// <summary>
        /// When pressure ref point can not be used (Option: false), a mean pressure field is calculated and substracted from the current pressure field.
        /// </summary>
        public bool UsePresRefPoint = true;

        /// <summary>
        /// Clones the NonLinearConfig
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            var clone = new NonLinearSolverConfig() {
                constantNewtonIterations = this.constantNewtonIterations,
                ConvergenceCriterion = this.ConvergenceCriterion,
                MaxSolverIterations = this.MaxSolverIterations,
                MinSolverIterations = this.MinSolverIterations,
                SolverCode = this.SolverCode,
                UnderRelax = this.UnderRelax,
                UsePresRefPoint = this.UsePresRefPoint,
                verbose = this.verbose
        };
            return clone;
        }

        /// <summary>
        /// Compares value not reference!
        /// </summary>
        /// <param name="compareto"></param>
        /// <returns></returns>
        public bool Equals(NonLinearSolverConfig compareto) {
            return this.constantNewtonIterations == compareto.constantNewtonIterations &&
                this.ConvergenceCriterion == compareto.ConvergenceCriterion &&
                this.MaxSolverIterations == compareto.MaxSolverIterations &&
                this.MinSolverIterations == compareto.MinSolverIterations &&
                this.SolverCode == compareto.SolverCode &&
                this.UnderRelax == compareto.UnderRelax &&
                this.UsePresRefPoint == compareto.UsePresRefPoint &&
                this.verbose == compareto.verbose;
        }
    }
}

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
using System.Threading.Tasks;
using BoSSS.Foundation;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Connectors.Matlab;
using ilPSP.Tracing;
using System.IO;
using System.Diagnostics;

namespace BoSSS.Solution.AdvancedSolvers {

    /*
    /// <summary>

    /// </summary>
    public class NLSolverSequence : NonlinearSolver {
        /// <summary>
        /// NonLinear solver sequence. Used for starting with Picard and later with Newton
        /// </summary>
        /// <param name="__AssembleMatrix"></param>
        /// <param name="__AggBasisSeq"></param>
        /// <param name="__MultigridOperatorConfig"></param>
        public NLSolverSequence(OperatorEvalOrLin __AssembleMatrix, IEnumerable<AggregationGridBasis[]> __AggBasisSeq, MultigridOperator.ChangeOfBasisConfig[][] __MultigridOperatorConfig) :
            base(__AssembleMatrix, __AggBasisSeq, __MultigridOperatorConfig) //
        {
        }
        /// <summary>
        /// sequence of NonLinear solvers
        /// </summary>
        public NonlinearSolver[] m_NLSequence;

        /// <summary>
        /// the linear solver
        /// </summary>
        public ISolverSmootherTemplate m_linsolver;

        /// <summary>
        /// Called at every iteration; the arguments are 
        ///  - iteration index 
        ///  - current solution 
        ///  - current residual 
        ///  - current multigrid operator
        /// </summary>
        //public event Action<int, double[], double[], MultigridOperator> IterationCallback;

        public override bool SolverDriver<S>(CoordinateVector SolutionVec, S RHS) {

            if (m_NLSequence == null)
                throw new FormatException("m_NLSequence is empty. You fool!");

            bool finalResult = false;
            foreach (NonlinearSolver NL in m_NLSequence) {

                //NL.IterationCallback += IterationCallback;
                finalResult = NL.SolverDriver(SolutionVec, RHS);
            }

            return finalResult;

        }
    }
    */
}




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

using NSE_SIMPLE.Multiphase;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace NSE_SIMPLE {

    /// <summary>
    /// Collection of matrix assemblies for level-set equation.
    /// </summary>
    public class MatrixFactoryLevelSet {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="LevelSetOperators"></param>
        /// <param name="multiphaseControl"></param>
        /// <param name="BDF"></param>
        public MatrixFactoryLevelSet(OperatorFactoryLevelSet LevelSetOperators, int LocalNoOfCells, SolverConfiguration solverConf, BDFScheme BDF) {
            LevelSet = new MatrixAssemblyLevelSet(LevelSetOperators.LevelSetAdvection);

            MultiphaseSIMPLEControl multiphaseControl = solverConf.Control as MultiphaseSIMPLEControl;
            if (multiphaseControl.LevelSetRelaxationType == RelaxationTypes.Implicit) {
                LevelSetApprox = new MatrixAssemblyApprox(
                    solverConf, LocalNoOfCells, LevelSet, BDF, 2 * multiphaseControl.PredictorApproximationUpdateCycle);
            }
        }

        /// <summary>
        /// Matrix for level-set equation, i.e. advection term.
        /// </summary>
        public SIMPLEMatrixAssembly LevelSet {
            get;
            private set;
        }

        /// <summary>
        /// Approximation of level-set matrix for under-relaxation.
        /// </summary>
        public SIMPLEMatrixAssembly LevelSetApprox {
            get;
            private set;
        }
    }
}

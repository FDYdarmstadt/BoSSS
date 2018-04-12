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
using System.Collections;
using BoSSS.Foundation;

namespace BoSSS.Solution.LevelSetTools.FastMarching {
    /// <summary>
    ///Interface for CellMarcher. The Local Solver reinitializes only one cell. 
    /// </summary>
    interface ILocalSolver {

        /// <summary>
        /// Initializes <paramref name="Phi"/> in cell <paramref name="jCell"/> using the cells from <paramref name="AcceptedMask"/>
        /// </summary>
        /// <param name="jCell"></param>
        /// <param name="AcceptedMask"></param>
        /// <param name="Phi"></param>
        void LocalSolve(int jCell, BitArray AcceptedMask, SinglePhaseField Phi);

    }
}

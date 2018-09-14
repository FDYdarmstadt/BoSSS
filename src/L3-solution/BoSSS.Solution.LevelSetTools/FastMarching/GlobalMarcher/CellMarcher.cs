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
using System.Collections;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.LevelSetTools.FastMarcher;
using BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher;

namespace BoSSS.Solution.LevelSetTools.FastMarching.GlobalMarcher {

    class CellMarcher {

        GridData gridDat;
        ILocalSolver localSolver;
        
        /// <summary>
        /// Fast marching solver. Initializes a Domain by fast marching. 
        /// Each cell must be initialized locally with a <paramref name="LocalSolver"/>.
        /// </summary>
        /// <param name="LevelSetBasis"></param>
        /// <param name="LocalSolver"> A solver that initializes only one cell</param>
        public CellMarcher(Basis LevelSetBasis, ILocalSolver LocalSolver) {
            gridDat = (GridData)(LevelSetBasis.GridDat);
            localSolver = LocalSolver;
        }

        /// <summary>
        /// Firstorder Reinit of <paramref name="Phi"/> on <paramref name="ReinitField"/> field.
        /// The order in which each cell is initialized is determined by fast marching. 
        /// Locally, each cell is initialized by the localSolver specified in the constructor.
        /// </summary>
        /// <param name="Phi">Field to reinitialize</param>
        /// <param name="Accepted">Start values</param>
        /// <param name="ReinitField">Specific Domain, e.g. whole domain or nearField</param>
        public void Reinit(SinglePhaseField Phi, CellMask Accepted, CellMask ReinitField) {

            //Build Marcher 
            IFastMarchingQueue<IMarchingNode> Heap = new MarchingHeap(this.gridDat.Cells.Count);
            Fastmarcher Solver = new Fastmarcher(Heap);

            //Initialize Graph for Marching and build initial accepted nodes
            MarchingCell.Initialize(localSolver, Phi, gridDat, ReinitField);
            MarchingCell[] AcceptedCells = MarchingCell.BuildInitialAcceptedCells(Accepted);

            //Solve
            Solver.march(AcceptedCells);

        }
    }
}

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
using BoSSS.Foundation.Grid;
using CNS.EquationSystem;
using System;

namespace CNS.LoadBalancing {

    /// <summary>
    /// AV cells are 1, all others are 0
    /// </summary>
    public class ArtificialViscosityCellClassifier : ICellClassifier {

        public (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program) {
            if (!program.Control.ActiveOperators.HasFlag(Operators.ArtificialViscosity)) {
                throw new Exception(
                    "The selected cell classifier is only sensible for runs with artificial viscosity");
            }

            int noOfClasses = 2;
            int J = program.gridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cellToPerformanceClassMap = new int[J];

            // old
            //foreach (Chunk chunk in program.Control.ArtificialViscosityLaw.GetShockedCellMask(program.GridData)) {
            //    foreach (int cell in chunk.Elements) {
            //        cellToPerformanceClassMap[cell] = 1;
            //    }
            //}

            // new
            DGField avField = program.WorkingSet.DerivedFields[CNSVariables.ArtificialViscosity];
            foreach (Chunk chunk in program.SpeciesMap.SubGrid.VolumeMask) {
                foreach (int cell in chunk.Elements) {
                    if (avField.GetMeanValue(cell) > 0) {
                        cellToPerformanceClassMap[cell] = 1;
                    }
                }
            }

            return (noOfClasses, cellToPerformanceClassMap);
        }
    }
}

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

using BoSSS.Foundation.Grid;
using BoSSS.Solution.CompressibleFlowCommon;
using System.Collections;

namespace CNS.ShockCapturing {

    /// <summary>
    /// A law to determine the appropriate amount of artificial viscosity in
    /// troubled cells
    /// </summary>
    public interface IArtificialViscosityLaw {

        /// <summary>
        /// Returns the appropriate viscosity.
        /// </summary>
        /// <param name="cellSize"></param>
        /// <param name="jCell"></param>
        /// <param name="state"></param>
        /// <returns></returns>
        double GetViscosity(int jCell, double cellSize, StateVector state);

        bool IsShocked(int jCell);
    }

    /// <summary>
    /// Extension methods for <see cref="IArtificialViscosityLaw"/>
    /// </summary>
    public static class IArtificialViscosityLawExtensions {

        /// <summary>
        /// Returns a cell mask containing all cells that are considered
        /// "shocked" according to the given
        /// <paramref name="artificalViscosityLaw"/> (i.e., where
        /// <see cref="IArtificialViscosityLaw.IsShocked(int)"/> returns true)
        /// </summary>
        /// <param name="artificalViscosityLaw">
        /// Some artificial viscosity law
        /// </param>
        /// <param name="gridData">
        /// Information about the grid.
        /// </param>
        /// <returns>
        /// A cell mask containing all shocked cells
        /// </returns>
        public static CellMask GetShockedCellMask(this IArtificialViscosityLaw artificalViscosityLaw, IGridData gridData) {
            BitArray mask = new BitArray(gridData.iLogicalCells.NoOfLocalUpdatedCells);
            for (int cell = 0; cell < mask.Length; cell++) {
                mask[cell] = artificalViscosityLaw.IsShocked(cell);
            }
            return new CellMask(gridData, mask);
        }
    }
}

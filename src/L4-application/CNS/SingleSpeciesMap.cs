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
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;

namespace CNS {

    /// <summary>
    /// A species map that always assumes the same species. Useful for
    /// single-phase calculations.
    /// </summary>
    public class SingleSpeciesMap : ISpeciesMap {

        /// <summary>
        /// The equation of state that applies everywhere.
        /// </summary>
        private Material material;

        /// <summary>
        /// Initializes a mock map that always returns the given equation of
        /// state.
        /// </summary>
        /// <param name="gridData">
        /// Information about the grid.
        /// </param>
        /// <param name="material">
        /// The material/fluid that applies everywhere.
        /// </param>
        public SingleSpeciesMap(IGridData gridData, Material material) {
            this.GridData = gridData;
            this.material = material;
        }

        #region ISpeciesMap Members

        /// <summary>
        /// Information about the grid
        /// </summary>
        public IGridData GridData {
            get;
            private set;
        }

        /// <summary>
        /// Returns the material/fluid specified in the constructor.
        /// </summary>
        /// <param name="levelSetValue">Not used.</param>
        /// <returns>
        /// <see cref="SingleSpeciesMap.SingleSpeciesMap"/>
        /// </returns>
        public Material GetMaterial(double levelSetValue) {
            return material;
        }

        /// <summary>
        /// Always returns the full domain
        /// </summary>
        public SubGrid SubGrid {
            get {
                return new SubGrid(CellMask.GetFullMask(GridData));
            }
        }

        #endregion
    }
}

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

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// Flags giving additional information on cells.
    /// </summary>
    [Flags]
    public enum CellInfo {

        /// <summary>
        /// initial value
        /// </summary>
        Undefined = 0,

        /// <summary>
        /// marks the region of the info in which the cell type index (index
        /// into <see cref="GridCommons.RefElements"/>) is encoded
        /// </summary>
        RefElementIndex_Mask = 0x7, // binary: 0111

        /// <summary>
        /// marks affine-linear cells
        /// </summary>
        CellIsAffineLinear = 0x8, // binary: 1000

        /// <summary>
        /// Marks cells which are aggregated from smaller parts.
        /// </summary>
        IsAggregate = 0x10, // binary: 10000

        /// <summary>
        /// All flags on, implementing this ensures that the enum 
        /// is compiled using 32 bits.
        /// </summary>
        AllOn = unchecked((int)~0)
    }
}

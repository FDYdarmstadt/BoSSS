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

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Provides basic information about a grid.
    /// This interface serves as a kind of filter, providing only the information
    /// that is necessary for database-related operations.
    /// </summary>
    public interface IGridInfo : IDatabaseEntityInfo<IGridInfo>, IEquatable<IGridInfo> {

        /// <summary>
        /// The number of grid cells.
        /// </summary>
        int NumberOfCells {
            get;
        }

        /// <summary>
        /// The simplex dimension in the sense of measure-theory.
        /// </summary>
        int SpatialDimension {
            get;
        }
    }
}

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
        /// a string to store some user-information about the grid;
        /// </summary>
        string Description {
            get;
            set;
        }

        /// <summary>
        /// The simplex dimension in the sense of measure-theory.
        /// </summary>
        int SpatialDimension {
            get;
        }

        /// <summary>
        /// collection of all (parallel) data vectors used in the grid
        /// </summary>
        IReadOnlyCollection<Guid> AllDataVectorIDs {
            get;
        }

        /// <summary>
        /// This is a mapping from each used <em>EdgeTag</em> to a string that
        /// provides a name and additional information about the EdgeTag. The
        /// intention for this member is to provide both, a name (e.g.
        /// 'Left wall') for different regions of the boundary as well as
        /// boundary condition type info (e.g. 'inlet' or 'wall' or 'outflow' ...).
        /// </summary>
        /// <remarks>
        /// The names have no impact on the application on this application
        /// layer (L2-layer of BoSSS). They may be used on a higher application
        /// layer; Usually, this member (as like mostly all other public
        /// variable of this class) should be initialized by grid generator
        /// programs.
        /// </remarks>
        IDictionary<byte, string> EdgeTagNames {
            get;
        }

    }
}

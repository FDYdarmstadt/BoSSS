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
    /// Information about a single time-step
    /// </summary>
    public interface ITimestepInfo : IDatabaseEntityInfo<ITimestepInfo> {

        /// <summary>
        /// The time-step index/number is represented by the first entry in the
        /// array. In order to handle inner iterations (and possibly inner
        /// iterations of inner iterations, and so on), the subsequent array
        /// elements can be used. E.g. time-step 12 with solver iteration 2323
        /// and sub-iteration 369 would be Indices == [12, 2323, 369].
        /// </summary>
        TimestepNumber TimeStepNumber {
            get;
        }

        /// <summary>
        /// The grid corresponding to this time-step.
        /// </summary>
        IGridInfo Grid {
            get;
        }

        /// <summary>
        /// The session this time-step belongs to.
        /// </summary>
        ISessionInfo Session {
            get;
        }

        /// <summary>
        /// Stores information which grid is assigned to this time-step without
        /// having to load the grid.
        /// </summary>
        Guid GridID {
            get;
        }

        /// <summary>
        /// Physical time in the simulation
        /// </summary>
        double PhysicalTime {
            get;
        }

        /// <summary>
        /// Unique identifier referring to the storage vector
        /// </summary>
        Guid StorageID {
            get;
        }

        /// <summary>
        /// Information about the DG fields in this time-step: DG basis, class
        /// type, identification, ...
        /// </summary>
        IEnumerable<DGField.FieldInitializer> FieldInitializers {
            get;
        }

        /// <summary>
        /// Contains information on the fields of this time-step.
        /// </summary>
        IEnumerable<DGField> Fields {
            get;
        }
    }
    
}

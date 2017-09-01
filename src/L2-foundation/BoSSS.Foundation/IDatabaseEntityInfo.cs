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
    /// Common interface for all entities that can be saved to a database.
    /// </summary>
    /// <typeparam name="T">
    /// The type of the database entity (required for type-safe copy
    /// operations)
    /// </typeparam>
    public interface IDatabaseEntityInfo<out T> where T : IDatabaseEntityInfo<T> {

        /// <summary>
        /// Unique identifier of the specific IDatabaseEntityInfo object.
        /// </summary>
        Guid ID {
            get;
        }

        /// <summary>
        /// The time when the represented entity has been created.
        /// </summary>
        DateTime CreationTime {
            get;
        }

        /// <summary>
        /// The time when this object has been written to disc.
        /// </summary>
        DateTime WriteTime {
            get;
            set;
        }

        /// <summary>
        /// The name of the entity.
        /// </summary>
        string Name {
            get;
            set;
        }

        /// <summary>
        /// The database where this info object is located.
        /// </summary>
        IDatabaseInfo Database {
            get;
        }

        /// <summary>
        /// Copies this IDatabaseObjectInfo object for storage in a different
        /// database.
        /// </summary>
        /// <param name="targetDatabase">The target database</param>
        /// <returns>
        /// A copy of this IDatabaseObjectInfo object with all the same
        /// information, except for the database field, which will be the one of
        /// the target database
        /// </returns>
        T CopyFor(IDatabaseInfo targetDatabase);

    }
}

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
    /// Relevant information of a database object
    /// </summary>
    public interface IDatabaseInfo : IEquatable<IDatabaseInfo> {

        /// <summary>
        /// Full path to the base directory of the database.
        /// </summary>
        string Path {
            get;
        }

        /// <summary>
        /// detects if some other path actually also points to this database
        /// </summary>
        bool PathMatch(string otherPath);
        
        /// <summary>
        /// Alternative paths to access the database, if <see cref="Path"/> is not present on a given machine.
        /// This allows to use the same control file or object on different machines, where the database is located in a different path.
        /// - 1st entry: path into the local file system
        /// - 2nd entry: optional machine name filter
        /// </summary>
        (string DbPath, string MachineFilter)[] AlternateDbPaths {
            get;
        }

        /// <summary>
        /// Provides functionality to copy/move/delete info objects stored in
        /// the database
        /// </summary>
        IDatabaseController Controller {
            get;
        }

        /// <summary>
        /// The sessions of this database.
        /// </summary>
        IList<ISessionInfo> Sessions {
            get;
        }

        /// <summary>
        /// The grids of this database.
        /// </summary>
        IList<IGridInfo> Grids {
            get;
        }

        /// <summary>
        /// Sessions sorted according to projects, see <see cref="ISessionInfo.ProjectName"/>.
        /// </summary>
        IDictionary<string, IEnumerable<ISessionInfo>> Projects {
            get;
        }
    }
}

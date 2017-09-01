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
    /// Null object for <see cref="IDatabaseInfo"/>
    /// </summary>
    public class NullDatabaseInfo : IDatabaseInfo {

        /// <summary>
        /// The one and only instance (reference equality)
        /// </summary>
        public static readonly NullDatabaseInfo Instance = new NullDatabaseInfo();

        /// <summary>
        /// Hides the public default constructor
        /// </summary>
        private NullDatabaseInfo() {
        }

        #region IDatabaseInfo Members

        /// <summary>
        /// Returns an empty path
        /// </summary>
        public string Path {
            get {
                return "";
            }
        }

        /// <summary>
        /// Returns <see cref="NullDatabaseController.Instance"/>
        /// </summary>
        public IDatabaseController Controller {
            get {
                return NullDatabaseController.Instance;
            }
        }

        /// <summary>
        /// An empty collection.
        /// </summary>
        public IEnumerable<ISessionInfo> Sessions {
            get {
                yield break;
            }
        }

        /// <summary>
        /// An empty collection.
        /// </summary>
        public IEnumerable<IGridInfo> Grids {
            get {
                yield break;
            }
        }

        /// <summary>
        /// An empty dictionary.
        /// </summary>
        public IDictionary<string, IEnumerable<ISessionInfo>> Projects {
            get {
                return (new Dictionary<string, IEnumerable<ISessionInfo>>());
            }
        }

        #endregion


        
    }
}

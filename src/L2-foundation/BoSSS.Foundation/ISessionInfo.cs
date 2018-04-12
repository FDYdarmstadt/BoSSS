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
    /// Provides basic information about a session.
    /// </summary>
    public interface ISessionInfo : IDatabaseEntityInfo<ISessionInfo>, IEquatable<ISessionInfo>, IDisposable {

        /// <summary>
        /// A description of the session. Can be more detailed than the name.
        /// </summary>
        string Description {
            get;
            set;
        }

        /// <summary>
        /// Associates this session with a project
        /// </summary>
        string ProjectName {
            get;
            set;
        }

        /// <summary>
        /// All the time-steps of this session.
        /// </summary>
        IList<ITimestepInfo> Timesteps {
            get;
        }

        /// <summary>
        /// The session ID this session has been restarted from.
        /// </summary>
        Guid RestartedFrom {
            get;
        }

        /// <summary>
        /// the git commit hash of the master HEAD
        /// </summary>
        string MasterGitCommit {
            get;
        }

        /// <summary>
        /// A collection of tags for this session.
        /// </summary>
        IEnumerable<string> Tags {
            get;
            set;
        }

        /// <summary>
        /// Names of compute nodes on which the session is running; Index: MPI rank index in the
        /// MPI_COMM_WORLD communicator;
        /// </summary>
        IList<string> ComputeNodeNames {
            get;
        }

        /// <summary>
        /// Returns all the grids used in this session
        /// </summary>
        IEnumerable<IGridInfo> GetGrids();

        /// <summary>
        /// Keys (DG degree, various settings like timestepping scheme, also the queries) and the corresponding values.
        /// </summary>
        IDictionary<string, object> KeysAndQueries {
            get;
        }
    }
}
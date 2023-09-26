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
using System.Linq;
using System.Text;

namespace BoSSS.Foundation.Grid {
    /// <summary>
    /// methods to compute/define a grid distribution
    /// </summary>
    public enum GridPartType {

        /// <summary>
        /// use <see cref="ParMETIS"/> (parallel partitioning)
        /// </summary>
        ParMETIS = 1,

        /// <summary>
        /// Predefined partition.
        /// </summary>
        Predefined = 2,

        /// <summary>
        /// Partitioning according to space-filling clusterHilbert curve considering Clusters, ...
        /// </summary>
        clusterHilbert = 3,

        /// <summary>
        /// Partitioning according to space-filling clusterHilbert curve, direct Costmapping ...
        /// </summary>
        Hilbert = 4,

        /// <summary>
        /// leave grid as it is; The first J/P cells will be on first processor, ...
        /// </summary>
        none = 0,

        /// <summary>
        /// Use <see cref="METIS"/> (serial partitioning)
        /// </summary>
        METIS = 5,

        /// <summary>
        /// Mostly for debugging:
        /// Using a partitioning stored in a field `MPIrank` in another session;
        /// the session ID is drawn from the partitioning options string, see <see cref="IGrid.Redistribute"/>
        /// </summary>
        OtherSession = 6
    }
}

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

using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Aggregation {
    partial class AggregationGridData {

        /// <summary>
        /// see <see cref="Parallelization"/>
        /// </summary>
        private Parallelization m_Parallel;

        
        /// <summary>
        /// see <see cref="Parallelization"/>
        /// </summary>
        public IParallelization iParallel {
            get {
                return m_Parallel;
            }
        }

        /// <summary>
        /// Contains information for MPI parallelization.
        /// </summary>
        public class Parallelization : IParallelization {

            /// <summary>
            /// pointer to owner object
            /// </summary>
            internal AggregationGridData m_owner;

            /// <summary>
            /// Conversion of global cell indices to local cell indices, 
            /// i.e. the inverse of <see cref="GlobalIndicesExternalCells"/>. 
            ///  - keys: global indices of external/ghost cells 
            ///  - values: local indices of external/ghost cells
            /// </summary>
            public Dictionary<long, int> Global2LocalIdx {
                get;
                internal set;
            }

            /// <summary>
            /// Conversion of global cell indices to local cell indices
            /// </summary>
            /// <param name="globalCellIndex"></param>
            /// <remarks>
            /// The conversion can only be done for cells which are stored on
            /// this MPI process;
            /// </remarks>
            public int GetLocalCellIndex(long globalCellIndex) {
                Partitioning part = m_owner.CellPartitioning;
                int i0 = part.i0;
                int J = part.LocalLength;
                int jLoc = (int)globalCellIndex - i0;
                if (jLoc >= 0 && jLoc < J) {
                    return jLoc;
                } else {
                    if (Global2LocalIdx.TryGetValue(globalCellIndex, out jLoc)) {
                        return jLoc;
                    } else {
                        return int.MinValue;
                    }
                }
            }

            /// <summary>
            /// Global indices of external cells (local indices j in the range
            /// <see cref="CellData.NoOfLocalUpdatedCells"/> &lt;= j &lt;
            /// <see cref="CellData.NoOfCells"/>); Note that there is an index
            /// offset, so the entry at index 0 is the global index of cell at
            /// local index <see cref="CellData.NoOfLocalUpdatedCells"/>;
            /// See also <see cref="GlobalIndicesExternalCells"/>;
            /// </summary>
            public long[] GlobalIndicesExternalCells {
                get;
                internal set;
            }

            /// <summary>
            /// local cell indices (only border cells) which must be send to
            /// other processors; for each processor, the communication list is
            /// stored in ascending order.<br/>
            /// </summary>
            /// <remarks>
            /// - 1st index: target processor; if the 'p'-th entry is null, there
            /// is no communication with processor 'p'
            /// - 2nd index: enumeration;
            /// - content: local cell index
            /// </remarks>
            public int[][] SendCommLists {
                get;
                internal set;
            }

            /// <summary>
            /// for each process, a local cell index at which items received
            /// from other processes should be inserted;
            /// - index: MPI process rank of process from which data is received;
            /// - content: local cell index, where the received eitems should be inserted.
            /// </summary>
            public int[] RcvCommListsInsertIndex {
                get;
                internal set;
            }

            /// <summary>
            /// for each process, the number of cells that are received from
            /// this process
            /// index: MPI process rank of process from which data is received;
            /// </summary>
            public int[] RcvCommListsNoOfItems {
                get;
                internal set;
            }

            /// <summary>
            /// List of processes (MPI ranks) which receive data from this process.
            /// index: enumeration.
            /// </summary>
            public int[] ProcessesToSendTo {
                get;
                internal set;
            }

            /// <summary>
            /// List of processes (MPI ranks) which send data to this process.
            /// </summary>
            public int[] ProcessesToReceiveFrom {
                get;
                internal set;
            }

        }
    }
}

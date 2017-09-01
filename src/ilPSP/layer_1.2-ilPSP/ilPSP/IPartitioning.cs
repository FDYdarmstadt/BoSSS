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

using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP {

    /// <summary>
    /// Defines a basis partitioning of an index over a range of MPI processors, 
    /// i.e. over the respective communicator <see cref="IPartitioning.MPI_Comm"/>.
    /// </summary>
    public interface IPartitioning {

        /// <summary>
        /// MPI Communicator on which this object lives on
        /// </summary>
        MPI_Comm MPI_Comm {
            get;
        }

        /// <summary>
        /// Number of MPI processors in the communicator <see cref="MPI_Comm"/>.
        /// </summary>
        int MpiSize {
            get;
        }

        /// <summary>
        /// MPI rank of actual processor within the communicator <see cref="MPI_Comm"/>.
        /// </summary>
        int MpiRank {
            get;
        }

        /// <summary>
        /// the first global index that is stored on the actual MPI process
        /// </summary>
        int i0 {
            get;
        }

        /// <summary>
        /// The first global index that is stored on the NEXT MPI process
        /// </summary>
        int iE {
            get;
        }

        /// <summary>
        /// Index offset for some processor
        /// </summary>
        /// <param name="proc">
        /// process index; can be equal to number of processes (i.e. highest process
        /// rank +1); In this case, the <see cref="TotalLength"/> is returned.
        /// </param>
        /// <returns>index of the first permutation entry stored by processor <paramref name="proc"/></returns>
        int GetI0Offest(int proc);

        /// <summary>
        /// returns the number of entries which are stored by
        /// this processor;
        /// </summary>
        int LocalLength {
            get;
        }

        /// <summary>
        /// Total length of the partition over all processes, i.e. sum of <see cref="LocalLength"/> over all MPI processors.
        /// </summary>
        int TotalLength {
            get;
        }

        /// <summary>
        /// True, if any of the properties of this interface resp. object can change in time
        /// (e.g. due to re-partitioning of a grid, dynamic load balancing, local refinement).
        /// </summary>
        bool IsMutable {
            get;
        }

        /// <summary>
        /// If this partition is mutable (see <see cref="IsMutable"/>), this method should 
        /// return a frozen version of this partitioning.
        /// </summary>
        IPartitioning GetImmutablePartition();


        /// <summary>
        /// Returns the process rank which stores the <paramref name="index"/>-th entry;
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        int FindProcess(int index);

        /// <summary>
        /// returns the process rank which stores the <paramref name="index"/>-th entry;
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        int FindProcess(long index);

        /// <summary>
        /// true, if some index <paramref name="i"/> is within the local range of this MPI Process,
        /// i.e <paramref name="i"/> is greater or equal to <see cref="i0"/> and smaller than <see cref="i0"/>+<see cref="LocalLength"/>;
        /// </summary>
        bool IsInLocalRange(int i);

        /// <summary>
        /// Returns the number of entries which are stored by
        /// processor <paramref name="proc"/>;
        /// </summary>
        /// <param name="proc"></param>
        /// <returns></returns>
        int GetLocalLength(int proc);
    }
}

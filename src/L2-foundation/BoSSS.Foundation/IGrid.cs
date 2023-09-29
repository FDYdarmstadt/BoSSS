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

using BoSSS.Foundation.IO;
using BoSSS.Platform.LinAlg;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// Common interface for all grids
    /// </summary>
    public interface IGrid : IGridInfo {

        /// <summary>
        /// Access to grid metrics
        /// </summary>
        IGridData iGridData {
            get;
        }

        /// <summary>
        /// Releases the object cached in <see cref="iGridData"/>; this is necessary if 
        /// the grid object was changed somehow.
        /// </summary>
        void InvalidateGridData();

        


        /// <summary>
        /// Driver method for grid redistribution; this includes 
        /// - computing a new partition 
        /// - application of the partition to this grid, i.e. invocation off <see cref="RedistributeGrid(int[])"/>
        /// </summary>
        void Redistribute(IDatabaseInfo iom, GridPartType method, string PartOptions);

        /// <summary>
        /// Low-level interface:
        /// redistributes this grid, i.e. sends cells to different processors
        /// </summary>
        /// <param name="part">
        /// - index: local cell index;
        /// - content: MPI processor rank for respective cell; 
        /// </param>
        void RedistributeGrid(int[] part);
               

        /// <summary>
        /// MPI process rank (within world communicator)
        /// </summary>
        int MyRank {
            get;
        }

        /// <summary>
        /// MPI world communicator size 
        /// </summary>
        int Size {
            get;
        }

        /// <summary>
        /// Gets the partition of cells over the MPI processes;
        /// </summary>
        Partitioning CellPartitioning {
            get;
        }

        /// <summary>
        /// %
        /// </summary>
        IGridSerializationHandler GridSerializationHandler {
            get;
        }

        /// <summary>
        /// list of transformations which describe how some edges should be transformed to other edges;
        /// </summary>
        /// <remarks>
        /// indices into this list are edge tags (<see cref="IGeometricalEdgeData.EdgeTags"/> minus <see cref="Grid.Classic.GridCommons.FIRST_PERIODIC_BC_TAG"/>)
        /// </remarks>
        IList<AffineTrafo> PeriodicTrafo {
            get;
        }
    }
}

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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution {

    
    /// <summary>
    /// Base-class to handle data persistence under any form of mesh update,
    /// this can be either MPI-repartitioning or adaptive mesh refinement.
    /// The purpose of this class is to 
    /// - backup/serialize objects before balancing and
    /// - restore/serialize objects after balancing.
    /// </summary>
    abstract public class GridUpdateData {

        /// <summary>
        /// Old Grid, before re-distribution.
        /// </summary>
        protected IGridData m_OldGrid;

        /// <summary>
        /// Number of locally updated cells in old grid, before load balancing. (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).
        /// </summary>
        protected int m_oldJ;

        /// <summary>
        /// Level-set tracker before redistribution.
        /// </summary>
        protected LevelSetTracker m_OldTracker;

        /// <summary>
        /// Number of locally updated cells in old grid, before load balancing. (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).
        /// </summary>
        protected int m_newJ = -1;

        /// <summary>
        /// New Grid, before re-distribution.
        /// </summary>
        internal IGridData m_NewGrid;

        /// <summary>
        /// New Level-Set tracker, after grid redistribution.
        /// </summary>
        protected LevelSetTracker m_NewTracker;

        /// <summary>
        /// Stores the backup data of DG Fields before re-distribution.
        /// - dictionary key: string reference to identify DG field
        /// - 1st value index: local cell index, old grid
        /// - 2nd value index: DG coordinate index within cell
        /// </summary>
        protected Dictionary<string, double[][]> m_oldDGFieldData = new Dictionary<string, double[][]>();


        /// <summary>
        /// Loads the DG coordinates after update of grid.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="Reference">
        /// Unique string reference under which data has been stored before grid-redistribution.
        /// </param>
        public abstract void RestoreDGField(DGField f, string Reference);
        
        /// <summary>
        /// Backup data for some DG field before update of grid.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="Reference">
        /// Unique string reference under which the re-distributed data can be accessed later, within <see cref="RestoreDGField(DGField, string)"/>.
        /// </param>
        abstract public void BackupField(DGField f, string Reference);

        /// <summary>
        /// See <see cref="BackupField(DGField, string)"/>.
        /// </summary>
        public void BackupField(DGField f) {
            BackupField(f, f.Identification);
        }

        /// <summary>
        /// See <see cref="RestoreDGField(DGField, string)"/>.
        /// </summary>
        public void RestoreDGField(DGField f) {
            RestoreDGField(f, f.Identification);
        }


        /// <summary>
        /// Set the new level-set-tracker after it is available.
        /// </summary>
        public void SetNewTracker(LevelSetTracker NewTracker) {
            if(NewTracker != null && !object.ReferenceEquals(NewTracker.GridDat, m_NewGrid))
                throw new ArgumentException();
            m_NewTracker = NewTracker;
        }


        /// <summary>
        /// Whether the grid is only re-distributed, or also adapted.
        /// - true: mesh adaptation is used; the grid/mesh may also be re-distributed.
        /// - false: only grid re-distribution.
        /// </summary>
        public bool GridAdaptation {
            get;
            protected set;
        }


    }
}

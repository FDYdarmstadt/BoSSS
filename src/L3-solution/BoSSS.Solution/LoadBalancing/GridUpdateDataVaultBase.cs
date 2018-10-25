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
    abstract public class GridUpdateDataVaultBase {

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
        public int SetNewTracker(LevelSetTracker NewTracker) {
            if(NewTracker != null && !object.ReferenceEquals(NewTracker.GridDat, m_NewGrid))
                throw new ArgumentException();
            m_NewTracker = NewTracker;
            return this.RestoreTracker();
        }


        /// <summary>
        /// Various parameter required for the restoration of the level-set tracker. 
        /// </summary>
        protected class LevelSetTrackerPrivateData {

            /// <summary>
            /// <see cref="LevelSetTracker.HistoryLength"/>
            /// </summary>
            public int HistoryLength;

            /// <summary>
            /// <see cref="LevelSetTracker.PopulatedHistoryLength"/>
            /// </summary>
            public int PopultatedHistoryLength;

            /// <summary>
            /// For each time level, <see cref="LevelSetTracker.LevelSetRegions.Version"/>.
            /// </summary>
            public int[] Versions;
        }

        /// <summary>
        /// the only instance of <see cref="LevelSetTrackerPrivateData"/>
        /// </summary>
        protected LevelSetTrackerPrivateData m_LsTrkPrivData;

        /// <summary>
        /// Backup reference for some level-set field.
        /// </summary>
        protected string GetLSbackupName(int iHistory, int iLevSet) {
            return "LevelSetTracker_LevelSet_#" + iLevSet + "__at_time_level_" + iHistory + "_somewordstomakeitauniquekey44312341";
        }

        /// <summary>
        /// Backup reference for some level-set region code vector.
        /// </summary>
        protected string GetLSregioncodeName(int iHistory) {
            return "LevelSetTracker_RegionCode_at_time_level_" + iHistory + "_somewordstomakeitauniquekey325176";
        }

        /// <summary>
        /// Saves the internal state of <see cref="m_OldTracker"/>.
        /// </summary>
        public void BackupTracker () {
            using(new FuncTrace()) {
                if(m_LsTrkPrivData != null)
                    throw new NotSupportedException("Can only be called once.");
                m_LsTrkPrivData = new LevelSetTrackerPrivateData();
                int NoOfLS = m_OldTracker.NoOfLevelSets;

                if(m_OldTracker != null) {
                    m_LsTrkPrivData.HistoryLength = m_OldTracker.HistoryLength;
                    m_LsTrkPrivData.PopultatedHistoryLength = m_OldTracker.PopulatedHistoryLength;
                    m_LsTrkPrivData.Versions = new int[m_LsTrkPrivData.PopultatedHistoryLength + 1];

                    for(int iH = 1; iH > -m_LsTrkPrivData.PopultatedHistoryLength; iH--) {
                        var TimeLevel = m_OldTracker.BackupTimeLevel(iH);
                        for(int iLs = 0; iLs < NoOfLS; iLs++) {
                            this.BackupField(TimeLevel.Item1[iLs], GetLSbackupName(iH, iLs));
                        }
                        this.BackupVector(TimeLevel.Item2, GetLSregioncodeName(iH));

                        m_LsTrkPrivData.Versions[1 - iH] = TimeLevel.Item3;
                    }

                }
            }
        }

        /// <summary>
        /// Restores the status of <see cref="m_NewTracker"/>.
        /// </summary>
        protected abstract int RestoreTracker();
        
        /// <summary>
        /// Whether the grid is only re-distributed, or also adapted.
        /// - true: mesh adaptation is used; the grid/mesh may also be re-distributed.
        /// - false: only grid re-distribution.
        /// </summary>
        public bool GridAdaptation {
            get;
            protected set;
        }

        /// <summary>
        /// Backup data for some vector before grid-redistribution.
        /// </summary>
        /// <param name="vec">
        /// Output, a vector which indices correlate with the logical cell indices,
        /// i.e. its length must be a multiple of the number of cells (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).    
        /// </param>
        /// <param name="Reference">
        /// Unique string reference under which the re-distributed data can be accessed later, within <see cref="RestoreDGField(DGField, string)"/>.
        /// </param>
        public void BackupVector<T, V>(V vec, string Reference)
            where V : IList<T> //
        {
            using(new FuncTrace()) {
                if(vec.Count != m_oldJ)
                    throw new ArgumentException();

                var SerializedData = new double[m_oldJ][];

                using(var ms = new MemoryStream()) {
                    var fmt = new BinaryFormatter();

                    for(int j = 0; j < m_oldJ; j++) {
                        ms.Position = 0;
                        fmt.Serialize(ms, vec[j]);

                        SerializedData[j] = CopyData(ms);
                    }
                }

                m_oldDGFieldData.Add(Reference, SerializedData);
            }
        }

        static double[] CopyData(MemoryStream ms) {
            long Pos = ms.Position;
            double[] R = new double[(Pos / sizeof(double)) + 1];
            byte[] buf = ms.GetBuffer();

            unsafe {
                fixed (void* _pbuf = buf, _bR = R) {
                    byte* pbuf = (byte*)_pbuf, pR = (byte*)_bR;

                    for(int i = 0; i < Pos; i++)
                        pR[i] = pbuf[i];
                }
            }
            return R;
        }

        /// <summary>
        /// Backup data for some floating-point vector before grid-redistribution.
        /// </summary>
        /// <param name="vec">
        /// Input, a vector which indices correlate with the logical cell indices,
        /// i.e. its length must be a multiple of the number of cells (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).        
        /// </param>
        /// <param name="Reference">
        /// Unique string reference under which the re-distributed data can be accessed later, within <see cref="RestoreDGField(DGField, string)"/>.
        /// </param>
        public void BackupVector<V>(V vec, string Reference) where V : IList<double> {
            int L = vec.Count;
            if(L % m_oldJ != 0)
                throw new ArgumentException();
            int BlockLen = L / m_oldJ;

            double[][] dataVec = new double[m_oldJ][];
            for(int j = 0; j < m_oldJ; j++) {
                double[] data_j = new double[BlockLen];
                for(int n = 0; n < BlockLen; n++)
                    data_j[n] = vec[n + BlockLen * j];
                dataVec[j] = data_j;
            }

            m_oldDGFieldData.Add(Reference, dataVec);
        }

        
        /// <summary>
        /// Backup data for some 16-bit-integer vector before grid-redistribution (<see cref="LevelSetTracker.LevelSetRegions.RegionsCode"/>).
        /// </summary>
        /// <param name="vec">
        /// Input, a vector which indices correlate with the logical cell indices,
        /// i.e. its length must be a multiple of the number of cells (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).        
        /// </param>
        /// <param name="Reference">
        /// Unique string reference under which the re-distributed data can be accessed later, within <see cref="RestoreDGField(DGField, string)"/>.
        /// </param>
        public void BackupVector(ushort[] vec, string Reference)  {
            int L = vec.Length;
            if(L % m_oldJ != 0)
                throw new ArgumentException();
            int BlockLen = L / m_oldJ;

            double[][] dataVec = new double[m_oldJ][];
            for(int j = 0; j < m_oldJ; j++) {
                double[] data_j = new double[BlockLen];
                for(int n = 0; n < BlockLen; n++)
                    data_j[n] = vec[n + BlockLen * j];
                dataVec[j] = data_j;
            }

            m_oldDGFieldData.Add(Reference, dataVec);
        }
        


    }
}

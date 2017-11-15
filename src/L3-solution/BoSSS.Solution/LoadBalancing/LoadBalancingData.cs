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
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;

namespace BoSSS.Solution {

    /// <summary>
    /// During dynamic load balancing <see cref="Application{T}.MpiRedistributeAndMeshAdapt(int, double)"/>, this is used to
    /// - backup/serialize objects before balancing and
    /// - restore/serialize objects after balancing.
    /// </summary>
    public class LoadBalancingData {

        /// <summary>
        /// ctor;
        /// </summary>
        /// <param name="oldGrid"></param>
        /// <param name="oldTracker"></param>
        internal LoadBalancingData(IGridData oldGrid, LevelSetTracker oldTracker) {
            if(oldTracker != null && !object.ReferenceEquals(oldTracker.GridDat, oldGrid))
                throw new ArgumentException();
            m_OldGrid = oldGrid;
            m_OldTracker = oldTracker;
            m_oldJ = m_OldGrid.CellPartitioning.LocalLength;
        }

        /// <summary>
        /// Old Grid, before re-distribution.
        /// </summary>
        IGridData m_OldGrid;

        /// <summary>
        /// Number of locally updated cells in old grid, before load balancing. (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).
        /// </summary>
        int m_oldJ;

        /// <summary>
        /// Level-set tracker before redistribution.
        /// </summary>
        LevelSetTracker m_OldTracker;

        /// <summary>
        /// Number of locally updated cells in old grid, before load balancing. (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).
        /// </summary>
        int m_newJ = -1;

        /// <summary>
        /// New Grid, before re-distribution.
        /// </summary>
        internal IGridData m_NewGrid;

        /// <summary>
        /// New Level-Set tracker, after grid redistribution.
        /// </summary>
        LevelSetTracker m_NewTracker;

        /// <summary>
        /// Stores the backup data of DG Fields before re-distribution.
        /// - dictionary key: string reference to identify DG field
        /// - 1st value index: local cell index, old grid
        /// - 2nd value index: DG coordinate index within cell
        /// </summary>
        Dictionary<string, double[][]> m_oldDGFieldData = new Dictionary<string, double[][]>();

        /// <summary>
        /// Stores the backup data of DG Fields after re-distribution (when no mesh adaptation is performed).
        /// - dictionary key: string reference to identify DG field
        /// - 1st value index: local cell index, new grid
        /// - 2nd value index: DG coordinate index within cell
        /// </summary>
        Dictionary<string, double[][]> m_newDGFieldData_OnlyRedist;

        /// <summary>
        /// Stores the backup data of DG Fields after mesh adaptation.
        /// - dictionary key: string reference to identify DG field
        /// - 1st value index: local cell index, new grid.
        /// - 2nd value index: enumeration of cells; for refined and conserved cells, this contains only one element; for cells which are coarsened, this corresponds to all cells of the coarsening cluster.
        /// - 3rd value index: DG coordinate index within cell.
        /// </summary>
        Dictionary<string, double[][][]> m_newDGFieldData_GridAdapta;

        /// <summary>
        /// Resorting permutation, points from old to new, i.e. New[<see cref="Resorting"/>[j]] = Old[j].
        /// </summary>
        Permutation Resorting;

        /// <summary>
        /// Resorting permutation, points from old to new, i.e. New[<see cref="Resorting"/>[j]] = Old[j].
        /// </summary>
        Permutation InverseResorting;

        /// <summary>
        /// Backup data for some DG field before grid-redistribution.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="Reference">
        /// Unique string reference under which the re-distributed data can be accessed later, within <see cref="RestoreDGField(DGField, string)"/>.
        /// </param>
        public void BackupField(DGField f, string Reference) {
            if(!object.ReferenceEquals(f.GridDat, m_OldGrid)) {
                throw new ArgumentException("DG field seems to be assigned to some other grid.");
            }

            Basis oldBasis = f.Basis;

            double[][] oldFieldsData = new double[m_oldJ][];
            m_oldDGFieldData.Add(Reference, oldFieldsData);

            for(int j = 0; j < m_oldJ; j++) {
                int Nj = oldBasis.GetLength(j);
                double[] data_j = new double[Nj];
                for(int n = 0; n < Nj; n++) {
                    data_j[n] = f.Coordinates[j, n];
                }

                oldFieldsData[j] = data_j;
            }
        }

        /// <summary>
        /// See <see cref="BackupField(DGField, string)"/>.
        /// </summary>
        public void BackupField(DGField f) {
            BackupField(f, f.Identification);
        }

        /// <summary>
        /// Loads the DG coordinates after grid-redistribution.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="Reference">
        /// Unique string reference under which data has been stored before grid-redistribution.
        /// </param>
        public void RestoreDGField(DGField f, string Reference) {
            int newJ = this.m_newJ;
            if(this.GridAdaptation) {
                

            } else {
                double[][] newFieldsData = m_newDGFieldData_OnlyRedist[Reference];
                Debug.Assert(newFieldsData.Length == newJ);
                Basis newBasis = f.Basis;
                if(!object.ReferenceEquals(f.GridDat, m_NewGrid)) {
                    throw new ArgumentException("DG field seems to be assigned to some other grid.");
                }

                for(int j = 0; j < newJ; j++) {
                    double[] data_j = newFieldsData[j];
                    int Nj = data_j.Length;
                    if(data_j.Length != newBasis.GetLength(j))
                        throw new ApplicationException(string.Format("Data length mismatch: cell {0} expecting {1} coordinates, got {2}.", j, newBasis.GetLength(j), Nj));

                    for(int n = 0; n < Nj; n++) {
                        f.Coordinates[j, n] = data_j[n];
                    }
                }
            }
        }

        /// <summary>
        /// See <see cref="RestoreDGField(DGField, string)"/>.
        /// </summary>
        public void RestoreDGField(DGField f) {
            RestoreDGField(f, f.Identification);
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
        /// Loads some data vector before grid-redistribution.
        /// </summary>
        /// <param name="vec">
        /// Output, a vector which indices correlate with the logical cell indices,
        /// i.e. its length must be a multiple of the number of cells (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).
        /// </param>
        /// <param name="Reference">
        /// Unique string reference under which data has been stored before grid-redistribution.
        /// </param>
        public void RestoreVector<T, V>(V vec, string Reference)
            where V : IList<T> //
        {
            using(new FuncTrace()) {
                if(vec.Count != m_newJ)
                    throw new ArgumentException();

                var SerializedData = m_newDGFieldData_OnlyRedist[Reference];
                Debug.Assert(SerializedData.Length == m_newJ);
                var fmt = new BinaryFormatter();

                for(int j = 0; j < m_newJ; j++) {
                    using(var ms = new MemoryStream(ConvertBuffer(SerializedData[j]))) {
                        vec[j] = (T)(fmt.Deserialize(ms));
                    }
                }
            }
        }

        static byte[] ConvertBuffer(double[] buf) {
            byte[] R = new byte[buf.Length * sizeof(double)];

            unsafe {
                fixed (void* _pbuf = buf, _bR = R) {
                    byte* pbuf = (byte*)_pbuf, pR = (byte*)_bR;

                    for(int i = 0; i < R.Length; i++)
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
        /// Loads some floating-point data vector after grid-redistribution.
        /// </summary>
        /// <param name="vec">
        /// Output, a vector which indices correlate with the logical cell indices,
        /// i.e. its length must be a multiple of the number of cells (<see cref="ILogicalCellData.NoOfLocalUpdatedCells"/>).
        /// </param>
        /// <param name="Reference">
        /// Unique string reference under which data has been stored before grid-redistribution.
        /// </param>
        public void RestoreVector<V>(V vec, string Reference) where V : IList<double> {
            int L = vec.Count;
            if(L % m_newJ != 0)
                throw new ArgumentException();
            int BlockLen = L / m_newJ;

            double[][] dataVec = m_newDGFieldData_OnlyRedist[Reference];
            Debug.Assert(dataVec.Length == m_newJ);
            for(int j = 0; j < m_newJ; j++) {
                double[] data_j = dataVec[j];
                for(int n = 0; n < BlockLen; n++)
                    vec[n + BlockLen * j] = data_j[n];
            }
        }

        /// <summary>
        /// Loads some floating-point data vector after grid-redistribution.
        /// </summary>
        /// <param name="OutputFunc">
        /// Output, which stores the output;
        /// </param>
        /// <param name="Reference">
        /// Unique string reference under which data has been stored before grid-redistribution.
        /// </param>
        public void RestoreVector(Action<int, double[]> OutputFunc, string Reference) {
            double[][] dataVec = m_newDGFieldData_OnlyRedist[Reference];
            Debug.Assert(dataVec.Length == m_newJ);
            for(int j = 0; j < m_newJ; j++) {
                double[] data_j = dataVec[j];
                OutputFunc(j, data_j);
            }
        }


        Dictionary<string, CcmData> m_CCMdata = new Dictionary<string, CcmData>();

        class CcmData {
            public string[] SpeciesNamesList;
            public int HMForder;
            public XQuadFactoryHelper.MomentFittingVariants HMFvariant;
        }


        /// <summary>
        /// Backup of a cut-cell metric.
        /// </summary>
        public void BackupCutCellMetrics(CutCellMetrics ccm, string Reference) {
            if(!object.ReferenceEquals(m_OldTracker, ccm.Tracker))
                throw new ArgumentException();

            CcmData da = new CcmData();

            SpeciesId[] speciesIdList = ccm.SpeciesList.ToArray();
            int NoOfSpc = speciesIdList.Length;
            da.SpeciesNamesList = speciesIdList.Select(spc => m_OldTracker.GetSpeciesName(spc)).ToArray();
            da.HMForder = ccm.HMForder;
            da.HMFvariant = ccm.HMFvariant;


            m_CCMdata.Add(Reference, da);

            for(int iSpc = 0; iSpc < speciesIdList.Length; iSpc++) {
                SpeciesId spc = speciesIdList[iSpc];
                string SpeciesName = m_OldTracker.GetSpeciesName(spc);

                double[] Avol = ccm.CutCellVolumes[spc].ExtractSubArrayShallow(new int[] { 0 }, new int[] { m_oldJ - 1 }).To1DArray();
                double[] Aint = ccm.InterfaceArea[spc].ExtractSubArrayShallow(new int[] { 0 }, new int[] { m_oldJ - 1 }).To1DArray();
                MultidimensionalArray Aedg = ccm.CutEdgeAreas[spc];

                Debug.Assert(Aedg.Dimension == 1 && Aedg.GetLength(0) == m_OldGrid.iGeomEdges.Count);

                Tuple<long, double>[][] AedgGathered = GatherEdge2Cell(m_OldGrid, iEdge => Aedg[iEdge]);
                Debug.Assert(AedgGathered.Length == m_oldJ);

#if DEBUG
                BitArray CheckTouch = new BitArray(Aedg.GetLength(0));
                ScatterCell2Edge(m_OldGrid, AedgGathered, delegate (int iEdge, double val) {
                    Debug.Assert(Aedg[iEdge] == val);
                    CheckTouch[iEdge] = true;
                });
                for(int i = 0; i < Aedg.GetLength(0); i++)
                    Debug.Assert(CheckTouch[i] == true);
#endif

                BackupVector(Avol, Reference + "::CutCellMetrics-Vol:" + SpeciesName);
                BackupVector(Aint, Reference + "::CutCellMetrics-Int:" + SpeciesName);
                BackupVector<Tuple<long, double>[], Tuple<long, double>[][]>(AedgGathered, Reference + "::CutCellMetrics-Edge:" + SpeciesName);

#if DEBUG
                Tuple<long, double>[][] check_AedgGathered = new Tuple<long, double>[m_oldJ][];

                int bkup_new_J = m_newJ;
                m_newJ = m_oldJ; // hack to make next line working
                RestoreVector<Tuple<long, double>[], Tuple<long, double>[][]>(check_AedgGathered, Reference + "::CutCellMetrics-Edge:" + SpeciesName);
                m_newJ = bkup_new_J;
                for(int j = 0; j < m_oldJ; j++) {
                    var Lst1 = AedgGathered[j];
                    var Lst2 = check_AedgGathered[j];
                    Debug.Assert(Lst1.Length == Lst2.Length);
                    int LL = Lst2.Length;
                    for(int l = 0; l < LL; l++) {
                        Debug.Assert(Lst1[l].Item1 == Lst2[l].Item1);
                        Debug.Assert(Lst1[l].Item2 == Lst2[l].Item2);
                    }
                }
#endif
            }
        }

        /// <summary>
        /// Reload cut-cell metrics.
        /// </summary>
        public void RestoreCutCellMetrics(out CutCellMetrics ccm, string Reference) {
            if(m_NewTracker == null)
                throw new NotSupportedException();

            CcmData da = m_CCMdata[Reference];

            SpeciesId[] speciesIdList = da.SpeciesNamesList.Select(spcName => m_NewTracker.GetSpeciesId(spcName)).ToArray();
            int NoOfSpc = speciesIdList.Length;

            Dictionary<SpeciesId, MultidimensionalArray> CcVol = new Dictionary<SpeciesId, MultidimensionalArray>();
            Dictionary<SpeciesId, MultidimensionalArray> CcInt = new Dictionary<SpeciesId, MultidimensionalArray>();
            Dictionary<SpeciesId, MultidimensionalArray> CcEdg = new Dictionary<SpeciesId, MultidimensionalArray>();

            int JE = m_newJ + m_NewGrid.iLogicalCells.NoOfExternalCells;
            double[] vec_cellMetrics = new double[JE * speciesIdList.Length * 2];

            // collect all per-cell-metrics in the same MultidimArry, for MPI-exchange
            // 1st index: cell
            // 2nd index: species
            // 3rd index: 0 for interface surface per cell, 1 for cut-cell-volume
            MultidimensionalArray cellMetrics = MultidimensionalArray.CreateWrapper(vec_cellMetrics, JE, speciesIdList.Length, 2);


            for(int iSpc = 0; iSpc < speciesIdList.Length; iSpc++) {
                SpeciesId spc = speciesIdList[iSpc];
                string SpeciesName = m_NewTracker.GetSpeciesName(spc);

                double[] Avol = new double[m_newJ];
                double[] Aint = new double[m_newJ];
                Tuple<long, double>[][] AedgGathered = new Tuple<long, double>[m_newJ][];

                RestoreVector(delegate (int j, double[] v) { cellMetrics[j, iSpc, 1] = v[0]; Debug.Assert(v.Length == 1); }, Reference + "::CutCellMetrics-Vol:" + SpeciesName);
                RestoreVector(delegate (int j, double[] v) { cellMetrics[j, iSpc, 0] = v[0]; Debug.Assert(v.Length == 1); }, Reference + "::CutCellMetrics-Int:" + SpeciesName);
                RestoreVector<Tuple<long, double>[], Tuple<long, double>[][]>(AedgGathered, Reference + "::CutCellMetrics-Edge:" + SpeciesName);

                CcVol.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1));
                CcInt.Add(spc, cellMetrics.ExtractSubArrayShallow(-1, iSpc, 1));

                MultidimensionalArray Aedg = MultidimensionalArray.Create(m_NewGrid.iGeomEdges.Count);
                ScatterCell2Edge(m_NewGrid, AedgGathered, delegate (int iEdge, double val) { Aedg[iEdge] = val; });

                CcEdg.Add(spc, Aedg);
            }

            vec_cellMetrics.MPIExchange(m_NewGrid);

            ccm = new CutCellMetrics(da.HMFvariant, da.HMForder, m_NewTracker,
                CcVol, CcInt, CcEdg);
        }



        /// <summary>
        /// Gathers data from edges to cells.
        /// </summary>
        /// <param name="InputFunc">
        /// Provides the input data, i.e. returns a data item for each edge index.
        /// </param>
        /// <param name="grd"></param>
        /// <returns></returns>
        static Tuple<long, T>[][] GatherEdge2Cell<T>(IGridData grd, Func<int, T> InputFunc) {
            int NoEdg = grd.iGeomEdges.Count;

            int[,] Edg2Cel = grd.iGeomEdges.CellIndices;
            int J = grd.iLogicalCells.NoOfLocalUpdatedCells;
            var FaceIdx = grd.iGeomEdges.FaceIndices;

            var Ret = new List<Tuple<long, T>>[J];

            for(int iEdge = 0; iEdge < NoEdg; iEdge++) {
                var Data = InputFunc(iEdge);

                int jCell0 = Edg2Cel[iEdge, 0];
                int jCell1 = Edg2Cel[iEdge, 1];
                int Face0 = FaceIdx[iEdge, 0];
                int Face1 = FaceIdx[iEdge, 1];
                if(Ret[jCell0] == null)
                    Ret[jCell0] = new List<Tuple<long, T>>();
                if(jCell1 >= 0 && jCell1 < J && Ret[jCell1] == null)
                    Ret[jCell1] = new List<Tuple<long, T>>();

                if(jCell1 >= 0) {
                    //
                    long Ident0 = grd.iLogicalCells.GetGlobalID(jCell1);
                    long Ident1 = grd.iLogicalCells.GetGlobalID(jCell0);
                    Debug.Assert(Ident0 >= 0);
                    Debug.Assert(Ident1 >= 0);
                    Ret[jCell0].Add(new Tuple<long, T>(Ident0, Data));
                    if(jCell1 < J)
                        Ret[jCell1].Add(new Tuple<long, T>(Ident1, Data));
                } else {
                    Face1 = 1234;
                    long Ident0 = (0xffffffff & (long)Face0 | (0xffffffff & (long)Face1) << 32);
                    Ident0 *= -1;
                    Debug.Assert(Ident0 < 0);
                    Ret[jCell0].Add(new Tuple<long, T>(Ident0, Data));
                }

                Debug.Assert(jCell0 < J);
            }

            // converting lists to arrays sucks...
            var _Ret = new Tuple<long, T>[J][];
            for(int j = 0; j < J; j++) {
                if(Ret[j] == null)
                    _Ret[j] = new Tuple<long, T>[0];
                else
                    _Ret[j] = Ret[j].ToArray();
            }

            return _Ret;
        }

        /// <summary>
        /// Scatters data from cells back to edges, i.e. the inverse to <see cref="GatherEdge2Cell{T}(Func{int, T})"/>
        /// </summary>
        /// <param name="Data">
        /// Function which provides the local edge index and the associated data.
        /// </param>
        /// <param name="OutputFunc"></param>
        /// <param name="grd"></param>
        static void ScatterCell2Edge<T>(IGridData grd, Tuple<long, T>[][] Data, Action<int, T> OutputFunc) {
            int NoEdg = grd.iGeomEdges.Count;

            int[,] Edg2Cel = grd.iGeomEdges.CellIndices;
            int J = grd.iLogicalCells.NoOfLocalUpdatedCells;
            int[][] C2E = grd.iLogicalCells.Cells2Edges;
            var FaceIdx = grd.iGeomEdges.FaceIndices;

            if(Data.Length != J)
                throw new ArgumentException();

            for(int j = 0; j < J; j++) { // loop over cells
                long Gid = grd.iLogicalCells.GetGlobalID(j);
                int[] Edges = C2E[j];

                foreach(var tt in Data[j]) { // for each data item in cell 'j'...
                    long Ident = tt.Item1;
                    T ed = tt.Item2;

                    bool boundaryEdge;
                    int IdentFace;
                    if(Ident < 0) {
                        boundaryEdge = true;
                        Ident *= -1;
                        IdentFace = (int)(0xffffffff & (Ident));

                        int MagicCheck = (int)(0xffffffff & (Ident >> 32));
                        Debug.Assert(MagicCheck == 1234); // set above

                        Debug.Assert(IdentFace >= 0);
                        Debug.Assert(IdentFace < byte.MaxValue);
                    } else {
                        boundaryEdge = false;
                        IdentFace = int.MinValue;
                    }

                    //int IdentFace0 = (int)(0xffffffff & (Ident));
                    //int IdentFace1 = (int)(0xffffffff & (Ident >> 32));

                    bool bFound = false;
                    for(int s = 0; s < Edges.Length; s++) {
                        int iEdge = Math.Abs(Edges[s]) - 1;
                        int Face0 = FaceIdx[iEdge, 0];
                        int jCell0 = Edg2Cel[iEdge, 0];
                        int jCell1 = Edg2Cel[iEdge, 1];

                        if(boundaryEdge) {
                            // searching Boundary Edge 
                            if(jCell1 < 0) {
                                Debug.Assert(j == jCell0);

                                if(IdentFace == Face0) {
                                    Debug.Assert(grd.iLogicalEdges.CellIndices[iEdge, 1] < 0);
                                    OutputFunc(iEdge, ed);
                                    bFound = true;
                                    continue;
                                }
                            }
                        } else {
                            if(jCell1 >= 0) {
                                long NeighGlobalId;
                                Debug.Assert((j == jCell0) || (j == jCell1));
                                Debug.Assert((j == jCell0) != (j == jCell1));
                                if(j == jCell0)
                                    NeighGlobalId = grd.iLogicalCells.GetGlobalID(jCell1);
                                else
                                    NeighGlobalId = grd.iLogicalCells.GetGlobalID(jCell0);
                                Debug.Assert(NeighGlobalId != grd.iLogicalCells.GetGlobalID(j));

                                // searching Inner Edge
                                if(Ident == NeighGlobalId) {
                                    OutputFunc(iEdge, ed);
                                    bFound = true;
                                    continue;
                                }
                            }
                        }
                    }
                    Debug.Assert(bFound);
                }
            }

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
        /// Apply the resorting (only mesh redistribution between MPI processors, no mesh adaptation)
        /// </summary>
        public void Resort(Permutation r, IGridData NewGrid) {
            using(new FuncTrace()) {
                this.GridAdaptation = false;
                this.Resorting = r;
                this.InverseResorting = r.Invert();
                m_OldGrid = null;
                m_OldTracker = null;
                m_NewGrid = NewGrid;
                m_newJ = m_NewGrid.iLogicalCells.NoOfLocalUpdatedCells;


                // fix the sequence in which we serialize/de-serialize fields
                string[] FieldNames = this.m_oldDGFieldData.Keys.ToArray();
                int NoFields = FieldNames.Length;

                // format data for serialization
                double[][] oldFieldsData = new double[m_oldJ][];
                {
                    double[][][] oldFields = FieldNames.Select(rstring => this.m_oldDGFieldData[rstring]).ToArray();

                    for(int j = 0; j < m_oldJ; j++) {
                        int Ntot = 0;
                        for(int iF = 0; iF < NoFields; iF++)
                            Ntot += oldFields[iF][j].Length;

                        double[] data_j = new double[Ntot + NoFields];
                        int cnt = 0;
                        for(int iF = 0; iF < NoFields; iF++) {
                            double[] data_iF_j = oldFields[iF][j];
                            int Nj = data_iF_j.Length;
                            data_j[cnt] = Nj;
                            cnt++;
                            for(int n = 0; n < Nj; n++) {
                                data_j[cnt] = data_iF_j[n];
                                cnt++;
                            }
                        }
                        oldFieldsData[j] = data_j;
                    }
                    oldFields = null;
                    m_oldDGFieldData = null;
                }

                // redistribute data
                double[][] newFieldsData = new double[m_newJ][];
                {
                    Resorting.ApplyToVector(oldFieldsData, newFieldsData, m_NewGrid.CellPartitioning);
                    oldFieldsData = null;
                    Debug.Assert(newFieldsData.Length == m_newJ);
                }

                // re-format re-distributed data
                m_newDGFieldData_OnlyRedist = new Dictionary<string, double[][]>();
                {
                    double[][][] newFields = new double[FieldNames.Length][][];
                    for(int iF = 0; iF < NoFields; iF++) {
                        newFields[iF] = new double[m_newJ][];
                    }
                    for(int j = 0; j < m_newJ; j++) {
                        double[] data_j = newFieldsData[j];

                        int cnt = 0;
                        for(int iF = 0; iF < NoFields; iF++) {
                            int Nj = (int)data_j[cnt];
                            cnt++;
                            double[] data_iF_j = new double[Nj];
                            for(int n = 0; n < Nj; n++) {
                                data_iF_j[n] = data_j[cnt];
                                cnt++;
                            }
                            newFields[iF][j] = data_iF_j;
                        }
                        Debug.Assert(cnt == data_j.Length);
                    }

                    for(int iF = 0; iF < NoFields; iF++) {
                        m_newDGFieldData_OnlyRedist.Add(FieldNames[iF], newFields[iF]);
                    }
                }
            }
        }

        /// <summary>
        /// Whether the grid is only re-distributed, or also adapted.
        /// - true: mesh adaptation is used; the grid/mesh may also be re-distributed.
        /// - false: only grid re-distribution.
        /// </summary>
        public bool GridAdaptation {
            get;
            private set;
        }


        /// <summary>
        /// Apply the resorting (including mesh adaptation).
        /// </summary>
        public void Resort(GridCorrelation r, GridData NewGrid) {
            using(new FuncTrace()) {
                this.GridAdaptation = false;
                this.Resorting = null;
                this.InverseResorting = null;

                m_OldGrid = null;
                m_OldTracker = null;
                m_NewGrid = NewGrid;
                m_newJ = m_NewGrid.iLogicalCells.NoOfLocalUpdatedCells;
                
                // fix the sequence in which we serialize/de-serialize fields
                string[] FieldNames = this.m_oldDGFieldData.Keys.ToArray();
                int NoFields = FieldNames.Length;

                // format data for serialization
                double[][] oldFieldsData = new double[m_oldJ][]; // 1st: cell; 2nd enum
                {

                    double[][][] oldFields = FieldNames.Select(rstring => this.m_oldDGFieldData[rstring]).ToArray(); // 1st: field; 2nd: cell; 3rd: DG coord

                    for(int j = 0; j < m_oldJ; j++) {
                        int Ntot = 0;
                        for(int iF = 0; iF < NoFields; iF++)
                            Ntot += oldFields[iF][j].Length;

                        // pack the DG coords for each field into one array of the form
                        // { Np_f1, DGcoords_f1, Np_f2, DGcoords_f2, .... };

                        double[] data_j = new double[Ntot + NoFields];
                        int cnt = 0;
                        for(int iF = 0; iF < NoFields; iF++) {
                            double[] data_iF_j = oldFields[iF][j];
                            int Nj = data_iF_j.Length;
                            data_j[cnt] = Nj;
                            cnt++;
                            for(int n = 0; n < Nj; n++) {
                                data_j[cnt] = data_iF_j[n];
                                cnt++;
                            }
                        }
                        oldFieldsData[j] = data_j;
                    }
                    oldFields = null;
                    m_oldDGFieldData = null;
                }

                // redistribute data
                double[][][] newFieldsData = new double[m_newJ][][];
                {
                    //Resorting.ApplyToVector(oldFieldsData, newFieldsData, m_NewGrid.CellPartitioning);
                    r.ApplyToVector(oldFieldsData, newFieldsData, m_NewGrid.CellPartitioning);
                    oldFieldsData = null;
                    Debug.Assert(newFieldsData.Length == m_newJ);
                }

                // re-format re-distributed data
                m_newDGFieldData_GridAdapta = new Dictionary<string, double[][][]>();
                {
                    double[][][][] newFields = new double[FieldNames.Length][][][]; // indices: [field, cell idx, enum over coarsening, DG mode]
                    for(int iF = 0; iF < NoFields; iF++) {
                        newFields[iF] = new double[m_newJ][][];
                    }
                    for(int j = 0; j < m_newJ; j++) {
                        double[][] data_j = newFieldsData[j];
                        int L = data_j.Length; // cell cluster size (equal 1 for refined or conserved cells)

                        for(int l = 0; l < L; l++) {
                            double[] data_jl = data_j[l];

                            int cnt = 0;
                            for(int iF = 0; iF < NoFields; iF++) {
                                int Nj = (int)data_jl[cnt];
                                cnt++;
                                double[] data_iF_jl = new double[Nj];
                                for(int n = 0; n < Nj; n++) {
                                    data_iF_jl[n] = data_jl[cnt];
                                    cnt++;
                                }
                                newFields[iF][j][l] = data_iF_jl;
                            }
                            Debug.Assert(cnt == data_jl.Length);
                        }
                    }

                    for(int iF = 0; iF < NoFields; iF++) {
                        m_newDGFieldData_GridAdapta.Add(FieldNames[iF], newFields[iF]);
                    }
                }
            }

        }



        /// <summary>
        /// Permutation matrix from an old to a new partitioning. 
        /// </summary>
        /// <param name="RowPart">Row partitioning, i.e. new data partitioning.</param>
        /// <param name="ColPart">Column partitioning, i.e. old data partitioning.</param>
        /// <param name="tau">
        /// Permutation from new to old Indices.
        /// </param>
        /// <returns></returns>
        static BlockMsrMatrix GetRowPermutationMatrix(IBlockPartitioning RowPart, IBlockPartitioning ColPart, Permutation tau) {
            BlockMsrMatrix P = new BlockMsrMatrix(RowPart, ColPart);
            //if (RowPart.LocalNoOfBlocks != tau.LocalLength)
            //    throw new ArgumentException();
            if(RowPart.TotalNoOfBlocks != tau.TotalLength)
                throw new ArgumentException();
            if(!RowPart.AllBlockSizesEqual)
                throw new NotSupportedException("unable to perform redistribution for variable size blocking (unable to compute offsets for variable size blocking).");
            if(!ColPart.AllBlockSizesEqual)
                throw new NotSupportedException("unable to perform redistribution for variable size blocking (unable to compute offsets for variable size blocking).");
            if(RowPart.TotalLength != ColPart.TotalLength)
                throw new ArgumentException();

            int IBlock = RowPart.GetBlockLen(RowPart.FirstBlock);
            if(ColPart.GetBlockLen(ColPart.FirstBlock) != IBlock)
                throw new ArgumentException();

            int J = RowPart.LocalNoOfBlocks;
            long FB = RowPart.FirstBlock;
            long[] LocalBlockIdxS = J.ForLoop(i => i + FB);
            long[] TargetBlockIdxS = new long[LocalBlockIdxS.Length];
            tau.EvaluatePermutation(LocalBlockIdxS, TargetBlockIdxS);

            MultidimensionalArray TempBlock = MultidimensionalArray.Create(IBlock, IBlock);

            for(int jSrc_Loc = 0; jSrc_Loc < J; jSrc_Loc++) { // loop over cells resp. local block-indices
                int jSrcGlob = jSrc_Loc + RowPart.FirstBlock; // block-row index
                int jDstGlob = (int)TargetBlockIdxS[jSrc_Loc]; // block-column index

                Debug.Assert(RowPart.IsLocalBlock(jSrcGlob));
                int i0 = RowPart.GetBlockI0(jSrcGlob);
                int BT = RowPart.GetBlockType(jSrcGlob);
                int[] _i0 = RowPart.GetSubblk_i0(BT);
                int[] Len = RowPart.GetSubblkLen(BT);
                Debug.Assert(IBlock == RowPart.GetBlockLen(jSrcGlob));


                int j0 = IBlock * jDstGlob; // this would not work for variable size blocking
#if DEBUG
                if(ColPart.IsLocalBlock(jDstGlob)) {
                    // column block corresponds to some cell 
                    Debug.Assert(IBlock == ColPart.GetBlockLen(jDstGlob));
                    Debug.Assert(j0 == ColPart.GetBlockI0(jDstGlob));
                    int CBT = ColPart.GetBlockType(jDstGlob);
                    Debug.Assert(ArrayTools.AreEqual(_i0, ColPart.GetSubblk_i0(CBT)));
                    Debug.Assert(ArrayTools.AreEqual(Len, ColPart.GetSubblkLen(CBT)));
                }
#endif
                Debug.Assert(_i0.Length == Len.Length);
                int K = _i0.Length;

                for(int i = 0; i < IBlock; i++)
                    TempBlock[i, i] = 0.0;

                for(int k = 0; k < K; k++) {
                    int A = _i0[k];
                    int E = Len[k] + A;
                    for(int i = A; i < E; i++)
                        TempBlock[i, i] = 1;
                }

                P.AccBlock(i0, j0, 1.0, TempBlock);
            }

            return P;
        }

        /// <summary>
        /// Old Matrices which will be redistributed.
        /// </summary>
        Dictionary<string, Tuple<int[], int[], BlockMsrMatrix>> m_Matrices = new Dictionary<string, Tuple<int[], int[], BlockMsrMatrix>>();


        static public int[] GetDGBasisHash(IList<Basis> BS) {
            int[] pS = new int[BS.Count];
            for(int i = 0; i < pS.Length; i++) {
                if(BS[i] is XDGBasis) {
                    pS[i] = -BS[i].Degree;
                } else {
                    pS[i] = -BS[i].Degree - 1;
                }
            }
            return pS;
        }

        static BlockPartitioning CloneBlockPartitioning(IBlockPartitioning part) {
            int J = part.LocalNoOfBlocks;
            int j0 = part.FirstBlock;
            if(!part.AllBlockSizesEqual)
                throw new NotSupportedException();

            int Sz = part.GetBlockI0(j0 + 1) - part.GetBlockI0(j0);
            int[] BlockType = new int[J];

            int maxBT = 0;
            for(int j = 0; j < J; j++) {
                int blockType = part.GetBlockType(j + j0);
                BlockType[j] = blockType;
                maxBT = Math.Max(maxBT, blockType);
            }
            int[][] _Subblk_i0 = new int[maxBT + 1][];
            int[][] _SubblkLen = new int[maxBT + 1][];
            for(int b = 0; b <= maxBT; b++) {
                _Subblk_i0[b] = part.GetSubblk_i0(b);
                _SubblkLen[b] = part.GetSubblkLen(b);
            }

            return new BlockPartitioning(
                part.LocalLength,
                Sz,
                _Subblk_i0, _SubblkLen,
                BlockType,
                part.MPI_Comm);
        }

        Dictionary<int[], BlockPartitioning> m_OldBlockings = new Dictionary<int[], BlockPartitioning>(new FuncEqualityComparer<int[]>((A, B) => ArrayTools.AreEqual(A, B)));

        /// <summary>
        /// Backup of some operator matrix before grid-redistribution.
        /// </summary>
        /// <param name="M">
        /// Matrix which should be saved.
        /// </param>
        /// <param name="Reference">
        /// Unique string reference under which data can be referenced after grid-redistribution.
        /// </param>
        /// <param name="RowMapping">
        /// Coordinate mapping to correlate the matrix rows with cells of the computational grid.
        /// </param>
        /// <param name="ColMapping">
        ///  Coordinate mapping to correlate the matrix columns with cells of the computational grid.
        /// </param>
        public void BackupMatrix(BlockMsrMatrix M, string Reference, UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping) {
            CheckMatrix(M, RowMapping, ColMapping, m_oldJ);

            int[] rowHash = GetDGBasisHash(RowMapping.BasisS);
            int[] colHash = GetDGBasisHash(ColMapping.BasisS);

            if(!m_OldBlockings.ContainsKey(rowHash))
                m_OldBlockings.Add(rowHash, CloneBlockPartitioning(M._RowPartitioning));
            if(!m_OldBlockings.ContainsKey(colHash))
                m_OldBlockings.Add(colHash, CloneBlockPartitioning(M._ColPartitioning));


            BlockMsrMatrix Msave = M.RecyclePermute(m_OldBlockings[rowHash], m_OldBlockings[colHash], new int[0, 2], new int[0, 2]);

            m_Matrices.Add(Reference, new Tuple<int[], int[], BlockMsrMatrix>(rowHash, colHash, Msave));
        }

        private static void CheckMatrix(BlockMsrMatrix M, UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping, int J) {
            if(M._RowPartitioning.LocalNoOfBlocks != J)
                throw new ArgumentException("Number of column blocks must be equal to number of cells.");
            if(M._ColPartitioning.LocalNoOfBlocks != J)
                throw new ArgumentException("Number of column blocks must be equal to number of cells.");
            if(!M._RowPartitioning.EqualsPartition(RowMapping))
                throw new ArgumentException("Mapping must match partitioning.");
            if(!M._ColPartitioning.EqualsPartition(ColMapping))
                throw new ArgumentException("Mapping must match partitioning.");
            if(!M._RowPartitioning.AllBlockSizesEqual)
                throw new NotSupportedException("Only mappings with constant frame blocks are supported.");
            if(!M._ColPartitioning.AllBlockSizesEqual)
                throw new NotSupportedException("Only mappings with constant frame blocks are supported.");
        }


        Dictionary<int[], BlockMsrMatrix> m_RowPermutationMatrices = new Dictionary<int[], BlockMsrMatrix>(new FuncEqualityComparer<int[]>((a, b) => ArrayTools.AreEqual(a, b)));
        Dictionary<int[], BlockMsrMatrix> m_ColPermutationMatrices = new Dictionary<int[], BlockMsrMatrix>(new FuncEqualityComparer<int[]>((a, b) => ArrayTools.AreEqual(a, b)));


        /// <summary>
        /// Reload of some operator after before grid-redistribution.
        /// </summary>
        /// <param name="Mtx_new">
        /// Output, matrix which should be restored.
        /// </param>
        /// <param name="Reference">
        /// Unique string reference under which data has been stored before grid-redistribution.
        /// </param>
        /// <param name="RowMapping">
        /// Coordinate mapping to correlate the matrix rows with cells of the computational grid.
        /// </param>
        /// <param name="ColMapping">
        ///  Coordinate mapping to correlate the matrix columns with cells of the computational grid.
        /// </param>
        public void RestoreMatrix(BlockMsrMatrix Mtx_new, string Reference, UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping) {
            int J = m_NewGrid.iLogicalCells.NoOfLocalUpdatedCells;
            CheckMatrix(Mtx_new, RowMapping, ColMapping, J);

            BlockMsrMatrix Mtx_old = m_Matrices[Reference].Item3;
            int[] RowHash = GetDGBasisHash(RowMapping.BasisS);
            int[] ColHash = GetDGBasisHash(ColMapping.BasisS);

            if(!ArrayTools.AreEqual(RowHash, m_Matrices[Reference].Item1))
                throw new ApplicationException();
            if(!ArrayTools.AreEqual(ColHash, m_Matrices[Reference].Item2))
                throw new ApplicationException();

            BlockMsrMatrix P_Row = GetRowPermutationMatrix(Mtx_new, Mtx_old, RowHash);
            BlockMsrMatrix P_Col = GetColPermutationMatrix(Mtx_new, Mtx_old, ColHash);

            Debug.Assert(Mtx_new._RowPartitioning.LocalNoOfBlocks == m_newJ);
            Debug.Assert(Mtx_new._ColPartitioning.LocalNoOfBlocks == m_newJ);
            Debug.Assert(Mtx_old._RowPartitioning.LocalNoOfBlocks == m_oldJ);
            Debug.Assert(Mtx_old._ColPartitioning.LocalNoOfBlocks == m_oldJ);

            Debug.Assert(P_Row._RowPartitioning.LocalNoOfBlocks == m_newJ);
            Debug.Assert(P_Row._ColPartitioning.LocalNoOfBlocks == m_oldJ);
            Debug.Assert(P_Col._RowPartitioning.LocalNoOfBlocks == m_oldJ);
            Debug.Assert(P_Col._ColPartitioning.LocalNoOfBlocks == m_newJ);

            BlockMsrMatrix.Multiply(Mtx_new, BlockMsrMatrix.Multiply(P_Row, Mtx_old), P_Col);
        }

        private BlockMsrMatrix GetRowPermutationMatrix(BlockMsrMatrix Mtx_new, BlockMsrMatrix Mtx_old, int[] DGHash) {
            BlockMsrMatrix P_Row;
            if(!m_RowPermutationMatrices.TryGetValue(DGHash, out P_Row)) {
                Debug.Assert(Mtx_new._RowPartitioning.LocalNoOfBlocks == m_newJ);
                Debug.Assert(Mtx_old._RowPartitioning.LocalNoOfBlocks == m_oldJ);

                P_Row = GetRowPermutationMatrix(Mtx_new._RowPartitioning, Mtx_old._RowPartitioning, InverseResorting);
                m_RowPermutationMatrices.Add(DGHash, P_Row);
            }

            return P_Row;
        }

        private BlockMsrMatrix GetColPermutationMatrix(BlockMsrMatrix Mtx_new, BlockMsrMatrix Mtx_old, int[] DGHash) {
            BlockMsrMatrix P_Col, P_Row;
            if(!m_ColPermutationMatrices.TryGetValue(DGHash, out P_Col)) {
                if(!m_RowPermutationMatrices.TryGetValue(DGHash, out P_Row)) {
                    Debug.Assert(Mtx_new._ColPartitioning.LocalNoOfBlocks == m_newJ);
                    Debug.Assert(Mtx_old._ColPartitioning.LocalNoOfBlocks == m_oldJ);
                    P_Row = GetRowPermutationMatrix(Mtx_new._ColPartitioning, Mtx_old._ColPartitioning, InverseResorting);
                    m_RowPermutationMatrices.Add(DGHash, P_Row);
                }

                P_Col = P_Row.Transpose();
                m_ColPermutationMatrices.Add(DGHash, P_Col);
            }

            return P_Col;
        }



    }
}

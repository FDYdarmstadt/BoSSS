﻿/* =======================================================================
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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Comm;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation.Caching;
using System.Collections;

namespace BoSSS.Foundation.Grid.Aggregation {
    public partial class AggregationGridData : IGridData {

        /// <summary>
        /// The grid from which this was coarsened.
        /// </summary>
        public IGridData ParentGrid {
            get;
            private set;
        }

        /// <summary>
        /// MPI process rank (within world communicator)
        /// </summary>
        public int MpiRank {
            get {
                return this.CellPartitioning.MpiRank;
            }
        }

        /// <summary>
        /// MPI world communicator size 
        /// </summary>
        public int MpiSize {
            get {
                return this.CellPartitioning.MpiSize;
            }
        }

        /// <summary>
        /// The ancestor grid, from which the aggregation sequence was derived.
        /// </summary>
        public GridData AncestorGrid {
            get {
                if(ParentGrid is AggregationGridData agd) {
                    GridData Ancestor =agd.AncestorGrid;
                    Debug.Assert(this.iGeomCells.Count == Ancestor.iGeomCells.Count);
                    return Ancestor;
                } else {
                    return (GridData) ParentGrid; // end of recursion
                }
            }
        }

        /// <summary>
        /// Multi-grid level index.
        /// </summary>
        public int MgLevel {
            get {
                if(ParentGrid is AggregationGridData) {
                    return ((AggregationGridData)ParentGrid).MgLevel + 1;
                } else {
                    return 0;
                }
            }
        }

        /// <summary>
        /// The spatial dimension, induced from the parent grid (<see cref="ParentGrid"/>).
        /// </summary>
        public int SpatialDimension {
            get {
                return ParentGrid.SpatialDimension;
            }
        }


        /// <summary>
        /// Gets the partitioning of cells over the MPI processes.
        /// </summary>
        public Partitioning CellPartitioning {
            get;
            private set;
        }

        /// <summary>
        /// The global ID for each cell
        /// </summary>
        public Permutation CurrentGlobalIdPermutation {
            get {
                return aggregationGrid.CurrentGlobalIdPermutation();
            }
        }

        AggregationGrid aggregationGrid;

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="aggregationGrid">
        /// Aggregation grid.
        /// </param>
        /// <param name="AggregationCells">
        /// Coarse cells which build up the fine cells.
        /// - 1st index: coarse (i.e. this) grid cell index
        /// - 2nd index: enumeration
        /// - content: local cell index into the parent grid <see cref="AggregationGrid.ParentGrid"/>.
        /// </param>
        public AggregationGridData(AggregationGrid aggregationGrid, int[][] AggregationCells) {

            this.aggregationGrid = aggregationGrid;
            IGridData pGrid = aggregationGrid.ParentGrid.iGridData;
            InitializeGridData(pGrid, AggregationCells);
        }

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="pGrid">
        /// Parent grid.
        /// </param>
        /// <param name="AggregationCells">
        /// Coarse cells which build up the fine cells.
        /// - 1st index: coarse (i.e. this) grid cell index
        /// - 2nd index: enumeration
        /// - content: local cell index into the parent grid <paramref name="pGrid"/>.
        /// </param>
        public AggregationGridData(IGridData pGrid, int[][] AggregationCells) {
            InitializeGridData(pGrid, AggregationCells);
            aggregationGrid = null;
        }

        void InitializeGridData(IGridData pGrid, int[][] AggregationCells) {
            ParentGrid = pGrid;

            int JlocFine = pGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int JElocFine = pGrid.iLogicalCells.Count;

            m_GeomCellData = new GeomCellData() { m_Owner = this };
            m_LogicalCellData = new LogicalCellData() { m_Owner = this };
            m_GeomEdgeData = new GeomEdgeData() { m_Owner = this };
            m_LogEdgeData = new LogEdgeData(this);
            m_VertexData = new VertexData();
            m_Parallel = new Parallelization() { m_owner = this };

            CellPartitioning = new Partitioning(AggregationCells.Length, pGrid.CellPartitioning.MPI_Comm);

            long j0Coarse = CellPartitioning.i0;

            BuildNeighborship(AggregationCells);
            DefineCellParts();
            CollectEdges();
            m_GeomEdgeData.CollectGeomEdges2logCells();

            m_ChefBasis = new _BasisData(this);
        }

        private void CollectEdges() {

            int[] F2C = this.jCellFine2jCellCoarse;
            
            int JlocCoarse = iLogicalCells.NoOfLocalUpdatedCells;
            int JElocCoarse = iLogicalCells.Count;

            // collect edges
            // =========================


            // temporary cells-to-edge map
            List<int>[] tmpC2E = new List<int>[JElocCoarse];
            for (int jC = 0; jC < JElocCoarse; jC++) {
                tmpC2E[jC] = new List<int>();
            }

            // temporary edges map
            var tmpEdges = new List<EdgeTmp>();


            int[,] E2Cfine = ParentGrid.iLogicalEdges.CellIndices;
            int NoOfEdgesFine = ParentGrid.iLogicalEdges.Count;
            int[][] FineLogicalToGeom = ParentGrid.iLogicalEdges.EdgeToParts;
            for (int iEdgeFine = 0; iEdgeFine < NoOfEdgesFine; iEdgeFine++) { // loop over all logical edges in the fine grid.
                int jCellFine1 = E2Cfine[iEdgeFine, 0];
                int jCellFine2 = E2Cfine[iEdgeFine, 1];

                int jCellCoarse1 = F2C[jCellFine1];
                int jCellCoarse2 = jCellFine2 >= 0 ? F2C[jCellFine2] : -7544890;

                if (jCellCoarse1 != jCellCoarse2) {


                    if (jCellCoarse2 >= 0) {
                        // +++++++++++++
                        // internal edge
                        // +++++++++++++

                        // does the edge exists already?
                        EdgeTmp edgTmp = null;
                        foreach (int i in tmpC2E[jCellCoarse1]) {
                            int iEdge = Math.Abs(i) - 1;
                            EdgeTmp _edgTmp = tmpEdges[iEdge];
                            if ((_edgTmp.jCell1 == jCellCoarse1 && _edgTmp.jCell2 == jCellCoarse2)
                                || (_edgTmp.jCell2 == jCellCoarse1 && _edgTmp.jCell1 == jCellCoarse2)) {
                                edgTmp = _edgTmp;
                                break;
                            }
                        }

                        // allocate edge structure 
                        if (edgTmp == null) {
                            edgTmp = new EdgeTmp();
                            edgTmp.jCell1 = jCellCoarse1;
                            edgTmp.jCell2 = jCellCoarse2;

                            Debug.Assert(edgTmp.jCell1 != edgTmp.jCell2);

                            tmpEdges.Add(edgTmp);
                            tmpC2E[jCellCoarse1].Add(+tmpEdges.Count); // rem.: edge index shifted by 1, 
                            tmpC2E[jCellCoarse2].Add(-tmpEdges.Count); // sign denotes in resp. out - cell.
                        }

                        // add edge parts
                        int[] FineL2G = null; // logical-to-geometrical edge translation on the fine grid.
                        if (FineLogicalToGeom != null) {
                            FineL2G = FineLogicalToGeom[iEdgeFine];
                        }
                        if (FineL2G == null) {
                            FineL2G = new int[] { iEdgeFine };
                        }
                        edgTmp.LogicalToGeometricalMap.AddRange(FineL2G);

                    } else {
                        // +++++++++++++
                        // boundary edge
                        // +++++++++++++


                        // does the edge exists already?
                        EdgeTmp edgTmp = null;
                        foreach (int i in tmpC2E[jCellCoarse1]) {
                            int iEdge = Math.Abs(i) - 1;
                            EdgeTmp _edgTmp = tmpEdges[iEdge];
                            if (_edgTmp.jCell1 == jCellCoarse1 && _edgTmp.jCell2 < 0) {
                                edgTmp = _edgTmp;
                                break;
                            }
                        }

                        // allocate edge structure 
                        if (edgTmp == null) {
                            edgTmp = new EdgeTmp();
                            edgTmp.jCell1 = jCellCoarse1;
                            edgTmp.jCell2 = jCellCoarse2;
                            tmpEdges.Add(edgTmp);
                            tmpC2E[jCellCoarse1].Add(+tmpEdges.Count); // rem.: edge index shifted by 1, 
                        }

                        // add edge parts
                        int[] FineL2G = null; // logical-to-geometrical edge translation on the fine grid.
                        if (FineLogicalToGeom != null) {
                            FineL2G = FineLogicalToGeom[iEdgeFine];
                        }
                        if (FineL2G == null) {
                            FineL2G = new int[] { iEdgeFine };
                        }
                        edgTmp.LogicalToGeometricalMap.AddRange(FineL2G);
                    }
                }
            }

            // convert temporary data structures to final ones
            // ===============================================

            int NoOfEdg = tmpEdges.Count;
            m_LogEdgeData.CellIndices = new int[NoOfEdg,2];
            m_LogicalCellData.Cells2Edges = new int[JElocCoarse][];
            m_LogEdgeData.EdgeToParts = new int[NoOfEdg][];

            int[,] E2C =  m_LogEdgeData.CellIndices;
            int[][] C2E = m_LogicalCellData.Cells2Edges;
            int[][] L2G = m_LogEdgeData.EdgeToParts;

            for(int j = 0; j < JElocCoarse; j++) {
                C2E[j] = tmpC2E[j].ToArray();
            }

            
            for(int e = 0; e < NoOfEdg; e++) {
                EdgeTmp eTmp = tmpEdges[e];
                L2G[e] = eTmp.LogicalToGeometricalMap.ToArray();
                Array.Sort<int>(L2G[e]);
                E2C[e, 0] = eTmp.jCell1;
                E2C[e, 1] = eTmp.jCell2;
            }


#if DEBUG
            {
                // test if the logical edges match the geometrical ones 

                bool[] geomMarker = new bool[iGeomEdges.Count];
                for (int iEdge = 0; iEdge < iLogicalEdges.Count; iEdge++) {
                    foreach (int iGE in iLogicalEdges.EdgeToParts[iEdge]) {
                        Debug.Assert(geomMarker[iGE] == false, "geometrical edge referenced by two different logical edges");
                        geomMarker[iGE] = true;
                    }
                }
                //for (int iGE = 0; iGE < iGeomEdges.Count; iGE++) {
                //    Debug.Assert(geomMarker[iGE] == true, "unreferenced geometrical edge");
                //}
            }
#endif
        }

        
        /// <summary>
        /// helper data structure
        /// </summary>
        class EdgeTmp {
            public int jCell1;
            public int jCell2;
            public List<int> LogicalToGeometricalMap = new List<int>();
        }

        /// <summary>
        /// Mapping from parent grid cell index to cell index within this aggregation grid.
        /// - index: local cell index within this aggregation grid.
        /// - content: local cell index within the parent grid (<see cref="ParentGrid"/>)
        /// </summary>
        public int[][] jCellCoarse2jCellFine {
            get;
            private set;
        }
        

        /// <summary>
        /// Mapping from parent grid cell index to cell index within this aggregation grid.
        /// - index: local cell index within the parent grid (<see cref="ParentGrid"/>)
        /// - content: local cell index within this aggregation grid.
        /// </summary>
        public int[] jCellFine2jCellCoarse {
            get;
            private set;
        }

        public void MergeWithPartentGrid(AggregationGridData parent) {
            var grid = this;
            Debug.Assert(parent.CellPartitioning.LocalLength >= grid.CellPartitioning.LocalLength, "target is smaller then this grid. Do you messed up things?");
        
            var JCoarse = grid.iLogicalCells.Count;
            var Jparentparent = parent.jCellFine2jCellCoarse.Length;

            var mergedC2F = new int[JCoarse][];
            var mergedF2C = new int[Jparentparent];

            int checkcounter = 0;
            for (int jC = 0; jC < JCoarse; jC++) {
                var tmp = new List<int>();
                foreach (int jP in grid.jCellCoarse2jCellFine[jC]) {
                    int[] jPPs = parent.jCellCoarse2jCellFine[jP];
                    tmp.AddRange(jPPs);
                }
                mergedC2F[jC] = tmp.ToArray();
                checkcounter += tmp.Count();
            }
            Debug.Assert(checkcounter== Jparentparent);

            for (int jF = 0; jF < Jparentparent; jF++) {
                int jC = grid.jCellFine2jCellCoarse[parent.jCellFine2jCellCoarse[jF]];
                mergedF2C[jF] = jC;
            }
#if DEBUG
            // Test 4 surjective mapping
            // test the coarse-to-fine map
            bool[] testMarker = new bool[Jparentparent];
            int[][] C2F = mergedC2F;
            Debug.Assert(C2F.Length == JCoarse);
            for (int jC = 0; jC < JCoarse; jC++) {
                foreach (int jF in C2F[jC]) {
                    Debug.Assert(testMarker[jF] == false, $"cell {jF} already appears in coarse grid.");
                    testMarker[jF] = true;
                }
            }
            for (int jF = 0; jF < Jparentparent; jF++) {
                Debug.Assert(testMarker[jF] == true, $"cell {jF} of fine grid was not agglomerated");
            }

            // test the fine-to-coarse map
            int[] F2C = mergedF2C;
            Debug.Assert(F2C.Length == Jparentparent);
            for (int jF = 0; jF < Jparentparent; jF++) {
                Debug.Assert(C2F[F2C[jF]].Contains(jF), $"mapping is not invertable!{jF} could not been mapped back");
            }
#endif
            grid.ParentGrid = parent.ParentGrid;
            grid.jCellCoarse2jCellFine = mergedC2F;
            grid.jCellFine2jCellCoarse = mergedF2C;
        }

        /// <summary>
        /// Sets up the decomposition of the aggregate cell into elementary parts, which can be mapped to reference elements.
        /// </summary>
        private void DefineCellParts() {
            IGridData pGrid = ParentGrid;
            int J = this.iLogicalCells.NoOfLocalUpdatedCells;
            int JE = this.iLogicalCells.Count;

            m_LogicalCellData.AggregateCellToParts = new int[JE][];
            m_GeomCellData.GeomCell2LogicalCell = new int[pGrid.iGeomCells.Count];

            int[][] AggPartC = m_LogicalCellData.AggregateCellToParts;
            int[][] AggPartF = pGrid.iLogicalCells.AggregateCellToParts; // aggregate-to-geometric on parent grid
            int[] Geom2Agg = m_GeomCellData.GeomCell2LogicalCell;
#if DEBUG
            ArrayTools.SetAll(Geom2Agg, -1);
#endif

            List<int> tmp = new List<int>();
            for (int j = 0; j < JE; j++) { // loop over coarse/this cells
                tmp.Clear();

                foreach(int jFine in jCellCoarse2jCellFine[j]) { // loop over all fine cells aggregated in 'j' ...
                    if(AggPartF == null) {
                        tmp.Add(jFine);
                    } else {
                        int[] FinePart = AggPartF[jFine];
                        if(FinePart == null) {
                            tmp.Add(jFine);
                        } else {
                            tmp.AddRange(FinePart);
                        }
                    }
                }

                AggPartC[j] = tmp.ToArray();
                foreach(int jGeom in AggPartC[j]) {
#if DEBUG
                    Debug.Assert(Geom2Agg[jGeom] < 0);
#endif
                    Geom2Agg[jGeom] = j;
                }
            }

#if DEBUG
            Debug.Assert(Geom2Agg.Where(i => i < 0).Count() == 0);
#endif
        }

        private void BuildNeighborship(int[][] AggregationCells) {
            IGridData pGrid = ParentGrid;
            long j0Coarse = CellPartitioning.i0;
            int JlocCoarse = CellPartitioning.LocalLength;
            int JElocFine = pGrid.iLogicalCells.Count;
            int JlocFine = pGrid.iLogicalCells.NoOfLocalUpdatedCells;

            // compute fine-to-coarse mapping
            // ==============================

            long[] Fine2CoarseGlobal = new long[JElocFine]; // index: local cell index on fine grid; maps to _global_ cell index on coarse grid.
            {
                ArrayTools.SetAll(Fine2CoarseGlobal, -111L);
                for (int jCellCoarse = 0; jCellCoarse < JlocCoarse; jCellCoarse++) {
                    foreach (int jCellFine in AggregationCells[jCellCoarse]) {
                        if (jCellFine < 0 || jCellFine >= JlocFine)
                            throw new ArgumentOutOfRangeException();
                        Fine2CoarseGlobal[jCellFine] = jCellCoarse + j0Coarse;
                    }
                }


                Fine2CoarseGlobal.MPIExchange<long[], long>(pGrid);
            }

            // define external coarse cells
            // ============================
            int JElocCoarse;
            {
                HashSet<long> tmpCoarseExternal = new HashSet<long>();
                for (int jF = JlocFine; jF < JElocFine; jF++) {
                    tmpCoarseExternal.Add(Fine2CoarseGlobal[jF]);
                }

                long[] ExtGlbIdx = tmpCoarseExternal.ToArray();
                Array.Sort(ExtGlbIdx); // since the send lists are ascending, the external cells also should be ascending;
                //                        since the MPI rank only increases with the global cell index, this sort-operation also guarantees 
                //                        that the external cells are sorted according to MPI rank.
                m_Parallel.GlobalIndicesExternalCells = ExtGlbIdx;
                JElocCoarse = JlocCoarse + ExtGlbIdx.Length;

                // global-to-local index mapping:
                Dictionary<long, int> GlobalIdx2Local = new Dictionary<long, int>();
                for(int jC = JlocCoarse; jC < JElocCoarse; jC++) {
                    GlobalIdx2Local.Add(ExtGlbIdx[jC - JlocCoarse], jC);
                }
                m_Parallel.Global2LocalIdx = GlobalIdx2Local;
            }

            // fine-to-coarse, coarse-to-fine mapping in local indices
            // ========================================================
            {
                // fine-to-coarse
                this.jCellFine2jCellCoarse = new int[JElocFine];
                var F2C = jCellFine2jCellCoarse;

                for (int jF = 0; jF < JElocFine; jF++) { // loop over fine cells...
                    int jCoarse;
                    if (jF < JlocFine) {
                        jCoarse = checked((int)(Fine2CoarseGlobal[jF] - j0Coarse));
                        Debug.Assert(jCoarse >= 0);
                        Debug.Assert(jCoarse < JlocCoarse);
                    } else {
                        jCoarse = m_Parallel.Global2LocalIdx[Fine2CoarseGlobal[jF]];
                        Debug.Assert(jCoarse >= JlocCoarse);
                        Debug.Assert(jCoarse < JElocCoarse);
                    }

                    F2C[jF] = jCoarse;
                }

                // coarse-to-fine mapping
                this.jCellCoarse2jCellFine = new int[JElocCoarse][];
                var C2F = jCellCoarse2jCellFine;
                var tmpC2F = new List<int>[JElocCoarse];
                for (int jC = 0; jC < JElocCoarse; jC++)
                    tmpC2F[jC] = new List<int>();
                for (int jFine = 0; jFine < JElocFine; jFine++) {
                    int jCoarse = F2C[jFine];
                    Debug.Assert(!tmpC2F[jCoarse].Contains(jFine));
                    tmpC2F[jCoarse].Add(jFine);
                }
                for (int jC = 0; jC < JElocCoarse; jC++)
                    C2F[jC] = tmpC2F[jC].ToArray();
            }

            // define neighborship
            // ===================
            {
                int[][] FineClNeigh = pGrid.iLogicalCells.CellNeighbours;

                HashSet<int>[] tmpClNeig = new HashSet<int>[JlocCoarse];
                Dictionary<long, int> ExtCells_GlobalIdx2Local = m_Parallel.Global2LocalIdx;


                for (int jCellCoarse = 0; jCellCoarse < JlocCoarse; jCellCoarse++) { // loop over all coarse grid cells...
                    foreach (int jCellFine in AggregationCells[jCellCoarse]) {
                        int[] FineNeighs = FineClNeigh[jCellFine];
                        foreach (int jCellFineNeigh in FineNeighs) { // loop over all neighbor cells in the fine grid...
                            long jCellCoarseNeighGlob = Fine2CoarseGlobal[jCellFineNeigh]; // map index of fine grid neighbor to coarse cell index

                            // convert global cell index to local 
                            int jCellCoarseNeighLoc;
                            if (jCellCoarseNeighGlob >= j0Coarse && jCellCoarseNeighGlob <= (j0Coarse + JlocCoarse)) {
                                // Neighbor is a locally updated cell
                                jCellCoarseNeighLoc = checked((int)(jCellCoarseNeighGlob - j0Coarse));

                            } else {
                                // Neighbor is an external cell

                                jCellCoarseNeighLoc = ExtCells_GlobalIdx2Local[jCellCoarseNeighGlob];
                            }

                            // add neighbor cell
                            if (jCellCoarse != jCellCoarseNeighLoc) {
                                if (tmpClNeig[jCellCoarse] == null)
                                    tmpClNeig[jCellCoarse] = new HashSet<int>();
                                tmpClNeig[jCellCoarse].Add(jCellCoarseNeighLoc);
                            }
                        }
                    }
                }

                m_LogicalCellData.CellNeighbours = new int[JlocCoarse][];
                var CN = m_LogicalCellData.CellNeighbours;
                for (int jCellCoarse = 0; jCellCoarse < JlocCoarse; jCellCoarse++) {
                    var hs = tmpClNeig[jCellCoarse];
                    Debug.Assert(hs == null || hs.Contains(jCellCoarse) == false);
                    if (hs != null) {
                        int[] xx = new int[hs.Count];
                        CN[jCellCoarse] = xx;
                        int i = 0;
                        foreach(int x in hs) {
                            xx[i] = x;
                            i++;
                        }
                    } else {
                        CN[jCellCoarse] = new int[0];
                    }
                    Debug.Assert(CN[jCellCoarse].Contains(jCellCoarse) == false);
                }
            }

            // MPI send lists
            // ==============
            int mpiSize = CellPartitioning.MpiSize;
            {
                

                m_Parallel.ProcessesToSendTo = ParentGrid.iParallel.ProcessesToSendTo.CloneAs();
                m_Parallel.ProcessesToReceiveFrom = ParentGrid.iParallel.ProcessesToReceiveFrom.CloneAs();

                int[] F2C = this.jCellFine2jCellCoarse;

                m_Parallel.SendCommLists = new int[mpiSize][];
                var tmpSendList = new HashSet<int>();
                for(int rnk = 0; rnk < mpiSize; rnk++) {
                    int[] ParrentSendList = ParentGrid.iParallel.SendCommLists[rnk];
                    if (ParrentSendList != null && ParrentSendList.Length > 0) {
                        tmpSendList.Clear();

                        foreach(int jFine in ParrentSendList) {
                            Debug.Assert(jFine >= 0);
                            Debug.Assert(jFine < ParentGrid.iLogicalCells.NoOfLocalUpdatedCells);

                            int jCoarse = F2C[jFine];
                            Debug.Assert(jCoarse >= 0);
                            Debug.Assert(jCoarse < this.iLogicalCells.NoOfLocalUpdatedCells);

                            tmpSendList.Add(jCoarse);
                        }

                        m_Parallel.SendCommLists[rnk] = tmpSendList.ToArray();
                        Array.Sort(m_Parallel.SendCommLists[rnk]);
                    }
                }
            }

            // MPI receive lists
            // =================
            {
                //ParentGrid.iParallel.ProcessesToReceiveFrom
                var GlobalIdx = this.m_Parallel.GlobalIndicesExternalCells;

                this.m_Parallel.RcvCommListsNoOfItems = new int[mpiSize];
                this.m_Parallel.RcvCommListsInsertIndex = new int[mpiSize];
                int[] NoOfItems = this.m_Parallel.RcvCommListsNoOfItems;
                int[] InsertIdx = this.m_Parallel.RcvCommListsInsertIndex;
                ArrayTools.SetAll(InsertIdx, -1);

#if DEBUG
                HashSet<int> ProcCheck = new HashSet<int>();
#endif

                for(int jC = JlocCoarse; jC < JElocCoarse; jC++) {
                    int jGlob = (int) GlobalIdx[jC - JlocCoarse];
                    int iProc = this.CellPartitioning.FindProcess(jGlob);
#if DEBUG
                    ProcCheck.Add(iProc);
#endif
                    if(InsertIdx[iProc] < 0) {
                        InsertIdx[iProc] = jC;
                    } else {
                        Debug.Assert(jC > InsertIdx[iProc]);
                    }
                    NoOfItems[iProc]++;
                }

#if DEBUG
                Debug.Assert(this.iParallel.ProcessesToReceiveFrom.SetEquals(ProcCheck));
#endif
            }
            
            // MPI check
            // =========
#if DEBUG
            {
                long[] TestData = new long[JElocCoarse];
                for (int jC = 0; jC < JlocCoarse; jC++) {
                    long GlobalIdx;
                    GlobalIdx = jC + j0Coarse;
                    TestData[jC] = GlobalIdx;
                }

                TestData.MPIExchange<long[], long>(this);

                for (int jC = 0; jC < JElocCoarse; jC++) {
                    long GlobalIdx;
                    if(jC < JlocCoarse)
                        GlobalIdx = jC + j0Coarse;
                    else 
                        GlobalIdx = (int)(m_Parallel.GlobalIndicesExternalCells[jC - JlocCoarse]);

                    Debug.Assert(TestData[jC] == GlobalIdx);
                }
            }
#endif
        }

        public void TransformLocal2Global(MultidimensionalArray LocalVerticesIn, MultidimensionalArray GlobalVerticesOut, int jCell) {
            ParentGrid.TransformLocal2Global(LocalVerticesIn, GlobalVerticesOut, jCell);
        }

        public void TransformLocal2Global(MultidimensionalArray LocalVerticesIn, int j0, int Len, MultidimensionalArray GlobalVerticesOut, int OutArrayOffset) {
            ParentGrid.TransformLocal2Global(LocalVerticesIn, j0, Len, GlobalVerticesOut, OutArrayOffset);
        }

        public void TransformLocal2Global(MultidimensionalArray NS, int j0, int Len, MultidimensionalArray Nodesglob) {
            ParentGrid.TransformLocal2Global(NS, j0, Len, Nodesglob);
        }

        public void TransformGlobal2Local(MultidimensionalArray GlobalVerticesIn, MultidimensionalArray LocalVerticesOut, int j0, int Len, int OutArrayOffset) {
            ParentGrid.TransformGlobal2Local(GlobalVerticesIn, LocalVerticesOut, j0, Len, OutArrayOffset);
        }

        public void TransformGlobal2Local(MultidimensionalArray GlobalVerticesIn, MultidimensionalArray LocalVerticesOut, int jCell, bool[] NewtonConvergence) {
            ParentGrid.TransformGlobal2Local(GlobalVerticesIn, LocalVerticesOut, jCell, NewtonConvergence);
        }

        
        public CacheLogicImplBy_CNs InverseJacobian {
            get {
                return ParentGrid.InverseJacobian;
            }
        }

        public CacheLogicImplBy_CNs JacobianDeterminat {
            get {
                return ParentGrid.JacobianDeterminat;
            }
        }

        public CacheLogic_CNs GlobalNodes {
            get {
                return ParentGrid.GlobalNodes;
            }
        }

        public CacheLogicImplBy_CNs Jacobian {
            get {
                return ParentGrid.Jacobian;
            }
        }

        public CacheLogicImplBy_CNs AdjungateJacobian {
            get {
                return ParentGrid.AdjungateJacobian;
            }
        }

        /// <summary>
        /// Identification of the grid in the BoSSS database, 
        /// equal to the <see cref="BoSSS.Foundation.IO.IDatabaseEntityInfo{T}.ID"/>.
        /// </summary>
        public Guid GridID {
            get {
                return Grid.ID;
            }
        }

        /// <summary>
        /// The grid for which information is provided
        /// </summary>
        public IGrid Grid {
            get {
                if (aggregationGrid == null) {
                    throw new Exception("Grid was never set.");
                }
                return aggregationGrid;
            }
        }

        /// <summary>
        /// returns the remaining multigrid hierarchy
        /// </summary>
        public AggregationGridData[] MultigridSequence {
            get {

                var baseMG = this.AncestorGrid.MultigridSequence;
                for(int iLevel = 0; iLevel < baseMG.Length; iLevel++) {
                    if(baseMG[iLevel].MgLevel == this.MgLevel + 1)
                        return baseMG.Skip(iLevel).ToArray();
                }

                return null;
            }
        }

        public IDictionary<byte, string> EdgeTagNames {
            get {
                return ParentGrid.EdgeTagNames;
            }
        }

        /// <summary>
        /// Clears (lots of) internal references for this object, to make sure that any attempt to use it leads to an exception.
        /// </summary>
        public void Invalidate() {
            m_IsAlive = false;
            this.m_GeomCellData = null;
            this.m_LogicalCellData = null;

            this.m_GeomEdgeData = null;
            this.m_LogEdgeData = null;

            this.m_ChefBasis = null;

            this.m_Parallel = null;
            this.m_VertexData = null;
            this.aggregationGrid = null;
            this.CellPartitioning = null;
        }

        bool m_IsAlive = true;

        /// <summary>
        /// indicates that <see cref="Invalidate"/> has been called
        /// </summary>
        public bool IsAlive() {
            return m_IsAlive;
        }
    }
}

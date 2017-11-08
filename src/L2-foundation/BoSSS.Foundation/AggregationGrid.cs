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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Comm;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation.Caching;

namespace BoSSS.Foundation.Grid.Aggregation {
    public partial class AggregationGrid : IGridData {

        /// <summary>
        /// The grid from which this was coarsened.
        /// </summary>
        public IGridData ParentGrid {
            get;
            private set;
        }

        /// <summary>
        /// The ancestor grid, from which the aggregation sequence was derived.
        /// </summary>
        public IGridData AncestorGrid {
            get {
                if(ParentGrid is AggregationGrid) {
                    IGridData Ancestor = ((AggregationGrid)ParentGrid).AncestorGrid;
                    Debug.Assert(this.iGeomCells.NoOfCells == Ancestor.iGeomCells.NoOfCells);
                    return Ancestor;
                } else {
                    return ParentGrid;
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
        /// Constructor.
        /// </summary>
        /// <param name="pGrid">
        /// Parrent grid.
        /// </param>
        /// <param name="AggregationCells">
        /// Coarse cells which build up the fine cells.
        /// - 1st index: coarse (i.e. this) grid cell index
        /// - 2nd index: enumeration
        /// - content: local cell index into the parrent grid <paramref name="pGrid"/>.
        /// </param>
        public AggregationGrid(IGridData pGrid, int[][] AggregationCells) {
            ParentGrid = pGrid;

            int JlocFine = pGrid.iLogicalCells.NoOfLocalUpdatedCells;
            int JElocFine = pGrid.iLogicalCells.NoOfCells;

            m_GeomCellData = new GeomCellData() { m_Owner = this };
            m_LogicalCellData = new LogicalCellData() { m_Owner = this };
            m_GeomEdgeData = new GeomEdgeData() { m_Owner = this };
            m_LogEdgeData = new LogEdgeData();
            m_VertexData = new VertexData();
            m_Parallel = new Parallelization() { m_owner = this };
            
            CellPartitioning = new Partitioning(AggregationCells.Length, pGrid.CellPartitioning.MPI_Comm);
            
            int j0Coarse = CellPartitioning.i0;

            BuildNeighborship(AggregationCells);
            DefineCellParts();
            CollectEdges();
        }


        private void CollectEdges() {

            int[] F2C = this.jCellFine2jCellCoarse;
            
            int JlocCoarse = iLogicalCells.NoOfLocalUpdatedCells;
            int JElocCoarse = iLogicalCells.NoOfCells;

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
            for(int iEdgeFine = 0; iEdgeFine < NoOfEdgesFine; iEdgeFine++) { // loop over all logical edges in the fine grid.
                int jCellFine1 = E2Cfine[iEdgeFine, 0];
                int jCellFine2 = E2Cfine[iEdgeFine, 1];

                int jCellCoarse1 = F2C[jCellFine1];
                int jCellCoarse2 = jCellFine2 >= 0 ? F2C[jCellFine2] : -7544890;


                if (jCellCoarse2 >= 0) {
                    // +++++++++++++
                    // internal edge
                    // +++++++++++++

                    // does the edge exists already?
                    EdgeTmp edgTmp = null;
                    foreach (int i in tmpC2E[jCellCoarse1]) {
                        int iEdge = Math.Abs(i) - 1;
                        EdgeTmp _edgTmp = tmpEdges[iEdge];
                        if (   (_edgTmp.jCell1 == jCellCoarse1 && _edgTmp.jCell2 == jCellCoarse2) 
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
        }

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

        /// <summary>
        /// Sets up the decomposition of the aggregate cell into elementary parts, which can be mapped to reference elements.
        /// </summary>
        private void DefineCellParts() {
            IGridData pGrid = ParentGrid;
            int J = this.iLogicalCells.NoOfLocalUpdatedCells;
            int JE = this.iLogicalCells.NoOfCells;

            m_LogicalCellData.AggregateCellToParts = new int[JE][];
            var AggPartC = m_LogicalCellData.AggregateCellToParts;
            var AggPartF = pGrid.iLogicalCells.AggregateCellToParts;
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
            }
        }

        private void BuildNeighborship(int[][] AggregationCells) {
            IGridData pGrid = ParentGrid;
            int j0Coarse = CellPartitioning.i0;
            int JlocCoarse = CellPartitioning.LocalLength;
            int JElocFine = pGrid.iLogicalCells.NoOfCells;
            int JlocFine = pGrid.iLogicalCells.NoOfLocalUpdatedCells;

            // compute fine-to-coarse mapping
            // ==============================

            int[] Fine2CoarseGlobal = new int[JElocFine]; // index: local cell index on fine grid; maps to _global_ cell index on coarse grid.
            {
                ArrayTools.SetAll(Fine2CoarseGlobal, -111);
                for (int jCellCoarse = 0; jCellCoarse < JlocCoarse; jCellCoarse++) {
                    foreach (int jCellFine in AggregationCells[jCellCoarse]) {
                        if (jCellFine < 0 || jCellFine >= JlocFine)
                            throw new ArgumentOutOfRangeException();
                        Fine2CoarseGlobal[jCellFine] = jCellCoarse + j0Coarse;
                    }
                }


                Fine2CoarseGlobal.MPIExchange<int[], int>(pGrid);
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
                        jCoarse = Fine2CoarseGlobal[jF] - j0Coarse;
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
                            int jCellCoarseNeighGlob = Fine2CoarseGlobal[jCellFineNeigh]; // map index of fine grid neighbor to coarse cell index

                            // convert global cell index to local 
                            int jCellCoarseNeighLoc;
                            if (jCellCoarseNeighGlob >= j0Coarse && jCellCoarseNeighGlob <= (j0Coarse + JlocCoarse)) {
                                // Neighbor is a locally updated cell
                                jCellCoarseNeighLoc = jCellCoarseNeighGlob - j0Coarse;

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
                int[] TestData = new int[JElocCoarse];
                for (int jC = 0; jC < JlocCoarse; jC++) {
                    int GlobalIdx;
                    GlobalIdx = jC + j0Coarse;
                    TestData[jC] = GlobalIdx;
                }

                TestData.MPIExchange<int[], int>(this);

                for (int jC = 0; jC < JElocCoarse; jC++) {
                    int GlobalIdx;
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

        public GridData.BasisData ChefBasis {
            get {
                throw new NotImplementedException();
            }
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

        public Guid GridID {
            get {
                throw new NotImplementedException();
            }
        }

        public IDictionary<byte, string> EdgeTagNames {
            get {
                return ParentGrid.EdgeTagNames;
            }
        }
    }
}

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

using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.Grid.Classic
{
    partial class GridData {

        /// <summary>
        /// Creates a new grid, which is an adaptive refinement (cell by cell) of this grid.
        /// </summary>
        public GridCommons Adapt(IEnumerable<int> CellsToRefine, IEnumerable<int[]> CellsToCoarsen, out GridCorrelation Old2New) {
            using(new FuncTrace()) {
                GridCommons oldGrid = this.m_Grid;
                GridCommons newGrid = new GridCommons(oldGrid.RefElements, oldGrid.EdgeRefElements);
                Partitioning cellPartitioning = this.CellPartitioning;
                Old2New = new GridCorrelation();
                                
                int J = this.Cells.NoOfLocalUpdatedCells;
                int JE = this.Cells.NoOfExternalCells + this.Cells.NoOfLocalUpdatedCells;

                BitArray CellsToRefineBitmask = new BitArray(J);
                BitArray CellsToCoarseBitmask = new BitArray(J);
                BitArray AdaptNeighborsBitmask = new BitArray(J);

                // templates for subdivision
                // =========================
                RefElement[] KrefS = oldGrid.RefElements; // all ref elements used
                RefElement.SubdivisionTreeNode[] KrefS_SubDiv = new RefElement.SubdivisionTreeNode[KrefS.Length]; // subdivision tree for each ref element
                RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves = new RefElement.SubdivisionTreeNode[KrefS.Length][]; // actual subdivision elements
                Tuple<int, int>[][,] KrefS_SubdivConnections = new Tuple<int, int>[KrefS.Length][,]; // connections between elements; 1st idx: ref elem; 2nd idx: subdiv elm; 3rd idx: face of subdiv elm; content: [idx of subdiv elm,idx of face]
                int[][][] KrefS_Faces2Subdiv = new int[KrefS.Length][][]; // mapping: [ref elm, face of ref elm] -> Subdivision elements which bound to this face.
                Old2New.GeometricMapping = new AffineTrafo[KrefS.Length][];

                for(int iKref = 0; iKref < KrefS.Length; iKref++) {
                    RefElement Kref = KrefS[iKref];
                    KrefS_SubDiv[iKref] = Kref.GetSubdivisionTree(1);
                    KrefS_SubdivLeaves[iKref] = KrefS_SubDiv[0].GetLeaves();
                    Debug.Assert(ArrayTools.ListEquals(KrefS_SubdivLeaves[iKref], KrefS_SubDiv[iKref].Children[0].GetLevel(), (a, b) => object.ReferenceEquals(a, b)));
                    KrefS_Faces2Subdiv[iKref] = new int[Kref.NoOfFaces][];

                    KrefS_SubdivConnections[iKref] = new Tuple<int, int>[KrefS_SubdivLeaves[iKref].Length, KrefS[iKref].NoOfFaces];
                    for(int iSubdiv = 0; iSubdiv < KrefS_SubdivConnections[iKref].GetLength(0); iSubdiv++) { // loop over subdivision elements
                        for(int iFace = 0; iFace < KrefS_SubdivConnections[iKref].GetLength(1); iFace++) { // loop over faces of subdivision elements
                            var t = KrefS_SubdivLeaves[iKref][iSubdiv].GetNeighbor(iFace);
                            if(t.Item1 < 0) {
                                // at the boundary of the subdivision
                                ArrayTools.AddToArray(iSubdiv, ref KrefS_Faces2Subdiv[iKref][t.Item2]);
                            }
                            KrefS_SubdivConnections[iKref][iSubdiv, iFace] = t;
                        }
                    }

                    Old2New.GeometricMapping[iKref] = new AffineTrafo[KrefS_SubdivLeaves[iKref].Length];
                    for(int iSubDiv = 0; iSubDiv < KrefS_SubdivLeaves[iKref].Length; iSubDiv++) {
                        Old2New.GeometricMapping[iKref][iSubDiv] = KrefS_SubdivLeaves[iKref][iSubDiv].TrafoFromRoot;
                    }
                }

                Old2New.KrefS_SubdivLeaves = KrefS_SubdivLeaves;

                BitArray AdaptNeighboursOtherProcess = new BitArray(JE);
                List<long>[] exchangeNeighbours = new List<long>[MpiSize];

                long[] externalCellsGlobalIndices = iParallel.GlobalIndicesExternalCells;
                int firstGlobalIndex = cellPartitioning.i0;
                // Check Input, set Bitmasks
                // =========================

                if (CellsToRefine != null)
                {
                    ThrowExceptionForDoubleEntryInEnumeration(CellsToRefine);
                    CellsToRefineBitmask = GetCellsToRefine(CellsToRefine);
                    GetLocalAndExternalNeighbourCells(CellsToRefine, ref AdaptNeighborsBitmask);
                }

                Debugger.Launch();
                if (CellsToCoarsen != null) {
                    foreach(int[] coarseningCluster in CellsToCoarsen)
                    {
                        Cell[] coarseningCellCluster = coarseningCluster.Select(j => Cells.GetCell(j)).ToArray();

                        int CoarseningClusterID = coarseningCellCluster[0].CoarseningClusterID;
                        int iKref = Cells.GetRefElementIndex(coarseningCluster[0]);
                        int RefinementLevel = coarseningCellCluster[0].RefinementLevel;
                        
                        if(coarseningCluster.Length != KrefS_SubdivLeaves[iKref].Length)
                            throw new ArgumentException("Number of elements in coarsening cluster does not match refinement template for respective element type.");
                        if(RefinementLevel <= 0 || CoarseningClusterID <= 0)
                            throw new ArgumentException("Coarsening not available for respective cell.");

                        if(coarseningCellCluster.Where(cl => cl.ParentCell != null).Count() != 1)
                            throw new ArgumentException("Coarsening cluster seems wrong, or internal data may be corrupted.");

                        for(int z = 0; z < coarseningCellCluster.Length; z++) {
                            int j = coarseningCluster[z];

                            if(CellsToRefineBitmask[j] == true)
                                throw new ArgumentException("Cannot refine and coarsen the same cell.");
                            if(CellsToCoarseBitmask[j] == true)
                                throw new ArgumentException("Double entry.", "CellsToCoarsen");
                            CellsToCoarseBitmask[j] = true;

                            Cell Cj = Cells.GetCell(j);
                            if(CoarseningClusterID != Cj.CoarseningClusterID)
                                throw new ArgumentException("Mismatch of 'CoarseningClusterID' within cluster.");

                            this.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] Neighs, out _);

                            foreach (int jNeigh in Neighs) {
                                if(Array.IndexOf(coarseningCluster, jNeigh) < 0) {
                                    AdaptNeighborsBitmask[jNeigh] = true;
                                }
                            }
                        }
                    }
                    GetLocalAndExternalNeighbourCells(CellsToCoarsen, ref AdaptNeighborsBitmask);
                }

                BitArray CellsToRefineBitmaskExchange = new BitArray(JE);
                BitArray CellsToCoarseBitmaskExchange = new BitArray(JE);
                for (int j = 0; j < J; j++)
                {
                    CellsToRefineBitmaskExchange[j] = CellsToRefineBitmask[j];
                    CellsToCoarseBitmaskExchange[j] = CellsToCoarseBitmask[j];
                }
                CellsToRefineBitmaskExchange.MPIExchange(this);
                CellsToCoarseBitmaskExchange.MPIExchange(this);


                int InsertCounter = J;
                int test = MpiRank;
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);


                // create new cells
                // ================

                //Debug.Assert(this.MpiSize == 1, "still need to adjust the following lines.");

                long GlobalIdCounter = oldGrid.NumberOfCells_l;
                Convert.ToDouble(GlobalIdCounter).MPISum();
                Convert.ToInt64(GlobalIdCounter);
                List<Cell> newCells = new List<Cell>();
                int newVertexCounter = (oldGrid.Cells.Max(cl => cl.NodeIndices.Max()) + 1).MPIMax();
                Cell[][] adaptedCells = new Cell[J][];

                int globalJ = J.MPISum();
                long[] tempGlobalIDs = CurrentGlobalIdPermutation.Values.CloneAs();
                long[][] exchangeGlobalIDs = tempGlobalIDs.MPIGatherO(0);
                long[] oldGlobalIDs = new long[globalJ];
                if(MpiRank == 0)
                {
                    int mpiOffset = 0;
                    for (int m = 0; m < MpiSize; m++)
                    {
                        for (int j = 0; j < exchangeGlobalIDs[m].Length; j++)
                        {
                            oldGlobalIDs[mpiOffset + j] = exchangeGlobalIDs[m][j];
                        }
                        mpiOffset += exchangeGlobalIDs[m].Length;
                    }
                }
                oldGlobalIDs = oldGlobalIDs.MPIBroadcast(0);

                Old2New.OldGlobalId = CurrentGlobalIdPermutation.Values.CloneAs();
                Old2New.MappingIndex = new int[J][];
                Old2New.DestGlobalId = new long[J][];

                // clone neighbors of refined/coarsened cells
                // ------------------------------------------
                for (int j = 0; j < J; j++) {
                    Debug.Assert(Old2New.OldGlobalId[j] == this.Cells.GetCell(j).GlobalID);
                    Debug.Assert(Old2New.OldGlobalId[j] == oldGrid.Cells[j].GlobalID);
                    Debug.Assert(object.ReferenceEquals(this.Cells.GetCell(j), oldGrid.Cells[j]));

                    Debug.Assert((CellsToRefineBitmask[j] && CellsToCoarseBitmask[j]) == false, "Cannot refine and coarsen the same cell.");

                    if((CellsToRefineBitmask[j] || CellsToCoarseBitmask[j]) == false) {
                        if(AdaptNeighborsBitmask[j]) {
                            // neighbor information needs to be updated

                            Cell oldCell = oldGrid.Cells[j];
                            Cell newCell = oldCell.CloneAs(); // data 
                            newCells.Add(newCell);
                            adaptedCells[j] = new Cell[] { newCell };
                            
                            // remove out-dated neighborship info
                            if(newCell.CellFaceTags != null && newCell.CellFaceTags.Length > 0) {
                                int[] oldNeighs = this.Cells.CellNeighbours[j];
                                foreach(int jNeigh in oldNeighs) {
                                    if(CellsToRefineBitmaskExchange[jNeigh] || CellsToCoarseBitmaskExchange[jNeigh]) {
                                        // one of the neighbors has changed, so _potentially_ the cell face tags have to be updated
                                        long gId_Neigh = this.Cells.GetGlobalID(jNeigh);

                                        for(int i = 0; i < newCell.CellFaceTags.Length; i++) {
                                            if(newCell.CellFaceTags[i].NeighCell_GlobalID == gId_Neigh) {
                                                Debug.Assert(newCell.CellFaceTags[i].EdgeTag == 0 || newCell.CellFaceTags[i].EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG);
                                                ArrayTools.RemoveAt(ref newCell.CellFaceTags, i);
                                                i--;
                                            }
                                        }
                                    }
                                }
                            }

                        } else {
                            // cell and neighbors remain unchanged
                            newCells.Add(oldGrid.Cells[j]);
                        }

                        Debug.Assert(Old2New.MappingIndex[j] == null);
                        Debug.Assert(Old2New.DestGlobalId[j] == null);
                        Old2New.MappingIndex[j] = null;
                        Old2New.DestGlobalId[j] = new long[] { newCells[newCells.Count - 1].GlobalID };

                    } else {
                        Debug.Assert(CellsToRefineBitmask[j] || CellsToCoarseBitmask[j]);
                    }
                }
                
                // coarsening
                // ----------
                int bCoarsened = 0;
                if(CellsToCoarsen != null) {
                    
                    foreach(int[] jCellS in CellsToCoarsen) {
                        bCoarsened = 0xFFFF;

                        // cluster of cells to coarsen
                        Cell[] CellS = jCellS.Select(j => this.Cells.GetCell(j)).ToArray();
                        Debug.Assert(jCellS.Length == CellS.Length);

                        int RefinementLevel = CellS[0].RefinementLevel - 1;
                        if (RefinementLevel < 0)
                            throw new ArgumentException("Refinement level out of range - corrupted data structure.");
                        foreach(var cl in CellS) {
                            if(cl.RefinementLevel != RefinementLevel+1)
                                throw new ArgumentException("Refinement varies within refinement cluster - corrupted data structure.");
                        }

                        Cell Cell0 = CellS.Single(cl => cl.ParentCell != null);
                        Cell Mother = Cell0.ParentCell;
                        Debug.Assert(CellS.Where(cl => cl.GlobalID == Mother.GlobalID).Count() == 1);
                        Debug.Assert(Mother.RefinementLevel == RefinementLevel);
                       
                        Cell restoredCell = new Cell();
                        restoredCell.Type = Mother.Type;
                        Debug.Assert(Mother.Type == Cell0.Type);
                        restoredCell.NodeIndices = Mother.NodeIndices;
                        restoredCell.CoarseningClusterID = Mother.CoarseningClusterID;
                        restoredCell.ParentCell = Mother.ParentCell;
                        restoredCell.GlobalID = Cell0.GlobalID;
                        restoredCell.TransformationParams = Mother.TransformationParams;
                        restoredCell.RefinementLevel = RefinementLevel;
                        restoredCell.CoarseningClusterSize = Mother.CoarseningClusterSize;
                        restoredCell.CoarseningLeafIndex = Mother.CoarseningLeafIndex;

                        // boundary conditions by cell face tags
                        if (Mother.CellFaceTags != null) {
                            restoredCell.CellFaceTags = Mother.CellFaceTags.Where(cftag => cftag.EdgeTag > 0 && cftag.EdgeTag < GridCommons.FIRST_PERIODIC_BC_TAG).ToArray();
                        }

                        for(int iSubDiv = 0; iSubDiv < jCellS.Length; iSubDiv++) {
                            int j = jCellS[iSubDiv];
                            Cell Cj = CellS[iSubDiv];
                            Debug.Assert(adaptedCells[j] == null);
                            adaptedCells[j] = new[] { restoredCell };

                            Debug.Assert(Old2New.MappingIndex[j] == null);
                            Debug.Assert(Old2New.DestGlobalId[j] == null);
                            Old2New.MappingIndex[j] = new int[] { Cj.CoarseningLeafIndex };
                            Old2New.DestGlobalId[j] = new long[] { restoredCell.GlobalID };
                        }
                        
                        newCells.Add(restoredCell);

                    }
                }

                // refinement
                // ----------

                double[] SendOffset = new double[1];
                SendOffset[0] = CellsToRefine.Count();
                double[] ReceiveOffset = new double[MpiSize];
                unsafe
                {
                    fixed (double* pCheckSend = SendOffset, pCheckReceive = ReceiveOffset)
                    {
                        csMPI.Raw.Allgather((IntPtr)pCheckSend, SendOffset.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, SendOffset.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                    }
                }
                int globalIDOffset = 0;
                for (int m = 0; m < MpiRank; m++)
                {
                    globalIDOffset += Convert.ToInt32(ReceiveOffset[m]) * 3;
                }

                //Debugger.Launch();


                long[][] cellMapping = new long[globalJ][];
                if (CellsToRefine != null)
                {
                    int NewCoarseningClusterId;
                    {
                        int[] locData = new int[]{
                           (this.m_Grid.Cells.Max(cl => cl.CoarseningClusterID)),
                           (CellsToRefine.Count() + 1)
                        };
                        int[] glbData = locData.MPIMax();

                        NewCoarseningClusterId = glbData[0] + 1 + glbData[1] * this.MpiRank;
                    }
                    foreach (int j in CellsToRefine)
                    {
                        long oldGlobalID = oldGrid.Cells[j].GlobalID;
                        
                        Cell oldCell = oldGrid.Cells[j];
                        int iKref = this.Cells.GetRefElementIndex(j);
                        RefElement Kref = KrefS[iKref];
                        RefElement.SubdivisionTreeNode[] Leaves = KrefS_SubdivLeaves[iKref];
                        Tuple<int, int>[,] Connections = KrefS_SubdivConnections[iKref];

                        NodeSet RefNodes = Kref.GetInterpolationNodes(oldCell.Type);

                        Debug.Assert(adaptedCells[j] == null);
                        Cell[] refinedCells = new Cell[Leaves.Length];
                        adaptedCells[j] = refinedCells;

                        Debug.Assert(Old2New.MappingIndex[j] == null);
                        Old2New.MappingIndex[j] = new int[Leaves.Length];
                        cellMapping[oldGlobalID] = new long[Leaves.Length];
                        for (int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 1: create new cells
                            // create new cell
                            Cell newCell = new Cell();
                            newCell.Type = oldCell.Type;
                            if(iSubDiv == 0) {
                                newCell.GlobalID = oldCell.GlobalID;
                                newCell.ParentCell = oldCell.CloneAs();
                            }
                            else
                            {
                                long currentID = GlobalIdCounter + globalIDOffset;
                                GlobalIdCounter++;
                                newCell.GlobalID = currentID;
                            }
                            newCell.RefinementLevel = oldCell.RefinementLevel + 1;
                            newCell.CoarseningClusterSize = Leaves.Length;
                            newCell.CoarseningClusterID = NewCoarseningClusterId;
                            newCell.CoarseningLeafIndex = iSubDiv;
                            refinedCells[iSubDiv] = newCell;

                            // Vertices
                            var RefNodesRoot = Leaves[iSubDiv].Trafo2Root.Transform(RefNodes);
                            newCell.TransformationParams = MultidimensionalArray.Create(RefNodes.Lengths);
                            this.TransformLocal2Global(RefNodesRoot, newCell.TransformationParams, j);

                            // node indices
                            newCell.NodeIndices = new int[Kref.NoOfVertices];
                            for(int i = 0; i < Kref.NoOfVertices; i++) {
                                newCell.NodeIndices[i] = newVertexCounter + i + Convert.ToInt32(newCell.GlobalID) * Kref.NoOfVertices;
                                //newVertexCounter++;
                            }

                            // correlation
                            cellMapping[oldGlobalID][iSubDiv] = newCell.GlobalID;
                            Old2New.MappingIndex[j][iSubDiv] = iSubDiv;
                        }
                        NewCoarseningClusterId++;

                        for(int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 2: do other things

                            // neighbors within
                            for(int iFace = 0; iFace < Kref.NoOfFaces; iFace++) {
                                int iSubDiv_Neigh = Connections[iSubDiv, iFace].Item1;
                                if(iSubDiv_Neigh >= 0) {
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        ConformalNeighborship = true,
                                        NeighCell_GlobalID = refinedCells[Connections[iSubDiv, iFace].Item1].GlobalID,
                                        FaceIndex = iFace
                                    }, ref refinedCells[iSubDiv].CellFaceTags);
                                }
                            }
                        }

                        newCells.AddRange(refinedCells);
                        Debug.Assert(Old2New.DestGlobalId[j] == null);
                        Old2New.DestGlobalId[j] = refinedCells.Select(Cl => Cl.GlobalID).ToArray();

                    }
                }
                
                newGrid.Cells = newCells.ToArray();
                var CafafaasfasffNglb = newGrid.GetCellNeighbourship(true);


                // fix neighborship
                // ================

                byte[,] Edge2Face = this.Edges.FaceIndices;
                int[,] Edge2Cell = this.Edges.CellIndices;
                MultidimensionalArray[] VerticesFor_KrefEdge = this.Edges.EdgeRefElements.Select(KrefEdge => KrefEdge.Vertices).ToArray();

                int[] ONE_NULL = new int[] { 0 };

                int NoOfEdges = this.Edges.Count;
                Debug.Assert(Edge2Face.GetLength(0) == NoOfEdges);
                Debug.Assert(Edge2Cell.GetLength(0) == NoOfEdges);

                List<Tuple<int, Cell[]>> cellsOnNeighbourProcess = ExchangeCellData(adaptedCells);

                for (int iEdge = 0; iEdge < NoOfEdges; iEdge++)
                {
                    int jCell1 = Edge2Cell[iEdge, 0];
                    int jCell2 = Edge2Cell[iEdge, 1];

                    if (jCell2 < 0)
                    {
                        Debug.Assert((CellsToRefineBitmask[jCell1] && CellsToCoarseBitmask[jCell1]) == false);
                        if ((CellsToRefineBitmask[jCell1] || CellsToCoarseBitmask[jCell1]) == false)
                            continue;

                        Cell[] adaptedBCells1 = adaptedCells[jCell1];
                        Debug.Assert(adaptedBCells1 != null);

                        int iBFace = Edge2Face[iEdge, 0];

                        foreach (Cell cl in adaptedBCells1)
                        {
                            if (cl.CellFaceTags.Where(cft => cft.FaceIndex == iBFace).Count() == 0 && this.Edges.EdgeTags[iEdge] > 0)
                            {
                                ArrayTools.AddToArray(new CellFaceTag()
                                {
                                    EdgeTag = this.Edges.EdgeTags[iEdge],
                                    ConformalNeighborship = false,
                                    NeighCell_GlobalID = long.MinValue,
                                    FaceIndex = iBFace
                                }, ref cl.CellFaceTags);
                            }
                        }
                        continue;
                    }

                    Debug.Assert((CellsToRefineBitmaskExchange[jCell1] && CellsToCoarseBitmaskExchange[jCell1]) == false);
                    Debug.Assert((CellsToRefineBitmaskExchange[jCell2] && CellsToCoarseBitmaskExchange[jCell2]) == false);

                    bool C1changed = CellsToRefineBitmaskExchange[jCell1] || CellsToCoarseBitmaskExchange[jCell1];
                    bool C2changed = CellsToRefineBitmaskExchange[jCell2] || CellsToCoarseBitmaskExchange[jCell2];

                    if ((C1changed || C2changed) == false)
                        // edge between two un-changed cells -- this neighborship remains the same.
                        continue;

                    Cell[] adaptedCells1 = new Cell[1];
                    Cell[] adaptedCells2 = new Cell[1];

                    if (IsPartOfLocalCells(J, jCell1) && IsPartOfLocalCells(J, jCell2))
                    {
                        adaptedCells1 = adaptedCells[jCell1];
                        adaptedCells2 = adaptedCells[jCell2];
                    }
                    else if (IsPartOfLocalCells(J, jCell1) && !IsPartOfLocalCells(J, jCell2))
                    {
                        adaptedCells1 = adaptedCells[jCell1];
                        adaptedCells2 = FindCellOnNeighbourProcess(cellsOnNeighbourProcess, jCell2);
                    }
                    else if (!IsPartOfLocalCells(J, jCell1) && IsPartOfLocalCells(J, jCell2))
                    {
                        adaptedCells1 = FindCellOnNeighbourProcess(cellsOnNeighbourProcess, jCell1);
                        adaptedCells2 = adaptedCells[jCell2];
                    }
                    else
                        throw new Exception("Both cells not on the current process");

                    ThrowExceptionIfCellIsMissing(jCell1, adaptedCells1);
                    ThrowExceptionIfCellIsMissing(jCell2, adaptedCells2);

                    

                    if (CellsToCoarseBitmaskExchange[jCell1] && CellsToCoarseBitmaskExchange[jCell2])
                    {
                        Debug.Assert(adaptedCells1.Length == 1);
                        Debug.Assert(adaptedCells2.Length == 1);
                        if (adaptedCells1[0].GlobalID == adaptedCells2[0].GlobalID)
                        {
                            // these two cells will be joint into one cell -> no new neighborship
                            Debug.Assert(ReferenceEquals(adaptedCells1[0], adaptedCells2[0]));
                            continue;
                        }
                    }

                    int iFace1 = Edge2Face[iEdge, 0];
                    int iFace2 = Edge2Face[iEdge, 1];

                    Debug.Assert((adaptedCells1.Length > 1) == (CellsToRefineBitmaskExchange[jCell1]));
                    Debug.Assert((adaptedCells2.Length > 1) == (CellsToRefineBitmaskExchange[jCell2]));

                    int iKref1 = this.Cells.GetRefElementIndex(jCell1);
                    int iKref2 = this.Cells.GetRefElementIndex(jCell2);
                    RefElement Kref1 = this.Cells.GetRefElement(jCell1);
                    RefElement Kref2 = this.Cells.GetRefElement(jCell2);

                    int[] idx1, idx2;
                    if (CellsToRefineBitmaskExchange[jCell1])
                    {
                        idx1 = KrefS_Faces2Subdiv[iKref1][iFace1];
                    }
                    else
                    {
                        Debug.Assert(adaptedCells1.Length == 1);
                        idx1 = ONE_NULL;
                    }

                    if (CellsToRefineBitmaskExchange[jCell2])
                    {
                        idx2 = KrefS_Faces2Subdiv[iKref2][iFace2];
                    }
                    else
                    {
                        Debug.Assert(adaptedCells2.Length == 1);
                        idx2 = ONE_NULL;
                    }

                    foreach (int i1 in idx1)
                    {

                        MultidimensionalArray VtxFace1;
                        if (CellsToRefineBitmaskExchange[jCell1])
                        {
                            VtxFace1 = KrefS_SubdivLeaves[iKref1][i1].GetFaceVertices(iFace1);
                        }
                        else
                        {
                            VtxFace1 = Kref1.GetFaceVertices(iFace1);
                        }

                        Cell Cl1 = adaptedCells1[i1];

                        foreach (int i2 in idx2)
                        {

                            Cell Cl2 = adaptedCells2[i2];
                            Debug.Assert(Cl1.GlobalID != Cl2.GlobalID);


                            int conCount1;
                            if (Cl1.CellFaceTags == null)
                            {
                                conCount1 = 0;
                            }
                            else
                            {
                                conCount1 = Cl1.CellFaceTags.Where(cfTag => cfTag.NeighCell_GlobalID == Cl2.GlobalID).Count();
                            }
                            Debug.Assert(conCount1 <= 1);
#if DEBUG
                            int conCount2;
                            if (Cl2.CellFaceTags == null)
                            {
                                conCount2 = 0;
                            }
                            else
                            {
                                conCount2 = Cl2.CellFaceTags.Where(cfTag => cfTag.NeighCell_GlobalID == Cl1.GlobalID).Count();
                            }
                            Debug.Assert(conCount1 == conCount2);
#endif                            
                            if (conCount1 > 0)
                                continue;

                            MultidimensionalArray VtxFace2;
                            {
                                MultidimensionalArray VtxFace2_L;
                                if (CellsToRefineBitmaskExchange[jCell2])
                                {
                                    VtxFace2_L = KrefS_SubdivLeaves[iKref2][i2].GetFaceVertices(iFace2);
                                }
                                else
                                {
                                    VtxFace2_L = Kref2.GetFaceVertices(iFace2);
                                }

                                MultidimensionalArray VtxFace2_G = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                                VtxFace2 = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                                this.TransformLocal2Global(VtxFace2_L, VtxFace2_G, jCell2);
                                bool[] Converged = new bool[VtxFace2_L.NoOfRows];
                                this.TransformGlobal2Local(VtxFace2_G, VtxFace2, jCell1, Converged);
                                if (Converged.Any(t => t == false))
                                    throw new ArithmeticException("Newton divergence");
                            }

                            bool bIntersect = GridData.EdgeData.FaceIntersect(VtxFace1, VtxFace2,
                                Kref1.GetFaceTrafo(iFace1), Kref1.GetInverseFaceTrafo(iFace1),
                                VerticesFor_KrefEdge,
                                out bool conformal1, out bool conformal2, out AffineTrafo newTrafo, out int Edg_idx);

                            if (bIntersect)
                            {
                                ArrayTools.AddToArray(new CellFaceTag()
                                {
                                    ConformalNeighborship = false,
                                    NeighCell_GlobalID = Cl2.GlobalID,
                                    FaceIndex = iFace1
                                }, ref Cl1.CellFaceTags);

                                ArrayTools.AddToArray(new CellFaceTag()
                                {
                                    ConformalNeighborship = false,
                                    NeighCell_GlobalID = Cl1.GlobalID,
                                    FaceIndex = iFace2
                                }, ref Cl2.CellFaceTags);
                            }
                        }
                    }
                }

                // add EdgeTagNames and periodic Transformations
                for (int etCnt = 1; etCnt < this.EdgeTagNames.Count; etCnt++)
                {
                    KeyValuePair<byte, string> etPair = this.EdgeTagNames.ElementAt(etCnt);
                    newGrid.EdgeTagNames.Add(etPair);
                }

                foreach (AffineTrafo trafo in this.Grid.PeriodicTrafo)
                {
                    newGrid.PeriodicTrafo.Add(trafo);
                }
                foreach (AffineTrafo itrafo in this.Grid.InversePeriodicTrafo)
                {
                    newGrid.InversePeriodicTrafo.Add(itrafo);
                }



                // finalize
                // ========
                int bCoarsenedGlobal = bCoarsened.MPIMax();
                if(bCoarsenedGlobal > 0) {
#if DEBUG
                    if(this.MpiSize == 1) {
                        List<int> allgids = new List<int>();
                        foreach(var cl in newGrid.Cells) {
                            allgids.Add((int)(cl.GlobalID));
                        }
                        bool[] markers = new bool[allgids.Max() + 1];
                        for(int i = 0; i < allgids.Count; i++) {
                            long gid = allgids[i];
                            Debug.Assert(markers[gid] == false, "Some GlobalID is used twice.");
                            markers[gid] = true;
                        }

                        foreach(var cl in newGrid.Cells) {
                            if(cl.CellFaceTags != null) {
                                for(int i = 0; i < cl.CellFaceTags.Length; i++) {
                                    long ngid = cl.CellFaceTags[i].NeighCell_GlobalID;
                                    if (ngid >= 0)
                                        Debug.Assert(markers[ngid] == true);
                                }
                            }
                        }
                    }
#endif

                    List<long> old2NewGid = new List<long>();
                    Debug.Assert(Old2New.DestGlobalId.Length == J);
                    for(int j = 0; j < J; j++) {
                        old2NewGid.AddRange(Old2New.DestGlobalId[j]);
                    }

                    newGrid.CompressGlobalID(old2NewGid);

                    int c2 = 0;
                    for(int j = 0; j < J; j++) {
                        long[] o2nj = Old2New.DestGlobalId[j];
                        int K = o2nj.Length;
                        for(int k = 0; k < K; k++) {
                            o2nj[k] = old2NewGid[c2];
                            c2++;
                        }
                    }
                    Debug.Assert(c2 == old2NewGid.Count);
                }

                var CNglb = newGrid.GetCellNeighbourship(true);
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
                return newGrid;
            }
        }

        private static void ThrowExceptionIfCellIsMissing(int jCell, Cell[] adaptedCells)
        {
            if (adaptedCells == null)
            {
                throw new Exception("Cell " + jCell + "does not exist!");
            }
        }

        private Cell[] FindCellOnNeighbourProcess(List<Tuple<int, Cell[]>> cellsOnNeighbourProcess, int localCellIndex)
        {
            int cellDivision2D = 4;
            Cell[] adaptedCell = new Cell[cellDivision2D];
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            long[] externalCellsGlobalIndices = iParallel.GlobalIndicesExternalCells;
            int jCell2GlobalIndex = (int)externalCellsGlobalIndices[localCellIndex - noOfLocalCells];
            for (int j = 0; j < cellsOnNeighbourProcess.Count(); j++)
            {
                if (jCell2GlobalIndex == cellsOnNeighbourProcess[j].Item1)
                {
                    adaptedCell = cellsOnNeighbourProcess[j].Item2;
                    break;
                }
            }

            return adaptedCell;
        }

        private List<Tuple<int, Cell[]>> ExchangeCellData(Cell[][] Cells)
        {
            List<Tuple<int, Cell[]>> exchangedCellData = new List<Tuple<int, Cell[]>>();
            Dictionary<int, List<Tuple<int, Cell[]>>> sendCellData = GetBoundaryCellsAndProcessToSend(Cells);
            IDictionary<int, List<Tuple<int, Cell[]>>> receiveCellData = SerialisationMessenger.ExchangeData(sendCellData);

            foreach (KeyValuePair<int, List<Tuple<int, Cell[]>>> kv in receiveCellData)
            {
                List<Tuple<int, Cell[]>> list = kv.Value;

                foreach (Tuple<int, Cell[]> t in list)
                { 
                    int globalCellIndex = t.Item1;
                    Cell[] cellOnOtherProc = t.Item2;

                    int localCellIndex = globalCellIndex + CellPartitioning.i0;

                    exchangedCellData.Add(new Tuple<int, Cell[]>(globalCellIndex, cellOnOtherProc));
                }
            }

            return exchangedCellData;
        }

        private Dictionary<int, List<Tuple<int, Cell[]>>> GetBoundaryCellsAndProcessToSend(Cell[][] Cells)
        {
            Dictionary<int, List<Tuple<int, Cell[]>>> boundaryCellsProcess = new Dictionary<int, List<Tuple<int, Cell[]>>>();
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;

            for (int j = 0; j < noOfLocalCells; j++)
            {
                this.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] neighbourCells, out _);

                for (int i = 0; i < neighbourCells.Length; i++)
                {
                    if (!IsPartOfLocalCells(noOfLocalCells, neighbourCells[i]))
                    {
                        int currentProcess = GetExternalCellProcess(neighbourCells[i]);
                        int globalCellIndex = CellPartitioning.i0 + j;
                        if (!boundaryCellsProcess.TryGetValue(currentProcess, out List<Tuple<int, Cell[]>> exchangeCellData))
                        {
                            exchangeCellData = new List<Tuple<int, Cell[]>>();
                            boundaryCellsProcess.Add(currentProcess, exchangeCellData);
                        }
                        foreach (KeyValuePair<int, List<Tuple<int, Cell[]>>> kv in boundaryCellsProcess)
                        {
                            List<Tuple<int, Cell[]>> list = kv.Value;
                            bool processAlreadyIncluded = false;
                            foreach (Tuple<int, Cell[]> t in list)
                            {
                                if (t.Item1 == j)
                                {
                                    processAlreadyIncluded = true;
                                    break;
                                }
                            }
                            if(!processAlreadyIncluded)
                                exchangeCellData.Add(new Tuple<int, Cell[]>(globalCellIndex, Cells[j]));
                        }
                    }
                }
            }

            return boundaryCellsProcess;
        }

        private int GetExternalCellProcess(int externalCellIndex)
        {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            long[] externalCellsGlobalIndices = iParallel.GlobalIndicesExternalCells;

            int globalCellIndex = (int)externalCellsGlobalIndices[externalCellIndex - noOfLocalCells];
            int currentProcess = CellPartitioning.FindProcess(globalCellIndex);

            return currentProcess;
        }
        
        private void ThrowExceptionForDoubleEntryInEnumeration(IEnumerable<int> enumeration)
        {
            for(int i = 0; i < enumeration.Count(); i++)
            {
                int currentEntry = enumeration.ElementAt(i);
                for (int j = i + 1; j < enumeration.Count(); j++)
                {
                    if(currentEntry == enumeration.ElementAt(j))
                        throw new ArgumentException("Double entry.", "CellsToRefine");
                }
            }
        }

        private BitArray GetCellsToRefine(IEnumerable<int> cellsToRefine)
        {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            BitArray cellsToRefineBitMask = new BitArray(noOfLocalCells);

            foreach (int currentCellLocalIndex in cellsToRefine)
            {
                cellsToRefineBitMask[currentCellLocalIndex] = true;
            }

            return cellsToRefineBitMask;
        }

        private void GetLocalAndExternalNeighbourCells(IEnumerable<int> cellsToRefine, ref BitArray AdaptNeighborsBitmask)
        {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            long[] externalCellsGlobalIndices = this.iParallel.GlobalIndicesExternalCells;
            List<long>[] exchangeNeighbours = new List<long>[MpiSize];

            foreach (int currentCellIndex in cellsToRefine)
            {
                this.GetCellNeighbours(currentCellIndex, GetCellNeighbours_Mode.ViaEdges, out int[] neighbourCells, out _);

                foreach (int neighbourCellIndex in neighbourCells)
                {
                    if (IsPartOfLocalCells(noOfLocalCells, neighbourCellIndex))
                        AdaptNeighborsBitmask[neighbourCellIndex] = true;
                    else
                    {
                        int neighbourGlobalIndex = (int)externalCellsGlobalIndices[neighbourCellIndex - noOfLocalCells];
                        int neighbourProcess = CellPartitioning.FindProcess(neighbourGlobalIndex);
                        if (exchangeNeighbours[neighbourProcess] == null)
                            exchangeNeighbours[neighbourProcess] = new List<long>();
                        exchangeNeighbours[neighbourProcess].Add(neighbourGlobalIndex);
                    }
                }
            }

            GetAndExchangeExternalNeighbours(exchangeNeighbours, ref AdaptNeighborsBitmask);
        }

        private void GetLocalAndExternalNeighbourCells(IEnumerable<int[]> coarseningClusters, ref BitArray AdaptNeighborsBitmask)
        {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            long[] externalCellsGlobalIndices = this.iParallel.GlobalIndicesExternalCells;
            List<long>[] exchangeNeighbours = new List<long>[MpiSize];

            foreach(int[] coarseningClusterID in coarseningClusters)
            {
                foreach (int currentCellIndex in coarseningClusterID)
                {
                    this.GetCellNeighbours(currentCellIndex, GetCellNeighbours_Mode.ViaVertices, out int[] neighbourCells, out _);

                    foreach (int neighbourCellIndex in neighbourCells)
                    {
                        if (IsPartOfLocalCells(noOfLocalCells, neighbourCellIndex))
                            AdaptNeighborsBitmask[neighbourCellIndex] = true;
                        else
                        {
                            int neighbourGlobalIndex = (int)externalCellsGlobalIndices[neighbourCellIndex - noOfLocalCells];
                            int neighbourProcess = CellPartitioning.FindProcess(neighbourGlobalIndex);
                            if (exchangeNeighbours[neighbourProcess] == null)
                                exchangeNeighbours[neighbourProcess] = new List<long>();
                            exchangeNeighbours[neighbourProcess].Add(neighbourGlobalIndex);
                        }
                    }
                }
            }
            GetAndExchangeExternalNeighbours(exchangeNeighbours, ref AdaptNeighborsBitmask);
        }

        private void GetAndExchangeExternalNeighbours(List<long>[] exchangeNeighbours, ref BitArray AdaptNeighborsBitmask)
        {
            List<long> externalNeighbours = ExchangeExternalNeighbours(exchangeNeighbours);
            GetExternalNeighbours(externalNeighbours, ref AdaptNeighborsBitmask);
        }

        private List<long> ExchangeExternalNeighbours(List<long>[] exchangeNeighbours)
        {
            List<long> externalNeighbours = new List<long>();

            List<long>[][] tempExchange = exchangeNeighbours.MPIGatherO(0);
            tempExchange = tempExchange.MPIBroadcast(0);

            for (int m = 0; m < MpiSize; m++)
            {
                if (m != MpiRank)
                {
                    if (tempExchange[m][MpiRank] != null)
                        externalNeighbours.AddRange(tempExchange[m][MpiRank]);
                }
            }

            return externalNeighbours;
        }

        private void GetExternalNeighbours(List<long> externalNeighbours, ref BitArray AdaptNeighborsBitmask)
        {
            int firstGlobalIndex = this.CellPartitioning.i0;
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;

            for (int j = 0; j < externalNeighbours.Count(); j++)
            {
                int externalNeighbour = (int)externalNeighbours[j] - firstGlobalIndex;
                AdaptNeighborsBitmask[externalNeighbour] = true;
            }
        }

        private static bool IsPartOfLocalCells(int noOfLocalCells, int currentCellIndex)
        {
            return currentCellIndex < noOfLocalCells;
        }

        static private void TestCellIntersection(Cell A, Cell B) {
            Debug.Assert(A.TransformationParams.GetLength(1) == B.TransformationParams.GetLength(1));

            if(A.GlobalID == 408 && B.GlobalID == 525)
                Debugger.Break();
            if(B.GlobalID == 408 && A.GlobalID == 525)
                Debugger.Break();

            int D = A.TransformationParams.GetLength(1);
            BoundingBox BBA = new BoundingBox(A.TransformationParams);
            BoundingBox BBB = new BoundingBox(B.TransformationParams);

            BBA.ExtendByFactor(1e-8);
            BBB.ExtendByFactor(1e-8);

            if(!BBA.Overlap(BBB))
                throw new ApplicationException("Internal error.");
        }

      
    }
}

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
        public GridCommons Adapt(IEnumerable<int> cellsToRefine, IEnumerable<int[]> cellsToCoarse, out GridCorrelation Old2New) {
            using(new FuncTrace())
            {
                GridCommons oldGrid = this.m_Grid;
                GridCommons newGrid = new GridCommons(oldGrid.RefElements, oldGrid.EdgeRefElements);
                Partitioning cellPartitioning = this.CellPartitioning;
                Old2New = new GridCorrelation();

                int J = this.Cells.NoOfLocalUpdatedCells;
                int JE = this.Cells.NoOfExternalCells + this.Cells.NoOfLocalUpdatedCells;

                BitArray cellsToRefineBitmask = new BitArray(JE);
                BitArray cellsToCoarseBitmask = new BitArray(JE);
                BitArray AdaptNeighborsBitmask = new BitArray(J);

                // templates for subdivision
                // =========================
                RefElement[] KrefS = oldGrid.RefElements; // all ref elements used
                RefElement.SubdivisionTreeNode[] KrefS_SubDiv = new RefElement.SubdivisionTreeNode[KrefS.Length]; // subdivision tree for each ref element
                RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves = new RefElement.SubdivisionTreeNode[KrefS.Length][]; // actual subdivision elements
                Tuple<int, int>[][,] KrefS_SubdivConnections = new Tuple<int, int>[KrefS.Length][,]; // connections between elements; 1st idx: ref elem; 2nd idx: subdiv elm; 3rd idx: face of subdiv elm; content: [idx of subdiv elm,idx of face]
                int[][][] KrefS_Faces2Subdiv = new int[KrefS.Length][][]; // mapping: [ref elm, face of ref elm] -> Subdivision elements which bound to this face.
                Old2New.GeometricMapping = new AffineTrafo[KrefS.Length][];

                for (int iKref = 0; iKref < KrefS.Length; iKref++)
                {
                    RefElement Kref = KrefS[iKref];
                    KrefS_SubDiv[iKref] = Kref.GetSubdivisionTree(1);
                    KrefS_SubdivLeaves[iKref] = KrefS_SubDiv[0].GetLeaves();
                    Debug.Assert(ArrayTools.ListEquals(KrefS_SubdivLeaves[iKref], KrefS_SubDiv[iKref].Children[0].GetLevel(), (a, b) => object.ReferenceEquals(a, b)));
                    KrefS_Faces2Subdiv[iKref] = new int[Kref.NoOfFaces][];

                    KrefS_SubdivConnections[iKref] = new Tuple<int, int>[KrefS_SubdivLeaves[iKref].Length, KrefS[iKref].NoOfFaces];
                    for (int iSubdiv = 0; iSubdiv < KrefS_SubdivConnections[iKref].GetLength(0); iSubdiv++)
                    { // loop over subdivision elements
                        for (int iFace = 0; iFace < KrefS_SubdivConnections[iKref].GetLength(1); iFace++)
                        { // loop over faces of subdivision elements
                            var t = KrefS_SubdivLeaves[iKref][iSubdiv].GetNeighbor(iFace);
                            if (t.Item1 < 0)
                            {
                                // at the boundary of the subdivision
                                ArrayTools.AddToArray(iSubdiv, ref KrefS_Faces2Subdiv[iKref][t.Item2]);
                            }
                            KrefS_SubdivConnections[iKref][iSubdiv, iFace] = t;
                        }
                    }

                    Old2New.GeometricMapping[iKref] = new AffineTrafo[KrefS_SubdivLeaves[iKref].Length];
                    for (int iSubDiv = 0; iSubDiv < KrefS_SubdivLeaves[iKref].Length; iSubDiv++)
                    {
                        Old2New.GeometricMapping[iKref][iSubDiv] = KrefS_SubdivLeaves[iKref][iSubDiv].TrafoFromRoot;
                    }
                }

                Old2New.KrefS_SubdivLeaves = KrefS_SubdivLeaves;
                long GlobalIdCounter = GetGlobalIDCounter(oldGrid);
                int globalIDOffset = GetGlobalIdMPIOffset(cellsToRefine, cellsToCoarse);
                int newVertexCounter = (oldGrid.Cells.Max(cl => cl.NodeIndices.Max()) + 1).MPIMax();

                Old2New.OldGlobalId = CurrentGlobalIdPermutation.Values.CloneAs();
                Old2New.MappingIndex = new int[J][];
                Old2New.DestGlobalId = new long[J][];
                // Check Input, set Bitmasks
                // =========================
                if (!cellsToRefine.IsNullOrEmpty())
                {
                    ThrowExceptionForDoubleEntryInEnumeration(cellsToRefine);
                    cellsToRefineBitmask = GetCellBitMask(cellsToRefine);
                    GetLocalAndExternalNeighbourCells(cellsToRefine, ref AdaptNeighborsBitmask);
                }

                List<Cell> cellsInNewGrid = new List<Cell>();
                Cell[][] adaptedCells = new Cell[J][];
                int bCoarsened = 0;
                if (!cellsToCoarse.IsNullOrEmpty())
                {
                    ThrowExceptionsForErrorsInCoarseningCluster(cellsToCoarse, cellsToRefineBitmask, cellsToCoarseBitmask, KrefS_SubdivLeaves);
                    cellsToCoarseBitmask = GetCellBitMask(cellsToCoarse);
                    GetLocalAndExternalNeighbourCells(cellsToCoarse, ref AdaptNeighborsBitmask);

                    foreach (int[] cellClusterID in cellsToCoarse)
                    {
                        bCoarsened = 0xFFFF;
                        CoarseCells(Old2New, cellsInNewGrid, adaptedCells, cellClusterID);
                    }
                }

                cellsToRefineBitmask.MPIExchange(this);
                cellsToCoarseBitmask.MPIExchange(this);


                // clone neighbors of refined/coarsened cells
                // ------------------------------------------
                for (int j = 0; j < J; j++)
                {
                    Debug.Assert(Old2New.OldGlobalId[j] == this.Cells.GetCell(j).GlobalID);
                    Debug.Assert(Old2New.OldGlobalId[j] == oldGrid.Cells[j].GlobalID);
                    Debug.Assert(ReferenceEquals(Cells.GetCell(j), oldGrid.Cells[j]));
                    Debug.Assert((cellsToRefineBitmask[j] && cellsToCoarseBitmask[j]) == false, "Cannot refine and coarsen the same cell.");

                    if ((cellsToRefineBitmask[j] || cellsToCoarseBitmask[j]) == false)
                    {
                        if (AdaptNeighborsBitmask[j])
                        {
                            AdaptNeighboursOfChangedCellsToNewGrid(Old2New, oldGrid, cellsToRefineBitmask, cellsToCoarseBitmask, cellsInNewGrid, adaptedCells, j);
                        }
                        else
                        {
                            CloneAllOtherCellsToNewGrid(Old2New, oldGrid, cellsInNewGrid, j);
                        }
                    }
                    else
                    {
                        Debug.Assert(cellsToRefineBitmask[j] || cellsToCoarseBitmask[j]);
                    }
                }



                // refinement
                // ----------

                if (cellsToRefine != null)
                {
                    int NewCoarseningClusterId = GetNewCoarseningClusterID(cellsToRefine, oldGrid);

                    foreach (int j in cellsToRefine)
                    {
                        CheckForDoubleEntries(Old2New, adaptedCells, j);

                        Cell oldCell = oldGrid.Cells[j];
                        int iKref = Cells.GetRefElementIndex(j);
                        RefElement Kref = KrefS[iKref];
                        NodeSet RefNodes = Kref.GetInterpolationNodes(oldCell.Type);
                        RefElement.SubdivisionTreeNode[] Leaves = KrefS_SubdivLeaves[iKref];
                        Tuple<int, int>[,] Connections = KrefS_SubdivConnections[iKref];

                        Cell[] refinedCells = new Cell[Leaves.Length];
                        adaptedCells[j] = refinedCells;

                        Old2New.MappingIndex[j] = new int[Leaves.Length];
                        for (int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++)
                        {
                            Cell newCell = RefineCell(ref GlobalIdCounter, globalIDOffset, NewCoarseningClusterId, oldCell, Leaves, iSubDiv);
                            refinedCells[iSubDiv] = newCell;

                            GetNewVerticesOfRefinedCells(j, RefNodes, Leaves, iSubDiv, newCell);

                            newCell.NodeIndices = GetNodeIndicesOfRefinedCells(newVertexCounter, Kref, newCell);

                            // correlation
                            Old2New.MappingIndex[j][iSubDiv] = iSubDiv;
                        }
                        NewCoarseningClusterId++;
                        AdaptNeighbourshipWithinOldCell(Kref, Leaves, Connections, refinedCells);

                        Debug.Assert(Old2New.DestGlobalId[j] == null);
                        Old2New.DestGlobalId[j] = refinedCells.Select(Cl => Cl.GlobalID).ToArray();
                        cellsInNewGrid.AddRange(refinedCells);
                    }
                }

                newGrid.Cells = cellsInNewGrid.ToArray();


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
                        Debug.Assert((cellsToRefineBitmask[jCell1] && cellsToCoarseBitmask[jCell1]) == false);
                        if ((cellsToRefineBitmask[jCell1] || cellsToCoarseBitmask[jCell1]) == false)
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

                    Debug.Assert((cellsToRefineBitmask[jCell1] && cellsToCoarseBitmask[jCell1]) == false);
                    Debug.Assert((cellsToRefineBitmask[jCell2] && cellsToCoarseBitmask[jCell2]) == false);

                    bool C1changed = cellsToRefineBitmask[jCell1] || cellsToCoarseBitmask[jCell1];
                    bool C2changed = cellsToRefineBitmask[jCell2] || cellsToCoarseBitmask[jCell2];

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



                    if (cellsToCoarseBitmask[jCell1] && cellsToCoarseBitmask[jCell2])
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

                    Debug.Assert((adaptedCells1.Length > 1) == (cellsToRefineBitmask[jCell1]));
                    Debug.Assert((adaptedCells2.Length > 1) == (cellsToRefineBitmask[jCell2]));

                    int iKref1 = this.Cells.GetRefElementIndex(jCell1);
                    int iKref2 = this.Cells.GetRefElementIndex(jCell2);
                    RefElement Kref1 = this.Cells.GetRefElement(jCell1);
                    RefElement Kref2 = this.Cells.GetRefElement(jCell2);

                    int[] idx1, idx2;
                    if (cellsToRefineBitmask[jCell1])
                    {
                        idx1 = KrefS_Faces2Subdiv[iKref1][iFace1];
                    }
                    else
                    {
                        Debug.Assert(adaptedCells1.Length == 1);
                        idx1 = ONE_NULL;
                    }

                    if (cellsToRefineBitmask[jCell2])
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
                        if (cellsToRefineBitmask[jCell1])
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
                                if (cellsToRefineBitmask[jCell2])
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
                if (bCoarsenedGlobal > 0)
                {
#if DEBUG
                    if (this.MpiSize == 1)
                    {
                        List<int> allgids = new List<int>();
                        foreach (var cl in newGrid.Cells)
                        {
                            allgids.Add((int)(cl.GlobalID));
                        }
                        bool[] markers = new bool[allgids.Max() + 1];
                        for (int i = 0; i < allgids.Count; i++)
                        {
                            long gid = allgids[i];
                            Debug.Assert(markers[gid] == false, "Some GlobalID is used twice.");
                            markers[gid] = true;
                        }

                        foreach (var cl in newGrid.Cells)
                        {
                            if (cl.CellFaceTags != null)
                            {
                                for (int i = 0; i < cl.CellFaceTags.Length; i++)
                                {
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
                    for (int j = 0; j < J; j++)
                    {
                        old2NewGid.AddRange(Old2New.DestGlobalId[j]);
                    }

                    newGrid.CompressGlobalID(old2NewGid);

                    int c2 = 0;
                    for (int j = 0; j < J; j++)
                    {
                        long[] o2nj = Old2New.DestGlobalId[j];
                        int K = o2nj.Length;
                        for (int k = 0; k < K; k++)
                        {
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

        private static void AdaptNeighbourshipWithinOldCell(RefElement Kref, RefElement.SubdivisionTreeNode[] Leaves, Tuple<int, int>[,] Connections, Cell[] refinedCells)
        {
            for (int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++)
            {
                // neighbors within
                for (int iFace = 0; iFace < Kref.NoOfFaces; iFace++)
                {
                    int iSubDiv_Neigh = Connections[iSubDiv, iFace].Item1;
                    if (iSubDiv_Neigh >= 0)
                    {
                        ArrayTools.AddToArray(new CellFaceTag()
                        {
                            ConformalNeighborship = true,
                            NeighCell_GlobalID = refinedCells[Connections[iSubDiv, iFace].Item1].GlobalID,
                            FaceIndex = iFace
                        }, ref refinedCells[iSubDiv].CellFaceTags);
                    }
                }
            }
        }

        private static void CheckForDoubleEntries(GridCorrelation Old2New, Cell[][] adaptedCells, int j)
        {
            if (adaptedCells[j] != null)
                throw new Exception("Error in refinement algorithm: Cell was already changend.");
            if (Old2New.MappingIndex[j] != null)
                throw new Exception("Error in refinement algorithm: Mapping index already exists.");
        }

        private static Cell RefineCell(ref long GlobalIdCounter, int globalIDOffset, int NewCoarseningClusterId, Cell oldCell, RefElement.SubdivisionTreeNode[] Leaves, int iSubDiv)
        {
            Cell newCell = new Cell
            {
                Type = oldCell.Type,
                RefinementLevel = oldCell.RefinementLevel + 1,
                CoarseningClusterSize = Leaves.Length,
                CoarseningClusterID = NewCoarseningClusterId,
                CoarseningLeafIndex = iSubDiv
            };
            if (iSubDiv == 0)
            {
                newCell.GlobalID = oldCell.GlobalID;
                newCell.ParentCell = oldCell.CloneAs();
            }
            else
            {
                newCell.GlobalID = GlobalIdCounter + globalIDOffset;
                GlobalIdCounter++;
            }

            return newCell;
        }

        private void GetNewVerticesOfRefinedCells(int j, NodeSet RefNodes, RefElement.SubdivisionTreeNode[] Leaves, int iSubDiv, Cell newCell)
        {
            MultidimensionalArray RefNodesRoot = Leaves[iSubDiv].Trafo2Root.Transform(RefNodes);
            newCell.TransformationParams = MultidimensionalArray.Create(RefNodes.Lengths);
            TransformLocal2Global(RefNodesRoot, newCell.TransformationParams, j);
        }

        private static int[] GetNodeIndicesOfRefinedCells(int newVertexCounter, RefElement Kref, Cell newCell)
        {
            int[] tempNodeIndices = new int[Kref.NoOfVertices];
            for (int i = 0; i < Kref.NoOfVertices; i++)
            {
                tempNodeIndices[i] = newVertexCounter + i + (int)newCell.GlobalID * Kref.NoOfVertices;
            }

            return tempNodeIndices;
        }

        private int GetNewCoarseningClusterID(IEnumerable<int> cellsToRefine, GridCommons oldGrid)
        {
            int NewCoarseningClusterId;
            {
                int[] locData = new int[]
                {
                           (oldGrid.Cells.Max(cl => cl.CoarseningClusterID)),
                           (cellsToRefine.Count() + 1)
                };
                int[] glbData = locData.MPIMax();

                NewCoarseningClusterId = glbData[0] + 1 + glbData[1] * MpiRank;
            }

            return NewCoarseningClusterId;
        }

        private int GetGlobalIdMPIOffset(IEnumerable<int> cellsToRefine, IEnumerable<int[]> cellsToCoarse)
        {
            int globalIDOffset = 0;
            double[] SendOffset = new double[1];
            SendOffset[0] = cellsToRefine.Count() - cellsToCoarse.Count();
            double[] ReceiveOffset = new double[MpiSize];
            unsafe
            {
                fixed (double* pCheckSend = SendOffset, pCheckReceive = ReceiveOffset)
                {
                    csMPI.Raw.Allgather((IntPtr)pCheckSend, SendOffset.Length, csMPI.Raw._DATATYPE.DOUBLE, (IntPtr)pCheckReceive, SendOffset.Length, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._COMM.WORLD);
                }
            }
            for (int m = 0; m < MpiRank; m++)
            {
                globalIDOffset += Convert.ToInt32(ReceiveOffset[m]) * 3;
            }

            return globalIDOffset;
        }

        private void CoarseCells(GridCorrelation Old2New, List<Cell> cellsInNewGrid, Cell[][] adaptedCells, int[] cellClusterID)
        {
            Cell[] currentCellClusterToCoarsen = cellClusterID.Select(j => Cells.GetCell(j)).ToArray();
            CheckCoarseningClusterLength(cellClusterID, currentCellClusterToCoarsen);

            int RefinementLevel = currentCellClusterToCoarsen[0].RefinementLevel - 1;
            CheckRefinementLevel(currentCellClusterToCoarsen, RefinementLevel);

            Cell Cell0 = currentCellClusterToCoarsen.Single(cl => cl.ParentCell != null);
            Cell Mother = Cell0.ParentCell;
            CheckMotherCell(currentCellClusterToCoarsen, RefinementLevel, Cell0, Mother);

            Cell restoredCell = RestoreCell(RefinementLevel, Cell0, Mother);
            AdaptBoundaryConditions(Mother, restoredCell);
            AdaptGlobalIDOfRestoredCell(Old2New, adaptedCells, cellClusterID, currentCellClusterToCoarsen, restoredCell);
            cellsInNewGrid.Add(restoredCell);
        }

        private static void CheckCoarseningClusterLength(int[] cellClusterID, Cell[] currentCellClusterToCoarsen)
        {
            if (cellClusterID.Length != currentCellClusterToCoarsen.Length)
                throw new Exception("Error in coarsening algorithm: The no of cells in the coarsening cluster differs from the length of the cluster.");
        }

        private static void CheckRefinementLevel(Cell[] currentCellCluster, int RefinementLevel)
        {
            if (RefinementLevel < 0)
                throw new ArgumentException("Refinement level out of range - corrupted data structure.");
            foreach (Cell cl in currentCellCluster)
            {
                if (cl.RefinementLevel != RefinementLevel + 1)
                    throw new ArgumentException("Refinement varies within refinement cluster - corrupted data structure.");
            }
        }

        private static void CheckMotherCell(Cell[] currentCellCluster, int RefinementLevel, Cell Cell0, Cell Mother)
        {
            if (Mother.Type != Cell0.Type)
                throw new Exception("Error in coarsening algorithm: Mother and child cell are of different types.");
            if (Mother.RefinementLevel != RefinementLevel)
                throw new Exception("Error in coarsening algorithm: Mother cell has a different refinement level.");
            if (currentCellCluster.Where(cl => cl.GlobalID == Mother.GlobalID).Count() != 1)
                throw new Exception("Error in coarsening algorithm: GlobalID of mother cell occurs multiple times in child cells");
        }

        private static Cell RestoreCell(int RefinementLevel, Cell Cell0, Cell Mother)
        {
            Cell restoredCell = new Cell
            {
                Type = Mother.Type,
                NodeIndices = Mother.NodeIndices,
                CoarseningClusterID = Mother.CoarseningClusterID,
                ParentCell = Mother.ParentCell,
                GlobalID = Cell0.GlobalID,
                TransformationParams = Mother.TransformationParams,
                RefinementLevel = RefinementLevel,
                CoarseningClusterSize = Mother.CoarseningClusterSize,
                CoarseningLeafIndex = Mother.CoarseningLeafIndex
            };
            return restoredCell;
        }
        private static void AdaptBoundaryConditions(Cell Mother, Cell restoredCell)
        {
            if (Mother.CellFaceTags != null)
            {
                restoredCell.CellFaceTags = Mother.CellFaceTags.Where(cftag => cftag.EdgeTag > 0 && cftag.EdgeTag < GridCommons.FIRST_PERIODIC_BC_TAG).ToArray();
            }
        }
        private static void AdaptGlobalIDOfRestoredCell(GridCorrelation Old2New, Cell[][] adaptedCells, int[] cellClusterID, Cell[] currentCellCluster, Cell restoredCell)
        {
            for (int iSubDiv = 0; iSubDiv < cellClusterID.Length; iSubDiv++)
            {
                int j = cellClusterID[iSubDiv];
                Cell Cj = currentCellCluster[iSubDiv];
                Debug.Assert(adaptedCells[j] == null);
                adaptedCells[j] = new[] { restoredCell };

                Debug.Assert(Old2New.MappingIndex[j] == null);
                Debug.Assert(Old2New.DestGlobalId[j] == null);
                Old2New.MappingIndex[j] = new int[] { Cj.CoarseningLeafIndex };
                Old2New.DestGlobalId[j] = new long[] { restoredCell.GlobalID };
            }
        }

        private static long GetGlobalIDCounter(GridCommons oldGrid)
        {
            long GlobalIdCounter = oldGrid.NumberOfCells_l;
            Convert.ToDouble(GlobalIdCounter).MPISum();
            Convert.ToInt64(GlobalIdCounter);
            return GlobalIdCounter;
        }

        private void AdaptNeighboursOfChangedCellsToNewGrid(GridCorrelation Old2New, GridCommons oldGrid, BitArray cellsToRefineBitmask, BitArray cellsToCoarseBitmask, List<Cell> cellsInNewGrid, Cell[][] adaptedCells, int j)
        {
            Cell oldCell = oldGrid.Cells[j];
            Cell newCell = oldCell.CloneAs();
            cellsInNewGrid.Add(newCell);
            adaptedCells[j] = new Cell[] { newCell };

            // remove out-dated neighborship info
            if (newCell.CellFaceTags != null && newCell.CellFaceTags.Length > 0)
            {
                int[] oldNeighs = this.Cells.CellNeighbours[j];
                foreach (int jNeigh in oldNeighs)
                {
                    if (cellsToRefineBitmask[jNeigh] || cellsToCoarseBitmask[jNeigh])
                    {
                        long neighbourCellGlobalID = this.Cells.GetGlobalID(jNeigh);

                        for (int i = 0; i < newCell.CellFaceTags.Length; i++)
                        {
                            if (newCell.CellFaceTags[i].NeighCell_GlobalID == neighbourCellGlobalID)
                            {
                                Debug.Assert(newCell.CellFaceTags[i].EdgeTag == 0 || newCell.CellFaceTags[i].EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG);
                                ArrayTools.RemoveAt(ref newCell.CellFaceTags, i);
                                i--;
                            }
                        }
                    }
                }
            }
            AdaptGlobalIDOfUnchangedCells(Old2New, cellsInNewGrid, j);
        }

        private static void CloneAllOtherCellsToNewGrid(GridCorrelation Old2New, GridCommons oldGrid, List<Cell> cellsInNewGrid, int j)
        {
            cellsInNewGrid.Add(oldGrid.Cells[j]);
            AdaptGlobalIDOfUnchangedCells(Old2New, cellsInNewGrid, j);
        }

        private static void AdaptGlobalIDOfUnchangedCells(GridCorrelation Old2New, List<Cell> cellsInNewGrid, int j)
        {
            Debug.Assert(Old2New.MappingIndex[j] == null);
            Debug.Assert(Old2New.DestGlobalId[j] == null);
            Old2New.MappingIndex[j] = null;
            Old2New.DestGlobalId[j] = new long[] { cellsInNewGrid[cellsInNewGrid.Count - 1].GlobalID };
        }

        private void ThrowExceptionsForErrorsInCoarseningCluster(IEnumerable<int[]> CellsToCoarsen, BitArray CellsToRefineBitmask, BitArray CellsToCoarseBitmask, RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves)
        {
            foreach (int[] coarseningCluster in CellsToCoarsen)
            {
                Cell[] coarseningCellCluster = coarseningCluster.Select(j => Cells.GetCell(j)).ToArray();

                int CoarseningClusterID = coarseningCellCluster[0].CoarseningClusterID;
                int iKref = Cells.GetRefElementIndex(coarseningCluster[0]);
                int RefinementLevel = coarseningCellCluster[0].RefinementLevel;

                if (coarseningCluster.Length != KrefS_SubdivLeaves[iKref].Length)
                    throw new ArgumentException("Number of elements in coarsening cluster does not match refinement template for respective element type.");
                if (RefinementLevel <= 0 || CoarseningClusterID <= 0)
                    throw new ArgumentException("Coarsening not available for respective cell.");

                if (coarseningCellCluster.Where(cl => cl.ParentCell != null).Count() != 1)
                    throw new ArgumentException("Coarsening cluster seems wrong, or internal data may be corrupted.");

                for (int z = 0; z < coarseningCellCluster.Length; z++)
                {
                    int j = coarseningCluster[z];

                    if (CellsToRefineBitmask[j] == true)
                        throw new ArgumentException("Cannot refine and coarsen the same cell.");
                    if (CellsToCoarseBitmask[j] == true)
                        throw new ArgumentException("Double entry.", "CellsToCoarsen");
                    CellsToCoarseBitmask[j] = true;

                    Cell Cj = Cells.GetCell(j);
                    if (CoarseningClusterID != Cj.CoarseningClusterID)
                        throw new ArgumentException("Mismatch of 'CoarseningClusterID' within cluster.");
                }
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

        private BitArray GetCellBitMask(IEnumerable<int> cellsToRefine)
        {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            int noOfExternalCells = this.Cells.NoOfExternalCells;
            BitArray cellsToRefineBitMask = new BitArray(noOfLocalCells + noOfExternalCells);

            foreach (int currentCellLocalIndex in cellsToRefine)
            {
                cellsToRefineBitMask[currentCellLocalIndex] = true;
            }

            return cellsToRefineBitMask;
        }

        private BitArray GetCellBitMask(IEnumerable<int[]> cellsToCoarsen)
        {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            int noOfExternalCells = this.Cells.NoOfExternalCells;
            BitArray cellsToRefineBitMask = new BitArray(noOfLocalCells + noOfExternalCells);

            foreach (int[] coarseningCluster in cellsToCoarsen)
            {
                Cell[] coarseningCellCluster = coarseningCluster.Select(j => Cells.GetCell(j)).ToArray();
                for (int z = 0; z < coarseningCellCluster.Length; z++)
                {
                    int currentCellLocalIndex = coarseningCluster[z];
                    cellsToRefineBitMask[currentCellLocalIndex] = true;
                }
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
                Cell[] coarseningCellCluster = coarseningClusterID.Select(j => Cells.GetCell(j)).ToArray();
                for (int z = 0; z < coarseningCellCluster.Length; z++)
                {
                    int currentCellIndex = coarseningClusterID[z];

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

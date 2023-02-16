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

namespace BoSSS.Foundation.Grid.Classic {
    partial class GridData {

        /// <summary>
        /// Creates a new grid, which is an adaptive refinement (cell by cell) of this grid.
        /// </summary>
        /// <param name="cellsToRefine">
        /// All cells to be refined (local indices)
        /// </param>
        /// <param name="cellsToCoarse">
        /// All coarsening clusters with their related cells;
        /// 
        /// </param>
        /// <param name="Old2New">
        /// The correlation between the old and the new grid.
        /// </param>
        public GridCommons Adapt(IEnumerable<int> cellsToRefine, IEnumerable<int[]> cellsToCoarse, out GridCorrelation Old2New) {
            using(new FuncTrace()) {
                void Debug_Assert(bool mustbeTrue, string msg) {
                    //Debug.Assert(mustbeTrue);
                    if(!mustbeTrue)
                        throw new ApplicationException("Internal error in mesh adaptation: " + msg);
                }
                


                // Check for refinement/coarsening on all processes.
                bool anyRefinement = !cellsToRefine.IsNullOrEmpty().MPIAnd();
                bool anyCoarsening = !cellsToCoarse.IsNullOrEmpty().MPIAnd();

                GridCommons oldGrid = m_Grid;
                GridCommons newGrid = new GridCommons(oldGrid.RefElements, oldGrid.EdgeRefElements);
                Partitioning cellPartitioning = CellPartitioning;
                Old2New = new GridCorrelation();

                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);// somehow this barrier is absolutly necessary, I have no idea why.
                int J = Cells.NoOfLocalUpdatedCells;
                int JE = Cells.NoOfExternalCells + Cells.NoOfLocalUpdatedCells;

                BitArray cellsToRefineBitmask = new BitArray(JE);
                BitArray cellsToCoarseBitmask = new BitArray(JE);
                BitArray AdaptNeighborsBitmask = new BitArray(J);

                // templates for subdivision
                // =========================
                RefElement[] KrefS = oldGrid.RefElements; // all ref elements used
                RefElement.SubdivisionTreeNode[] KrefS_SubDiv = new RefElement.SubdivisionTreeNode[KrefS.Length]; // subdivision tree for each ref element
                RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves = new RefElement.SubdivisionTreeNode[KrefS.Length][]; // actual subdivision elements
                (int iSubDiv, int iFace)[][,] KrefS_SubdivConnections = new (int, int)[KrefS.Length][,]; // connections between elements; 1st idx: ref elem; 2nd idx: subdiv elm; 3rd idx: face of subdiv elm; content: [idx of subdiv elm,idx of face]
                int[][][] KrefS_Faces2Subdiv = new int[KrefS.Length][][]; // mapping: [ref elm, face of ref elm] -> Subdivision elements which bound to this face.
                Old2New.GeometricMapping = new AffineTrafo[KrefS.Length][];

                for(int iKref = 0; iKref < KrefS.Length; iKref++) {
                    RefElement Kref = KrefS[iKref];
                    KrefS_SubDiv[iKref] = Kref.GetSubdivisionTree(1);
                    KrefS_SubdivLeaves[iKref] = KrefS_SubDiv[0].GetLeaves();
                    Debug.Assert(ArrayTools.ListEquals(KrefS_SubdivLeaves[iKref], KrefS_SubDiv[iKref].Children[0].GetLevel(), (a, b) => object.ReferenceEquals(a, b)));
                    KrefS_Faces2Subdiv[iKref] = new int[Kref.NoOfFaces][];

                    KrefS_SubdivConnections[iKref] = new (int, int)[KrefS_SubdivLeaves[iKref].Length, KrefS[iKref].NoOfFaces];
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
                Old2New.OldGlobalId = CurrentGlobalIdPermutation.Values.CloneAs();
                Old2New.MappingIndex = new int[J][];
                Old2New.DestGlobalId = new long[J][];

                // define counters to get the correct global id and vertex id
                long noOfGlobalCells = oldGrid.NumberOfCells_l;
                long globalCellIDOffsetLocalProcess = GetGlobalIdOffset(cellsToRefine, KrefS_SubdivLeaves);
                long newVertexCounter = (oldGrid.Cells.Max(cl => cl.NodeIndices.Max()) + 1).MPIMax();

                // all locally adapted cells (on this mpi process)
                Cell[][] adaptedCells = new Cell[J][]; // 1st index: cell index of current mesh; 2nd index: enumeration of subdivisions
                List<Cell> cellsInNewGrid = new List<Cell>();

                // Create and exchange bit-masks
                // =============================
                if(anyRefinement) {
                    CheckForDoubleEntryInEnumeration(cellsToRefine);
                    foreach(int currentCellLocalIndex in cellsToRefine) {
                        cellsToRefineBitmask[currentCellLocalIndex] = true;
                    }
                    GetLocalAndExternalNeighbourCells(cellsToRefine, ref AdaptNeighborsBitmask);
                }
                if(anyCoarsening) {
                    CheckCoarseningCluster(cellsToCoarse, cellsToRefineBitmask, cellsToCoarseBitmask, KrefS_SubdivLeaves);
                    cellsToCoarseBitmask = GetCellBitMask(cellsToCoarse);
                    GetLocalAndExternalNeighbourCells(cellsToCoarse, ref AdaptNeighborsBitmask);
                }
                cellsToRefineBitmask.MPIExchange(this);
                cellsToCoarseBitmask.MPIExchange(this);


                List<Tuple<int, byte, bool>> CellFaceTagsWithPeriodicInverse = new List<Tuple<int, byte, bool>>();
                int noOfEdges = this.Edges.Count;
                for(int j = 0; j < J; j++) {
                    var cellFaceTags = this.Grid.GridData.Cells.GetCell(j).CellFaceTags;
                    if(cellFaceTags.IsNullOrEmpty())
                        continue;
                    for(int c = 0; c < cellFaceTags.Length; c++) {
                        if(cellFaceTags[c].EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG)
                            CellFaceTagsWithPeriodicInverse.Add(new Tuple<int, byte, bool>(j, cellFaceTags[c].EdgeTag, cellFaceTags[c].PeriodicInverse));
                    }
                }

                // Handle unchanged cells
                // ========================= 
                for(int j = 0; j < J; j++) {
                    Debug_Assert(Old2New.OldGlobalId[j] == this.Cells.GetCell(j).GlobalID, "Old to new Gid mismatch");
                    Debug_Assert(Old2New.OldGlobalId[j] == oldGrid.Cells[j].GlobalID, "Old to new Gid mismatch (2)");
                    Debug_Assert(object.ReferenceEquals(this.Cells.GetCell(j), oldGrid.Cells[j]), "missing ref eq");
                    Debug_Assert((cellsToRefineBitmask[j] && cellsToCoarseBitmask[j]) == false, "Cannot refine and coarsen the same cell.");

                    if((cellsToRefineBitmask[j] || cellsToCoarseBitmask[j]) == false) {
                        if(AdaptNeighborsBitmask[j]) {
                            // neighbor information needs to be updated
                            Cell oldCell = oldGrid.Cells[j];
                            Cell newCell = oldCell.CloneAs(); // data 

                            // remove out-dated neighborship info
                            if(newCell.CellFaceTags != null && newCell.CellFaceTags.Length > 0) {
                                int[] oldNeighs = GetNeighboursViaEdgesAndVertices(j);
                                foreach(int jNeigh in oldNeighs) {
                                    if(cellsToRefineBitmask[jNeigh] || cellsToCoarseBitmask[jNeigh]) {
                                        // one of the neighbors has changed, so _potentially_ the cell face tags have to be updated
                                        long gId_Neigh = this.Cells.GetGlobalID(jNeigh);

                                        for(int i = 0; i < newCell.CellFaceTags.Length; i++) {
                                            if(newCell.CellFaceTags[i].NeighCell_GlobalID == gId_Neigh) {
                                                Debug_Assert(newCell.CellFaceTags[i].EdgeTag == 0 || newCell.CellFaceTags[i].EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG, "Illegal face tag");
                                                ArrayTools.RemoveAt(ref newCell.CellFaceTags, i);
                                                i--;
                                            }
                                        }
                                    }
                                }
                            }
                            cellsInNewGrid.Add(newCell);
                            adaptedCells[j] = new Cell[] { newCell };
                        } else {
                            // cell and neighbors remain unchanged
                            var oldCell = oldGrid.Cells[j];
                            var newCell = oldCell.CloneAs(); // data 
                            long gId_Neigh = newCell.GlobalID;
                            int[] oldNeighs = this.Cells.CellNeighbours[j];
                            cellsInNewGrid.Add(newCell);
                            adaptedCells[j] = new Cell[] { newCell };
                        }

                        Debug_Assert(Old2New.MappingIndex[j] == null, "should be null!");
                        Debug_Assert(Old2New.DestGlobalId[j] == null, "should be null (2)!");
                        Old2New.MappingIndex[j] = null;
                        Old2New.DestGlobalId[j] = new long[] { cellsInNewGrid[cellsInNewGrid.Count - 1].GlobalID };

                    } else {
                        Debug_Assert(cellsToRefineBitmask[j] || cellsToCoarseBitmask[j], "either refine or coarsen");
                    }
                }

                // Coarse cells
                // =========================
                if(anyCoarsening) {
                    foreach(int[] cellClusterID in cellsToCoarse) {
                        CoarseCells(Old2New, cellsInNewGrid, adaptedCells, cellClusterID);
                    }
                }

                // Refine cells
                // =========================
                if(anyRefinement) {
                    int NewCoarseningClusterId = GetNewCoarseningClusterID(cellsToRefine, oldGrid);

                    foreach(int j in cellsToRefine) {
                        CheckForDoubleEntries(Old2New, adaptedCells, j);

                        Cell oldCell = oldGrid.Cells[j];
                        int iKref = Cells.GetRefElementIndex(j);
                        RefElement Kref = KrefS[iKref];
                        RefElement.SubdivisionTreeNode[] Leaves = KrefS_SubdivLeaves[iKref];
                        Cell[] refinedCells = new Cell[Leaves.Length];

                        for(int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) {
                            Cell newCell = CreateRefinedCell(ref noOfGlobalCells, globalCellIDOffsetLocalProcess, NewCoarseningClusterId, oldCell, Leaves, iSubDiv);
                            refinedCells[iSubDiv] = newCell;

                            NodeSet RefNodes = Kref.GetInterpolationNodes(oldCell.Type);
                            GetNewNodesOfRefinedCells(j, RefNodes, Leaves, iSubDiv, newCell);
                            newCell.NodeIndices = GetNodeIndicesOfRefinedCells(newVertexCounter, Kref.NoOfVertices, newCell);
                        }
                        NewCoarseningClusterId++;

                        Old2New.MappingIndex[j] = GetOldToNewCorrelation(Leaves);

                        (int, int)[,] Connections = KrefS_SubdivConnections[iKref];
                        AdaptNeighbourshipWithinOldCell(Kref.NoOfFaces, Leaves.Length, Connections, refinedCells);

                        Debug_Assert(Old2New.DestGlobalId[j] == null, "should also be null");
                        Old2New.DestGlobalId[j] = refinedCells.Select(Cl => Cl.GlobalID).ToArray();
                        cellsInNewGrid.AddRange(refinedCells);
                        adaptedCells[j] = refinedCells;
                    }
                }

                // fix neighbor-ship
                // ================
                byte[,] Edge2Face = Edges.FaceIndices;
                int[,] Edge2Cell = Edges.CellIndices;
                MultidimensionalArray[] VerticesFor_KrefEdge = Edges.EdgeRefElements.Select(KrefEdge => KrefEdge.Vertices).ToArray();

                int[] ONE_NULL = new int[] { 0 };

                int NoOfEdges = Edges.Count;
                if(Edge2Face.GetLength(0) != NoOfEdges)
                    throw new Exception("Edge2Face to long");
                if(Edge2Cell.GetLength(0) != NoOfEdges)
                    throw new Exception("Edge2Cell to long");

                // exchange cell data between processes
                List<(long, Cell[])> cellsOnNeighbourProcess = SerialExchangeCellData(adaptedCells);
                for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) { // loop over edges in actual grid...
                    int localCellIndex1 = Edge2Cell[iEdge, 0];
                    int localCellIndex2 = Edge2Cell[iEdge, 1];
                    int iFace1 = Edge2Face[iEdge, 0];
                    int iFace2 = Edge2Face[iEdge, 1];

                    if(localCellIndex2 < 0) { // negative cell index -> boundary cell!
                        if((cellsToRefineBitmask[localCellIndex1] || cellsToCoarseBitmask[localCellIndex1]) == false)
                            continue;
                        AdaptBoundaryCellFaces(adaptedCells, Edge2Face, iEdge, localCellIndex1);
                        continue;
                    }

                    int iKref1 = Cells.GetRefElementIndex(localCellIndex1);
                    int iKref2 = Cells.GetRefElementIndex(localCellIndex2);
                    RefElement Kref1 = Cells.GetRefElement(localCellIndex1);
                    RefElement Kref2 = Cells.GetRefElement(localCellIndex2);

                    Debug_Assert((cellsToRefineBitmask[localCellIndex1] && cellsToCoarseBitmask[localCellIndex1]) == false, "Error in input data: cannot refine and coarsen at the same time (cell1)");
                    Debug_Assert((cellsToRefineBitmask[localCellIndex2] && cellsToCoarseBitmask[localCellIndex2]) == false, "Error in input data: cannot refine and coarsen at the same time (cell2)");

                    bool cell1Changed = cellsToRefineBitmask[localCellIndex1] || cellsToCoarseBitmask[localCellIndex1];
                    bool cell2Changed = cellsToRefineBitmask[localCellIndex2] || cellsToCoarseBitmask[localCellIndex2];

                    if((cell1Changed || cell2Changed) == false)
                        continue;

                    Cell[] adaptedCells1;
                    Cell[] adaptedCells2;

                    bool periodicInverse1 = true;

                    long i0 = CellPartitioning.i0;
                    if(IsPartOfLocalCells(J, localCellIndex1) && IsPartOfLocalCells(J, localCellIndex2)) {
                        // both local
                        adaptedCells1 = adaptedCells[localCellIndex1];
                        adaptedCells2 = adaptedCells[localCellIndex2];
                    } else if(IsPartOfLocalCells(J, localCellIndex1) && !IsPartOfLocalCells(J, localCellIndex2)) {
                        adaptedCells1 = adaptedCells[localCellIndex1];
                        int iKref = Cells.GetRefElementIndex(localCellIndex2);
                        RefElement Kref = KrefS[iKref];
                        RefElement.SubdivisionTreeNode[] Leaves = KrefS_SubdivLeaves[iKref];
                        adaptedCells2 = FindCellOnNeighbourProcess(cellsOnNeighbourProcess, localCellIndex2, Leaves.Length, out int neighbourProcess);

                        periodicInverse1 = neighbourProcess > MpiRank;
                        CheckIfCellIsMissing(localCellIndex2, adaptedCells2, i0);
                    } else if(!IsPartOfLocalCells(J, localCellIndex1) && IsPartOfLocalCells(J, localCellIndex2)) {
                        int iKref = Cells.GetRefElementIndex(localCellIndex1);
                        RefElement Kref = KrefS[iKref];
                        RefElement.SubdivisionTreeNode[] Leaves = KrefS_SubdivLeaves[iKref];
                        adaptedCells1 = FindCellOnNeighbourProcess(cellsOnNeighbourProcess, localCellIndex1, Leaves.Length, out int neighbourProcess);
                        adaptedCells2 = adaptedCells[localCellIndex2];
                        periodicInverse1 = neighbourProcess < MpiRank;
                        CheckIfCellIsMissing(localCellIndex1, adaptedCells1, i0);
                    } else {
                        throw new Exception("Error in refinement and coarsening algorithm: Both cells not on the current process");
                    }
                    bool periodicInverse2 = !periodicInverse1;

                    CheckIfCellIsMissing(localCellIndex1, adaptedCells1, i0);
                    CheckIfCellIsMissing(localCellIndex2, adaptedCells2, i0);

                    if(cellsToCoarseBitmask[localCellIndex1] && cellsToCoarseBitmask[localCellIndex2]) {
                        Debug_Assert(adaptedCells1.Length == 1, "data wrong (3)");
                        Debug_Assert(adaptedCells2.Length == 1, "data wrong (4)");
                        if(adaptedCells1[0] == null)
                            throw new Exception("Cell is missing! Cell with global index " + (i0 + localCellIndex1));
                        if(adaptedCells2[0] == null)
                            throw new Exception("Cell is missing! Cell with global index " + (i0 + localCellIndex2));
                        if(adaptedCells1[0].GlobalID == adaptedCells2[0].GlobalID) {
                            // these two cells will be joint into one cell -> no new neighborship
                            Debug.Assert(ReferenceEquals(adaptedCells1[0], adaptedCells2[0]));
                            continue;
                        }
                    }

                    Debug.Assert((adaptedCells1.Length > 1) == (cellsToRefineBitmask[localCellIndex1]));
                    Debug.Assert((adaptedCells2.Length > 1) == (cellsToRefineBitmask[localCellIndex2]));

                    int[] idx1, idx2; // mapping: Subdivision elements which bound to respective face.
                    if(cellsToRefineBitmask[localCellIndex1]) {
                        idx1 = KrefS_Faces2Subdiv[iKref1][iFace1];
                    } else {
                        Debug_Assert(adaptedCells1.Length == 1, "data wrong (5)");
                        idx1 = ONE_NULL;
                    }

                    if(cellsToRefineBitmask[localCellIndex2]) {
                        idx2 = KrefS_Faces2Subdiv[iKref2][iFace2];
                    } else {
                        Debug_Assert(adaptedCells2.Length == 1, "data wrong (6)");
                        idx2 = ONE_NULL;
                    }

                    foreach(int i1 in idx1) { // compare all candidates in first original cell ...
                        MultidimensionalArray VtxFace1;
                        if(cellsToRefineBitmask[localCellIndex1]) {
                            VtxFace1 = KrefS_SubdivLeaves[iKref1][i1].GetFaceVertices(iFace1);
                        } else {
                            VtxFace1 = Kref1.GetFaceVertices(iFace1);
                        }

                        Cell Cl1 = adaptedCells1[i1];
                        foreach(int i2 in idx2) { // ... with all candidates in second original cell ...

                            Cell Cl2 = adaptedCells2[i2];
                            Debug_Assert(Cl1.GlobalID != Cl2.GlobalID, "The same GlobalID was assigned to to different cells.");

                            int conCount1; // number of connections from Cell1 to Cell2
                            if(Cl1.CellFaceTags == null) {
                                conCount1 = 0;
                            } else {
                                conCount1 = Cl1.CellFaceTags.Where(cfTag => cfTag.NeighCell_GlobalID == Cl2.GlobalID).Count();
                            }
                            Debug_Assert(conCount1 <= 1, "data wrong (8)");
                            //#if DEBUG
/*                            {
                                int conCount2; // number of connections from Cell2 to Cell1
                                if(Cl2.CellFaceTags == null) {
                                    conCount2 = 0;
                                } else {
                                    conCount2 = Cl2.CellFaceTags.Where(cfTag => cfTag.NeighCell_GlobalID == Cl1.GlobalID).Count();
                                }

                                ilPSP.Environment.StdoutOnlyOnRank0 = false;
                                if(conCount1 != conCount2) {
                                    //int D = Cl2.TransformationParams.;
                                    int D = this.SpatialDimension;
                                    Vector centercoordinates = new Vector(D);
                                    int NoOfNodes = Cl2.TransformationParams.NoOfRows;
                                    //Compute Barycenter of rectangular cells
                                    for(int d = 0; d < D; d++) {
                                        double center = 0;
                                        for(int k = 0; k < NoOfNodes; k++) {
                                            center += Cl2.TransformationParams[k, d];
                                        }

                                        centercoordinates[d] = center / ((double)NoOfNodes);
                                    }
                                    Console.Error.WriteLine($"proc{ilPSP.Environment.MPIEnv.MPI_Rank} reporting: coord of {Cl2.GlobalID}: {centercoordinates}");
                                }
                                ilPSP.Environment.StdoutOnlyOnRank0 = true;

                                Debug_Assert(conCount1 == conCount2, $"unsymmetric connectivity: #Con cell1->Cell2: {conCount1}, but vice-versa: {conCount2}; localCellIndex1 = {localCellIndex1}, localCellIndex2 = {localCellIndex2}, J = {J}");
                            }*/
                            //#endif                          

                            if(conCount1 > 0)
                                continue;

                            byte EdgeTag = 0;
                            MultidimensionalArray VtxFace2;
                            {
                                MultidimensionalArray VtxFace2_L;
                                if(cellsToRefineBitmask[localCellIndex2]) {
                                    VtxFace2_L = KrefS_SubdivLeaves[iKref2][i2].GetFaceVertices(iFace2);
                                } else {
                                    VtxFace2_L = Kref2.GetFaceVertices(iFace2);
                                }

                                MultidimensionalArray VtxFace2_G = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                                VtxFace2 = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));

                                this.TransformLocal2Global(VtxFace2_L, VtxFace2_G, localCellIndex2);

                                if(this.Grid.GridData.Edges.EdgeTags[iEdge] >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                                    var perTrf = this.Grid.GridData.Edges.GetPeriodicTrafo(iEdge, false);
                                    MultidimensionalArray VtxFace2_Gtmp = VtxFace2_G.CloneAs();
                                    perTrf.Transform(VtxFace2_Gtmp, VtxFace2_G);
                                    EdgeTag = this.Grid.GridData.Edges.EdgeTags[iEdge];
                                }

                                bool[] Converged = new bool[VtxFace2_L.NoOfRows];
                                this.TransformGlobal2Local(VtxFace2_G, VtxFace2, localCellIndex1, Converged);
                                if(Converged.Any(t => t == false))
                                    throw new ArithmeticException("Newton divergence");
                            }

                            bool bIntersect = false;
                            try {
                                bIntersect = EdgeData.FaceIntersect(VtxFace1, VtxFace2,
                                        Kref1.GetFaceTrafo(iFace1), Kref1.GetInverseFaceTrafo(iFace1),
                                        VerticesFor_KrefEdge,
                                        out bool conformal1, out bool conformal2, out AffineTrafo newTrafo, out int Edg_idx);
                            } catch (Exception e) {
                                var rnd = new Random();
                                
                                using (var stw = new System.IO.StreamWriter("AdaptCrash.rank" + MpiRank + ".txt")) {
                                    stw.WriteLine(e.GetType().Name + ": " + e.Message);
                                    
                                    stw.WriteLine("Face1: " + VtxFace1.NoOfRows + " vertices");
                                    for (int i = 0; i < VtxFace1.NoOfRows; i++) {
                                        stw.WriteLine("#" + i + ": " + VtxFace1.GetRowPt(i).ToString());
                                    }
                                    stw.WriteLine("Face2: " + VtxFace2.NoOfRows + " vertices");
                                    for (int i = 0; i < VtxFace2.NoOfRows; i++) {
                                        stw.WriteLine("#" + i + ": " + VtxFace2.GetRowPt(i).ToString());
                                    }
                                    stw.WriteLine("Face1 Trafo:");
                                    stw.WriteLine(Kref1.GetFaceTrafo(iFace1).ToString());
                                    stw.WriteLine("Face1 Inverse Trafo:");
                                    stw.WriteLine(Kref1.GetInverseFaceTrafo(iFace1).ToString());
                                    stw.Flush();

                                    stw.WriteLine("Cell 1:");
                                    stw.WriteLine(Cl1.ToString());
                                    stw.WriteLine("Cell 2:");
                                    stw.WriteLine(Cl2.ToString());

                                    stw.Close();
                                }

                                throw new AggregateException(e);
                            }
                            bool periodicInverse = false;
                            for(int j = 0; j < CellFaceTagsWithPeriodicInverse.Count(); j++) {
                                if(CellFaceTagsWithPeriodicInverse[j].Item1 == localCellIndex1 && CellFaceTagsWithPeriodicInverse[j].Item2 == EdgeTag) {
                                    periodicInverse = CellFaceTagsWithPeriodicInverse[j].Item3;
                                    break;
                                }
                            }

                            if(bIntersect) {
                                // add connections to the refined cells
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    ConformalNeighborship = false,
                                    NeighCell_GlobalID = Cl2.GlobalID,
                                    FaceIndex = iFace1,
                                    EdgeTag = EdgeTag,
                                    PeriodicInverse = (EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) && periodicInverse
                                }, ref Cl1.CellFaceTags);

                                ArrayTools.AddToArray(new CellFaceTag() {
                                    ConformalNeighborship = false,
                                    NeighCell_GlobalID = Cl1.GlobalID,
                                    FaceIndex = iFace2,
                                    EdgeTag = EdgeTag,
                                    PeriodicInverse = (EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) && (!periodicInverse)
                                }, ref Cl2.CellFaceTags);
                            }
                        }
                    }
                }

                newGrid.Cells = cellsInNewGrid.ToArray();

                AddEdgeTagNames(newGrid);

                AddPeriodicTransformations(newGrid);

                // finalize
                // ========
                if(anyCoarsening) {
#if DEBUG
                    if (MpiSize == 1) {
                        List<int> allgids = new List<int>();
                        foreach (Cell cl in newGrid.Cells) {
                            allgids.Add((int)(cl.GlobalID));
                        }
                        bool[] markers = new bool[allgids.Max() + 1];
                        for (int i = 0; i < allgids.Count; i++) {
                            long gid = allgids[i];
                            Debug.Assert(markers[gid] == false, "Some GlobalID is used twice.");
                            markers[gid] = true;
                        }

                        foreach (Cell cl in newGrid.Cells) {
                            if (cl.CellFaceTags != null) {
                                for (int i = 0; i < cl.CellFaceTags.Length; i++) {
                                    long ngid = cl.CellFaceTags[i].NeighCell_GlobalID;
                                    //if (ngid >= 0)
                                    //Debug.Assert(markers[ngid] == true);
                                }
                            }
                        }
                    }
#endif
                    if(Old2New.DestGlobalId.Length != J)
                        throw new Exception("Error in coarsening algorithm: No of global ID does not fit to no of locally updated cells");

                    List<long> old2NewGlobalId = new List<long>();
                    for(int j = 0; j < J; j++) {
                        old2NewGlobalId.AddRange(Old2New.DestGlobalId[j]);
                    }
                    newGrid.CompressGlobalID(old2NewGlobalId);

                    int c2 = 0;
                    for(int j = 0; j < J; j++) {
                        long[] o2nj = Old2New.DestGlobalId[j];
                        int K = o2nj.Length;
                        for(int k = 0; k < K; k++) {
                            o2nj[k] = old2NewGlobalId[c2];
                            c2++;
                        }
                    }
                    Debug.Assert(c2 == old2NewGlobalId.Count);
                }
                return newGrid;
            }
        }

        // <summary>
        // Calculates the global id offset 
        // </summary>
        // <param name="CellsToRefine">
        // </param>
        //private int GetGlobalIdOffset(IEnumerable<int> cellsToRefine) {
        //    int globalIDOffset = 0;
        //    int noOfRefinedCellsPerCell = 4;
        //    int noOfNewCellsIDPerCell = noOfRefinedCellsPerCell - 1;
        //    int sendOffset = cellsToRefine.Count();
        //    int[] receiveOffset = sendOffset.MPIAllGatherO();
        //    for (int m = 0; m < MpiRank; m++) {
        //        globalIDOffset += receiveOffset[m] * noOfNewCellsIDPerCell;
        //    }
        //    return globalIDOffset;
        //}
        private long GetGlobalIdOffset(IEnumerable<int> CellsToRefine, RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves) {
            long countLeaves = 0;
            for (int j = 0; j < CellsToRefine.Count(); j++) {
                int iKref = Cells.GetRefElementIndex(j);
                countLeaves += KrefS_SubdivLeaves[iKref].Length - 1;
            }
            long sendOffset = countLeaves;
            long[] receiveOffset = sendOffset.MPIAllGatherO();
            long globalIDOffset = 0;
            for (int m = 0; m < MpiRank; m++) {
                globalIDOffset += receiveOffset[m];
            }
            return globalIDOffset;
        }

        /// <summary>
        /// </summary>
        /// <param name="enumeration">
        /// </param>
        private void CheckForDoubleEntryInEnumeration(IEnumerable<int> enumeration) {
            for (int i = 0; i < enumeration.Count(); i++) {
                int currentEntry = enumeration.ElementAt(i);
                for (int j = i + 1; j < enumeration.Count(); j++) {
                    if (currentEntry == enumeration.ElementAt(j))
                        throw new ArgumentException("Double entry.", "CellsToRefine");
                }
            }
        }

        /// <summary>
        /// </summary>
        /// <param name="cellsToRefine">
        /// </param>
        /// <param name="AdaptNeighborsBitmask">
        /// </param>
        private void GetLocalAndExternalNeighbourCells(IEnumerable<int> cellsToRefine, ref BitArray AdaptNeighborsBitmask) {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            long[] externalCellsGlobalIndices = this.iParallel.GlobalIndicesExternalCells;
            List<long>[] exchangeNeighbours = new List<long>[MpiSize];

            foreach (int currentCellIndex in cellsToRefine) {
                int[] neighbourCells = GetNeighboursViaEdgesAndVertices(currentCellIndex);

                foreach (int neighbourCellIndex in neighbourCells) {
                    if (IsPartOfLocalCells(noOfLocalCells, neighbourCellIndex))
                        AdaptNeighborsBitmask[neighbourCellIndex] = true;
                    else {
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

        private int[] GetNeighboursViaEdgesAndVertices(int currentCellIndex) {
            this.GetCellNeighbours(currentCellIndex, GetCellNeighbours_Mode.ViaEdges, out int[] neighbourCellsEdges, out _);
            var edgeNeighbours = this.GetCellNeighboursViaEdges(currentCellIndex);
            List<int> testSomething = new List<int>();
            foreach(var i in edgeNeighbours) {
                testSomething.Add(i.Item1);
            }
            neighbourCellsEdges = testSomething.ToArray();
            this.GetCellNeighbours(currentCellIndex, GetCellNeighbours_Mode.ViaVertices, out int[] neighbourCellsVertices, out _);
            int[] neighbourCells = new int[neighbourCellsEdges.Length + neighbourCellsVertices.Length];
            for (int i = 0; i < neighbourCellsEdges.Length; i++) {
                neighbourCells[i] = neighbourCellsEdges[i];
            }
            for (int i = 0; i < neighbourCellsVertices.Length; i++) {
                neighbourCells[i + neighbourCellsEdges.Length] = neighbourCellsVertices[i];
            }
            return neighbourCells;
        }

        /// <summary>
        /// </summary>
        /// <param name="coarseningClusters">
        /// </param>
        /// <param name="AdaptNeighborsBitmask">
        /// </param>
        private void GetLocalAndExternalNeighbourCells(IEnumerable<int[]> coarseningClusters, ref BitArray AdaptNeighborsBitmask) {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            long[] externalCellsGlobalIndices = this.iParallel.GlobalIndicesExternalCells;
            List<long>[] exchangeNeighbours = new List<long>[MpiSize];

            foreach (int[] coarseningClusterID in coarseningClusters) {
                Cell[] coarseningCellCluster = coarseningClusterID.Select(j => Cells.GetCell(j)).ToArray();
                for (int z = 0; z < coarseningCellCluster.Length; z++) {
                    int currentCellIndex = coarseningClusterID[z];

                    int[] neighbourCells = GetNeighboursViaEdgesAndVertices(currentCellIndex);

                    foreach (int neighbourCellIndex in neighbourCells) {

                        if (IsPartOfLocalCells(noOfLocalCells, neighbourCellIndex))
                            AdaptNeighborsBitmask[neighbourCellIndex] = true;

                        else {
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

        /// <summary>
        /// Checks whether a certain cell is part of the local cells
        /// </summary>
        /// <param name="noOfLocalCells">
        /// </param>
        /// <param name="currentCellIndex">
        /// </param>
        private static bool IsPartOfLocalCells(int noOfLocalCells, int currentCellIndex) {
            return currentCellIndex < noOfLocalCells;
        }

        /// <summary>
        /// Exchange method of the process boundary cells.
        /// </summary>
        /// <param name="exchangeNeighbours">
        /// </param>
        /// <param name="AdaptNeighborsBitmask">
        /// </param>
        private void GetAndExchangeExternalNeighbours(List<long>[] exchangeNeighbours, ref BitArray AdaptNeighborsBitmask) {
            long firstGlobalIndex = this.CellPartitioning.i0;
            List<long> externalNeighbours = new List<long>();

            List<long>[][] tempExchange = exchangeNeighbours.MPIGatherO(0);
            tempExchange = tempExchange.MPIBroadcast(0);
            for (int m = 0; m < MpiSize; m++) {
                if (m != MpiRank) {
                    if (tempExchange[m][MpiRank] != null)
                        externalNeighbours.AddRange(tempExchange[m][MpiRank]);
                }
            }

            for (int j = 0; j < externalNeighbours.Count(); j++) {
                int externalNeighbour = checked((int)(externalNeighbours[j] - firstGlobalIndex));
                AdaptNeighborsBitmask[externalNeighbour] = true;
            }
        }

        /// <summary>
        /// Creates a bit mask from the coarsening clusters
        /// </summary>
        /// <param name="cellsToCoarsen">
        /// </param>
        private BitArray GetCellBitMask(IEnumerable<int[]> cellsToCoarsen) {
            BitArray cellsToRefineBitMask = new BitArray(this.Cells.NoOfLocalUpdatedCells + this.Cells.NoOfExternalCells);

            foreach (int[] coarseningCluster in cellsToCoarsen) {
                Cell[] coarseningCellCluster = coarseningCluster.Select(j => Cells.GetCell(j)).ToArray();
                for (int z = 0; z < coarseningCellCluster.Length; z++) {
                    int currentCellLocalIndex = coarseningCluster[z];
                    cellsToRefineBitMask[currentCellLocalIndex] = true;
                }
            }

            return cellsToRefineBitMask;
        }

        /// <summary>
        /// Gets the new coarsening cluster id 
        /// </summary>
        /// <param name="cellsToRefine">
        /// </param>
        /// <param name="oldGrid">
        /// </param>
        private int GetNewCoarseningClusterID(IEnumerable<int> cellsToRefine, GridCommons oldGrid) {
            int[] locData = new int[] { oldGrid.Cells.Max(cl => cl.CoarseningClusterID), cellsToRefine.Count() + 1 };
            int[] glbData = locData.MPIMax();
            int NewCoarseningClusterId = glbData[0] + 1 + glbData[1] * MpiRank;
            return NewCoarseningClusterId;
        }

        /// <summary>
        /// Checks whether a cell was already refined
        /// </summary>
        /// <param name="Old2New">
        /// </param>
        /// <param name="adaptedCells">
        /// - 1st index: cell index of current mesh
        /// - 2nd index: enumeration of subdivisions
        /// </param>
        /// <param name="j"></param>
        private static void CheckForDoubleEntries(GridCorrelation Old2New, Cell[][] adaptedCells, int j) {
            if (adaptedCells[j] != null)
                throw new Exception("Error in refinement algorithm: Cell was already changend.");
            if (Old2New.MappingIndex[j] != null)
                throw new Exception("Error in refinement algorithm: Mapping index already exists.");
        }

        /// <summary>
        /// Creates a new refined cell
        /// </summary>
        /// <param name="GlobalIdCounter">
        /// </param>
        /// <param name="globalIDOffset">
        /// </param>
        /// <param name="NewCoarseningClusterId"></param>
        /// <param name="oldCell">
        /// </param>
        /// <param name="Leaves">
        /// </param>
        /// <param name="iSubDiv">
        /// </param>
        private static Cell CreateRefinedCell(ref long GlobalIdCounter, long globalIDOffset, int NewCoarseningClusterId, Cell oldCell, RefElement.SubdivisionTreeNode[] Leaves, int iSubDiv) {
            Cell newCell = new Cell {
                Type = oldCell.Type,
                RefinementLevel = oldCell.RefinementLevel + 1,
                CoarseningClusterSize = Leaves.Length,
                CoarseningClusterID = NewCoarseningClusterId,
                CoarseningLeafIndex = iSubDiv
            };
            if (iSubDiv == 0) {
                newCell.GlobalID = oldCell.GlobalID;
                newCell.ParentCell = oldCell.CloneAs();
            }
            else {
                newCell.GlobalID = GlobalIdCounter + globalIDOffset;
                GlobalIdCounter++;
            }
            return newCell;
        }

        /// <summary>
        /// Gets the nodes of the refined cell.
        /// </summary>
        /// <param name="j">
        /// </param>
        /// <param name="RefNodes">
        /// </param>
        /// <param name="Leaves"></param>
        /// <param name="iSubDiv">
        /// </param>
        /// <param name="newCell">
        /// </param>
        private void GetNewNodesOfRefinedCells(int j, NodeSet RefNodes, RefElement.SubdivisionTreeNode[] Leaves, int iSubDiv, Cell newCell) {
            MultidimensionalArray RefNodesRoot = Leaves[iSubDiv].Trafo2Root.Transform(RefNodes);
            newCell.TransformationParams = MultidimensionalArray.Create(RefNodes.Lengths);
            TransformLocal2Global(RefNodesRoot, newCell.TransformationParams, j);
        }

        /// <summary>
        /// Gives the new cell the correct node indices's
        /// </summary>
        /// <param name="newVertexCounter">
        /// </param>
        /// <param name="noOfVertices">
        /// </param>
        /// <param name="newCell"></param>
        private long[] GetNodeIndicesOfRefinedCells(long newVertexCounter, int noOfVertices, Cell newCell) {
            long[] tempNodeIndices = new long[noOfVertices];
            for (int i = 0; i < noOfVertices; i++) {
                tempNodeIndices[i] = newVertexCounter + i + newCell.GlobalID * noOfVertices;
            }
            return tempNodeIndices;
        }

        /// <summary>
        /// Correlates the new cells to the old.
        /// </summary>
        /// <param name="Leaves">
        /// </param>
        private static int[] GetOldToNewCorrelation(RefElement.SubdivisionTreeNode[] Leaves) {
            int[] mappingIndex = new int[Leaves.Length];
            for (int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) {
                mappingIndex[iSubDiv] = iSubDiv;
            }
            return mappingIndex;
        }

        /// <summary>
        /// Updates the neighbourship info between all new cells within a single old cell
        /// </summary>
        /// <param name="noOfFaces">
        /// </param>
        /// <param name="noOfLeaves">
        /// </param>
        /// <param name="Connections">
        /// </param>
        /// <param name="refinedCells">
        /// </param>
        private static void AdaptNeighbourshipWithinOldCell(int noOfFaces, int noOfLeaves, (int, int)[,] Connections, Cell[] refinedCells) {
            for (int iSubDiv = 0; iSubDiv < noOfLeaves; iSubDiv++) {
                for (int iFace = 0; iFace < noOfFaces; iFace++) {
                    int iSubDiv_Neigh = Connections[iSubDiv, iFace].Item1;
                    if (iSubDiv_Neigh >= 0) {
                        ArrayTools.AddToArray(new CellFaceTag() {
                            ConformalNeighborship = true,
                            NeighCell_GlobalID = refinedCells[Connections[iSubDiv, iFace].Item1].GlobalID,
                            FaceIndex = iFace
                        }, ref refinedCells[iSubDiv].CellFaceTags);
                    }
                }
            }
        }

        /// <summary>
        /// Error check of the coaresening cluster.
        /// </summary>
        /// <param name="CellsToCoarsen">
        /// </param>
        /// <param name="CellsToRefineBitmask">
        /// </param>
        /// <param name="CellsToCoarseBitmask">
        /// </param>
        /// <param name="KrefS_SubdivLeaves">
        /// </param>
        private void CheckCoarseningCluster(IEnumerable<int[]> CellsToCoarsen, BitArray CellsToRefineBitmask, BitArray CellsToCoarseBitmask, RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves) {
            foreach (int[] coarseningCluster in CellsToCoarsen) {
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

                for (int z = 0; z < coarseningCellCluster.Length; z++) {
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

        /// <summary>
        /// Restores the coarsed cell from the finer cells in the coarsening cluster.
        /// </summary>
        private void CoarseCells(GridCorrelation Old2New, List<Cell> cellsInNewGrid, Cell[][] adaptedCells, int[] cellClusterID) {
            Cell[] currentCellClusterToCoarsen = cellClusterID.Select(j => Cells.GetCell(j)).ToArray();
            if (cellClusterID.Length != currentCellClusterToCoarsen.Length)
                throw new Exception("Error in coarsening algorithm: The no of cells in the coarsening cluster differs from the length of the cluster.");

            int RefinementLevel = currentCellClusterToCoarsen[0].RefinementLevel - 1;
            if (RefinementLevel < 0)
                throw new ArgumentException("Refinement level out of range - corrupted data structure.");
            foreach (Cell cl in currentCellClusterToCoarsen) {
                if (cl.RefinementLevel != RefinementLevel + 1)
                    throw new ArgumentException("Refinement varies within refinement cluster - corrupted data structure.");
            }

            Cell Cell0 = currentCellClusterToCoarsen.Single(cl => cl.ParentCell != null);
            Cell Mother = Cell0.ParentCell;
            if (Mother.Type != Cell0.Type)
                throw new Exception("Error in coarsening algorithm: Mother and child cell are of different types.");
            if (Mother.RefinementLevel != RefinementLevel)
                throw new Exception("Error in coarsening algorithm: Mother cell has a different refinement level.");
            if (currentCellClusterToCoarsen.Where(cl => cl.GlobalID == Mother.GlobalID).Count() > 1)
                throw new Exception("Error in coarsening algorithm: GlobalID of mother cell occurs multiple times in child cells");

            Cell restoredCell = Mother.CloneAs();
            restoredCell.GlobalID = Cell0.GlobalID;
            restoredCell.RefinementLevel = RefinementLevel;
            if (Mother.CellFaceTags != null)
                restoredCell.CellFaceTags = Mother.CellFaceTags.Where(cftag => cftag.EdgeTag > 0 && cftag.EdgeTag < GridCommons.FIRST_PERIODIC_BC_TAG).ToArray();

            for (int iSubDiv = 0; iSubDiv < cellClusterID.Length; iSubDiv++) {
                int j = cellClusterID[iSubDiv];
                Cell Cj = currentCellClusterToCoarsen[iSubDiv];
                Debug.Assert(adaptedCells[j] == null);
                adaptedCells[j] = new[] { restoredCell };

                Debug.Assert(Old2New.MappingIndex[j] == null);
                Debug.Assert(Old2New.DestGlobalId[j] == null);
                Old2New.MappingIndex[j] = new int[] { Cj.CoarseningLeafIndex };
                Old2New.DestGlobalId[j] = new long[] { restoredCell.GlobalID };
            }
            cellsInNewGrid.Add(restoredCell);
        }

        /// <summary>
        /// Gets the cell data for all process boundary cells on the neighboring process of the current process.
        /// </summary>
        /// <remarks>
        /// Copies the entire cell data to the neighboring process in order to update neighboring information between cells at the process boundary.
        /// </remarks>
        /// - 1st index: cell index of current mesh
        /// - 2nd index: enumeration of subdivisions
        /// <param name="Cells">
        /// All adapted cells, i.e. every cell which was either refined or coarsened.
        /// </param>
        /// <returns>
        /// A collection of refinements/aggregations of cells in the current mesh.
        /// - global index: with respect to current grid, the global index (external range) of the original cell which is refined/coarsened
        /// - cell cluster: refined, resp. coarsened cells in respective cluster (Refinement: multiple entries, Coarsening: a single entry)
        /// </returns>
        private List<(long Gidx, Cell[] AdaptCluster)> SerialExchangeCellData(Cell[][] Cells) {
            List<(long Gidx, Cell[] AdaptCluster)> exchangedCellData = new List<(long Gidx, Cell[] AdaptCluster)>(); 
            Dictionary<int, List<(long, Cell[])>> sendCellData = GetBoundaryCellsAndProcessToSend(Cells);
            IDictionary<int, List<(long, Cell[])>> receiveCellData = SerialisationMessenger.ExchangeData(sendCellData);

            foreach (KeyValuePair<int, List<(long, Cell[])>> kv in receiveCellData) {
                List<(long, Cell[])> list = kv.Value;
                exchangedCellData.AddRange(list);
            }
            return exchangedCellData;
        }

        /// <summary>
        /// Creates a dict. where the first entry is the process to be send to and the second entry is a tuple with the global cell index and the cell data of all process boundary cells of the current process.
        /// </summary>
        /// <param name="Cells">
        /// </param>
        /// <returns>
        /// - key: MPI processor rank (destination)
        /// - value: a set of pairs, containing the global index and the respective cell data
        /// </returns>
        private Dictionary<int, List<(long Gidx, Cell[] AdaptCluster)>> GetBoundaryCellsAndProcessToSend(Cell[][] Cells) {
            Dictionary<int, List<(long, Cell[])>> boundaryCellsProcess = new Dictionary<int, List<(long, Cell[])>>();
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;

            for (int j = 0; j < noOfLocalCells; j++) {
                int[] neighbourCells = GetNeighboursViaEdgesAndVertices(j);

                for (int i = 0; i < neighbourCells.Length; i++) {
                    if (!IsPartOfLocalCells(noOfLocalCells, neighbourCells[i])) {
                        int processorToBeSent = GetExternalCellProcess(neighbourCells[i]);
                        long globalCellIndex = CellPartitioning.i0 + j; // not globalId

                        if (!boundaryCellsProcess.TryGetValue(processorToBeSent, out List<(long, Cell[])> exchangeCellData)) {
                            exchangeCellData = new List<(long, Cell[])>();
                            boundaryCellsProcess.Add(processorToBeSent, exchangeCellData);
                        }

                        bool IsCellAlreadyAdded = false; // there might be repeated neighbour cell info coming from both vertices and edges, so this is to avoid it.
                        foreach ((long, Cell[]) CellsToBeSent in boundaryCellsProcess[processorToBeSent]) {
                            if (CellsToBeSent.Item1 == globalCellIndex) {
                                IsCellAlreadyAdded = true;
                                break;
                            }
                        }

                        if (!IsCellAlreadyAdded)
                            exchangeCellData.Add((globalCellIndex, Cells[j])); // add cell to the list to be exchanged 
                    }
                }
            }
            return boundaryCellsProcess;
        }

        /// <summary>
        /// Gets the process of an external cell.
        /// </summary>
        /// <param name="externalCellIndex">
        /// </param>
        private int GetExternalCellProcess(int externalCellIndex) {
            int noOfLocalCells = this.Cells.NoOfLocalUpdatedCells;
            long[] externalCellsGlobalIndices = iParallel.GlobalIndicesExternalCells;
            int globalCellIndex = (int)externalCellsGlobalIndices[externalCellIndex - noOfLocalCells];
            int process = CellPartitioning.FindProcess(globalCellIndex);
            return process;
        }

        /// <summary>
        /// Searches for a cell on the neighbour process and creates a temporary cell to create the neighbourship info.
        /// </summary>
        /// <param name="cellsOnNeighbourProcess">
        /// </param>
        /// <param name="localCellIndex">
        /// </param>
        /// <param name="leavesLength">
        /// No of cell subdivisions
        /// </param>
        /// <param name="neighbourProcess">
        /// MPI rank of respective cell cells
        /// </param>
        private Cell[] FindCellOnNeighbourProcess(List<(long GlobalIdx, Cell[] AdaptCluster)> cellsOnNeighbourProcess, int localCellIndex, int leavesLength, out int neighbourProcess) {
            Cell[] adaptedCell = null;
            int noOfLocalCells = Cells.NoOfLocalUpdatedCells;
            long[] externalCellsGlobalIndices = iParallel.GlobalIndicesExternalCells;
            long globalIndex = externalCellsGlobalIndices[localCellIndex - noOfLocalCells];
            bool foundNothing = true;
            neighbourProcess = 0;
            
            for (int j = 0; j < cellsOnNeighbourProcess.Count; j++) {
                if (globalIndex == cellsOnNeighbourProcess[j].Item1) {
                    adaptedCell = cellsOnNeighbourProcess[j].Item2;
                    neighbourProcess = this.CellPartitioning.FindProcess(cellsOnNeighbourProcess[j].Item1);
                    foundNothing = false;
                    break;
                }
            }
            if (foundNothing)
                throw new Exception("Found no cell with localCellIndex " + localCellIndex + " global index: " + globalIndex + " on MPIRank " + MpiRank + " no of local cells: " + noOfLocalCells);
            return adaptedCell;
        }

        /// <summary>
        /// Recreates the boundary info on adapted cell faces
        /// </summary>
        private void AdaptBoundaryCellFaces(Cell[][] adaptedCells, byte[,] Edge2Face, int iEdge, int localCellIndex1) {
            Cell[] adaptedBCells1 = adaptedCells[localCellIndex1];
            Debug.Assert(adaptedBCells1 != null);
            int iBFace = Edge2Face[iEdge, 0];
            foreach (Cell cl in adaptedBCells1) {
                if (cl.CellFaceTags.Where(cft => cft.FaceIndex == iBFace).Count() == 0 && this.Edges.EdgeTags[iEdge] > 0) {
                    ArrayTools.AddToArray(new CellFaceTag() {
                        EdgeTag = Edges.EdgeTags[iEdge],
                        ConformalNeighborship = false,
                        NeighCell_GlobalID = long.MinValue,
                        FaceIndex = iBFace
                    }, ref cl.CellFaceTags);
                }
            }
        }

        /// <summary>
        /// Error check for missing cell.
        /// </summary>
        private static void CheckIfCellIsMissing(int jCell, Cell[] adaptedCell, long i0) {
            if (adaptedCell[0] == null) {
                throw new Exception("Cell with global index " + (i0 + jCell) + ", i0 = " + i0 + ", does not exist!");
            }
        }

        private void AddEdgeTagNames(GridCommons newGrid) {
            for (int etCnt = 1; etCnt < this.EdgeTagNames.Count; etCnt++) {
                KeyValuePair<byte, string> etPair = this.EdgeTagNames.ElementAt(etCnt);
                newGrid.EdgeTagNames.Add(etPair);
            }
        }

        private void AddPeriodicTransformations(GridCommons newGrid) {
            foreach (AffineTrafo trafo in this.Grid.PeriodicTrafo) {
                newGrid.PeriodicTrafo.Add(trafo);
            }
            foreach (AffineTrafo itrafo in this.Grid.InversePeriodicTrafo) {
                newGrid.InversePeriodicTrafo.Add(itrafo);
            }
        }

       
    }
}

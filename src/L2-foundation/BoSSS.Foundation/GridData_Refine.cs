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
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Classic {
    partial class GridData {

        /// <summary>
        /// Creates a new grid, which is an adaptive refinement (cell by cell) of this grid.
        /// </summary>
        public GridCommons Adapt(IEnumerable<int> CellsToRefine, int[][] CellsToCoarsen) {
            using(new FuncTrace()) {
                GridCommons oldGrid = this.m_Grid;
                GridCommons newGrid = new GridCommons(oldGrid.RefElements, oldGrid.EdgeRefElements);
                
                int J = this.Cells.NoOfLocalUpdatedCells;

                BitArray CellsToRefineBitmask = new BitArray(J);
                BitArray CellsToCoarseBitmask = new BitArray(J);
                BitArray AdaptNeighborsBitmask = new BitArray(J);
                if(CellsToRefine != null) {
                    foreach(int jCell in CellsToRefine) {
                        if(CellsToRefineBitmask[jCell] == true)
                            throw new ArgumentException("Double entry.", "CellsToRefine");

                        CellsToRefineBitmask[jCell] = true;

                        int[] Neighs, dummy;
                        this.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaEdges, out Neighs, out dummy);

                        foreach(int jNeigh in Neighs) {
                            AdaptNeighborsBitmask[jNeigh] = true;
                        }
                    }
                }

                if(CellsToCoarsen != null) {
                    foreach(int[] jCellS in CellsToCoarsen) {

                        // cluster of cells to coarsen
                        Cell[] CellS = jCellS.Select(j => this.Cells.GetCell(j)).ToArray();

                        for(int z = 0; z < CellS.Length; z++) {
                            int j = jCellS[z];

                            if(CellsToRefineBitmask[j] == true)
                                throw new ArgumentException("Cannot refine and coarsen the same cell.");
                            if(CellsToCoarseBitmask[j] == true)
                                throw new ArgumentException("Double entry.", "CellsToCoarsen");
                            CellsToCoarseBitmask[j] = true;

                            Cell Cj = this.Cells.GetCell(j);
                            if(Cj.CoarseningPeers == null)
                                throw new ArgumentException("Coarsening not available for respective cell.");
                            if(Cj.CoarseningPeers.Length != jCellS.Length - 1)
                                throw new ArgumentException("Coarsening cluster seems incomplete.");

                            foreach(long gid in Cj.CoarseningPeers) {
                                if(CellS.Where(cl => cl.GlobalID == gid).Count() != 1)
                                    throw new ArgumentException("Coarsening cluster seems incomplete.");
                            }

                            int[] Neighs, dummy;
                            this.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out Neighs, out dummy);

                            foreach(int jNeigh in Neighs) {
                                if(Array.IndexOf(jCellS, jNeigh) < 0) {
                                    AdaptNeighborsBitmask[jNeigh] = true;
                                }
                            }
                        }
                    }
                }


                int InsertCounter = J;

                // templates for subdivision
                // =========================

                RefElement[] KrefS = oldGrid.RefElements; // all ref elements used
                RefElement.SubdivisionTreeNode[] KrefS_SubDiv = new RefElement.SubdivisionTreeNode[KrefS.Length]; // subdivision tree for each ref element
                RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves = new RefElement.SubdivisionTreeNode[KrefS.Length][]; // actual subdivision elements
                Tuple<int, int>[][,] KrefS_SubdivConnections = new Tuple<int, int>[KrefS.Length][,]; // connections between elements; 1st idx: ref elem; 2nd idx: subdiv elm; 3rd idx: face of subdiv elm; content: [idx of subdiv elm,idx of face]
                int[][][] KrefS_Faces2Subdiv = new int[KrefS.Length][][]; // mapping: [ref elm, face of ref elm] -> Subdivision elements which bound to this face.
                for(int iKref = 0; iKref < KrefS.Length; iKref++) {
                    RefElement Kref = KrefS[iKref];
                    KrefS_SubDiv[iKref] = Kref.GetSubdivisionTree(1);
                    KrefS_SubdivLeaves[iKref] = KrefS_SubDiv[0].GetLeaves();
                    Debug.Assert(ArrayTools.ListEquals(KrefS_SubdivLeaves[iKref], KrefS_SubDiv[iKref].Children[0].GetLevel(), (a, b) => object.ReferenceEquals(a, b)));
                    KrefS_Faces2Subdiv[iKref] = new int[Kref.NoOfFaces][];

                    KrefS_SubdivConnections[iKref] = new Tuple<int, int>[KrefS[iKref].NoOfFaces, KrefS[iKref].NoOfFaces];
                    for(int iSubdiv = 0; iSubdiv < KrefS_SubdivConnections[iKref].GetLength(0); iSubdiv++) { // loop over subdivision elements
                        for(int iFace = 0; iFace < KrefS_SubdivConnections[iKref].GetLength(0); iFace++) { // loop over faces of subdivision elements
                            var t = KrefS_SubdivLeaves[iKref][iSubdiv].GetNeighbor(iFace);
                            if(t.Item1 < 0) {
                                // at the boundary of the subdivision
                                ArrayTools.AddToArray(iSubdiv, ref KrefS_Faces2Subdiv[iKref][t.Item2]);
                            }
                            KrefS_SubdivConnections[iKref][iSubdiv, iFace] = t;
                        }
                    }
                }

                // create new cells
                // ================

                Debug.Assert(this.MpiSize == 1, "still need to adjust the following lines.");

                long GlobalIdCounter = oldGrid.NumberOfCells_l;
                //int PtrNewCells = oldGrid.NoOfUpdateCells;
                //newGrid.Cells = new Cell[NewNoOfCells];
                List<Cell> newCells = new List<Cell>();
                int newVertexCounter = oldGrid.Cells.Max(cl => cl.NodeIndices.Max()) + 1;
                Cell[][] adaptedCells = new Cell[J][];

                // clone neighbors of refined/coarsened cells
                // ------------------------------------------
                for (int j = 0; j < J; j++) {
                    Debug.Assert((CellsToRefineBitmask[j] && CellsToCoarseBitmask[j]) == false, "Cannot refine and coarsen the same cell.");

                    if((CellsToRefineBitmask[j] || CellsToCoarseBitmask[j]) == false) {
                        if(AdaptNeighborsBitmask[j]) {
                            // neighbor information needs to be updated

                            var oldCell = oldGrid.Cells[j];
                            var newCell = oldCell.CloneAs(); // data 
                            newCells.Add(newCell);
                            adaptedCells[j] = new Cell[] { newCell };
                            
                            // remove out-dated neighborship info
                            if(newCell.CellFaceTags != null && newCell.CellFaceTags.Length > 0) {
                                int[] oldNeighs = this.Cells.CellNeighbours[j];
                                foreach(int jNeigh in oldNeighs) {
                                    if(CellsToRefineBitmask[jNeigh] || CellsToCoarseBitmask[jNeigh]) {
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

                        int RefinementLevel = CellS[0].RefinementLevel - 1;
                        if (RefinementLevel < 0)
                            throw new ArgumentException("Refinement level out of range - corrupted data structure.");
                        foreach(var cl in CellS) {
                            if(cl.RefinementLevel != RefinementLevel+1)
                                throw new ArgumentException("Refinement varies within refinement cluster - corrupted data structure.");
                        }

                        Cell Cell0 = CellS.Single(cl => cl.ParentCell != null);
                        Cell Mother = Cell0.ParentCell;

                       
                        Cell restoredCell = new Cell();
                        restoredCell.Type = Mother.Type;
                        Debug.Assert(Mother.Type == Cell0.Type);
                        restoredCell.NodeIndices = Mother.NodeIndices;
                        restoredCell.CoarseningPeers = Mother.CoarseningPeers;
                        restoredCell.ParentCell = Mother.ParentCell;
                        restoredCell.GlobalID = Mother.GlobalID;
                        restoredCell.TransformationParams = Mother.TransformationParams;
                        restoredCell.RefinementLevel = RefinementLevel;

                        // boundary conditions by cell face tags
                        restoredCell.CellFaceTags = Mother.CellFaceTags.Where(cftag => cftag.EdgeTag > 0 && cftag.EdgeTag < GridCommons.FIRST_PERIODIC_BC_TAG).ToArray();
                        
                        foreach(var j in jCellS) {
                            Debug.Assert(adaptedCells[j] == null);
                            adaptedCells[j] = new[] { restoredCell };
                        }


                        newCells.Add(restoredCell);
                    }
                }

                // refinement
                // ----------
                if(CellsToRefine != null) {
                    foreach(int j in CellsToRefine) {
                        var oldCell = oldGrid.Cells[j];
                        int iKref = this.Cells.GetRefElementIndex(j);
                        var Kref = KrefS[iKref];
                        var Leaves = KrefS_SubdivLeaves[iKref];
                        Tuple<int, int>[,] Connections = KrefS_SubdivConnections[iKref];

                        NodeSet RefNodes = Kref.GetInterpolationNodes(oldCell.Type);

                        Debug.Assert(adaptedCells[j] == null);
                        Cell[] refinedCells = new Cell[Leaves.Length];
                        adaptedCells[j] = refinedCells;
                        for(int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 1: create new cells

                            // create new cell
                            Cell newCell = new Cell();
                            newCell.Type = oldCell.Type;
                            if(iSubDiv == 0) {
                                newCell.GlobalID = oldCell.GlobalID;
                                newCell.ParentCell = oldCell.CloneAs();
                            } else {
                                newCell.GlobalID = GlobalIdCounter;
                                GlobalIdCounter++;
                            }
                            newCell.RefinementLevel = oldCell.RefinementLevel + 1;
                            refinedCells[iSubDiv] = newCell;
                            

                            // Vertices
                            var RefNodesRoot = Leaves[iSubDiv].Trafo2Root.Transform(RefNodes);
                            newCell.TransformationParams = MultidimensionalArray.Create(RefNodes.Lengths);
                            this.TransformLocal2Global(RefNodesRoot, newCell.TransformationParams, j);

                            // node indices
                            newCell.NodeIndices = new int[Kref.NoOfVertices];
                            for(int i = 0; i < Kref.NoOfVertices; i++) {
                                newCell.NodeIndices[i] = newVertexCounter;
                                newVertexCounter++;
                            }
                        }


                        for(int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 2: do other things
                            // record information for (later) coarsening
                            refinedCells[iSubDiv].CoarseningPeers = refinedCells
                                .Where(cell => cell.GlobalID != refinedCells[iSubDiv].GlobalID)
                                .Select(cell => cell.GlobalID)
                                .ToArray();

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
                    }
                }


                newGrid.Cells = newCells.ToArray();

                // fix neighborship
                // ================

                byte[,] Edge2Face = this.Edges.FaceIndices;
                //int[][] Cells2Edges = this.Cells.Cells2Edges;
                int[,] Edge2Cell = this.Edges.CellIndices;
                //byte[] EdgeTags = this.Edges.EdgeTags;
                MultidimensionalArray[] VerticesFor_KrefEdge = this.Edges.EdgeRefElements.Select(KrefEdge => KrefEdge.Vertices).ToArray();

                int[] ONE_NULL = new int[] { 0 };

                int NoOfEdges = this.Edges.Count;
                Debug.Assert(Edge2Face.GetLength(0) == NoOfEdges);
                Debug.Assert(Edge2Cell.GetLength(0) == NoOfEdges);

                for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                    int jCell1 = Edge2Cell[iEdge, 0];
                    int jCell2 = Edge2Cell[iEdge, 1];
                    if(jCell2 < 0)
                        continue;

                    Debug.Assert((CellsToRefineBitmask[jCell1] && CellsToCoarseBitmask[jCell1]) == false);
                    Debug.Assert((CellsToRefineBitmask[jCell2] && CellsToCoarseBitmask[jCell2]) == false);

                    bool C1changed = CellsToRefineBitmask[jCell1] || CellsToCoarseBitmask[jCell1];
                    bool C2changed = CellsToRefineBitmask[jCell2] || CellsToCoarseBitmask[jCell2];
                   
                    if((C1changed || C2changed) == false)
                        // edge between two un-changed cells -- this neighborship remains the same.
                        continue;


                    Cell[] adaptedCells1 = adaptedCells[jCell1];
                    Cell[] adaptedCells2 = adaptedCells[jCell2];

                    if(CellsToCoarseBitmask[jCell1] && CellsToCoarseBitmask[jCell2]) {
                        Debug.Assert(adaptedCells1.Length == 1);
                        Debug.Assert(adaptedCells2.Length == 1);
                        if(adaptedCells1[0].GlobalID == adaptedCells2[0].GlobalID) {
                            // these two cells will be joint into one cell -> no new neighborship
                            Debug.Assert(ReferenceEquals(adaptedCells1[0], adaptedCells2[0]));
                            continue;
                        }

                    }

                    Debug.Assert(adaptedCells1 != null);
                    Debug.Assert(adaptedCells2 != null);

                    int iFace1 = Edge2Face[iEdge, 0];
                    int iFace2 = Edge2Face[iEdge, 1];

                    Debug.Assert((adaptedCells1.Length > 1) == (CellsToRefineBitmask[jCell1]));
                    Debug.Assert((adaptedCells2.Length > 1) == (CellsToRefineBitmask[jCell2]));

                    int iKref1 = this.Cells.GetRefElementIndex(jCell1);
                    int iKref2 = this.Cells.GetRefElementIndex(jCell2);
                    var Kref1 = this.Cells.GetRefElement(jCell1);
                    var Kref2 = this.Cells.GetRefElement(jCell2);

                    int[] idx1, idx2;
                    if(CellsToRefineBitmask[jCell1]) {
                        idx1 = KrefS_Faces2Subdiv[iKref1][iFace1];
                    } else {
                        Debug.Assert(adaptedCells1.Length == 1);
                        idx1 = ONE_NULL;
                    }

                    if(CellsToRefineBitmask[jCell2]) {
                        idx2 = KrefS_Faces2Subdiv[iKref2][iFace2];
                    } else {
                        Debug.Assert(adaptedCells2.Length == 1);
                        idx2 = ONE_NULL;
                    }

                    foreach(int i1 in idx1) {

                        MultidimensionalArray VtxFace1;
                        if(CellsToRefineBitmask[jCell1]) {
                            VtxFace1 = KrefS_SubdivLeaves[iKref1][i1].GetFaceVertices(iFace1);
                        } else {
                            VtxFace1 = Kref1.GetFaceVertices(iFace1);
                        }

                        Cell Cl1 = adaptedCells1[i1];

                        foreach(int i2 in idx2) {

                            Cell Cl2 = adaptedCells2[i2];
                            Debug.Assert(Cl1.GlobalID != Cl2.GlobalID);


                                int conCount1;
                            if(Cl1.CellFaceTags == null) {
                                conCount1 = 0;
                            } else {
                                conCount1 = Cl1.CellFaceTags.Where(cfTag => cfTag.NeighCell_GlobalID == Cl2.GlobalID).Count();
                            }
                            Debug.Assert(conCount1 <= 1);
#if DEBUG
                            int conCount2;
                            if(Cl2.CellFaceTags == null) {
                                conCount2 = 0;
                            } else {
                                conCount2 = Cl2.CellFaceTags.Where(cfTag => cfTag.NeighCell_GlobalID == Cl1.GlobalID).Count();
                            }
                            Debug.Assert(conCount1 == conCount2);
#endif                            
                            if(conCount1 > 0)
                                continue;

                            MultidimensionalArray VtxFace2;
                            {
                                MultidimensionalArray VtxFace2_L;
                                if(CellsToRefineBitmask[jCell2]) {
                                    VtxFace2_L = KrefS_SubdivLeaves[iKref2][i2].GetFaceVertices(iFace2);
                                } else {
                                    VtxFace2_L = Kref2.GetFaceVertices(iFace2);
                                }

                                MultidimensionalArray VtxFace2_G = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                                VtxFace2 = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                                this.TransformLocal2Global(VtxFace2_L, VtxFace2_G, jCell2);
                                bool[] Converged = new bool[VtxFace2_L.NoOfRows];
                                this.TransformGlobal2Local(VtxFace2_G, VtxFace2, jCell1, Converged);
                                if(Converged.Any(t => t == false))
                                    throw new ArithmeticException("Newton divergence");
                            }

                            bool bIntersect = GridData.EdgeData.FaceIntersect(VtxFace1, VtxFace2,
                                Kref1.GetFaceTrafo(iFace1), Kref1.GetInverseFaceTrafo(iFace1),
                                VerticesFor_KrefEdge,
                                out bool conformal1, out bool conformal2, out AffineTrafo newTrafo, out int Edg_idx);

                            if(bIntersect) {

                                ArrayTools.AddToArray(new CellFaceTag() {
                                    ConformalNeighborship = false,
                                    NeighCell_GlobalID = Cl2.GlobalID,
                                    FaceIndex = iFace1
                                }, ref Cl1.CellFaceTags);

                                ArrayTools.AddToArray(new CellFaceTag() {
                                    ConformalNeighborship = false,
                                    NeighCell_GlobalID = Cl1.GlobalID,
                                    FaceIndex = iFace2
                                }, ref Cl2.CellFaceTags);
                            }
                        }
                    }
                }



                /*
                for(int j = 0; j < J; j++) { // loop over all original cells...
                    if(CellsToRefineBitmask[j]) {
                        // +++++++++++++++
                        // cell is refined 
                        // +++++++++++++++
                        int iKref = this.Cells.GetRefElementIndex(j);
                        var Kref = this.Cells.GetRefElement(j);

                        foreach(int i in Cells2Edges[j]) { // loop over old edges 
                                                           // edge index, in and out
                            int iEdge = Math.Abs(i) - 1;
                            int Other = i > 0 ? 1 : 0;
                            int Me = i < 0 ? 1 : 0;
                            Debug.Assert(Edge2Cell[iEdge, Me] == j);

                            // neighbor index
                            int jNeigh = Edge2Cell[iEdge, Other];
                            if(jNeigh >= 0) {
                                int iKrefNeigh = this.Cells.GetRefElementIndex(jNeigh);
                                var KrefNeigh = this.Cells.GetRefElement(jNeigh);

                                // face indices
                                int iFace_Neigh = Edge2Face[iEdge, Other]; // face at neighbor cell
                                int iFace = Edge2Face[iEdge, Me]; // face at cell j
                                
                                // connect all affected cells on 'Me'-side of the edge to original or refined cell on the other side of the edge
                                foreach(int SubdivIdx in KrefS_Faces2Subdiv[iKref][iFace]) { // loop over all affected cells on 'Me'-side...
                                    Cell hC = adaptedCells[j][SubdivIdx]; // cell on 'Me'-side

                                    AddCFT(newGrid, CellsToRefineBitmask, CellsToCoarseBitmask, AdaptNeighborsBitmask, KrefS_SubdivLeaves, KrefS_Faces2Subdiv, adaptedCells, VerticesFor_KrefEdge, 
                                        iEdge, j, hC, jNeigh, SubdivIdx);
                                }
                            } else {
                                // boundary edge: nothing to do
                            }
                        }

                        //for(int iSubdiv = 0; iSubdiv < KrefS_SubdivLeaves[iKref].Length; iSubdiv++) {
                        //    var newCell = refinedOnes[j][iSubdiv];
                        //}

                    } else if(CellsToCoarseBitmask[j]) {
                        // ++++++++++++++++++
                        // cell was coarsened
                        // ++++++++++++++++++

                        Debug.Assert(adaptedCells[j].Length == 1);
                        Cell hC = adaptedCells[j][0];

                        foreach(int i in Cells2Edges[j]) { // loop over old edges 
                                                           // edge index, in and out
                            int iEdge = Math.Abs(i) - 1;
                            int Other = i > 0 ? 1 : 0;
                            int Me = i < 0 ? 1 : 0;
                            Debug.Assert(Edge2Cell[iEdge, Me] == j);

                            // neighbor index
                            int jNeigh = Edge2Cell[iEdge, Other];
                            if(jNeigh >= 0) {
                                if(!(adaptedCells[jNeigh].Length == 1 && object.ReferenceEquals(adaptedCells[jNeigh][0], hC))) { // skip edges that were removed due to coarsening (internal edges of the subdivision cluster)

                                    AddCFT(newGrid, CellsToRefineBitmask, CellsToCoarseBitmask, AdaptNeighborsBitmask, KrefS_SubdivLeaves, KrefS_Faces2Subdiv, adaptedCells, VerticesFor_KrefEdge, 
                                        iEdge, j, hC, jNeigh, int.MinValue);
#if DEBUG
                                    foreach(var cfTag in hC.CellFaceTags) {
                                        long ngid = cfTag.NeighCell_GlobalID;
                                        if(ngid >= 0) {
                                            Debug.Assert(hC.CellFaceTags.Where(cfTag2 => cfTag2.NeighCell_GlobalID == ngid).Count() == 1);
                                        }
                                    }
#endif
                                }
                            } else {
                                // boundary edge: nothing to do
                            }
                        }

                    } else if(AdaptNeighborsBitmask[j]) {
                        // +++++++++++++++++++++++++++++++++++++++++++++
                        // some neighbor cell was refined *or* coarsened
                        // +++++++++++++++++++++++++++++++++++++++++++++

                        Debug.Assert(CellsToRefineBitmask[j] == false); // only for cells which are not refined themselves
                        Debug.Assert(CellsToCoarseBitmask[j] == false); // only for cells which are not refined themselves

                        Debug.Assert(adaptedCells[j].Length == 1);
                        var newCell = adaptedCells[j][0];
                        Debug.Assert(object.ReferenceEquals(newCell, oldGrid.Cells[j]) == false); // we are allowed to change the cell object
                        int[] oldNeighs = this.Cells.CellNeighbours[j];


                        // remove out-dated neighborship info
                        if(newCell.CellFaceTags != null && newCell.CellFaceTags.Length > 0) {
                            foreach(int jNeigh in oldNeighs) {
                                if(CellsToRefineBitmask[jNeigh] || CellsToCoarseBitmask[jNeigh]) {
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
                                                
                        // add new neighborship info
                        foreach(int jNeigh in oldNeighs) {
                            // find edge
                            int iEdge = FindEdge(j, jNeigh, out int Me, out int Other);
                            int iFace = Edge2Face[iEdge, Me]; // face at cell j

                            if(CellsToRefineBitmask[jNeigh]) { // neighbor cell was refined.
                                // we have to add the new neighbors
                                Cell[] refined_Neighs = adaptedCells[jNeigh];
                                Debug.Assert(refined_Neighs != null);
                                int iKrefNeigh = this.Cells.GetRefElementIndex(jNeigh);
                                                                
                                // face of old neighbor cell
                                int iFace_Neigh = Edge2Face[iEdge, Other]; // face at neighbor cell

                                foreach(int iSubdiv_Neigh in KrefS_Faces2Subdiv[iKrefNeigh][iFace_Neigh]) {
                                    Cell newNeigh = refined_Neighs[iSubdiv_Neigh];

                                    // add new cell-face-tag
                                    CellFaceTag cft = default(CellFaceTag);
                                    cft.EdgeTag = EdgeTags[iEdge];
                                    Debug.Assert(cft.EdgeTag == 0 || cft.EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG);
                                    cft.FaceIndex = iFace;
                                    cft.NeighCell_GlobalID = newNeigh.GlobalID;

                                    ArrayTools.AddToArray(cft, ref newCell.CellFaceTags);

                                    TestCellIntersection(newCell, newNeigh);

                                }
                            } else if (CellsToCoarseBitmask[jNeigh]) { // neighbor cell was coarsened.
                                Debug.Assert(adaptedCells[jNeigh].Length == 1);
                                Cell coarsenedNeigh = adaptedCells[jNeigh][0];

                                if(newCell.CellFaceTags.Where(cftag => cftag.NeighCell_GlobalID == coarsenedNeigh.GlobalID && cftag.FaceIndex == iFace).Count() <= 0) {
                                    ArrayTools.AddToArray(
                                        new CellFaceTag() {
                                            NeighCell_GlobalID = coarsenedNeigh.GlobalID,
                                            FaceIndex = iFace
                                        },
                                        ref newCell.CellFaceTags);
                                    TestCellIntersection(newCell, coarsenedNeigh);
                                }
                            }
                        }

                    } else {
                        Debug.Assert(adaptedCells[j] == null);
                        //Debug.Assert(object.ReferenceEquals(newGrid.Cells[j], oldGrid.Cells[j]) == true);
                    }

                    

                }
                */

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
                        foreach(int gid in allgids) {
                            Debug.Assert(markers[gid] == false);
                            markers[gid] = true;
                        }

                        foreach(var cl in newGrid.Cells) {
                            if(cl.CellFaceTags != null) {
                                for(int i = 0; i < cl.CellFaceTags.Length; i++) {
                                    long ngid = cl.CellFaceTags[i].NeighCell_GlobalID;
                                    Debug.Assert(markers[ngid] == true);
                                }
                            }
                        }
                    }
#endif

                    newGrid.CompressGlobalID();

                }

                return newGrid;
            }
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

        /*
        /// <summary>
        /// Subroutine of <see cref="Adapt(int[], int[][])"/>.
        /// </summary>
        private void AddCFT(GridCommons newGrid, 
            BitArray CellsToRefineBitmask, BitArray CellsToCoarseBitmask, BitArray AdaptNeighborsBitmask, 
            RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves, int[][][] KrefS_Faces2Subdiv, Cell[][] adaptedCells, 
            MultidimensionalArray[] VerticesFor_KrefEdge, 
            int iEdge, int j, Cell hC, int jNeigh, int idx) {

            Debug.Assert(j != jNeigh);
            Debug.Assert(adaptedCells[j].IndexOf(hC, (C1, C2) => object.ReferenceEquals(C1, C2)) >= 0);

            int iKref = this.Cells.GetRefElementIndex(j);
            var Kref = this.Cells.GetRefElement(j);
            int iKrefNeigh = this.Cells.GetRefElementIndex(jNeigh);
            var KrefNeigh = this.Cells.GetRefElement(jNeigh);

            // face indices
            byte[,] Edge2Face = this.Edges.FaceIndices;
            int[,] Edge2Cell = this.Edges.CellIndices;
            Debug.Assert((Edge2Cell[iEdge, 0] == j && Edge2Cell[iEdge, 1] == jNeigh) || (Edge2Cell[iEdge, 1] == j && Edge2Cell[iEdge, 0] == jNeigh));
            int Me, Other;
            if (Edge2Cell[iEdge, 0] == j) {
                Me = 0;
                Other = 1;
            } else {
                Me = 1;
                Other = 0;
            }
            Debug.Assert(Edge2Cell[iEdge, Me] == j);
            Debug.Assert(Edge2Cell[iEdge, Other] == jNeigh);


            int iFace = Edge2Face[iEdge, Me]; // face at cell j
            int iFace_Neigh = Edge2Face[iEdge, Other]; // face at neighbor cell



            if (CellsToRefineBitmask[jNeigh]) {
                // + + + + + + + + + + + + + 
                // neighbor cell is refined
                // + + + + + + + + + + + + + 

                Debug.Assert(adaptedCells[jNeigh] != null);
                Debug.Assert(adaptedCells[jNeigh].Length == KrefS_SubdivLeaves[iKrefNeigh].Length);

                foreach(int nidx in KrefS_Faces2Subdiv[iKrefNeigh][iFace_Neigh]) {
                    var pC = adaptedCells[jNeigh][nidx]; // cell on 'Other' side
                    
                    bool bIntersect;
                    MultidimensionalArray VtxFace1;
                    if(CellsToRefineBitmask[j]) {
                        VtxFace1 = KrefS_SubdivLeaves[iKref][idx].GetFaceVertices(iFace);
                    } else {
                        VtxFace1 = Kref.GetFaceVertices(iFace);
                    }

                    MultidimensionalArray VtxFace2;
                    {
                        MultidimensionalArray VtxFace2_L = KrefS_SubdivLeaves[iKrefNeigh][nidx].GetFaceVertices(iFace_Neigh);
                        MultidimensionalArray VtxFace2_G = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                        VtxFace2 = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                        this.TransformLocal2Global(VtxFace2_L, VtxFace2_G, jNeigh);
                        bool[] Converged = new bool[VtxFace2_L.NoOfRows];
                        this.TransformGlobal2Local(VtxFace2_G, VtxFace2, j, Converged);
                        if(Converged.Any(t => t == false))
                            throw new ArithmeticException("Newton divergence");
                    }
                    
                    bIntersect = GridData.EdgeData.FaceIntersect(VtxFace1, VtxFace2,
                        Kref.GetFaceTrafo(iFace), Kref.GetInverseFaceTrafo(iFace),
                        VerticesFor_KrefEdge,
                        out bool conformal1, out bool conformal2, out AffineTrafo newTrafo, out int Edg_idx);

                    if(bIntersect) {
                        ArrayTools.AddToArray(new CellFaceTag() {
                            ConformalNeighborship = false,
                            NeighCell_GlobalID = pC.GlobalID,
                            FaceIndex = iFace
                        }, ref hC.CellFaceTags);

                        TestCellIntersection(hC, pC);
                    }
                }
            } else if(CellsToCoarseBitmask[jNeigh]) {
                // + + + + + + + + + + + + + + 
                // neighbor cell is coarsened
                // + + + + + + + + + + + + + +

                Debug.Assert(adaptedCells[jNeigh] != null);
                Debug.Assert(adaptedCells[jNeigh].Length == 1);

                var pC = adaptedCells[jNeigh][0]; // cell on 'Other' side

                if (hC.CellFaceTags.Where(cfTag => cfTag.NeighCell_GlobalID == pC.GlobalID && cfTag.FaceIndex == iFace).Count() <= 0) {
                    ArrayTools.AddToArray(new CellFaceTag() {
                        NeighCell_GlobalID = pC.GlobalID,
                        FaceIndex = iFace
                    }, ref hC.CellFaceTags);
                    TestCellIntersection(hC, pC);
                }

            } else {
                // + + + + + + + + + + + + + + + +
                // neighbor cell is *not* changed
                // + + + + + + + + + + + + + + + +

                Debug.Assert(adaptedCells[jNeigh] != null);
                Debug.Assert(adaptedCells[jNeigh].Length == 1);
                Debug.Assert(AdaptNeighborsBitmask[jNeigh] == true);
                Debug.Assert(CellsToRefineBitmask[jNeigh] == false);
                Debug.Assert(CellsToCoarseBitmask[jNeigh] == false);

                var pC = adaptedCells[jNeigh][0];

                if(hC.CellFaceTags.Where(cfTag => cfTag.NeighCell_GlobalID == pC.GlobalID && cfTag.FaceIndex == iFace).Count() <= 0) {

                    bool bIntersect;
                    MultidimensionalArray VtxFace1;
                    if(CellsToRefineBitmask[j]) {
                        VtxFace1 = KrefS_SubdivLeaves[iKref][idx].GetFaceVertices(iFace);
                    } else {
                        VtxFace1 = Kref.GetFaceVertices(iFace);
                    }

                    MultidimensionalArray VtxFace2;
                    {
                        MultidimensionalArray VtxFace2_L = KrefNeigh.GetFaceVertices(iFace_Neigh);
                        MultidimensionalArray VtxFace2_G = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                        VtxFace2 = MultidimensionalArray.Create(VtxFace2_L.GetLength(0), VtxFace2_L.GetLength(1));
                        this.TransformLocal2Global(VtxFace2_L, VtxFace2_G, jNeigh);
                        bool[] Converged = new bool[VtxFace2_L.NoOfRows];
                        this.TransformGlobal2Local(VtxFace2_G, VtxFace2, j, Converged);
                        if(Converged.Any(t => t == false))
                            throw new ArithmeticException("Newton divergence");
                    }

                    bIntersect = GridData.EdgeData.FaceIntersect(VtxFace1, VtxFace2,
                        Kref.GetFaceTrafo(iFace), Kref.GetInverseFaceTrafo(iFace),
                        VerticesFor_KrefEdge,
                        out bool conformal1, out bool conformal2, out AffineTrafo newTrafo, out int Edg_idx);

                    if(bIntersect) {
                        ArrayTools.AddToArray(new CellFaceTag() {
                            ConformalNeighborship = false,
                            NeighCell_GlobalID = pC.GlobalID,
                            FaceIndex = iFace
                        }, ref hC.CellFaceTags);
                        TestCellIntersection(hC, pC);
                    }
                }
            }
        }

        /// <summary>
        /// Finds the edge between cell <paramref name="j"/> and cell <paramref name="jNeigh"/>.
        /// </summary>
        /// <param name="j">First cell index.</param>
        /// <param name="jNeigh">Second cell index.</param>
        /// <param name="Me">
        /// 0, if cell <paramref name="j"/> is the in-cell of the respective edge; 1 if it is the out-cell.
        /// </param>
        /// <param name="Other"></param>
        /// <returns>
        /// Edge index.
        /// </returns>
        private int FindEdge(int j, int jNeigh, out int Me, out int Other) {
            int iEdge = -1;
            Other = -1;
            Me = -1;

            int[] C2E_j = this.Cells.Cells2Edges[j];
            int[,] Edge2Cell = this.Edges.CellIndices;

            foreach(int i in C2E_j) {
                Debug.Assert(i != 0);
                int _iEdge = Math.Abs(i) - 1;
                int _Other = i > 0 ? 1 : 0;
                int _Me = i < 0 ? 1 : 0;
                Debug.Assert(Edge2Cell[_iEdge, _Me] == j);
                if(Edge2Cell[_iEdge, _Other] == jNeigh) {
                    Me = _Me;
                    Other = _Other;
                    Debug.Assert(iEdge < 0); // Edge found twice
                    iEdge = _iEdge;
                }
            }
            Debug.Assert(iEdge >= 0); // unable to find edge
            Debug.Assert(Edge2Cell[iEdge, Me] == j);
            Debug.Assert(Edge2Cell[iEdge, Other] == jNeigh);

            return iEdge;
        }
        */
    }
}

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
using ilPSP;
using ilPSP.Utils;
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
        public GridCommons Adapt(int[] CellsToRefine, int[][] CellsToCoarsen) {
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
            Tuple<int, int>[][,] KrefS_SubdivConnections = new Tuple<int, int>[KrefS.Length][,]; // connections between elements; 1st idx: ref elem; 2nd idx: subdiv elm; 3rd idx: face of subdiv elm; content: idx of subdiv elm
            int[][][] KrefS_Faces2Subdiv = new int[KrefS.Length][][]; // mapping: [ref elm, face of ref elm] -> Subdivision elements which bound to this face.
            for(int iKref = 0; iKref < KrefS.Length; iKref++) {
                RefElement Kref = KrefS[iKref];
                KrefS_SubDiv[iKref] = Kref.GetSubdivisionTree(1);
                KrefS_SubdivLeaves[iKref] = KrefS_SubDiv[0].GetLeaves();
                Debug.Assert(ArrayTools.ListEquals(KrefS_SubdivLeaves[iKref], KrefS_SubDiv[iKref].Children[0].GetLevel(), (a, b) => object.ReferenceEquals(a, b)));
                KrefS_Faces2Subdiv[iKref] = new int[Kref.NoOfFaces][];

                KrefS_SubdivConnections[iKref] = new Tuple<int, int>[KrefS[iKref].NoOfFaces, KrefS[iKref].NoOfFaces];
                for (int iSubdiv = 0; iSubdiv < KrefS_SubdivConnections[iKref].GetLength(0); iSubdiv++) { // loop over subdivision elements
                    for (int iFace = 0; iFace < KrefS_SubdivConnections[iKref].GetLength(0); iFace++) { // loop over faces of subdivision elements
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
            int PtrNewCells = oldGrid.NoOfUpdateCells;
            //newGrid.Cells = new Cell[NewNoOfCells];
            List<Cell> newCells = new List<Cell>();

            int newVertexCounter = oldGrid.Cells.Max(cl => cl.NodeIndices.Max()) + 1;
            Cell[][] adaptedCells = new Cell[J][];
            for(int j = 0; j < J; j++) {
                Debug.Assert(CellsToCoarseBitmask[j] && CellsToCoarseBitmask[j] == false, "Cannot refine and coarsen the same cell.");

                if(CellsToRefineBitmask[j] && CellsToCoarseBitmask[j] == false) {
                    if(AdaptNeighborsBitmask[j]) {
                        // neighbor information needs to be updated

                        var oldCell = oldGrid.Cells[j];
                        var newCell = oldCell.CloneAs(); // data 
                        newCells.Add(newCell);
                        adaptedCells[j] = new Cell[] { newCell };

                    } else {
                        // cell and neighbors remain unchanged
                        newCells.Add(oldGrid.Cells[j]);
                    }
                }
            }

            

            // coarsening
            // ----------
            if(CellsToCoarsen != null) {
                foreach(int[] jCellS in CellsToCoarsen) {

                    // cluster of cells to coarsen
                    Cell[] CellS = jCellS.Select(j => this.Cells.GetCell(j)).ToArray();

                    Cell Cell0 = CellS.Single(cl => cl.ParentCell != null);
                    Cell Mother = Cell0.ParentCell;

                    /*
                    Mother.CellFaceTags
                    Mother.CoarseningPeers
                    Mother.NodeIndices
                    Mother.ParentCell
                    Mother.Type
                    Mother.GlobalID
                    Mother.TransformationParams
                    */

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

                    Cell[] refinedCells = new Cell[Leaves.Length];
                    adaptedCells[j] = refinedCells;
                    for(int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 1: create new cells

                        // create new cell
                        Cell newCell = new Cell();
                        newCell.Type = oldCell.Type;
                        if(iSubDiv == 0) {
                            newCell.GlobalID = oldCell.GlobalID;
                            newGrid.Cells[j] = newCell;
                            newCell.ParentCell = oldCell.CloneAs();
                        } else {
                            newCell.GlobalID = GlobalIdCounter;
                            GlobalIdCounter++;

                            newGrid.Cells[PtrNewCells] = newCell;
                            PtrNewCells++;
                        }
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
            int[][] Cells2Edges = this.Cells.Cells2Edges;
            int[,] Edge2Cell = this.Edges.CellIndices;
            byte[] EdgeTags = this.Edges.EdgeTags;
            MultidimensionalArray[] VerticesFor_KrefEdge = this.Edges.EdgeRefElements.Select(KrefEdge => KrefEdge.Vertices).ToArray();


            for(int j = 0; j < J; j++) { // loop over all original cells...
                if (CellsToRefineBitmask[j]) { 
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
                        if (jNeigh >= 0) {
                            int iKrefNeigh = this.Cells.GetRefElementIndex(jNeigh);
                            var KrefNeigh = this.Cells.GetRefElement(jNeigh);

                            // face indices
                            int iFace_Neigh = Edge2Face[iEdge, Other]; // face at neighbor cell
                            int iFace = Edge2Face[iEdge, Me]; // face at cell j


                            // connect all affected cells on 'Me'-side of the edge to original or refined cell on the other side of the edge
                            foreach(int idx in KrefS_Faces2Subdiv[iKref][iFace]) { // loop over all affected cells on 'Me'-side...
                                Cell hC = adaptedCells[j][idx]; // cell on 'Me'-side
                                
                                if( CellsToRefineBitmask[jNeigh]) {
                                    // +++++++++++++++++++++++++++++
                                    // neighbor cell is also refined
                                    // +++++++++++++++++++++++++++++

                                    foreach(int nidx in KrefS_Faces2Subdiv[iKrefNeigh][iFace_Neigh]) {
                                        var pC = adaptedCells[jNeigh][nidx]; // cell on 'Other' side

                                        MultidimensionalArray VtxFace1 = KrefS_SubdivLeaves[iKref][idx].GetFaceVertices(iFace);

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

                                        
                                        bool bIntersect = GridData.EdgeData.FaceIntersect(VtxFace1, VtxFace2,
                                            Kref.GetFaceTrafo(iFace), Kref.GetInverseFaceTrafo(iFace),
                                            VerticesFor_KrefEdge,
                                            out bool conformal1, out bool conformal2, out AffineTrafo newTrafo, out int Edg_idx);
                                            
                                        if(bIntersect) {
                                            ArrayTools.AddToArray(new CellFaceTag() {
                                                ConformalNeighborship = false,
                                                NeighCell_GlobalID = pC.GlobalID,
                                                FaceIndex = iFace
                                            }, ref hC.CellFaceTags);
                                        }
                                    }
                                } else {
                                    // ++++++++++++++++++++++++++++++
                                    // neighbor cell is *not* refined
                                    // ++++++++++++++++++++++++++++++
                                    ArrayTools.AddToArray(new CellFaceTag() {
                                        ConformalNeighborship = false,
                                        NeighCell_GlobalID = newGrid.Cells[jNeigh].GlobalID,
                                        FaceIndex = iFace
                                    }, ref hC.CellFaceTags);
                                }
                            }
                        } else {

                        }
                    }

                    //for(int iSubdiv = 0; iSubdiv < KrefS_SubdivLeaves[iKref].Length; iSubdiv++) {
                    //    var newCell = refinedOnes[j][iSubdiv];
                    //}

                } else if (AdaptNeighborsBitmask[j]) {
                    // ++++++++++++++++++++++++++
                    // neighbor cell was refined
                    // ++++++++++++++++++++++++++

                    Debug.Assert(CellsToRefineBitmask[j] == false); // only for cells which are not refined themselves
                    Debug.Assert(object.ReferenceEquals(newGrid.Cells[j], oldGrid.Cells[j]) == false); // we are allowed to change the cell object
                    var newCell = newGrid.Cells[j];
                    int[] oldNeighs = this.Cells.CellNeighbours[j];


                    // remove out-dated neighborship info
                    if (newCell.CellFaceTags != null && newCell.CellFaceTags.Length > 0) {
                        foreach (int jNeigh in oldNeighs) {
                            if (CellsToRefineBitmask[jNeigh]) {
                                // one of the neighbors has changed, so _potentially_ the cell face tags have to be updated
                                long gId_Neigh = this.Cells.GetGlobalID(jNeigh);
                                                                
                                for (int i = 0; i < newCell.CellFaceTags.Length; i++) {
                                    if (newCell.CellFaceTags[i].NeighCell_GlobalID == gId_Neigh) {
                                        Debug.Assert(newCell.CellFaceTags[i].EdgeTag == 0 || newCell.CellFaceTags[i].EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG);
                                        ArrayTools.RemoveAt(ref newCell.CellFaceTags, i);
                                        i--;
                                    }
                                }
                            }
                        }
                    }

                    // add new neighborship info
                    int[] C2E_j = Cells2Edges[j];
                    foreach (int jNeigh in oldNeighs) {
                        if (CellsToRefineBitmask[jNeigh]) {
                            // we have to add the new neighbors
                            Cell[] refined_Neighs = adaptedCells[jNeigh];
                            Debug.Assert(refined_Neighs != null);
                            int iKrefNeigh = this.Cells.GetRefElementIndex(jNeigh);


                            // find edge
                            int iEdge = -1;
                            int Other = -1, Me = -1;
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


                            // know face of old neighbor cell
                            int iFace_Neigh = Edge2Face[iEdge, Other]; // face at neighbor cell
                            int iFace = Edge2Face[iEdge, Me]; // face at cell j

                            foreach (int iSubdiv_Neigh in KrefS_Faces2Subdiv[iKrefNeigh][iFace_Neigh]) {
                                Cell newNeigh = refined_Neighs[iSubdiv_Neigh];

                                // add new cell-face-tag
                                CellFaceTag cft = default(CellFaceTag);
                                cft.EdgeTag = EdgeTags[iEdge];
                                Debug.Assert(cft.EdgeTag == 0 || cft.EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG);
                                cft.FaceIndex = iFace;
                                cft.NeighCell_GlobalID = newNeigh.GlobalID;

                                ArrayTools.AddToArray(cft, ref newCell.CellFaceTags);
                            }
                        }
                    }


                } else {
                    Debug.Assert(object.ReferenceEquals(newGrid.Cells[j], oldGrid.Cells[j]) == true);
                }
            }



            return newGrid;
        }



         



    }
}

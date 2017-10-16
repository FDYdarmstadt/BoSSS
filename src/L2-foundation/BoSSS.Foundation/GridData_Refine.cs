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
        public GridCommons Refine() {
            GridCommons oldGrid = this.m_Grid;
            GridCommons newGrid = new GridCommons(oldGrid.RefElements, oldGrid.EdgeRefElements);

            int[] CellsToRefine = new[] { 4 };


            int NewNoOfCells = CellsToRefine.Length * 3 + this.Cells.NoOfLocalUpdatedCells;
            int J = this.Cells.NoOfLocalUpdatedCells;

            BitArray CellsToRefineBitmask = new BitArray(J);
            BitArray RefineNeighborsBitmask = new BitArray(J);
            foreach(int jCell in CellsToRefine) {
                CellsToRefineBitmask[jCell] = true;

                int[] Neighs, dummy;
                this.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaEdges, out Neighs, out dummy);

                foreach(int jNeigh in Neighs) {
                    RefineNeighborsBitmask[jNeigh] = true;
                }
            }


            int InsertCounter = J;

            // templates for subdivision
            // =========================

            RefElement[] KrefS = oldGrid.RefElements; // all ref elements used
            RefElement.SubdivisionTreeNode[] KrefS_SubDiv = new RefElement.SubdivisionTreeNode[KrefS.Length]; // subdivision tree for each ref element
            RefElement.SubdivisionTreeNode[][] KrefS_SubdivLeaves = new RefElement.SubdivisionTreeNode[KrefS.Length][]; // actual subdivision elements
            Tuple<int, int>[][,] KrefS_SubdivConnections = new Tuple<int, int>[KrefS.Length][,]; // connections between elements; 1st idx: ref elem; 2nd idx: subdiv elm; 3rd idx: face of subdiv elm; content: idx od subdiv elm
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
            newGrid.Cells = new Cell[NewNoOfCells];
            int newVertexCounter = oldGrid.Cells.Max(cl => cl.NodeIndices.Max()) + 1;
            Cell[][] refinedOnes = new Cell[J][];
            for(int j = 0; j < J; j++) {
                if (CellsToRefineBitmask[j]) {
                    var oldCell = oldGrid.Cells[j];
                    int iKref = this.Cells.GetRefElementIndex(j);
                    var Kref = KrefS[iKref];
                    var Leaves = KrefS_SubdivLeaves[iKref];
                    Tuple<int, int>[,] Connections = KrefS_SubdivConnections[iKref];

                    NodeSet RefNodes = Kref.GetInterpolationNodes(oldCell.Type);

                    Cell[] newCellS = new Cell[Leaves.Length];
                    refinedOnes[j] = newCellS;
                    for (int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 1: create new cells

                        // create new cell
                        Cell newCell = new Cell();
                        newCell.Type = oldCell.Type;
                        if (iSubDiv == 0) {
                            newCell.GlobalID = oldCell.GlobalID;
                            newGrid.Cells[j] = newCell;
                        } else {
                            newCell.GlobalID = GlobalIdCounter;
                            GlobalIdCounter++;

                            newGrid.Cells[PtrNewCells] = newCell;
                            PtrNewCells++;
                        }
                        newCellS[iSubDiv] = newCell;

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
                    for (int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 2: do other things
                        // neighbors within
                        for (int iFace = 0; iFace < Kref.NoOfFaces; iFace++) {
                            int iSubDiv_Neigh = Connections[iSubDiv, iFace].Item1;
                            if (iSubDiv_Neigh >= 0) {
                                CellFaceTag cft;
                                cft.ConformalNeighborship = true;
                                cft.NeighCell_GlobalID = newCellS[Connections[iSubDiv, iFace].Item1].GlobalID;
                                cft.FaceIndex = iFace;
                            }
                        }
                    }
                } else if(RefineNeighborsBitmask[j]) {
                    // neighbor information needs to be updated

                    var oldCell = oldGrid.Cells[j];
                    var newCell = oldCell.CloneAs(); // data 
                    newGrid.Cells[j] = newCell;
                    
                } else {
                    // cell and neighbors remain unchanged
                    newGrid.Cells[j] = oldGrid.Cells[j];
                }
            }

            // fix neighborship
            // ================

            byte[,] Edge2Face = this.Edges.FaceIndices;
            int[][] Cells2Edges = this.Cells.Cells2Edges;
            int[,] Edge2Cell = this.Edges.CellIndices;
            byte[] EdgeTags = this.Edges.EdgeTags;

            for (int j = 0; j < J; j++) {
                if (CellsToRefineBitmask[j]) {
                    int iKref = this.Cells.GetRefElementIndex(j);

                    foreach(int i in Cells2Edges[j]) { // loop over old edges 
                        // edge index, in and out
                        int iEdge = Math.Abs(i) - 1;
                        int Other = i > 0 ? 1 : 0;
                        int Me = i < 0 ? 1 : 0;
                        Debug.Assert(Edge2Cell[iEdge, Me] == j);

                        // neighbor index
                        int jNeigh = Edge2Cell[iEdge, Other];
                        int iKrefNeigh = this.Cells.GetRefElementIndex(jNeigh);

                        // face indices
                        int iFace_Neigh = Edge2Face[iEdge, Other]; // face at neighbor cell
                        int iFace = Edge2Face[iEdge, Me]; // face at cell j

                        // affected cells on 'Me'-side of the edge
                        IEnumerable<Cell> hereCells;
                        hereCells = KrefS_Faces2Subdiv[iKref][iFace].Select(idx => refinedOnes[j][idx]);

                        // peer cells
                        IEnumerable<Cell> peerCells;
                        if(CellsToRefineBitmask[jNeigh]) {
                            peerCells = KrefS_Faces2Subdiv[iKrefNeigh][iFace_Neigh].Select(idx => refinedOnes[jNeigh][idx]);
                        } else {
                            peerCells = new[] { newGrid.Cells[jNeigh] };
                        }

                        // connect all 'hereCells' to the 'peerCells'
                        todo

                    }


                    for(int iSubdiv = 0; iSubdiv < KrefS_SubdivLeaves[iKref].Length; iSubdiv++) {
                        var newCell = refinedOnes[j][iSubdiv];




                    }

                } else if (RefineNeighborsBitmask[j]) {
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
                            Cell[] refined_Neighs = refinedOnes[jNeigh];
                            Debug.Assert(refined_Neighs != null);
                            int iKrefNeigh = this.Cells.GetRefElementIndex(jNeigh);


                            // find edge
                            int iEdge = -1;
                            int Other = -1, Me = -1;
                            foreach(int i in C2E_j) {
                                Debug.Assert(i != 0);
                                int _iEdge = Math.Abs(i) - 1;
                                Other = i > 0 ? 1 : 0;
                                Me = i < 0 ? 1 : 0;
                                Debug.Assert(Edge2Cell[_iEdge, Me] == j);
                                if(Edge2Cell[_iEdge, Other] == jNeigh) {
                                    Debug.Assert(iEdge < 0); // Edge found twice
                                    iEdge = _iEdge;
                                }
                            }
                            Debug.Assert(iEdge >= 0); // unable to find edge

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

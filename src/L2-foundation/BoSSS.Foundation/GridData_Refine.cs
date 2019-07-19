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
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Classic {
    partial class GridData {

        



        /// <summary>
        /// Creates a new grid, which is an adaptive refinement (cell by cell) of this grid.
        /// </summary>
        public GridCommons Adapt(IEnumerable<int> CellsToRefine, IEnumerable<int[]> CellsToCoarsen, out GridCorrelation Old2New) {
            using(new FuncTrace()) {
                GridCommons oldGrid = this.m_Grid;
                GridCommons newGrid = new GridCommons(oldGrid.RefElements, oldGrid.EdgeRefElements);


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
                //List<AffineTrafo> InterCellTrafos = new List<AffineTrafo>();
                //int[][] RefinementIctIdx = new int[KrefS.Length][];
                //int[][] CoarseningIctIdx = new int[KrefS.Length][];

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

                    //RefinementIctIdx[iKref] = new int[KrefS_SubdivLeaves[iKref].Length];
                    //CoarseningIctIdx[iKref] = new int[KrefS_SubdivLeaves[iKref].Length];
                    Old2New.GeometricMapping[iKref] = new AffineTrafo[KrefS_SubdivLeaves[iKref].Length];
                    for(int iSubDiv = 0; iSubDiv < KrefS_SubdivLeaves[iKref].Length; iSubDiv++) {
                        Old2New.GeometricMapping[iKref][iSubDiv] = KrefS_SubdivLeaves[iKref][iSubDiv].TrafoFromRoot;
                        //InterCellTrafos.Add(KrefS_SubdivLeaves[iKref][iSubDiv].TrafoFromRoot);
                        //RefinementIctIdx[iKref][iSubDiv] = InterCellTrafos.Count - 1;
                        //InterCellTrafos.Add(KrefS_SubdivLeaves[iKref][iSubDiv].Trafo2Root);
                        //CoarseningIctIdx[iKref][iSubDiv] = InterCellTrafos.Count - 1;
                    }
                }

                Old2New.KrefS_SubdivLeaves = KrefS_SubdivLeaves;

                BitArray AdaptNeighboursOtherProcess = new BitArray(JE);

                // Check Input, set Bitmasks
                // =========================

                if (CellsToRefine != null) {
                    foreach(int jCell in CellsToRefine) {
                        if(CellsToRefineBitmask[jCell] == true)
                            throw new ArgumentException("Double entry.", "CellsToRefine");

                        CellsToRefineBitmask[jCell] = true;

                        int[] Neighs, dummy;
                        this.GetCellNeighbours(jCell, GetCellNeighbours_Mode.ViaEdges, out Neighs, out dummy);

                        foreach(int jNeigh in Neighs) {
                            if (jNeigh < J)
                                AdaptNeighborsBitmask[jNeigh] = true;
                            else
                                AdaptNeighboursOtherProcess[jNeigh] = true;
                        }
                    }
                }
                AdaptNeighboursOtherProcess.MPIExchange(this);

                for (int j = 0; j < J; j++)
                {
                    if (AdaptNeighboursOtherProcess[j])
                        AdaptNeighborsBitmask[j] = true;
                }

                if (CellsToCoarsen != null) {
                    foreach(int[] jCellS in CellsToCoarsen) { // loop over all coarsening clusters...

                        // cluster of cells to coarsen
                        Cell[] CellS = jCellS.Select(j => this.Cells.GetCell(j)).ToArray();

                        int CoarseningClusterID = CellS[0].CoarseningClusterID;
                        int iKref = this.Cells.GetRefElementIndex(jCellS[0]);
                        int RefinementLevel = CellS[0].RefinementLevel;
                        
                        if(jCellS.Length != KrefS_SubdivLeaves[iKref].Length)
                            throw new ArgumentException("Number of elements in coarsening cluster does not match refinement template for respective element type.");
                        if(RefinementLevel <= 0 || CoarseningClusterID <= 0)
                            throw new ArgumentException("Coarsening not available for respective cell.");

                        if(CellS.Where(cl => cl.ParentCell != null).Count() != 1)
                            throw new ArgumentException("Coarsening cluster seems wrong, or internal data may be corrupted.");

                        for(int z = 0; z < CellS.Length; z++) {
                            int j = jCellS[z];

                            if(CellsToRefineBitmask[j] == true)
                                throw new ArgumentException("Cannot refine and coarsen the same cell.");
                            if(CellsToCoarseBitmask[j] == true)
                                throw new ArgumentException("Double entry.", "CellsToCoarsen");
                            CellsToCoarseBitmask[j] = true;

                            Cell Cj = this.Cells.GetCell(j);
                            //if(Cj.CoarseningPeers == null)
                            //    throw new ArgumentException("Coarsening not available for respective cell.");
                            //if(Cj.CoarseningPeers.Length != jCellS.Length - 1)
                            //    throw new ArgumentException("Coarsening cluster seems incomplete.");
                            //if(Cj.CoarseningPeers.Length != jCellS.Length - 1)
                            //    throw new ArgumentException("Coarsening cluster seems incomplete.");
                            //foreach(long gid in jCellS) {
                            //    if(CellS.Where(cl => cl.GlobalID == gid).Count() != 1)
                            //        throw new ArgumentException("Coarsening cluster seems incomplete.");
                            //}
                            if(CoarseningClusterID != Cj.CoarseningClusterID)
                                throw new ArgumentException("Mismatch of 'CoarseningClusterID' within cluster.");

                            int[] Neighs, dummy;
                            this.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out Neighs, out dummy);

                            foreach(int jNeigh in Neighs) {
                                if(Array.IndexOf(jCellS, jNeigh) < 0) {
                                    AdaptNeighborsBitmask[jNeigh] = true;
                                }
                            }
                        }
                    }
                }


                int InsertCounter = J;
                int test = MpiRank;
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);


                // create new cells
                // ================

                //Debug.Assert(this.MpiSize == 1, "still need to adjust the following lines.");

                long GlobalIdCounter = oldGrid.NumberOfCells_l;
                Convert.ToDouble(GlobalIdCounter).MPISum();
                Convert.ToInt64(GlobalIdCounter);
                //int PtrNewCells = oldGrid.NoOfUpdateCells;
                //newGrid.Cells = new Cell[NewNoOfCells];
                List<Cell> newCells = new List<Cell>();
                int newVertexCounter = oldGrid.Cells.Max(cl => cl.NodeIndices.Max()) + 1;
                Cell[][] adaptedCells = new Cell[J][];

                Old2New.OldGlobalId = this.CurrentGlobalIdPermutation.Values.CloneAs();
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
                        //Debug.Assert(CellS.Where(cl => cl.GlobalID == Mother.GlobalID).Count() == 1);
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

                Debugger.Launch();
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
                Console.WriteLine("jgadjflg " + ReceiveOffset[0] + "ishgosdfg " + ReceiveOffset[1]);
                int Offset = 0;
                for (int m = 0; m < MpiRank; m++)
                {
                    Offset += Convert.ToInt32(ReceiveOffset[m]) * 3;
                }

                if (CellsToRefine != null) {

                    int NewCoarseningClusterId;
                    {
                        int[] locData = new int[]{
                           (this.m_Grid.Cells.Max(cl => cl.CoarseningClusterID)),
                           (CellsToRefine.Count() + 1)
                        };
                        int[] glbData = locData.MPIMax();

                        NewCoarseningClusterId = glbData[0] + 1 + glbData[1] * this.MpiRank;
                    }

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

                        Debug.Assert(Old2New.MappingIndex[j] == null);
                        Old2New.MappingIndex[j] = new int[Leaves.Length];
                        for(int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 1: create new cells
                            // create new cell
                            Cell newCell = new Cell();
                            newCell.Type = oldCell.Type;
                            if(iSubDiv == 0) {
                                newCell.GlobalID = oldCell.GlobalID;
                                newCell.ParentCell = oldCell.CloneAs();
                            }
                            else
                            {
                                long currentID = GlobalIdCounter + Offset;
                                GlobalIdCounter++;
                                newCell.GlobalID = currentID;
                                Debug.Assert(currentID != 0);
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
                                newCell.NodeIndices[i] = newVertexCounter;
                                newVertexCounter++;
                            }

                            // correlation
                            Old2New.MappingIndex[j][iSubDiv] = iSubDiv;
                        }
                        NewCoarseningClusterId++;

                        for(int iSubDiv = 0; iSubDiv < Leaves.Length; iSubDiv++) { // pass 2: do other things
                            //// record information for (later) coarsening
                            //refinedCells[iSubDiv].CoarseningPeers = refinedCells
                            //    .Where(cell => cell.GlobalID != refinedCells[iSubDiv].GlobalID)
                            //    .Select(cell => cell.GlobalID)
                            //    .ToArray();

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
                    if (jCell2 < 0) {
                        Debug.Assert((CellsToRefineBitmask[jCell1] && CellsToCoarseBitmask[jCell1]) == false);
                        if((CellsToRefineBitmask[jCell1] || CellsToCoarseBitmask[jCell1]) == false)
                            continue;

                        Cell[] adaptedBCells1 = adaptedCells[jCell1];
                        Debug.Assert(adaptedBCells1 != null);

                        int iBFace = Edge2Face[iEdge, 0];

                        foreach(Cell cl in adaptedBCells1) {
                            if(cl.CellFaceTags.Where(cft => cft.FaceIndex == iBFace).Count() == 0 && this.Edges.EdgeTags[iEdge] > 0) {
                                ArrayTools.AddToArray(new CellFaceTag() {
                                    EdgeTag = this.Edges.EdgeTags[iEdge],
                                    ConformalNeighborship = false,
                                    NeighCell_GlobalID = long.MinValue,
                                    FaceIndex = iBFace
                                }, ref cl.CellFaceTags);
                            }
                        }

                        continue;
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

                    Debug.Assert((CellsToRefineBitmaskExchange[jCell1] && CellsToCoarseBitmaskExchange[jCell1]) == false);
                    Debug.Assert((CellsToRefineBitmaskExchange[jCell2] && CellsToCoarseBitmaskExchange[jCell2]) == false);

                    bool C1changed = CellsToRefineBitmaskExchange[jCell1] || CellsToCoarseBitmaskExchange[jCell1];
                    bool C2changed = CellsToRefineBitmaskExchange[jCell2] || CellsToCoarseBitmaskExchange[jCell2];
                   
                    if((C1changed || C2changed) == false)
                        // edge between two un-changed cells -- this neighborship remains the same.
                        continue;

                    Cell[] adaptedCells1 = adaptedCells[jCell1];
                    Cell[] adaptedCells2 = adaptedCells[jCell2];

                    if(CellsToCoarseBitmaskExchange[jCell1] && CellsToCoarseBitmaskExchange[jCell2]) {
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
                    if(CellsToRefineBitmaskExchange[jCell1]) {
                        idx1 = KrefS_Faces2Subdiv[iKref1][iFace1];
                    } else {
                        Debug.Assert(adaptedCells1.Length == 1);
                        idx1 = ONE_NULL;
                    }

                    if(CellsToRefineBitmaskExchange[jCell2]) {
                        idx2 = KrefS_Faces2Subdiv[iKref2][iFace2];
                    } else {
                        Debug.Assert(adaptedCells2.Length == 1);
                        idx2 = ONE_NULL;
                    }

                    foreach(int i1 in idx1) {

                        MultidimensionalArray VtxFace1;
                        if(CellsToRefineBitmaskExchange[jCell1]) {
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

                // add EdgeTagNames and periodic Transformations
                for (int etCnt = 1; etCnt < this.EdgeTagNames.Count; etCnt++) {
                    var etPair = this.EdgeTagNames.ElementAt(etCnt);
                    newGrid.EdgeTagNames.Add(etPair);
                }

                foreach (AffineTrafo trafo in this.Grid.PeriodicTrafo) {
                    newGrid.PeriodicTrafo.Add(trafo);
                }
                foreach (AffineTrafo itrafo in this.Grid.InversePeriodicTrafo) {
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

      
    }
}

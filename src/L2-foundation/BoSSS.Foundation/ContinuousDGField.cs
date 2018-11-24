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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using ilPSP.LinSolvers;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;


namespace BoSSS.Foundation {

    /// <summary>
    /// continuous DG field via L2-projection with continuity constraints
    /// </summary>
    public class ContinuousDGField {


        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="b"></param>
        public ContinuousDGField(Basis b) {
            m_Basis = b;
            m_grd = (GridData)b.GridDat;
            m_Mapping = new UnsetteledCoordinateMapping(b);
            m_Coordinates = MultidimensionalArray.Create(m_Mapping.LocalLength);
        }

        Basis m_Basis;

        public Basis Basis {
            get {
                return m_Basis;
            }
        }

        GridData m_grd;

        UnsetteledCoordinateMapping m_Mapping;

        MultidimensionalArray m_Coordinates;

        public MultidimensionalArray Coordinates {
            get {
                return m_Coordinates;
            }
        }

        ISparseSolver m_OpSolver;

        //int[] CellMask2Coord;

        int[] GlobalVert2Local;


        /// <summary>
        /// linear solver for the quadratic optimization problem, matrix A has to be defined! 
        /// </summary>
        public ilPSP.LinSolvers.ISparseSolver OpSolver {
            get {
                if (m_OpSolver == null) {
                    //var solver = new ilPSP.LinSolvers.monkey.CG();
                    //solver.MatrixType = ilPSP.LinSolvers.monkey.MatrixType.Auto;
                    //solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.CPU;
                    //solver.Tolerance = 1.0e-12;
                    //var solver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
                    var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
                    m_OpSolver = solver;
                }

                return m_OpSolver;
            }
        }



        public void ProjectDGField(double alpha, ConventionalDGField DGField, CellMask mask = null) {
            if (DGField.Basis.Degree > this.m_Basis.Degree)
                throw new ArgumentException("continuous projection on a lower degree basis is not recommended");

            if (mask == null) {
                mask = CellMask.GetFullMask(m_grd);
            }

            // hack
            //CellMask2Coord = new int[m_grd.Cells.NoOfLocalUpdatedCells];
            List<int> maskedVert = new List<int>();
            //int numCoord = 0;
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    //CellMask2Coord[j] = numCoord;
                    //numCoord += DGField.Basis.GetLength(j);
                    int[] vertAtCell = m_grd.Cells.CellVertices[j];
                    foreach (int vert in vertAtCell) {
                        if (!maskedVert.Contains(vert)) {
                            maskedVert.Add(vert);
                        }
                    }
                }
            }
            //m_Coordinates = MultidimensionalArray.Create(numCoord);
            GlobalVert2Local = new int[m_grd.Vertices.NoOfNodes4LocallyUpdatedCells];
            int localVertInd = 0;
            foreach (int vert in maskedVert) {
                GlobalVert2Local[vert] = localVertInd;
                localVertInd++;
            }


            int degree = m_Basis.Degree;

            // get DG-coordinates (change of basis for projection on a higher polynomial degree)
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int N = DGField.Basis.GetLength(j);
                    for (int n = 0; n < N; n++) {
                        m_Coordinates[m_Mapping.LocalUniqueCoordinateIndex(0, j, n)] = DGField.Coordinates[j, n];
                        //m_Coordinates[CellMask2Coord[j] + n] = DGField.Coordinates[j, n];
                    }
                }
            }

            // construction of constraints matrix A
            // ====================================

            SubGrid maskSG = new SubGrid(mask);
            EdgeMask innerEM = maskSG.InnerEdgesMask;


            // get interpolation points for the continuity constraints
            NodeSet qNodes = getEdgeInterpolationNodes(0, 0);


            //MultidimensionalArray B = MultidimensionalArray.Create(innerEM.NoOfItemsLocally * qNodes.NoOfNodes, (int)m_Mapping.GlobalCount);
            //MultidimensionalArray B = MultidimensionalArray.Create(innerEM.NoOfItemsLocally * qNodes.NoOfNodes, m_Coordinates.Length);
            List<int> AcceptedEdges = new List<int>();
            List<NodeSet> AcceptedNodes = new List<NodeSet>();



            int nodeCount = 0;
            //int NumVert = m_grd.Vertices.NoOfNodes4LocallyUpdatedCells;
            //MultidimensionalArray CondAtVert = MultidimensionalArray.Create(NumVert, 3);    // 1. index: local vertice index / 2. index: 0 = actual cond, 1 = potential own cond, 2 = potential cond
            //MultidimensionalArray CondIncidenceMatrix = MultidimensionalArray.Create(NumVert, NumVert, 3);  // upper triangular matrix in the first two indices
            int NumVert = localVertInd;
            MultidimensionalArray CondAtVert = MultidimensionalArray.Create(NumVert, 3);
            MultidimensionalArray CondIncidenceMatrix = MultidimensionalArray.Create(NumVert, NumVert, 3);

            int[][] EdgeSend = m_grd.Edges.EdgeSendLists;
            int[][] EdgeInsert = m_grd.Edges.EdgeInsertLists;
            List<int> InterprocEdges = new List<int>();
            List<List<int>> VertAtInterprocEdges = new List<List<int>>();
            bool isInterprocEdge;
            int neighbourProc;
            int[] local2Interproc = new int[NumVert];
            local2Interproc.SetAll(-1);
            List<List<int>> ProcsAtInterprocVert = new List<List<int>>();
            //bool nonConformEdge;
            //BitArray CellsAtNonConform = new BitArray(m_grd.Cells.NoOfLocalUpdatedCells);
            //List<int> nonConformEdges = new List<int>();
            //List<List<int>> VertAtNonConformEdges = new List<List<int>>();

            // local edges per process
            foreach (var chunk in innerEM) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {

                    var edgeInfo = m_grd.Edges.Info[j];

                    // check for interprocess edges to be processed in a second step
                    isInterprocEdge = false;
                    neighbourProc = -1;
                    if (edgeInfo.HasFlag(EdgeInfo.Interprocess)) {
                        isInterprocEdge = true;
                        InterprocEdges.Add(j);
                        //foreach (int[] list in externalEdges) {
                        //    if (list != null && list.Contains(j)) {
                        //        InterprocEdges.Remove(j);
                        //    }
                        //}
                        for (int iL = 0; iL < EdgeSend.Count(); iL++) {
                            int[] sendList = EdgeSend[iL];
                            int[] insertList = EdgeInsert[iL];
                            if ((sendList != null && sendList.Contains(j)) || (insertList != null && insertList.Contains(j))) {     
                                if (sendList != null && sendList.Contains(j)) {     // find owned edge
                                    InterprocEdges.Remove(j);
                                }
                                neighbourProc = iL;
                                break;
                            } 
                        }
                    }


                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];


                    // check for cells which are located at a non-conformal cell (hanging nodes)
                    //nonConformEdge = false;
                    //if (edgeInfo.HasFlag(EdgeInfo.Cell1_Nonconformal)) {
                    //    CellsAtNonConform[cell2] = true;
                    //    nonConformEdges.Add(j);
                    //    nonConformEdge = true;
                    //} else if (edgeInfo.HasFlag(EdgeInfo.Cell2_Nonconformal)) {
                    //    CellsAtNonConform[cell1] = true;
                    //    nonConformEdges.Add(j);
                    //    nonConformEdge = true;
                    //};


                    // check the conditions at vertices (actual, own, potential)
                    int[] vertAtCell1 = m_grd.Cells.CellVertices[cell1];
                    int[] vertAtCell2 = m_grd.Cells.CellVertices[cell2];
                    int numVCond = 0;
                    List<int> VertAtEdge = new List<int>();
                    for (int i = 0; i < vertAtCell1.Length; i++) {
                        int vert = vertAtCell1[i];
                        if (vertAtCell2.Contains(vert)) {
                            VertAtEdge.Add(vert);
                            CondAtVert[GlobalVert2Local[vert], 2] += 1;   // potential condition at vert
                            if (!isInterprocEdge) { // && !nonConformEdge) {
                                CondAtVert[GlobalVert2Local[vert], 0] += 1;   // actual condition at vert 
                                CondAtVert[GlobalVert2Local[vert], 1] += 1;   // potential own cond at vert (local edge)
                            } else {
                                if (InterprocEdges.Contains(j)) {
                                    CondAtVert[GlobalVert2Local[vert], 1] += 1;   // potential own cond at vert (interproc edge)
                                }
                                if (local2Interproc[GlobalVert2Local[vert]] == -1) {
                                    local2Interproc[GlobalVert2Local[vert]] = ProcsAtInterprocVert.Count();
                                    List<int> procsAtVert = new List<int>();
                                    procsAtVert.Add(neighbourProc);
                                    ProcsAtInterprocVert.Add(procsAtVert);
                                } else {
                                    List<int> procsAtVert = ProcsAtInterprocVert.ElementAt(local2Interproc[GlobalVert2Local[vert]]);
                                    if (!procsAtVert.Contains(neighbourProc)) {
                                        procsAtVert.Add(neighbourProc);
                                    }
                                }
                            }
                            if (CondAtVert[GlobalVert2Local[vert], 0] == 4) {
                                numVCond++;
                            }
                        }
                    }
                    //if (nonConformEdge) {
                    //    VertAtNonConformEdges.Add(VertAtEdge);
                    //}
                    if (isInterprocEdge && InterprocEdges.Contains(j)) {
                        VertAtInterprocEdges.Add(VertAtEdge);
                    }


                    // check for overdetermined edges (additional for 3D)
                    int numECond = 0;
                    int[] edgeOrientation = new int[2];
                    if (m_grd.SpatialDimension == 3) {
                        int i0 = 0;
                        for (int indV1 = 0; indV1 < 4; indV1++) {
                            i0++;
                            for (int indV2 = i0; indV2 < 4; indV2++) {
                                int m, n;
                                if (VertAtEdge.ElementAt(indV1) <= VertAtEdge.ElementAt(indV2)) {
                                    m = VertAtEdge.ElementAt(indV1);
                                    n = VertAtEdge.ElementAt(indV2);
                                } else {
                                    m = VertAtEdge.ElementAt(indV2);
                                    n = VertAtEdge.ElementAt(indV1);
                                }
                                int m_loc = GlobalVert2Local[m];
                                int n_loc = GlobalVert2Local[n];
                                CondIncidenceMatrix[m_loc, n_loc, 2] += 1;  // potential condition at edge
                                if (!isInterprocEdge) {
                                    CondIncidenceMatrix[m_loc, n_loc, 0] += 1;   // actual condition at edge 
                                    CondIncidenceMatrix[m_loc, n_loc, 1] += 1;   // potential own cond at edge (local edge/face)
                                } else {
                                    if (InterprocEdges.Contains(j)) {
                                        CondIncidenceMatrix[m_loc, n_loc, 1] += 1;   // potential own cond at vert (interproc edge)
                                    }
                                }
                                if (CondIncidenceMatrix[m_loc, n_loc, 0] == 4) {
                                    numECond++;
                                    if (m_grd.Vertices.Coordinates[m, 0] == m_grd.Vertices.Coordinates[n, 0]) {
                                        edgeOrientation[0] = 1;
                                    } else if (m_grd.Vertices.Coordinates[m, 1] == m_grd.Vertices.Coordinates[n, 1]) {
                                        edgeOrientation[1] = 1;
                                    } else {
                                        throw new ApplicationException("should not occur");
                                    }

                                }
                            }
                        }
                    }

                   
                    if (!isInterprocEdge) { // && !nonConformEdge) {

                        AcceptedEdges.Add(j);

                        qNodes = getEdgeInterpolationNodes(numVCond, numECond, edgeOrientation);
                        AcceptedNodes.Add(qNodes);

                        //if (qNodes != null) {
                        //    Console.WriteLine("proc {0}: numECond = {1}, No of qNodes = {2}", m_grd.MpiRank, numECond, qNodes.NoOfNodes);
                        //} else {
                        //    Console.WriteLine("proc {0}: numECond = {1}, No qNodes", m_grd.MpiRank, numECond);
                        //}

                        if (qNodes != null) {

                            //// set continuity constraints
                            //var results = m_Basis.EdgeEval(qNodes, j, 1);

                            //for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                            //    // Cell1
                            //    for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                            //        //B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                            //        B[nodeCount + qN, CellMask2Coord[cell1]] = results.Item1[0, qN, p];
                            //    }
                            //    // Cell2
                            //    for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                            //        //B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
                            //        B[nodeCount + qN, CellMask2Coord[cell2]] = results.Item1[0, qN, p];
                            //    }
                            //}
                            nodeCount += qNodes.NoOfNodes;
                        }
                    }

                }
            }

            #region hanging nodes

            //BitArray ishangingNode = new BitArray(m_grd.Vertices.NoOfNodes4LocallyUpdatedCells);

            //// get edges at hanging nodes
            //CellMask CaNCmsk = new CellMask(m_grd, CellsAtNonConform);
            //SubGrid CaNCsgrd = new SubGrid(CaNCmsk);
            //EdgeMask hangingEdges = CaNCsgrd.InnerEdgesMask;

            //// find corresponding hanging nodes
            //int[] hangVert2NonConformCell = new int[m_grd.Vertices.NoOfNodes4LocallyUpdatedCells];
            //hangVert2NonConformCell.SetAll(-1);
            //foreach (var chunk in hangingEdges) {
            //    int j0 = chunk.i0;
            //    int jE = chunk.JE;
            //    for (int j = j0; j < jE; j++) {

            //        int cell1 = m_grd.Edges.CellIndices[j, 0];
            //        int cell2 = m_grd.Edges.CellIndices[j, 1];
            //        int[] cell1_neighbours = m_grd.Cells.CellNeighbours[cell1];
            //        int[] cell2_neighbours = m_grd.Cells.CellNeighbours[cell2];
            //        int NonConformCell = -1;
            //        for (int i = 0; i < cell1_neighbours.Length; i++) {
            //            int cell = cell1_neighbours[i];
            //            if (cell2_neighbours.Contains(cell)) {
            //                NonConformCell = cell;
            //                break;
            //            }
            //        }
            //        if (NonConformCell > -1) {

            //            int[] vertAtCell1 = m_grd.Cells.CellVertices[cell1];
            //            int[] vertAtCell2 = m_grd.Cells.CellVertices[cell2];
            //            for (int i = 0; i < vertAtCell1.Length; i++) {
            //                int vert = vertAtCell1[i];
            //                if (vertAtCell2.Contains(vert) && m_grd.Vertices.VerticeToCell[vert].Length == 2) {
            //                    ishangingNode[vert] = true;
            //                    hangVert2NonConformCell[vert] = NonConformCell;
            //                    break;
            //                }
            //            }

            //        }
            //    }
            //}


            //// non-confrom edges (only for single core in 2D so far!)
            //foreach (int j in nonConformEdges) {

            //    int cell1 = m_grd.Edges.CellIndices[j, 0];
            //    int cell2 = m_grd.Edges.CellIndices[j, 1];
            //    int hangCell;
            //    int NonConformCell;
            //    if (CellsAtNonConform[cell1]) {
            //        hangCell = cell1;
            //        NonConformCell = cell2;
            //    } else {
            //        hangCell = cell2;
            //        NonConformCell = cell1;
            //    }

            //    int[] vertAtHCell = m_grd.Cells.CellVertices[hangCell];
            //    int[] vertAtNcCell = m_grd.Cells.CellVertices[NonConformCell];
            //    int numVCond = 0;
            //    foreach (int vert in vertAtHCell) {
            //        if (vertAtNcCell.Contains(vert)) {
            //            CondAtVert[vert, 0] += 1;
            //            CondAtVert[vert, 1] += 1;
            //            CondAtVert[vert, 2] += 1;
            //            if (CondAtVert[vert, 0] > maxVCond) {
            //                numVCond++;
            //            }
            //        }
            //        if (ishangingNode[vert] && hangVert2NonConformCell[vert] == NonConformCell) {
            //            CondAtVert[vert, 0] += 1;
            //            CondAtVert[vert, 1] += 1;
            //            CondAtVert[vert, 2] += 1;
            //            if (CondAtVert[vert, 0] == maxVCond) {
            //                numVCond++;
            //            }
            //        }
            //    }


            //    qNodes = getEdgeInterpolationNodes(numVCond, 0);

            //    if (qNodes != null) {

            //        // set continuity constraints
            //        var results = m_Basis.EdgeEval(qNodes, j, 1);

            //        for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
            //            // Cell1
            //            for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
            //                B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
            //            }
            //            // Cell2
            //            for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
            //                B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
            //            }
            //        }
            //        nodeCount += qNodes.NoOfNodes;
            //    }


            //}



            // additional conditions for hanging nodes
            //ishangingNode = new BitArray(m_grd.Vertices.NoOfNodes4LocallyUpdatedCells);

            // get edges at hanging nodes
            //CellMask CaNCmsk = new CellMask(m_grd, CellsAtNonConform);
            //SubGrid CaNCsgrd = new SubGrid(CaNCmsk);
            //EdgeMask hangingEdges = CaNCsgrd.InnerEdgesMask;

            //// find corresponding hanging nodes
            //hangVert2hangEdge = new int[m_grd.Vertices.NoOfNodes4LocallyUpdatedCells];
            //hangVert2hangEdge.SetAll(-1);         
            //hangVert2NonConformCell = new int[m_grd.Vertices.NoOfNodes4LocallyUpdatedCells];
            //hangVert2NonConformCell.SetAll(-1);
            //int[] hangEdge2hangVert = new int[hangingEdges.NoOfItemsLocally];
            //int Eind = 0;
            //int HNcount = 0;
            //foreach (var chunk in hangingEdges) {
            //    int j0 = chunk.i0;
            //    int jE = chunk.JE;
            //    for (int j = j0; j < jE; j++) {

            //        int cell1 = m_grd.Edges.CellIndices[j, 0];
            //        int cell2 = m_grd.Edges.CellIndices[j, 1];
            //        int[] cell1_neighbours = m_grd.Cells.CellNeighbours[cell1];
            //        int[] cell2_neighbours = m_grd.Cells.CellNeighbours[cell2];
            //        int NonConformCell = -1;
            //        for (int i = 0; i < cell1_neighbours.Length; i++) {
            //            int cell = cell1_neighbours[i];
            //            if (cell2_neighbours.Contains(cell)) {
            //                NonConformCell = cell;
            //                break;
            //            }
            //        }
            //        if (NonConformCell > -1) {

            //            int[] vertAtCell1 = m_grd.Cells.CellVertices[cell1];
            //            int[] vertAtCell2 = m_grd.Cells.CellVertices[cell2];
            //            for (int i = 0; i < vertAtCell1.Length; i++) {
            //                int vert = vertAtCell1[i];
            //                if (vertAtCell2.Contains(vert) && m_grd.Vertices.VerticeToCell[vert].Length == 2) {
            //                    ishangingNode[vert] = true;
            //                    hangVert2hangEdge[vert] = j;
            //                    hangVert2NonConformCell[vert] = NonConformCell;
            //                    hangEdge2hangVert[Eind] = vert;
            //                    HNcount++;
            //                    break;
            //                }
            //            }

            //        } else {
            //            hangEdge2hangVert[Eind] = -1;
            //        }
            //        Eind++;
            //    }
            //}

            // compute derivatives of the basis polynomials along on spatial direction
            //PolynomialList[,] edgeDeriv = ComputePartialDerivatives(m_Basis);

            // compute additional constraints at hanging edges
            //var Trafo = m_grd.ChefBasis.Scaling;
            //MultidimensionalArray B2 = MultidimensionalArray.Create(HNcount * m_Basis.Degree, (int)m_Mapping.GlobalCount);
            //Eind = 0;
            //HNcount = 0;
            //foreach (var chunk in hangingEdges) {
            //    int j0 = chunk.i0;
            //    int jE = chunk.JE;
            //    for (int j = j0; j < jE; j++) {

            //        if (hangEdge2hangVert[Eind] != -1) {

            //            int cell1 = m_grd.Edges.CellIndices[j, 0];
            //            int cell2 = m_grd.Edges.CellIndices[j, 1];

            //            int trf1 = m_grd.Edges.Edge2CellTrafoIndex[j, 0];
            //            int trf2 = m_grd.Edges.Edge2CellTrafoIndex[j, 1];

            //            MultidimensionalArray hNd_global = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(new int[] { hangEdge2hangVert[Eind], 0 },
            //                                                                        new int[] { hangEdge2hangVert[Eind], m_grd.SpatialDimension - 1 });

            //            MultidimensionalArray hND_local1 = MultidimensionalArray.Create(1, 1, m_grd.SpatialDimension);
            //            MultidimensionalArray hND_local2 = MultidimensionalArray.Create(1, 1, m_grd.SpatialDimension);

            //            m_grd.TransformGlobal2Local(hNd_global, hND_local1, cell1, 1, 0);
            //            m_grd.TransformGlobal2Local(hNd_global, hND_local2, cell2, 1, 0);

            //            NodeSet hangNode1 = new NodeSet(m_grd.Grid.GetRefElement(0), hND_local1.ExtractSubArrayShallow(0, -1, -1));
            //            NodeSet hangNode2 = new NodeSet(m_grd.Grid.GetRefElement(0), hND_local2.ExtractSubArrayShallow(0, -1, -1));

            //            int derivInd = -1;
            //            for (int d = 0; d < m_grd.SpatialDimension; d++) {
            //                if (m_grd.Edges.NormalsForAffine[j, d] != 0)
            //                    derivInd = d;
            //            }

            //            for (int drv = 0; drv < m_Basis.Degree; drv++) {
            //                MultidimensionalArray Res = MultidimensionalArray.Create(1, m_Basis.Polynomials[0].Count);
            //                edgeDeriv[drv, derivInd].Evaluate(hangNode1, Res);
            //                for (int p = 0; p < m_Basis.Polynomials[0].Count; p++) {
            //                    B2[HNcount * m_Basis.Degree + drv, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = Res[0, p] * Trafo[cell1];
            //                }
            //                Res.Clear();
            //                edgeDeriv[drv, derivInd].Evaluate(hangNode2, Res);
            //                for (int p = 0; p < m_Basis.Polynomials[0].Count; p++) {
            //                    B2[HNcount * m_Basis.Degree + drv, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -Res[0, p] * Trafo[cell2];
            //                }
            //            }
            //            HNcount++;
            //        }
            //        Eind++;
            //    }
            //}


            //// non-confrom edges (only for single core in 2D so far!)
            //int edge = 0;
            //foreach (int j in nonConformEdges) {

            //    List<int> vertAtEdge = VertAtNonConformEdges.ElementAt(edge);
            //    bool internalNonConformEdge = true;
            //    if (vertAtEdge.Count > 0) {
            //        internalNonConformEdge = false;
            //    }

            //    if (!internalNonConformEdge) {

            //        int cell1 = m_grd.Edges.CellIndices[j, 0];
            //        int cell2 = m_grd.Edges.CellIndices[j, 1];
            //        int hangCell;
            //        int NonConformCell;
            //        if (CellsAtNonConform[cell1]) {
            //            hangCell = cell1;
            //            NonConformCell = cell2;
            //        } else {
            //            hangCell = cell2;
            //            NonConformCell = cell1;
            //        }

            //        int hangNode = -1;
            //        foreach (var vert in m_grd.Cells.CellVertices[hangCell]) {
            //            if (ishangingNode[vert] && hangVert2NonConformCell[vert] == NonConformCell) {
            //                hangNode = vert;
            //            }
            //        }

            //        if (hangNode == -1) {
            //            throw new ArgumentOutOfRangeException("should not happen");
            //        }

            //        bool NonConfromEdgeIsProcessed = false;
            //        if (CondAtVert[hangNode, 0] == maxVCond) {
            //            NonConfromEdgeIsProcessed = true;
            //        }


            //        if (!NonConfromEdgeIsProcessed) {

            //            int numVCond = 0;
            //            foreach (int vert in vertAtEdge) {
            //                CondAtVert[vert, 0] += 1;
            //                CondAtVert[vert, 1] += 1;
            //                CondAtVert[vert, 2] += 1;
            //                if (CondAtVert[vert, 0] > maxVCond) {
            //                    numVCond++;
            //                }
            //            }

            //            CondAtVert[hangNode, 0] += 1;
            //            CondAtVert[hangNode, 1] += 1;
            //            CondAtVert[hangNode, 2] += 1;

            //            ProcessNonConformEdge(hangNode, hangCell, NonConformCell,  ref numVCond);


            //            qNodes = getEdgeInterpolationNodes(numVCond, 0);

            //            if (qNodes != null) {

            //                // set continuity constraints
            //                var results = m_Basis.EdgeEval(qNodes, j, 1);

            //                for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
            //                    // Cell1
            //                    for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
            //                        B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
            //                    }
            //                    // Cell2
            //                    for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
            //                        B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
            //                    }
            //                }
            //                nodeCount += qNodes.NoOfNodes;
            //            }

            //        }
            //    }
            //    edge++;
            //}
#endregion

            //int nodeCount_OnProc = nodeCount;


            // interprocess edges 
            int edge = 0;
            foreach (int j in InterprocEdges) {


                int cell1 = m_grd.Edges.CellIndices[j, 0];
                int cell2 = m_grd.Edges.CellIndices[j, 1];

                int numVCond = 0;
                int numECond = 0;
                List<int> VertAtEdge = VertAtInterprocEdges.ElementAt(edge);

                switch (m_grd.SpatialDimension) {
                    case 1: {
                            throw new NotImplementedException("TODO");
                        }
                    case 2: {
                            // condtions at vertices
                            foreach (int vert in VertAtEdge) {
                                numVCond += NegotiateNumVCond(vert, CondAtVert, ProcsAtInterprocVert[local2Interproc[vert]]);
                                CondAtVert[vert, 0] += 1;
                            }
                            break;
                        }
                    case 3: {
                            // conditions at edges 
                            //Console.WriteLine("proc {0}: edge {1}", m_grd.MpiRank, j);
                            int i0 = 0;
                            for (int indV1 = 0; indV1 < VertAtEdge.Count; indV1++) {
                                i0++;
                                for (int indV2 = i0; indV2 < VertAtEdge.Count; indV2++) {
                                    int m, n;
                                    if (VertAtEdge.ElementAt(indV1) <= VertAtEdge.ElementAt(indV2)) {
                                        m = VertAtEdge.ElementAt(indV1);
                                        n = VertAtEdge.ElementAt(indV2);
                                    } else {
                                        m = VertAtEdge.ElementAt(indV2);
                                        n = VertAtEdge.ElementAt(indV1);
                                    }

                                    numECond += NegotiateNumECond(m, n, CondIncidenceMatrix, ProcsAtInterprocVert[local2Interproc[m]]);
                                    CondIncidenceMatrix[m, n, 0] += 1;
                                }
                            }
                            break;
                        }
                }


                //if (VertAtEdge.Count == 4) {
                //    int i0 = 0;
                //    for (int indV1 = 0; indV1 < VertAtEdge.Count; indV1++) {
                //        i0++;
                //        for (int indV2 = i0; indV2 < VertAtEdge.Count; indV2++) {
                //            int m, n;
                //            if (VertAtEdge.ElementAt(indV1) <= VertAtEdge.ElementAt(indV2)) {
                //                m = VertAtEdge.ElementAt(indV1);
                //                n = VertAtEdge.ElementAt(indV2);
                //            } else {
                //                m = VertAtEdge.ElementAt(indV2);
                //                n = VertAtEdge.ElementAt(indV1);
                //            }
                //            CondIncidenceMatrix[m, n, 0] += 1;
                //        }
                //    }
                //}

                AcceptedEdges.Add(j);

                qNodes = getEdgeInterpolationNodes(numVCond, numECond);
                AcceptedNodes.Add(qNodes);

                //if (qNodes != null) {
                //    Console.WriteLine("proc {0}: numECond = {1}, No of qNodes = {2} - interproc Edge", m_grd.MpiRank, numECond, qNodes.NoOfNodes);
                //} else {
                //    Console.WriteLine("proc {0}: numECond = {1}, No qNodes - interproc Edge", m_grd.MpiRank, numECond);
                //}

                if (qNodes != null) {

                    // set continuity constraints
                    //var results = m_Basis.EdgeEval(qNodes, j, 1);

                    //for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                    //    // Cell1
                    //    for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                    //        B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                    //    }
                    //    // Cell2
                    //    for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                    //        B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
                    //    }
                    //}
                    nodeCount += qNodes.NoOfNodes;
                }
                edge++;

            }


            Partitioning rowPart = new Partitioning(nodeCount);
            MsrMatrix A = new MsrMatrix(rowPart, m_Mapping);

            int count = 0;
            nodeCount = 0;
            foreach (int j in AcceptedEdges) {

                int cell1 = m_grd.Edges.CellIndices[j, 0];
                int cell2 = m_grd.Edges.CellIndices[j, 1];


                // set continuity constraints
                qNodes = AcceptedNodes.ElementAt(count);

                var results = m_Basis.EdgeEval(qNodes, j, 1);

                for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                    // Cell1
                    for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                        A[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                    }
                    // Cell2
                    for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                        A[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
                    }
                }
                count++;
                nodeCount += qNodes.NoOfNodes;
                    

         
            }


            //Partitioning rowPart = new Partitioning(nodeCount);
            //MsrMatrix A = new MsrMatrix(rowPart, m_Mapping);
            //A.AccBlock(rowPart.i0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, (int)m_Mapping.GlobalCount - 1 }));
            //MsrMatrix A = new MsrMatrix(nodeCount, m_Coordinates.Length, 1, 1);
            //A.AccBlock(0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, m_Coordinates.Length - 1 }));
            //A.AccBlock(rowPart.i0 + nodeCount, 0, 1.0, B2);

            //A.SaveToTextFile("C:\\tmp\\AMatrix.txt");

            // test with matlab
            //MultidimensionalArray output = MultidimensionalArray.Create(1, 2);
            //Console.WriteLine("Calling MATLAB/Octave...");
            //using (BatchmodeConnector bmc = new BatchmodeConnector()) {
            //    bmc.PutSparseMatrix(A, "A");
            //    bmc.Cmd("rank_A = rank(full(A))");
            //    bmc.Cmd("rank_AT = rank(full(A'))");
            //    bmc.GetMatrix(output, "[rank_A, rank_AT]");

            //    bmc.Execute(false);
            //}

            //Console.WriteLine("A: No of Rows = {0}; rank = {1}", A.NoOfRows, output[0, 0]);
            //Console.WriteLine("AT: No of Rows = {0}; rank = {1}", A.Transpose().NoOfRows, output[0, 1]);


            // solve
            MsrMatrix AAT = A * A.Transpose();

            double[] RHS = new double[rowPart.LocalLength];
            A.SpMVpara(1.0, m_Coordinates.To1DArray(), 0.0, RHS);

            double[] v = new double[rowPart.LocalLength];
            double[] x = new double[m_Coordinates.Length];

            //OpSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            OpSolver.DefineMatrix(AAT);
            OpSolver.Solve(v, RHS);

            A.Transpose().SpMVpara(-1.0, v, 0.0, x);

            m_Coordinates.AccVector(1.0, x);

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="numVCond"></param>
        /// <param name="numECond"></param>
        /// <returns></returns>
        private NodeSet getEdgeInterpolationNodes(int numVCond, int numECond, int[] edgeOrientation = null) {

            int degree = m_Basis.Degree;

            switch (m_grd.SpatialDimension) {
                case 1: {
                        throw new NotImplementedException("TODO");
                    }
                case 2: {
                        if (numVCond > degree) {
                            return null;
                        } else {
                            QuadRule quad = m_grd.Edges.EdgeRefElements[0].GetQuadratureRule((degree - numVCond) * 2);
                            return quad.Nodes;
                        }
                    }
                case 3: {
                        if (numECond > degree) {
                            return null;
                        } else {
                            QuadRule quad = m_grd.Edges.EdgeRefElements[0].FaceRefElement.GetQuadratureRule((degree - numECond) * 2);
                            NodeSet qNodes = quad.Nodes;
                            int Nnds = (((degree - numECond) + 1) * ((degree - numECond) + 2) / 2);
                            MultidimensionalArray nds = MultidimensionalArray.Create(Nnds, 2);
                            int node = 0;
                            for (int n1 = 0; n1 < quad.NoOfNodes; n1++) {
                                for (int n2 = 0; n2 < quad.NoOfNodes - n1; n2++) {
                                    nds[node, 0] = qNodes[n1, 0];
                                    nds[node, 1] = qNodes[n2, 0];
                                    node++;
                                }
                            }
                            return new NodeSet(m_grd.Edges.EdgeRefElements[0], nds);
                            
                        }
                    }
                default:
                    throw new NotSupportedException("spatial dimension not supported");
            }

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="basis"></param>
        /// <returns></returns>
        //private PolynomialList[,] ComputePartialDerivatives(Basis basis) {

        //    int deg = basis.Degree;
        //    int D = basis.GridDat.SpatialDimension;

        //    PolynomialList[,] edgeDeriv = new PolynomialList[deg, D];

        //    for (int drv = 0; drv < deg; drv++) {
        //        for (int d = 0; d < D; d++) {
        //            List<Polynomial> polyList = new List<Polynomial>();
        //            int[] deriv = new int[D];
        //            deriv[d] = drv + 1;
        //            foreach (Polynomial poly in basis.Polynomials[0]) {
        //                Polynomial polyDeriv = poly.Derive(deriv);
        //                polyList.Add(polyDeriv);
        //            }
        //            edgeDeriv[drv, d] = new PolynomialList(polyList);
        //        }
        //    }

        //    return edgeDeriv;

        //}


        //BitArray ishangingNode;

        //int[] hangVert2hangEdge;

        //int[] hangVert2NonConformCell;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="hangNode"></param>
        /// <param name="hangCell"></param>
        /// <param name="numVCond"></param>
        //private void ProcessNonConformEdge(int currentHNode, int currentHCell, int NonConformCell, ref int numVCond) {

        //    // set conditions at current hanging node
        //    CondAtVert[currentHNode, 0] += 1;
        //    CondAtVert[currentHNode, 1] += 1;
        //    CondAtVert[currentHNode, 2] += 1;

        //    // search for the next hanging node and cell / or terminate
        //    int currentHEdge = hangVert2hangEdge[currentHNode];

        //    int cell1 = m_grd.Edges.CellIndices[currentHEdge, 0];
        //    int cell2 = m_grd.Edges.CellIndices[currentHEdge, 1];

        //    int nextHCell;
        //    if (cell1 == currentHCell) {
        //        nextHCell = cell2;
        //    } else {
        //        nextHCell = cell1;
        //    }

        //    int[] vertAtHCell = m_grd.Cells.CellVertices[nextHCell];
        //    int[] vertAtNcCell = m_grd.Cells.CellVertices[NonConformCell];
        //    foreach (int vert in vertAtHCell) {
        //        if (ishangingNode[vert] && vert != currentHNode && hangVert2NonConformCell[vert] == NonConformCell) {   // next recursion
        //            // set conditions at next hanging node
        //            CondAtVert[vert, 0] += 1;
        //            CondAtVert[vert, 1] += 1;
        //            CondAtVert[vert, 2] += 1;
        //            ProcessNonConformEdge(vert, nextHCell, NonConformCell, ref numVCond);
        //        }
        //        if (vertAtNcCell.Contains(vert)) {      // terminate
        //            // set conditions at non conform edge
        //            CondAtVert[vert, 0] += 1;
        //            CondAtVert[vert, 1] += 1;
        //            CondAtVert[vert, 2] += 1;
        //            if (CondAtVert[vert, 0] > maxVCond) {
        //                numVCond++;
        //            }
        //        }
        //    }


        //}


        /// <summary>
        /// 
        /// </summary>
        /// <param name="vertAtEdge"></param>
        /// <param name="CondAtVert"></param>
        /// <returns></returns>
        private int NegotiateNumVCond(int vert, MultidimensionalArray CondAtVert, List<int> procsAtVert = null) {


            if (CondAtVert[vert, 2] == 2) {

                if (procsAtVert.Count == 1) {
                    if (CondAtVert[vert, 1] == 2) {
                        if (CondAtVert[vert, 0] == 1) {
                            return 1;
                        }
                    }
                } else if (procsAtVert.Count == 2) {
                    if (procsAtVert.Min() > m_grd.MpiRank) {
                        if (CondAtVert[vert, 0] == CondAtVert[vert, 1] - 1) {
                            return 1;
                        }
                    }
                }
            }

            if (CondAtVert[vert, 2] == 3) {

                if (procsAtVert.Count == 1) {
                    if (CondAtVert[vert, 1] == 3) {
                        if (CondAtVert[vert, 0] == 2) {
                            return 1;
                        }
                    } else if (CondAtVert[vert, 1] == 2) {
                        if (procsAtVert.Min() > m_grd.MpiRank) {
                            return 1;
                        }
                    }

                } else if (procsAtVert.Count == 2) {
                    if (CondAtVert[vert, 0] == CondAtVert[vert, 1] - 1) {
                        return 1;
                    }
                }

            }

            if (CondAtVert[vert, 2] == 4) {
                if (CondAtVert[vert, 1] == 4) {
                    if (CondAtVert[vert, 0] == 3) {
                        return 1;
                    }
                }
                if (CondAtVert[vert, 1] == 3) {
                    if (CondAtVert[vert, 0] == 2) {
                        return 1;
                    }
                }
            }


            return 0;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="vertAtEdge"></param>
        /// <param name="CondIncidenceMatrix"></param>
        /// <returns></returns>
        //private int NegotiateNumECond(List<int> VertAtEdge, MultidimensionalArray CondIncidenceMatrix) {
        private int NegotiateNumECond(int m, int n, MultidimensionalArray CondIncidenceMatrix, List<int> procsAtVertm = null, List<int> procsAtVertn = null) {

            //int numECond = 0;

            //int i0 = 0;
            //for (int indV1 = 0; indV1 < VertAtEdge.Count; indV1++) {
            //    i0++;
            //    for (int indV2 = i0; indV2 < VertAtEdge.Count; indV2++) {
            //        int m, n;
            //        if (VertAtEdge.ElementAt(indV1) <= VertAtEdge.ElementAt(indV2)) {
            //            m = VertAtEdge.ElementAt(indV1);
            //            n = VertAtEdge.ElementAt(indV2);
            //        } else {
            //            m = VertAtEdge.ElementAt(indV2);
            //            n = VertAtEdge.ElementAt(indV1);
            //        }

            if (CondIncidenceMatrix[m, n, 2] > 1) {
                //Console.WriteLine("proc {0}: condAtEdge ({1}/{2}/{3})", m_grd.MpiRank, CondIncidenceMatrix[m, n, 0], CondIncidenceMatrix[m, n, 1], CondIncidenceMatrix[m, n, 2]);
            }

            if (CondIncidenceMatrix[m, n, 2] == 2) {
                if (CondIncidenceMatrix[m, n, 1] == 2) {
                    if (CondIncidenceMatrix[m, n, 0] == 1) {
                        //Console.WriteLine("numVCond++");
                        return 1;
                    }
                }
            }


            if (CondIncidenceMatrix[m, n, 2] == 3) {
                if (CondIncidenceMatrix[m, n, 1] == 3) {
                    if (CondIncidenceMatrix[m, n, 0] == 2) {
                        //Console.WriteLine("numVCond++");
                        return 1;
                    }
                } else if (CondIncidenceMatrix[m, n, 1] == 2) {
                    if (procsAtVertm.Min() > m_grd.MpiRank) {
                        //Console.WriteLine("numVCond++");
                        return 1;
                    }
                    //if (m < m_grd.Vertices.NodePartitioning.LocalLength && n < m_grd.Vertices.NodePartitioning.LocalLength) {   // proc owns both vert
                    //    return 1;
                    //} else if (m < m_grd.Vertices.NodePartitioning.LocalLength && n >= m_grd.Vertices.NodePartitioning.LocalLength) {
                    //    int proc = 0;
                    //    foreach (int[] list in m_grd.Vertices.VertexInsertLists) {
                    //        if (list != null && list.Contains(n)) {
                    //            if (m_grd.MpiRank < proc) {
                    //                return 1;
                    //            }
                    //        }
                    //        proc++;
                    //    }
                    //} else if (m >= m_grd.Vertices.NodePartitioning.LocalLength && n < m_grd.Vertices.NodePartitioning.LocalLength) {
                    //    int proc = 0;
                    //    foreach (int[] list in m_grd.Vertices.VertexInsertLists) {
                    //        if (list != null && list.Contains(m)) {
                    //            if (m_grd.MpiRank < proc) {
                    //                return 1;
                    //            }
                    //        }
                    //        proc++;
                    //    }
                    //}
                }
            }


            if (CondIncidenceMatrix[m, n, 2] == 4) {
                if (CondIncidenceMatrix[m, n, 1] == 4) {
                    if (CondIncidenceMatrix[m, n, 0] == 3) {
                        //Console.WriteLine("numVCond++");
                        return 1;
                    }
                }
                if (CondIncidenceMatrix[m, n, 1] == 3) {
                    if (CondIncidenceMatrix[m, n, 0] == 2) {
                        //Console.WriteLine("numVCond++");
                        return 1;
                    }
                }
            }

            //    }
            //}

            return 0;
        }


        /// <summary>
        /// Accumulate this field to a DG Field
        /// </summary>
        /// <param name="alpha">Scaling factor</param>
        /// <param name="DGField"></param>
        /// <param name="mask"></param>
        public void AccToDGField(double alpha, ConventionalDGField DGField, CellMask mask = null) {
            if (!DGField.Basis.Equals(this.m_Basis))
                throw new ArgumentException("Basis does not match.");

            if (mask == null) {
                mask = CellMask.GetFullMask(m_grd);
            }

            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int JE = chunk.JE;
                for (int j = j0; j < JE; j++) {
                    // collect coordinates for cell j
                    int N = DGField.Basis.GetLength(j);
                    double[] CDGcoord = new double[N];
                    for (int n = 0; n < N; n++) {
                        CDGcoord[n] = m_Coordinates[m_Mapping.LocalUniqueCoordinateIndex(0, j, n)];
                        //CDGcoord[n] = m_Coordinates[CellMask2Coord[j] + n];
                    }

                    double[] DGcoord = new double[N];
                    DGField.Coordinates.GetRow(j, DGcoord);

                    DGcoord.AccV(alpha, CDGcoord);

                    DGField.Coordinates.SetRow(j, DGcoord);

                }
            }

        }
    }

}

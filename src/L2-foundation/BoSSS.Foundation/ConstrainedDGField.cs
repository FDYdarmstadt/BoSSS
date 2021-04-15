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
using System.Diagnostics;
using MPI.Wrappers;

namespace BoSSS.Foundation {

    /// <summary>
    /// continuous DG field via L2-projection with continuity constraints
    /// </summary>
    public class ConstrainedDGField {


        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="b"></param>
        public ConstrainedDGField(Basis b) {
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
                    //var solver = new ilPSP.LinSolvers.monkey.PCG();
                    //solver.MatrixType = ilPSP.LinSolvers.monkey.MatrixType.Auto;
                    //solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.CPU;
                    //solver.Tolerance = 1.0e-12;
                    //var solver = new ilPSP.LinSolvers.HYPRE.PCG();
                    //var precond = new ilPSP.LinSolvers.HYPRE.Euclid();
                    //solver.NestedPrecond = precond;
                    var solver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
                    //var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
                    m_OpSolver = solver;
                }

                return m_OpSolver;
            }
        }


        /// <summary>
        /// Projects some DG field <paramref name="DGField"/> onto the internal, continuous representation
        /// </summary>
        /// <param name="DGField">
        /// input; unchanged on exit
        /// </param>
        /// <param name="mask"></param>
        public void ProjectDGField2(ConventionalDGField DGField, CellMask mask = null) {
            if (DGField.Basis.Degree > this.m_Basis.Degree)
                throw new ArgumentException("continuous projection on a lower degree basis is not recommended");
            this.Coordinates.Clear(); // clear internal state, to get the same result for the same input every time

            if (mask == null) {
                mask = CellMask.GetFullMask(m_grd);
            }
            if(mask.NoOfItemsLocally.MPISum() <= 0) {
                throw new ArgumentOutOfRangeException("Domain mask cannot be empty.");
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
                                    if ((m_grd.Vertices.Coordinates[m, 0] - m_grd.Vertices.Coordinates[n, 0]) < 1.0e-15) {
                                        edgeOrientation[0] = 1;
                                    } else if ((m_grd.Vertices.Coordinates[m, 1] - m_grd.Vertices.Coordinates[n, 1]) < 1.0e-15) {
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

                        //qNodes = getEdgeInterpolationNodes(numVCond, numECond, edgeOrientation);
                        qNodes = getEdgeInterpolationNodes(numVCond, numECond);
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
                                int lvert = GlobalVert2Local[vert];// = vert, vorher IndexOutOfRange exeption
                                numVCond += NegotiateNumVCond(lvert, CondAtVert, ProcsAtInterprocVert[local2Interproc[lvert]]);
                                CondAtVert[lvert, 0] += 1;
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
                                    int v1 = GlobalVert2Local[VertAtEdge.ElementAt(indV1)];
                                    int v2 = GlobalVert2Local[VertAtEdge.ElementAt(indV2)];
                                    int m, n;
                                    if (local2Interproc[v1] >= 0 && local2Interproc[v2] >= 0) {
                                        if (v1 <= v2) {
                                            m = v1;
                                            n = v2;
                                        } else {
                                            m = v2;
                                            n = v1;
                                        }
                                    } else {
                                        throw new ArgumentException("should not occur");
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
            long _nodeCount = A.RowPartitioning.i0; // start at global index
            foreach (int j in AcceptedEdges) {

                int cell1 = m_grd.Edges.CellIndices[j, 0];
                int cell2 = m_grd.Edges.CellIndices[j, 1];


                // set continuity constraints
                qNodes = AcceptedNodes.ElementAt(count);

                var results = m_Basis.EdgeEval(qNodes, j, 1);

                for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                    // Cell1       
                    for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                        A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                    }
                    // Cell2
                    for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                        A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
                    }
                }
                count++;
                _nodeCount += qNodes.NoOfNodes;
                    

         
            }


            //Partitioning rowPart = new Partitioning(nodeCount);
            //MsrMatrix A = new MsrMatrix(rowPart, m_Mapping);
            //A.AccBlock(rowPart.i0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, (int)m_Mapping.GlobalCount - 1 }));
            //MsrMatrix A = new MsrMatrix(nodeCount, m_Coordinates.Length, 1, 1);
            //A.AccBlock(0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, m_Coordinates.Length - 1 }));
            //A.AccBlock(rowPart.i0 + nodeCount, 0, 1.0, B2);

            //A.SaveToTextFile("C:\\tmp\\AMatrix.txt");

            // test with matlab
            MultidimensionalArray output = MultidimensionalArray.Create(1, 2);
            Console.WriteLine("Calling MATLAB/Octave...");
            using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                bmc.PutSparseMatrix(A, "A");
                bmc.Cmd("rank_A = rank(full(A))");
                bmc.Cmd("rank_AT = rank(full(A'))");
                bmc.GetMatrix(output, "[rank_A, rank_AT]");

                bmc.Execute(false);
            }

            Console.WriteLine("A: No of Rows = {0}; rank = {1}", A.NoOfRows, output[0, 0]);
            Console.WriteLine("AT: No of Rows = {0}; rank = {1}", A.Transpose().NoOfRows, output[0, 1]);


            //var r = IMatrixExtensions.ReducedRowEchelonForm(Amtx);
            //Console.WriteLine("rank of Amtx = {0}", r.Item4);
            //Console.WriteLine("No of rows of reduced row-echelon form = {0}", r.Item1.Lengths[0]);
            //Console.WriteLine("No of columns of reduced row-echelon form = {0}", r.Item1.Lengths[1]);


            // solve
            MsrMatrix AAT = A * A.Transpose();

            double condNum = AAT.condest();
            Console.WriteLine("Condition Number of AAT is " + condNum);

            double[] RHS = new double[rowPart.LocalLength];
            A.SpMVpara(1.0, m_Coordinates.To1DArray(), 0.0, RHS);

            double[] v = new double[rowPart.LocalLength];
            double[] x = new double[m_Coordinates.Length];

            //OpSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            OpSolver.DefineMatrix(AAT);
            OpSolver.Solve(v, RHS);
            OpSolver.Dispose();

            A.Transpose().SpMVpara(-1.0, v, 0.0, x);

            m_Coordinates.AccVector(1.0, x);

        }


        /// <summary>
        /// Projects some DG field <paramref name="DGField"/> onto the internal, continuous representation
        /// </summary>
        /// <param name="DGField">
        /// input; unchanged on exit
        /// </param>
        /// <param name="mask"></param>
        public void ProjectDGField(ConventionalDGField DGField, CellMask mask = null) {

            bool useRowEchelonForm = false;

            if (useRowEchelonForm)
                ProjectDGField_rowEchelonForm(DGField, mask);
            else
                ProjectDGField_geometric(DGField, mask);

        }


        public void ProjectDGField_geometric(ConventionalDGField DGField, CellMask mask = null) {
            if (DGField.Basis.Degree > this.m_Basis.Degree)
                throw new ArgumentException("continuous projection on a lower degree basis is not recommended");
            this.Coordinates.Clear(); // clear internal state, to get the same result for the same input every time

            int D = this.Basis.GridDat.SpatialDimension;

            if (mask == null) {
                mask = CellMask.GetFullMask(m_grd);
            }
            if (mask.NoOfItemsLocally.MPISum() <= 0) {
                throw new ArgumentOutOfRangeException("Domain mask cannot be empty.");
            }

            // list of masked vertices inside domain mask
            List<GeometricVerticeForProjection> maskedVert = new List<GeometricVerticeForProjection>();
            List<GeometricEdgeForProjection> maskedEdges = new List<GeometricEdgeForProjection>();
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int[] vertAtCell = m_grd.Cells.CellVertices[j];
                    foreach (int vert in vertAtCell) {
                        GeometricVerticeForProjection gVert = new GeometricVerticeForProjection(vert);
                        if (!maskedVert.Contains(gVert)) {
                            maskedVert.Add(gVert);
                        }
                    }
                    List<GeometricEdgeForProjection> edgesAtCell = GetGeometricEdgesForCell(vertAtCell);
                    foreach (var gEdge in edgesAtCell) {
                        if (!maskedEdges.Contains(gEdge)) {
                            maskedEdges.Add(gEdge);
                        }
                    }
                }
            }


            // get DG-coordinates (change of basis for projection on a higher polynomial degree)
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int N = DGField.Basis.GetLength(j);
                    for (int n = 0; n < N; n++) {
                        m_Coordinates[m_Mapping.LocalUniqueCoordinateIndex(0, j, n)] = DGField.Coordinates[j, n];
                    }
                }
            }


            // construction of constraints matrix A
            // ====================================

            SubGrid maskSG = new SubGrid(mask);
            EdgeMask innerEM = maskSG.InnerEdgesMask;

            if (innerEM.NoOfItemsLocally.MPISum() <= 0) {
                Console.WriteLine("no inner edges: return without changes");
                return;
            }

            //List<int> AcceptedEdges = new List<int>();
            List<NodeSet> AcceptedNodes = new List<NodeSet>();

            int nodeCount = 0;


            // local edges per process
            foreach (var chunk in innerEM) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {

                    var edgeInfo = m_grd.Edges.Info[j];

                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];

                    int trafoIdx1 = m_grd.Edges.Edge2CellTrafoIndex[j, 0];
                    int trafoIdx2 = m_grd.Edges.Edge2CellTrafoIndex[j, 1];

                    int[] vertAtCell1 = m_grd.Cells.CellVertices[cell1];
                    int[] vertAtCell2 = m_grd.Cells.CellVertices[cell2];

                    // get geometric vertices/edges(3d) at considered edge/(face)
                    //List<GeometricVerticeForProjection> geomVertAtEdge = new List<GeometricVerticeForProjection>();
                    //List<GeometricEdgeForProjection> geomEdgeAtEdge = new List<GeometricEdgeForProjection>();


                    // 
                    int maxCondAtVert = (D == 2) ? 4 : 12;
                    int OverdeterminedCondAtVertice = 0;
                    //for (int i = 0; i < vertAtCell1.Length; i++) {
                    //    int vert = vertAtCell1[i];
                    //    if (vertAtCell2.Contains(vert)) {
                    //        GeometricVerticeForProjection gVert = maskedVert.Find(vrt => vrt.Equals(vert));
                    //        //geomVertAtEdge.Add(gVert);
                    //        gVert.IncreaseNoOfConditions();
                    //        int condAtVert = gVert.GetNoOfConditions();
                    //        //Console.WriteLine("conditions at vertice {0}: {1}", vert, condAtVert);
                    //        Debug.Assert(condAtVert <= maxCondAtVert);
                    //        if (condAtVert == maxCondAtVert)
                    //            OverdeterminedCondAtVertice++;
                    //    }
                    //}
                    ////Debug.Assert(geomVertAtEdge.Count() == (D - 1) * 2);

                    // 
                    int OverdeterminedCondAtGeomEdge = 0;
                    //int[] OverdeterminedEdgeDirection = new int[m_Basis.Degree + 1];
                    //if (D == 3) {
                    //    List<GeometricEdgeForProjection> edgesAtFace1 = GetGeometricEdgesForCell(vertAtCell1);
                    //    List<GeometricEdgeForProjection> edgesAtFace2 = GetGeometricEdgesForCell(vertAtCell2);
                    //    foreach (var gEdge1 in edgesAtFace1) {
                    //        if (edgesAtFace2.Contains(gEdge1)) {
                    //            GeometricEdgeForProjection gEdge = maskedEdges.Find(edg => edg.Equals(gEdge1));
                    //            //geomEdgeAtEdge.Add(gEdge);
                    //            gEdge.IncreaseNoOfConditions();
                    //            int condAtEdge = gEdge.GetNoOfConditions();
                    //            //Console.WriteLine("conditions at edge ({0}/{1}): {2}", gEdge1.VerticeInd1, gEdge1.VerticeInd2, condAtEdge);
                    //            Debug.Assert(condAtEdge <= 4);
                    //            if (condAtEdge == 4) {
                    //                //int dir1 = gEdge.GetRefDirection(m_grd, cell1, trafoIdx1);
                    //                //int dir2 = gEdge.GetRefDirection(m_grd, cell2, trafoIdx2);
                    //                //if (dir1 != dir2)
                    //                //    throw new ArgumentException("constrainedDG field: dir1 != dir2");
                    //                OverdeterminedEdgeDirection[OverdeterminedCondAtGeomEdge] = gEdge.GetRefDirection(m_grd, cell1, trafoIdx1);
                    //                OverdeterminedCondAtGeomEdge++;
                    //            }
                    //        }
                    //    }
                    //}


                    //if (!isInterprocEdge) 
                    {
                        //AcceptedEdges.Add(j);
                        //NodeSet qNds = getEdgeInterpolationNodes(OverdeterminedCondAtVertice, OverdeterminedCondAtGeomEdge, OverdeterminedEdgeDirection);
                        NodeSet qNds = getEdgeInterpolationNodes(0, 0);
                        AcceptedNodes.Add(qNds);

                        if (qNds != null) {
                            nodeCount += qNds.NoOfNodes;

                            //Console.WriteLine("continuity at edge {0}: numVcond = {1}, numEcond = {2}, nodeCount = {3}",
                            //    j, OverdeterminedCondAtVertice, OverdeterminedCondAtGeomEdge, qNds.NoOfNodes);
                        } else {
                            //Console.WriteLine("no continuity at edge {0}: numVcond = {1}, numEcond = {2}",
                            //    j, OverdeterminedCondAtVertice, OverdeterminedCondAtGeomEdge);
                        }

                    }

                }
            }
            //Console.WriteLine("node count: {0}", nodeCount);

            //NodeSet qNodes = getEdgeInterpolationNodes(0,0);
            //int nodeCount = innerEM.NoOfItemsLocally * qNodes.Lengths[0];

            List<long> BlockI0 = new List<long>();
            List<int> BlockLen = new List<int>();
            long i0 = 0;
            //for (int i = 0; i < innerEM.NoOfItemsLocally; i++) {
            foreach (NodeSet ns in AcceptedNodes) {
                if (ns != null) {
                    BlockI0.Add(i0);
                    BlockLen.Add(ns.Lengths[0]);
                    i0 += ns.Lengths[0];
                }
            }

            BlockPartitioning rowBlockPart = new BlockPartitioning(nodeCount, BlockI0, BlockLen, m_Mapping.MPI_Comm, true);
            BlockMsrMatrix A = new BlockMsrMatrix(rowBlockPart, m_Mapping);

            //MultidimensionalArray Amtx = MultidimensionalArray.Create(rowBlockPart.LocalLength, m_Mapping.LocalLength);

            //MsrMatrix AAT = new MsrMatrix(rowPart, rowPart);

            //Console.WriteLine("rank {0}: construct constraint matrix A", this.m_grd.MpiRank);
            //Console.WriteLine("No of nodes on rank {0}: {1}", this.m_grd.MpiRank, nodeCount);

            int count = 0;
            long _nodeCount = A.RowPartitioning.i0; // start at global index
            foreach (var chunk in innerEM) {
                foreach (int j in chunk.Elements) {

                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];

                    // set continuity constraints
                    NodeSet qNodes = AcceptedNodes.ElementAt(count);
                    if (qNodes == null)
                        break;

                    var results = m_Basis.EdgeEval(qNodes, j, 1);

                    //int eToCell1 = m_grd.Edges.Edge2CellTrafoIndex[j, 0];
                    //int eToCell2 = m_grd.Edges.Edge2CellTrafoIndex[j, 1];

                    //// results cell 1
                    //var trafo = m_grd.Edges.Edge2CellTrafos[eToCell1];
                    //MultidimensionalArray qCellNodes = trafo.Transform(qNodes);
                    //NodeSet qCellNdset = new NodeSet(m_grd.Cells.RefElements[0], qCellNodes);
                    //MultidimensionalArray results1 = m_Basis.Evaluate(qCellNdset);
                    //results1.Scale(m_grd.Cells.JacobiDet[cell1]);

                    //// results cell 2
                    //trafo = m_grd.Edges.Edge2CellTrafos[eToCell2];
                    //qCellNodes = trafo.Transform(qNodes);
                    //qCellNdset = new NodeSet(m_grd.Cells.RefElements[0], qCellNodes);
                    //MultidimensionalArray results2 = m_Basis.Evaluate(qCellNdset);
                    //results2.Scale(m_grd.Cells.JacobiDet[cell2]);

                    for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                        // Cell1   
                        for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                            A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                            //A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results1[0, qN, p];
                        }
                        // Cell2
                        for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                            A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
                            //A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results2[0, qN, p];
                        }
                    }
                    count++;
                    _nodeCount += qNodes.NoOfNodes;
                }
            }

            //Console.WriteLine("rank {0}: finished assembly", this.m_grd.MpiRank);

            //Partitioning rowPart = new Partitioning(nodeCount);
            //MsrMatrix A = new MsrMatrix(rowPart, m_Mapping);
            //A.AccBlock(rowPart.i0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, (int)m_Mapping.GlobalCount - 1 }));
            //MsrMatrix A = new MsrMatrix(nodeCount, m_Coordinates.Length, 1, 1);
            //A.AccBlock(0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, m_Coordinates.Length - 1 }));
            //A.AccBlock(rowPart.i0 + nodeCount, 0, 1.0, B2);

            //A.SaveToTextFileSparse("C:\\tmp\\AMatrix.txt");
            //Console.WriteLine("Writing A-matrix to text file!!!");

            //// test with matlab
            //MultidimensionalArray output = MultidimensionalArray.Create(1, 2);
            //Console.WriteLine("Calling MATLAB/Octave...");
            //using (BatchmodeConnector bmc = new BatchmodeConnector()) {
            //    bmc.PutSparseMatrix(A, "A");
            //    bmc.Cmd("rank_A = rank(full(A))");
            //    //bmc.Cmd("rank_AT = rank(full(A'))");
            //    //bmc.GetMatrix(output, "[rank_A, rank_AT]");
            //    bmc.GetMatrix(output, "[rank_A, 0]");

            //    bmc.Execute(false);
            //}

            //Console.WriteLine("A: No of Rows = {0}; rank = {1} (MATLAB)", A.NoOfRows, output[0, 0]);
            ////Console.WriteLine("AT: No of Rows = {0}; rank = {1}", A.Transpose().NoOfRows, output[0, 1]);

            // solve
            BlockMsrMatrix AT = A.Transpose();

            BlockMsrMatrix AAT = new BlockMsrMatrix(rowBlockPart, rowBlockPart);
            BlockMsrMatrix.Multiply(AAT, A, AT);

            //Console.WriteLine("A: No of Rows = {0}", A.NoOfRows);
            //double condNum = AAT.condest();
            //Console.WriteLine("Condition Number of AAT is " + condNum);

            double[] RHS = new double[rowBlockPart.LocalLength];
            //A.SpMVpara(1.0, m_Coordinates.To1DArray(), 0.0, RHS);
            A.SpMV(1.0, m_Coordinates.To1DArray(), 0.0, RHS);

            //m_Coordinates.To1DArray().SaveToTextFile("C:\\tmp\\Coord.txt");
            //Console.WriteLine("Writing coordinates to text file!!!");
            ////AAT.SaveToTextFileSparse("C:\\tmp\\AATMatrix.txt");
            ////Console.WriteLine("Writing AAT-matrix to text file!!!");
            //RHS.SaveToTextFile("C:\\tmp\\RHS.txt");
            //Console.WriteLine("Writing RHS to text file!!!");

            double[] v = new double[rowBlockPart.LocalLength];
            double[] x = new double[m_Coordinates.Length];

            //OpSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            OpSolver.DefineMatrix(AAT);
            //Console.WriteLine("rank {0}: solve constraint variables", this.m_grd.MpiRank);
            OpSolver.Solve(v, RHS);
            //Console.WriteLine("rank {0}: done", this.m_grd.MpiRank);
            OpSolver.Dispose();

            //A.Transpose().SpMVpara(-1.0, v, 0.0, x);
            AT.SpMV(-1.0, v, 0.0, x);

            m_Coordinates.AccVector(1.0, x);


            //// solve with matlab
            //double[] RHS = new double[rowBlockPart.LocalLength];
            //A.SpMV(1.0, m_Coordinates.To1DArray(), 0.0, RHS);
            //MultidimensionalArray v = MultidimensionalArray.Create(rowBlockPart.LocalLength, 1);

            //Console.WriteLine("Calling MATLAB/Octave...");
            //using (BatchmodeConnector bmc = new BatchmodeConnector()) {
            //    bmc.PutSparseMatrix(A, "A");
            //    bmc.PutVector(RHS, "b");
            //    bmc.Cmd("AAT = A * (A.');");
            //    bmc.Cmd("[L, U] = ilu(AAT);");
            //    bmc.Cmd("v = pcg(AAT, b, 1e-8, 500, L, U);");
            //    bmc.GetMatrix(v, "v");

            //    bmc.Execute(false);
            //}
            //Console.WriteLine("Closing MATLAB/Octave.");

            //double[] x = new double[m_Coordinates.Length];

            //BlockMsrMatrix AT = A.Transpose();
            //AT.SpMV(-1.0, v.ExtractSubArrayShallow(new int[] { -1, 0 }).To1DArray(), 0.0, x);

            //m_Coordinates.AccVector(1.0, x);


            //// row echelon form
            //var r = IMatrixExtensions.ReducedRowEchelonForm(Amtx);
            //Console.WriteLine("rank of A = {0} (reduced row-echelon form)", r.Item4);

            //MultidimensionalArray rMtx = MultidimensionalArray.Create(r.Item2.Length, m_Mapping.LocalLength);
            //rMtx.SetMatrix(r.Item1.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { r.Item4 - 1, m_Mapping.LocalLength - 1 }));

            //MsrMatrix R = rMtx.ToMsrMatrix();
            //MsrMatrix RT = R.Transpose();
            //MsrMatrix RRT = R * RT;

            //double condNum = RRT.condest();
            //Console.WriteLine("Condition Number of RRT is " + condNum);

            //double[] RHS = new double[r.Item2.Length];
            ////A.SpMVpara(1.0, m_Coordinates.To1DArray(), 0.0, RHS);
            //R.SpMV(1.0, m_Coordinates.To1DArray(), 0.0, RHS);

            //double[] v = new double[r.Item2.Length];
            //double[] x = new double[m_Coordinates.Length];

            ////OpSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            //OpSolver.DefineMatrix(RRT);
            ////Console.WriteLine("rank {0}: solve constraint variables", this.m_grd.MpiRank);
            //OpSolver.Solve(v, RHS);
            ////Console.WriteLine("rank {0}: done", this.m_grd.MpiRank);
            //OpSolver.Dispose();

            ////A.Transpose().SpMVpara(-1.0, v, 0.0, x);
            //RT.SpMV(-1.0, v, 0.0, x);

            //m_Coordinates.AccVector(1.0, x);

        }


        public void ProjectDGField_rowEchelonForm(ConventionalDGField DGField, CellMask mask = null) {
            if (DGField.Basis.Degree > this.m_Basis.Degree)
                throw new ArgumentException("continuous projection on a lower degree basis is not recommended");
            this.Coordinates.Clear(); // clear internal state, to get the same result for the same input every time

            int D = this.Basis.GridDat.SpatialDimension;

            if (mask == null) {
                mask = CellMask.GetFullMask(m_grd);
            }
            if (mask.NoOfItemsLocally.MPISum() <= 0) {
                throw new ArgumentOutOfRangeException("Domain mask cannot be empty.");
            }

            // list of masked vertices inside domain mask
            List<GeometricVerticeForProjection> maskedVert = new List<GeometricVerticeForProjection>();
            List<GeometricEdgeForProjection> maskedEdges = new List<GeometricEdgeForProjection>();
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int[] vertAtCell = m_grd.Cells.CellVertices[j];
                    foreach (int vert in vertAtCell) {
                        GeometricVerticeForProjection gVert = new GeometricVerticeForProjection(vert);
                        if (!maskedVert.Contains(gVert)) {
                            maskedVert.Add(gVert);
                        }
                    }
                    List<GeometricEdgeForProjection> edgesAtCell = GetGeometricEdgesForCell(vertAtCell);
                    foreach (var gEdge in edgesAtCell) {
                        if (!maskedEdges.Contains(gEdge)) {
                            maskedEdges.Add(gEdge);
                        }
                    }
                }
            }


            // get DG-coordinates (change of basis for projection on a higher polynomial degree)
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int N = DGField.Basis.GetLength(j);
                    for (int n = 0; n < N; n++) {
                        m_Coordinates[m_Mapping.LocalUniqueCoordinateIndex(0, j, n)] = DGField.Coordinates[j, n];
                    }
                }
            }


            // construction of constraints matrix A
            // ====================================

            SubGrid maskSG = new SubGrid(mask);
            EdgeMask innerEM = maskSG.InnerEdgesMask;

            if (innerEM.NoOfItemsLocally.MPISum() <= 0) {
                Console.WriteLine("no inner edges: return without changes");
                return;
            }

            List<int> edgeList = new List<int>();
            List<NodeSet> AcceptedNodes = new List<NodeSet>();

            int edgeCount = 0;
            int nodeCount = 0;

            int OverdeterminedCondAtGeomEdge = 0;

            // local edges per process
            foreach (var chunk in innerEM) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {

                    var edgeInfo = m_grd.Edges.Info[j];

                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];

                    int[] vertAtCell1 = m_grd.Cells.CellVertices[cell1];
                    int[] vertAtCell2 = m_grd.Cells.CellVertices[cell2];

                    // get geometric vertices/edges(3d) at considered edge/(face)
                    //List<GeometricVerticeForProjection> geomVertAtEdge = new List<GeometricVerticeForProjection>();
                    //List<GeometricEdgeForProjection> geomEdgeAtEdge = new List<GeometricEdgeForProjection>();

                    edgeList.Add(j);
                    edgeCount++;

                    // 
                    int maxCondAtVert = (D == 2) ? 4 : 12;
                    int OverdeterminedCondAtVertice = 0;
                    for (int i = 0; i < vertAtCell1.Length; i++) {
                        int vert = vertAtCell1[i];
                        if (vertAtCell2.Contains(vert)) {
                            GeometricVerticeForProjection gVert = maskedVert.Find(vrt => vrt.Equals(vert));
                            //geomVertAtEdge.Add(gVert);
                            gVert.IncreaseNoOfConditions();
                            int condAtVert = gVert.GetNoOfConditions();
                            //Console.WriteLine("conditions at vertice {0}: {1}", vert, condAtVert);
                            Debug.Assert(condAtVert <= maxCondAtVert);
                            if (condAtVert == maxCondAtVert)
                                OverdeterminedCondAtVertice++;
                        }
                    }
                    //Debug.Assert(geomVertAtEdge.Count() == (D - 1) * 2);

                    // 
                    //int OverdeterminedCondAtGeomEdge = 0;
                    if (D == 3) {
                        List<GeometricEdgeForProjection> edgesAtFace1 = GetGeometricEdgesForCell(vertAtCell1);
                        List<GeometricEdgeForProjection> edgesAtFace2 = GetGeometricEdgesForCell(vertAtCell2);
                        foreach (var gEdge1 in edgesAtFace1) {
                            if (edgesAtFace2.Contains(gEdge1)) {
                                GeometricEdgeForProjection gEdge = maskedEdges.Find(edg => edg.Equals(gEdge1));
                                gEdge.AddEdge(j);
                                gEdge.IncreaseNoOfConditions();
                                int condAtEdge = gEdge.GetNoOfConditions();
                                //Console.WriteLine("conditions at edge ({0}/{1}): {2}", gEdge1.VerticeInd1, gEdge1.VerticeInd2, condAtEdge);
                                Debug.Assert(condAtEdge <= 4);
                                if (condAtEdge == 4)
                                    OverdeterminedCondAtGeomEdge++;
                            }
                        }
                    }


                    //if (!isInterprocEdge) 
                    {
                        //AcceptedEdges.Add(j);
                        //NodeSet qNds = getEdgeInterpolationNodes(OverdeterminedCondAtVertice, OverdeterminedCondAtGeomEdge);
                        NodeSet qNds = getEdgeInterpolationNodes(0, 0);
                        AcceptedNodes.Add(qNds);

                        if (qNds != null) {
                            nodeCount += qNds.NoOfNodes;

                            //Console.WriteLine("continuity at edge {0}: numVcond = {1}, numEcond = {2}, nodeCount = {3}",
                            //    j, OverdeterminedCondAtVertice, OverdeterminedCondAtGeomEdge, qNds.NoOfNodes);
                        } else {
                            //Console.WriteLine("no continuity at edge {0}: numVcond = {1}, numEcond = {2}",
                            //    j, OverdeterminedCondAtVertice, OverdeterminedCondAtGeomEdge);
                        }

                    }

                }
            }

            NodeSet qNodes = getEdgeInterpolationNodes(0,0);
            //int nodeCount = innerEM.NoOfItemsLocally * qNodes.Lengths[0];

            int NoConditions = OverdeterminedCondAtGeomEdge * 4 * qNodes.NoOfNodes;
            MultidimensionalArray rMtx = MultidimensionalArray.Create(NoConditions, m_Mapping.LocalLength);

            int _nodeCount = 0;
            foreach (var gEdge in maskedEdges) {
                List<int> blockList = gEdge.GetEdgeList();
                if (blockList.Count == 4) {
                    MultidimensionalArray Ablock = MultidimensionalArray.Create(4 * qNodes.NoOfNodes, m_Mapping.LocalLength);
                    int _nodeBlockCount = 0;
                    foreach (int jEdge in blockList) {

                        int cell1 = m_grd.Edges.CellIndices[jEdge, 0];
                        int cell2 = m_grd.Edges.CellIndices[jEdge, 1];

                        // set continuity constraints
                        var results = m_Basis.EdgeEval(qNodes, jEdge, 1);

                        for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                            // Cell1       
                            for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                                Ablock[(int)(_nodeBlockCount + qN), (int)m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                            }
                            // Cell2
                            for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                                Ablock[(int)(_nodeBlockCount + qN), (int)(m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p))] = -results.Item2[0, qN, p];
                            }
                        }
                        _nodeBlockCount += qNodes.NoOfNodes;
                    }

                    var rblock = IMatrixExtensions.ReducedRowEchelonForm(Ablock);
                    //MultidimensionalArray rblockMtx = MultidimensionalArray.Create(rblock.Item2.Length, m_Mapping.LocalLength);
                    //rblockMtx.SetMatrix(rblock.Item1.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { rblock.Item4 - 1, m_Mapping.LocalLength - 1 }));
                    rMtx.SetSubMatrix(rblock.Item1, new int[] { _nodeCount, 0 }, new int[] { _nodeCount + _nodeBlockCount - 1, m_Mapping.LocalLength - 1 });

                    _nodeCount += _nodeBlockCount;
                }
            }
            rMtx.SaveToTextFile("C:\\tmp\\rMtx.txt");


            List<long> BlockI0 = new List<long>();
            List<int> BlockLen = new List<int>();
            long i0 = 0;
            //for (int i = 0; i < innerEM.NoOfItemsLocally; i++) {
            foreach (NodeSet ns in AcceptedNodes) {
                if (ns != null) {
                    BlockI0.Add(i0);
                    BlockLen.Add(ns.Lengths[0]);
                    i0 += ns.Lengths[0];
                }
            }

            BlockPartitioning rowBlockPart = new BlockPartitioning(nodeCount, BlockI0, BlockLen, m_Mapping.MPI_Comm, true);
            BlockMsrMatrix A = new BlockMsrMatrix(rowBlockPart, m_Mapping);

            MultidimensionalArray Amtx = MultidimensionalArray.Create(rowBlockPart.LocalLength, m_Mapping.LocalLength);

            int count = 0;
            long _nodeCountLong = A.RowPartitioning.i0; // start at global index
            foreach (var chunk in innerEM) {
                foreach (int j in chunk.Elements) {

                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];

                    // set continuity constraints
                    //NodeSet qNodes = AcceptedNodes.ElementAt(count);
                    //if (qNodes == null)
                    //    break;

                    var results = m_Basis.EdgeEval(qNodes, j, 1);

                    for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                        // Cell1       
                        for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                            //A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                            Amtx[(int)(_nodeCountLong + qN), (int)m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                        }
                        // Cell2
                        for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                            //A[_nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
                            Amtx[(int)(_nodeCountLong + qN), (int)(m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p))] = -results.Item2[0, qN, p];
                        }
                    }
                    count++;
                    _nodeCountLong += qNodes.NoOfNodes;
                }
            }

            //Console.WriteLine("rank {0}: finished assembly", this.m_grd.MpiRank);

            //Partitioning rowPart = new Partitioning(nodeCount);
            //MsrMatrix A = new MsrMatrix(rowPart, m_Mapping);
            //A.AccBlock(rowPart.i0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, (int)m_Mapping.GlobalCount - 1 }));
            //MsrMatrix A = new MsrMatrix(nodeCount, m_Coordinates.Length, 1, 1);
            //A.AccBlock(0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, m_Coordinates.Length - 1 }));
            //A.AccBlock(rowPart.i0 + nodeCount, 0, 1.0, B2);

            //A.SaveToTextFile("C:\\tmp\\AMatrix.txt");

            //// test with matlab
            //MultidimensionalArray output = MultidimensionalArray.Create(1, 2);
            //Console.WriteLine("Calling MATLAB/Octave...");
            //using (BatchmodeConnector bmc = new BatchmodeConnector()) {
            //    bmc.PutSparseMatrix(A, "A");
            //    bmc.Cmd("rank_A = rank(full(A))");
            //    //bmc.Cmd("rank_AT = rank(full(A'))");
            //    //bmc.GetMatrix(output, "[rank_A, rank_AT]");
            //    bmc.GetMatrix(output, "[rank_A, 0]");

            //    bmc.Execute(false);
            //}

            //Console.WriteLine("A: No of Rows = {0}; rank = {1} (MATLAB)", A.NoOfRows, output[0, 0]);
            ////Console.WriteLine("AT: No of Rows = {0}; rank = {1}", A.Transpose().NoOfRows, output[0, 1]);


            //// solve
            //BlockMsrMatrix AT = A.Transpose();

            //BlockMsrMatrix AAT = new BlockMsrMatrix(rowBlockPart, rowBlockPart);
            //BlockMsrMatrix.Multiply(AAT, A, AT);

            //double condNum = AAT.condest();
            //Console.WriteLine("Condition Number of AAT is " + condNum);

            //double[] RHS = new double[rowBlockPart.LocalLength];
            ////A.SpMVpara(1.0, m_Coordinates.To1DArray(), 0.0, RHS);
            //A.SpMV(1.0, m_Coordinates.To1DArray(), 0.0, RHS);

            ////AAT.SaveToTextFileSparse("C:\\tmp\\AATMatrix.txt");
            ////Console.WriteLine("Writing AAT-matrix to text file!!!");
            ////RHS.SaveToTextFile("C:\\tmp\\RHS.txt");
            ////Console.WriteLine("Writing RHS to text file!!!");

            //double[] v = new double[rowBlockPart.LocalLength];
            //double[] x = new double[m_Coordinates.Length];

            ////OpSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            //OpSolver.DefineMatrix(AAT);
            ////Console.WriteLine("rank {0}: solve constraint variables", this.m_grd.MpiRank);
            //OpSolver.Solve(v, RHS);
            ////Console.WriteLine("rank {0}: done", this.m_grd.MpiRank);
            //OpSolver.Dispose();

            ////A.Transpose().SpMVpara(-1.0, v, 0.0, x);
            //AT.SpMV(-1.0, v, 0.0, x);

            //m_Coordinates.AccVector(1.0, x);


            // row echelon form
            var r = IMatrixExtensions.ReducedRowEchelonForm(Amtx);
            Console.WriteLine("rank of A = {0} (reduced row-echelon form)", r.Item4);

            rMtx = MultidimensionalArray.Create(r.Item2.Length, m_Mapping.LocalLength);
            rMtx.SetMatrix(r.Item1.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { r.Item4 - 1, m_Mapping.LocalLength - 1 }));

            MsrMatrix R = rMtx.ToMsrMatrix();
            MsrMatrix RT = R.Transpose();
            MsrMatrix RRT = R * RT;

            rMtx.SaveToTextFile("C:\\tmp\\RMatrix.txt");
            double condNum = RRT.condest();
            Console.WriteLine("Condition Number of RRT is " + condNum);

            double[] RHS = new double[NoConditions];
            //A.SpMVpara(1.0, m_Coordinates.To1DArray(), 0.0, RHS);
            R.SpMV(1.0, m_Coordinates.To1DArray(), 0.0, RHS);

            double[] v = new double[NoConditions];
            double[] x = new double[m_Coordinates.Length];

            //OpSolver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            OpSolver.DefineMatrix(RRT);
            //Console.WriteLine("rank {0}: solve constraint variables", this.m_grd.MpiRank);
            OpSolver.Solve(v, RHS);
            //Console.WriteLine("rank {0}: done", this.m_grd.MpiRank);
            OpSolver.Dispose();

            //A.Transpose().SpMVpara(-1.0, v, 0.0, x);
            RT.SpMV(-1.0, v, 0.0, x);

            m_Coordinates.AccVector(1.0, x);

        }




        class GeometricVerticeForProjection {

            int VerticeIndex;

            int NoOfConditions;

            public GeometricVerticeForProjection(int vertInd) {
                VerticeIndex = vertInd;
                NoOfConditions = 0;
            }

            public void IncreaseNoOfConditions() {
                this.NoOfConditions++;
            }

            public int GetNoOfConditions() {
                return this.NoOfConditions;
            }

            public override bool Equals(Object obj) {

                if (obj is GeometricVerticeForProjection) {
                    return (this.VerticeIndex - ((GeometricVerticeForProjection)obj).VerticeIndex) == 0;
                } else if (obj is int) {
                    return (this.VerticeIndex - (int)obj) == 0;
                } else
                    throw new ArgumentException("wrong type of object");
             
            }

        }

        class GeometricEdgeForProjection {

            public int VerticeInd1;
            public int VerticeInd2;

            List<int> edgeList;
            int NoOfConditions;

            public GeometricEdgeForProjection(int vertInd1, int vertInd2) {
                this.VerticeInd1 = vertInd1;
                this.VerticeInd2 = vertInd2;
                edgeList = new List<int>();
                NoOfConditions = 0;
            }

            public override bool Equals(object obj) {

                if (obj is GeometricEdgeForProjection) {

                    if ((this.VerticeInd1 == ((GeometricEdgeForProjection)obj).VerticeInd1 
                        && this.VerticeInd2 == ((GeometricEdgeForProjection)obj).VerticeInd2)
                        || (this.VerticeInd1 == ((GeometricEdgeForProjection)obj).VerticeInd2 
                        && this.VerticeInd2 == ((GeometricEdgeForProjection)obj).VerticeInd1))
                        return true;
                    else
                        return false;

                } else
                    throw new ArgumentException("wrong type of object");
            }

            public void AddEdge(int jEdge) {
                if (edgeList.Count == 4)
                    throw new InvalidOperationException("edgeList: maximum count of 4 elements reached");
                edgeList.Add(jEdge);
            }

            public List<int> GetEdgeList() {
                return edgeList;
            }

            public void IncreaseNoOfConditions() {
                this.NoOfConditions++;
            }

            public int GetNoOfConditions() {
                return this.NoOfConditions;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <returns></returns>
            public int GetRefDirection(GridData m_grd, int jCell, int iFace) {

                double[] vertCoord1 = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(this.VerticeInd1, -1).To1DArray();
                double[] vertCoord2 = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(this.VerticeInd2, -1).To1DArray();

                int D = vertCoord1.Length;
                MultidimensionalArray vertCoord1_glb = MultidimensionalArray.Create(1, D);
                MultidimensionalArray vertCoord2_glb = MultidimensionalArray.Create(1, D);
                vertCoord1_glb.SetRow<double[]>(0, vertCoord1);
                vertCoord2_glb.SetRow<double[]>(0, vertCoord2);

                MultidimensionalArray vertCoord1_loc = MultidimensionalArray.Create(1, 1, D);
                MultidimensionalArray vertCoord2_loc = MultidimensionalArray.Create(1, 1, D);

                m_grd.TransformGlobal2Local(vertCoord1_glb, vertCoord1_loc, jCell, 1, 0);
                m_grd.TransformGlobal2Local(vertCoord2_glb, vertCoord2_loc, jCell, 1, 0);
                double[] vertCellCoord1 = vertCoord1_loc.ExtractSubArrayShallow(0, 0, -1).To1DArray();
                double[] vertCellCoord2 = vertCoord2_loc.ExtractSubArrayShallow(0, 0, -1).To1DArray();

                // identify face ref vertices
                NodeSet refV = m_grd.Edges.EdgeRefElements[0].Vertices;
                NodeSet refVvol = m_grd.Edges.EdgeRefElements[0].Vertices.GetVolumeNodeSet(m_grd, iFace);
                //var trafo = m_grd.Edges.Edge2CellTrafos[iFace];
                int idx1 = -1; int idx2 = -1;
                for (int i = 0; i < refV.Lengths[0]; i++) {
                    double dist1 = 0.0; double dist2 = 0.0;
                    for (int d = 0; d < vertCoord1.Length; d++) {
                        dist1 += (vertCellCoord1[d] - refVvol[i,d]).Pow2();
                        dist2 += (vertCellCoord2[d] - refVvol[i,d]).Pow2();
                    }
                    dist1 = dist1.Sqrt();
                    if (dist1 < 1.0e-12)
                        idx1 = i;
                    dist2 = dist2.Sqrt();
                    if (dist2 < 1.0e-12)
                        idx2 = i;            
                }
                Debug.Assert(idx1 != idx2);

                // identify direction of edge in ref element
                double dist = 0.0;
                for (int d = 0; d < refV.Lengths[1]; d++) {
                    dist += (refV[idx1,d] - refV[idx2,d]).Pow2();
                }
                dist = dist.Sqrt();
                int dir = -1;
                for (int d = 0; d < refV.Lengths[1]; d++) {
                    if (Math.Abs(Math.Abs(refV[idx1, d] - refV[idx2, d]) - dist) < 1.0e-15) {
                        dir = d;
                    }
                }

                return dir;

            }

        }


        List<GeometricEdgeForProjection> GetGeometricEdgesForCell(int[] vertAtCell) {

            List<GeometricEdgeForProjection> geomEdges = new List<GeometricEdgeForProjection>();

            //int[] vertAtCell = m_grd.Cells.CellVertices[jCell];
            for (int v1 = 0; v1 < vertAtCell.Length; v1++) {
                double[] vertCoord1 = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(vertAtCell[v1], -1).To1DArray();
                for (int v2 = v1+1; v2 < vertAtCell.Length; v2++) {
                    double[] vertCoord2 = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(vertAtCell[v2], -1).To1DArray();
                    double dist = 0.0;
                    for (int d = 0; d < vertCoord1.Length; d++) {
                        dist += (vertCoord2[d] - vertCoord1[d]).Pow2();
                    }
                    dist = dist.Sqrt();
                    bool isEdge = false;
                    //int dir = -1;
                    for (int d = 0; d < vertCoord1.Length; d++) {
                        if (Math.Abs(Math.Abs(vertCoord2[d] - vertCoord1[d]) - dist) < 1.0e-15) {
                            isEdge = true;
                            //dir = d;
                        }
                    }
                    if (isEdge) {
                        GeometricEdgeForProjection gEdge = new GeometricEdgeForProjection(vertAtCell[v1], vertAtCell[v2]);
                        //Console.WriteLine("define edge ({0},{1}) in direction {2}", gEdge.VerticeInd1, gEdge.VerticeInd2, dir);
                        geomEdges.Add(gEdge);
                    }
                }
            }

            return geomEdges;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="numVCond"></param>
        /// <param name="numECond"></param>
        /// <returns></returns>
        private NodeSet getEdgeInterpolationNodes(int numVcond, int numEcond, int[] edgeOrientation = null) {

            int degree = m_Basis.Degree;

            switch (m_grd.SpatialDimension) {
                case 2: {
                    if (numVcond > degree) {
                        return null;
                    } else {
                        QuadRule quad = m_grd.Edges.EdgeRefElements[0].GetQuadratureRule((degree - numVcond) * 2);
                        return quad.Nodes;
                    }
                }
                case 3: {

                    QuadRule quad1D = m_grd.Edges.EdgeRefElements[0].FaceRefElement.GetQuadratureRule(degree * 2);
                    NodeSet qNodes = quad1D.Nodes;
                    int Nnds = ((degree + 1) * (degree + 2) / 2); 

                    int degreeR = degree - numEcond; 
                    int NoNdsR = ((degreeR + 1) * (degreeR + 2) / 2);
                    if (NoNdsR <= 0) {
                        return null;

                    } else {
                        if (edgeOrientation == null && numEcond > 0)
                            throw new ArgumentException();
                        if (edgeOrientation == null && numEcond < 1)
                            edgeOrientation = new int[degree + 1];

                        MultidimensionalArray nds = MultidimensionalArray.Create(Nnds, 2);
                        int node = 0;
                        int[] dirCount = new int[2];
                        for (int dirIdx = 0; dirIdx < edgeOrientation.Length; dirIdx++) {
                            int dir = edgeOrientation[dirIdx];
                            int m = dirCount[dir];
                            int n0 = dirCount[(dir == 0) ? 1 : 0];
                            int nL = quad1D.NoOfNodes - dirIdx;
                            for (int n = n0; n < n0 + nL; n++) {
                                if (dir == 0) {
                                    nds[node, 0] = qNodes[n, 0];
                                    nds[node, 1] = qNodes[m, 0];
                                }
                                if (dir == 1) {
                                    nds[node, 0] = qNodes[m, 0];
                                    nds[node, 1] = qNodes[n, 0];
                                }
                                node++;
                            }
                            dirCount[dir]++;
                        }
                        MultidimensionalArray ndsR = nds.ExtractSubArrayShallow(new int[] {Nnds - NoNdsR, 0 }, new int[] { Nnds - 1, 1 });
                        //Console.WriteLine("No ndsR = {0}", ndsR.Lengths[0]);
                        return new NodeSet(m_grd.Edges.EdgeRefElements[0], ndsR);
                    }
                }
                default:
                    throw new NotSupportedException("spatial dimension not supported");
            }

        }

        ///// <summary>
        ///// 
        ///// </summary>
        ///// <param name="basis"></param>
        ///// <returns></returns>
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
        /// <param name="DGField">output</param>
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


    public class myCG : IDisposable {
        public void Init(BlockMsrMatrix M) {
            m_Matrix = M;
            //var M_test = M.CloneAs();
            //M_test.Acc(-1.0, M.Transpose());
            //Console.WriteLine("Symm-test: " + M_test.InfNorm());
            //Console.WriteLine("Inf-Norm: " + M.InfNorm());
            //m_Matrix.SaveToTextFileSparse("M");
            BlockJacInit();
            ILUInit();
        }

        private BlockMsrMatrix ILU_M;
        public void ILUInit() {
            ILU_M = m_Matrix.CloneAs();
        }
        public void ILUSolve<P, Q>(P X, Q B)
            where P : IList<double>
            where Q : IList<double> {

            long n = ILU_M.RowPartitioning.LocalLength;
            double sum = 0;

            var tempMtx = ILU_M;

            // Zeros on diagonal elements because of saddle point structure
            //for (int bla = 0; bla < n; bla++) {
            //    if (ILU_M.GetDiagonalElement(bla) == 0)
            //        throw new Exception("One or more diagonal elements are zero, ILU cannot work");
            //    ILU_M.SetDiagonalElement(bla, 1);
            //}

            // ILU decomposition of matrix
            for (long k = 0; k < n - 1; k++) {
                for (long i = k; i < n; i++) {
                    if (tempMtx[i, k] == 0) { i = n; } else {
                        ILU_M[i, k] = ILU_M[i, k] / ILU_M[k, k];
                        for (long j = k + 1; j < n; j++) {
                            if (tempMtx[i, j] == 0) {
                                j = n;
                            } else {
                                ILU_M[i, j] = ILU_M[i, j] - ILU_M[i, k] * ILU_M[k, j];
                            }
                        }
                    }
                }
            }

            if (this.ILU_M.RowPartitioning.MpiSize != 1)
                throw new NotSupportedException();

            // find solution of Ly = b
            double[] y = new double[n];
            for (int i = 0; i < n; i++) {
                sum = 0;
                for (int k = 0; k < i; k++)
                    sum += ILU_M[i, k] * y[k];
                y[i] = B[i] - sum;
            }
            // find solution of Ux = y
            for (long i = n - 1; i >= 0; i--) {
                sum = 0;
                for (long k = i + 1; k < n; k++)
                    sum += ILU_M[i, k] * X[(int)k]; // index into Mtx and X must be different for more than 1 MPI process.
                X[(int)i] = (1 / ILU_M[i, i]) * (y[i] - sum);
            }

        }
        public void BlockJacInit() {
            BlockMsrMatrix M = m_Matrix;

            int L = M.RowPartitioning.LocalLength;

            Diag = new BlockMsrMatrix(M._RowPartitioning, M._ColPartitioning);
            invDiag = new BlockMsrMatrix(M._RowPartitioning, M._ColPartitioning);
            int Jloc = M._RowPartitioning.LocalNoOfBlocks;
            long j0 = M._RowPartitioning.FirstBlock;
            MultidimensionalArray temp = null;
            for (int j = 0; j < Jloc; j++) {
                long jBlock = j + j0;
                int Nblk = M._RowPartitioning.GetBlockLen(jBlock);
                long i0 = M._RowPartitioning.GetBlockI0(jBlock);

                if (temp == null || temp.NoOfCols != Nblk)
                    temp = MultidimensionalArray.Create(Nblk, Nblk);

                M.ReadBlock(i0, i0, temp);
                Diag.AccBlock(i0, i0, 1.0, temp, 0.0);
                temp.InvertInPlace();
                invDiag.AccBlock(i0, i0, 1.0, temp, 0.0);
            }
        }
        private double omega = 0.5;

        public void BlockJacSolve<U, V>(U xl, V bl)
            where U : IList<double>
            where V : IList<double> {
            int L = xl.Count;
            double[] ql = new double[L];

            double iter0_ResNorm = 0;

            for (int iIter = 0; iIter < 10; iIter++) {
                ql.SetV(bl);
                m_Matrix.SpMV(-1.0, xl, 1.0, ql);
                double ResNorm = ql.L2NormPow2().MPISum().Sqrt();

                if (iIter == 0) {
                    iter0_ResNorm = ResNorm;
                }

                Diag.SpMV(1.0, xl, 1.0, ql);
                invDiag.SpMV(omega, ql, 1.0 - omega, xl);

            }
        }

        BlockMsrMatrix m_Matrix;
        BlockMsrMatrix Diag;
        BlockMsrMatrix invDiag;

        public void Solve<Vec1, Vec2>(Vec1 _x, Vec2 _R)
            where Vec1 : IList<double>
            where Vec2 : IList<double> //
        {

            double[] x, R;
            if (_x is double[]) {
                x = _x as double[];
            } else {
                x = _x.ToArray();
            }
            if (_R is double[]) {
                R = _R as double[];
            } else {
                R = _R.ToArray();
            }

            int L = x.Length;

            double[] P = new double[L];

            //double[] R = rhs; // rhs is only needed once, so we can use it to store residuals
            double[] V = new double[L];
            double[] Z = new double[L];


            // compute P0, R0
            // ==============
            GenericBlas.dswap(L, x, 1, P, 1);
            m_Matrix.SpMV(-1.0, P, 1.0, R);

            GenericBlas.dswap(L, x, 1, P, 1);
            //if (Precond != null) {
            //    Precond.Solve(Z, R);
            //    P.SetV(Z);
            //} else {
            //P.SetV(R);
            //}
            BlockJacSolve(Z, R);
            P.SetV(Z);

            double alpha = R.InnerProd(P).MPISum();
            double alpha_0 = alpha;
            double ResNorm;
            var ResReal = new double[R.Length];
            var Xdummy = new double[R.Length];

            ResNorm = Math.Sqrt(alpha);
            double ResNorm0 = ResNorm;

            // iterate
            // =======
            for (int n = 1; true; n++) {

                if (n % 1 == 0) {
                    Xdummy.SetV(x);
                    ResReal.SetV(R);
                    m_Matrix.SpMV(-1.0, Xdummy, 1.0, ResReal);
                    Console.WriteLine("Res real at n" + n + ":" + ResReal.MPI_L2Norm());
                }

                //Console.WriteLine("ResNorm at n"+n+":"+ResNorm);
                if (ResNorm / ResNorm0 + ResNorm < 1E-8 || ResNorm < 1E-16 || n > 100) {
                    if (n > 1000) Console.WriteLine("maximum number of iterations reached. Solution maybe not converged.");
                    break;
                }

                if (Math.Abs(alpha) <= double.Epsilon) {
                    // numerical breakdown
                    break;
                }


                m_Matrix.SpMV(1.0, P, 0, V);
                double VxP = V.InnerProd(P).MPISum();
                if (double.IsNaN(VxP) || double.IsInfinity(VxP))
                    throw new ArithmeticException();
                double lambda = alpha / VxP;
                if (double.IsNaN(lambda) || double.IsInfinity(lambda))
                    throw new ArithmeticException();


                x.AccV(lambda, P);

                R.AccV(-lambda, V);

                //if (Precond != null) {
                //    Z.Clear();
                //    Precond.Solve(Z, R);
                //} else {
                //Z.SetV(R);
                //}
                Z.Clear();
                BlockJacSolve(Z, R);

                double alpha_neu = R.InnerProd(Z).MPISum();

                // compute residual norm
                ResNorm = R.L2NormPow2().MPISum().Sqrt();

                P.ScaleV(alpha_neu / alpha);
                P.AccV(1.0, Z);

                alpha = alpha_neu;
            }

            if (!object.ReferenceEquals(_x, x))
                _x.SetV(x);
            if (!object.ReferenceEquals(_R, R))
                _R.SetV(R);


            return;

        }

        public void Dispose() {
            m_Matrix.Clear();
        }

    }


}

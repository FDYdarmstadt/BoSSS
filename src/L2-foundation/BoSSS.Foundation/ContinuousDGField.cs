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


        int[] Cell2Coord;

        int[] mask2OpCoord;


        /// <summary>
        /// linear solver for the quadratic optimization problem, Opmatrix has to be defined! 
        /// </summary>
        //public ilPSP.LinSolvers.ISparseSolver OpSolver {
        //    get {
        //        if (m_OpSolver == null) {
        //            //var solver = new ilPSP.LinSolvers.monkey.CG();
        //            //solver.MatrixType = ilPSP.LinSolvers.monkey.MatrixType.Auto;
        //            //solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.CPU;
        //            //solver.Tolerance = 1.0e-12;
        //            //var solver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();
        //            var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
        //            m_OpSolver = solver;
        //        }

        //        return m_OpSolver;
        //    }
        //}


        public void ProjectDGField(double alpha, ConventionalDGField DGField, CellMask mask = null) {
            if (DGField.Basis.Degree > this.m_Basis.Degree)
                throw new ArgumentException("continuous projection on a lower degree basis is not recommended");

            if (mask == null) {
                mask = CellMask.GetFullMask(m_grd);
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
                    }
                }
            }

            // construction of constraints matrix A
            // ====================================

            SubGrid maskSG = new SubGrid(mask);
            EdgeMask innerEM = maskSG.InnerEdgesMask;

            // get interpolation points for the continuity constraints
            NodeSet qNodes = getEdgeInterpolationNodes();


            MultidimensionalArray B = MultidimensionalArray.Create(innerEM.NoOfItemsLocally * qNodes.NoOfNodes, (int)m_Mapping.GlobalCount);

            int nodeCount = 0;
            int NumVert = m_grd.Vertices.NoOfNodes4LocallyUpdatedCells;
            MultidimensionalArray CondAtVert = MultidimensionalArray.Create(NumVert, 3);    // 1. index: local vertice index / 2. index: 0 = actual cond, 1 = potential own cond, 2 = potential cond
            MultidimensionalArray CondIncidenceMatrix = MultidimensionalArray.Create(NumVert, NumVert, 3);  // upper triangular matrix in the first two indices
            int maxVCond = ((2 * m_grd.SpatialDimension) - 1) * (m_grd.SpatialDimension - 1) - (m_grd.SpatialDimension - 2);
            int[][] externalEdges = m_grd.Edges.EdgeSendLists;
            List<int> InterprocEdges = new List<int>();
            List<List<int>> VertAtInterprocEdges = new List<List<int>>();
            bool isInterprocEdge;

            // local edges per process
            foreach (var chunk in innerEM) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {


                    // check for interprocess edges to be processed in a second step
                    isInterprocEdge = false;
                    if (m_grd.Edges.Info[j].HasFlag(EdgeInfo.Interprocess)) {
                        isInterprocEdge = true;
                        Console.WriteLine("interproc edge");
                        InterprocEdges.Add(j);
                        foreach (int[] list in externalEdges) {
                            if (list != null && list.Contains(j)) {
                                InterprocEdges.Remove(j);
                            }
                        }
                    }



                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];

                    // check the conditions at vertices (actual, own, potential)
                    int[] vertAtCell1 = m_grd.Cells.CellVertices[cell1];
                    int[] vertAtCell2 = m_grd.Cells.CellVertices[cell2];
                    int numVCond = 0;
                    List<int> VertAtEdge = new List<int>();
                    for (int i = 0; i < vertAtCell1.Length; i++) {
                        int vert = vertAtCell1[i];
                        if (vertAtCell2.Contains(vert)) {
                            VertAtEdge.Add(vert);
                            CondAtVert[vert, 2] += 1;   // potential condition at vert
                            if (!isInterprocEdge) {
                                CondAtVert[vert, 0] += 1;   // actual condition at vert 
                                CondAtVert[vert, 1] += 1;   // potential own cond at vert (local edge)
                            } else {
                                if (InterprocEdges.Contains(j)) {
                                    CondAtVert[vert, 1] += 1;   // potential own cond at vert (interproc edge)
                                } 
                            }                           
                            if (CondAtVert[vert,0] > maxVCond) {
                                numVCond++;
                            }
                        }
                    }
                    if (isInterprocEdge && InterprocEdges.Contains(j)) {
                        VertAtInterprocEdges.Add(VertAtEdge);
                    }


                    // check for overdetermined edges (additional for 3D)
                    int numECond = 0;
                    int VertC = VertAtEdge.Count;
                    if (VertC == 4) {
                        int i0 = 0;
                        for (int indV1 = 0; indV1 < VertC; indV1++) {
                            i0++;
                            for (int indV2 = i0; indV2 < VertC; indV2++) {
                                int m, n;
                                if (VertAtEdge.ElementAt(indV1) <= VertAtEdge.ElementAt(indV2)) {
                                    m = VertAtEdge.ElementAt(indV1);
                                    n = VertAtEdge.ElementAt(indV2);
                                } else {
                                    m = VertAtEdge.ElementAt(indV2);
                                    n = VertAtEdge.ElementAt(indV1);
                                }
                                Console.WriteLine("m = {0}, n = {1}", m, n);
                                CondIncidenceMatrix[m, n, 2] += 1;  // potential condition at edge
                                Console.WriteLine("potential conditions = {0}", CondIncidenceMatrix[m, n, 2]);
                                if (!isInterprocEdge) {
                                    CondIncidenceMatrix[m, n, 0] += 1;   // actual condition at edge 
                                    CondIncidenceMatrix[m, n, 1] += 1;   // potential own cond at edge (local edge/face)
                                } else {
                                    if (InterprocEdges.Contains(j)) {
                                        CondIncidenceMatrix[m, n, 1] += 1;   // potential own cond at vert (interproc edge)
                                    }
                                }
                                if (CondIncidenceMatrix[m, n, 0] == 4) {
                                    numECond++;
                                }

                            }
                        }
                    }


                    if (!isInterprocEdge) {

                        qNodes = getEdgeInterpolationNodes(numVCond, numECond);

                        if (qNodes != null) {

                            // set continuity constraints
                            var results = m_Basis.EdgeEval(qNodes, j, 1);

                            for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                                // Cell1
                                for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                                    B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                                }
                                // Cell2
                                for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                                    B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
                                }
                            }
                            nodeCount += qNodes.NoOfNodes;
                        }
                    }

                }
            }


            // interprocess edges 
            int edge = 0;
            foreach (int j in InterprocEdges) {


                int cell1 = m_grd.Edges.CellIndices[j, 0];
                int cell2 = m_grd.Edges.CellIndices[j, 1];

                int numVCond = 0;
                int numECond = 0;
                List<int> VertAtEdge = VertAtInterprocEdges.ElementAt(edge);
                if (VertAtInterprocEdges.First().Count == 2) {
                    // condtions at vertices 
                    numVCond = NegotiateNumVCond(VertAtEdge, CondAtVert);

                    foreach (int vert in VertAtEdge) {
                        CondAtVert[vert, 0] += 1;
                    }

                } else {
                    // conditions at edges 
                    numECond = NegotiateNumECond(VertAtEdge, CondIncidenceMatrix);

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
                            CondIncidenceMatrix[m, n, 0] += 1;
                        }
                    }
                }

                qNodes = getEdgeInterpolationNodes(numVCond, numECond);

                if (qNodes != null) {

                    // set continuity constraints
                    var results = m_Basis.EdgeEval(qNodes, j, 1);

                    for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                        // Cell1
                        for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                            B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell1, p)] = results.Item1[0, qN, p];
                        }
                        // Cell2
                        for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                            B[nodeCount + qN, m_Mapping.GlobalUniqueCoordinateIndex(0, cell2, p)] = -results.Item2[0, qN, p];
                        }
                    }
                    nodeCount += qNodes.NoOfNodes;
                }
                edge++;

            }


            Partitioning rowPart = new Partitioning(nodeCount);
            MsrMatrix A = new MsrMatrix(rowPart, m_Mapping);
            A.AccBlock(rowPart.i0, 0, 1.0, B.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { nodeCount - 1, (int)m_Mapping.GlobalCount - 1 }));


            // test with matlab
            MultidimensionalArray output = MultidimensionalArray.Create(1, 2);
            Console.WriteLine("Calling MATLAB/Octave...");
            using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                bmc.PutSparseMatrix(A, "A");
                bmc.Cmd("[m,n] = size(full(A))");
                bmc.Cmd("rank_A = rank(full(A))");
                bmc.Cmd("rank_AT = rank(full(A'))");
                bmc.GetMatrix(output, "[rank_A, rank_AT]");

                bmc.Execute(false);
            }

            Console.WriteLine("A: No of Rows = {0}; rank = {1}", A.NoOfRows, output[0, 0]);
            Console.WriteLine("AT: No of Rows = {0}; rank = {1}", A.Transpose().NoOfRows, output[0, 1]);


            // solve
            MsrMatrix AAT = A * A.Transpose();

            double[] RHS = new double[rowPart.LocalLength];
            A.SpMVpara(1.0, m_Coordinates.To1DArray(), 0.0, RHS);

            double[] v = new double[rowPart.LocalLength];
            double[] x = new double[m_Coordinates.Length];

            var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();

            solver.DefineMatrix(AAT);
            solver.Solve(v, RHS);

            A.Transpose().SpMVpara(-1.0, v, 0.0, x);

            m_Coordinates.AccVector(1.0, x);

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="numVCond"></param>
        /// <param name="numECond"></param>
        /// <returns></returns>
        private NodeSet getEdgeInterpolationNodes(int numVCond = 0, int numECond = 0) {

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
                                for (int n2 = 0; n2 <= n1; n2++) {
                                    nds[node, 0] = qNodes[n1, 0];
                                    nds[node, 1] = qNodes[n2, 0];
                                    node++;
                                }
                            }
                            return new NodeSet(m_grd.Edges.EdgeRefElements[0], nds); //.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Nnds - (1 + numVCond), 1 }));
                        }
                    }
                default:
                    throw new NotSupportedException("spatial dimension not supported");
            }

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="vertAtEdge"></param>
        /// <param name="CondAtVert"></param>
        /// <returns></returns>
        private int NegotiateNumVCond(List<int> VertAtEdge, MultidimensionalArray CondAtVert) {

            int numVCond = 0;

            foreach (int vert in VertAtEdge) {
                if (CondAtVert[vert, 2] == 3) {
                    if (CondAtVert[vert, 1] == 3) {
                        if (CondAtVert[vert, 0] == 2) {
                            numVCond++;
                        }
                    } else if (CondAtVert[vert, 1] == 2) {
                        if (vert < m_grd.Vertices.NodePartitioning.LocalLength) {   // proc owns vert
                            numVCond++;
                        }
                    }
                } 
            }

            return numVCond;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="vertAtEdge"></param>
        /// <param name="CondIncidenceMatrix"></param>
        /// <returns></returns>
        private int NegotiateNumECond(List<int> VertAtEdge, MultidimensionalArray CondIncidenceMatrix) {

            Console.WriteLine("negotiate number of edge conditons");

            int numECond = 0;

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
                    Console.WriteLine("m = {0}, n = {1}", m, n);
                    if (CondIncidenceMatrix[m, n, 2] == 3) {
                        Console.WriteLine("potential conditions = 3");
                        if (CondIncidenceMatrix[m, n, 1] == 3) {
                            Console.WriteLine("potential own conditions = 3");
                            if (CondIncidenceMatrix[m, n, 0] == 2) {
                                Console.WriteLine("actual conditions = 2");
                                numECond++;
                            }
                        } else if (CondIncidenceMatrix[m, n, 1] == 2) {
                            Console.WriteLine("potential own conditions = 2");
                            if (m < m_grd.Vertices.NodePartitioning.LocalLength && n < m_grd.Vertices.NodePartitioning.LocalLength) {   // proc owns both vert
                                Console.WriteLine("proc owns both vertices");
                                numECond++;
                            } else if (m < m_grd.Vertices.NodePartitioning.LocalLength && n >= m_grd.Vertices.NodePartitioning.LocalLength) {
                                int proc = 0;
                                foreach (int[] list in m_grd.Vertices.VertexInsertLists) {
                                    if (list != null && list.Contains(n)) {
                                        if (m_grd.MpiRank < proc) {
                                            Console.WriteLine("edge negotiated to this proc 0");
                                            numECond++;
                                        }
                                    }
                                    proc++;
                                }
                            } else if (m >= m_grd.Vertices.NodePartitioning.LocalLength && n < m_grd.Vertices.NodePartitioning.LocalLength) {
                                int proc = 0;
                                foreach (int[] list in m_grd.Vertices.VertexInsertLists) {
                                    if (list != null && list.Contains(m)) {
                                        if (m_grd.MpiRank < proc) {
                                            Console.WriteLine("edge negotiated to this proc 0");
                                            numECond++;
                                        }
                                    }
                                    proc++;
                                }
                            }
                        }
                    }

                }
            }

            return numECond;
        }


        /// <summary>
        /// project some conventional DG field onto this, this field will be continuous 
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="DGField"></param>
        /// <param name="mask"></param>
        public void ProjectDGField_serial(double alpha, ConventionalDGField DGField, CellMask mask = null) {
            if (DGField.Basis.Degree > this.m_Basis.Degree)
                throw new ArgumentException("continuous projection on a lower degree basis is not recommended");

            if (mask == null) {
                mask = CellMask.GetFullMask(m_grd);
            }


            // init (change of basis for projection on a higher degree)
            int J = m_grd.Cells.NoOfLocalUpdatedCells;
            Cell2Coord = new int[J];
            int sum = 0;
            for (int j = 0; j < J; j++) {
                Cell2Coord[j] = sum;
                sum += m_Basis.GetLength(j);
            }
            m_Coordinates = MultidimensionalArray.Create(sum);

            int Jmask = 0;
            mask2OpCoord = new int[J];
            mask2OpCoord.SetAll(-1);
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    mask2OpCoord[j] = Jmask;
                    int K = this.m_Basis.GetLength(j);
                    Jmask += K;
                    for (int k = 0; k < K; k++) {
                        if (k < DGField.Basis.GetLength(j)) {
                            m_Coordinates[Cell2Coord[j] + k] = DGField.Coordinates[j, k];
                        }
                    }
                }
            }

            // construction of G (constraints matrix for the optimization problem) only on mask cells
            // ======================================================================================

            var Trafo = m_grd.ChefBasis.Scaling;


            SubGrid maskSG = new SubGrid(mask);
            EdgeMask innerEM = maskSG.InnerEdgesMask;

            // get interpolation points for the continuity constraints
            NodeSet qNodes;
            int degree = m_Basis.Degree;
            int VCond;
            switch (m_grd.SpatialDimension) {
                case 1: {
                        throw new NotImplementedException("TODO");
                    }
                case 2: {
                        QuadRule quad = m_grd.Edges.EdgeRefElements[0].GetQuadratureRule(degree * 2);
                        qNodes = quad.Nodes;
                        VCond = 3;
                        break;
                    }
                case 3: {
                        QuadRule quad = m_grd.Edges.EdgeRefElements[0].FaceRefElement.GetQuadratureRule(degree * 2);
                        qNodes = quad.Nodes;
                        int Nnds = (degree + 1) * (degree + 2) / 2;
                        MultidimensionalArray nds = MultidimensionalArray.Create(Nnds, 2);
                        int node = 0;
                        for (int n1 = 0; n1 < quad.NoOfNodes; n1++) {
                            for (int n2 = 0; n2 <= n1; n2++) {
                                nds[node, 0] = qNodes[n1, 0];
                                nds[node, 1] = qNodes[n2, 0];
                                node++;
                            }
                        }
                        qNodes = new NodeSet(m_grd.Edges.EdgeRefElements[0], nds);
                        VCond = 9;
                        break;
                    }
                default:
                    throw new NotSupportedException("spatial dimension not supported");
            }


            //QuadRule quad = m_grd.Edges.EdgeRefElements[0].GetQuadratureRule(m_Basis.Degree * 2);
            //NodeSet qNodes = quad.Nodes;

            //NodeSet qNodes = m_grd.Edges.EdgeRefElements[0].GetInterpolationNodes(m_grd.Edges.EdgeRefElements[0].SupportedCellTypes.ElementAt(m_Basis.Degree));

            //int NofNds = m_Basis.Degree + 1;
            //MultidimensionalArray edgeNds = MultidimensionalArray.Create(new int[] { NofNds, 1});
            //for (int q = 0; q < NofNds; q++) {
            //    edgeNds[q, 0] = -1.0 + (2.0 / (double)m_Basis.Degree)*q;
            //}
            //NodeSet qNodes = new NodeSet(m_grd.Edges.EdgeRefElements[0], edgeNds);


            // construction of A (constraints matrix for the optimization problem) only on mask cells
            MultidimensionalArray A = MultidimensionalArray.Create(innerEM.NoOfItemsLocally * qNodes.NoOfNodes, Jmask);
            //int edgeCount = 0;
            int nodeCount = 0;

            int NoVert = m_grd.Vertices.Count;
            int[] CondAtVertice = new int[NoVert];
            int[,] IncidenceMatrix = new int[NoVert, NoVert];

            BitArray CellsAtNonConform = new BitArray(J);
            BitArray CondAtNonConform = new BitArray(J);
            List<int> NonConformCells = new List<int>();
            bool nonConformEdge;
            foreach (var chunk in innerEM) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {

                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];


                    // check for cells which are located at a non-conformal cell (hanging nodes)
                    nonConformEdge = false;
                    var edgeInfo = m_grd.Edges.Info[j];
                    if (edgeInfo.HasFlag(EdgeInfo.Cell1_Nonconformal)) {                      
                        CellsAtNonConform[cell2] = true;
                        if (NonConformCells.Contains(cell1)) {
                            nonConformEdge = true;
                        } else {
                            NonConformCells.Add(cell1);
                            //CondAtNonConform[cell2] = true;
                        }
                        //if (CellsAtNonConform[cell2]) {
                        //    //nonConformEdge = true;
                        //} else {
                        //    CellsAtNonConform[cell2] = true;
                        //}
                    } else if (edgeInfo.HasFlag(EdgeInfo.Cell2_Nonconformal)) {
                        CellsAtNonConform[cell1] = true;
                        if (NonConformCells.Contains(cell2)) {
                            nonConformEdge = true;
                        } else {
                            NonConformCells.Add(cell2);
                            //CondAtNonConform[cell1] = true;
                        }
                        //if (CellsAtNonConform[cell1]) {
                        //    //nonConformEdge = true;
                        //} else {
                        //    CellsAtNonConform[cell1] = true;
                        //}
                    };


                    // check for overdetermined vertices
                    int[] vertAtCell1 = m_grd.Cells.CellVertices[cell1];
                    int[] vertAtCell2 = m_grd.Cells.CellVertices[cell2];
                    int numVCond = 0;
                    List<int> VertAtEdge = new List<int>();
                    for (int i = 0; i < vertAtCell1.Length; i++) {
                        int vert = vertAtCell1[i];
                        if (vertAtCell2.Contains(vert)) {
                            VertAtEdge.Add(vert);
                            CondAtVertice[vert] += 1;
                            if (CondAtVertice[vert] > VCond) {
                                numVCond++;
                            }
                        }
                    }
                    // check for overdetermined edges (additional for 3D)
                    int numECond = 0;
                    int VertC = VertAtEdge.Count;
                    if (VertC > 2) {
                        int i0 = 0;
                        for (int indV1 = 0; indV1 < VertAtEdge.Count; indV1++) {
                            i0++;
                            for (int indV2 = i0; indV2 < VertAtEdge.Count; indV2++) {
                                int V1 = VertAtEdge.ElementAt(indV1);
                                int V2 = VertAtEdge.ElementAt(indV2);
                                IncidenceMatrix[V1, V2] += 1;
                                IncidenceMatrix[V2, V1] += 1;
                                if (IncidenceMatrix[V1, V2] == 4) {
                                    numECond++;
                                }

                            }
                        }
                    }

                    //QuadRule quad = m_grd.Edges.EdgeRefElements[0].GetQuadratureRule((m_Basis.Degree - numCond) * 2);
                    //qNodes = quad.Nodes;

                    // construct points
                    switch (m_grd.SpatialDimension) {
                        case 1: {
                                throw new NotImplementedException("TODO");
                            }
                        case 2: {
                                QuadRule quad = m_grd.Edges.EdgeRefElements[0].GetQuadratureRule((degree - numVCond) * 2);
                                qNodes = quad.Nodes;
                                break;
                            }
                        case 3: {
                                QuadRule quad = m_grd.Edges.EdgeRefElements[0].FaceRefElement.GetQuadratureRule((degree - numECond) * 2);
                                qNodes = quad.Nodes;
                                int Nnds = (((degree - numECond)+ 1) * ((degree - numECond) + 2) / 2);
                                MultidimensionalArray nds = MultidimensionalArray.Create(Nnds, 2);
                                int node = 0;
                                for (int n1 = 0; n1 < quad.NoOfNodes; n1++) {
                                    for (int n2 = 0; n2 <= n1; n2++) {
                                        nds[node, 0] = qNodes[n1, 0];
                                        nds[node, 1] = qNodes[n2, 0];
                                        node++;
                                    }
                                }
                                qNodes = new NodeSet(m_grd.Edges.EdgeRefElements[0], nds); //.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Nnds - (1 + numVCond), 1 }));
                                break;
                            }
                        default:
                            throw new NotSupportedException("spatial dimension not supported");
                    }


 


                    if (!nonConformEdge) {

                        int trf1 = m_grd.Edges.Edge2CellTrafoIndex[j, 0];
                        int trf2 = m_grd.Edges.Edge2CellTrafoIndex[j, 1];

                        NodeSet VolNds1 = qNodes.GetVolumeNodeSet(m_grd, trf1);
                        NodeSet VolNds2 = qNodes.GetVolumeNodeSet(m_grd, trf2);

                        //MultidimensionalArray VolNds1_global = NodeSet.Create(VolNds1.Lengths);
                        //MultidimensionalArray VolNds2_global = NodeSet.Create(VolNds2.Lengths);

                        //m_grd.TransformLocal2Global(VolNds1, VolNds1_global, cell1);
                        //m_grd.TransformLocal2Global(VolNds2, VolNds2_global, cell2);

                        //VolNds1 = new NodeSet(m_grd.Grid.GetRefElement(0), VolNds1_global);
                        //VolNds2 = new NodeSet(m_grd.Grid.GetRefElement(0), VolNds2_global);

                       
                        MultidimensionalArray NdsValues1 = m_Basis.Evaluate(VolNds1);
                        MultidimensionalArray NdsValues2 = m_Basis.Evaluate(VolNds2);

                        //NdsValues1.Scale(cdgTrafo[cell1]);
                        //NdsValues2.Scale(cdgTrafo[cell2]);

                        var results = m_Basis.EdgeEval(qNodes, j, 1);


                        for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                            // Cell1
                            for (int p = 0; p < this.m_Basis.GetLength(cell1); p++) {
                                //G[edgeCount * qNodes.NoOfNodes + qN, mask2OpCoord[cell1] + p] = NdsValues1[qN, p] * Trafo[cell1];
                                A[nodeCount + qN, mask2OpCoord[cell1] + p] = results.Item1[0, qN, p];
                            }   
                            // Cell2
                            for (int p = 0; p < this.m_Basis.GetLength(cell2); p++) {
                                //G[edgeCount * qNodes.NoOfNodes + qN, mask2OpCoord[cell2] + p] = -NdsValues2[qN, p] * Trafo[cell2];
                                A[nodeCount + qN, mask2OpCoord[cell2] + p] = -results.Item2[0, qN, p];
                            }
                        }
                        //edgeCount++;
                        nodeCount += qNodes.NoOfNodes;
                    }
                }
            }

            // resize G incase of hanging nodes
            MultidimensionalArray A_temp = A.CloneAs();
            A = A_temp.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { (nodeCount) - 1, Jmask - 1 });


            // get edges at hanging nodes
            CellMask CaNCmsk = new CellMask(m_grd, CellsAtNonConform);
            SubGrid CaNCsgrd = new SubGrid(CaNCmsk);
            EdgeMask hangingEdges = CaNCsgrd.InnerEdgesMask;

            // find corresponding hanging nodes
            int[] hangingVertices = new int[hangingEdges.NoOfItemsLocally];
            int Vind = 0;
            foreach (var chunk in hangingEdges) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];
                    int[] vertAtCell1 = m_grd.Cells.CellVertices[cell1];
                    int[] vertAtCell2 = m_grd.Cells.CellVertices[cell2];
                    for (int i = 0; i < vertAtCell1.Length; i++) {
                        if (vertAtCell2.Contains(vertAtCell1[i]) && m_grd.Vertices.VerticeToCell[vertAtCell1[i]].Length == 2) {
                            hangingVertices[Vind] = vertAtCell1[i];
                            Vind++;
                            break;
                        }
                    }
                }
            }

            // compute derivatives of the basis polynomials along the edges
            PolynomialList[,] edgeDeriv = new PolynomialList[m_Basis.Degree, m_grd.SpatialDimension];
            for (int drv = 0; drv < m_Basis.Degree; drv++) {
                for (int d = 0; d < m_grd.SpatialDimension; d++) {
                    List<Polynomial> polyList = new List<Polynomial>();
                    int[] deriv = new int[m_grd.SpatialDimension];
                    deriv[d] = drv + 1;
                    foreach (Polynomial poly in m_Basis.Polynomials[0]) {
                        Polynomial polyDeriv = poly.Derive(deriv);
                        polyList.Add(polyDeriv);
                    }
                    edgeDeriv[drv, d] = new PolynomialList(polyList);
                }
            }

            // setup additional constraints for hanging nodes
            MultidimensionalArray A2 = MultidimensionalArray.Create(hangingEdges.NoOfItemsLocally * m_Basis.Degree, Jmask);
            Vind = 0;
            foreach (var chunk in hangingEdges) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int cell1 = m_grd.Edges.CellIndices[j, 0];
                    int cell2 = m_grd.Edges.CellIndices[j, 1];

                    int trf1 = m_grd.Edges.Edge2CellTrafoIndex[j, 0];
                    int trf2 = m_grd.Edges.Edge2CellTrafoIndex[j, 1];

                    MultidimensionalArray hNd_global = m_grd.Vertices.Coordinates.ExtractSubArrayShallow(new int[] { hangingVertices[Vind], 0 },
                                                                                new int[] { hangingVertices[Vind], m_grd.SpatialDimension - 1 });

                    MultidimensionalArray hND_local1 = MultidimensionalArray.Create(1, 1, m_grd.SpatialDimension);
                    MultidimensionalArray hND_local2 = MultidimensionalArray.Create(1, 1, m_grd.SpatialDimension);

                    m_grd.TransformGlobal2Local(hNd_global, hND_local1, cell1, 1, 0);
                    m_grd.TransformGlobal2Local(hNd_global, hND_local2, cell2, 1, 0);

                    NodeSet hangNode1 = new NodeSet(m_grd.Grid.GetRefElement(0), hND_local1.ExtractSubArrayShallow(0, -1, -1));
                    NodeSet hangNode2 = new NodeSet(m_grd.Grid.GetRefElement(0), hND_local2.ExtractSubArrayShallow(0, -1, -1));

                    //NodeSet VolNds1 = qNodes.GetVolumeNodeSet(m_grd, trf1);
                    //NodeSet VolNds2 = qNodes.GetVolumeNodeSet(m_grd, trf2);

                    int derivInd = -1;
                    for (int d = 0; d < m_grd.SpatialDimension; d++) {
                        if (m_grd.Edges.NormalsForAffine[j, d] != 0)
                            derivInd = d;
                    }

                    //for (int qN = 0; qN < qNodes.NoOfNodes; qN++) {
                        for (int drv = 0; drv < m_Basis.Degree; drv++) {
                            MultidimensionalArray Res = MultidimensionalArray.Create(1, m_Basis.Polynomials[0].Count);
                            edgeDeriv[drv, derivInd].Evaluate(hangNode1, Res);
                            for (int p = 0; p < m_Basis.Polynomials[0].Count; p++) {
                                A2[Vind * m_Basis.Degree + drv, mask2OpCoord[cell1] + p] = Res[0, p] * Trafo[cell1];
                            }
                            Res.Clear();
                            edgeDeriv[drv, derivInd].Evaluate(hangNode2, Res);
                            for (int p = 0; p < m_Basis.Polynomials[0].Count; p++) {
                                A2[Vind * m_Basis.Degree + drv, mask2OpCoord[cell2] + p] = -Res[0, p] * Trafo[cell2];
                            }
                        }
                    //}
                    Vind++;
                }
            }


            // construct optimization matrix
            MsrMatrix OpMatrix = new MsrMatrix(Jmask + A.Lengths[0] + A2.Lengths[0]);
            // mass matrix
            int ind = 0;
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int K = this.m_Basis.GetLength(j);
                    for (int k = 0; k < K; k++) {
                        OpMatrix[ind, ind] = 1;
                        ind++;
                    }
                }
            }
            // A and A_Transpose
            MsrMatrix B = new MsrMatrix(A.Lengths[0], A.Lengths[1], 1, 1);
            for (int i = Jmask; i < Jmask + A.Lengths[0]; i++) {
                for (int j = 0; j < Jmask; j++) {
                    B[i - Jmask, j] = A[i - Jmask, j];
                    OpMatrix[i, j] = A[i - Jmask, j];
                    OpMatrix[j, i] = A[i - Jmask, j];
                }
            }
            // A2 and A2_Transpose
            for (int i = Jmask + A.Lengths[0]; i < Jmask + A.Lengths[0] + A2.Lengths[0]; i++) {
                for (int j = 0; j < Jmask; j++) {
                    OpMatrix[i, j] = A2[i - (Jmask + A.Lengths[0]), j];
                    OpMatrix[j, i] = A2[i - (Jmask + A.Lengths[0]), j];
                }
            }


            // construct RHS
            double[] RHS = new double[Jmask]; // + A.Lengths[0] + A2.Lengths[0]];
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int K = this.m_Basis.GetLength(j);
                    for (int k = 0; k < K; k++) {
                        RHS[mask2OpCoord[j] + k] = this.m_Coordinates[Cell2Coord[j] + k];
                    }
                }
            }

            // check
            //double[] res = new double[Jmask + G.Lengths[0]];
            //OpMatrix.SpMV(1, RHS, 0, res);

            //OpMatrix.SaveToTextFile("C:\\tmp\\OpMatrix.txt");
            //RHS.SaveToTextFile("C:\\tmp\\RHS.txt");

            // solve linear system
            //OpSolver.DefineMatrix(OpMatrix);

            //double[] x = new double[Jmask + G.Lengths[0] + G2.Lengths[0]];
            //var OpSol = OpSolver.Solve(x, RHS);


            // test with matlab
            //MultidimensionalArray xWrapper = MultidimensionalArray.Create(x.Length, 1);
            MultidimensionalArray output = MultidimensionalArray.Create(1, 3);
            MultidimensionalArray Q2 = MultidimensionalArray.Create(A.Lengths[1], A.Lengths[1] - A.Lengths[0]);
            MultidimensionalArray H = MultidimensionalArray.Create(Q2.Lengths[1], Q2.Lengths[1]);
            Console.WriteLine("Calling MATLAB/Octave...");
            using (BatchmodeConnector bmc = new BatchmodeConnector()) {
                //bmc.PutMatrix(OpMatrix.ToFullMatrixOnProc0(), "OpMatrix");
                //bmc.PutVector<double[]>(RHS, "RHS");
                //bmc.Cmd("x = linsolve(OpMatrix, RHS);");
                //bmc.GetMatrix(xWrapper, "x");

                bmc.PutSparseMatrix(A.ToMsrMatrix(), "A");
                bmc.Cmd("[m,n] = size(A)");
                bmc.Cmd("rank_A = rank(full(A))");
                bmc.Cmd("rank_AT = rank(full(A'))");
                //bmc.GetMatrix(output, "[rank_A, rank_AT]");

                bmc.Cmd("[Q,R] = qr(full(A)')");
                bmc.Cmd("[k,l] = size(Q)");
                bmc.Cmd("Q2 = Q(:,m+1:end)");
                bmc.GetMatrix(Q2, "Q2");

                bmc.Cmd("H = Q2' * Q2");
                bmc.Cmd("[V,r] = chol(H)");
                bmc.GetMatrix(H, "H");

                bmc.GetMatrix(output, "[rank_A, rank_AT, r]");

                bmc.Execute(false);

                //x.Clear();
                //xWrapper.GetColumn(0, x);
            }

            Console.WriteLine("A: No of Rows = {0}; rank = {1}", A.NoOfRows, output[0, 0]);
            Console.WriteLine("AT: No of Rows = {0}; rank = {1}", A.Transpose().NoOfRows, output[0, 1]);
            Console.WriteLine("ZT * Z positive definite = {0}", output[0, 2]);


            // null space

            //MultidimensionalArray Q;
            //MultidimensionalArray R;
            //A.Transpose().GetQRFactorization(out Q, out R);

            //int n = Q.Lengths[0];
            //int m = R.Lengths[0];
            //MsrMatrix Z = Q.ExtractSubArrayShallow(new int[] { 0, m }, new int[] { n, n }).ToMsrMatrix();
            //MsrMatrix ZT = Z.Transpose();

            //MsrMatrix ZTZ = ZT * Z;


            //MsrMatrix Z = Q2.ToMsrMatrix();
            //MsrMatrix ZT = Z.Transpose();

            ////MsrMatrix ZTZ = H.ToMsrMatrix();

            ////var solver = new ilPSP.LinSolvers.monkey.CG();
            ////solver.MatrixType = ilPSP.LinSolvers.monkey.MatrixType.Auto;
            ////solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.CPU;
            ////solver.Tolerance = 1.0e-12;

            ////double[] ZTf = new double[ZT.NoOfCols];
            //double[] v = new double[ZT.NoOfRows];
            //double[] x = new double[Z.NoOfRows];

            ////solver.DefineMatrix(ZTZ);
            //ZT.SpMV(1.0, RHS, 0.0, v);
            ////solver.Solve(v, ZTf);

            //Z.SpMV(1.0, v, 0.0, x);


            // Schur
            MsrMatrix AAT = B * B.Transpose();

            RHS = new double[B.NoOfRows];
            B.SpMV(1.0, m_Coordinates.To1DArray(), 0.0, RHS);

            double[] v = new double[B.NoOfRows];
            double[] x = new double[m_Coordinates.Length];

            var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();

            solver.DefineMatrix(AAT);
            solver.Solve(v, RHS);

            x.AccV(1.0, m_Coordinates.To1DArray());
            B.Transpose().SpMV(-1.0, v, 1.0, x);


            // reconstruct coordinate vector
            foreach (var chunk in mask) {
                int j0 = chunk.i0;
                int jE = chunk.JE;
                for (int j = j0; j < jE; j++) {
                    int K = this.m_Basis.GetLength(j);
                    for (int k = 0; k < K; k++) {
                        this.m_Coordinates[Cell2Coord[j] + k] = x[mask2OpCoord[j] + k];
                    }
                }
            }

        }


        /// <summary>
        /// accumulate this field to a DG Field
        /// </summary>
        /// <param name="alpha"></param>
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

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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.LevelSetTools.FastMarcher;
using ilPSP.Utils;

namespace BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher {

    class LocalMarcher_Structured : ILocalSolver {

        Basis solutionBasis;
        GridData gridDat;
        Dictionary<Position, Node> edgeNodes;
        List<Node> nodeGrid;
        MultidimensionalArray localCoordinates;
        int resolution;
        int iKref;
        QuadRule[] daRuleS;

        /// <summary>
        /// 
        /// ToDo: Validate for a different basis than Levelset basis.
        /// </summary>
        /// <param name="SolutionBasis"></param>
        public LocalMarcher_Structured(Basis SolutionBasis) {
            solutionBasis = SolutionBasis;
            resolution = SolutionBasis.Degree + 3;
            gridDat = (GridData)SolutionBasis.GridDat;
            daRuleS = gridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(2 * (resolution - 3))).ToArray();

            NodeSet QuadNodes = daRuleS[0].Nodes;
            double[] coordinates1D = new double[resolution];
            coordinates1D[0] = -1;
            coordinates1D[resolution - 1] = 1;
            for (int i = 0; i < resolution - 2; ++i) {
                coordinates1D[i + 1] = QuadNodes[i, 0];
            }

            localCoordinates = BuildCoordinatesAndQuadRule(coordinates1D);

            //Build Graph for Fast Marching
            SetupNodeGraph(coordinates1D);
            SetupEdgeNodes();
        }

        MultidimensionalArray BuildCoordinatesAndQuadRule(double[] coordinates1D) {

            int D = gridDat.SpatialDimension;
            int NoOfNodes = (int)Math.Pow(coordinates1D.Length, D - 1);
            MultidimensionalArray coordinates = MultidimensionalArray.Create(NoOfNodes, D - 1);
            int ndIdx = 0;
            if (D == 2) {
                for (int nd0 = 0; nd0 < coordinates1D.Length; nd0++) {
                    coordinates[ndIdx, 0] = coordinates1D[nd0];
                    ndIdx++;
                }
            } else if (D == 3) {
                for (int nd0 = 0; nd0 < coordinates1D.Length; nd0++) {
                    for (int nd1 = 0; nd1 < coordinates1D.Length; nd1++) {
                        coordinates[ndIdx, 0] = coordinates1D[nd0];
                        coordinates[ndIdx, 1] = coordinates1D[nd1];
                        ndIdx++;
                    }
                }
            } else
                throw new ArgumentOutOfRangeException("not supported spaial dimension");


            return coordinates;
        }

        public void LocalSolve(int jCell, BitArray AcceptedMask, SinglePhaseField Phi) {
            //Reset
            ResetNodeGraphValues();

            //Build Fast Marcher
            IFastMarchingQueue<IMarchingNode> Heap = new MarchingHeap((int)Math.Pow(resolution, gridDat.SpatialDimension));
            Fastmarcher marcher = new Fastmarcher(Heap);

            //Fill NodeGraph with information, i.e. pos coordinates
            FillInNodeGraph(jCell);

            //Find Accepted Nodes on Cell Edges and Set them
            Node[] Accepted = SetupAccepted(jCell, AcceptedMask, Phi);

            //Solve
            //Console.WriteLine("LocalMarcher.LocalSolver: macher.march()");
            marcher.march(Accepted);

            //Project solution to DGSpace
            MapToDGSpace(jCell, Phi);
        }

        void FillInNodeGraph(int jCell) {
            iKref = gridDat.Cells.GetRefElementIndex(jCell);
            int D = gridDat.SpatialDimension;
            MultidimensionalArray globalCoordinates = MultidimensionalArray.Create(1, D);
            MultidimensionalArray localCoordinates = MultidimensionalArray.Create(1, D);
            foreach (Node node in nodeGrid) {
                localCoordinates.ExtractSubArrayShallow(0,-1).SetVector(node.Pos_local);
                gridDat.TransformLocal2Global(localCoordinates, globalCoordinates, jCell);
                node.SetGlobalPosition(globalCoordinates.ExtractSubArrayShallow(0,-1).To1DArray());
            }
        }

        Node[] SetupAccepted(int jCell, BitArray AcceptedMask, SinglePhaseField Phi) {

            List<Node> Accepted = new List<Node>();

            //Find accepted edges
            var neighbors = gridDat.GetCellNeighboursViaEdges(jCell);

            //Set values for nodes on edges in NodeGraph
            //-------------------------------------------------------------------------------------------------------------------------------
            int D = gridDat.SpatialDimension;
            int[,] TrafoIdx = gridDat.Edges.Edge2CellTrafoIndex;
            MultidimensionalArray PhiEdge = MultidimensionalArray.Create(1, (int)Math.Pow(resolution, D - 1));

            NodeSet EdgeNodes = new NodeSet(gridDat.Edges.EdgeRefElements[0], (int)Math.Pow(resolution, D - 1), D - 1, false);
            EdgeNodes.Set(localCoordinates);
            EdgeNodes.LockForever();

            MultidimensionalArray EdgeNodesGlobal = MultidimensionalArray.Create((int)Math.Pow(resolution, D - 1), D);
            MultidimensionalArray EdgeNodesLocal = MultidimensionalArray.Create((int)Math.Pow(resolution, D - 1), D);

            //For each accepted edge set value for nodes in NodeGraph
            foreach (var neighbor in neighbors) {
                if (AcceptedMask[neighbor.Item1]) {

                    //Build Nodeset for evaluation
                    int iTrafo = TrafoIdx[neighbor.Item2, neighbor.Item3];
                    NodeSet CellNodes = EdgeNodes.GetVolumeNodeSet(gridDat, iTrafo, false);

                    //Evaluate and find global position
                    Phi.Evaluate(neighbor.Item1, 1, CellNodes, PhiEdge);
                    gridDat.TransformLocal2Global(EdgeNodes.GetVolumeNodeSet(this.gridDat, iTrafo, false), EdgeNodesGlobal, neighbor.Item1);
                    gridDat.TransformGlobal2Local(EdgeNodesGlobal, EdgeNodesLocal, jCell, null);

                    //Find corresponding nodes in NodeGraph and Set them. 
                    Node[] acceptedNodes = FindAndSetCorrespondingNodes(EdgeNodesLocal, PhiEdge);

                    //Add to Accepted  
                    Accepted.AddRange(acceptedNodes);
                }
            }
            return Accepted.ToArray();
        }

        Node[] FindAndSetCorrespondingNodes(MultidimensionalArray EdgeNode, MultidimensionalArray EdgeValue) {
            //To find corresponding nodes, the method looks up the edgenodes in the dictionary EdgeNodes using the position as the key.
            Node[] nodes = new Node[EdgeNode.GetLength(0)];
            for (int i = 0; i < EdgeNode.GetLength(0); ++i) {
                Position edgeNodePosition = new Position(EdgeNode.ExtractSubArrayShallow(i,-1).To1DArray());
                //Find correspoding Node in Edgenodes and write into nodes,
                if (!this.edgeNodes.TryGetValue(edgeNodePosition, out nodes[i])) {
                    throw new DataMisalignedException("Cannot align Edge with NodeGraph");
                }
                nodes[i].Phi = EdgeValue[0, i];
            }
            return nodes;
        }

        void SetupNodeGraph(double[] coordinates1D) {

            //Build matrix of nodes and initialize each node
            nodeGrid = new List<Node>();
            for (int nd = 0; nd < localCoordinates.Lengths[0]; ++nd) {
                for (int nd1D = 0; nd1D < coordinates1D.Length; nd1D++) {
                    double[] pos = ArrayTools.Cat(localCoordinates.ExtractSubArrayShallow(nd, -1).To1DArray(), coordinates1D[nd1D]);
                    nodeGrid.Add(new Node(pos, double.MaxValue));
                }
            }

            //Link nodes
            if (gridDat.SpatialDimension == 2)
                SetupNeighbors2D();
            else if (gridDat.SpatialDimension == 3)
                SetupNeighbors3D();
            else
                throw new ArgumentOutOfRangeException("not supported spatial dimension");

        }

        void SetupNeighbors2D() {

            for (int nd = 0; nd < nodeGrid.Count(); nd++) {
                Node node = nodeGrid.ElementAt(nd);
                //int[] potentialNeighbors = new int[] { nd - 1, nd + 1, nd - resolution, nd + resolution };

                int neighbor = nd - 1;
                if (neighbor > -1 && neighbor % resolution < resolution - 1)
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));

                neighbor = nd + 1;
                if (neighbor < nodeGrid.Count && neighbor % resolution > 0)
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));

                neighbor = nd - resolution;
                if (neighbor > -1)
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));

                neighbor = nd + resolution;
                if (neighbor < nodeGrid.Count)
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));
            }
        }

        void SetupNeighbors3D() {

            for (int nd = 0; nd < nodeGrid.Count(); nd++) {
                Node node = nodeGrid.ElementAt(nd);
                //int[] potentialNeighbors = new int[] { nd - 1, nd + 1, nd - resolution, nd + resolution, nd - resolution * resolution, nd + resolution * resolution };

                int neighbor = nd - 1;
                if (neighbor > -1 && neighbor % resolution < resolution - 1)
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));

                neighbor = nd + 1;
                if (neighbor < nodeGrid.Count() && neighbor % resolution > 0)
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));

                neighbor = nd - resolution;
                if (neighbor > -1 && neighbor % (resolution * resolution) < resolution * (resolution - 1))
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));

                neighbor = nd + resolution;
                if (neighbor < nodeGrid.Count() && neighbor % (resolution * resolution) > resolution - 1)
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));

                neighbor = nd - resolution * resolution;
                if (neighbor > -1)
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));

                neighbor = nd + resolution * resolution;
                if (neighbor < nodeGrid.Count())
                    node.neighbors.AddLast(nodeGrid.ElementAt(neighbor));
            }
        }


        void ResetNodeGraphValues() {
            foreach (Node node in nodeGrid) {
                node.Phi = double.MaxValue;
            }
        }

        void SetupEdgeNodes() {
            edgeNodes = new Dictionary<Position, Node>(new PositionComparer());

            BitArray edgeNodesBit;
            if (gridDat.SpatialDimension == 2)
                edgeNodesBit = EdgeNodes2D();
            else if (gridDat.SpatialDimension == 3)
                edgeNodesBit = EdgeNodes3D();
            else
                throw new ArgumentOutOfRangeException("not supported spatial dimension");

            for (int nd = 0; nd < edgeNodesBit.Length; nd++) {
                if (edgeNodesBit[nd]) {
                    Node EdgeNode = nodeGrid.ElementAt(nd);
                    Position positionOfEdgeNode = new Position(EdgeNode);
                    edgeNodes.Add(positionOfEdgeNode, EdgeNode);
                }
            }

        }


        BitArray EdgeNodes2D() {

            BitArray edgeNodesBit = new BitArray(nodeGrid.Count());
            for (int nd = 0; nd < edgeNodesBit.Length; nd++) {
                if (nd < resolution)
                    edgeNodesBit[nd] = true;
                if (nd % resolution == 0)
                    edgeNodesBit[nd] = true;
                if (nd % resolution == resolution - 1)
                    edgeNodesBit[nd] = true;
                if (nd >= edgeNodesBit.Length - resolution)
                    edgeNodesBit[nd] = true;
            }

            return edgeNodesBit;
        }

        BitArray EdgeNodes3D() {

            BitArray edgeNodesBit = new BitArray(nodeGrid.Count());
            for (int nd = 0; nd < edgeNodesBit.Length; nd++) {
                if (nd < resolution * resolution)
                    edgeNodesBit[nd] = true;
                if (nd % resolution == 0)
                    edgeNodesBit[nd] = true;
                if (nd % resolution == resolution - 1)
                    edgeNodesBit[nd] = true;
                if (nd >= edgeNodesBit.Length - resolution * resolution)
                    edgeNodesBit[nd] = true;
                if (nd % (resolution * resolution) > resolution * (resolution - 1))
                    edgeNodesBit[nd] = true;
                if (nd % (resolution * resolution) < (resolution - 1))
                    edgeNodesBit[nd] = true;
            }

            return edgeNodesBit;
        }


        MultidimensionalArray ConvertToQuadNodes(List<Node> Grid) {
            //int length = Grid.Count - 2;
            int ConvertedArray_Length = Grid.Count - edgeNodes.Count;
            MultidimensionalArray ConvertedArray = MultidimensionalArray.Create(ConvertedArray_Length);

            int D = gridDat.SpatialDimension;
            BitArray edgeNodesBit;
            if (D == 2)
                edgeNodesBit = EdgeNodes2D();
            else if (D == 3)
                edgeNodesBit = EdgeNodes3D();
            else
                throw new ArgumentOutOfRangeException("not supported spatial dimension");


            //Remove Edges and write to Multidimensionalarray. Make sure that everything is kept in the correct order!
            if (D == 2) {
                int i_ConvertedArray = 0;
                int i_increment = resolution - 2;
                int iShift_ConvertedArray = 1;
                for (int i = 0; i < Grid.Count; i++) {
                    if (!edgeNodesBit[i]) {
                        ConvertedArray[i_ConvertedArray] = Grid[i].Phi;
                        i_ConvertedArray += i_increment;
                        if (i_ConvertedArray >= ConvertedArray_Length) {
                            i_ConvertedArray = iShift_ConvertedArray;
                            iShift_ConvertedArray++;
                        }
                    }
                }
            }

            if (D == 3) {
                int i_ConvertedArray = 0;
                int i_increment = (resolution - 2) * (resolution - 2);
                int iShift_ConvertedArray = (resolution - 2);
                int iShift_increment = (resolution - 2);
                int iShift_ConvertedArray2 = 1;
                for (int i = 0; i < Grid.Count; i++) {
                    if (!edgeNodesBit[i]) {
                        //Console.WriteLine("i_ConvertedArray = {0}", i_ConvertedArray);
                        ConvertedArray[i_ConvertedArray] = Grid[i].Phi;
                        i_ConvertedArray += i_increment;
                        if (i_ConvertedArray >= ConvertedArray_Length) {
                            i_ConvertedArray = iShift_ConvertedArray;
                            iShift_ConvertedArray += iShift_increment;
                            if (iShift_ConvertedArray >= i_increment) {
                                iShift_ConvertedArray = iShift_ConvertedArray2;
                                iShift_ConvertedArray2++;
                            }
                        }
                    }
                }
            }



            return ConvertedArray;
        }

        void MapToDGSpace(int jCell, SinglePhaseField Phi) {

            MultidimensionalArray PhiAtQuadNodes = ConvertToQuadNodes(nodeGrid);
            int numberOfQuadNodes = PhiAtQuadNodes.Length;
            MultidimensionalArray weighted_PhiAtQuadNodes = MultidimensionalArray.Create(numberOfQuadNodes);

            //NodeSet quadNodes = this.daRuleS[0].Nodes;

            for (int i = 0; i < numberOfQuadNodes; i++) { // loop over all quadrature nodes
                weighted_PhiAtQuadNodes[i] = PhiAtQuadNodes[i] * this.daRuleS[iKref].Weights[i];
            }

            var BasisValues = this.solutionBasis.CellEval(this.daRuleS[iKref].Nodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);

            if (this.gridDat.Cells.IsCellAffineLinear(jCell)) {
                int N = this.solutionBasis.GetLength(jCell);
                int N2 = Phi.Basis.GetLength(jCell);

                MultidimensionalArray Phi_1 = MultidimensionalArray.Create(N);
                double scale = this.gridDat.Cells.JacobiDet[jCell];
                Phi_1.Multiply(scale, BasisValues, weighted_PhiAtQuadNodes, 0.0, "m", "km", "k");
                for (int n = 0; n < N; n++)
                    Phi.Coordinates[jCell, n] = Phi_1[n];
                for (int n = N; n < N2; n++) {
                    Phi.Coordinates[jCell, n] = 0;
                }
            } else {
                throw new NotImplementedException("not implemented for curved cells");
            }
        }
    }
}

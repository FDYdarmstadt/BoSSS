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

namespace BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher {

    class LocalMarcher_2DStructured : ILocalSolver {

        Basis solutionBasis;
        GridData gridDat;
        Dictionary<Position,Node> edgeNodes;
        Node[,] nodeGrid;
        double[] localCoordinates;
        int resolution;
        int iKref;
        QuadRule[] daRuleS;

        /// <summary>
        /// 
        /// ToDo: Validate for a different basis than Levelset basis.
        /// </summary>
        /// <param name="SolutionBasis"></param>
        public LocalMarcher_2DStructured(Basis SolutionBasis) {
            solutionBasis = SolutionBasis;
            resolution = SolutionBasis.Degree + 3; 
            gridDat = (GridData)SolutionBasis.GridDat;
            localCoordinates = BuildCoordinatesAndQuadRule(resolution);

            //Build Graph for Fast Marching
            SetupNodeGraph();
            SetupEdgeNodes();
        }

        double[] BuildCoordinatesAndQuadRule(int resolution) {
            daRuleS = gridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(2* (resolution - 3))).ToArray();
            NodeSet QuadNodes = daRuleS[0].Nodes;
            double[] coordinates = new double[resolution];
            coordinates[0] = -1;
            coordinates[resolution - 1] = 1;
            for (int i = 0; i < resolution - 2; ++i) {
                coordinates[i + 1] = QuadNodes[i,0];
            }
            return( coordinates);
        }

        public void LocalSolve(int jCell, BitArray AcceptedMask, SinglePhaseField Phi) {
            //Reset
            ResetNodeGraphValues();

            //Build Fast Marcher
            IFastMarchingQueue < IMarchingNode > Heap = new MarchingHeap(resolution * resolution);
            Fastmarcher marcher = new Fastmarcher(Heap);

            //Fill NodeGraph with information, i.e. x and y coordinates
            FillInNodeGraph(jCell);

            //Find Accepted Nodes on Cell Edges and Set them
            Node[] Accepted = SetupAccepted(jCell, AcceptedMask, Phi);
            
            //Solve
            marcher.march(Accepted);
            
            //Project solution to DGSpace
            MapToDGSpace(jCell, Phi);
        }

        void FillInNodeGraph(int jCell) {
            iKref = gridDat.Cells.GetRefElementIndex(jCell);
            MultidimensionalArray globalCoordinates = MultidimensionalArray.Create(1,2);
            MultidimensionalArray localCoordinates = MultidimensionalArray.Create(1,2);
            foreach (Node node in nodeGrid) {
                localCoordinates[0, 0] = node.X_local;
                localCoordinates[0, 1] = node.Y_local;
                gridDat.TransformLocal2Global(localCoordinates, globalCoordinates, jCell);
                node.SetGlobalPosition(globalCoordinates[0, 0], globalCoordinates[0, 1]);
            }
        }

        Node[] SetupAccepted(int jCell, BitArray AcceptedMask, SinglePhaseField Phi) {

            List<Node> Accepted = new List<Node>();
            
            //Find accepted edges
            Tuple<int, int, int>[] neighbors = gridDat.GetCellNeighboursViaEdges(jCell);

            //Set values for nodes on edges in NodeGraph
            //-------------------------------------------------------------------------------------------------------------------------------
            int[,] TrafoIdx = gridDat.Edges.Edge2CellTrafoIndex;
            MultidimensionalArray PhiEdge = MultidimensionalArray.Create(1,resolution);

            NodeSet EdgeNodes = new NodeSet(gridDat.Edges.EdgeRefElements[0], resolution, 1);
            EdgeNodes.ExtractSubArrayShallow(-1, 0).SetVector(localCoordinates);
            EdgeNodes.LockForever();

            MultidimensionalArray EdgeNodesGlobal = MultidimensionalArray.Create(resolution, gridDat.SpatialDimension);
            MultidimensionalArray EdgeNodesLocal = MultidimensionalArray.Create(resolution, gridDat.SpatialDimension );

            //For each accepted edge set value for nodes in NodeGraph
            foreach (Tuple<int, int , int> neighbor in neighbors) {
                if (AcceptedMask[neighbor.Item1]) {

                    //Build Nodeset for evaluation
                    int iTrafo = TrafoIdx[neighbor.Item2, neighbor.Item3];
                    NodeSet CellNodes = EdgeNodes.GetVolumeNodeSet(gridDat, iTrafo);

                    //Evaluate and find global position
                    Phi.Evaluate(neighbor.Item1, 1, CellNodes, PhiEdge);
                    gridDat.TransformLocal2Global(EdgeNodes.GetVolumeNodeSet(this.gridDat, iTrafo), EdgeNodesGlobal, neighbor.Item1);
                    gridDat.TransformGlobal2Local( EdgeNodesGlobal, EdgeNodesLocal, jCell, null);

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
                Position edgeNodePosition = new Position(EdgeNode[i, 0], EdgeNode[i, 1]);
                //Find correspoding Node in Edgenodes and write into nodes,
                if (!this.edgeNodes.TryGetValue(edgeNodePosition, out nodes[i])) {
                    throw new DataMisalignedException("Cannot align Edge with NodeGraph");
                }
                nodes[i].Phi = EdgeValue[0, i];
            }
            return nodes;   
        }

        void SetupNodeGraph() {
            //Build matrix of nodes and initialize each node
            nodeGrid = new Node[resolution,resolution];
            for (int i = 0; i < resolution; ++i) {
                for (int j = 0; j < resolution; ++j) {
                    nodeGrid[i, j] = new Node(localCoordinates[i], localCoordinates[j], double.MaxValue);
                }
            }
            //Link nodes
            for (int i = 0; i < resolution; ++i) {
                for (int j = 0; j < resolution; ++j) {
                    Node node = nodeGrid[i, j];
                    for (int k = -1; k < 2; k += 2) {
                        int neighbor_x = i + k;
                        int neighbor_y = j ;
                        if (neighbor_x > -1 && neighbor_x < resolution && neighbor_y > -1 && neighbor_y < resolution)
                            node.neighbors.AddLast(nodeGrid[neighbor_x, neighbor_y]);
                    }
                    for (int l = -1; l < 2; l += 2) {
                        int neighbor_x = i;
                        int neighbor_y = j + l;
                        if (neighbor_x > -1 && neighbor_x < resolution && neighbor_y > -1 && neighbor_y < resolution)
                            node.neighbors.AddLast(nodeGrid[neighbor_x, neighbor_y]);
                    }
                }
            }
        }

        void ResetNodeGraphValues() {
            foreach (Node node in nodeGrid) {
                node.Phi = double.MaxValue;
            }
        }

        void SetupEdgeNodes() {
            edgeNodes = new Dictionary<Position, Node>(new PositionComparer());
            
            //Go through all 4 edges.
            for (int j = 0; j < resolution; ++j) {
                int i_x = 0;
                int i_y = j; 
                Node EdgeNode = nodeGrid[i_x, i_y];
                Position positionOfEdgeNode = new Position(EdgeNode);
                edgeNodes.Add(positionOfEdgeNode, EdgeNode);
            }
            for (int j = 0; j < resolution; ++j) {
                int i_x = resolution - 1;
                int i_y = j;
                Node EdgeNode = nodeGrid[i_x, i_y];
                Position positionOfEdgeNode = new Position(EdgeNode);
                edgeNodes.Add(positionOfEdgeNode, EdgeNode);
            }
            for (int j = 1; j < resolution - 1; ++j) {
                int i_x = j ;
                int i_y = 0;
                Node EdgeNode = nodeGrid[i_x, i_y];
                Position positionOfEdgeNode = new Position(EdgeNode);
                edgeNodes.Add(positionOfEdgeNode, EdgeNode);
            }
            for (int j = 1; j < resolution - 1; ++j) {
                int i_x = j;
                int i_y = resolution - 1;
                Node EdgeNode = nodeGrid[i_x, i_y];
                Position positionOfEdgeNode = new Position(EdgeNode);
                edgeNodes.Add(positionOfEdgeNode, EdgeNode);
            }

        }

        MultidimensionalArray ConvertToQuadNodes(Node[,] Grid) {
            int length = Grid.GetLength(0) - 2;
            int ConvertetdArray_Length = length * length;
            MultidimensionalArray ConvertedArray = MultidimensionalArray.Create(ConvertetdArray_Length);

            //Remove Edges and write to Multidimensionalarray. Make sure that everything is kept in the correct order!
            for (int i = 0; i < length; ++i) {
                for (int j = 0; j < length; ++j) {
                    int i_ConvertedArray = i * length + j;
                    ConvertedArray[i_ConvertedArray] = Grid[j + 1, i + 1].Phi;
                }
            }
            return ConvertedArray;  
        }

        void MapToDGSpace(int jCell, SinglePhaseField Phi) {

            int numberOfQuadNodes = (resolution - 2) * (resolution - 2);
            MultidimensionalArray PhiAtQuadNodes = ConvertToQuadNodes(nodeGrid);
            MultidimensionalArray weighted_PhiAtQuadNodes = MultidimensionalArray.Create(numberOfQuadNodes);


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
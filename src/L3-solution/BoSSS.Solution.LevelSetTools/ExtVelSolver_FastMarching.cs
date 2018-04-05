using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Solution.LevelSetTools.Reinit.FastMarch {
    
    /// <summary>
    /// Finds the extension velocity with a geometric approach. Only for 2d-LevelSet.
    /// Phi is calculated at a defined set of nodes and subsequently projected onto the original Phi field.
    /// This is necessary, because Phi is usually not reperesented by a nodal basis.
    /// Input: CellBoundaries with corresponding Phi values. 
    /// Output: Updated Phi. 
    /// </summary>
    /// <param name="Phi"></param>
    /// <param name="GradPhi"></param>
    /// <param name="ExtProperty"></param>
    /// <param name="ExtPropertyMin"></param>
    /// <param name="ExtPropertyMax"></param>
    /// <param name="jCell"></param>
    /// <param name="Accepted"></param>
    /// <param name="signMod"></param>
    public class ExtVelSolver_FastMarching {

        Basis SolverBasis;
        GridData GridDat;
        Edge[] edgesOfCell;
        MultidimensionalArray[] QuadNodesGlobal;
        QuadRule[] DaRuleS;
        int numberOfUsedEdges;

        /// <summary>
        /// Initializes the Solver. Creates the objects that are used in the solver-method.
        /// Outlines the QuadNodesGlobal which are used as nodes for the geometric approach.
        /// Creates the edges on which the values will be defined in the solver-method.
        /// Creates a Cell which holds all the information of the parameters in the inside of the cell
        /// </summary>
        /// <param name="ExtPropertyBasis"></param>
        public ExtVelSolver_FastMarching(Basis ExtPropertyBasis) {
            this.SolverBasis = ExtPropertyBasis;
            this.GridDat = ExtPropertyBasis.GridDat; 

            //Build Edges.
            double[] NodeGrid = Grid1D.Linspace(-1, 1, Math.Max(20, (this.SolverBasis.Degree + 1) * 2));
            int maxNumberOfEdges = 4;
            edgesOfCell = new Edge[maxNumberOfEdges]; 
            for (int i = 0; i < maxNumberOfEdges; i++) {
                edgesOfCell[i] = new Edge(NodeGrid, GridDat); 
            }

            DaRuleS = this.GridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(this.SolverBasis.Degree * 2)).ToArray();
            QuadNodesGlobal = DaRuleS.Select(rule => MultidimensionalArray.Create(rule.NoOfNodes, GridDat.SpatialDimension)).ToArray();

        }

        /// <summary>
        /// Finds the new Extension Property and overwrites the respective field. 
        /// </summary>
        /// <param name="Phi"></param> The LevelSetFunction /field
        /// <param name="ExtProperty"></param> Extension property that will be updated
        /// <param name="Accepted_Mutuable"></param> Mapping of accepted Cells
        /// <param name="jCell"></param> Name of the cell, in which we find the extension Property
        /// <param name="signMod"></param> Information on the sign of LevelSetFunction Phi
        public void ExtVelSolve_FastMarching(SinglePhaseField Phi, ConventionalDGField ExtProperty, BitArray Accepted_Mutuable, int jCell, double signMod) {

            GridData grid = Phi.GridDat;
            numberOfUsedEdges = 0;

            //Find neighbors that define the cells edges
            //Scheme: find used edges, associate a neighbouring edge to each used edge, count the number of used edges.
            Tuple<int, int, int>[] neighbors = this.GridDat.Cells.GetCellNeighboursViaEdges(jCell);
            foreach (Tuple<int, int, int> neighbor in neighbors) {
                if (Accepted_Mutuable[neighbor.Item1] == true) {
                    edgesOfCell[numberOfUsedEdges].AssociatedNeighbour = neighbor; //NumberOfEdges is incremented in class Edge
                    numberOfUsedEdges++;
                }
            }

            //Create nodes on the accepted neighbor edges, clear the buffered values.
            int[,] TrafoIdx = this.GridDat.Edges.Edge2CellTrafoIndex;
            for (int i = 0; i < numberOfUsedEdges; i++) { // loop over accepted neighbours
                edgesOfCell[i].clearBuffer();
                edgesOfCell[i].createNodes(TrafoIdx, Phi, ExtProperty);
            }
        
            //Create cell which holds information about:
            //The quadraturenodes, i.e. value and position
            Cell cell = new Cell(jCell, QuadNodesGlobal, GridDat, DaRuleS);

            //Find new values for the nodes in the cell 
            foreach (Node quadNode in cell) {  
                int[] selectedNode = findCorrelatingNode(quadNode, edgesOfCell, signMod);
                //Get Edge value at selectedNode (node with smallest difference in Phi) 
                double valueOfSelectedNode = edgesOfCell[selectedNode[0]].getExtEdgeValue(selectedNode[1]); 
                //Write the found value for Ext into the array, which holds the value of each Node. This array will be used to to update ExtProperty
                cell.setExtAtQuadNode(valueOfSelectedNode, quadNode);
            }

            //update ExtProperty by projection
            updateExtProperty(ExtProperty, cell);
        }

        /// <summary>
        /// Holds the information on each edge of the cell. The values of Phi/Ext at the neighbouring cells are stored in the PhiEdgeEvalBuffer,ExtEdgeEvalBuffer, respectively.
        /// The numberOfUsedEdges indicates how many neighbouring edges are used to find ExtProperty. 
        /// </summary>
        class Edge { 

            static int numberOfUsedEdges;
            public MultidimensionalArray ExtEdgeEvalBuffer;
            public MultidimensionalArray PhiEdgeEvalBuffer;
            public MultidimensionalArray EdgeNodesGlobal;
            public NodeSet EdgeNodes;
            public int noOfEdgeNodes;
            GridData Solver_Grid;
            Tuple<int, int, int> associatedNeighbour; 
            
            public Edge(double[] Grid1D, GridData Solver_Grid) {
                int noOfEdgeNodes = Grid1D.Length;
                this.Solver_Grid = Solver_Grid;

                initializeNodeSet(Grid1D);
                initializeBuffer(); 
            }

            private void initializeNodeSet(double[] Grid1D) {
                noOfEdgeNodes = Grid1D.Length;
                this.EdgeNodes = new NodeSet(Solver_Grid.Edges.EdgeRefElements[0], noOfEdgeNodes, 1);
                this.EdgeNodes.ExtractSubArrayShallow(-1, 0).SetVector(Grid1D);
                this.EdgeNodes.LockForever();
            }

            private void initializeBuffer() {
                PhiEdgeEvalBuffer = MultidimensionalArray.Create(1, noOfEdgeNodes);
                ExtEdgeEvalBuffer = MultidimensionalArray.Create(1, noOfEdgeNodes);
                EdgeNodesGlobal = MultidimensionalArray.Create(noOfEdgeNodes, Solver_Grid.SpatialDimension);
            }

            public void clearBuffer() {
                this.PhiEdgeEvalBuffer.Clear();
                this.ExtEdgeEvalBuffer.Clear();
            } 

            public Tuple<int, int, int> AssociatedNeighbour {
                set {
                    this.associatedNeighbour = value;
                    Edge.numberOfUsedEdges++;
                }
                get {
                    return associatedNeighbour;
                }
            }

            public void createNodes(int[,] TrafoIdx, SinglePhaseField Phi, ConventionalDGField ExtProperty) {
                int iTrafo = TrafoIdx[associatedNeighbour.Item2, associatedNeighbour.Item3];
                NodeSet CellNodes = this.EdgeNodes.GetVolumeNodeSet(this.Solver_Grid, iTrafo);

                //Writes Phi and ExtProperty values at edge nodes into the respective Buffer
                Phi.Evaluate(associatedNeighbour.Item1, 1, CellNodes, this.PhiEdgeEvalBuffer);
                ExtProperty.Evaluate(associatedNeighbour.Item1, 1, CellNodes, this.ExtEdgeEvalBuffer);

                //Writes the corresponding nodes into CellNodesGlobalBuffer
                this.Solver_Grid.TransformLocal2Global(this.EdgeNodes.GetVolumeNodeSet(this.Solver_Grid, iTrafo), this.EdgeNodesGlobal, associatedNeighbour.Item1);
            }

            public static int NumberOfEdges {
                get {
                    return numberOfUsedEdges;
                }
            }

            public static void clearStatic() {
                Edge.numberOfUsedEdges = 0;
            }

            public double getExtEdgeValue(int numberOfQuadNode) {
                return ExtEdgeEvalBuffer[0, numberOfQuadNode]; 
            }



        }

        /// <summary>
        /// Holds the local information of the cell itself. 
        /// </summary>
        class Cell : IEnumerable {
            public int jCell;
            public int iKref; //RefNumber of each cell
            public int numberOfQuadNodes;
            MultidimensionalArray QuadNodesGlobal; 
            public MultidimensionalArray ExtAtQuadNodes;

            public Cell(int JCell, MultidimensionalArray[] AllQuadNodesGlobal, GridData Solver_Grid, QuadRule[] DaRuleS) {
                this.jCell = JCell;
                this.iKref = Solver_Grid.Cells.GetRefElementIndex(JCell);
                this.QuadNodesGlobal = AllQuadNodesGlobal[iKref];
                Solver_Grid.TransformLocal2Global(DaRuleS[iKref].Nodes, this.QuadNodesGlobal, jCell);
                numberOfQuadNodes = QuadNodesGlobal.GetLength(0);
                this.ExtAtQuadNodes = MultidimensionalArray.Create(numberOfQuadNodes);

            }

            public void setExtAtQuadNodes(double value, int i) {
                ExtAtQuadNodes[i] = value;
            }

            public IEnumerator GetEnumerator() {
                for (int i = 0; i < numberOfQuadNodes; i++) {
                    double[] node = { QuadNodesGlobal[i, 0], QuadNodesGlobal[i, 1] };
                    Node quadNode = new Node(node, i);

                    yield return quadNode;
                }
            }

            public void setExtAtQuadNode(double value, Node quadNode) {
                ExtAtQuadNodes[quadNode.I] = value;
            }

        }

        /// <summary>
        /// Holds the information of one node, i.e. the position in double[] node , the number of the node in int i
        /// </summary>
        class Node{
            int i;
            double[] node;

            public Node(double[] Node, int I) {
                this.i = I;
                this.node = Node;
            }

            public int I {
                get {
                    return i; 
                }
            }

            public double[] NodePosition {
                get {
                    return node; 
                }
            }
        }

        void FastMarch(Edge[] edgesOfCell, Cell cell, double signMod) {
            //To do
            
        }

        /// <summary>
        /// Finds the node on the edge that corresponds to the node in the cell (InputNode).
        /// </summary>
        /// <param name="InputNode"></param>Node in cell.
        /// <param name="edgesOfCell"></param> Array of surrounding edges
        /// <param name="signMod"></param>sign of Phi/LevelSet
        /// <returns></returns>
        int[] findCorrelatingNode(Node InputNode, Edge[] edgesOfCell, double signMod) {
            double[] inputNode = InputNode.NodePosition; 
            double[] EdgeNode = new double[2];
            int[] selectedNode = new int[2]; // {edge, No of node}
            double dist_min = double.MaxValue;
            double dist;

            // find closest point: brute force approach
            // loop over all edges with known values
            for (int i = 0; i < numberOfUsedEdges; i++ ) {
                // loop over all nodes on this edge
                for (int i_EdgeNode = 0; i_EdgeNode < edgesOfCell[i].noOfEdgeNodes; i_EdgeNode++) { 
                    EdgeNode[0] = edgesOfCell[i].EdgeNodesGlobal[i_EdgeNode, 0];
                    EdgeNode[1] = edgesOfCell[i].EdgeNodesGlobal[i_EdgeNode, 1];
                    double phi = edgesOfCell[i].PhiEdgeEvalBuffer[0,i_EdgeNode];

                    phi *= signMod;
                    dist = GenericBlas.L2Dist(EdgeNode, inputNode) + phi;

                    if (dist < dist_min) {
                        dist_min = dist;
                        selectedNode[0] = i;
                        selectedNode[1] = i_EdgeNode;
                    }
                }
            }
            return (selectedNode); 
        }

        /// <summary>
        /// Updates the ExtProperty field via projection. Needs more work. 
        /// </summary>
        /// <param name="ExtProperty"></param> Field that will be updated
        /// <param name="cell"></param> The cell that will be updated. Cell holds the new values etc. 
        void updateExtProperty(ConventionalDGField ExtProperty, Cell cell) {
            int jCell = cell.jCell;
            int iKref = cell.iKref;
            int numberOfQuadNodes = cell.numberOfQuadNodes;
            MultidimensionalArray ExtAtQuadNodes = cell.ExtAtQuadNodes; 

            MultidimensionalArray weighted_ExtAtQuadNodes = MultidimensionalArray.Create(numberOfQuadNodes);

            for (int i = 0; i < numberOfQuadNodes; i++) { // loop over all quadrature nodes
                weighted_ExtAtQuadNodes[i] = ExtAtQuadNodes[i] * this.DaRuleS[iKref].Weights[i];
            }

            var BasisValues = this.SolverBasis.CellEval(this.DaRuleS[iKref].Nodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);

            if (this.GridDat.Cells.IsCellAffineLinear(jCell)) {
                int N = this.SolverBasis.GetLength(jCell);
                int N2 = ExtProperty.Basis.GetLength(jCell);

                MultidimensionalArray Phi_1 = MultidimensionalArray.Create(N);
                double scale = this.GridDat.Cells.JacobiDet[jCell];
                Phi_1.Multiply(scale, BasisValues, weighted_ExtAtQuadNodes, 0.0, "m", "km", "k");
                for (int n = 0; n < N; n++)
                    ExtProperty.Coordinates[jCell, n] = Phi_1[n];
                for (int n = N; n < N2; n++) {
                    ExtProperty.Coordinates[jCell, n] = 0;
                }
            } else {
                throw new NotImplementedException("not implemented for curved cells");
            }
        }
    }
}




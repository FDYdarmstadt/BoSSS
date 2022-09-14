using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Connectors.Matlab;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading;

namespace BoSSS.Solution.AdvancedSolvers {
    public class AggregationGridCurvedInjector {
        /// <summary>
        /// The function computes a Np*Np operator for each cell of the grid.
        /// This operator transforms the basis from aggregation level 0 to the
        /// maximum aggregation level. It preserves continuity along inner aggregation cell edges
        /// and the orthonormality of the basis on the aggregation cell. 
        /// Furthermore it seeks a solution that minimizes the gradient jump along inner edges.
        /// </summary>
        /// <param name="_agGrd">Grid Data for the coarse aggregation level</param>
        /// <param name="_maxDgBasis">underlying DG basis</param>
        /// <param name="_InjectorCoarse"></param>
        /// <returns> 
        /// - index corresponds to the index of the geometric cell
        /// </returns>
        public static void AggregateCurvedCells(AggregationGridData _agGrd, Basis _maxDgBasis, MultidimensionalArray[] _InjectorCoarse) {
            AggregationGridData agGrd = _agGrd;
            Basis dgBasis = _maxDgBasis;

            // check type of underlying reference elements
            foreach(var refel in agGrd.AncestorGrid.iGeomCells.RefElements) {
                if(refel.GetType() != typeof(BoSSS.Foundation.Grid.RefElements.Square)) { throw new NotSupportedException("currently only square elements are supported"); }
            }

            // get dimensions of Injection operator
            int Jagg = agGrd.iLogicalCells.NoOfLocalUpdatedCells;
            int Np = dgBasis.Length;

            MultidimensionalArray[] InjectorCoarseLogicalCell = new MultidimensionalArray[Jagg];
            MultidimensionalArray[][] RotationCoarseLogicalCell = new MultidimensionalArray[Jagg][];

            Console.WriteLine($"Building direct Injector to coarsest level for {Jagg} cells on MPI-Rank {agGrd.MpiRank}: ");

            // get a list of all inner edges
            int[] iedges;
            List<int> l_iedges = new List<int>();

            for(int i = 0; i < Jagg; i++) {

                Console.WriteLine($"aggregating cell {i}");
                // get parts of current aggregation cell
                int[] parts = agGrd.iLogicalCells.AggregateCellToParts[i];

                l_iedges.Clear();

                for(int j = 0; j < agGrd.iGeomEdges.Count; j++) {
                    // check for outer edges and if both cells on the edge belong to the current aggregation cell     
                    if(agGrd.iGeomEdges.LogicalCellIndices[j, 0] == agGrd.iGeomEdges.LogicalCellIndices[j, 1] && agGrd.iGeomEdges.LogicalCellIndices[j, 0] == i) {
                        l_iedges.Add(j);
                    }
                }
                // convert to array
                iedges = l_iedges.ToArray();

                InjectorCoarseLogicalCell[i] = MultidimensionalArray.Create(parts.Length, Np, Np);
                RotationCoarseLogicalCell[i] = new MultidimensionalArray[parts.Length];

                // retrieve rotation matrices

                RotationCoarseLogicalCell[i] = RotateBasis(agGrd, dgBasis, parts);


                // special case for the first basis function
                InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(-1, 0, 0).Acc(1, GetCoefficients0Degree(agGrd, dgBasis, parts));

                // compute the columns in the Injection operator for subsequent basis functions
                for(int n = 1; n < Np; n++) {
                    //InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(new int[] { 0, 0, n }, new int[] { parts.Length - 1, n, -1 }).Acc(1, GetCoefficients(agGrd, dgBasis, InjectorCoarseLogicalCell[i], RotationCoarseLogicalCell[i], parts, iedges, n));
                    //InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(new int[] { 0, 0, n }, new int[] { parts.Length - 1, n, -1 }).Acc(1, GetCoefficientsIntegral(agGrd, dgBasis, InjectorCoarseLogicalCell[i], RotationCoarseLogicalCell[i], parts, iedges, n));
                    InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(new int[] { 0, 0, n }, new int[] { parts.Length - 1, n, -1 }).Acc(1, GetCoefficientsGlobalSF(agGrd, dgBasis, InjectorCoarseLogicalCell[i], RotationCoarseLogicalCell[i], parts, iedges, n));
                }

                // initial orthonormalization
                #region reorthonormalize

                MultidimensionalArray MM = MultidimensionalArray.Create(Np, Np);
                MultidimensionalArray ortho = MultidimensionalArray.Create(Np, Np);
                MultidimensionalArray orthoInv = MultidimensionalArray.Create(Np, Np);

                foreach(int k in parts) {
                    MM.GEMM(1.0, InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(Array.IndexOf(parts, k), -1, -1), InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(Array.IndexOf(parts, k), -1, -1), 1.0, true);
                }

                MM.Cholesky();
                MM.TransposeTo(ortho);
                ortho.InvertTo(orthoInv);

                foreach(int k in parts) {
                    InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(Array.IndexOf(parts, k), -1, -1).GEMM(1.0, InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(Array.IndexOf(parts, k), -1, -1).CloneAs(), orthoInv, 0.0);
                }

                #endregion

                // alternative algorithm using a hierarchical cell to cell procedure with final reorthonormalization
                //InjectorCoarseLogicalCell[i] = GetInjector(agGrd, dgBasis, RotationCoarseLogicalCell[i], parts, iedges);


                // accumulate the direct operator
                foreach(int k in parts) {
                    // apply rotation to raw injector and accumulate in direct injector array
                    _InjectorCoarse[k] = MultidimensionalArray.Create(Np, Np);
                    _InjectorCoarse[k].GEMM(1.0, RotationCoarseLogicalCell[i][Array.IndexOf(parts, k)], InjectorCoarseLogicalCell[i].ExtractSubArrayShallow(Array.IndexOf(parts, k), -1, -1), 0.0);
                }

            }

            Console.WriteLine("Done. ");

        }



        // Successive build-up of the injectors (unlikely to work :'( , somethings wrong the plot is not continuous)
        private static MultidimensionalArray GetInjector(AggregationGridData _agGrd, Basis _dgBasis, MultidimensionalArray[] _rot, int[] _parts, int[] _iedges) {
            int partCount = _parts.Length;
            int Np = _dgBasis.Polynomials[0].Count;
            var Injector = MultidimensionalArray.Create(partCount, Np, Np);

            // extract cell pairs
            int[,] pairs = ExtractCellPairs(_agGrd, _parts);

            // leave root cell unchanged
            Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, pairs[0, 0]), -1, -1).AccEye(1.0);

            // iterate over each cell pair
            for(int i = 0; i < pairs.GetLength(0); i++) {
                var thisInjector = Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, pairs[i, 1]), -1, -1);
                int cell_i = pairs[i, 0];
                int cell_j = pairs[i, 1];

                int edge = pairs[i, 2];
                int[] thisParts = new int[] { cell_i, cell_j };
                MultidimensionalArray[] thisRots = new MultidimensionalArray[] { _rot[Array.IndexOf(_parts, pairs[i, 0])], _rot[Array.IndexOf(_parts, pairs[i, 1])] };

                thisInjector[0, 0] = GetCoefficients0DegreeSuccessive(_agGrd, _dgBasis, Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, pairs[i, 0]), -1, -1), thisParts);

                // compute both variants with nonzero diagonal (+,-) and select the better
                for(int n = 1; n < Np; n++) {
                    //(var injPos, var qPos) = GetCoefficientsSuccessive(_agGrd, _dgBasis, Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, pairs[i, 0]), -1, -1), thisRots, thisParts, edge, n, true);
                    //(var injNeg, var qNeg) = GetCoefficientsSuccessive(_agGrd, _dgBasis, Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, pairs[i, 0]), -1, -1), thisRots, thisParts, edge, n, false);

                    (var injPos, var qPos) = GetCoefficientsScaleFit(_agGrd, _dgBasis, Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, pairs[i, 0]), -1, -1), thisRots, thisParts, edge, n, true);
                    (var injNeg, var qNeg) = GetCoefficientsScaleFit(_agGrd, _dgBasis, Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, pairs[i, 0]), -1, -1), thisRots, thisParts, edge, n, true);

                    var inj = qPos <= qNeg ? injPos : injNeg;

                    thisInjector.ExtractSubArrayShallow(new int[] { 0, n }, new int[] { n, -1 }).Acc(1, inj);

                }



            }

            #region reorthonormalize

            MultidimensionalArray MM = MultidimensionalArray.Create(Np, Np);
            MultidimensionalArray ortho = MultidimensionalArray.Create(Np, Np);
            MultidimensionalArray orthoInv = MultidimensionalArray.Create(Np, Np);

            foreach(int k in _parts) {
                MM.GEMM(1.0, Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, k), -1, -1), Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, k), -1, -1), 1.0, true);
            }

            MM.Cholesky();
            MM.TransposeTo(ortho);
            ortho.InvertTo(orthoInv);

            foreach(int k in _parts) {
                Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, k), -1, -1).GEMM(1.0, Injector.ExtractSubArrayShallow(Array.IndexOf(_parts, k), -1, -1).CloneAs(), orthoInv, 0.0);
            }

            #endregion

            return Injector;
        }

        // Construct rotation Matrices to orient all cells in one aggregated cell in the same way
        private static MultidimensionalArray[] RotateBasis(AggregationGridData _agGrd, Basis _maxDgBasis, int[] parts) {
            int D = _maxDgBasis.GridDat.SpatialDimension;
            MultidimensionalArray[] RotationMatrices = new MultidimensionalArray[parts.Length];
            int Np = _maxDgBasis.Length;

            // extract cell pairs
            int[,] pairs = ExtractCellPairs(_agGrd, parts);

            if(D == 3) {
                throw new NotImplementedException();
            } else if(D == 2) {



                RotationMatrices[Array.IndexOf(parts, pairs[0, 0])] = MultidimensionalArray.Create(Np, Np);
                // First cell serves as reference without rotation
                RotationMatrices[Array.IndexOf(parts, pairs[0, 0])].AccEye(1.0);

                for(int i = 0; i < pairs.GetLength(0); i++) {
                    int iCell = Array.IndexOf(parts, pairs[i, 0]);
                    int jCell = Array.IndexOf(parts, pairs[i, 1]);

                    RotationMatrices[jCell] = MultidimensionalArray.Create(Np, Np);

                    // extract the normals in reference coordinates for the edge of the current pair
                    int face_i = _agGrd.AncestorGrid.Edges.FaceIndices[pairs[i, 2], Array.IndexOf(_agGrd.AncestorGrid.Edges.CellIndices.GetRow(pairs[i, 2]), pairs[i, 0])];
                    double[] normal_i = _maxDgBasis.GridDat.iGeomCells.GetRefElement(pairs[i, 0]).FaceNormals.ExtractSubArrayShallow(face_i, -1).To1DArray();


                    int face_j = _agGrd.AncestorGrid.Edges.FaceIndices[pairs[i, 2], Array.IndexOf(_agGrd.AncestorGrid.Edges.CellIndices.GetRow(pairs[i, 2]), pairs[i, 1])];
                    double[] normal_j = _maxDgBasis.GridDat.iGeomCells.GetRefElement(pairs[i, 1]).FaceNormals.ExtractSubArrayShallow(face_j, -1).To1DArray();

                    // compute the angle between the two cells
                    double alpha = Math.Acos(GenericBlas.InnerProd(normal_i, normal_j) / (normal_i.L2Norm() * normal_j.L2Norm()));
                    // direction???
                    alpha = Math.PI - alpha;

                    // compute the total rotation Matrix through multiplication
                    double k1 = Math.Cos(alpha);
                    double k2 = Math.Sin(alpha);

                    // Pattern of R: Note for alpha % 2PI == 0 this is a Eye matrix
                    // Polynomials: 1, y, x, y^2, xy, x^2 ...
                    // 1    0   0   0   0   0
                    // 0    k1 -k2  0   0   0
                    // 0    k2  k1  0   0   0
                    // 0    0   0   k1  0  -k2
                    // 0    0   0   0   1   0
                    // 0    0   0   k2  0   k1
                    // and so on
                    int start = 0;
                    for(int n = 0; n <= _maxDgBasis.Degree; n++) {
                        start += n;
                        MultidimensionalArray rot = RotationMatrices[jCell].ExtractSubArrayShallow(new int[] { start, start }, new int[] { start + n, start + n });
                        if(n % 2 == 0) {
                            rot[n / 2, n / 2] = 1;
                            for(int m = 0; m < n - 1; m++) {
                                rot[m, m] = k1;
                                rot[n - m, n - m] = k1;

                                rot[m, n - m] = -k2;
                                rot[n - m, m] = k2;
                            }
                        } else {
                            for(int m = 0; m < n; m++) {
                                rot[m, m] = k1;
                                rot[n - m, n - m] = k1;

                                rot[m, n - m] = -k2;
                                rot[n - m, m] = k2;
                            }
                        }
                    }

                    RotationMatrices[jCell].GEMM(1.0, RotationMatrices[iCell], RotationMatrices[jCell].CloneAs(), 0.0);
                }

            }

            return RotationMatrices;
        }

        // starts with a reference cell and extracts a graph of neigbour cells for the aggregated cell
        private static int[,] ExtractCellPairs(AggregationGridData agGrd, int[] parts) {
            int[,] pairs = new int[parts.Length - 1, 3];

            int k = 0;
            // create a list to build the graph
            List<int> excluded = new List<int>();
            excluded.Add(parts[0]);

            for(int i = 0; i < parts.Length - 1; i++) {
                // get inner neighbours that are not yet part of the graph
                int[] neighbours = agGrd.AncestorGrid.iLogicalCells.CellNeighbours[excluded[i]].Where(t => parts.Except(excluded).Contains(t)).ToArray();

                foreach(int n in neighbours) {
                    // write cellpair
                    pairs[k, 0] = excluded[i];
                    pairs[k, 1] = n;
                    // add the connecting edge (for one cell the edge is marked by a negative index) note -1 for correct array address
                    pairs[k, 2] = Math.Abs(agGrd.AncestorGrid.iLogicalCells.Cells2Edges[excluded[i]].Where(t => agGrd.AncestorGrid.iLogicalCells.Cells2Edges[n].Contains(-t)).First()) - 1;
                    // Now n is part of the graph
                    excluded.Add(n);
                    k++;
                }
            }

            return pairs;
        }

        // computes the injector coefficients for the first basis functions
        private static MultidimensionalArray GetCoefficients0Degree(AggregationGridData _agGrd, Basis _maxDgBasis, int[] parts) {
            var Injector0Degree = MultidimensionalArray.Create(parts.Length);

            // calculate total aggregation cell volume
            double aggVol = 0;
            foreach(int i in parts) {
                aggVol += _agGrd.AncestorGrid.Cells.GetCellVolume(i);
            }

            // For the first basis functions with constant value this breaks down to a simple scaling
            double scale = 1 / Math.Sqrt(aggVol);

            int count = 0;
            foreach(int i in parts) {
                Injector0Degree[count] = Math.Sqrt(_agGrd.AncestorGrid.Cells.GetCellVolume(i)) * scale;
                count++;
            }

            return Injector0Degree;
        }

        // computes the injector coefficients higher order basis functions
        private static MultidimensionalArray GetCoefficients(AggregationGridData _agGrd, Basis _maxDgBasis, MultidimensionalArray _inj, MultidimensionalArray[] _rot, int[] parts, int[] iedges, int basisIndex) {

            #region startup

            int partCount = parts.Length;
            int varCount = partCount * (basisIndex + 1);
            int currentDegree = _maxDgBasis.Polynomials[0][basisIndex].AbsoluteDegree;
            var InjectorNDegree = MultidimensionalArray.Create(partCount, basisIndex + 1);

            // equation matrix
            // basisIndex + 1 equations for othogonality
            // (currentdegree + 1) * iedges.Length equations for continuity along inneredges
            // (basisIndex + 1) * iedges.Length equations for edge-normal gradient jump minimization (smootheness)
            MultidimensionalArray aggE1 = MultidimensionalArray.Create(basisIndex + 1 + (currentDegree + 1) * iedges.Length, varCount);
            MultidimensionalArray aggE2 = MultidimensionalArray.Create((basisIndex * 2) * iedges.Length, varCount);
            MultidimensionalArray aggE2T = MultidimensionalArray.Create(varCount, (basisIndex * 2) * iedges.Length);


            // rhs
            MultidimensionalArray aggRhs1 = MultidimensionalArray.Create(basisIndex + 1 + (currentDegree + 1) * iedges.Length);
            MultidimensionalArray aggRhs2 = MultidimensionalArray.Create((basisIndex * 2) * iedges.Length);

            // construct the system
            // row[basisIndex] is a placeholder for the normality condition, to ensure full rank of the system
            //aggE1[basisIndex, (partCount/2) * (basisIndex + 1) - 1] = 1;
            aggE1[basisIndex, basisIndex] = 1;
            aggRhs1[basisIndex] = 1;

            #endregion

            #region orthogonality conditions

            // orthogonality
            for(int row = 0; row < basisIndex; row++) {
                for(int i = 0; i < partCount; i++) {
                    for(int column = 0; column <= row; column++) {
                        aggE1[row, i * (basisIndex + 1) + column] = _inj[i, column, row];
                        //E[row, i * (basisIndex + 1) + column] = _inj[i, column, row];
                    }
                }
            }

            #endregion

            #region continuity conditions

            // continuity
            // TODO harder part: Select currentDegree + 1 points per shared edge and ensure continuity there
            // 1. per edge construct currentDegree + 1 points
            // 2. evaluate the basisIndex -basisfunction on the adjacent cells for these points
            // 3. write the values in aggEq
            int rowOffset = basisIndex + 1;
            List<int> usableRowsE = Enumerable.Range(0, basisIndex + 1).ToList();


            for(int row = 0; row < iedges.Length; row++) {

                // FOR THE WHOLE SEGMENT, SUPER CLUNKY RIGHT NOW NEEDS TO BE WORKED OVER IF CODE IS FUNCTIONAL
                // select current edge
                int edge = iedges[row];

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[iedges[row], 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[iedges[row], 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // get basis values on the edge
                // construct evaluation points, currently only viable for 1D edges                
                // create NodeSet
                int edgeRef = _agGrd.iGeomEdges.GetRefElementIndex(edge);
                var refEl = _agGrd.iGeomEdges.EdgeRefElements[edgeRef];
                int edgeVertCount = refEl.NoOfVertices;

                MultidimensionalArray edgeEvalPoints = MultidimensionalArray.Create(currentDegree + 1, refEl.SpatialDimension);
                refEl.Vertices.To2DArray();
                if(edgeVertCount < currentDegree + 1) {
                    for(int i = 0; i < edgeVertCount; i++) {
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(1.0, refEl.Vertices.ExtractSubArrayShallow(i, -1));
                    }

                    for(int i = edgeVertCount; i < currentDegree + 1; i++) {
                        // construct a point as middle point between the edge vertices
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(0.5, edgeEvalPoints.ExtractSubArrayShallow(i - edgeVertCount, -1));
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(0.5, edgeEvalPoints.ExtractSubArrayShallow(i - edgeVertCount + 1, -1));
                    }
                } else {
                    for(int i = 0; i < currentDegree + 1; i++) {
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(1.0, refEl.Vertices.ExtractSubArrayShallow(i, -1));
                    }
                }


                NodeSet edge_coords = new NodeSet(refEl, edgeEvalPoints, false);
                edge_coords.LockForever();

                // evaluate the polynomials on the edge for both cells
                Tuple<MultidimensionalArray, MultidimensionalArray> edge_values = _maxDgBasis.EdgeEval(edge_coords, edge, 1);

                // apply the rotation
                edge_values.Item1.Multiply(1.0, edge_values.Item1.CloneAs(), _rot[index_i], 0.0, "cij", "cik", "kj");
                edge_values.Item2.Multiply(1.0, edge_values.Item2.CloneAs(), _rot[index_j], 0.0, "cij", "cik", "kj");

                for(int i = 0; i < (currentDegree + 1); i++) {

                    for(int column = 0; column < basisIndex + 1; column++) {
                        aggE1[row * (currentDegree + 1) + i + rowOffset, index_i * (basisIndex + 1) + column] = edge_values.Item1[0, i, column];
                        aggE1[row * (currentDegree + 1) + i + rowOffset, index_j * (basisIndex + 1) + column] = -edge_values.Item2[0, i, column];
                        //E[row * (currentDegree + 1) + i + rowOffset, index_i * (basisIndex + 1) + column] = edge_values.Item1[0, i, column];
                        //E[row * (currentDegree + 1) + i + rowOffset, index_j * (basisIndex + 1) + column] = -edge_values.Item2[0, i, column];
                    }
                }

                // test if rows are linear dependent
                usableRowsE.AddRange(CheckLinDependency(aggE1, currentDegree + 1, row * (currentDegree + 1) + rowOffset));
            }

            #endregion

            #region gradient jump conditions

            // smootheness
            // TODO hardest part: Evaluate the gradients on the edge in normal direction then compute the L2-norm of the jump and enforce a minimization condition
            // 1. compute the gradient parts for the relevant (0-basisIndex) basis functions in normal direction on each edge
            // 2. compute the Gram determinant for the edge
            // 3. integrate the basisfunction couplings (originating from the L2-norm) over the edge (Caching would be highly advisable, because the same integrals will be needed for the next basis function to aggregate)
            // 4. write values in aggEq

            // create array containing the edge-normal gradient coupling integrals
            MultidimensionalArray gradient = MultidimensionalArray.Create(2, basisIndex, 2, basisIndex);
            List<int> usableRowsQ = new List<int>();
            // iterate over edges
            for(int row = 0; row < iedges.Length; row++) {
                // FOR THE WHOLE SEGMENT, SUPER CLUNKY RIGHT NOW NEEDS TO BE WORKED OVER IF CODE IS FUNCTIONAL
                // select current edge
                int edge = iedges[row];

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[iedges[row], 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[iedges[row], 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // create array containing the edge-normal gradient coupling integrals
                gradient.Clear();
                //gradient = GetEdgeGradient(_maxDgBasis, _agGrd, edge, cell_i, cell_j, basisIndex);

                // compute the gradient coupled integrals
                var Context = _maxDgBasis.GridDat;
                var edgeMask = new EdgeMask(Context, new int[] { edge }, MaskType.Geometrical);

                EdgeQuadrature.GetQuadrature(new int[4] { 2, basisIndex, 2, basisIndex }, Context,
                    (new EdgeQuadratureScheme(true, edgeMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes_ref_edge = QR.Nodes;
                        //Debug.Assert(Length == 1);

                        MultidimensionalArray edge_normals = Context.iGeomEdges.NormalsCache.GetNormals_Edge(nodes_ref_edge, edge, 1);
                        Tuple<MultidimensionalArray, MultidimensionalArray> edge_gradient = _maxDgBasis.EdgeEvalGradient(nodes_ref_edge, edge, 1);

                        // apply the rotation
                        edge_gradient.Item1.Multiply(1.0, edge_gradient.Item1.CloneAs(), _rot[index_i], 0.0, "cijn", "cikn", "kj");
                        edge_gradient.Item2.Multiply(1.0, edge_gradient.Item2.CloneAs(), _rot[index_j], 0.0, "cijn", "cikn", "kj");

                        MultidimensionalArray edge_normal_gradient = MultidimensionalArray.Create(2, nodes_ref_edge.Length, basisIndex);

                        // compute the gradients in normal direction at the quadnodes reverse the normal for one cell
                        edge_normal_gradient.ExtractSubArrayShallow(0, -1, -1).Multiply(1.0, edge_gradient.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");
                        edge_normal_gradient.ExtractSubArrayShallow(1, -1, -1).Multiply(-1.0, edge_gradient.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");

                        // couple the gradients at each node (out "kinjm" "node_k,cell_i,polynom_n,cell_j,polynom_m")
                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1, -1, -1);
                        EvalResult.Multiply(1.0, edge_normal_gradient, edge_normal_gradient, 0.0, "kinjm", "ikn", "jkm");

                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1, -1, -1);
                                                     gradient.Clear();
                                                     gradient.Acc(1.0, res);
                                                 }).Execute();



                // iterate over both cells (leading coeffcient)
                for(int i = 0; i <= 1; i++) {

                    // iterate over all conditions for the current cell
                    for(int k = 1; k < basisIndex + 1; k++) {

                        // iterate over the variables
                        for(int column = 1; column < basisIndex + 1; column++) {
                            aggE2[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_i * (basisIndex + 1) + column] = gradient[i, k - 1, 0, column - 1];
                            aggE2[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_j * (basisIndex + 1) + column] = gradient[i, k - 1, 1, column - 1];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_i * (basisIndex + 1) + column] = gradient[i, k - 1, 0, column - 1];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_j * (basisIndex + 1) + column] = gradient[i, k - 1, 1, column - 1];
                        }
                    }
                }
                // test if rows are linear dependent
                usableRowsQ.AddRange(CheckLinDependency(aggE2, 2 * basisIndex, row * (basisIndex * 2)));
            }

            #endregion

            #region solve system

            //Partitioning rowPartA = new Partitioning((basisIndex * 2) * iedges.Length);
            //Partitioning rowPartE = new Partitioning(basisIndex + 1 + (currentDegree + 1) * iedges.Length);
            //Partitioning colPart = new Partitioning(varCount);

            //// Matrix storing the optimization goal
            //MsrMatrix A = new MsrMatrix(rowPartA, colPart);
            //// Matrix storing the side conditions
            ////MsrMatrix E = new MsrMatrix(rowPartE, colPart);

            // construct the complete matrix for the quadratic program
            MultidimensionalArray E = MultidimensionalArray.Create(usableRowsE.Count(), varCount);
            MultidimensionalArray ET = MultidimensionalArray.Create(varCount, usableRowsE.Count());

            for(int i = 0; i < usableRowsE.Count(); i++) {
                E.SetRow(i, aggE1.ExtractSubArrayShallow(usableRowsE.ElementAt(i), -1).To1DArray());
                ET.SetColumn(i, aggE1.ExtractSubArrayShallow(usableRowsE.ElementAt(i), -1).To1DArray());
            }

            MultidimensionalArray hQ = MultidimensionalArray.Create(usableRowsQ.Count(), varCount);

            for(int i = 0; i < usableRowsQ.Count(); i++) {
                hQ.SetRow(i, aggE2.ExtractSubArrayShallow(usableRowsQ.ElementAt(i), -1).To1DArray());
            }

            // solution using quadratic program and sparse matrices
            MultidimensionalArray Q = MultidimensionalArray.Create(hQ.NoOfCols, hQ.NoOfCols);
            double Qscale = 1.0; // 1 / Math.Pow(aggE2.InfNorm(), 2.0);
            Q.GEMM(Qscale, hQ, hQ, 0.0, true);

            Partitioning qRowPart = new Partitioning(Q.NoOfRows + E.NoOfRows, MPI.Wrappers.csMPI.Raw._COMM.SELF);
            Partitioning qColPart = new Partitioning(Q.NoOfCols + E.NoOfRows, MPI.Wrappers.csMPI.Raw._COMM.SELF);
            MsrMatrix QP = new MsrMatrix(qRowPart, qColPart);

            double[] c = new double[varCount];
            double[] d = new double[E.NoOfRows];

            d[basisIndex] = 1;
            QP.AccBlock(0, 0, 1.0, Q);
            QP.AccBlock(0, Q.NoOfCols, 1.0, ET);
            QP.AccBlock(Q.NoOfRows, 0, 1.0, E);



            double[] RHS = c.Concat(d).ToArray();
            double[] v = new double[QP.NoOfCols];

            var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
            //var solver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            solver.DefineMatrix(QP);
            solver.Solve(v, RHS);

            // compute the solution that fullfills the orthogonality and continuity exactly and the minimizes the gradient jumps
            double[] sol = new double[varCount];
            //aggE1.SolveWithCondition(sol, aggRhs1.To1DArray(), aggE2, aggRhs2.To1DArray());
            sol = v.GetSubVector(0, varCount);


            // check sparsity
            //int nonzeroCount1 = 0;
            //foreach (double t in aggE1.To2DArray()) {
            //    nonzeroCount1 += t != 0.0 ? 1 : 0; 
            //}
            //int nonzeroCount2 = 0;
            //foreach (double t in aggE2.To2DArray())
            //{
            //    nonzeroCount2 += t != 0.0 ? 1 : 0;
            //}

            //double sparsity1 = (double)nonzeroCount1 / (aggE1.NoOfRows * aggE1.NoOfCols);
            //double sparsity2 =  (double)nonzeroCount2 / (aggE2.NoOfRows * aggE2.NoOfCols);
            //Console.WriteLine($"Sparsity of first equation set: {sparsity1}; Sparsity of second equation set: {sparsity2}");

            // compute the diagonal entry of the Injector for the first cell in parts and then solve the dependencies to the other cells
            // normality condition
            double scale = 1 / sol.L2Norm();
            sol.ScaleV(scale);
            MultidimensionalArray vector_InjectorNDegree = MultidimensionalArray.CreateWrapper(sol, varCount);
            //vector_InjectorNDegree.Scale(scale);

            #endregion


            #region fmincon

            // check the diagonal entries
            bool isPD = true;
            for(int i = 0; i < parts.Length; i++) {
                if(Math.Abs(vector_InjectorNDegree[(i + 1) * (basisIndex + 1) - 1]) <= 1E-5) { isPD = false; }
            }

            // if diagonal entries are too small use fmincon
            if(!isPD) {
                Console.WriteLine("Using fmincon to ensure Positive Definiteness");
                // solve with fmincon as optimization problem, note the 3 additional scripts already present in the working directory
                BatchmodeConnector bmc = new BatchmodeConnector(@"D:\bosss_db_masterthesis\Matlab");
                bmc.Cmd("clear all");

                // gather the necessary matrices
                bmc.PutMatrix(Q, "Q");

                // quadratic constraints (note these are nonconvex) 0.5*x'H_ix + d_i <= 0
                MultidimensionalArray[] H = new MultidimensionalArray[parts.Length];
                double lowerlimit = Math.Sqrt((double)1 / parts.Length * 1E-8);
                string infimum = $@"1/{parts.Length}*1E-8";
                bmc.Cmd($"H = cell({parts.Length},1);");
                bmc.Cmd($"d = cell({parts.Length},1);");
                for(int i = 0; i < parts.Length; i++) {
                    H[i] = MultidimensionalArray.Create(varCount, varCount);
                    H[i][(i + 1) * (basisIndex + 1) - 1, (i + 1) * (basisIndex + 1) - 1] = -1;
                    bmc.PutMatrix(H[i], $"H{{{i + 1}}}");
                    bmc.Cmd($"d{{{i + 1}}} = {infimum};");
                }

                // equality constraint Gx = 0
                bmc.PutMatrix(E, "G");
                // remove the fix condition
                bmc.Cmd($"G({basisIndex + 1},:) = [];");
                bmc.Cmd($"G = sparse(G);");
                bmc.Cmd("b = sparse(size(G,1),1);");

                // ensure normality 1/2x'Heqx + d = 0
                bmc.Cmd("Heq = cell(1,1);");
                bmc.Cmd("deq = cell(1,1);");

                bmc.Cmd("Heq{1} = 2*eye(size(G,2));");
                bmc.Cmd("deq{1} = -1;");



                bmc.Cmd(@"options = optimoptions(@fmincon,'Algorithm','interior-point',...
                    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
                    'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H,Heq));");

                bmc.Cmd("fun = @(x)quadobj(x,Q);");
                bmc.Cmd("nonlconstr = @(x)quadconstr(x,H,d,Heq,deq);");
                bmc.Cmd("x0 = sparse(size(G,2),1);");

                bmc.Cmd("flin = sparse(size(x0,1),1);");
                bmc.Cmd("fnlin = @(x)quadobj(x,sparse(size(Q,1),size(Q,2)));");

                bmc.Cmd("x0 = linprog(flin,[],[],G,b,[],[]);");
                bmc.Cmd(@"x0 = fmincon(fnlin,x0,...
                    [],[], G, b,[],[], nonlconstr, options); ");
                bmc.Cmd(@"[x,fval,eflag,output,lambda] = fmincon(fun,x0,...
                    [],[], G, b,[],[], nonlconstr, options); ");


                MultidimensionalArray vector_InjectorNDegree_Matlab = MultidimensionalArray.Create(varCount, 1);
                bmc.GetMatrix(vector_InjectorNDegree_Matlab, "x");

                bmc.Execute(false);

                vector_InjectorNDegree = vector_InjectorNDegree_Matlab.ExtractSubArrayShallow(-1, 0);

                // simple test to at least get an idea if somethings gone horribly wrong
                if(Math.Abs(vector_InjectorNDegree_Matlab.L2Norm() - 1) >= 1E-3) Console.WriteLine("Warning converged to infeasable point");

                isPD = true;
                for(int i = 0; i < parts.Length; i++) {
                    if(Math.Abs(vector_InjectorNDegree[(i + 1) * (basisIndex + 1) - 1]) - lowerlimit <= 0) { isPD = false; }
                }

                if(!isPD) Console.WriteLine($"WARNING: Basisfunction {basisIndex} violates pd condition");
            }
            #endregion

            // restructure to Injector operator form
            for(int i = 0; i < parts.Length; i++) {
                InjectorNDegree.ExtractSubArrayShallow(i, -1).Acc(1.0, vector_InjectorNDegree.ExtractSubArrayShallow(new int[] { i * (basisIndex + 1) }, new int[] { (i + 1) * (basisIndex + 1) - 1 }));
            }

            return InjectorNDegree;
        }

        // computes the injector coefficients higher order basis functions, using integral conditions for jump and gradient
        private static MultidimensionalArray GetCoefficientsIntegral(AggregationGridData _agGrd, Basis _maxDgBasis, MultidimensionalArray _inj, MultidimensionalArray[] _rot, int[] parts, int[] iedges, int basisIndex) {

            #region startup

            int partCount = parts.Length;
            int varCount = partCount * (basisIndex + 1);
            int currentDegree = _maxDgBasis.Polynomials[0][basisIndex].AbsoluteDegree;
            var InjectorNDegree = MultidimensionalArray.Create(partCount, basisIndex + 1);

            // equation matrix
            // basisIndex + 1 equations for othogonality
            // (currentdegree + 1) * iedges.Length equations for continuity along inneredges
            // (basisIndex + 1) * iedges.Length equations for edge-normal gradient jump minimization (smootheness)
            MultidimensionalArray aggE1 = MultidimensionalArray.Create(basisIndex + 1 + 2 * (basisIndex + 1) * iedges.Length, varCount);
            MultidimensionalArray aggE2 = MultidimensionalArray.Create((basisIndex * 2) * iedges.Length, varCount);
            MultidimensionalArray aggE2T = MultidimensionalArray.Create(varCount, (basisIndex * 2) * iedges.Length);


            // rhs
            MultidimensionalArray aggRhs1 = MultidimensionalArray.Create(basisIndex + 1 + (currentDegree + 1) * iedges.Length);
            MultidimensionalArray aggRhs2 = MultidimensionalArray.Create((basisIndex * 2) * iedges.Length);

            // construct the system
            // row[basisIndex] is a placeholder for the normality condition, to ensure full rank of the system
            //aggE1[basisIndex, (partCount/2) * (basisIndex + 1) - 1] = 1;
            int cellOffset = basisIndex % parts.Length;
            aggE1[basisIndex, basisIndex] = 1;
            //aggE1[basisIndex, (cellOffset + 1) * (basisIndex + 1) - 1] = 1;
            aggRhs1[basisIndex] = 1;

            #endregion

            #region orthogonality conditions

            // orthogonality
            for(int row = 0; row < basisIndex; row++) {
                for(int i = 0; i < partCount; i++) {
                    for(int column = 0; column <= row; column++) {
                        aggE1[row, i * (basisIndex + 1) + column] = _inj[i, column, row];
                        //E[row, i * (basisIndex + 1) + column] = _inj[i, column, row];
                    }
                }
            }

            #endregion

            #region continuity conditions

            // continuity
            // TODO harder part: Select currentDegree + 1 points per shared edge and ensure continuity there
            // 1. per edge construct currentDegree + 1 points
            // 2. evaluate the basisIndex -basisfunction on the adjacent cells for these points
            // 3. write the values in aggEq
            int rowOffset = basisIndex + 1;


            // create array containing the jump coupling integrals
            MultidimensionalArray jump = MultidimensionalArray.Create(2, basisIndex + 1, 2, basisIndex + 1);
            List<int> usableRowsE = Enumerable.Range(0, basisIndex + 1).ToList();
            // iterate over edges
            for(int row = 0; row < iedges.Length; row++) {
                // FOR THE WHOLE SEGMENT, SUPER CLUNKY RIGHT NOW NEEDS TO BE WORKED OVER IF CODE IS FUNCTIONAL
                // select current edge
                int edge = iedges[row];

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[iedges[row], 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[iedges[row], 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // create array containing the edge-normal gradient coupling integrals
                jump.Clear();

                // compute the jump coupled integrals
                var Context = _maxDgBasis.GridDat;
                var edgeMask = new EdgeMask(Context, new int[] { edge }, MaskType.Geometrical);

                EdgeQuadrature.GetQuadrature(new int[4] { 2, basisIndex + 1, 2, basisIndex + 1 }, Context,
                    (new EdgeQuadratureScheme(true, edgeMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes_ref_edge = QR.Nodes;
                        //Debug.Assert(Length == 1);

                        Tuple<MultidimensionalArray, MultidimensionalArray> edge_jump = _maxDgBasis.EdgeEval(nodes_ref_edge, edge, 1);

                        // apply the rotation
                        edge_jump.Item1.Multiply(1.0, edge_jump.Item1.CloneAs(), _rot[index_i], 0.0, "cij", "cik", "kj");
                        edge_jump.Item2.Multiply(1.0, edge_jump.Item2.CloneAs(), _rot[index_j], 0.0, "cij", "cik", "kj");

                        MultidimensionalArray edge_coupled_jump = MultidimensionalArray.Create(2, nodes_ref_edge.Length, basisIndex + 1);

                        // compute the gradients in normal direction at the quadnodes reverse the normal for one cell
                        edge_coupled_jump.ExtractSubArrayShallow(0, -1, -1).Acc(1.0, edge_jump.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex }));
                        edge_coupled_jump.ExtractSubArrayShallow(1, -1, -1).Acc(1.0, edge_jump.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex }));

                        // couple the gradients at each node (out "kinjm" "node_k,cell_i,polynom_n,cell_j,polynom_m")
                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1, -1, -1);
                        EvalResult.Multiply(1.0, edge_coupled_jump, edge_coupled_jump, 0.0, "kinjm", "ikn", "jkm");

                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1, -1, -1);
                                                     jump.Clear();
                                                     jump.Acc(1.0, res);
                                                 }).Execute();



                // iterate over both cells (leading coeffcient)
                for(int i = 0; i <= 1; i++) {

                    // iterate over all conditions for the current cell
                    for(int k = 0; k < basisIndex + 1; k++) {

                        // iterate over the variables
                        for(int column = 0; column < basisIndex + 1; column++) {

                            aggE1[rowOffset + row * (basisIndex + 1) * 2 + k + i * (basisIndex + 1), index_i * (basisIndex + 1) + column] = jump[i, k, 0, column];
                            aggE1[rowOffset + row * (basisIndex + 1) * 2 + k + i * (basisIndex + 1), index_j * (basisIndex + 1) + column] = -jump[i, k, 1, column];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_i * (basisIndex + 1) + column] = gradient[i, k - 1, 0, column - 1];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_j * (basisIndex + 1) + column] = gradient[i, k - 1, 1, column - 1];
                        }
                    }
                }
                // test if rows are linear dependent
                usableRowsE.AddRange(CheckLinDependency(aggE1, 2 * (basisIndex + 1), rowOffset + row * (basisIndex + 1) * 2));
            }

            #endregion

            #region gradient jump conditions

            // smootheness
            // TODO hardest part: Evaluate the gradients on the edge in normal direction then compute the L2-norm of the jump and enforce a minimization condition
            // 1. compute the gradient parts for the relevant (0-basisIndex) basis functions in normal direction on each edge
            // 2. compute the Gram determinant for the edge
            // 3. integrate the basisfunction couplings (originating from the L2-norm) over the edge (Caching would be highly advisable, because the same integrals will be needed for the next basis function to aggregate)
            // 4. write values in aggEq

            // create array containing the edge-normal gradient coupling integrals
            MultidimensionalArray gradient = MultidimensionalArray.Create(2, basisIndex, 2, basisIndex);
            List<int> usableRowsQ = new List<int>();
            // iterate over edges
            for(int row = 0; row < iedges.Length; row++) {
                // FOR THE WHOLE SEGMENT, SUPER CLUNKY RIGHT NOW NEEDS TO BE WORKED OVER IF CODE IS FUNCTIONAL
                // select current edge
                int edge = iedges[row];

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[iedges[row], 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[iedges[row], 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // create array containing the edge-normal gradient coupling integrals
                gradient.Clear();
                //gradient = GetEdgeGradient(_maxDgBasis, _agGrd, edge, cell_i, cell_j, basisIndex);

                // compute the gradient coupled integrals
                var Context = _maxDgBasis.GridDat;
                var edgeMask = new EdgeMask(Context, new int[] { edge }, MaskType.Geometrical);

                EdgeQuadrature.GetQuadrature(new int[4] { 2, basisIndex, 2, basisIndex }, Context,
                    (new EdgeQuadratureScheme(true, edgeMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes_ref_edge = QR.Nodes;
                        //Debug.Assert(Length == 1);

                        MultidimensionalArray edge_normals = Context.iGeomEdges.NormalsCache.GetNormals_Edge(nodes_ref_edge, edge, 1);
                        Tuple<MultidimensionalArray, MultidimensionalArray> edge_gradient = _maxDgBasis.EdgeEvalGradient(nodes_ref_edge, edge, 1);

                        // apply the rotation
                        edge_gradient.Item1.Multiply(1.0, edge_gradient.Item1.CloneAs(), _rot[index_i], 0.0, "cijn", "cikn", "kj");
                        edge_gradient.Item2.Multiply(1.0, edge_gradient.Item2.CloneAs(), _rot[index_j], 0.0, "cijn", "cikn", "kj");

                        MultidimensionalArray edge_normal_gradient = MultidimensionalArray.Create(2, nodes_ref_edge.Length, basisIndex);

                        // compute the gradients in normal direction at the quadnodes reverse the normal for one cell
                        edge_normal_gradient.ExtractSubArrayShallow(0, -1, -1).Multiply(1.0, edge_gradient.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");
                        edge_normal_gradient.ExtractSubArrayShallow(1, -1, -1).Multiply(-1.0, edge_gradient.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");

                        // couple the gradients at each node (out "kinjm" "node_k,cell_i,polynom_n,cell_j,polynom_m")
                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1, -1, -1);
                        EvalResult.Multiply(1.0, edge_normal_gradient, edge_normal_gradient, 0.0, "kinjm", "ikn", "jkm");

                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1, -1, -1);
                                                     gradient.Clear();
                                                     gradient.Acc(1.0, res);
                                                 }).Execute();



                // iterate over both cells (leading coeffcient)
                for(int i = 0; i <= 1; i++) {

                    // iterate over all conditions for the current cell
                    for(int k = 1; k < basisIndex + 1; k++) {

                        // iterate over the variables
                        for(int column = 1; column < basisIndex + 1; column++) {
                            aggE2[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_i * (basisIndex + 1) + column] = gradient[i, k - 1, 0, column - 1];
                            aggE2[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_j * (basisIndex + 1) + column] = gradient[i, k - 1, 1, column - 1];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_i * (basisIndex + 1) + column] = gradient[i, k - 1, 0, column - 1];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_j * (basisIndex + 1) + column] = gradient[i, k - 1, 1, column - 1];
                        }
                    }
                }
                // test if rows are linear dependent
                usableRowsQ.AddRange(CheckLinDependency(aggE2, 2 * basisIndex, row * (basisIndex * 2)));
            }

            #endregion

            #region solve system

            //Partitioning rowPartA = new Partitioning((basisIndex * 2) * iedges.Length);
            //Partitioning rowPartE = new Partitioning(basisIndex + 1 + (currentDegree + 1) * iedges.Length);
            //Partitioning colPart = new Partitioning(varCount);

            //// Matrix storing the optimization goal
            //MsrMatrix A = new MsrMatrix(rowPartA, colPart);
            //// Matrix storing the side conditions
            ////MsrMatrix E = new MsrMatrix(rowPartE, colPart);

            // construct the complete matrix for the quadratic program
            MultidimensionalArray E = MultidimensionalArray.Create(usableRowsE.Count(), varCount);
            MultidimensionalArray ET = MultidimensionalArray.Create(varCount, usableRowsE.Count());

            for(int i = 0; i < usableRowsE.Count(); i++) {
                E.SetRow(i, aggE1.ExtractSubArrayShallow(usableRowsE.ElementAt(i), -1).To1DArray());
                ET.SetColumn(i, aggE1.ExtractSubArrayShallow(usableRowsE.ElementAt(i), -1).To1DArray());
            }

            MultidimensionalArray hQ = MultidimensionalArray.Create(usableRowsQ.Count(), varCount);

            for(int i = 0; i < usableRowsQ.Count(); i++) {
                hQ.SetRow(i, aggE2.ExtractSubArrayShallow(usableRowsQ.ElementAt(i), -1).To1DArray());
            }

            // solution using quadratic program and sparse matrices
            MultidimensionalArray Q = MultidimensionalArray.Create(hQ.NoOfCols, hQ.NoOfCols);
            double Qscale = 1.0; // 1 / Math.Pow(aggE2.InfNorm(), 2.0);
            Q.GEMM(Qscale, hQ, hQ, 0.0, true);

            Partitioning qRowPart = new Partitioning(Q.NoOfRows + E.NoOfRows, MPI.Wrappers.csMPI.Raw._COMM.SELF);
            Partitioning qColPart = new Partitioning(Q.NoOfCols + E.NoOfRows, MPI.Wrappers.csMPI.Raw._COMM.SELF);
            MsrMatrix QP = new MsrMatrix(qRowPart, qColPart);

            double[] c = new double[varCount];
            double[] d = new double[E.NoOfRows];

            d[basisIndex] = 1;

            QP.AccBlock(0, 0, 1.0, Q);
            QP.AccBlock(0, Q.NoOfCols, 1.0, ET);
            QP.AccBlock(Q.NoOfRows, 0, 1.0, E);
            QP.AssumeSymmetric = true;


            double[] RHS = c.Concat(d).ToArray();
            double[] v = new double[QP.NoOfCols];


            #region Sparse Solve
            //var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
            var solver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            solver.DefineMatrix(QP);
            solver.Solve(v, RHS);

            // compute the solution that fullfills the orthogonality and continuity exactly and the minimizes the gradient jumps
            double[] sol = new double[varCount];
            //aggE1.SolveWithCondition(sol, aggRhs1.To1DArray(), aggE2, aggRhs2.To1DArray());
            sol = v.GetSubVector(0, varCount);


            // check sparsity
            //int nonzeroCount1 = 0;
            //foreach (double t in aggE1.To2DArray()) {
            //    nonzeroCount1 += t != 0.0 ? 1 : 0; 
            //}
            //int nonzeroCount2 = 0;
            //foreach (double t in aggE2.To2DArray())
            //{
            //    nonzeroCount2 += t != 0.0 ? 1 : 0;
            //}

            //double sparsity1 = (double)nonzeroCount1 / (aggE1.NoOfRows * aggE1.NoOfCols);
            //double sparsity2 =  (double)nonzeroCount2 / (aggE2.NoOfRows * aggE2.NoOfCols);
            //Console.WriteLine($"Sparsity of first equation set: {sparsity1}; Sparsity of second equation set: {sparsity2}");

            // compute the diagonal entry of the Injector for the first cell in parts and then solve the dependencies to the other cells
            // normality condition
            double scale = 1 / sol.L2Norm();
            sol.ScaleV(scale);
            MultidimensionalArray vector_InjectorNDegree = MultidimensionalArray.CreateWrapper(sol, varCount);
            #endregion


            #region fmincon

            // check the diagonal entries
            bool isPD = true;
            for(int i = 0; i < parts.Length; i++) {
                if(Math.Abs(vector_InjectorNDegree[(i + 1) * (basisIndex + 1) - 1]) <= 1E-5) { isPD = false; }
            }

            // if diagonal entries are too small use fmincon
            if(!isPD) {
                Console.WriteLine("Using fmincon to ensure Positive Definiteness");
                // solve with fmincon as optimization problem, note the 3 additional scripts already present in the working directory
                BatchmodeConnector bmc = new BatchmodeConnector(@"D:\bosss_db_masterthesis\Matlab");
                bmc.Cmd("clear all");

                // gather the necessary matrices
                bmc.PutMatrix(Q, "Q");

                // quadratic constraints (note these are nonconvex) 0.5*x'H_ix + d_i <= 0
                MultidimensionalArray[] H = new MultidimensionalArray[parts.Length];
                double lowerlimit = Math.Sqrt((double)1 / parts.Length * 1E-8);
                string infimum = $@"1/{parts.Length}*1E-8";
                bmc.Cmd($"H = cell({parts.Length},1);");
                bmc.Cmd($"d = cell({parts.Length},1);");
                for(int i = 0; i < parts.Length; i++) {
                    H[i] = MultidimensionalArray.Create(varCount, varCount);
                    H[i][(i + 1) * (basisIndex + 1) - 1, (i + 1) * (basisIndex + 1) - 1] = -1;
                    bmc.PutMatrix(H[i], $"H{{{i + 1}}}");
                    bmc.Cmd($"d{{{i + 1}}} = {infimum};");
                }

                // equality constraint Gx = 0
                bmc.PutMatrix(E, "G");
                // remove the fix condition
                bmc.Cmd($"G({basisIndex + 1},:) = [];");
                bmc.Cmd($"G = sparse(G);");
                bmc.Cmd("b = sparse(size(G,1),1);");

                // ensure normality 1/2x'Heqx + d = 0
                bmc.Cmd("Heq = cell(1,1);");
                bmc.Cmd("deq = cell(1,1);");

                bmc.Cmd("Heq{1} = 2*eye(size(G,2));");
                bmc.Cmd("deq{1} = -1;");



                bmc.Cmd(@"options = optimoptions(@fmincon,'Algorithm','interior-point',...
                    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
                    'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H,Heq));");

                bmc.Cmd("fun = @(x)quadobj(x,Q);");
                bmc.Cmd("nonlconstr = @(x)quadconstr(x,H,d,Heq,deq);");
                bmc.Cmd("x0 = sparse(size(G,2),1);");

                bmc.Cmd("flin = sparse(size(x0,1),1);");
                bmc.Cmd("fnlin = @(x)quadobj(x,sparse(size(Q,1),size(Q,2)));");

                bmc.Cmd("x0 = linprog(flin,[],[],G,b,[],[]);");
                bmc.Cmd(@"x0 = fmincon(fnlin,x0,...
                    [],[], G, b,[],[], nonlconstr, options); ");
                bmc.Cmd(@"[x,fval,eflag,output,lambda] = fmincon(fun,x0,...
                    [],[], G, b,[],[], nonlconstr, options); ");


                MultidimensionalArray vector_InjectorNDegree_Matlab = MultidimensionalArray.Create(varCount, 1);
                bmc.GetMatrix(vector_InjectorNDegree_Matlab, "x");

                bmc.Execute(false);

                vector_InjectorNDegree = vector_InjectorNDegree_Matlab.ExtractSubArrayShallow(-1, 0);

                // simple test to at least get an idea if somethings gone horribly wrong
                if(Math.Abs(vector_InjectorNDegree_Matlab.L2Norm() - 1) >= 1E-3) Console.WriteLine("Warning converged to infeasable point");

                isPD = true;
                for(int i = 0; i < parts.Length; i++) {
                    if(Math.Abs(vector_InjectorNDegree[(i + 1) * (basisIndex + 1) - 1]) - lowerlimit <= 0) { isPD = false; }
                }

                if(!isPD) Console.WriteLine($"WARNING: Basisfunction {basisIndex} violates pd condition");
            }
            #endregion

            #endregion



            // restructure to Injector operator form
            for(int i = 0; i < parts.Length; i++) {
                InjectorNDegree.ExtractSubArrayShallow(i, -1).Acc(1.0, vector_InjectorNDegree.ExtractSubArrayShallow(new int[] { i * (basisIndex + 1) }, new int[] { (i + 1) * (basisIndex + 1) - 1 }));
            }

            return InjectorNDegree;
        }

        // looks for a solution with only 0 degree and diagonal part in the injectors
        private static MultidimensionalArray GetCoefficientsGlobalSF(AggregationGridData _agGrd, Basis _maxDgBasis, MultidimensionalArray _inj, MultidimensionalArray[] _rot, int[] parts, int[] iedges, int basisIndex) {

            #region startup

            int partCount = parts.Length;
            int varCount = partCount * 2;
            int currentDegree = _maxDgBasis.Polynomials[0][basisIndex].AbsoluteDegree;
            var InjectorNDegree = MultidimensionalArray.Create(partCount, basisIndex + 1);

            // equation matrix
            // basisIndex + 1 equations for othogonality
            // (currentdegree + 1) * iedges.Length equations for continuity along inneredges
            // (basisIndex + 1) * iedges.Length equations for edge-normal gradient jump minimization (smootheness)
            MultidimensionalArray aggE1 = MultidimensionalArray.Create(2 + 4 * iedges.Length, varCount);
            MultidimensionalArray aggE2 = MultidimensionalArray.Create(2 * iedges.Length, varCount);
            MultidimensionalArray aggE2T = MultidimensionalArray.Create(varCount, 2 * iedges.Length);


            // rhs
            MultidimensionalArray aggRhs1 = MultidimensionalArray.Create(2 + 4 * iedges.Length);
            MultidimensionalArray aggRhs2 = MultidimensionalArray.Create(2 * iedges.Length);

            // construct the system
            // row[basisIndex] is a placeholder for the normality condition, to ensure full rank of the system
            //aggE1[basisIndex, (partCount/2) * (basisIndex + 1) - 1] = 1;
            //int cellOffset = basisIndex % parts.Length;
            aggE1[0, 0] = 1;
            aggE1[1, 1] = 1;
            //aggE1[basisIndex, (cellOffset + 1) * (basisIndex + 1) - 1] = 1;
            aggRhs1[1] = 1;

            #endregion

            #region orthogonality conditions

            //// orthogonality
            //for (int row = 0; row < basisIndex; row++)
            //{
            //    for (int i = 0; i < partCount; i++)
            //    {
            //        for (int column = 0; column <= row; column++)
            //        {
            //            aggE1[row, i * (basisIndex + 1) + column] = _inj[i, column, row];
            //            //E[row, i * (basisIndex + 1) + column] = _inj[i, column, row];
            //        }
            //    }
            //}

            #endregion

            #region continuity conditions

            // continuity
            // TODO harder part: Select currentDegree + 1 points per shared edge and ensure continuity there
            // 1. per edge construct currentDegree + 1 points
            // 2. evaluate the basisIndex -basisfunction on the adjacent cells for these points
            // 3. write the values in aggEq
            int rowOffset = 2;


            // create array containing the jump coupling integrals
            MultidimensionalArray jump = MultidimensionalArray.Create(2, basisIndex + 1, 2, basisIndex + 1);
            List<int> usableRowsE = Enumerable.Range(0, 2).ToList();
            // iterate over edges
            for(int row = 0; row < iedges.Length; row++) {
                // FOR THE WHOLE SEGMENT, SUPER CLUNKY RIGHT NOW NEEDS TO BE WORKED OVER IF CODE IS FUNCTIONAL
                // select current edge
                int edge = iedges[row];

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[iedges[row], 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[iedges[row], 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // create array containing the edge-normal gradient coupling integrals
                jump.Clear();

                // compute the jump coupled integrals
                var Context = _maxDgBasis.GridDat;
                var edgeMask = new EdgeMask(Context, new int[] { edge }, MaskType.Geometrical);

                EdgeQuadrature.GetQuadrature(new int[4] { 2, basisIndex + 1, 2, basisIndex + 1 }, Context,
                    (new EdgeQuadratureScheme(true, edgeMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes_ref_edge = QR.Nodes;
                        //Debug.Assert(Length == 1);

                        Tuple<MultidimensionalArray, MultidimensionalArray> edge_jump = _maxDgBasis.EdgeEval(nodes_ref_edge, edge, 1);

                        // apply the rotation
                        edge_jump.Item1.Multiply(1.0, edge_jump.Item1.CloneAs(), _rot[index_i], 0.0, "cij", "cik", "kj");
                        edge_jump.Item2.Multiply(1.0, edge_jump.Item2.CloneAs(), _rot[index_j], 0.0, "cij", "cik", "kj");

                        MultidimensionalArray edge_coupled_jump = MultidimensionalArray.Create(2, nodes_ref_edge.Length, basisIndex + 1);

                        // compute the gradients in normal direction at the quadnodes reverse the normal for one cell
                        edge_coupled_jump.ExtractSubArrayShallow(0, -1, -1).Acc(1.0, edge_jump.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex }));
                        edge_coupled_jump.ExtractSubArrayShallow(1, -1, -1).Acc(1.0, edge_jump.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex }));

                        // couple the gradients at each node (out "kinjm" "node_k,cell_i,polynom_n,cell_j,polynom_m")
                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1, -1, -1);
                        EvalResult.Multiply(1.0, edge_coupled_jump, edge_coupled_jump, 0.0, "kinjm", "ikn", "jkm");

                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1, -1, -1);
                                                     jump.Clear();
                                                     jump.Acc(1.0, res);
                                                 }).Execute();



                // iterate over both cells (leading coeffcient)
                for(int i = 0; i <= 1; i++) {

                    // iterate over all conditions for the current cell
                    for(int k = 0; k * basisIndex < basisIndex + 1; k++) {

                        // iterate over the variables
                        for(int column = 0; column * basisIndex < basisIndex + 1; column++) {

                            aggE1[rowOffset + row * 4 + k + i * 2, index_i * 2 + column] = jump[i, k * basisIndex, 0, column * basisIndex];
                            aggE1[rowOffset + row * 4 + k + i * 2, index_j * 2 + column] = -jump[i, k * basisIndex, 1, column * basisIndex];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_i * (basisIndex + 1) + column] = gradient[i, k - 1, 0, column - 1];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_j * (basisIndex + 1) + column] = gradient[i, k - 1, 1, column - 1];
                        }
                    }
                }
                // test if rows are linear dependent
                usableRowsE.AddRange(CheckLinDependency(aggE1, 4, rowOffset + row * 4));
            }

            #endregion

            #region gradient jump conditions

            // smootheness
            // TODO hardest part: Evaluate the gradients on the edge in normal direction then compute the L2-norm of the jump and enforce a minimization condition
            // 1. compute the gradient parts for the relevant (0-basisIndex) basis functions in normal direction on each edge
            // 2. compute the Gram determinant for the edge
            // 3. integrate the basisfunction couplings (originating from the L2-norm) over the edge (Caching would be highly advisable, because the same integrals will be needed for the next basis function to aggregate)
            // 4. write values in aggEq

            // create array containing the edge-normal gradient coupling integrals
            MultidimensionalArray gradient = MultidimensionalArray.Create(2, basisIndex, 2, basisIndex);
            List<int> usableRowsQ = new List<int>();
            // iterate over edges
            for(int row = 0; row < iedges.Length; row++) {
                // FOR THE WHOLE SEGMENT, SUPER CLUNKY RIGHT NOW NEEDS TO BE WORKED OVER IF CODE IS FUNCTIONAL
                // select current edge
                int edge = iedges[row];

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[iedges[row], 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[iedges[row], 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // create array containing the edge-normal gradient coupling integrals
                gradient.Clear();
                //gradient = GetEdgeGradient(_maxDgBasis, _agGrd, edge, cell_i, cell_j, basisIndex);

                // compute the gradient coupled integrals
                var Context = _maxDgBasis.GridDat;
                var edgeMask = new EdgeMask(Context, new int[] { edge }, MaskType.Geometrical);

                EdgeQuadrature.GetQuadrature(new int[4] { 2, basisIndex, 2, basisIndex }, Context,
                    (new EdgeQuadratureScheme(true, edgeMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes_ref_edge = QR.Nodes;
                        //Debug.Assert(Length == 1);

                        MultidimensionalArray edge_normals = Context.iGeomEdges.NormalsCache.GetNormals_Edge(nodes_ref_edge, edge, 1);
                        Tuple<MultidimensionalArray, MultidimensionalArray> edge_gradient = _maxDgBasis.EdgeEvalGradient(nodes_ref_edge, edge, 1);

                        // apply the rotation
                        edge_gradient.Item1.Multiply(1.0, edge_gradient.Item1.CloneAs(), _rot[index_i], 0.0, "cijn", "cikn", "kj");
                        edge_gradient.Item2.Multiply(1.0, edge_gradient.Item2.CloneAs(), _rot[index_j], 0.0, "cijn", "cikn", "kj");

                        MultidimensionalArray edge_normal_gradient = MultidimensionalArray.Create(2, nodes_ref_edge.Length, basisIndex);

                        // compute the gradients in normal direction at the quadnodes reverse the normal for one cell
                        edge_normal_gradient.ExtractSubArrayShallow(0, -1, -1).Multiply(1.0, edge_gradient.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");
                        edge_normal_gradient.ExtractSubArrayShallow(1, -1, -1).Multiply(-1.0, edge_gradient.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");

                        // couple the gradients at each node (out "kinjm" "node_k,cell_i,polynom_n,cell_j,polynom_m")
                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1, -1, -1);
                        EvalResult.Multiply(1.0, edge_normal_gradient, edge_normal_gradient, 0.0, "kinjm", "ikn", "jkm");

                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1, -1, -1);
                                                     gradient.Clear();
                                                     gradient.Acc(1.0, res);
                                                 }).Execute();



                // iterate over both cells (leading coeffcient)
                for(int i = 0; i <= 1; i++) {

                    // iterate over all conditions for the current cell
                    for(int k = 1; k * basisIndex < basisIndex + 1; k++) {

                        // iterate over the variables
                        for(int column = 1; column * basisIndex < basisIndex + 1; column++) {
                            aggE2[row * 2 + (k - 1) + i * 1, index_i * 2 + column] = gradient[i, k * basisIndex - 1, 0, column * basisIndex - 1];
                            aggE2[row * 2 + (k - 1) + i * 1, index_j * 2 + column] = gradient[i, k * basisIndex - 1, 1, column * basisIndex - 1];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_i * (basisIndex + 1) + column] = gradient[i, k - 1, 0, column - 1];
                            //A[row * (basisIndex * 2) + (k - 1) + i * basisIndex, index_j * (basisIndex + 1) + column] = gradient[i, k - 1, 1, column - 1];
                        }
                    }
                }
                // test if rows are linear dependent
                usableRowsQ.AddRange(CheckLinDependency(aggE2, 2, row * 2));
            }

            #endregion

            #region solve system

            //Partitioning rowPartA = new Partitioning((basisIndex * 2) * iedges.Length);
            //Partitioning rowPartE = new Partitioning(basisIndex + 1 + (currentDegree + 1) * iedges.Length);
            //Partitioning colPart = new Partitioning(varCount);

            //// Matrix storing the optimization goal
            //MsrMatrix A = new MsrMatrix(rowPartA, colPart);
            //// Matrix storing the side conditions
            ////MsrMatrix E = new MsrMatrix(rowPartE, colPart);

            // construct the complete matrix for the quadratic program
            MultidimensionalArray E = MultidimensionalArray.Create(usableRowsE.Count(), varCount);
            MultidimensionalArray ET = MultidimensionalArray.Create(varCount, usableRowsE.Count());

            for(int i = 0; i < usableRowsE.Count(); i++) {
                E.SetRow(i, aggE1.ExtractSubArrayShallow(usableRowsE.ElementAt(i), -1).To1DArray());
                ET.SetColumn(i, aggE1.ExtractSubArrayShallow(usableRowsE.ElementAt(i), -1).To1DArray());
            }

            MultidimensionalArray hQ = MultidimensionalArray.Create(usableRowsQ.Count(), varCount);

            for(int i = 0; i < usableRowsQ.Count(); i++) {
                hQ.SetRow(i, aggE2.ExtractSubArrayShallow(usableRowsQ.ElementAt(i), -1).To1DArray());
            }

            // solution using quadratic program and sparse matrices
            MultidimensionalArray Q = MultidimensionalArray.Create(hQ.NoOfCols, hQ.NoOfCols);
            double Qscale = 1.0; // 1 / Math.Pow(aggE2.InfNorm(), 2.0);
            Q.GEMM(Qscale, hQ, hQ, 0.0, true);

            Partitioning qRowPart = new Partitioning(Q.NoOfRows + E.NoOfRows, MPI.Wrappers.csMPI.Raw._COMM.SELF);
            Partitioning qColPart = new Partitioning(Q.NoOfCols + E.NoOfRows, MPI.Wrappers.csMPI.Raw._COMM.SELF);
            MsrMatrix QP = new MsrMatrix(qRowPart, qColPart);

            double[] c = new double[varCount];
            double[] d = new double[E.NoOfRows];

            d[1] = 1;

            QP.AccBlock(0, 0, 1.0, Q);
            QP.AccBlock(0, Q.NoOfCols, 1.0, ET);
            QP.AccBlock(Q.NoOfRows, 0, 1.0, E);
            QP.AssumeSymmetric = true;


            double[] RHS = c.Concat(d).ToArray();
            double[] v = new double[QP.NoOfCols];


            #region Sparse Solve
            ////var solver = new ilPSP.LinSolvers.MUMPS.MUMPSSolver();
            //var solver = new ilPSP.LinSolvers.PARDISO.PARDISOSolver();

            //solver.DefineMatrix(QP);
            //solver.Solve(v, RHS);

            //MultidimensionalArray system = MultidimensionalArray.Create(usableRowsE.Count + usableRowsQ.Count, varCount);
            //system.ExtractSubArrayShallow(new int[] {0,0 }, new int[] { usableRowsE.Count - 1, varCount - 1 }).Acc(1.0, E);
            //system.ExtractSubArrayShallow(new int[] { usableRowsE.Count, 0 }, new int[] { usableRowsQ.Count + usableRowsE.Count - 1, varCount - 1 }).Acc(1.0,hQ);

            //MultidimensionalArray expsystem = MultidimensionalArray.Create(usableRowsQ.Count + usableRowsE.Count, usableRowsQ.Count + usableRowsE.Count);
            //expsystem.DGEMM(1.0, system, system, 0.0, false, true);

            //double[] systemrhs = new double[usableRowsE.Count + usableRowsQ.Count];
            //systemrhs[1] = 1;

            MultidimensionalArray system = MultidimensionalArray.Create(usableRowsE.Count + varCount, usableRowsE.Count + varCount);
            system.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { varCount - 1, varCount - 1 }).GEMM(1.0, hQ, hQ, 0.0, true);
            system.ExtractSubArrayShallow(new int[] { varCount, 0 }, new int[] { varCount + usableRowsE.Count - 1, varCount - 1 }).Acc(1.0, E);
            system.ExtractSubArrayShallow(new int[] { 0, varCount }, new int[] { varCount - 1, varCount + usableRowsE.Count - 1 }).Acc(1.0, E.TransposeTo());


            double[] systemrhs = new double[usableRowsE.Count + varCount];
            systemrhs[varCount + 1] = 1;

            // compute the solution that fullfills the orthogonality and continuity exactly and the minimizes the gradient jumps
            double[] sol = new double[varCount];
            //double[] expsol = new double[usableRowsQ.Count + usableRowsE.Count];

            //// classic solve
            //expsystem.Solve(expsol, systemrhs);
            //// sol = A'*expsol
            //for (int j = 0; j < varCount; j++)
            //{
            //    for (int i = 0; i < usableRowsE.Count + usableRowsQ.Count; i++)
            //    {
            //        sol[j] += system[i, j] * expsol[i];
            //    }
            //}

            //system.LeastSquareSolve(sol, systemrhs);
            sol = new double[varCount + usableRowsE.Count];
            system.LeastSquareSolve(sol, systemrhs);

            //aggE1.SolveWithCondition(sol, aggRhs1.To1DArray(), aggE2, aggRhs2.To1DArray());
            //sol = v.GetSubVector(0, varCount);


            // check sparsity
            //int nonzeroCount1 = 0;
            //foreach (double t in aggE1.To2DArray()) {
            //    nonzeroCount1 += t != 0.0 ? 1 : 0; 
            //}
            //int nonzeroCount2 = 0;
            //foreach (double t in aggE2.To2DArray())
            //{
            //    nonzeroCount2 += t != 0.0 ? 1 : 0;
            //}

            //double sparsity1 = (double)nonzeroCount1 / (aggE1.NoOfRows * aggE1.NoOfCols);
            //double sparsity2 =  (double)nonzeroCount2 / (aggE2.NoOfRows * aggE2.NoOfCols);
            //Console.WriteLine($"Sparsity of first equation set: {sparsity1}; Sparsity of second equation set: {sparsity2}");

            // compute the diagonal entry of the Injector for the first cell in parts and then solve the dependencies to the other cells
            // normality condition
            //double scale = 1 / sol.L2Norm();
            //sol.ScaleV(scale);
            MultidimensionalArray vector_InjectorNDegree = MultidimensionalArray.CreateWrapper(sol, varCount);
            #endregion


            //#region fmincon

            //// check the diagonal entries
            //bool isPD = true;
            //for (int i = 0; i < parts.Length; i++)
            //{
            //    if (Math.Abs(vector_InjectorNDegree[(i + 1) * (basisIndex + 1) - 1]) <= 1E-5) { isPD = false; }
            //}

            //// if diagonal entries are too small use fmincon
            //if (!isPD)
            //{
            //    Console.WriteLine("Using fmincon to ensure Positive Definiteness");
            //    // solve with fmincon as optimization problem, note the 3 additional scripts already present in the working directory
            //    BatchmodeConnector bmc = new BatchmodeConnector(@"D:\bosss_db_masterthesis\Matlab");
            //    bmc.Cmd("clear all");

            //    // gather the necessary matrices
            //    bmc.PutMatrix(Q, "Q");

            //    // quadratic constraints (note these are nonconvex) 0.5*x'H_ix + d_i <= 0
            //    MultidimensionalArray[] H = new MultidimensionalArray[parts.Length];
            //    double lowerlimit = Math.Sqrt((double)1 / parts.Length * 1E-8);
            //    string infimum = $@"1/{parts.Length}*1E-8";
            //    bmc.Cmd($"H = cell({parts.Length},1);");
            //    bmc.Cmd($"d = cell({parts.Length},1);");
            //    for (int i = 0; i < parts.Length; i++)
            //    {
            //        H[i] = MultidimensionalArray.Create(varCount, varCount);
            //        H[i][(i + 1) * (basisIndex + 1) - 1, (i + 1) * (basisIndex + 1) - 1] = -1;
            //        bmc.PutMatrix(H[i], $"H{{{i + 1}}}");
            //        bmc.Cmd($"d{{{i + 1}}} = {infimum};");
            //    }

            //    // equality constraint Gx = 0
            //    bmc.PutMatrix(E, "G");
            //    // remove the fix condition
            //    bmc.Cmd($"G({basisIndex + 1},:) = [];");
            //    bmc.Cmd($"G = sparse(G);");
            //    bmc.Cmd("b = sparse(size(G,1),1);");

            //    // ensure normality 1/2x'Heqx + d = 0
            //    bmc.Cmd("Heq = cell(1,1);");
            //    bmc.Cmd("deq = cell(1,1);");

            //    bmc.Cmd("Heq{1} = 2*eye(size(G,2));");
            //    bmc.Cmd("deq{1} = -1;");



            //    bmc.Cmd(@"options = optimoptions(@fmincon,'Algorithm','interior-point',...
            //        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
            //        'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H,Heq));");

            //    bmc.Cmd("fun = @(x)quadobj(x,Q);");
            //    bmc.Cmd("nonlconstr = @(x)quadconstr(x,H,d,Heq,deq);");
            //    bmc.Cmd("x0 = sparse(size(G,2),1);");

            //    bmc.Cmd("flin = sparse(size(x0,1),1);");
            //    bmc.Cmd("fnlin = @(x)quadobj(x,sparse(size(Q,1),size(Q,2)));");

            //    bmc.Cmd("x0 = linprog(flin,[],[],G,b,[],[]);");
            //    bmc.Cmd(@"x0 = fmincon(fnlin,x0,...
            //        [],[], G, b,[],[], nonlconstr, options); ");
            //    bmc.Cmd(@"[x,fval,eflag,output,lambda] = fmincon(fun,x0,...
            //        [],[], G, b,[],[], nonlconstr, options); ");


            //    MultidimensionalArray vector_InjectorNDegree_Matlab = MultidimensionalArray.Create(varCount, 1);
            //    bmc.GetMatrix(vector_InjectorNDegree_Matlab, "x");

            //    bmc.Execute(false);

            //    vector_InjectorNDegree = vector_InjectorNDegree_Matlab.ExtractSubArrayShallow(-1, 0);

            //    // simple test to at least get an idea if somethings gone horribly wrong
            //    if (Math.Abs(vector_InjectorNDegree_Matlab.L2Norm() - 1) >= 1E-3) Console.WriteLine("Warning converged to infeasable point");

            //    isPD = true;
            //    for (int i = 0; i < parts.Length; i++)
            //    {
            //        if (Math.Abs(vector_InjectorNDegree[(i + 1) * (basisIndex + 1) - 1]) - lowerlimit <= 0) { isPD = false; }
            //    }

            //    if (!isPD) Console.WriteLine($"WARNING: Basisfunction {basisIndex} violates pd condition");
            //}
            //#endregion

            #endregion



            // restructure to Injector operator form
            for(int i = 0; i < parts.Length; i++) {
                InjectorNDegree[i, 0] = vector_InjectorNDegree[i * 2];
                InjectorNDegree[i, basisIndex] = vector_InjectorNDegree[i * 2 + 1];
            }

            return InjectorNDegree;
        }

        // computes the injector coefficients for the first basis functions, successive variant
        private static double GetCoefficients0DegreeSuccessive(AggregationGridData _agGrd, Basis _maxDgBasis, MultidimensionalArray _inj_i, int[] parts) {
            if(parts.Length != 2) throw new ArgumentException();

            // calculate the aggregation cell volumes
            double aggVol_i = _agGrd.AncestorGrid.Cells.GetCellVolume(parts[0]);
            double aggVol_j = _agGrd.AncestorGrid.Cells.GetCellVolume(parts[1]);

            double Injector0Degree = Math.Sqrt(aggVol_j / aggVol_i) * _inj_i[0, 0];

            return Injector0Degree;
        }

        private static (MultidimensionalArray, double) GetCoefficientsSuccessive(AggregationGridData _agGrd, Basis _maxDgBasis, MultidimensionalArray _inj_i, MultidimensionalArray[] _rot, int[] parts, int edge, int basisIndex, bool direction) {

            if(parts.Length != 2) throw new ArgumentException();

            var InjectorNDegree = MultidimensionalArray.Create(1, basisIndex + 1);
            double quality;
            quality = Double.PositiveInfinity;

            #region startup

            int partCount = parts.Length;
            int currentDegree = _maxDgBasis.Polynomials[0][basisIndex].AbsoluteDegree;

            // equation matrix
            // basisIndex + 1 equations for othogonality
            // (currentdegree + 1) * iedges.Length equations for continuity along inneredges
            // (basisIndex + 1) * iedges.Length equations for edge-normal gradient jump minimization (smootheness)

            // stores continuity conditions
            MultidimensionalArray aggE1 = MultidimensionalArray.Create(currentDegree + 1, basisIndex + 1);
            // rhs
            MultidimensionalArray aggRhs1 = MultidimensionalArray.Create(currentDegree + 1);


            // stores gradient jump conditions
            MultidimensionalArray aggE2 = MultidimensionalArray.Create(basisIndex * 2, basisIndex + 1);
            MultidimensionalArray aggRhs2 = MultidimensionalArray.Create(basisIndex * 2);

            // array to ensure positive definiteness
            MultidimensionalArray f = MultidimensionalArray.Create(basisIndex + 1);


            #endregion            

            #region continuity conditions
            {
                // continuity
                // TODO harder part: Select currentDegree + 1 points per shared edge and ensure continuity there
                // 1. per edge construct currentDegree + 1 points
                // 2. evaluate the basisIndex -basisfunction on the adjacent cells for these points
                // 3. write the values in aggEq

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[edge, 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[edge, 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // get basis values on the edge
                // construct evaluation points, currently only viable for 1D edges                
                // create NodeSet
                int edgeRef = _agGrd.iGeomEdges.GetRefElementIndex(edge);
                var refEl = _agGrd.iGeomEdges.EdgeRefElements[edgeRef];
                int edgeVertCount = refEl.NoOfVertices;

                MultidimensionalArray edgeEvalPoints = MultidimensionalArray.Create(currentDegree + 1, refEl.SpatialDimension);
                refEl.Vertices.To2DArray();
                if(edgeVertCount < currentDegree + 1) {
                    for(int i = 0; i < edgeVertCount; i++) {
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(1.0, refEl.Vertices.ExtractSubArrayShallow(i, -1));
                    }

                    for(int i = edgeVertCount; i < currentDegree + 1; i++) {
                        // construct a point as middle point between the edge vertices
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(0.5, edgeEvalPoints.ExtractSubArrayShallow(i - edgeVertCount, -1));
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(0.5, edgeEvalPoints.ExtractSubArrayShallow(i - edgeVertCount + 1, -1));
                    }
                } else {
                    for(int i = 0; i < currentDegree + 1; i++) {
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(1.0, refEl.Vertices.ExtractSubArrayShallow(i, -1));
                    }
                }


                NodeSet edge_coords = new NodeSet(refEl, edgeEvalPoints, false);
                edge_coords.LockForever();

                // evaluate the polynomials on the edge for both cells
                Tuple<MultidimensionalArray, MultidimensionalArray> edge_values = _maxDgBasis.EdgeEval(edge_coords, edge, 1);

                // apply the rotation
                edge_values.Item1.Multiply(1.0, edge_values.Item1.CloneAs(), _rot[index_i], 0.0, "cij", "cik", "kj");
                edge_values.Item2.Multiply(1.0, edge_values.Item2.CloneAs(), _rot[index_j], 0.0, "cij", "cik", "kj");

                var edge_values_0 = index_i == 0 ? edge_values.Item1 : edge_values.Item2;
                var edge_values_1 = index_i == 0 ? edge_values.Item2 : edge_values.Item1;

                // also apply the injector for the currently fixed cell
                //edge_values_0.Multiply(1.0, edge_values_0.CloneAs(), _inj_i, 0.0, "cij", "cik", "kj");

                for(int i = 0; i < (currentDegree + 1); i++) {

                    for(int column = 0; column < basisIndex + 1; column++) {
                        aggRhs1[i] += edge_values_0[0, i, column] * _inj_i[column, basisIndex];
                        aggE1[i, column] = edge_values_1[0, i, column];
                    }
                }
            }
            #endregion

            #region gradient jump conditions
            {
                // smootheness
                // TODO hardest part: Evaluate the gradients on the edge in normal direction then compute the L2-norm of the jump and enforce a minimization condition
                // 1. compute the gradient parts for the relevant (0-basisIndex) basis functions in normal direction on each edge
                // 2. compute the Gram determinant for the edge
                // 3. integrate the basisfunction couplings (originating from the L2-norm) over the edge (Caching would be highly advisable, because the same integrals will be needed for the next basis function to aggregate)
                // 4. write values in aggEq

                // create array containing the edge-normal gradient coupling integrals
                MultidimensionalArray gradient = MultidimensionalArray.Create(2, basisIndex, 2, basisIndex);

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[edge, 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[edge, 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // compute the gradient coupled integrals
                var Context = _maxDgBasis.GridDat;
                var edgeMask = new EdgeMask(Context, new int[] { edge }, MaskType.Geometrical);

                EdgeQuadrature.GetQuadrature(new int[4] { 2, basisIndex, 2, basisIndex }, Context,
                    (new EdgeQuadratureScheme(true, edgeMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes_ref_edge = QR.Nodes;
                        //Debug.Assert(Length == 1);

                        MultidimensionalArray edge_normals = Context.iGeomEdges.NormalsCache.GetNormals_Edge(nodes_ref_edge, edge, 1);
                        Tuple<MultidimensionalArray, MultidimensionalArray> edge_gradient = _maxDgBasis.EdgeEvalGradient(nodes_ref_edge, edge, 1);

                        // apply the rotation
                        edge_gradient.Item1.Multiply(1.0, edge_gradient.Item1.CloneAs(), _rot[index_i], 0.0, "cijn", "cikn", "kj");
                        edge_gradient.Item2.Multiply(1.0, edge_gradient.Item2.CloneAs(), _rot[index_j], 0.0, "cijn", "cikn", "kj");

                        //// apply the injector for the fixed cell
                        //if (index_i == 0)
                        //{
                        //    edge_gradient.Item1.Multiply(1.0, edge_gradient.Item1.CloneAs(), _inj_i, 0.0, "cijn", "cikn", "kj");
                        //}
                        //else
                        //{
                        //    edge_gradient.Item2.Multiply(1.0, edge_gradient.Item2.CloneAs(), _inj_i, 0.0, "cijn", "cikn", "kj");
                        //}

                        MultidimensionalArray edge_normal_gradient = MultidimensionalArray.Create(2, nodes_ref_edge.Length, basisIndex);

                        // compute the gradients in normal direction at the quadnodes reverse the normal for one cell
                        edge_normal_gradient.ExtractSubArrayShallow(0, -1, -1).Multiply(1.0, edge_gradient.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");
                        edge_normal_gradient.ExtractSubArrayShallow(1, -1, -1).Multiply(-1.0, edge_gradient.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");

                        // couple the gradients at each node (out "kinjm" "node_k,cell_i,polynom_n,cell_j,polynom_m")
                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1, -1, -1);
                        EvalResult.Multiply(1.0, edge_normal_gradient, edge_normal_gradient, 0.0, "kinjm", "ikn", "jkm");

                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1, -1, -1);
                                                     gradient.Clear();
                                                     gradient.Acc(1.0, res);
                                                 }).Execute();



                var gradient_values_0 = index_i == 0 ? gradient.ExtractSubArrayShallow(-1, -1, 0, -1) : gradient.ExtractSubArrayShallow(-1, -1, 1, -1);
                var gradient_values_1 = index_i == 0 ? gradient.ExtractSubArrayShallow(-1, -1, 1, -1) : gradient.ExtractSubArrayShallow(-1, -1, 0, -1);

                // iterate over both cells (leading coeffcient)
                for(int i = 0; i <= 1; i++) {

                    // iterate over all conditions for the current cell
                    for(int k = 1; k < basisIndex + 1; k++) {

                        // iterate over the variables
                        for(int column = 1; column < basisIndex + 1; column++) {
                            aggE2[(k - 1) + i * basisIndex, column] = gradient_values_1[i, k - 1, column - 1];
                            aggRhs2[(k - 1) + i * basisIndex] -= gradient_values_0[i, k - 1, column - 1] * _inj_i[column, basisIndex];

                        }
                    }
                }

            }
            #endregion

            #region solve system

            // we solve the problem
            // minimize     F = x'*(aggE2'*aggE2)*x + f'*x - aggRhs2
            // s.t.         aggE1*x - aggRhs1 = 0

            MultidimensionalArray Q = MultidimensionalArray.Create(basisIndex + 1, basisIndex + 1);
            Q.GEMM(1.0, aggE2, aggE2, 0.0, true);
            //Q[basisIndex, basisIndex] += fscale;

            // accumulate system matrix
            MultidimensionalArray S = MultidimensionalArray.Create(basisIndex + 1 + currentDegree + 1, basisIndex + 1 + currentDegree + 1);
            S.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { basisIndex, basisIndex }).Acc(1.0, Q);
            S.ExtractSubArrayShallow(new int[] { 0, basisIndex + 1 }, new int[] { basisIndex, basisIndex + 1 + currentDegree }).Acc(1.0, aggE1.TransposeTo());
            S.ExtractSubArrayShallow(new int[] { basisIndex + 1, 0 }, new int[] { basisIndex + 1 + currentDegree, basisIndex }).Acc(1.0, aggE1);

            // accumulate rhs
            MultidimensionalArray rS = MultidimensionalArray.Create(basisIndex + 1 + currentDegree + 1);
            rS.ExtractSubArrayShallow(new int[] { 0 }, new int[] { basisIndex }).Acc(-1.0, f);
            rS.ExtractSubArrayShallow(new int[] { basisIndex + 1 }, new int[] { basisIndex + 1 + currentDegree }).Acc(1.0, aggRhs1);

            // solution vector
            double[] x = new double[rS.Length];
            double factor = Math.Sqrt(_agGrd.AncestorGrid.Cells.GetCellVolume(parts[1]) / _agGrd.AncestorGrid.Cells.GetCellVolume(parts[0]));
            rS[basisIndex] += direction ? -factor * S[basisIndex, basisIndex] : factor * S[basisIndex, basisIndex];
            S.LeastSquareSolve(x, rS.To1DArray());

            double[] sol = x.GetSubVector(0, basisIndex + 1);


            #region ensure positive definiteness

            while(sol[basisIndex].Abs() <= 1E-1 * factor || sol[basisIndex].Abs() >= 1E+4) {
                double fscale = 1E+2 * Math.Sign(S[basisIndex, basisIndex]);//aggRhs2.L2Norm();
                rS[basisIndex] += direction ? -factor * fscale : factor * fscale;
                S[basisIndex, basisIndex] += fscale;

                Console.WriteLine($"Ensuring spd");

                S.LeastSquareSolve(x, rS.To1DArray());

                sol = x.GetSubVector(0, basisIndex + 1);
            }

            #endregion


            MultidimensionalArray vector_InjectorNDegree = MultidimensionalArray.CreateWrapper(sol, basisIndex + 1);

            #endregion
            MultidimensionalArray temp = MultidimensionalArray.Create(basisIndex + 1);
            temp.Multiply(1.0, Q, vector_InjectorNDegree, 0.0, "i", "ij", "j");

            quality = vector_InjectorNDegree.InnerProduct(temp) - rS.ExtractSubArrayShallow(new int[] { 0 }, new int[] { basisIndex }).InnerProduct(vector_InjectorNDegree) - aggRhs2.InnerProduct(aggRhs2.CloneAs());

            return (vector_InjectorNDegree, quality);
        }

        private static (MultidimensionalArray, double) GetCoefficientsScaleFit(AggregationGridData _agGrd, Basis _maxDgBasis, MultidimensionalArray _inj_i, MultidimensionalArray[] _rot, int[] parts, int edge, int basisIndex, bool direction) {

            if(parts.Length != 2) throw new ArgumentException();

            var InjectorNDegree = MultidimensionalArray.Create(1, basisIndex + 1);
            double quality;
            quality = Double.PositiveInfinity;

            #region startup

            int partCount = parts.Length;
            int currentDegree = _maxDgBasis.Polynomials[0][basisIndex].AbsoluteDegree;

            // equation matrix
            // basisIndex + 1 equations for othogonality
            // (currentdegree + 1) * iedges.Length equations for continuity along inneredges
            // (basisIndex + 1) * iedges.Length equations for edge-normal gradient jump minimization (smootheness)

            // stores continuity conditions
            MultidimensionalArray aggE1 = MultidimensionalArray.Create(currentDegree + 1, basisIndex + 1);
            // rhs
            MultidimensionalArray aggRhs1 = MultidimensionalArray.Create(currentDegree + 1);


            // stores gradient jump conditions
            MultidimensionalArray aggE2 = MultidimensionalArray.Create(basisIndex * 2, basisIndex + 1);
            MultidimensionalArray aggRhs2 = MultidimensionalArray.Create(basisIndex * 2);

            // array to ensure positive definiteness
            MultidimensionalArray f = MultidimensionalArray.Create(basisIndex + 1);


            #endregion            

            #region continuity conditions
            {
                // continuity
                // TODO harder part: Select currentDegree + 1 points per shared edge and ensure continuity there
                // 1. per edge construct currentDegree + 1 points
                // 2. evaluate the basisIndex -basisfunction on the adjacent cells for these points
                // 3. write the values in aggEq

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[edge, 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[edge, 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // get basis values on the edge
                // construct evaluation points, currently only viable for 1D edges                
                // create NodeSet
                int edgeRef = _agGrd.iGeomEdges.GetRefElementIndex(edge);
                var refEl = _agGrd.iGeomEdges.EdgeRefElements[edgeRef];
                int edgeVertCount = refEl.NoOfVertices;

                MultidimensionalArray edgeEvalPoints = MultidimensionalArray.Create(currentDegree + 1, refEl.SpatialDimension);
                refEl.Vertices.To2DArray();
                if(edgeVertCount < currentDegree + 1) {
                    for(int i = 0; i < edgeVertCount; i++) {
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(1.0, refEl.Vertices.ExtractSubArrayShallow(i, -1));
                    }

                    for(int i = edgeVertCount; i < currentDegree + 1; i++) {
                        // construct a point as middle point between the edge vertices
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(0.5, edgeEvalPoints.ExtractSubArrayShallow(i - edgeVertCount, -1));
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(0.5, edgeEvalPoints.ExtractSubArrayShallow(i - edgeVertCount + 1, -1));
                    }
                } else {
                    for(int i = 0; i < currentDegree + 1; i++) {
                        edgeEvalPoints.ExtractSubArrayShallow(i, -1).Acc(1.0, refEl.Vertices.ExtractSubArrayShallow(i, -1));
                    }
                }


                NodeSet edge_coords = new NodeSet(refEl, edgeEvalPoints, false);
                edge_coords.LockForever();

                // evaluate the polynomials on the edge for both cells
                Tuple<MultidimensionalArray, MultidimensionalArray> edge_values = _maxDgBasis.EdgeEval(edge_coords, edge, 1);

                // apply the rotation
                edge_values.Item1.Multiply(1.0, edge_values.Item1.CloneAs(), _rot[index_i], 0.0, "cij", "cik", "kj");
                edge_values.Item2.Multiply(1.0, edge_values.Item2.CloneAs(), _rot[index_j], 0.0, "cij", "cik", "kj");

                var edge_values_0 = index_i == 0 ? edge_values.Item1 : edge_values.Item2;
                var edge_values_1 = index_i == 0 ? edge_values.Item2 : edge_values.Item1;

                // also apply the injector for the currently fixed cell
                //edge_values_0.Multiply(1.0, edge_values_0.CloneAs(), _inj_i, 0.0, "cij", "cik", "kj");

                for(int i = 0; i < (currentDegree + 1); i++) {

                    for(int column = 0; column < basisIndex + 1; column++) {
                        aggRhs1[i] += edge_values_0[0, i, column] * _inj_i[column, basisIndex];
                        aggE1[i, column] = edge_values_1[0, i, column];
                    }
                }
            }
            #endregion

            #region gradient jump conditions
            {
                // smootheness
                // TODO hardest part: Evaluate the gradients on the edge in normal direction then compute the L2-norm of the jump and enforce a minimization condition
                // 1. compute the gradient parts for the relevant (0-basisIndex) basis functions in normal direction on each edge
                // 2. compute the Gram determinant for the edge
                // 3. integrate the basisfunction couplings (originating from the L2-norm) over the edge (Caching would be highly advisable, because the same integrals will be needed for the next basis function to aggregate)
                // 4. write values in aggEq

                // create array containing the edge-normal gradient coupling integrals
                MultidimensionalArray gradient = MultidimensionalArray.Create(2, basisIndex, 2, basisIndex);

                // get the cells for the edge
                int cell_i = _agGrd.iGeomEdges.CellIndices[edge, 0];
                int cell_j = _agGrd.iGeomEdges.CellIndices[edge, 1];

                // get indices of left and right hand cell on the respective edge
                int index_i = Array.IndexOf(parts, cell_i);
                int index_j = Array.IndexOf(parts, cell_j);

                // compute the gradient coupled integrals
                var Context = _maxDgBasis.GridDat;
                var edgeMask = new EdgeMask(Context, new int[] { edge }, MaskType.Geometrical);

                EdgeQuadrature.GetQuadrature(new int[4] { 2, basisIndex, 2, basisIndex }, Context,
                    (new EdgeQuadratureScheme(true, edgeMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                        NodeSet nodes_ref_edge = QR.Nodes;
                        //Debug.Assert(Length == 1);

                        MultidimensionalArray edge_normals = Context.iGeomEdges.NormalsCache.GetNormals_Edge(nodes_ref_edge, edge, 1);
                        Tuple<MultidimensionalArray, MultidimensionalArray> edge_gradient = _maxDgBasis.EdgeEvalGradient(nodes_ref_edge, edge, 1);

                        // apply the rotation
                        edge_gradient.Item1.Multiply(1.0, edge_gradient.Item1.CloneAs(), _rot[index_i], 0.0, "cijn", "cikn", "kj");
                        edge_gradient.Item2.Multiply(1.0, edge_gradient.Item2.CloneAs(), _rot[index_j], 0.0, "cijn", "cikn", "kj");

                        //// apply the injector for the fixed cell
                        //if (index_i == 0)
                        //{
                        //    edge_gradient.Item1.Multiply(1.0, edge_gradient.Item1.CloneAs(), _inj_i, 0.0, "cijn", "cikn", "kj");
                        //}
                        //else
                        //{
                        //    edge_gradient.Item2.Multiply(1.0, edge_gradient.Item2.CloneAs(), _inj_i, 0.0, "cijn", "cikn", "kj");
                        //}

                        MultidimensionalArray edge_normal_gradient = MultidimensionalArray.Create(2, nodes_ref_edge.Length, basisIndex);

                        // compute the gradients in normal direction at the quadnodes reverse the normal for one cell
                        edge_normal_gradient.ExtractSubArrayShallow(0, -1, -1).Multiply(1.0, edge_gradient.Item1.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");
                        edge_normal_gradient.ExtractSubArrayShallow(1, -1, -1).Multiply(-1.0, edge_gradient.Item2.ExtractSubArrayShallow(new int[] { 0, 0, 1, 0 }, new int[] { -1, nodes_ref_edge.Length - 1, basisIndex, Context.SpatialDimension - 1 }), edge_normals.ExtractSubArrayShallow(0, -1, -1), 0.0, "ik", "ikn", "in");

                        // couple the gradients at each node (out "kinjm" "node_k,cell_i,polynom_n,cell_j,polynom_m")
                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1, -1, -1);
                        EvalResult.Multiply(1.0, edge_normal_gradient, edge_normal_gradient, 0.0, "kinjm", "ikn", "jkm");

                    },
                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                     Debug.Assert(Length == 1);

                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1, -1, -1);
                                                     gradient.Clear();
                                                     gradient.Acc(1.0, res);
                                                 }).Execute();



                var gradient_values_0 = index_i == 0 ? gradient.ExtractSubArrayShallow(-1, -1, 0, -1) : gradient.ExtractSubArrayShallow(-1, -1, 1, -1);
                var gradient_values_1 = index_i == 0 ? gradient.ExtractSubArrayShallow(-1, -1, 1, -1) : gradient.ExtractSubArrayShallow(-1, -1, 0, -1);

                // iterate over both cells (leading coeffcient)
                for(int i = 0; i <= 1; i++) {

                    // iterate over all conditions for the current cell
                    for(int k = 1; k < basisIndex + 1; k++) {

                        // iterate over the variables
                        for(int column = 1; column < basisIndex + 1; column++) {
                            aggE2[(k - 1) + i * basisIndex, column] = gradient_values_1[i, k - 1, column - 1];
                            aggRhs2[(k - 1) + i * basisIndex] -= gradient_values_0[i, k - 1, column - 1] * _inj_i[column, basisIndex];

                        }
                    }
                }

            }
            #endregion

            #region solve system


            // calculate the aggregation cell volumes
            double aggVol_i = _agGrd.AncestorGrid.Cells.GetCellVolume(parts[0]);
            double aggVol_j = _agGrd.AncestorGrid.Cells.GetCellVolume(parts[1]);

            double Injector0Degree = Math.Sqrt(aggVol_j / aggVol_i) * _inj_i[0, 0];

            MultidimensionalArray Q = MultidimensionalArray.Create(2, 2);
            Q[0, 1] = direction ? -1 : 1;
            Q[1, 0] = aggE1[0, 0];
            Q[1, 1] = aggE1[0, basisIndex];

            // accumulate rhs
            MultidimensionalArray rS = MultidimensionalArray.Create(2);
            rS[0] = Injector0Degree;
            rS[1] = aggRhs1[0];

            // solution vector
            double[] x = new double[rS.Length];
            Q.Solve(x, rS.To1DArray());

            MultidimensionalArray vector_InjectorNDegree = MultidimensionalArray.Create(basisIndex + 1);
            vector_InjectorNDegree[0] = x[0];
            vector_InjectorNDegree[basisIndex] = x[1];

            #endregion
            MultidimensionalArray S = MultidimensionalArray.Create(basisIndex + 1, basisIndex + 1);
            S.GEMM(1.0, aggE2, aggE2, 0.0, true);

            MultidimensionalArray temp = MultidimensionalArray.Create(basisIndex + 1);
            temp.Multiply(1.0, S, vector_InjectorNDegree, 0.0, "i", "ij", "j");

            quality = vector_InjectorNDegree.InnerProduct(temp) - aggRhs2.InnerProduct(aggRhs2.CloneAs());

            return (vector_InjectorNDegree, quality);
        }

        // recursively extracts the Injection operator from ilevel-1 to  ilevel from the direct operator to ilevel
        public static void ExtractInjectorCurved(AggregationGridData[] _agGrd, Basis _maxDgBasis, MultidimensionalArray[][] _Injectors, MultidimensionalArray[] _injectorCoarse, int ilevel) {
            if(ilevel >= 1) {

                Console.WriteLine($"Extracting Injector from level {ilevel - 1} to {ilevel}");

                // get dimensions of Injection operator
                int Jagg = _agGrd[ilevel].iLogicalCells.NoOfLocalUpdatedCells;
                int Np = _maxDgBasis.Length;
                var Context = _maxDgBasis.GridDat;


                _Injectors[ilevel] = new MultidimensionalArray[Jagg];

                MultidimensionalArray ortho = MultidimensionalArray.Create(Np, Np);
                MultidimensionalArray orthoInv = MultidimensionalArray.Create(Np, Np);

                // compute the mass matrix for lCell on ilevel - 1 when using the Injector to ilevel
                MultidimensionalArray MMtemp = MultidimensionalArray.Create(Np, Np);
                MultidimensionalArray MMcontrol = MultidimensionalArray.Create(Np, Np);

                // iterate over all aggregate (logical) cells of ilevel
                for(int aggCell = 0; aggCell < Jagg; aggCell++) {

                    // logical cells on ilevel-1 that build the aggregated cell on ilevel
                    int[] aggParts = _agGrd[ilevel].jCellCoarse2jCellFine[aggCell];
                    _Injectors[ilevel][aggCell] = MultidimensionalArray.Create(aggParts.Length, Np, Np);

                    // iterate over all aggregate (logical) cells of ilevel - 1
                    foreach(int lCell in aggParts) {
                        int[] parts = _agGrd[ilevel - 1].iLogicalCells.AggregateCellToParts[lCell];

                        ortho.Clear();
                        orthoInv.Clear();
                        MMtemp.Clear();

                        // based on the cell ordering we have to distinguish different cases
                        // case cell on level k consists of more than one k-1 cells
                        if(aggParts.Length > 1) {

                            // case cell on level k-1 consists of more than one level 0 cells
                            if(parts.Length > 1) {

                                // assemble the Massmatrix of the finer level logical cell
                                foreach(int gCell in parts) {
                                    var cellMask = new CellMask(Context, new[] { new Chunk() { i0 = gCell, Len = 1 } }, MaskType.Geometrical);

                                    CellQuadrature.GetQuadrature(new int[2] { Np, Np }, Context,
                                    (new CellQuadratureScheme(true, cellMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                                        // get a set of quad points on the reference element
                                        NodeSet nodes = QR.Nodes;
                                        Debug.Assert(Length == 1);

                                        var phi_0 = _maxDgBasis.CellEval(nodes, gCell, 1).ExtractSubArrayShallow(0, -1, -1);
                                        // apply the coarse injector
                                        MultidimensionalArray phi = phi_0.CloneAs();
                                        phi.Multiply(1.0, phi_0, _injectorCoarse[gCell], 0.0, "ki", "kj", "ji");

                                        var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1);

                                        EvalResult.Multiply(1.0, phi, phi, 0.0, "kmn", "kn", "km");
                                    },
                                    /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                                                                     Debug.Assert(Length == 1);

                                                                     var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1);
                                                                     MMtemp.Acc(1.0, res);
                                                                 }).Execute();

                                    //MMtemp.DGEMM(1.0, _injectorCoarse[gCell], _injectorCoarse[gCell], 1.0, true);
                                }

                                // orthonormalize
                                //MMtemp.SymmetricLDLInversion(orthoInv, null);
                                MMtemp.Cholesky();
                                MMtemp.TransposeTo(ortho);
                                ortho.InvertTo(orthoInv);

                                // form the new direct injector to ilevel - 1
                                foreach(int gCell in parts) {
                                    MultidimensionalArray coarseClone = _injectorCoarse[gCell].CloneAs();
                                    _injectorCoarse[gCell].Multiply(1.0, coarseClone, orthoInv, 0.0, "ij", "ik", "kj");
                                }
                            }
                            // case cell on level k-1 consists of one level 0 cell
                            else {
                                ortho.Acc(1.0, _injectorCoarse[parts[0]]);
                                _injectorCoarse[parts[0]].Clear();
                                _injectorCoarse[parts[0]].AccEye(1.0);
                            }
                        }
                        // case cell on level k consists of one k-1 cell
                        else {
                            ortho.AccEye(1.0);
                        }

                        //orthoInv.InvertTo(ortho);
                        _Injectors[ilevel][aggCell].ExtractSubArrayShallow(Array.IndexOf(aggParts, lCell), -1, -1).Acc(1.0, ortho);



#if DEBUG
                        MMtemp.Clear();
                        MMcontrol.Clear();

                        //check orthonormality(move to a more suitable location later on)
                        foreach (int gCell in parts)
                        {

                            var cellMask = new CellMask(Context, new[] { new Chunk() { i0 = gCell, Len = 1 } }, MaskType.Geometrical);

                            CellQuadrature.GetQuadrature(new int[2] { Np, Np }, Context,
                            (new CellQuadratureScheme(true, cellMask)).Compile(Context, _maxDgBasis.Degree * 2), // integrate over target cell
                            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray _EvalResult)
                            {
                                // get a set of quad points on the reference element
                                NodeSet nodes = QR.Nodes;
                                Debug.Assert(Length == 1);

                                var phi_0 = _maxDgBasis.CellEval(nodes, gCell, 1).ExtractSubArrayShallow(0, -1, -1);
                                // apply the coarse injector
                                MultidimensionalArray phi = phi_0.CloneAs();
                                phi.Multiply(1.0, phi_0, _injectorCoarse[gCell], 0.0, "ki", "kj", "ji");

                                var EvalResult = _EvalResult.ExtractSubArrayShallow(0, -1, -1, -1);

                                EvalResult.Multiply(1.0, phi, phi, 0.0, "kmn", "kn", "km");
                            },
                            /*_SaveIntegrationResults:*/ delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration)
                                                         {
                                                             Debug.Assert(Length == 1);

                                                             var res = ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1);
                                                             MMtemp.Acc(1.0, res);
                                                         }).Execute();
                        }

                        foreach (int gCell in parts)
                        {
                            MMcontrol.GEMM(1.0, _injectorCoarse[gCell], _injectorCoarse[gCell], 1.0, true);
                        }


                        // calculate distance to expected unity matrix
                        MultidimensionalArray Eye = MultidimensionalArray.Create(Np, Np);
                        Eye.AccEye(1.0);
                        Eye.Acc(-1.0, MMtemp);
                        double dist = Eye.L2Norm();

                        MMcontrol.Acc(-1.0, MMtemp);
                        double distC = MMcontrol.L2Norm();
                        //Debug.Assert(dist < 1E-3);
                        Console.ForegroundColor = dist < 1E-3 ? ConsoleColor.White : ConsoleColor.Red;
                        Console.WriteLine($"distance between Eye and MM for level {ilevel - 1} and lCell {lCell}: {dist}", Console.ForegroundColor);
                        Console.ForegroundColor = ConsoleColor.White;
                        Console.ForegroundColor = distC < 1E-3 ? ConsoleColor.White : ConsoleColor.Red;
                        Console.WriteLine($"distance between MM inj and MM int for level {ilevel - 1} and lCell {lCell}: {distC}", Console.ForegroundColor);
                        Console.ForegroundColor = ConsoleColor.White;
#endif

                    }
                }

                Console.WriteLine("Done. ");

                // recursive repeat until ilevel = 1
                ExtractInjectorCurved(_agGrd, _maxDgBasis, _Injectors, _injectorCoarse, ilevel - 1);
            }

        }

        /// returns the inices of rows in Matrix <paramref name="E"/>
        /// that are linearly independent
        private static List<int> CheckLinDependency(MultidimensionalArray E, int rows, int offset) {
            int cols = E.NoOfCols;
            List<int> usableRows = new List<int>();

            int[] pivot = IMatrixExtensions.ReducedRowEchelonForm(E.ExtractSubArrayShallow(new int[] { offset, 0 }, new int[] { offset + rows - 1, cols - 1 }).TransposeTo()).Item2;

            for(int i = 0; i < pivot.Length; i++) pivot[i] += offset;

            usableRows.AddRange(pivot);
            return usableRows;
        }

        /// <summary>
        /// Calculate an Injection operator by projection of the basis functions from the curved cells on <paramref name="ilevel"/>-1
        /// to (linear) bounding boxes around the aggregated cells on <paramref name="ilevel"/>. Followed by an reorthonormalization on the aggregated cells.
        /// </summary>
        /// <param name="_agGrd"></param>
        /// <param name="_maxDgBasis"></param>
        /// <param name="_Injectors"></param> the Injection operator from (level-1 to level) with 1st index level; 2nd index logical (aggregate) cell
        /// <param name="_injectorCoarse"></param> the direct injector from level 0 to ilevel
        /// <param name="ilevel"></param>
        public static void ProjectBasis(AggregationGridData[] _agGrd, Basis _maxDgBasis, MultidimensionalArray[][] _Injectors, MultidimensionalArray[] _injectorCoarse, int ilevel) {
            // get dimensions of Injection operator
            int Jagg = _agGrd[ilevel].iLogicalCells.NoOfLocalUpdatedCells;
            int Np = _maxDgBasis.Length;
            var Context = _maxDgBasis.GridDat;
            // special treatment of 0-th level
            if(ilevel == 0) {
                for(int i = 0; i < Jagg; i++) {
                    //_Injectors[ilevel][i] = MultidimensionalArray.Create(Np, Np);
                    //_Injectors[ilevel][i].AccEye(1.0);
                    _injectorCoarse[i] = MultidimensionalArray.Create(Np, Np);
                    _injectorCoarse[i].AccEye(1.0);
                }
            } else {
                _Injectors[ilevel] = new MultidimensionalArray[Jagg];

                for(int i = 0; i < Jagg; i++) {


                    // number of logical cells from previous level
                    int[] lparts = _agGrd[ilevel].jCellCoarse2jCellFine[i];
                    int Jlparts = lparts.Length;

                    // extract Bounding Box
                    Platform.Utils.Geom.BoundingBox BB = new Platform.Utils.Geom.BoundingBox(Context.SpatialDimension);
                    _agGrd[ilevel].iLogicalCells.GetCellBoundingBox(i, BB);

                    // get Polynomials
                    PolynomialList polyList = _maxDgBasis.Polynomials[0];

                    MultidimensionalArray a = MultidimensionalArray.Create(Jlparts, Np, Np);
                    _Injectors[ilevel][i] = MultidimensionalArray.Create(Jlparts, Np, Np);
                    // iterate over logical parts
                    for(int k = 0; k < Jlparts; k++) {
                        // number of parts
                        int[] parts = _agGrd[ilevel - 1].iLogicalCells.AggregateCellToParts[lparts[k]];
                        int Jparts = parts.Length;

                        // project onto the aggregated cell (consisting of aggregated cells from the previous level)
                        for(int j = 0; j < Jparts; j++) {
                            var cellMask = new CellMask(Context, new[] { new Chunk() { i0 = parts[j], Len = 1 } }, MaskType.Geometrical);

                            CellQuadrature.GetQuadrature(new int[2] { Np, Np }, Context,
                                (new CellQuadratureScheme(true, cellMask)).Compile(Context, _maxDgBasis.Degree * 2),
                                delegate (int j0, int Length, QuadRule QR, MultidimensionalArray _EvalResult) {
                                    NodeSet nodes = QR.Nodes;
                                    int D = nodes.SpatialDimension;
                                    NodeSet GlobalNodes = new NodeSet(Context.iGeomCells.GetRefElement(parts[j]), Length * nodes.NoOfNodes, D, false);
                                    Context.TransformLocal2Global(nodes, j0, Length, GlobalNodes.ResizeShallow(Length, nodes.NoOfNodes, D));

                                    for(int d = 0; d < D; d++) {
                                        var Cd = GlobalNodes.ExtractSubArrayShallow(-1, d);
                                        Cd.Scale(2 / (BB.Max[d] - BB.Min[d]));
                                        Cd.AccConstant((BB.Max[d] + BB.Min[d]) / (BB.Min[d] - BB.Max[d]));
                                    }
#if DEBUG
                                Debug.Assert(GlobalNodes.Min() >= -1.00001);
                                Debug.Assert(GlobalNodes.Max() <= +1.00001);
#endif
                                    GlobalNodes.LockForever();

                                    var BasisValuesRaw = _maxDgBasis.CellEval(nodes, j0, Length);
                                    // apply previous injection
                                    MultidimensionalArray BasisValues = MultidimensionalArray.Create(Length, nodes.NoOfNodes, Np);
                                    BasisValues.Multiply(1.0, BasisValuesRaw, _injectorCoarse[parts[j]], 0.0, "jkn", "jkm", "mn");
                                    var PolyVals = MultidimensionalArray.Create(GlobalNodes.NoOfNodes, Np);
                                    polyList.Evaluate(GlobalNodes, PolyVals);
                                    PolyVals = PolyVals.ResizeShallow(Length, nodes.NoOfNodes, Np);

                                    _EvalResult.Multiply(1.0, BasisValues, PolyVals, 0.0, "jknm", "jkn", "jkm");
                                },
                                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // _SaveIntegrationResults:
                                    a.ExtractSubArrayShallow(k, -1, -1).Acc(1.0, ResultsOfIntegration.ExtractSubArrayShallow(0, -1, -1));
                                }).Execute();
                        }
                    }

                    // reorthonormalize
                    MultidimensionalArray MM = MultidimensionalArray.Create(Np, Np);
                    MultidimensionalArray U = MultidimensionalArray.Create(Np, Np);
                    MultidimensionalArray orthoInv = MultidimensionalArray.Create(Np, Np);
                    for(int k = 0; k < Jlparts; k++) {
                        MM.GEMM(1.0, a.ExtractSubArrayShallow(k, -1, -1), a.ExtractSubArrayShallow(k, -1, -1), 1.0, true);
                    }

                    MM.Cholesky();
                    MM.TransposeTo(U);
                    U.InvertTo(orthoInv);

                    // populate Injectors
                    _Injectors[ilevel][i].Multiply(1.0, a, orthoInv, 0.0, "jkn", "jkm", "mn");
                    // next direct Injector
                    for(int k = 0; k < _Injectors[ilevel][i].Lengths[0]; k++) {
                        // number of parts
                        int[] parts = _agGrd[ilevel - 1].iLogicalCells.AggregateCellToParts[lparts[k]];
                        int Jparts = parts.Length;

                        for(int j = 0; j < Jparts; j++) {
                            int kCell = parts[j];
                            _injectorCoarse[kCell].GEMM(1.0, _injectorCoarse[kCell].CloneAs(), _Injectors[ilevel][i].ExtractSubArrayShallow(k, -1, -1), 0.0);
                        }

                    }
                }
            }

            // procede to next grid level
            ilevel++;
            if(ilevel < _agGrd.Length) {
                ProjectBasis(_agGrd, _maxDgBasis, _Injectors, _injectorCoarse, ilevel);
            }
        }

    }
}
